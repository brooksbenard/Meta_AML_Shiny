# =============================================================================
# Load and process BeatAML2 data for drug sensitivity analysis
# Data from https://biodev.github.io/BeatAML2/
# =============================================================================

load_beataml2 <- function(data_dir = file.path(getwd(), "beataml2_data")) {
  mutations_path <- file.path(data_dir, "mutations.txt")
  auc_path <- file.path(data_dir, "inhibitor_auc.txt")
  clinical_path <- file.path(data_dir, "clinical.xlsx")

  if (!file.exists(mutations_path) || !file.exists(auc_path)) {
    return(list(ok = FALSE, msg = "BeatAML2 data not found. Run setup_beataml2.R to download."))
  }

  # Mutations
  mut <- utils::read.delim(mutations_path, stringsAsFactors = FALSE, check.names = FALSE)
  mut$Sample <- as.character(mut$dbgap_sample_id)
  mut$VAF <- as.numeric(mut$t_vaf) * 100  # convert to 0-100
  mut$symbol <- as.character(mut$symbol)

  # FLT3 annotation (ITD vs TKD)
  vc <- tolower(as.character(mut$variant_classification))
  mut$Gene <- mut$symbol
  mut$Gene[mut$symbol == "FLT3" & vc == "inframe_insertion"] <- "FLT3-ITD"
  mut$Gene[mut$symbol == "FLT3" & vc %in% c("missense_variant", "inframe_deletion")] <- "FLT3-TKD"
  mut$Gene[mut$symbol == "FLT3" & mut$Gene == "FLT3"] <- "FLT3-TKD"  # default remaining

  # Keep one VAF per sample-gene (max if multiple)
  agg <- aggregate(VAF ~ Sample + Gene, data = mut, FUN = max, na.rm = TRUE)
  mut_wide <- agg

  # AUC data - use dbgap_dnaseq_sample for matching mutations
  auc <- utils::read.delim(auc_path, stringsAsFactors = FALSE, check.names = FALSE)
  auc$Sample <- ifelse(auc$dbgap_dnaseq_sample != "" & !is.na(auc$dbgap_dnaseq_sample),
    as.character(auc$dbgap_dnaseq_sample), as.character(auc$dbgap_rnaseq_sample))
  auc <- auc[auc$paper_inclusion == "TRUE" | auc$paper_inclusion == TRUE, , drop = FALSE]
  auc$auc <- as.numeric(auc$auc)
  auc <- auc[!is.na(auc$auc) & !is.na(auc$Sample) & auc$Sample != "", , drop = FALSE]

  # Clinical for subset
  clin <- NULL
  if (file.exists(clinical_path) && requireNamespace("readxl", quietly = TRUE)) {
    clin <- as.data.frame(readxl::read_excel(clinical_path), stringsAsFactors = FALSE)
    clin$Sample <- as.character(clin$dbgap_dnaseq_sample)
    na_idx <- clin$Sample == "" | is.na(clin$Sample)
    clin$Sample[na_idx] <- as.character(clin$dbgap_rnaseq_sample[na_idx])
  }

  # Merge: drug-gene pairs with both mutation and AUC
  samples_with_mut <- unique(mut_wide$Sample)
  samples_with_auc <- unique(auc$Sample)
  overlap_samples <- intersect(samples_with_mut, samples_with_auc)

  if (length(overlap_samples) < 5) {
    return(list(ok = FALSE, msg = paste("Only", length(overlap_samples), "samples with both mutation and drug data.")))
  }

  list(
    ok = TRUE,
    mutations = mut_wide,
    auc = auc,
    clinical = clin,
    overlap_samples = overlap_samples,
    data_dir = data_dir
  )
}

# Return sample IDs allowed for a given subset (same logic as correlation/LOOCV).
# Used so scatter plot and correlation table use identical data.
get_beataml_allowed_samples <- function(beataml, subset = "de_novo") {
  if (is.null(beataml) || !beataml$ok) return(character(0))
  if (!is.null(beataml$clinical) && subset != "All") {
    clin <- beataml$clinical
    clin <- clin[!is.na(clin$Sample) & clin$Sample != "", , drop = FALSE]
    if (nrow(clin) == 0) return(beataml$overlap_samples)
    if (subset == "de_novo") {
      keep <- as.character(clin$isDenovo) == "TRUE" | clin$isDenovo == TRUE
    } else if (subset == "secondary") {
      keep <- as.character(clin$isTransformed) == "TRUE" | clin$isTransformed == TRUE
    } else {
      keep <- rep(TRUE, nrow(clin))
    }
    allowed <- unique(clin$Sample[keep & !is.na(keep)])
    return(intersect(allowed, beataml$overlap_samples))
  }
  beataml$overlap_samples
}

# Compute drug-gene VAF-AUC correlations with LOOCV for overfitting (following Benard et al. methods)
# LOOCV: bestglm::LOOCV when available; otherwise manual LOO MSE. Same model auc ~ VAF.
compute_drug_vaf_correlations <- function(beataml, subset = "de_novo", min_n = 5, min_delta_vaf = 25, min_delta_auc = 75) {
  mut <- beataml$mutations
  auc <- beataml$auc

  if (!is.null(beataml$clinical) && subset != "All") {
    clin <- beataml$clinical
    clin <- clin[!is.na(clin$Sample) & clin$Sample != "", , drop = FALSE]
    if (subset == "de_novo") {
      keep <- as.character(clin$isDenovo) == "TRUE" | clin$isDenovo == TRUE
    } else if (subset == "secondary") {
      keep <- as.character(clin$isTransformed) == "TRUE" | clin$isTransformed == TRUE
    } else {
      keep <- rep(TRUE, nrow(clin))
    }
    allowed_samples <- unique(clin$Sample[keep & !is.na(keep)])
    allowed_samples <- intersect(allowed_samples, beataml$overlap_samples)
  } else {
    allowed_samples <- beataml$overlap_samples
  }

  mut <- mut[mut$Sample %in% allowed_samples, , drop = FALSE]
  auc <- auc[auc$Sample %in% allowed_samples, , drop = FALSE]

  drug_n <- sort(table(auc$inhibitor), decreasing = TRUE)
  drugs <- names(drug_n)[drug_n >= min_n]
  drugs <- head(drugs, 80)
  gene_n <- sort(table(mut$Gene), decreasing = TRUE)
  genes <- names(gene_n)[gene_n >= min_n]
  genes <- head(genes, 60)
  if (length(drugs) == 0 || length(genes) == 0) return(data.frame())

  use_bestglm <- requireNamespace("bestglm", quietly = TRUE)

  res <- list()
  for (d in drugs) {
    auc_d <- auc[auc$inhibitor == d, c("Sample", "auc"), drop = FALSE]
    if (nrow(auc_d) < min_n) next
    merged_d <- merge(auc_d, mut[, c("Sample", "Gene", "VAF")], by = "Sample")
    for (g in genes) {
      sub <- merged_d[merged_d$Gene == g, , drop = FALSE]
      if (nrow(sub) < min_n) next
      dv <- diff(range(sub$VAF, na.rm = TRUE))
      da <- diff(range(sub$auc, na.rm = TRUE))
      if (dv < min_delta_vaf || da < min_delta_auc) next
      fit <- tryCatch(lm(auc ~ VAF, data = sub), error = function(e) NULL)
      if (is.null(fit)) next
      cf <- summary(fit)$coefficients
      if (nrow(cf) < 2) next
      slope <- cf[2, 1]
      pval <- cf[2, 4]
      r2 <- summary(fit)$r.squared
      direction <- if (slope > 0) "Resistant" else "Sensitive"
      delta_auc <- slope * dv
      n_obs <- nrow(sub)

      if (use_bestglm) {
        X <- cbind(1, sub$VAF)
        y <- sub$auc
        loocv_out <- tryCatch(bestglm::LOOCV(X, y), error = function(e) NULL)
        if (!is.null(loocv_out) && length(loocv_out) >= 1) {
          loocv_mse <- as.numeric(loocv_out[1])
          loocv_mse_sd <- if (length(loocv_out) >= 2) as.numeric(loocv_out[2]) else NA
        } else {
          loocv_mse <- NA
          loocv_mse_sd <- NA
        }
      } else {
        loocv_mse <- NA
        loocv_mse_sd <- NA
      }
      if (!use_bestglm || is.na(loocv_mse)) {
        sq_errors <- vapply(seq_len(n_obs), function(i) {
          fit_loo <- tryCatch(lm(auc ~ VAF, data = sub[-i, ]), error = function(e) NULL)
          if (!is.null(fit_loo) && length(coef(fit_loo)) >= 2) {
            pred_i <- predict(fit_loo, newdata = sub[i, , drop = FALSE])
            (sub$auc[i] - pred_i)^2
          } else NA
        }, numeric(1))
        loocv_mse <- mean(sq_errors, na.rm = TRUE)
        loocv_mse_sd <- NA
      }

      # RMSE = sqrt(MSE) so value is in same units as AUC (0-100 scale)
      loocv_rmse <- if (!is.na(loocv_mse) && loocv_mse >= 0) sqrt(loocv_mse) else NA
      res[[length(res) + 1]] <- data.frame(
        Inhibitor = d, Gene = g, Slope = slope, R_squared = r2, p_value = pval,
        Direction = direction, n = n_obs, VAF_range = dv, AUC_range = da,
        delta_AUC = delta_auc, LOOCV_MSE = loocv_mse, LOOCV_MSE_sd = loocv_mse_sd,
        LOOCV_RMSE = loocv_rmse,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(res) == 0) return(data.frame())
  out <- do.call(rbind, res)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out
}

# Leave-one-out cross-validation for VAF-AUC correlations (MSE via bestglm::LOOCV)
# See https://www.rdocumentation.org/packages/bestglm/versions/0.37.3/topics/LOOCV
compute_drug_vaf_loo <- function(beataml, subset = "de_novo", min_n = 8, min_delta_vaf = 25, min_delta_auc = 75) {
  mut <- beataml$mutations
  auc <- beataml$auc

  if (!is.null(beataml$clinical) && subset != "All") {
    clin <- beataml$clinical
    clin <- clin[!is.na(clin$Sample) & clin$Sample != "", , drop = FALSE]
    if (subset == "de_novo") {
      keep <- as.character(clin$isDenovo) == "TRUE" | clin$isDenovo == TRUE
    } else if (subset == "secondary") {
      keep <- as.character(clin$isTransformed) == "TRUE" | clin$isTransformed == TRUE
    } else {
      keep <- rep(TRUE, nrow(clin))
    }
    allowed_samples <- unique(clin$Sample[keep & !is.na(keep)])
    allowed_samples <- intersect(allowed_samples, beataml$overlap_samples)
  } else {
    allowed_samples <- beataml$overlap_samples
  }

  mut <- mut[mut$Sample %in% allowed_samples, , drop = FALSE]
  auc <- auc[auc$Sample %in% allowed_samples, , drop = FALSE]

  drug_n <- sort(table(auc$inhibitor), decreasing = TRUE)
  drugs <- head(names(drug_n)[drug_n >= min_n], 40)
  gene_n <- sort(table(mut$Gene), decreasing = TRUE)
  genes <- head(names(gene_n)[gene_n >= min_n], 30)
  if (length(drugs) == 0 || length(genes) == 0) return(data.frame())

  use_bestglm <- requireNamespace("bestglm", quietly = TRUE)

  res <- list()
  for (d in drugs) {
    auc_d <- auc[auc$inhibitor == d, c("Sample", "auc"), drop = FALSE]
    if (nrow(auc_d) < min_n) next
    merged_d <- merge(auc_d, mut[, c("Sample", "Gene", "VAF")], by = "Sample")
    for (g in genes) {
      sub <- merged_d[merged_d$Gene == g, , drop = FALSE]
      if (nrow(sub) < min_n) next
      dv <- diff(range(sub$VAF, na.rm = TRUE))
      da <- diff(range(sub$auc, na.rm = TRUE))
      if (dv < min_delta_vaf || da < min_delta_auc) next

      fit_full <- tryCatch(lm(auc ~ VAF, data = sub), error = function(e) NULL)
      if (is.null(fit_full)) next
      slope_full <- coef(fit_full)[2]
      n_obs <- nrow(sub)

      if (use_bestglm) {
        X <- cbind(1, sub$VAF)
        y <- sub$auc
        loocv_out <- tryCatch(bestglm::LOOCV(X, y), error = function(e) NULL)
        if (!is.null(loocv_out) && length(loocv_out) >= 1) {
          mse <- as.numeric(loocv_out[1])
          mse_sd <- if (length(loocv_out) >= 2) as.numeric(loocv_out[2]) else NA
        } else {
          mse <- NA
          mse_sd <- NA
        }
      } else {
        mse <- NA
        mse_sd <- NA
      }

      slopes_loo <- numeric(n_obs)
      for (i in seq_len(n_obs)) {
        fit_loo <- tryCatch(lm(auc ~ VAF, data = sub[-i, ]), error = function(e) NULL)
        if (!is.null(fit_loo) && length(coef(fit_loo)) >= 2) {
          slopes_loo[i] <- coef(fit_loo)[2]
        } else {
          slopes_loo[i] <- NA
        }
      }
      slopes_loo_valid <- slopes_loo[!is.na(slopes_loo)]
      if (length(slopes_loo_valid) < n_obs * 0.8) next

      sign_flip <- sum((slopes_loo_valid > 0) != (slope_full > 0), na.rm = TRUE)
      pct_stable <- 100 * (1 - sign_flip / length(slopes_loo_valid))
      if (!use_bestglm || is.na(mse)) {
        sq_errors <- numeric(n_obs)
        for (i in seq_len(n_obs)) {
          fit_loo <- tryCatch(lm(auc ~ VAF, data = sub[-i, ]), error = function(e) NULL)
          if (!is.null(fit_loo) && length(coef(fit_loo)) >= 2) {
            pred_i <- predict(fit_loo, newdata = sub[i, , drop = FALSE])
            sq_errors[i] <- (sub$auc[i] - pred_i)^2
          } else {
            sq_errors[i] <- NA
          }
        }
        mse <- mean(sq_errors, na.rm = TRUE)
        mse_sd <- NA
      }
      rmse <- if (!is.na(mse) && mse >= 0) sqrt(mse) else NA
      res[[length(res) + 1]] <- data.frame(
        Inhibitor = d, Gene = g, Slope_full = slope_full,
        Slope_loo_mean = mean(slopes_loo_valid, na.rm = TRUE),
        Slope_loo_sd = sd(slopes_loo_valid, na.rm = TRUE),
        MSE = mse, MSE_sd = mse_sd, RMSE = rmse, n_sign_flip = sign_flip, pct_stable = pct_stable,
        n = n_obs, stringsAsFactors = FALSE
      )
    }
  }
  if (length(res) == 0) return(data.frame())
  do.call(rbind, res)
}
