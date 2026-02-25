# =============================================================================
# Precompute Drug Sensitivity results for Drug Sensitivity tab.
# Run once (e.g. before deploy) to generate drug_sensitivity_precomputed.rds.
# The Shiny app will load this file and use it instead of computing on the fly,
# reducing memory use and avoiding OOM on shinyapps.io free tier.
# Computes both VAF vs AUC and CCF vs AUC results for the same subsets.
# =============================================================================

# Run from the Meta_AML_Shiny app directory:
#   setwd("/path/to/Meta_AML_Shiny"); source("precompute_drug_sensitivity.R")
# Or: cd /path/to/Meta_AML_Shiny && Rscript precompute_drug_sensitivity.R
app_dir <- if (exists("ofile")) dirname(ofile) else getwd()
if (nzchar(app_dir)) setwd(app_dir)
if (!file.exists("load_beataml2.R")) stop("Run this script from the Meta_AML_Shiny app directory (where load_beataml2.R lives).")
source("load_beataml2.R")

b <- load_beataml2()
if (!b$ok) stop("BeatAML2 data not loaded: ", b$msg)

# Enrich Beat AML mutations with CCF from the meta cohort (for CCF vs AUC precomputation)
b_with_ccf <- b
cohort_path <- file.path(getwd(), "AML_Meta_Cohort_v2.rds")
if (!file.exists(cohort_path)) cohort_path <- file.path(getwd(), "AML_Meta_Cohort.rds")
if (file.exists(cohort_path)) {
  cohort <- tryCatch(as.data.frame(readRDS(cohort_path), stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(cohort) && nrow(cohort) > 0 &&
      "Cohort" %in% colnames(cohort) && "CCF" %in% colnames(cohort) &&
      "Sample" %in% colnames(cohort) && "Gene" %in% colnames(cohort)) {
    beat <- cohort[as.character(cohort$Cohort) %in% c("Beat AML", "Beat_AML"),
                   c("Sample", "Gene", "CCF"), drop = FALSE]
    if (nrow(beat) > 0) {
      beat$CCF <- as.numeric(beat$CCF)
      beat <- aggregate(CCF ~ Sample + Gene, data = beat, FUN = max, na.rm = TRUE)
      b_with_ccf$mutations <- merge(b$mutations, beat, by = c("Sample", "Gene"), all.x = TRUE)
      message("Merged CCF from cohort (Beat AML) for ", nrow(beat), " sample-gene pairs.")
    } else {
      message("No Beat AML rows in cohort; CCF precomputation will be skipped.")
      b_with_ccf <- NULL
    }
  } else {
    message("Cohort missing required columns (Cohort, CCF, Sample, Gene); CCF precomputation skipped.")
    b_with_ccf <- NULL
  }
} else {
  message("Cohort RDS not found; CCF precomputation skipped.")
  b_with_ccf <- NULL
}
if (is.null(b_with_ccf) || !"CCF" %in% colnames(b_with_ccf$mutations)) b_with_ccf <- NULL

subsets <- c("All", "de_novo", "secondary")
precomputed <- list(meta_aml4 = list())

message("Precomputing for Meta AML4 (all waves)...")
for (sub in subsets) {
  message("  subset: ", sub)
  precomputed$meta_aml4[[sub]] <- list(
    mut_wt = compute_mut_wt_all(b, sub),
    correlations = compute_drug_vaf_correlations(b, subset = sub, metric = "VAF"),
    loo = compute_drug_vaf_loo(b, subset = sub, metric = "VAF")
  )
  if (!is.null(b_with_ccf)) {
    precomputed$meta_aml4[[sub]]$correlations_ccf <- compute_drug_vaf_correlations(b_with_ccf, subset = sub, metric = "CCF")
    precomputed$meta_aml4[[sub]]$loo_ccf <- compute_drug_vaf_loo(b_with_ccf, subset = sub, metric = "CCF")
    message("    CCF vs AUC done.")
  } else {
    precomputed$meta_aml4[[sub]]$correlations_ccf <- NULL
    precomputed$meta_aml4[[sub]]$loo_ccf <- NULL
  }
}

out_path <- "drug_sensitivity_precomputed.rds"
saveRDS(precomputed, out_path)
message("Saved to ", out_path)
