# =============================================================================
# Build a combined clinical annotation file for all samples in the four datasets.
# Run from Meta_AML_Shiny: Rscript build_aggregated_clinical.R
# Output: Meta_AML_aggregated_clinical.tsv
# =============================================================================

setwd(if (exists("ofile")) dirname(ofile) else getwd())
if (!dir.exists("laml_tcga") || !dir.exists("beataml2_data") || !dir.exists("amlsg_data") || !dir.exists("UK_NCRI_data")) {
  stop("Run from Meta_AML_Shiny root.")
}

# Common output columns (one row per sample)
CLIN_COLS <- c(
  "Patient_ID", "Cohort", "Age", "Sex", "AML_type",
  "ELN_2017_Risk", "ELN_2022_Risk", "Time_to_OS_years", "OS_status", "WBC", "Hemoglobin", "Platelet",
  "BM_blasts", "PB_blasts", "Karyotype", "source_file"
)

fill_clin <- function(n) {
  x <- as.data.frame(matrix(NA_character_, nrow = n, ncol = length(CLIN_COLS)))
  colnames(x) <- CLIN_COLS
  x
}

# ---------- 1. TCGA LAML ----------
tcga_clin <- NULL
pat_path <- "laml_tcga/data_clinical_patient.txt"
samp_path <- "laml_tcga/data_clinical_sample.txt"
if (file.exists(pat_path) && file.exists(samp_path)) {
  pat <- read.delim(pat_path, skip = 4, stringsAsFactors = FALSE, check.names = FALSE)
  samp <- read.delim(samp_path, skip = 4, stringsAsFactors = FALSE, check.names = FALSE)
  pat_id <- if ("PATIENT_ID" %in% colnames(pat)) "PATIENT_ID" else "Patient Identifier"
  samp_id <- if ("SAMPLE_ID" %in% colnames(samp)) "SAMPLE_ID" else "Sample Identifier"
  pid_samp <- if ("PATIENT_ID" %in% colnames(samp)) "PATIENT_ID" else "Patient Identifier"
  m <- merge(samp[, c(pid_samp, samp_id)], pat, by.x = pid_samp, by.y = pat_id, all.x = TRUE)
  tcga_clin <- fill_clin(nrow(m))
  tcga_clin$Patient_ID <- as.character(m[[samp_id]])
  tcga_clin$Cohort <- "TCGA"
  tcga_clin$Age <- if ("AGE" %in% colnames(m)) as.character(m$AGE) else if ("Diagnosis Age" %in% colnames(m)) as.character(m$`Diagnosis Age`) else NA_character_
  tcga_clin$Sex <- if ("SEX" %in% colnames(m)) as.character(m$SEX) else if ("Sex" %in% colnames(m)) as.character(m$Sex) else NA_character_
  tcga_clin$AML_type <- "De novo"
  if ("OS_MONTHS" %in% colnames(m)) {
    os_mo <- as.numeric(m$OS_MONTHS)
    tcga_clin$Time_to_OS_years <- as.character(round(os_mo / 12, 2))
  } else if ("Overall Survival (Months)" %in% colnames(m)) {
    os_mo <- as.numeric(m$`Overall Survival (Months)`)
    tcga_clin$Time_to_OS_years <- as.character(round(os_mo / 12, 2))
  }
  if ("OS_STATUS" %in% colnames(m)) {
    os <- as.character(m$OS_STATUS)
    tcga_clin$OS_status <- ifelse(grepl("DECEASED|1:", os, ignore.case = TRUE), "1", "0")
  } else if ("Overall Survival Status" %in% colnames(m)) {
    os <- as.character(m$`Overall Survival Status`)
    tcga_clin$OS_status <- ifelse(grepl("DECEASED|1:", os, ignore.case = TRUE), "1", "0")
  }
  if ("PLATELET_COUNT_PRERESECTION" %in% colnames(m)) tcga_clin$Platelet <- as.character(m$PLATELET_COUNT_PRERESECTION)
  else if ("Platelet count preresection" %in% colnames(m)) tcga_clin$Platelet <- as.character(m$`Platelet count preresection`)
  if ("BLAST_COUNT" %in% colnames(m)) tcga_clin$BM_blasts <- as.character(m$BLAST_COUNT)
  else if ("Blast Count" %in% colnames(m)) tcga_clin$BM_blasts <- as.character(m$`Blast Count`)
  # Simplify TCGA cytogenetics to Normal / Complex / Other/Unknown
  raw_karyo <- NA_character_
  if ("CYTOGENETIC_ABNORMALITY_TYPE" %in% colnames(m)) raw_karyo <- as.character(m$CYTOGENETIC_ABNORMALITY_TYPE)
  else if ("Cytogenetic abnormality type" %in% colnames(m)) raw_karyo <- as.character(m$`Cytogenetic abnormality type`)
  raw_karyo[is.na(raw_karyo) | raw_karyo == "" | raw_karyo == "[Not Available]"] <- NA_character_
  simp_karyo <- rep("Other/Unknown", length(raw_karyo))
  is_complex <- !is.na(raw_karyo) & (grepl("complex", raw_karyo, ignore.case = TRUE) | grepl(">=\\s*3", raw_karyo))
  simp_karyo[is_complex] <- "Complex"
  # Normal only if the string is exactly 'Normal' (case/whitespace-insensitive)
  is_normal <- !is.na(raw_karyo) & (tolower(trimws(raw_karyo)) == "normal")
  simp_karyo[!is_complex & is_normal] <- "Normal"
  tcga_clin$Karyotype <- simp_karyo
  # Derive ELN 2017 risk for TCGA from cytogenetics + mutations (no ELN in source clinical)
  mut_path <- "laml_tcga/data_mutations.txt"
  if (file.exists(mut_path)) {
    mut <- read.delim(mut_path, stringsAsFactors = FALSE, check.names = FALSE)
    samp_col <- if ("Tumor_Sample_Barcode" %in% colnames(mut)) "Tumor_Sample_Barcode" else grep("Sample|Tumor", colnames(mut), value = TRUE)[1]
    gene_col <- if ("Hugo_Symbol" %in% colnames(mut)) "Hugo_Symbol" else "Gene"
    vc_col <- "Variant_Classification"
    if (!vc_col %in% colnames(mut)) vc_col <- "Variant_Classification"
    samples <- tcga_clin$Patient_ID
    mut$Sample <- as.character(mut[[samp_col]])
    mut$Gene <- as.character(mut[[gene_col]])
    mut <- mut[mut$Sample %in% samples, , drop = FALSE]
    vc <- tolower(as.character(mut[[vc_col]]))
    mut$is_ITD <- vc %in% c("in_frame_ins", "in_frame_insertion")
    has_gene <- function(g) samples %in% unique(mut$Sample[mut$Gene == g])
    n_cebpa <- setNames(rep(0L, length(samples)), samples)
    if (any(mut$Gene == "CEBPA")) {
      n_cebpa[names(which(table(mut$Sample[mut$Gene == "CEBPA"]) >= 2))] <- 2L
      n_cebpa[names(which(table(mut$Sample[mut$Gene == "CEBPA"]) == 1))] <- 1L
    }
    NPM1_mut <- has_gene("NPM1")
    FLT3_ITD <- samples %in% mut$Sample[mut$Gene == "FLT3" & mut$is_ITD]
    FLT3_TKD <- samples %in% mut$Sample[mut$Gene == "FLT3" & !mut$is_ITD]
    CEBPA_bi <- n_cebpa[samples] >= 2
    CEBPA_bi[is.na(CEBPA_bi)] <- FALSE
    RUNX1_mut <- has_gene("RUNX1")
    ASXL1_mut <- has_gene("ASXL1")
    TP53_mut <- has_gene("TP53")
    cyto_fav <- !is.na(raw_karyo) & (grepl("t\\(8;21\\)|t\\(8; 21\\)|inv\\(16\\)|t\\(16;16\\)|t\\(15;17\\)", raw_karyo, ignore.case = TRUE))
    cyto_adv <- !is.na(raw_karyo) & (is_complex | grepl("inv\\(3\\)|t\\(3;3\\)|t\\(6;9\\)|del\\(5q\\)|-5|del\\(7q\\)|-7|17p|t\\(v;11\\)|t\\(9;22\\)", raw_karyo, ignore.case = TRUE))
    adverse <- cyto_adv | RUNX1_mut | ASXL1_mut | TP53_mut | (!NPM1_mut & FLT3_ITD)
    favorable <- (cyto_fav | (NPM1_mut & !FLT3_ITD & !adverse) | (CEBPA_bi & !adverse)) & !adverse
    tcga_clin$ELN_2017_Risk <- "Intermediate"
    tcga_clin$ELN_2017_Risk[adverse] <- "Adverse"
    tcga_clin$ELN_2017_Risk[favorable] <- "Favorable"
    # ELN 2022: FLT3-ITD all intermediate; MDS-related genes adverse; NPM1+adverse cyto = adverse
    MDS_genes <- c("BCOR", "EZH2", "SF3B1", "SRSF2", "STAG2", "U2AF1", "ZRSR2")
    MDS_adv <- Reduce(`|`, lapply(MDS_genes, function(g) has_gene(g)))
    CEBPA_any <- n_cebpa[samples] >= 1
    CEBPA_any[is.na(CEBPA_any)] <- FALSE
    adverse_2022 <- cyto_adv | RUNX1_mut | ASXL1_mut | TP53_mut | MDS_adv | (NPM1_mut & cyto_adv)
    favorable_2022 <- (cyto_fav | CEBPA_any | (NPM1_mut & !FLT3_ITD & !adverse_2022)) & !adverse_2022
    tcga_clin$ELN_2022_Risk <- "Intermediate"
    tcga_clin$ELN_2022_Risk[adverse_2022] <- "Adverse"
    tcga_clin$ELN_2022_Risk[favorable_2022] <- "Favorable"
  }
  tcga_clin$source_file <- "laml_tcga/data_clinical_patient.txt + data_clinical_sample.txt"
}

# ---------- 2. Beat AML2 (patient-level: one row per patient, ~800 pts) ----------
beat_clin <- NULL
beat_path <- "beataml2_data/beataml_wv1to4_clinical.xlsx"
if (file.exists(beat_path) && requireNamespace("readxl", quietly = TRUE)) {
  b <- as.data.frame(readxl::read_excel(beat_path), stringsAsFactors = FALSE)
  pid_col <- "dbgap_subject_id"
  if (!pid_col %in% colnames(b)) pid_col <- grep("subject|patient", colnames(b), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(pid_col) || length(pid_col) == 0) pid_col <- colnames(b)[1]
  pids <- as.character(b[[pid_col]])
  pids <- pids[!is.na(pids) & pids != ""]
  if (length(pids) > 0) {
    uids <- unique(pids)
    beat_clin <- fill_clin(length(uids))
    beat_clin$Patient_ID <- uids
    beat_clin$Cohort <- "Beat AML"
    b_sub <- b[match(uids, pids), , drop = FALSE]
    if ("ageAtDiagnosis" %in% colnames(b_sub)) beat_clin$Age <- as.character(round(suppressWarnings(as.numeric(b_sub$ageAtDiagnosis)), 1))
    if ("consensus_sex" %in% colnames(b_sub)) beat_clin$Sex <- as.character(b_sub$consensus_sex)
    beat_clin$AML_type <- "Other"
    if ("isDenovo" %in% colnames(b_sub)) beat_clin$AML_type[as.character(b_sub$isDenovo) %in% c("TRUE", "True", "1")] <- "De novo"
    if ("isTransformed" %in% colnames(b_sub)) beat_clin$AML_type[as.character(b_sub$isTransformed) %in% c("TRUE", "True", "1")] <- "Secondary"
    if ("ELN2017" %in% colnames(b_sub)) beat_clin$ELN_2017_Risk <- as.character(b_sub$ELN2017)
    if ("overallSurvival" %in% colnames(b_sub)) {
      os_days <- suppressWarnings(as.numeric(b_sub$overallSurvival))
      beat_clin$Time_to_OS_years <- as.character(round(os_days / 365.25, 2))
    }
    if ("vitalStatus" %in% colnames(b_sub)) beat_clin$OS_status <- ifelse(grepl("dead|deceased", as.character(b_sub$vitalStatus), ignore.case = TRUE), "1", "0")
    if ("wbcCount" %in% colnames(b_sub)) beat_clin$WBC <- as.character(b_sub$wbcCount)
    if ("hemoglobin" %in% colnames(b_sub)) beat_clin$Hemoglobin <- as.character(b_sub$hemoglobin)
    if ("plateletCount" %in% colnames(b_sub)) beat_clin$Platelet <- as.character(b_sub$plateletCount)
    if ("%.Blasts.in.BM" %in% colnames(b_sub)) beat_clin$BM_blasts <- as.character(b_sub$`%.Blasts.in.BM`)
    if ("%.Blasts.in.PB" %in% colnames(b_sub)) beat_clin$PB_blasts <- as.character(b_sub$`%.Blasts.in.PB`)
    # Simplify Beat AML cytogenetics to Normal / Complex / Other/Unknown
    if ("karyotype" %in% colnames(b_sub)) {
      raw_karyo <- as.character(b_sub$karyotype)
      raw_karyo[is.na(raw_karyo) | raw_karyo == ""] <- NA_character_
      simp_karyo <- rep("Other/Unknown", length(raw_karyo))
      is_normal <- !is.na(raw_karyo) & grepl("^46,(XX|XY)(\\[[0-9]+\\])?$", trimws(raw_karyo), ignore.case = TRUE)
      simp_karyo[is_normal] <- "Normal"

      # Approximate complex as >=3 abnormality markers in the ISCN string
      abnormality_count <- rep(NA_integer_, length(raw_karyo))
      has_val <- !is.na(raw_karyo)
      if (any(has_val)) {
        s <- tolower(raw_karyo[has_val])
        markers <- gregexpr("t\\(|inv\\(|del\\(|dup\\(|add\\(|der\\(|dic\\(|ins\\(|i\\(|r\\(|\\+[0-9xy]+|-([0-9xy]+)|mar|idem|trp", s, perl = TRUE)
        abnormality_count[has_val] <- vapply(markers, function(m) if (length(m) == 1 && m[1] == -1) 0L else length(m), integer(1))
      }
      is_complex <- !is_normal & !is.na(abnormality_count) & abnormality_count >= 3

      # Also treat explicit "complex" text as complex if available
      if ("otherCytogenetics" %in% colnames(b_sub)) {
        oc <- as.character(b_sub$otherCytogenetics)
        oc[is.na(oc) | oc == ""] <- NA_character_
        is_complex <- is_complex | (!is.na(oc) & grepl("complex", oc, ignore.case = TRUE))
      }
      simp_karyo[is_complex] <- "Complex"
      beat_clin$Karyotype <- simp_karyo
    }
    beat_clin$source_file <- "beataml2_data/beataml_wv1to4_clinical.xlsx"
  }
}

# ---------- 3. AML-SG ----------
amlsg_clin <- NULL
if (file.exists("amlsg_data/AMLSG_Clinical_Anon.RData")) {
  load("amlsg_data/AMLSG_Clinical_Anon.RData")
  cl <- clinicalData
  amlsg_clin <- fill_clin(nrow(cl))
  amlsg_clin$Patient_ID <- as.character(cl$PDID)
  amlsg_clin$Cohort <- "AML-SG"
  if ("AOD" %in% colnames(cl)) amlsg_clin$Age <- as.character(cl$AOD)
  if ("gender" %in% colnames(cl)) {
    g <- as.character(cl$gender)
    amlsg_clin$Sex <- ifelse(g == "1", "Male", ifelse(g == "2", "Female", NA_character_))
  }
  if ("HB" %in% colnames(cl)) amlsg_clin$Hemoglobin <- as.character(cl$HB)
  if ("platelet" %in% colnames(cl)) amlsg_clin$Platelet <- as.character(cl$platelet)
  if ("wbc" %in% colnames(cl)) amlsg_clin$WBC <- as.character(cl$wbc)
  if ("BM_Blasts" %in% colnames(cl)) amlsg_clin$BM_blasts <- as.character(cl$BM_Blasts)
  if ("PB_Blasts" %in% colnames(cl)) amlsg_clin$PB_blasts <- as.character(cl$PB_Blasts)
  if ("Status" %in% colnames(cl)) amlsg_clin$OS_status <- as.character(cl$Status)
  if ("OS" %in% colnames(cl)) amlsg_clin$Time_to_OS_years <- as.character(round(as.numeric(cl$OS) / 365.25, 2))
  if ("TypeAML" %in% colnames(cl)) amlsg_clin$AML_type <- as.character(cl$TypeAML)
  # ELN risk (Papaemmanuil NEJM cohort) from df_nejm_personalization.tsv
  nejm_path <- "UK_NCRI_data/data/df_nejm_personalization.tsv"
  if (file.exists(nejm_path)) {
    nejm <- tryCatch(
      read.table(nejm_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (!is.null(nejm) && nrow(nejm) > 0 && "PDID" %in% colnames(nejm) && "eln_2017" %in% colnames(nejm)) {
      midx <- match(amlsg_clin$Patient_ID, as.character(nejm$PDID))
      eln <- tolower(as.character(nejm$eln_2017[midx]))
      amlsg_clin$ELN_2017_Risk <- ifelse(eln == "favorable", "Favorable",
        ifelse(eln == "intermediate", "Intermediate",
          ifelse(eln == "adverse", "Adverse", NA_character_)
        )
      )
    }
  }
  if ("NK" %in% colnames(cl) && "complex" %in% colnames(cl)) {
    nk <- as.numeric(cl$NK)
    comp <- as.numeric(cl$complex)
    amlsg_clin$Karyotype <- "Other/Unknown"
    amlsg_clin$Karyotype[!is.na(nk) & nk == 1] <- "Normal"
    amlsg_clin$Karyotype[!is.na(comp) & comp == 1] <- "Complex"
  }
  amlsg_clin$source_file <- "amlsg_data/AMLSG_Clinical_Anon.RData"
}

# ---------- 4. UK-NCRI ----------
# Pull clinical and karyotype from both aml_prognosis_updated.tsv and UK_NCRI_Clinical_data.csv
ncri_clin <- NULL
ncri_path <- "UK_NCRI_data/UK_NCRI_Clinical_data.csv"
prog_path <- "UK_NCRI_data/data/aml_prognosis_updated.tsv"
prog <- NULL
if (file.exists(prog_path)) {
  prog <- read.delim(prog_path, stringsAsFactors = FALSE, check.names = FALSE)
  keep <- !is.na(prog[[1]]) & as.character(prog[[1]]) != ""
  if (ncol(prog) >= 2) keep <- keep & (as.character(prog[[2]]) != "eln_2017_adverse")
  prog <- prog[keep, , drop = FALSE]
  prog_ids <- rownames(prog)
  if (length(prog_ids) > 0 && all(grepl("^PD[0-9]", prog_ids))) {
    prog$Patient_ID <- prog_ids
  } else if ("eln_2017" %in% colnames(prog) && is.character(prog[[1]])) {
    prog$Patient_ID <- as.character(prog[[1]])
  } else {
    prog$Patient_ID <- rownames(prog)
  }
}
n <- NULL
if (file.exists(ncri_path)) {
  raw <- readLines(ncri_path, n = 3)
  skip <- if (any(grepl("^S\\.Data|^#", raw))) 1 else 0
  n <- read.csv(ncri_path, skip = skip, stringsAsFactors = FALSE, check.names = FALSE)
  keep <- !is.na(n[[1]]) & n[[1]] != ""
  if (ncol(n) >= 2) keep <- keep & n[[2]] != "ahd"
  n <- n[keep, , drop = FALSE]
  n[[1]] <- as.character(n[[1]])
}
if (!is.null(prog) && nrow(prog) > 0) {
  ncri_ids <- prog$Patient_ID
} else if (!is.null(n) && nrow(n) > 0) {
  ncri_ids <- n[[1]]
} else {
  ncri_ids <- character(0)
}
if (length(ncri_ids) > 0) {
  ncri_clin <- fill_clin(length(ncri_ids))
  ncri_clin$Patient_ID <- ncri_ids
  ncri_clin$Cohort <- "UK-NCRI"
  ncri_clin$source_file <- "UK_NCRI_data/UK_NCRI_Clinical_data.csv + data/aml_prognosis_updated.tsv"
  midx_prog <- if (!is.null(prog)) match(ncri_ids, prog$Patient_ID) else rep(NA_integer_, length(ncri_ids))
  midx_csv  <- if (!is.null(n)) match(ncri_ids, n[[1]]) else rep(NA_integer_, length(ncri_ids))
  if (!is.null(prog) && nrow(prog) > 0) {
    idx <- midx_prog
    if ("age" %in% colnames(prog)) ncri_clin$Age <- as.character(round(suppressWarnings(as.numeric(prog$age[idx])), 1))
    if ("gender" %in% colnames(prog)) {
      g <- as.character(prog$gender[idx])
      ncri_clin$Sex <- ifelse(g == "0", "Male", ifelse(g == "1", "Female", NA_character_))
    }
    if ("secondary" %in% colnames(prog)) {
      sec <- as.character(prog$secondary[idx])
      ncri_clin$AML_type <- ifelse(sec == "1", "De novo", ifelse(sec == "2", "Secondary", "Other"))
    }
    if ("wbc" %in% colnames(prog)) ncri_clin$WBC <- as.character(prog$wbc[idx])
    if ("hb" %in% colnames(prog)) ncri_clin$Hemoglobin <- as.character(prog$hb[idx])
    if ("plt" %in% colnames(prog)) ncri_clin$Platelet <- as.character(prog$plt[idx])
    if ("bm_blasts" %in% colnames(prog)) ncri_clin$BM_blasts <- as.character(prog$bm_blasts[idx])
    if ("eln_2017_favorable" %in% colnames(prog)) {
      fav <- as.numeric(prog$eln_2017_favorable[idx])
      int <- as.numeric(prog$eln_2017_intermediate[idx])
      adv <- as.numeric(prog$eln_2017_adverse[idx])
      ncri_clin$ELN_2017_Risk <- "Other/Unknown"
      ncri_clin$ELN_2017_Risk[fav == 1] <- "Favorable"
      ncri_clin$ELN_2017_Risk[int == 1] <- "Intermediate"
      ncri_clin$ELN_2017_Risk[adv == 1] <- "Adverse"
    }
    if ("os" %in% colnames(prog)) ncri_clin$Time_to_OS_years <- as.character(prog$os[idx])
    if ("os_status" %in% colnames(prog)) ncri_clin$OS_status <- as.character(prog$os_status[idx])
    # Karyotype assignment from cytogenetic flags in aml_prognosis_updated.tsv:
    # - Complex if complex==1
    # - Normal if no cytogenetic abnormality flags are present
    # - Other/Unknown otherwise
    cyto_cols <- grep("^(add_|del_|t_|inv_|minus)|^(others_transloc|complex)$", colnames(prog), value = TRUE)
    if (length(cyto_cols) > 0) {
      cyto_mat <- suppressWarnings(sapply(prog[idx, cyto_cols, drop = FALSE], function(x) as.numeric(as.character(x))))
      cyto_mat <- as.matrix(cyto_mat)
      abnormal <- rowSums(cyto_mat == 1, na.rm = TRUE) > 0
      ncri_clin$Karyotype <- "Other/Unknown"
      ncri_clin$Karyotype[!abnormal] <- "Normal"
      if ("complex" %in% colnames(prog)) {
        comp <- suppressWarnings(as.numeric(as.character(prog$complex[idx])))
        ncri_clin$Karyotype[!is.na(comp) & comp == 1] <- "Complex"
      }
    } else if ("complex" %in% colnames(prog)) {
      comp <- suppressWarnings(as.numeric(as.character(prog$complex[idx])))
      ncri_clin$Karyotype <- "Other/Unknown"
      ncri_clin$Karyotype[!is.na(comp) & comp == 1] <- "Complex"
    }
  }
  if (!is.null(n) && nrow(n) > 0) {
    idx <- midx_csv
    nas <- is.na(idx)
    if (any(!nas)) {
      fill_age    <- is.na(ncri_clin$Age) | ncri_clin$Age == ""
      if (ncol(n) >= 10 && any(fill_age & !nas)) ncri_clin$Age[fill_age & !nas] <- as.character(round(suppressWarnings(as.numeric(n[[10]][idx[fill_age & !nas]])), 1))
      fill_sex    <- is.na(ncri_clin$Sex) | ncri_clin$Sex == ""
      if (ncol(n) >= 9 && any(fill_sex & !nas)) { g <- as.character(n[[9]][idx[fill_sex & !nas]]); ncri_clin$Sex[fill_sex & !nas] <- ifelse(g == "0", "Male", ifelse(g == "1", "Female", NA_character_)) }
      fill_aml    <- is.na(ncri_clin$AML_type) | ncri_clin$AML_type == ""
      if (ncol(n) >= 5 && any(fill_aml & !nas)) { sec <- as.character(n[[5]][idx[fill_aml & !nas]]); ncri_clin$AML_type[fill_aml & !nas] <- ifelse(sec == "1", "De novo", ifelse(sec == "2", "Secondary", "Other")) }
      fill_wbc    <- is.na(ncri_clin$WBC) | ncri_clin$WBC == ""
      if (ncol(n) >= 6 && any(fill_wbc & !nas)) ncri_clin$WBC[fill_wbc & !nas] <- as.character(n[[6]][idx[fill_wbc & !nas]])
      fill_hb     <- is.na(ncri_clin$Hemoglobin) | ncri_clin$Hemoglobin == ""
      if (ncol(n) >= 7 && any(fill_hb & !nas)) ncri_clin$Hemoglobin[fill_hb & !nas] <- as.character(n[[7]][idx[fill_hb & !nas]])
      fill_plt    <- is.na(ncri_clin$Platelet) | ncri_clin$Platelet == ""
      if (ncol(n) >= 8 && any(fill_plt & !nas)) ncri_clin$Platelet[fill_plt & !nas] <- as.character(n[[8]][idx[fill_plt & !nas]])
      fill_bm     <- is.na(ncri_clin$BM_blasts) | ncri_clin$BM_blasts == ""
      if (ncol(n) >= 4 && any(fill_bm & !nas)) ncri_clin$BM_blasts[fill_bm & !nas] <- as.character(n[[4]][idx[fill_bm & !nas]])
    }
  }
  if (is.null(prog) || nrow(prog) == 0) {
    if (!is.null(n) && nrow(n) > 0 && nrow(ncri_clin) == nrow(n)) {
      ncri_clin$Patient_ID <- as.character(n[[1]])
      if (ncol(n) >= 10) ncri_clin$Age <- as.character(round(suppressWarnings(as.numeric(n[[10]])), 1))
      if (ncol(n) >= 9) { g <- as.character(n[[9]]); ncri_clin$Sex <- ifelse(g == "0", "Male", ifelse(g == "1", "Female", NA_character_)) }
      if (ncol(n) >= 5) { sec <- as.character(n[[5]]); ncri_clin$AML_type <- ifelse(sec == "1", "De novo", ifelse(sec == "2", "Secondary", "Other")) }
      if (ncol(n) >= 6) ncri_clin$WBC <- as.character(n[[6]])
      if (ncol(n) >= 7) ncri_clin$Hemoglobin <- as.character(n[[7]])
      if (ncol(n) >= 8) ncri_clin$Platelet <- as.character(n[[8]])
      if (ncol(n) >= 4) ncri_clin$BM_blasts <- as.character(n[[4]])
      ncri_clin$source_file <- "UK_NCRI_data/UK_NCRI_Clinical_data.csv"
    }
  }
}

# ---------- Bind ----------
out_list <- list()
if (!is.null(tcga_clin) && nrow(tcga_clin) > 0) out_list$tcga <- tcga_clin
if (!is.null(beat_clin) && nrow(beat_clin) > 0) out_list$beat <- beat_clin
if (!is.null(amlsg_clin) && nrow(amlsg_clin) > 0) out_list$amlsg <- amlsg_clin
if (!is.null(ncri_clin) && nrow(ncri_clin) > 0) out_list$ncri <- ncri_clin

if (length(out_list) == 0) stop("No clinical files found.")

all_cols <- unique(c(CLIN_COLS, unlist(lapply(out_list, colnames))))
for (i in seq_along(out_list)) {
  for (c in all_cols) if (!c %in% colnames(out_list[[i]])) out_list[[i]][[c]] <- NA_character_
  out_list[[i]] <- out_list[[i]][, CLIN_COLS, drop = FALSE]
}
combined <- do.call(rbind, out_list)
rownames(combined) <- NULL

out_file <- "Meta_AML_aggregated_clinical.tsv"
write.table(combined, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Written ", nrow(combined), " rows to ", out_file)
message("Cohort counts:")
print(table(combined$Cohort, useNA = "ifany"))
