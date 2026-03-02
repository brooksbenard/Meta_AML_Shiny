# =============================================================================
# Build a summary of cytogenetic variables across the four cohorts.
# Sources: NCRI = aml_molecular_bdp.tsv; AML-SG = paper_full_data_validation.tsv
#          + AMLSG_Karyotypes.txt (ISCN); Beat AML = beataml_wv1to4_clinical.xlsx
#          karyotype; TCGA = data_clinical_patient + sample.
# Run from Meta_AML_Shiny: Rscript build_aggregated_cytogenetics.R
# Output: Meta_AML_aggregated_cytogenetics.tsv (one row per sample)
# =============================================================================

setwd(if (exists("ofile")) dirname(ofile) else getwd())
if (!dir.exists("laml_tcga") || !dir.exists("beataml2_data") || !dir.exists("amlsg_data") || !dir.exists("UK_NCRI_data")) {
  stop("Run from Meta_AML_Shiny root.")
}

# Standard cytogenetic flag columns (0/1 or NA) used across cohorts where available
CYTO_FLAGS <- c(
  "complex", "inv_3", "inv_16", "t_15_17", "t_8_21", "t_6_9",
  "del_5", "del_7", "t_v_11", "t_9_22", "add_8", "minus_y",
  "add_11", "add_13", "add_21", "add_22", "del_9", "del_12", "del_17", "del_20"
)
OUT_COLS <- c("Patient_ID", "Cohort", "Karyotype", "Karyotype_raw", CYTO_FLAGS)

safe_num <- function(x) {
  out <- suppressWarnings(as.numeric(as.character(x)))
  out[is.na(out)] <- 0
  out
}
`%||%` <- function(x, y) if (length(x) > 0 && !is.na(x)) x else y

fill_cyto <- function(n, cohort) {
  x <- as.data.frame(matrix(NA_character_, nrow = n, ncol = length(OUT_COLS)))
  colnames(x) <- OUT_COLS
  x$Patient_ID <- NA_character_
  x$Cohort <- cohort
  x$Karyotype <- NA_character_
  x$Karyotype_raw <- NA_character_
  for (c in CYTO_FLAGS) x[[c]] <- NA_character_
  x
}

# ---------- 1. TCGA LAML ----------
# TCGA: use data_cna_hg19.seg for arm-level CNA to assign Karyotype_raw and Karyotype when available
tcga_cyto <- NULL
pat_path <- "laml_tcga/data_clinical_patient.txt"
samp_path <- "laml_tcga/data_clinical_sample.txt"
seg_path <- "laml_tcga/data_cna_hg19.seg"
if (file.exists(pat_path) && file.exists(samp_path)) {
  pat <- read.delim(pat_path, skip = 4, stringsAsFactors = FALSE, check.names = FALSE)
  samp <- read.delim(samp_path, skip = 4, stringsAsFactors = FALSE, check.names = FALSE)
  pat_id <- if ("PATIENT_ID" %in% colnames(pat)) "PATIENT_ID" else "Patient Identifier"
  samp_id <- if ("SAMPLE_ID" %in% colnames(samp)) "SAMPLE_ID" else "Sample Identifier"
  pid_samp <- if ("PATIENT_ID" %in% colnames(samp)) "PATIENT_ID" else "Patient Identifier"
  m <- merge(samp[, c(pid_samp, samp_id)], pat, by.x = pid_samp, by.y = pat_id, all.x = TRUE)
  raw_karyo <- NA_character_
  if ("CYTOGENETIC_ABNORMALITY_TYPE" %in% colnames(m)) raw_karyo <- as.character(m$CYTOGENETIC_ABNORMALITY_TYPE)
  else if ("Cytogenetic abnormality type" %in% colnames(m)) raw_karyo <- as.character(m$`Cytogenetic abnormality type`)
  raw_karyo[is.na(raw_karyo) | raw_karyo == "" | raw_karyo == "[Not Available]"] <- NA_character_
  n <- nrow(m)
  tcga_cyto <- fill_cyto(n, "TCGA")
  tcga_cyto$Patient_ID <- as.character(m[[samp_id]])
  tcga_cyto$Karyotype_raw <- raw_karyo
  simp <- rep("Other/Unknown", n)
  is_complex <- !is.na(raw_karyo) & (grepl("complex", raw_karyo, ignore.case = TRUE) | grepl(">=\\s*3", raw_karyo))
  is_normal <- !is.na(raw_karyo) & (tolower(trimws(raw_karyo)) == "normal")
  simp[is_complex] <- "Complex"
  simp[!is_complex & is_normal] <- "Normal"
  tcga_cyto$Karyotype <- simp
  # Arm-level CNA from data_cna_hg19.seg (hg19): derive Karyotype_raw and Karyotype when clinical is missing or to supplement
  if (file.exists(seg_path)) {
    seg <- read.delim(seg_path, stringsAsFactors = FALSE)
    id_col <- if ("ID" %in% colnames(seg)) "ID" else colnames(seg)[1]
    seg$chrom <- as.character(seg$chrom)
    seg$chrom[seg$chrom == "23"] <- "X"
    seg$chrom[seg$chrom == "24"] <- "Y"
    # hg19 p-arm end (approximate centromere start) in bp
    hg19_p_end <- c("1" = 125000000, "2" = 93300000, "3" = 91000000, "4" = 49700000, "5" = 46900000,
      "6" = 58900000, "7" = 58100000, "8" = 43200000, "9" = 47300000, "10" = 39800000,
      "11" = 51000000, "12" = 33200000, "13" = 16500000, "14" = 16100000, "15" = 17000000,
      "16" = 35300000, "17" = 22200000, "18" = 15400000, "19" = 24400000, "20" = 26100000,
      "21" = 10900000, "22" = 13700000, "X" = 58600000, "Y" = 10100000)
    seg_mean_col <- if ("seg.mean" %in% colnames(seg)) "seg.mean" else grep("mean|log2", colnames(seg), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(seg_mean_col)) seg_mean_col <- colnames(seg)[ncol(seg)]
    seg$mean <- suppressWarnings(as.numeric(seg[[seg_mean_col]]))
    seg$len <- seg$loc.end - seg$loc.start + 1
    p_end <- as.numeric(hg19_p_end[seg$chrom])
    p_end[is.na(p_end)] <- 0
    seg$arm <- ifelse(seg$chrom %in% c("13", "14", "15", "21", "22"), "q",
      ifelse(seg$loc.end <= p_end, "p", "q"))
    seg$arm_key <- paste0(seg$chrom, seg$arm)
    seg$weighted <- seg$mean * seg$len
    seg$seg_id <- seg[[id_col]]
    arm_agg <- aggregate(cbind(weighted = weighted, len = len) ~ seg_id + arm_key, data = seg, sum)
    arm_agg$wmean <- arm_agg$weighted / arm_agg$len
    thresh_gain <- 0.25
    thresh_loss <- -0.25
    arm_agg$state <- ifelse(arm_agg$wmean >= thresh_gain, "gain", ifelse(arm_agg$wmean <= thresh_loss, "loss", "neutral"))
    samples_seg <- unique(seg[[id_col]])
    for (i in seq_len(n)) {
      sid <- tcga_cyto$Patient_ID[i]
      if (!sid %in% samples_seg) next
      sub <- arm_agg[arm_agg$seg_id == sid & arm_agg$state != "neutral", , drop = FALSE]
      if (nrow(sub) == 0) {
        if (is.na(tcga_cyto$Karyotype_raw[i]) || tcga_cyto$Karyotype_raw[i] == "")
          tcga_cyto$Karyotype_raw[i] <- "Normal (arm-level CNA)"
        if (is.na(simp[i]) || simp[i] == "Other/Unknown") { simp[i] <- "Normal"; tcga_cyto$Karyotype[i] <- "Normal" }
        next
      }
      parts <- character(0)
      for (k in seq_len(nrow(sub))) {
        ak <- sub$arm_key[k]
        st <- sub$state[k]
        chr <- sub("p$|q$", "", ak)
        arm <- sub("^[0-9XY]+", "", ak)
        if (st == "gain") parts <- c(parts, paste0("+", chr, arm))
        else parts <- c(parts, paste0("-", chr, if (arm == "q" && chr %in% c("5", "7")) "q" else arm))
      }
      seg_karyo <- paste(parts, collapse = ",")
      if (is.na(tcga_cyto$Karyotype_raw[i]) || tcga_cyto$Karyotype_raw[i] == "")
        tcga_cyto$Karyotype_raw[i] <- paste0("arm-level CNA: ", seg_karyo)
      n_abn <- nrow(sub)
      if (n_abn >= 3) { simp[i] <- "Complex"; tcga_cyto$Karyotype[i] <- "Complex" }
      else if ((is.na(simp[i]) || simp[i] == "Other/Unknown") && n_abn == 0) { simp[i] <- "Normal"; tcga_cyto$Karyotype[i] <- "Normal" }
      else if (is.na(simp[i]) || simp[i] == "Other/Unknown") tcga_cyto$Karyotype[i] <- "Other/Unknown"
      if (n_abn >= 1) {
        if (any(grepl("5q|5p", parts))) tcga_cyto$del_5[i] <- "1"
        if (any(grepl("7q|7p|-7", parts))) tcga_cyto$del_7[i] <- "1"
        if (any(grepl("\\+8", parts))) tcga_cyto$add_8[i] <- "1"
      }
    }
  }
  # Flags from clinical raw string (when not set by seg)
  r <- tolower(raw_karyo)
  for (i in seq_len(n)) {
    if (is.na(tcga_cyto$t_8_21[i]))   tcga_cyto$t_8_21[i]   <- ifelse(!is.na(r[i]) & grepl("t\\(8;21\\)|t\\(8; 21\\)", r[i]), "1", NA_character_)
    if (is.na(tcga_cyto$inv_16[i]))   tcga_cyto$inv_16[i]   <- ifelse(!is.na(r[i]) & grepl("inv\\(16\\)|t\\(16;16\\)", r[i]), "1", NA_character_)
    if (is.na(tcga_cyto$t_15_17[i]))  tcga_cyto$t_15_17[i]  <- ifelse(!is.na(r[i]) & grepl("t\\(15;17\\)|t\\(15; 17\\)", r[i]), "1", NA_character_)
    if (is.na(tcga_cyto$inv_3[i]))   tcga_cyto$inv_3[i]    <- ifelse(!is.na(r[i]) & grepl("inv\\(3\\)|t\\(3;3\\)", r[i]), "1", NA_character_)
    if (is.na(tcga_cyto$t_6_9[i]))   tcga_cyto$t_6_9[i]    <- ifelse(!is.na(r[i]) & grepl("t\\(6;9\\)", r[i]), "1", NA_character_)
    if (is.na(tcga_cyto$del_5[i]))   tcga_cyto$del_5[i]   <- ifelse(!is.na(r[i]) & grepl("del\\(5q\\)|-5|del\\(5\\)", r[i]), "1", NA_character_)
    if (is.na(tcga_cyto$del_7[i]))   tcga_cyto$del_7[i]   <- ifelse(!is.na(r[i]) & grepl("del\\(7q\\)|-7|del\\(7\\)", r[i]), "1", NA_character_)
    if (is.na(tcga_cyto$t_v_11[i]))  tcga_cyto$t_v_11[i]  <- ifelse(!is.na(r[i]) & grepl("t\\(v;11\\)|11q23", r[i]), "1", NA_character_)
    if (is.na(tcga_cyto$t_9_22[i]))  tcga_cyto$t_9_22[i]  <- ifelse(!is.na(r[i]) & grepl("t\\(9;22\\)", r[i]), "1", NA_character_)
    if (is.na(tcga_cyto$complex[i])) tcga_cyto$complex[i]  <- ifelse(is_complex[i], "1", NA_character_)
  }
}
# %||% helper
`%||%` <- function(x, y) if (length(x) > 0 && !is.na(x)) x else y

# ---------- 2. Beat AML2 ----------
beat_cyto <- NULL
beat_path <- "beataml2_data/beataml_wv1to4_clinical.xlsx"
if (file.exists(beat_path) && requireNamespace("readxl", quietly = TRUE)) {
  b <- as.data.frame(readxl::read_excel(beat_path), stringsAsFactors = FALSE)
  sample_ids <- as.character(b$dbgap_dnaseq_sample)
  if (!"dbgap_dnaseq_sample" %in% colnames(b)) sample_ids <- rep(NA_character_, nrow(b))
  na_samp <- is.na(sample_ids) | sample_ids == ""
  if (any(na_samp) && "dbgap_rnaseq_sample" %in% colnames(b))
    sample_ids[na_samp] <- as.character(b$dbgap_rnaseq_sample[na_samp])
  keep <- !na_samp & sample_ids != ""
  b_sub <- b[keep, , drop = FALSE]
  sample_ids <- sample_ids[keep]
  uidx <- match(unique(sample_ids), sample_ids)
  b_sub <- b_sub[uidx, , drop = FALSE]
  sample_ids <- unique(sample_ids)
  n <- length(sample_ids)
  if (n > 0) {
    beat_cyto <- fill_cyto(n, "Beat AML")
    beat_cyto$Patient_ID <- sample_ids
    raw_karyo <- if ("karyotype" %in% colnames(b_sub)) as.character(b_sub$karyotype) else rep(NA_character_, n)
    raw_karyo[is.na(raw_karyo) | raw_karyo == ""] <- NA_character_
    beat_cyto$Karyotype_raw <- raw_karyo
    simp <- rep("Other/Unknown", n)
    is_normal <- !is.na(raw_karyo) & grepl("^46,(XX|XY)(\\[[0-9]+\\])?$", trimws(raw_karyo), ignore.case = TRUE)
    abnormality_count <- rep(NA_integer_, n)
    has_val <- !is.na(raw_karyo)
    if (any(has_val)) {
      s <- tolower(raw_karyo[has_val])
      markers <- gregexpr("t\\(|inv\\(|del\\(|dup\\(|add\\(|der\\(|dic\\(|ins\\(|i\\(|r\\(|\\+[0-9xy]+|-([0-9xy]+)|mar|idem|trp", s, perl = TRUE)
      abnormality_count[has_val] <- vapply(markers, function(m) if (length(m) == 1 && m[1] == -1) 0L else length(m), integer(1))
    }
    is_complex <- !is_normal & !is.na(abnormality_count) & abnormality_count >= 3
    if ("otherCytogenetics" %in% colnames(b_sub)) {
      oc <- as.character(b_sub$otherCytogenetics)
      oc[is.na(oc) | oc == ""] <- NA_character_
      is_complex <- is_complex | (!is.na(oc) & grepl("complex", oc, ignore.case = TRUE))
    }
    simp[is_complex] <- "Complex"
    simp[is_normal] <- "Normal"
    beat_cyto$Karyotype <- simp
    r <- tolower(raw_karyo)
    beat_cyto$t_8_21   <- ifelse(!is.na(r) & grepl("t\\(8;21\\)", r), "1", NA_character_)
    beat_cyto$inv_16   <- ifelse(!is.na(r) & grepl("inv\\(16\\)|t\\(16;16\\)", r), "1", NA_character_)
    beat_cyto$t_15_17  <- ifelse(!is.na(r) & grepl("t\\(15;17\\)", r), "1", NA_character_)
    beat_cyto$inv_3    <- ifelse(!is.na(r) & grepl("inv\\(3\\)|t\\(3;3\\)", r), "1", NA_character_)
    beat_cyto$t_6_9    <- ifelse(!is.na(r) & grepl("t\\(6;9\\)", r), "1", NA_character_)
    beat_cyto$del_5    <- ifelse(!is.na(r) & grepl("del\\(5|-5|5q", r), "1", NA_character_)
    beat_cyto$del_7    <- ifelse(!is.na(r) & grepl("del\\(7|-7|7q", r), "1", NA_character_)
    beat_cyto$t_v_11   <- ifelse(!is.na(r) & grepl("t\\(v;11\\)|11q23|t\\(9;11\\)", r), "1", NA_character_)
    beat_cyto$t_9_22   <- ifelse(!is.na(r) & grepl("t\\(9;22\\)", r), "1", NA_character_)
    beat_cyto$complex  <- ifelse(is_complex, "1", NA_character_)
    beat_cyto$add_8    <- ifelse(!is.na(r) & grepl("\\+8|+8", r), "1", NA_character_)
  }
}

# ---------- 3. AML-SG (paper_full_data_validation.tsv + AMLSG_Karyotypes.txt) ----------
amlsg_cyto <- NULL
val_path <- "UK_NCRI_data/data/paper_full_data_validation.tsv"
karyo_path <- "amlsg_data/AMLSG_Karyotypes.txt"
if (file.exists(val_path)) {
  lns <- readLines(val_path)
  pdids <- character(max(0, length(lns) - 1))
  for (i in seq_along(pdids)) {
    tok <- strsplit(lns[i + 1], " +", fixed = FALSE)[[1]][1]
    pdids[i] <- gsub("^\"|\"$", "", tok)
  }
  val <- read.table(val_path, header = TRUE, quote = "\"", stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, comment.char = "")
  if (nrow(val) == length(pdids)) val[[1]] <- pdids
  n <- length(pdids)
  amlsg_cyto <- fill_cyto(n, "AML-SG")
  amlsg_cyto$Patient_ID <- pdids
  # Cytogenetic columns in paper_full_data_validation (names may have parentheses)
  g <- function(nm) {
    if (!nm %in% colnames(val)) return(rep(NA_character_, n))
    x <- safe_num(val[[nm]])
    ifelse(x == 1, "1", NA_character_)
  }
  amlsg_cyto$add_8     <- g("+8")
  amlsg_cyto$add_11    <- g("+11")
  amlsg_cyto$add_13    <- g("+13")
  amlsg_cyto$add_21    <- g("+21")
  amlsg_cyto$add_22    <- g("+22")
  amlsg_cyto$complex   <- g("complex")
  amlsg_cyto$del_12    <- g("-12")
  amlsg_cyto$del_17    <- g("-17")
  amlsg_cyto$del_20    <- g("-20")
  amlsg_cyto$del_5     <- g("-5")
  amlsg_cyto$del_7     <- g("-7")
  amlsg_cyto$del_9     <- g("-9")
  amlsg_cyto$inv_16    <- g("inv(16)")
  amlsg_cyto$inv_3     <- g("inv(3)")
  amlsg_cyto$minus_y   <- g("-Y")
  amlsg_cyto$t_15_17   <- g("t(15;17)")
  amlsg_cyto$t_6_9     <- g("t(6;9)")
  amlsg_cyto$t_8_21    <- g("t(8;21)")
  amlsg_cyto$t_9_22    <- g("t(9;22)")
  amlsg_cyto$t_v_11    <- g("t(v;11)")
  # t(9;11) -> map to t_v_11 for ELN (v;11) or keep separate; paper has "t(9;11)" - leave as t_v_11 for v;11 and add t_9_11 if needed; here we only have t(v;11)
  if ("t(9;11)" %in% colnames(val)) {
    t911 <- safe_num(val[["t(9;11)"]]) == 1
    amlsg_cyto$t_v_11[!is.na(amlsg_cyto$t_v_11) & amlsg_cyto$t_v_11 == "1"] <- "1"
    amlsg_cyto$t_v_11[t911 & (is.na(amlsg_cyto$t_v_11) | amlsg_cyto$t_v_11 != "1")] <- "1"
  }
  # Karyotype simplified: Complex if complex==1, Normal if no cyto flags, else Other
  cyto_flag_cols <- c("+8", "+11", "+13", "+21", "+22", "complex", "-12", "-17", "-18", "-20", "-3", "-4", "-5", "-7", "-9", "inv(16)", "inv(3)", "-Y", "t(15;17)", "t(6;9)", "t(8;21)", "t(9;11)", "t(9;22)", "t(v;11)")
  present <- rep(FALSE, n)
  for (col in cyto_flag_cols) {
    if (col %in% colnames(val)) present <- present | (safe_num(val[[col]]) == 1)
  }
  amlsg_cyto$Karyotype <- "Other/Unknown"
  amlsg_cyto$Karyotype[!present] <- "Normal"
  if ("complex" %in% colnames(val))
    amlsg_cyto$Karyotype[safe_num(val$complex) == 1] <- "Complex"
  # Karyotype_raw from AMLSG_Karyotypes.txt (ISCN string by PDID)
  if (file.exists(karyo_path)) {
    ktab <- read.delim(karyo_path, stringsAsFactors = FALSE, check.names = FALSE)
    karo_col <- which(tolower(colnames(ktab)) == "karyotype")[1]
    pdid_col_ktab <- which(tolower(colnames(ktab)) == "pdid")[1]
    if (length(pdid_col_ktab) > 0 && length(karo_col) > 0) {
      kidx <- match(pdids, as.character(ktab[[pdid_col_ktab]]))
      amlsg_cyto$Karyotype_raw <- as.character(ktab[[karo_col]][kidx])
    } else if (ncol(ktab) >= 3) {
      kidx <- match(pdids, as.character(ktab[[2]]))
      amlsg_cyto$Karyotype_raw <- as.character(ktab[[3]][kidx])
    }
  }
  # If Karyotype_raw is no metaphases / na / no analysis, build from cytogenetic flags in AMLSG_Clinical_Anon.RData
  missing_karyo <- tolower(trimws(amlsg_cyto$Karyotype_raw)) %in% c("no metaphases", "na", "no analysis", "no material", "")
  if (any(missing_karyo) && file.exists("amlsg_data/AMLSG_Clinical_Anon.RData")) {
    load_env <- new.env()
    load("amlsg_data/AMLSG_Clinical_Anon.RData", envir = load_env)
    cl_names <- ls(load_env)
    cl <- if ("cl" %in% cl_names) get("cl", load_env) else get(cl_names[1], load_env)
    if (!is.data.frame(cl)) cl <- as.data.frame(cl)
    pdid_rdata <- as.character(cl$PDID)
    if (!"PDID" %in% colnames(cl)) pdid_rdata <- as.character(cl[[1]])
    flag_to_iscn <- list(
      t_15_17 = "t(15;17)", t_8_21 = "t(8;21)", inv16_t16_16 = "inv(16)", inv_16 = "inv(16)",
      inv3_t3_3 = "inv(3)", inv_3 = "inv(3)", minus5_5q = "del(5q)", minus7 = "-7",
      t_6_9 = "t(6;9)", t_v_11 = "t(v;11)", t_9_22 = "t(9;22)", complex = "complex"
    )
    cyto_cols_rdata <- unique(c(names(flag_to_iscn), "NK", "complex"))
    for (i in which(missing_karyo)) {
      pid <- amlsg_cyto$Patient_ID[i]
      j <- match(pid, pdid_rdata)
      if (is.na(j)) next
      parts <- character(0)
      for (col in names(flag_to_iscn)) {
        if (col %in% colnames(cl) && safe_num(cl[j, col]) == 1)
          parts <- c(parts, flag_to_iscn[[col]])
      }
      if (length(parts) > 0)
        amlsg_cyto$Karyotype_raw[i] <- paste0("derived from flags: ", paste(parts, collapse = ", "))
    }
  }
}

# ---------- 4. UK-NCRI (aml_molecular_bdp.tsv) ----------
ncri_cyto <- NULL
bdp_path <- "UK_NCRI_data/data/aml_molecular_bdp.tsv"
if (file.exists(bdp_path)) {
  bdp <- read.table(bdp_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"")
  # First column is sample ID (no header name or first col name)
  id_col <- colnames(bdp)[1]
  ncri_ids <- as.character(bdp[[id_col]])
  n <- length(ncri_ids)
  ncri_cyto <- fill_cyto(n, "UK-NCRI")
  ncri_cyto$Patient_ID <- ncri_ids
  # Cytogenetic columns in aml_molecular_bdp: add_*, del_*, minusy, t_*, inv_*, complex, others_transloc
  cyto_cols <- grep("^(add_|del_|minusy|t_|inv_|complex|others_transloc)$", colnames(bdp), value = TRUE)
  for (col in cyto_cols) {
    std <- col
    if (col == "minusy") std <- "minus_y"
    if (std %in% colnames(ncri_cyto)) {
      ncri_cyto[[std]] <- ifelse(safe_num(bdp[[col]]) == 1, "1", NA_character_)
    }
  }
  # NCRI uses t_15_17, t_8_21, etc. - same names
  for (c in c("complex", "inv_3", "inv_16", "t_15_17", "t_8_21", "t_6_9", "del_5", "del_7", "t_v_11", "t_9_22", "add_8")) {
    if (c %in% colnames(bdp)) ncri_cyto[[c]] <- ifelse(safe_num(bdp[[c]]) == 1, "1", NA_character_)
  }
  if ("minusy" %in% colnames(bdp)) ncri_cyto$minus_y <- ifelse(safe_num(bdp$minusy) == 1, "1", NA_character_)
  for (c in c("add_11", "add_13", "add_21", "add_22", "del_9", "del_12", "del_17", "del_20")) {
    if (c %in% colnames(bdp)) ncri_cyto[[c]] <- ifelse(safe_num(bdp[[c]]) == 1, "1", NA_character_)
  }
  # Karyotype: Complex if complex==1, Normal if no cyto abnormality flags, else Other
  all_cyto <- grep("^(add_|del_|minusy|t_|inv_|complex|others_transloc)$", colnames(bdp), value = TRUE)
  abnormal <- rep(FALSE, n)
  for (col in all_cyto) abnormal <- abnormal | (safe_num(bdp[[col]]) == 1)
  ncri_cyto$Karyotype <- "Other/Unknown"
  ncri_cyto$Karyotype[!abnormal] <- "Normal"
  if ("complex" %in% colnames(bdp)) ncri_cyto$Karyotype[!is.na(bdp$complex) & safe_num(bdp$complex) == 1] <- "Complex"
  ncri_cyto$Karyotype_raw <- NA_character_
  # NCRI: build Karyotype_raw and Karyotype from UK_NCRI_Cytogenetics_data.csv when available
  ncri_cyto_path <- "UK_NCRI_data/UK_NCRI_Cytogenetics_data.csv"
  if (file.exists(ncri_cyto_path)) {
    cyto_csv <- tryCatch(read.csv(ncri_cyto_path, skip = 1, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (!is.null(cyto_csv) && nrow(cyto_csv) > 0) {
      id_cyto <- if ("data_pd" %in% colnames(cyto_csv)) cyto_csv$data_pd else cyto_csv[[1]]
      id_cyto <- as.character(id_cyto)
      # Map CSV columns to ISCN-like tokens: plus8 -> +8, minus7 -> -7, del_5q -> del(5q), inv_16 -> inv(16), t_8_21 -> t(8;21)
      plus_cols <- grep("^plus[0-9]+$|^plusx$|^plusy$", colnames(cyto_csv), value = TRUE)
      minus_cols <- grep("^minus[0-9]+$|^minusx$|^minusy$", colnames(cyto_csv), value = TRUE)
      del_cols <- grep("^del_[0-9]+[pq]?$", colnames(cyto_csv), value = TRUE)
      inv_cols <- grep("^inv_[0-9]+$", colnames(cyto_csv), value = TRUE)
      t_cols <- grep("^t_[0-9]+_[0-9]+$|^t_v_11$", colnames(cyto_csv), value = TRUE)
      for (ii in seq_len(nrow(cyto_csv))) {
        pid <- id_cyto[ii]
        midx <- match(pid, ncri_ids)
        if (is.na(midx)) next
        parts <- character(0)
        for (col in plus_cols) {
          if (safe_num(cyto_csv[ii, col]) == 1) {
            num <- sub("^plus", "", col)
            parts <- c(parts, paste0("+", num))
          }
        }
        for (col in minus_cols) {
          if (safe_num(cyto_csv[ii, col]) == 1) {
            num <- sub("^minus", "", col)
            parts <- c(parts, paste0("-", num))
          }
        }
        for (col in del_cols) {
          if (safe_num(cyto_csv[ii, col]) == 1) {
            m <- regexpr("([0-9]+)([pq]?)", col)
            chr <- regmatches(col, m)
            arm <- if (grepl("q", col)) "q" else "p"
            parts <- c(parts, paste0("del(", sub("^del_", "", col), ")"))
          }
        }
        for (col in inv_cols) {
          if (safe_num(cyto_csv[ii, col]) == 1)
            parts <- c(parts, paste0("inv(", sub("^inv_", "", col), ")"))
        }
        for (col in t_cols) {
          if (safe_num(cyto_csv[ii, col]) == 1) {
            if (col == "t_v_11") parts <- c(parts, "t(v;11)")
            else parts <- c(parts, paste0("t(", gsub("_", ";", sub("^t_", "", col)), ")"))
          }
        }
        n_abn <- length(parts)
        if (n_abn > 0) {
          ncri_cyto$Karyotype_raw[midx] <- paste(parts, collapse = ",")
          if (n_abn >= 3) ncri_cyto$Karyotype[midx] <- "Complex"
          else ncri_cyto$Karyotype[midx] <- "Other/Unknown"
        } else {
          ncri_cyto$Karyotype_raw[midx] <- "Normal"
          ncri_cyto$Karyotype[midx] <- "Normal"
        }
      }
    }
  }
}

# ---------- Bind and align to aggregated clinical (one row per clinical sample) ----------
# Option A: rbind cohort tables -> one row per sample with cyto data
# Option B: start from Meta_AML_aggregated_clinical and fill cyto columns (ensures same rows)
# We do Option B so the cytogenetics file has exactly the same samples as clinical.
clin_path <- "Meta_AML_aggregated_clinical.tsv"
if (!file.exists(clin_path)) {
  message("Meta_AML_aggregated_clinical.tsv not found; outputting cytogenetics for samples with data only.")
  out_list <- list()
  if (!is.null(tcga_cyto) && nrow(tcga_cyto) > 0) out_list$tcga <- tcga_cyto
  if (!is.null(beat_cyto) && nrow(beat_cyto) > 0) out_list$beat <- beat_cyto
  if (!is.null(amlsg_cyto) && nrow(amlsg_cyto) > 0) out_list$amlsg <- amlsg_cyto
  if (!is.null(ncri_cyto) && nrow(ncri_cyto) > 0) out_list$ncri <- ncri_cyto
  combined <- if (length(out_list) > 0) do.call(rbind, out_list) else NULL
} else {
  clin <- read.delim(clin_path, stringsAsFactors = FALSE)
  combined <- fill_cyto(nrow(clin), "")
  combined$Patient_ID <- clin$Patient_ID
  combined$Cohort <- clin$Cohort
  combined$Karyotype <- clin$Karyotype
  combined$Karyotype_raw <- NA_character_
  for (co in c("TCGA", "Beat AML", "AML-SG", "UK-NCRI")) {
    idx <- which(clin$Cohort == co)
    if (length(idx) == 0) next
    if (co == "TCGA" && !is.null(tcga_cyto)) {
      midx <- match(clin$Patient_ID[idx], tcga_cyto$Patient_ID)
      for (col in c("Karyotype_raw", CYTO_FLAGS)) {
        if (col %in% colnames(tcga_cyto))
          combined[[col]][idx] <- tcga_cyto[[col]][midx]
      }
      if ("Karyotype" %in% colnames(tcga_cyto)) combined$Karyotype[idx] <- tcga_cyto$Karyotype[midx]
    }
    if (co == "Beat AML" && !is.null(beat_cyto)) {
      midx <- match(clin$Patient_ID[idx], beat_cyto$Patient_ID)
      for (col in c("Karyotype_raw", CYTO_FLAGS)) {
        if (col %in% colnames(beat_cyto))
          combined[[col]][idx] <- beat_cyto[[col]][midx]
      }
      if ("Karyotype" %in% colnames(beat_cyto)) combined$Karyotype[idx] <- beat_cyto$Karyotype[midx]
    }
    if (co == "AML-SG" && !is.null(amlsg_cyto)) {
      midx <- match(clin$Patient_ID[idx], amlsg_cyto$Patient_ID)
      for (col in c("Karyotype_raw", CYTO_FLAGS)) {
        if (col %in% colnames(amlsg_cyto))
          combined[[col]][idx] <- amlsg_cyto[[col]][midx]
      }
      if ("Karyotype" %in% colnames(amlsg_cyto)) combined$Karyotype[idx] <- amlsg_cyto$Karyotype[midx]
    }
    if (co == "UK-NCRI" && !is.null(ncri_cyto)) {
      midx <- match(clin$Patient_ID[idx], ncri_cyto$Patient_ID)
      for (col in c("Karyotype_raw", CYTO_FLAGS)) {
        if (col %in% colnames(ncri_cyto))
          combined[[col]][idx] <- ncri_cyto[[col]][midx]
      }
      if ("Karyotype" %in% colnames(ncri_cyto)) combined$Karyotype[idx] <- ncri_cyto$Karyotype[midx]
    }
  }
}

if (is.null(combined) || nrow(combined) == 0) stop("No cytogenetic data to write.")

combined <- combined[, OUT_COLS, drop = FALSE]
out_file <- "Meta_AML_aggregated_cytogenetics.tsv"
write.table(combined, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Written ", nrow(combined), " rows to ", out_file)
message("Cohort counts:")
print(table(combined$Cohort, useNA = "ifany"))
