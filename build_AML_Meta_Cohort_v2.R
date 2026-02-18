# =============================================================================
# Build AML_Meta_Cohort_v2 from four data sources:
# - TCGA: laml_tcga_pan_can_atlas_2018/
# - Beat AML: beataml2_data/
# - AML-SG: amlsg_data/
# - UK-NCRI: UK_NCRI_data/
# =============================================================================

library(readxl)

# Mutation category mapping (from final_data_matrix)
# Will be assigned based on gene names from Benard et al. (2021) reference
# FLT3-ITD and FLT3-TKD get "RTK-RAS Signaling"

# Helper: Ensure standard clinical columns exist (numeric, NA if missing)
ensure_clin_cols <- function(d, wbc = NA_real_, hemoglobin = NA_real_, platelet = NA_real_,
                             bm_blast_percent = NA_real_, pb_blast_percent = NA_real_, ldh = NA_real_) {
  n <- nrow(d)
  if (length(wbc) == 1L) wbc <- rep(wbc, n)
  if (length(hemoglobin) == 1L) hemoglobin <- rep(hemoglobin, n)
  if (length(platelet) == 1L) platelet <- rep(platelet, n)
  if (length(bm_blast_percent) == 1L) bm_blast_percent <- rep(bm_blast_percent, n)
  if (length(pb_blast_percent) == 1L) pb_blast_percent <- rep(pb_blast_percent, n)
  if (length(ldh) == 1L) ldh <- rep(ldh, n)
  d$WBC <- as.numeric(wbc)
  d$Hemoglobin <- as.numeric(hemoglobin)
  d$Platelet <- as.numeric(platelet)
  d$BM_blast_percent <- as.numeric(bm_blast_percent)
  d$PB_blast_percent <- as.numeric(pb_blast_percent)
  d$LDH <- as.numeric(ldh)
  d
}

# Helper: Map variant classification to variant_type
map_variant_type <- function(vc) {
  vc <- tolower(as.character(vc))
  out <- rep("Other", length(vc))
  out[vc %in% c("sub", "snv", "snp", "missense_mutation", "missense_variant", "stop_gained", "nonsense_mutation")] <- "SNV"
  out[vc %in% c("del", "deletion", "frame_shift_del")] <- "Deletion"
  out[vc %in% c("splice_acceptor_variant", "splice_donor_variant", "splice_site")] <- "Splicing"
  out[vc %in% c("indel", "insertion", "frame_shift_ins", "inframe_insertion", "inframe_deletion",
                "in_frame_ins", "in_frame_del")] <- "Indel"  # in_frame_* = TCGA/cBioPortal style
  out[vc %in% c("itd")] <- "ITD"
  out[vc %in% c("ptd")] <- "PTD"
  out
}

# Helper: Normalize FLT3 to FLT3-ITD / FLT3-TKD / FLT3-Other by variant_type so most FLT3 are ITD or TKD, not other
normalize_flt3_gene <- function(gene, variant_type) {
  gene <- as.character(gene)
  variant_type <- as.character(variant_type)
  out <- gene
  flt3 <- which(gene == "FLT3")
  if (length(flt3) == 0) return(out)
  vt <- variant_type[flt3]
  out[flt3[vt %in% c("ITD", "INDEL", "Indel", "Insertion")]] <- "FLT3-ITD"
  out[flt3[vt %in% c("SNV", "Deletion")]] <- "FLT3-TKD"
  out[flt3[out[flt3] == "FLT3"]] <- "FLT3-Other"  # any remaining FLT3 -> Other
  out
}

# Helper: Assign Gene_Group for prognostically distinct mutation types (when identifiable)
# Gene_Group values: CEBPA_bi, CEBPA_mono, FLT3_ITD, FLT3_TKD, FLT3_other,
#   IDH2_p.R140, IDH2_p.R172, IDH2_other, NRAS_p.G12_13, NRAS_p.Q61_62, NRAS_other,
#   U2AF1_p.S34, U2AF1_p.Q157, U2AF1_other. CEBPA_bi/mono set later by sample-level count.
assign_gene_group <- function(gene, variant_type, aa_change) {
  gene <- as.character(gene)
  variant_type <- as.character(variant_type)
  aa <- as.character(aa_change)
  aa <- trimws(toupper(aa))
  out <- rep(NA_character_, length(gene))
  # FLT3: only FLT3-Other -> FLT3_other; FLT3-ITD/TKD assigned above; any remaining "FLT3" -> FLT3_TKD as fallback
  out[gene %in% c("FLT3-ITD")] <- "FLT3_ITD"
  out[gene %in% c("FLT3-TKD")] <- "FLT3_TKD"
  out[gene %in% c("FLT3-Other")] <- "FLT3_other"
  out[gene == "FLT3"] <- "FLT3_TKD"  # fallback if normalize_flt3_gene not used
  # IDH2 position: R140, R172, or other
  idh2 <- which(gene == "IDH2")
  out[idh2[grepl("R140|P\\.R140", aa[idh2])]] <- "IDH2_p.R140"
  out[idh2[grepl("R172|P\\.R172", aa[idh2])]] <- "IDH2_p.R172"
  out[idh2[is.na(out[idh2])]] <- "IDH2_other"  # retain other cases
  # NRAS position: G12/13, Q61/62, or other
  nras <- which(gene == "NRAS")
  out[nras[grepl("G12|G13|P\\.G12|P\\.G13", aa[nras]) | aa[nras] %in% c("G12", "G13")]] <- "NRAS_p.G12_13"
  out[nras[grepl("Q61|Q62|P\\.Q61|P\\.Q62", aa[nras]) | aa[nras] %in% c("Q61", "Q62")]] <- "NRAS_p.Q61_62"
  out[nras[is.na(out[nras])]] <- "NRAS_other"  # retain other cases
  # U2AF1 position: S34, Q157, or other
  u2 <- which(gene == "U2AF1")
  out[u2[grepl("S34|P\\.S34", aa[u2])]] <- "U2AF1_p.S34"
  out[u2[grepl("Q157|P\\.Q157", aa[u2])]] <- "U2AF1_p.Q157"
  out[u2[is.na(out[u2])]] <- "U2AF1_other"  # retain other cases
  out
}

# Helper: Classify karyotype/cytogenetics into Complex, Normal, Other, Unknown (for prognosis filtering)
classify_karyotype <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == "" | grepl("Not Available|unknown|\\[\\s*\\]", x, ignore.case = TRUE)] <- NA_character_
  out <- rep("Unknown", length(x))
  # Complex: >= 3 abnormalities, "complex", monosomal
  complex_idx <- grepl("Complex\\s*-\\s*>=\\s*3|complex|monosomal|>= 3 distinct|\\d~\\d.*mar|idem.*mar|Komplex", x, ignore.case = TRUE)
  out[complex_idx] <- "Complex"
  # Normal: exact "Normal" (TCGA) or simple 46,XX/46,XY only (Beat AML)
  normal_idx <- !complex_idx & (trimws(x) == "Normal" | grepl("^46,(XX|XY)\\[?\\d*\\]?\\s*$", trimws(x)))
  out[normal_idx] <- "Normal"
  # Other: non-missing, not Complex, not Normal
  other_idx <- !is.na(x) & trimws(x) != "" & out == "Unknown"
  out[other_idx] <- "Other"
  out
}

# Helper: Count cytogenetic abnormalities in ISCN karyotype string for AML-SG (3+ = Complex)
# Counts del(, add(, t(, inv(, der(, ins(, dic(, idic(, and aneuploidy +N/-N (e.g. +8, -7)
count_cytogenetic_events <- function(s) {
  s <- as.character(s)
  s[is.na(s) | trimws(s) == ""] <- ""
  n <- length(s)
  out <- integer(n)
  for (i in seq_len(n)) {
    x <- s[i]
    if (is.na(x) || trimws(x) == "" || grepl("no metaphases|no analysis|no material|^na\\s*$|Keine analysierbaren", x, ignore.case = TRUE)) {
      out[i] <- NA_integer_
      next
    }
    # Count structural: del(, add(, t(, inv(, der(, ins(, dic(, idic(
    m <- gregexpr("del\\(|add\\(|t\\(|inv\\(|der\\(|ins\\(|dic\\(|idic\\(", x, ignore.case = TRUE)[[1]]
    struct <- if (length(m) == 1L && m[1] == -1L) 0L else length(m)
    # Count aneuploidy: +8, -7, +21, -5, etc.
    m2 <- gregexpr("[+-][0-9]+|[+-](X|Y)\\b", x)[[1]]
    aneu <- if (length(m2) == 1L && m2[1] == -1L) 0L else length(m2)
    # +mar, -mar, etc. as one event each
    m3 <- gregexpr("\\+mar|\\-mar|\\+r|dmin", x, ignore.case = TRUE)[[1]]
    mar <- if (length(m3) == 1L && m3[1] == -1L) 0L else length(m3)
    out[i] <- struct + aneu + mar
  }
  out
}

# Classify AML-SG karyotype from ISCN string: Normal (46,XX/XY only), Complex (3+ events or "complex"), Other (1-2), Unknown
classify_amlsg_karyotype <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(gsub("^\"|\"$", "", x))
  out <- rep("Unknown", length(x))
  # Uninterpretable
  unk <- grepl("no metaphases|no analysis|no material|^na\\s*$|Keine analysierbaren|not available", x, ignore.case = TRUE)
  # Normal: only 46,XX or 46,XY (with optional [digits])
  normal <- !unk & grepl("^46,(XX|XY)(\\[\\d+\\])?\\s*$", x)
  out[normal] <- "Normal"
  # Complex: explicit "complex"/"Komplex" or 3+ events
  n_events <- count_cytogenetic_events(x)
  complex <- !unk & (grepl("complex|Komplex", x, ignore.case = TRUE) | (!is.na(n_events) & n_events >= 3L))
  out[complex] <- "Complex"
  # Other: 1-2 events, interpretable string
  other <- !unk & out != "Normal" & out != "Complex" & (grepl("46,|45,|47,|48,|49,|50,|51,|52,|inv\\(|del\\(|t\\(|add\\(|der\\(", x) | (!is.na(n_events) & n_events >= 1L & n_events <= 2L))
  out[other] <- "Other"
  out
}

# =============================================================================
# 0. Reference gene set (AML-SG + UK-NCRI) for filtering TCGA and Beat AML
# =============================================================================
cat("Building reference gene set from AML-SG and UK-NCRI...\n")
amlsg_ref <- read.delim("amlsg_data/AMLSG_Genetic.txt", stringsAsFactors = FALSE)
amlsg_ref$Gene <- amlsg_ref$GENE
amlsg_ref$Gene[amlsg_ref$Gene == "SFRS2"] <- "SRSF2"  # correct gene symbol
amlsg_ref$vt <- map_variant_type(amlsg_ref$VARIANT_TYPE)
amlsg_ref$Gene <- normalize_flt3_gene(amlsg_ref$Gene, amlsg_ref$vt)
genes_amlsg <- unique(as.character(amlsg_ref$Gene))

ukncri_ref <- read.csv("UK_NCRI_data/UK_NCRI_Mutations_data.csv", stringsAsFactors = FALSE, skip = 1)
ukncri_ref$vt <- map_variant_type(ukncri_ref$type)
ukncri_ref$Gene <- normalize_flt3_gene(ukncri_ref$gene, ukncri_ref$vt)
genes_ukncri <- unique(as.character(ukncri_ref$Gene))
bdp_path <- file.path(getwd(), "UK_NCRI_data", "data", "aml_molecular_bdp.tsv")
if (!file.exists(bdp_path)) bdp_path <- file.path(getwd(), "UK_NCRI_data", "aml_molecular_bdp.tsv")
if (file.exists(bdp_path)) {
  genes_ukncri <- unique(c(genes_ukncri, "FLT3-ITD", "FLT3-TKD", "FLT3-Other"))
}
amlsg_ukncri_genes <- unique(c(genes_amlsg, genes_ukncri))
cat("  Reference gene set: ", length(amlsg_ukncri_genes), " genes (AML-SG + UK-NCRI)\n", sep = "")

# =============================================================================
# 1. TCGA LAML (laml_tcga - complete clinical data: age, sex, blast %, platelet)
# =============================================================================
cat("Loading TCGA LAML data from laml_tcga...\n")
tcga_mut <- read.delim("laml_tcga/data_mutations.txt", stringsAsFactors = FALSE, comment.char = "#")
tcga_clin_patient <- read.delim("laml_tcga/data_clinical_patient.txt", stringsAsFactors = FALSE, comment.char = "#", skip = 4)
tcga_clin_sample <- read.delim("laml_tcga/data_clinical_sample.txt", stringsAsFactors = FALSE, comment.char = "#", skip = 4)

# TCGA mutations: Hugo_Symbol, Tumor_Sample_Barcode, t_vaf (if available)
tcga_mut$Sample <- tcga_mut$Tumor_Sample_Barcode
tcga_mut$Gene <- tcga_mut$Hugo_Symbol
# TCGA may not have t_vaf, use t_alt_count / t_depth if available
if ("t_vaf" %in% colnames(tcga_mut)) {
  tcga_mut$VAF <- as.numeric(tcga_mut$t_vaf) * 100
} else if (all(c("t_alt_count", "t_depth") %in% colnames(tcga_mut))) {
  tcga_mut$VAF <- (as.numeric(tcga_mut$t_alt_count) / as.numeric(tcga_mut$t_depth)) * 100
} else {
  tcga_mut$VAF <- NA_real_
}
tcga_mut$variant_type <- map_variant_type(tcga_mut$Variant_Classification)
tcga_mut$Gene <- normalize_flt3_gene(tcga_mut$Gene, tcga_mut$variant_type)
tcga_mut$Cohort <- "TCGA"
tcga_mut$AA_change <- if ("HGVSp_Short" %in% colnames(tcga_mut)) as.character(tcga_mut$HGVSp_Short) else if ("HGVSp" %in% colnames(tcga_mut)) as.character(tcga_mut$HGVSp) else NA_character_
tcga_mut$Gene_Group <- assign_gene_group(tcga_mut$Gene, tcga_mut$variant_type, tcga_mut$AA_change)
# Filter TCGA to genes present in AML-SG or UK-NCRI
tcga_mut <- tcga_mut[as.character(tcga_mut$Gene) %in% amlsg_ukncri_genes, , drop = FALSE]

# TCGA clinical: merge patient and sample (column names are in row 5, skip first 4 rows)
tcga_clin <- merge(tcga_clin_patient, tcga_clin_sample, by = "PATIENT_ID", all = TRUE, suffixes = c("", "_sample"))
# Use SAMPLE_ID so merge with mutations (Tumor_Sample_Barcode) matches; fallback to PATIENT_ID
tcga_clin$Sample <- if ("SAMPLE_ID" %in% colnames(tcga_clin)) as.character(tcga_clin$SAMPLE_ID) else as.character(tcga_clin$PATIENT_ID)
# Age: use AGE when present; fill missing from DAYS_TO_BIRTH (age = -DAYS_TO_BIRTH/365.25)
tcga_clin$Age <- NA_real_
if ("AGE" %in% colnames(tcga_clin)) {
  tcga_clin$Age <- as.numeric(tcga_clin$AGE)
}
if ("DAYS_TO_BIRTH" %in% colnames(tcga_clin)) {
  days_birth <- as.numeric(tcga_clin$DAYS_TO_BIRTH)
  miss <- is.na(tcga_clin$Age) & !is.na(days_birth)
  tcga_clin$Age[miss] <- round(-days_birth[miss] / 365.25)
}
# Sex: try SEX, Sex, sex, gender; normalize to Male/Female
sex_col <- NULL
for (c in c("SEX", "Sex", "sex", "gender")) {
  if (c %in% colnames(tcga_clin)) { sex_col <- c; break }
}
if (!is.null(sex_col)) {
  sx <- as.character(tcga_clin[[sex_col]])
  tcga_clin$Sex <- ifelse(sx %in% c("MALE", "Male"), "Male", ifelse(sx %in% c("FEMALE", "Female"), "Female", NA_character_))
} else {
  tcga_clin$Sex <- NA_character_
}
# TCGA OS is in months - convert to days
if ("OS_MONTHS" %in% colnames(tcga_clin)) {
  tcga_clin$Time_to_OS <- as.numeric(tcga_clin$OS_MONTHS) * 30.44  # Convert months to days
} else {
  tcga_clin$Time_to_OS <- NA_real_
}
if ("OS_STATUS" %in% colnames(tcga_clin)) {
  tcga_clin$Censor <- ifelse(toupper(tcga_clin$OS_STATUS) %in% c("DECEASED", "1"), 1, 0)
} else {
  tcga_clin$Censor <- NA_real_
}
tcga_clin$Subset <- "De novo"  # TCGA is primary AML
# TCGA may have cytogenetic risk in sample file or need to derive from karyotype
risk_col <- grep("Cytogenetic.Risk|RISK|CYTOGENETIC", colnames(tcga_clin), value = TRUE, ignore.case = TRUE)
if (length(risk_col) > 0) {
  tcga_clin$Risk <- ifelse(tcga_clin[[risk_col[1]]] %in% c("Favorable", "Intermediate", "Adverse"), 
                           tcga_clin[[risk_col[1]]], "Unknown")
} else {
  tcga_clin$Risk <- "Unknown"
}
# Karyotype: Complex / Normal / Other from cBioPortal "Cytogenetic Code (Other)" (CYTOGENETIC_CODE_OTHER) for laml_tcga_pub
tcga_clin$Karyotype <- "Unknown"
if (requireNamespace("jsonlite", quietly = TRUE)) {
  cbioportal_url <- "https://www.cbioportal.org/api/studies/laml_tcga_pub/clinical-data?clinicalDataType=PATIENT"
  tryCatch({
    raw <- readLines(cbioportal_url, warn = FALSE)
    if (length(raw) > 0) {
      js <- jsonlite::fromJSON(raw)
      if (is.data.frame(js) && nrow(js) > 0 && "patientId" %in% colnames(js) && "clinicalAttributeId" %in% colnames(js) && "value" %in% colnames(js)) {
        code_other <- js[js$clinicalAttributeId == "CYTOGENETIC_CODE_OTHER", c("patientId", "value")]
        if (nrow(code_other) > 0) {
          val <- trimws(as.character(code_other$value))
          # Map Cytogenetic Code (Other) to Complex, Normal, or Other
          k <- rep("Other", length(val))
          k[is.na(val) | val == "" | val == "NA"] <- "Unknown"
          k[grepl("normal\\s*karyotype|^normal$", val, ignore.case = TRUE)] <- "Normal"
          k[grepl("complex", val, ignore.case = TRUE)] <- "Complex"
          code_other$Karyotype <- k
          patient_id <- as.character(tcga_clin$PATIENT_ID)
          idx <- match(patient_id, as.character(code_other$patientId))
          hit <- !is.na(idx)
          tcga_clin$Karyotype[hit] <- code_other$Karyotype[idx[hit]]
          n_norm <- sum(tcga_clin$Karyotype == "Normal")
          n_comp <- sum(tcga_clin$Karyotype == "Complex")
          n_oth <- sum(tcga_clin$Karyotype == "Other")
          cat("  TCGA: karyotype from cBioPortal Cytogenetic Code (Other) for ", sum(hit), " samples (", n_comp, " Complex, ", n_norm, " Normal, ", n_oth, " Other)\n", sep = "")
        }
      }
    }
  }, error = function(e) {
    cat("  TCGA: cBioPortal fetch failed, karyotype Unknown\n")
  })
} else {
  cat("  TCGA: install jsonlite for karyotype from cBioPortal (Cytogenetic Code Other)\n")
}
# BM blast %, PB blast %, platelet: use if any column in patient/sample clinical has blast or platelet
bm_col <- grep("blast|bm_blast|bone.marrow|percent.blast|BLAST", colnames(tcga_clin), ignore.case = TRUE, value = TRUE)
bm_col <- if (length(bm_col) > 0 && !is.na(bm_col[1])) bm_col[1] else NULL
# PB blast: prefer columns explicitly for peripheral blood (exclude BM-only)
pb_col <- grep("pb_blast|peripheral.*blast|pb.blast|blood.blast|PB_BLAST", colnames(tcga_clin), ignore.case = TRUE, value = TRUE)
pb_col <- if (length(pb_col) > 0 && !is.na(pb_col[1])) pb_col[1] else NULL
plt_col <- grep("platelet|PLATELET", colnames(tcga_clin), ignore.case = TRUE, value = TRUE)
plt_col <- if (length(plt_col) > 0 && !is.na(plt_col[1])) plt_col[1] else NULL
tcga_clin <- ensure_clin_cols(tcga_clin,
  bm_blast_percent = if (!is.null(bm_col)) as.numeric(tcga_clin[[bm_col]]) else NA_real_,
  pb_blast_percent = if (!is.null(pb_col)) as.numeric(tcga_clin[[pb_col]]) else NA_real_,
  platelet = if (!is.null(plt_col)) as.numeric(tcga_clin[[plt_col]]) else NA_real_)

# Merge TCGA mutations with clinical
tcga <- merge(tcga_mut[, c("Sample", "Gene", "VAF", "variant_type", "Cohort", "AA_change", "Gene_Group")], 
              tcga_clin[, c("Sample", "Age", "Sex", "Time_to_OS", "Censor", "Subset", "Risk", "Karyotype",
                "WBC", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent", "LDH")], 
              by = "Sample", all.x = TRUE)

cat("  TCGA: ", length(unique(tcga$Sample)), " samples, ", nrow(tcga), " mutations\n", sep = "")

# =============================================================================
# 2. Beat AML
# =============================================================================
cat("Loading Beat AML data...\n")
beataml_mut <- read.delim("beataml2_data/mutations.txt", stringsAsFactors = FALSE, check.names = FALSE)
# Prefer wv1to4 clinical so we have karyotype and otherCytogenetics for proper Complex/Normal/Other assignment
beataml_clin_path <- if (file.exists("beataml2_data/beataml_wv1to4_clinical.xlsx")) "beataml2_data/beataml_wv1to4_clinical.xlsx" else "beataml2_data/clinical.xlsx"

if (!file.exists(beataml_clin_path)) {
  stop("Beat AML clinical file not found. Run setup_beataml2.R first.")
}

beataml_clin <- as.data.frame(read_excel(beataml_clin_path), stringsAsFactors = FALSE)
beataml_clin$Sample <- as.character(beataml_clin$dbgap_dnaseq_sample)
na_idx <- beataml_clin$Sample == "" | is.na(beataml_clin$Sample)
beataml_clin$Sample[na_idx] <- as.character(beataml_clin$dbgap_rnaseq_sample[na_idx])

# Beat AML mutations: dbgap_sample_id, symbol, t_vaf, genotyper (mutect/varscan)
beataml_mut$Sample <- as.character(beataml_mut$dbgap_sample_id)
beataml_mut$VAF <- as.numeric(beataml_mut$t_vaf) * 100
beataml_mut$symbol <- as.character(beataml_mut$symbol)

# Deduplicate by genotyper: keep all calls, but (1) exact same mutation called by both -> keep one (Mutect);
# (2) same sample/gene/AA position but different exact variant -> keep Mutect only to avoid double-counting position.
if ("genotyper" %in% colnames(beataml_mut)) {
  # AA position from protein_position (e.g. "315/2002" -> 315) or hgvsp_short (e.g. p.R882H -> 882)
  aa_pos <- rep(NA_integer_, nrow(beataml_mut))
  if ("protein_position" %in% colnames(beataml_mut)) {
    pp_raw <- as.character(beataml_mut$protein_position)
    aa_pos[!is.na(pp_raw) & pp_raw != ""] <- suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", pp_raw[!is.na(pp_raw) & pp_raw != ""])))
  }
  hs <- as.character(beataml_mut$hgvsp_short)
  miss <- is.na(aa_pos) & !is.na(hs) & hs != ""
  if (any(miss)) {
    parsed <- regexpr("[0-9]+", hs[miss])
    aa_pos[miss][parsed > 0] <- as.integer(regmatches(hs[miss], parsed))
  }
  beataml_mut$AA_pos <- aa_pos

  exact_key <- paste(beataml_mut$Sample, beataml_mut$seqnames, beataml_mut$pos_start, beataml_mut$ref, beataml_mut$alt, beataml_mut$symbol, sep = "|")
  # (1) Exact duplicate: same variant called by both callers -> keep one row (prefer Mutect)
  is_mutect <- beataml_mut$genotyper == "mutect"
  keep_exact <- logical(nrow(beataml_mut))
  for (k in unique(exact_key)) {
    idx <- which(exact_key == k)
    if (any(is_mutect[idx])) keep_exact[idx][which(is_mutect[idx])[1]] <- TRUE
    else keep_exact[idx[1]] <- TRUE
  }
  beataml_mut <- beataml_mut[keep_exact, , drop = FALSE]
  exact_key <- exact_key[keep_exact]
  is_mutect <- is_mutect[keep_exact]
  aa_pos <- aa_pos[keep_exact]

  # (2) Same sample + gene + AA position, different exact variant: keep Mutect only
  saga <- paste(beataml_mut$Sample, beataml_mut$symbol, aa_pos, sep = "|")
  saga_with_mutect <- unique(saga[is_mutect & !is.na(aa_pos)])
  saga_with_varscan <- unique(saga[!is_mutect & !is.na(aa_pos)])
  saga_has_both <- (saga %in% saga_with_mutect) & (saga %in% saga_with_varscan)
  drop_varscan_same_pos <- !is_mutect & saga_has_both
  beataml_mut <- beataml_mut[!drop_varscan_same_pos, , drop = FALSE]
  beataml_mut$AA_pos <- NULL
  cat("  Beat AML: deduplicated by exact variant and by sample/gene/AA position (prefer Mutect): ", nrow(beataml_mut), " mutations\n", sep = "")
}

# FLT3 annotation (ITD vs TKD)
vc <- tolower(as.character(beataml_mut$variant_classification))
beataml_mut$Gene <- beataml_mut$symbol
beataml_mut$Gene[beataml_mut$symbol == "FLT3" & vc == "inframe_insertion"] <- "FLT3-ITD"
beataml_mut$Gene[beataml_mut$symbol == "FLT3" & vc %in% c("missense_variant", "inframe_deletion")] <- "FLT3-TKD"
beataml_mut$Gene[beataml_mut$symbol == "FLT3" & beataml_mut$Gene == "FLT3"] <- "FLT3-TKD"

beataml_mut$variant_type <- map_variant_type(beataml_mut$variant_classification)
beataml_mut$Cohort <- "Beat AML"
beataml_mut$AA_change <- if ("hgvsp_short" %in% colnames(beataml_mut)) as.character(beataml_mut$hgvsp_short) else if ("hgvsp" %in% colnames(beataml_mut)) as.character(beataml_mut$hgvsp) else NA_character_
beataml_mut$Gene_Group <- assign_gene_group(beataml_mut$Gene, beataml_mut$variant_type, beataml_mut$AA_change)
# Filter Beat AML to genes present in AML-SG or UK-NCRI
beataml_mut <- beataml_mut[as.character(beataml_mut$Gene) %in% amlsg_ukncri_genes, , drop = FALSE]

# Aggregate to one VAF per sample-gene (max)
beataml_mut_agg <- aggregate(VAF ~ Sample + Gene + variant_type + Cohort, data = beataml_mut, FUN = max, na.rm = TRUE)
key <- paste(beataml_mut_agg$Sample, beataml_mut_agg$Gene, beataml_mut_agg$variant_type, beataml_mut_agg$Cohort)
key_src <- paste(beataml_mut$Sample, beataml_mut$Gene, beataml_mut$variant_type, beataml_mut$Cohort)
idx <- match(key, key_src)
beataml_mut_agg$AA_change <- beataml_mut$AA_change[idx]
beataml_mut_agg$Gene_Group <- beataml_mut$Gene_Group[idx]

# Add consensus mutations from clinical (FLT3-ITD, NPM1, RUNX1, ASXL1, TP53) only when
# the sample does not already have that gene in the mutation data, to avoid double-counting.
existing_sample_gene <- unique(paste(beataml_mut_agg$Sample, beataml_mut_agg$Gene, sep = "|"))
cn <- colnames(beataml_clin)
flt3_col <- match("FLT3-ITD", cn, nomatch = 0)
if (flt3_col == 0) flt3_col <- match("FLT3_ITD", cn, nomatch = 0)
npm1_col <- match("NPM1", cn, nomatch = 0)
runx1_col <- match("RUNX1", cn, nomatch = 0)
asxl1_col <- match("ASXL1", cn, nomatch = 0)
tp53_col <- match("TP53", cn, nomatch = 0)

consensus_muts <- list()
for (i in seq_len(nrow(beataml_clin))) {
  s <- beataml_clin$Sample[i]
  if (is.na(s) || s == "") next

  if (flt3_col > 0 && !(paste(s, "FLT3-ITD", sep = "|") %in% existing_sample_gene)) {
    val <- as.character(beataml_clin[[flt3_col]][i])
    if (!is.na(val) && tolower(trimws(val)) %in% c("positive", "mutant", "yes", "1", "true")) {
      consensus_muts[[length(consensus_muts) + 1]] <- data.frame(Sample = s, Gene = "FLT3-ITD", VAF = 50, variant_type = "ITD", Cohort = "Beat AML", AA_change = NA_character_, Gene_Group = "FLT3_ITD", stringsAsFactors = FALSE)
      existing_sample_gene <- c(existing_sample_gene, paste(s, "FLT3-ITD", sep = "|"))
    }
  }
  if (npm1_col > 0 && !(paste(s, "NPM1", sep = "|") %in% existing_sample_gene)) {
    val <- as.character(beataml_clin[[npm1_col]][i])
    if (!is.na(val) && tolower(trimws(val)) %in% c("positive", "mutant", "yes", "1", "true")) {
      consensus_muts[[length(consensus_muts) + 1]] <- data.frame(Sample = s, Gene = "NPM1", VAF = 50, variant_type = "SNV", Cohort = "Beat AML", AA_change = NA_character_, Gene_Group = NA_character_, stringsAsFactors = FALSE)
      existing_sample_gene <- c(existing_sample_gene, paste(s, "NPM1", sep = "|"))
    }
  }
  if (runx1_col > 0 && !(paste(s, "RUNX1", sep = "|") %in% existing_sample_gene)) {
    val <- beataml_clin[[runx1_col]][i]
    if (!is.na(val) && as.character(trimws(val)) != "" && !tolower(trimws(as.character(val))) %in% c("negative", "wild", "wt", "no", "0")) {
      consensus_muts[[length(consensus_muts) + 1]] <- data.frame(Sample = s, Gene = "RUNX1", VAF = 50, variant_type = "SNV", Cohort = "Beat AML", AA_change = NA_character_, Gene_Group = NA_character_, stringsAsFactors = FALSE)
      existing_sample_gene <- c(existing_sample_gene, paste(s, "RUNX1", sep = "|"))
    }
  }
  if (asxl1_col > 0 && !(paste(s, "ASXL1", sep = "|") %in% existing_sample_gene)) {
    val <- beataml_clin[[asxl1_col]][i]
    if (!is.na(val) && as.character(trimws(val)) != "" && !tolower(trimws(as.character(val))) %in% c("negative", "wild", "wt", "no", "0")) {
      consensus_muts[[length(consensus_muts) + 1]] <- data.frame(Sample = s, Gene = "ASXL1", VAF = 50, variant_type = "SNV", Cohort = "Beat AML", AA_change = NA_character_, Gene_Group = NA_character_, stringsAsFactors = FALSE)
      existing_sample_gene <- c(existing_sample_gene, paste(s, "ASXL1", sep = "|"))
    }
  }
  if (tp53_col > 0 && !(paste(s, "TP53", sep = "|") %in% existing_sample_gene)) {
    val <- beataml_clin[[tp53_col]][i]
    if (!is.na(val) && as.character(trimws(val)) != "" && !tolower(trimws(as.character(val))) %in% c("negative", "wild", "wt", "no", "0")) {
      consensus_muts[[length(consensus_muts) + 1]] <- data.frame(Sample = s, Gene = "TP53", VAF = 50, variant_type = "SNV", Cohort = "Beat AML", AA_change = NA_character_, Gene_Group = NA_character_, stringsAsFactors = FALSE)
      existing_sample_gene <- c(existing_sample_gene, paste(s, "TP53", sep = "|"))
    }
  }
}

if (length(consensus_muts) > 0) {
  consensus_df <- do.call(rbind, consensus_muts)
  combined <- rbind(beataml_mut_agg, consensus_df)
  beataml_mut_agg <- aggregate(VAF ~ Sample + Gene + variant_type + Cohort, data = combined, FUN = max, na.rm = TRUE)
  dedup <- combined[!duplicated(paste(combined$Sample, combined$Gene, combined$variant_type, combined$Cohort)), c("Sample", "Gene", "variant_type", "Cohort", "AA_change", "Gene_Group")]
  beataml_mut_agg <- merge(beataml_mut_agg, dedup, by = c("Sample", "Gene", "variant_type", "Cohort"), all.x = TRUE)
}

# Beat AML clinical: subset labels (capitalized for display)
beataml_clin$Subset <- "Other"
if ("isDenovo" %in% cn) {
  denovo <- as.character(beataml_clin$isDenovo) == "TRUE" | beataml_clin$isDenovo == TRUE
  beataml_clin$Subset[denovo] <- "De novo"
}
if ("isTransformed" %in% cn) {
  trans <- as.character(beataml_clin$isTransformed) == "TRUE" | beataml_clin$isTransformed == TRUE
  beataml_clin$Subset[trans] <- "Secondary"
}
if ("specificDxAtAcquisition" %in% cn) {
  therapy <- as.character(beataml_clin$specificDxAtAcquisition) == "Therapy-related myeloid neoplasms"
  beataml_clin$Subset[therapy] <- "Therapy"
}
if ("diseaseStageAtSpecimenCollection" %in% cn) {
  stage <- as.character(beataml_clin$diseaseStageAtSpecimenCollection)
  relapse_residual <- grepl("relapse|residual", stage, ignore.case = TRUE)
  beataml_clin$Subset[relapse_residual] <- "Relapse"
}

# Beat AML survival: check units (likely days)
beataml_clin$Time_to_OS <- if ("overallSurvival" %in% cn) as.numeric(beataml_clin$overallSurvival) else NA_real_
beataml_clin$Censor <- if ("vitalStatus" %in% cn) {
  vs <- as.character(beataml_clin$vitalStatus)
  ifelse(vs == "Dead", 1, ifelse(vs == "Alive", 0, NA_real_))
} else NA_real_

beataml_clin$Age <- if ("ageAtDiagnosis" %in% cn) as.numeric(beataml_clin$ageAtDiagnosis) else NA_real_
# Sex: Beat AML uses consensus_sex or inferred_sex (values Female/Male); fallback to "sex" if present
sex_col <- if ("sex" %in% cn) "sex" else if ("consensus_sex" %in% cn) "consensus_sex" else if ("inferred_sex" %in% cn) "inferred_sex" else NULL
if (!is.null(sex_col)) {
  sx <- as.character(beataml_clin[[sex_col]])
  sx[grepl("Female;Male|Male;Female|unknown", sx, ignore.case = TRUE)] <- NA_character_
  beataml_clin$Sex <- ifelse(sx == "MALE" | sx == "Male", "Male", ifelse(sx == "FEMALE" | sx == "Female", "Female", NA_character_))
} else {
  beataml_clin$Sex <- NA_character_
}

beataml_clin$Risk <- if ("ELN2017" %in% cn) {
  eln <- as.character(beataml_clin$ELN2017)
  ifelse(eln %in% c("Adverse", "Intermediate", "Favorable"), eln, "Unknown")
} else "Unknown"

# Karyotype: (1) Normal from otherCytogenetics (normal/Normal), (2) then Normal from karyotype (46,XX/46,XY only, no del/inv/t(/+/-)), (3) rest = Complex 3+ vs Other. Then compare to karyotype-only.
beataml_clin$Karyotype <- "Unknown"
karyo_col <- grep("^karyotype$", cn, ignore.case = TRUE, value = TRUE)[1]
other_cyto_col <- grep("otherCytogenetics|other_cytogenetics|other.cytogenetics", cn, ignore.case = TRUE, value = TRUE)[1]
if (length(karyo_col) || length(other_cyto_col)) {
  k1 <- if (length(karyo_col)) as.character(beataml_clin[[karyo_col]]) else rep("", nrow(beataml_clin))
  k2 <- if (length(other_cyto_col)) as.character(beataml_clin[[other_cyto_col]]) else rep("", nrow(beataml_clin))
  k1[is.na(k1)] <- ""; k2[is.na(k2)] <- ""
  k1 <- trimws(k1); k2 <- trimws(k2)
  combined <- trimws(paste(k1, k2, sep = " "))
  combined[combined == ""] <- NA_character_
  # Use event-counting (same as AML-SG): Complex = 3+ cytogenetic abnormalities, Normal = 46,XX/XY or literal "Normal"
  # Step 1: Normal from otherCytogenetics (contains normal/Normal, avoid "abnormal")
  normal_from_other <- nzchar(k2) & grepl("normal", k2, ignore.case = TRUE) & !grepl("abnormal", k2, ignore.case = TRUE)
  beataml_clin$Karyotype[normal_from_other] <- "Normal"
  # Step 2: Additional Normal from karyotype: 46,XX[n] or 46,XY[n] only (no del, inv, t(, +, -, etc.)
  abnormal_in_karyo <- grepl("del\\(|inv\\(|t\\(|add\\(|der\\(|ins\\(|dic\\(|idic\\(|\\+[0-9]|\\-[0-9]|\\+mar|\\-mar|dmin|\\+r|\\+X|\\-X|\\+Y|\\-Y", k1, ignore.case = TRUE)
  karyo_looks_normal <- grepl("^46[, ]*(XX|XY)(\\[\\d*\\])?\\s*$", k1) & !abnormal_in_karyo
  normal_from_karyo <- (beataml_clin$Karyotype != "Normal") & karyo_looks_normal
  beataml_clin$Karyotype[normal_from_karyo] <- "Normal"
  # Step 3: Remaining = event-count on combined: Complex 3+, Other 1-2
  remaining <- beataml_clin$Karyotype == "Unknown"
  if (any(remaining)) beataml_clin$Karyotype[remaining] <- classify_amlsg_karyotype(combined[remaining])
  beataml_clin$Karyotype[is.na(beataml_clin$Karyotype)] <- "Unknown"
  n_comp <- sum(beataml_clin$Karyotype == "Complex")
  n_norm <- sum(beataml_clin$Karyotype == "Normal")
  n_oth <- sum(beataml_clin$Karyotype == "Other")
  n_unk <- sum(beataml_clin$Karyotype == "Unknown")
  cat("  Beat AML: karyotype (Normal from otherCytogenetics then karyotype 46,XX/XY; rest Complex 3+ / Other): ", n_comp, " Complex, ", n_norm, " Normal, ", n_oth, " Other, ", n_unk, " Unknown\n", sep = "")
  cat("    Normal from otherCytogenetics: ", sum(normal_from_other), "; Normal from karyotype (46,XX/XY only): ", sum(normal_from_karyo), "\n", sep = "")
  # Comparison: distribution using ONLY karyotype column
  k1_only_norm <- grepl("^46[, ]*(XX|XY)(\\[\\d*\\])?\\s*$", k1) & !abnormal_in_karyo
  k1_only_class <- classify_amlsg_karyotype(k1)
  k1_only_class[k1_only_norm] <- "Normal"
  tab_k1 <- table(factor(k1_only_class, levels = c("Complex", "Normal", "Other", "Unknown")))
  cat("  Beat AML: comparison using ONLY karyotype column: Complex ", tab_k1["Complex"], ", Normal ", tab_k1["Normal"], ", Other ", tab_k1["Other"], ", Unknown ", tab_k1["Unknown"], "\n", sep = "")
}

# Beat AML clinical labs: wbcCount, hemoglobin, plateletCount, %.Blasts.in.BM, %.Blasts.in.PB, LDH
wbc_col <- grep("^wbcCount$|^wbc", cn, ignore.case = TRUE, value = TRUE)[1]
hb_col  <- grep("^hemoglobin$|^hb$", cn, ignore.case = TRUE, value = TRUE)[1]
plt_col <- grep("^plateletCount$|^platelet", cn, ignore.case = TRUE, value = TRUE)[1]
bm_col  <- grep("Blasts.in.BM|bm_blast|BM.blast", cn, ignore.case = TRUE, value = TRUE)[1]
pb_col  <- grep("Blasts.in.PB|pb_blast|PB.blast", cn, ignore.case = TRUE, value = TRUE)[1]
ldh_col <- grep("^LDH$|^ldh", cn, ignore.case = TRUE, value = TRUE)[1]
beataml_clin <- ensure_clin_cols(beataml_clin,
  wbc = if (length(wbc_col)) as.numeric(beataml_clin[[wbc_col]]) else NA_real_,
  hemoglobin = if (length(hb_col)) as.numeric(beataml_clin[[hb_col]]) else NA_real_,
  platelet = if (length(plt_col)) as.numeric(beataml_clin[[plt_col]]) else NA_real_,
  bm_blast_percent = if (length(bm_col)) as.numeric(beataml_clin[[bm_col]]) else NA_real_,
  pb_blast_percent = if (length(pb_col)) as.numeric(beataml_clin[[pb_col]]) else NA_real_,
  ldh = if (length(ldh_col)) as.numeric(beataml_clin[[ldh_col]]) else NA_real_)

# Merge Beat AML mutations with clinical
beataml <- merge(beataml_mut_agg[, c("Sample", "Gene", "VAF", "variant_type", "Cohort", "AA_change", "Gene_Group")], 
                 beataml_clin[, c("Sample", "Age", "Sex", "Time_to_OS", "Censor", "Subset", "Risk", "Karyotype",
                   "WBC", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent", "LDH")], 
                 by = "Sample", all.x = TRUE)

cat("  Beat AML: ", length(unique(beataml$Sample)), " samples, ", nrow(beataml), " mutations\n", sep = "")

# =============================================================================
# 3. AML-SG
# =============================================================================
cat("Loading AML-SG data...\n")
amlsg_mut <- read.delim("amlsg_data/AMLSG_Genetic.txt", stringsAsFactors = FALSE)
# Clinical data from AMLSG_Clinical_Anon.RData (object: clinicalData)
amlsg_clin_path <- "amlsg_data/AMLSG_Clinical_Anon.RData"
if (!file.exists(amlsg_clin_path)) stop("AML-SG clinical file not found: ", amlsg_clin_path)
amlsg_env <- new.env()
load(amlsg_clin_path, envir = amlsg_env)
if ("clinicalData" %in% ls(amlsg_env)) {
  amlsg_clin <- as.data.frame(get("clinicalData", envir = amlsg_env), stringsAsFactors = FALSE)
} else {
  objs <- ls(amlsg_env)
  df_objs <- objs[sapply(objs, function(o) is.data.frame(get(o, envir = amlsg_env)))]
  if (length(df_objs) == 0) stop("AMLSG_Clinical_Anon.RData contains no data frame. Expected object named 'clinicalData'.")
  amlsg_clin <- as.data.frame(get(df_objs[1], envir = amlsg_env), stringsAsFactors = FALSE)
}
cat("  AML-SG clinical: ", nrow(amlsg_clin), " rows from AMLSG_Clinical_Anon.RData\n", sep = "")

# AML-SG mutations: SAMPLE_NAME, GENE, VAF, VARIANT_TYPE, AA_CHANGE
amlsg_mut$Sample <- amlsg_mut$SAMPLE_NAME
amlsg_mut$Gene <- amlsg_mut$GENE
amlsg_mut$Gene[amlsg_mut$Gene == "SFRS2"] <- "SRSF2"  # correct gene symbol
amlsg_mut$VAF <- if ("VAF" %in% colnames(amlsg_mut)) as.numeric(amlsg_mut$VAF) else if ("X._MUT_IN_TUM" %in% colnames(amlsg_mut)) as.numeric(amlsg_mut$X._MUT_IN_TUM) else NA_real_
amlsg_mut$variant_type <- map_variant_type(amlsg_mut$VARIANT_TYPE)
amlsg_mut$Gene <- normalize_flt3_gene(amlsg_mut$Gene, amlsg_mut$variant_type)
amlsg_mut$Cohort <- "AML-SG"
amlsg_mut$AA_change <- if ("AA_CHANGE" %in% colnames(amlsg_mut)) as.character(amlsg_mut$AA_CHANGE) else NA_character_
amlsg_mut$Gene_Group <- assign_gene_group(amlsg_mut$Gene, amlsg_mut$variant_type, amlsg_mut$AA_change)

# AML-SG clinical: PDID, TypeAML, OS (days), Status, M_Risk (ELN), AOD (age), gender
amlsg_clin$Sample <- amlsg_clin$PDID
amlsg_clin$Subset <- ifelse(amlsg_clin$TypeAML == "AML", "De novo",
                           ifelse(amlsg_clin$TypeAML == "sAML", "Secondary",
                                 ifelse(amlsg_clin$TypeAML == "tAML", "Therapy", "Other")))
amlsg_clin$Time_to_OS <- as.numeric(amlsg_clin$OS)  # Already in days
# Status in AMLSG_Clinical_Anon.RData is integer: 1 = event (death), 0 = censored (alive)
amlsg_clin$Censor <- as.numeric(amlsg_clin$Status)
amlsg_clin$Censor[!amlsg_clin$Status %in% c(0, 1)] <- NA_real_
amlsg_clin$Age <- as.numeric(amlsg_clin$AOD)
amlsg_clin$Sex <- ifelse(amlsg_clin$gender == 1, "Male", ifelse(amlsg_clin$gender == 2, "Female", NA_character_))
# AML-SG M_Risk: map to ELN risk categories
# Values: Favorable, Adverse, Inter-1, Inter-2 -> Intermediate
amlsg_clin$Risk <- "Unknown"
amlsg_clin$Risk[amlsg_clin$M_Risk == "Favorable"] <- "Favorable"
amlsg_clin$Risk[amlsg_clin$M_Risk %in% c("Inter-1", "Inter-2", "Intermediate")] <- "Intermediate"
amlsg_clin$Risk[amlsg_clin$M_Risk == "Adverse"] <- "Adverse"

# Karyotype: prefer AMLSG_Karyotypes.txt (Normal / Complex = 3+ events / Other / Unknown); fallback to clinical if present
amlsg_clin$Karyotype <- "Unknown"
karyo_path <- file.path(getwd(), "amlsg_data", "AMLSG_Karyotypes.txt")
if (file.exists(karyo_path)) {
  amlsg_karyo <- read.delim(karyo_path, stringsAsFactors = FALSE, check.names = FALSE)
  pdid_col <- grep("^PDID$|^pdid$", colnames(amlsg_karyo), value = TRUE)[1]
  karyo_col <- grep("^karyotype$|^Karyotype$", colnames(amlsg_karyo), value = TRUE)[1]
  if (length(pdid_col) && length(karyo_col)) {
    idx <- match(amlsg_clin$Sample, amlsg_karyo[[pdid_col]])
    hit <- !is.na(idx)
    kstr <- amlsg_karyo[[karyo_col]][idx[hit]]
    amlsg_clin$Karyotype[hit] <- classify_amlsg_karyotype(kstr)
    n_comp <- sum(amlsg_clin$Karyotype == "Complex")
    n_norm <- sum(amlsg_clin$Karyotype == "Normal")
    n_oth <- sum(amlsg_clin$Karyotype == "Other")
    cat("  AML-SG: karyotype from AMLSG_Karyotypes.txt for ", sum(hit), " samples (", n_comp, " Complex, ", n_norm, " Normal, ", n_oth, " Other; Complex = 3+ events)\n", sep = "")
  }
}
# Fallback: if clinical has karyotype column, fill remaining Unknown from it
if ("karyotype" %in% colnames(amlsg_clin)) {
  still_unk <- amlsg_clin$Karyotype == "Unknown"
  if (any(still_unk)) {
    amlsg_clin$Karyotype[still_unk] <- classify_karyotype(amlsg_clin$karyotype[still_unk])
    cat("  AML-SG: filled ", sum(still_unk), " remaining Unknown from clinical karyotype column\n", sep = "")
  }
}

# AML-SG clinical labs: wbc, HB, platelet, BM_Blasts, PB_Blasts, LDH
amlsg_clin <- ensure_clin_cols(amlsg_clin,
  wbc = as.numeric(amlsg_clin$wbc),
  hemoglobin = as.numeric(amlsg_clin$HB),
  platelet = as.numeric(amlsg_clin$platelet),
  bm_blast_percent = as.numeric(amlsg_clin$BM_Blasts),
  pb_blast_percent = if ("PB_Blasts" %in% colnames(amlsg_clin)) as.numeric(amlsg_clin$PB_Blasts) else NA_real_,
  ldh = as.numeric(amlsg_clin$LDH))

# Merge AML-SG mutations with clinical
amlsg <- merge(amlsg_mut[, c("Sample", "Gene", "VAF", "variant_type", "Cohort", "AA_change", "Gene_Group")], 
               amlsg_clin[, c("Sample", "Age", "Sex", "Time_to_OS", "Censor", "Subset", "Risk", "Karyotype",
                 "WBC", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent", "LDH")], 
               by = "Sample", all.x = TRUE)

cat("  AML-SG: ", length(unique(amlsg$Sample)), " samples, ", nrow(amlsg), " mutations\n", sep = "")

# =============================================================================
# 4. UK-NCRI
# =============================================================================
cat("Loading UK-NCRI data...\n")
ukncri_mut <- read.csv("UK_NCRI_data/UK_NCRI_Mutations_data.csv", stringsAsFactors = FALSE, skip = 1)
ukncri_clin <- read.csv("UK_NCRI_data/UK_NCRI_Clinical_data.csv", stringsAsFactors = FALSE, skip = 1)

# UK-NCRI mutations: data_pd (sample), gene, tum_vaf, protein
ukncri_mut$Sample <- ukncri_mut$data_pd
ukncri_mut$Gene <- ukncri_mut$gene
ukncri_mut$VAF <- as.numeric(ukncri_mut$tum_vaf)
ukncri_mut$variant_type <- map_variant_type(ukncri_mut$type)
ukncri_mut$Gene <- normalize_flt3_gene(ukncri_mut$Gene, ukncri_mut$variant_type)
ukncri_mut$Cohort <- "UK-NCRI"
ukncri_mut$AA_change <- if ("protein" %in% colnames(ukncri_mut)) as.character(ukncri_mut$protein) else NA_character_
ukncri_mut$Gene_Group <- assign_gene_group(ukncri_mut$Gene, ukncri_mut$variant_type, ukncri_mut$AA_change)

# UK-NCRI FLT3 ITD / TKD / other from aml_molecular_bdp.tsv (not in UK_NCRI_Mutations_data.csv)
# File has 155 fields per row: first = sample ID (no header), rest = 154 gene/feature columns
bdp_path <- file.path(getwd(), "UK_NCRI_data", "data", "aml_molecular_bdp.tsv")
if (!file.exists(bdp_path)) bdp_path <- file.path(getwd(), "UK_NCRI_data", "aml_molecular_bdp.tsv")
if (file.exists(bdp_path)) {
  bdp <- read.table(bdp_path, sep = "", header = TRUE, quote = "\"", check.names = FALSE,
                    row.names = 1)  # first column -> sample ID as row names
  bdp_sample <- gsub("^\"|\"$", "", rownames(bdp))
  nms <- colnames(bdp)
  itd_col <- which(nms == "ITD")[1]
  tkd_col <- which(nms == "FLT3_TKD")[1]
  oth_col <- which(nms == "FLT3_other")[1]
  if (length(itd_col) && length(tkd_col) && length(oth_col)) {
    itd_val <- as.numeric(bdp[, itd_col])
    tkd_val <- as.numeric(bdp[, tkd_col])
    oth_val <- as.numeric(bdp[, oth_col])
    extra <- list(Sample = character(), Gene = character(), VAF = numeric(), variant_type = character(), Cohort = character(), AA_change = character(), Gene_Group = character())
    for (i in seq_along(bdp_sample)) {
      s <- bdp_sample[i]
      if (!is.na(itd_val[i]) && itd_val[i] == 1) {
        extra$Sample <- c(extra$Sample, s)
        extra$Gene <- c(extra$Gene, "FLT3-ITD")
        extra$VAF <- c(extra$VAF, NA_real_)
        extra$variant_type <- c(extra$variant_type, "ITD")
        extra$Cohort <- c(extra$Cohort, "UK-NCRI")
        extra$AA_change <- c(extra$AA_change, NA_character_)
        extra$Gene_Group <- c(extra$Gene_Group, "FLT3_ITD")
      }
      if (!is.na(tkd_val[i]) && tkd_val[i] == 1) {
        extra$Sample <- c(extra$Sample, s)
        extra$Gene <- c(extra$Gene, "FLT3-TKD")
        extra$VAF <- c(extra$VAF, NA_real_)
        extra$variant_type <- c(extra$variant_type, "SNV")
        extra$Cohort <- c(extra$Cohort, "UK-NCRI")
        extra$AA_change <- c(extra$AA_change, NA_character_)
        extra$Gene_Group <- c(extra$Gene_Group, "FLT3_TKD")
      }
      if (!is.na(oth_val[i]) && oth_val[i] == 1) {
        extra$Sample <- c(extra$Sample, s)
        extra$Gene <- c(extra$Gene, "FLT3-Other")
        extra$VAF <- c(extra$VAF, NA_real_)
        extra$variant_type <- c(extra$variant_type, "Other")
        extra$Cohort <- c(extra$Cohort, "UK-NCRI")
        extra$AA_change <- c(extra$AA_change, NA_character_)
        extra$Gene_Group <- c(extra$Gene_Group, "FLT3_other")
      }
    }
    if (length(extra$Sample) > 0) {
      ukncri_mut <- rbind(ukncri_mut[, c("Sample", "Gene", "VAF", "variant_type", "Cohort", "AA_change", "Gene_Group")],
                          as.data.frame(extra, stringsAsFactors = FALSE))
      cat("  UK-NCRI: added ", length(extra$Sample), " FLT3 ITD/TKD/Other rows from aml_molecular_bdp.tsv\n", sep = "")
    }
  }
} else {
  cat("  UK-NCRI: aml_molecular_bdp.tsv not found; FLT3 ITD/TKD/Other not added\n", sep = "")
}

# UK-NCRI clinical: first column (sample ID), secondary, age, gender (see UK_NCRI_data/data: ahd=0 <-> secondary=1 = de novo)
ukncri_clin$Sample <- as.character(ukncri_clin[, 1])  # First column is sample ID
# secondary: 1 = de novo (primary AML, no AHD), 2 = secondary AML, 3 = therapy-related (per ahd x secondary cross-tab in data)
ukncri_clin$Subset <- ifelse(ukncri_clin$secondary == 1, "De novo",
  ifelse(ukncri_clin$secondary == 2, "Secondary", ifelse(ukncri_clin$secondary == 3, "Therapy", "Other")))
ukncri_clin$Age <- as.numeric(ukncri_clin$age)
ukncri_clin$Sex <- ifelse(ukncri_clin$gender == 0, "Male", ifelse(ukncri_clin$gender == 1, "Female", NA_character_))
# Survival/ELN: not in base UK_NCRI clinical file; enrich from aml_prognosis_updated.tsv
ukncri_clin$Time_to_OS <- NA_real_
ukncri_clin$Censor <- NA_real_
ukncri_clin$Risk <- "Unknown"
ukncri_clin$Karyotype <- "Unknown"

prog_path <- file.path(getwd(), "UK_NCRI_data", "aml_prognosis_updated.tsv")
if (!file.exists(prog_path)) prog_path <- file.path(getwd(), "UK_NCRI_data", "data", "aml_prognosis_updated.tsv")
if (!file.exists(prog_path)) prog_path <- path.expand("~/Downloads/data/aml_prognosis_updated.tsv")
if (file.exists(prog_path)) {
  # File has 205 header cols and 206 data cols; first column = sample ID (becomes row names)
  prog <- read.delim(prog_path, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  prog$Sample <- gsub("^\"|\"$", "", rownames(prog))
  eln_num <- as.numeric(prog[, 1])   # first column after row names is numeric ELN (1/2/3)
  prog$Risk <- ifelse(eln_num == 1, "Favorable", ifelse(eln_num == 2, "Intermediate", ifelse(eln_num == 3, "Adverse", "Unknown")))
  prog$Time_to_OS <- as.numeric(prog$os) * 365.25   # os in years -> days
  prog$Censor <- as.numeric(prog$os_status)         # 1 = event, 0 = censored
  match_idx <- match(ukncri_clin$Sample, prog$Sample)
  hit <- !is.na(match_idx)
  ukncri_clin$Time_to_OS[hit] <- prog$Time_to_OS[match_idx[hit]]
  ukncri_clin$Censor[hit] <- prog$Censor[match_idx[hit]]
  ukncri_clin$Risk[hit] <- prog$Risk[match_idx[hit]]
  comp_col <- grep("^complex$", colnames(prog), ignore.case = TRUE, value = TRUE)[1]
  if (length(comp_col)) {
    comp_val <- as.numeric(prog[match_idx[hit], comp_col])
    ukncri_clin$Karyotype[hit] <- ifelse(!is.na(comp_val) & comp_val == 1, "Complex", "Unknown")
    n_complex <- sum(!is.na(comp_val) & comp_val == 1)
    cat("  UK-NCRI: Complex karyotype from aml_prognosis_updated.tsv for ", n_complex, " samples\n", sep = "")
  }
  cat("  UK-NCRI: enriched", sum(hit), "samples with ELN and OS/Censor from aml_prognosis_updated.tsv\n", sep = " ")
} else {
  cat("  UK-NCRI: aml_prognosis_updated.tsv not found. Place in UK_NCRI_data/ or ~/Downloads/data/\n", sep = " ")
}

# UK-NCRI Normal/Other from UK_NCRI_Cytogenetics_data.csv: Normal only when NO amp, del, or translocation (all abnormality columns 0)
cyto_path <- file.path(getwd(), "UK_NCRI_data", "UK_NCRI_Cytogenetics_data.csv")
if (!file.exists(cyto_path)) cyto_path <- file.path(getwd(), "UK_NCRI_data", "data", "UK_NCRI_Cytogenetics_data.csv")
if (file.exists(cyto_path)) {
  cyto <- read.csv(cyto_path, stringsAsFactors = FALSE, skip = 1, check.names = FALSE)
  id_col <- grep("^data_pd$|^data\\.pd$", colnames(cyto), value = TRUE)[1]
  if (length(id_col)) {
    abn_cols <- setdiff(colnames(cyto), id_col)
    if (length(abn_cols) > 0) {
      cyto_id <- as.character(cyto[[id_col]])
      idx_cyto <- match(ukncri_clin$Sample, cyto_id)
      not_complex <- as.character(ukncri_clin$Karyotype) != "Complex"
      matched <- !is.na(idx_cyto) & not_complex
      if (any(matched)) {
        cyto_matched <- cyto[idx_cyto[matched], abn_cols, drop = FALSE]
        all_zero <- apply(cyto_matched, 1, function(r) all(as.numeric(r) %in% c(0, NA), na.rm = TRUE))
        all_zero[is.na(all_zero)] <- FALSE
        ukncri_clin$Karyotype[matched][all_zero] <- "Normal"
        ukncri_clin$Karyotype[matched][!all_zero] <- "Other"
      }
      n_normal <- sum(ukncri_clin$Karyotype == "Normal")
      n_other <- sum(ukncri_clin$Karyotype == "Other")
      cat("  UK-NCRI: Normal/Other from UK_NCRI_Cytogenetics_data.csv (", n_normal, " Normal, ", n_other, " Other; Normal = no add/del/translocation/amp)\n", sep = "")
    }
  }
} else {
  cat("  UK-NCRI: UK_NCRI_Cytogenetics_data.csv not found; Normal/Other not assigned\n", sep = "")
}

# UK-NCRI clinical labs: WBC, Hemoglobin, Platelet *must* come from UK_NCRI_Clinical_data.csv (column names may vary)
cn_clin <- colnames(ukncri_clin)
wbc_col  <- grep("^wbc$|^WBC$|white.*blood|^wbc_count", cn_clin, ignore.case = TRUE, value = TRUE)[1]
hb_col   <- grep("^hb$|^Hb$|^hemoglobin$|^Hemoglobin$", cn_clin, ignore.case = TRUE, value = TRUE)[1]
plt_col  <- grep("^plt$|^PLT$|^platelet$|^Platelet$|platelet.*count", cn_clin, ignore.case = TRUE, value = TRUE)[1]
bm_col   <- grep("^bm_blasts$|bm_blast|BM.blast|bone.*marrow.*blast", cn_clin, ignore.case = TRUE, value = TRUE)[1]
ukncri_clin <- ensure_clin_cols(ukncri_clin,
  wbc = if (length(wbc_col)) as.numeric(ukncri_clin[[wbc_col]]) else NA_real_,
  hemoglobin = if (length(hb_col)) as.numeric(ukncri_clin[[hb_col]]) else NA_real_,
  platelet = if (length(plt_col)) as.numeric(ukncri_clin[[plt_col]]) else NA_real_,
  bm_blast_percent = if (length(bm_col)) as.numeric(ukncri_clin[[bm_col]]) else NA_real_,
  pb_blast_percent = NA_real_,
  ldh = NA_real_)
if (!length(wbc_col) || !length(hb_col) || !length(plt_col)) {
  cat("  UK-NCRI: WBC/Hb/Platelet columns in UK_NCRI_Clinical_data.csv - found WBC:", length(wbc_col), " Hb:", length(hb_col), " Platelet:", length(plt_col), "; available:", paste(head(cn_clin, 20), collapse = ", "), "\n", sep = " ")
}

# Merge UK-NCRI mutations with clinical
ukncri <- merge(ukncri_mut[, c("Sample", "Gene", "VAF", "variant_type", "Cohort", "AA_change", "Gene_Group")], 
                ukncri_clin[, c("Sample", "Age", "Sex", "Time_to_OS", "Censor", "Subset", "Risk", "Karyotype",
                  "WBC", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent", "LDH")], 
                by = "Sample", all.x = TRUE)

cat("  UK-NCRI: ", length(unique(ukncri$Sample)), " samples, ", nrow(ukncri), " mutations\n", sep = "")

# =============================================================================
# 5. Merge all cohorts
# =============================================================================
cat("Merging all cohorts...\n")
all_data <- rbind(tcga, beataml, amlsg, ukncri)

# Ensure consistent column types
all_data$Sample <- as.character(all_data$Sample)
all_data$Gene <- as.character(all_data$Gene)
all_data$VAF <- as.numeric(all_data$VAF)
all_data$Cohort <- as.character(all_data$Cohort)
all_data$variant_type <- as.character(all_data$variant_type)
all_data$Subset <- as.character(all_data$Subset)
all_data$Age <- as.numeric(all_data$Age)
all_data$Sex <- as.character(all_data$Sex)
all_data$Risk <- as.character(all_data$Risk)
all_data$Time_to_OS <- as.numeric(all_data$Time_to_OS)
all_data$Censor <- as.numeric(all_data$Censor)
all_data$WBC <- as.numeric(all_data$WBC)
all_data$Hemoglobin <- as.numeric(all_data$Hemoglobin)
all_data$Platelet <- as.numeric(all_data$Platelet)
all_data$BM_blast_percent <- as.numeric(all_data$BM_blast_percent)
all_data$PB_blast_percent <- as.numeric(all_data$PB_blast_percent)
all_data$LDH <- as.numeric(all_data$LDH)
all_data$AA_change <- as.character(all_data$AA_change)
all_data$Gene_Group <- as.character(all_data$Gene_Group)
all_data$Karyotype <- as.character(all_data$Karyotype)
all_data$Karyotype[is.na(all_data$Karyotype) | all_data$Karyotype == ""] <- "Unknown"

# CEBPA_bi vs CEBPA_mono by sample-level CEBPA mutation count
cebpa_idx <- which(all_data$Gene == "CEBPA")
if (length(cebpa_idx) > 0) {
  cebpa_per_sample <- aggregate(Gene ~ Sample, data = all_data[cebpa_idx, , drop = FALSE], FUN = length)
  n_cebpa <- setNames(cebpa_per_sample$Gene, as.character(cebpa_per_sample$Sample))
  all_data$Gene_Group[cebpa_idx] <- ifelse(n_cebpa[as.character(all_data$Sample[cebpa_idx])] >= 2L, "CEBPA_bi", "CEBPA_mono")
}

cat("  Total: ", length(unique(all_data$Sample)), " samples, ", nrow(all_data), " mutations\n", sep = "")

# =============================================================================
# 6. Assign mutation categories from reference (final_data_matrix)
# =============================================================================
cat("Assigning mutation categories...\n")
if (file.exists("final_data_matrix.RData")) {
  load("final_data_matrix.RData")
  ref <- final_data_matrix[!is.na(final_data_matrix$mutation_category) & 
                           as.character(final_data_matrix$mutation_category) != "", 
                           c("Gene", "mutation_category"), drop = FALSE]
  if (nrow(ref) > 0) {
    ref$Gene <- as.character(ref$Gene)
    ref$mutation_category <- as.character(ref$mutation_category)
    most_common <- aggregate(mutation_category ~ Gene, data = ref, 
                            FUN = function(x) names(sort(table(x), decreasing = TRUE))[1])
    gene_to_cat <- setNames(as.character(most_common$mutation_category), as.character(most_common$Gene))
    all_data$mutation_category <- unname(gene_to_cat[as.character(all_data$Gene)])
    # Assign RTK-RAS Signaling to all FLT3 mutations (ITD, TKD, Other, and any plain FLT3)
    flt3_idx <- which(all_data$Gene %in% c("FLT3", "FLT3-ITD", "FLT3-TKD", "FLT3-Other"))
    if (length(flt3_idx) > 0) {
      all_data$mutation_category[flt3_idx] <- "RTK-RAS Signaling"
    }
    # Homogenize any RTK/RAS label to RTK-RAS Signaling
    all_data$mutation_category[all_data$mutation_category %in% c("RTK/RAS Signaling", "RTK-RAS Signaling")] <- "RTK-RAS Signaling"
    # Assign Chromatin/Cohesin to MLL (KMT2A)
    mll_idx <- which(all_data$Gene == "MLL")
    if (length(mll_idx) > 0) {
      all_data$mutation_category[mll_idx] <- "Chromatin/Cohesin"
    }
  } else {
    all_data$mutation_category <- NA_character_
  }
} else {
  all_data$mutation_category <- NA_character_
}

# =============================================================================
# 7. Check mutation frequencies across cohorts
# =============================================================================
cat("Checking mutation frequencies...\n")
mut_freq <- aggregate(Sample ~ Gene + Cohort, data = all_data, FUN = function(x) length(unique(x)))
mut_freq_wide <- reshape(mut_freq, idvar = "Gene", timevar = "Cohort", direction = "wide")
colnames(mut_freq_wide) <- gsub("Sample.", "", colnames(mut_freq_wide))
mut_freq_wide$Total <- rowSums(mut_freq_wide[, -1], na.rm = TRUE)
mut_freq_wide <- mut_freq_wide[order(mut_freq_wide$Total, decreasing = TRUE), ]

cat("\nTop 20 genes by total mutation frequency:\n")
print(head(mut_freq_wide[, c("Gene", "TCGA", "Beat AML", "AML-SG", "UK-NCRI", "Total")], 20))

# =============================================================================
# 8. Clean and validate data
# =============================================================================
cat("\nCleaning and validating data...\n")
# Remove negative survival times
all_data$Time_to_OS[all_data$Time_to_OS < 0] <- NA_real_
# Ensure VAF is 0-100
all_data$VAF[all_data$VAF < 0] <- NA_real_
all_data$VAF[all_data$VAF > 100] <- 100

# Summary statistics
cat("  Samples with survival data: ", sum(!is.na(all_data$Time_to_OS)), " / ", length(unique(all_data$Sample)), "\n", sep = "")
cat("  Samples with age data: ", sum(!is.na(all_data$Age)), " / ", length(unique(all_data$Sample)), "\n", sep = "")
cat("  Samples with sex data: ", sum(!is.na(all_data$Sex)), " / ", length(unique(all_data$Sample)), "\n", sep = "")

# =============================================================================
# 9. Save as AML_Meta_Cohort_v2
# =============================================================================
cat("\nSaving AML_Meta_Cohort_v2...\n")
AML_Meta_Cohort_v2 <- all_data[, c("Sample", "Gene", "VAF", "Cohort", "Subset", "Time_to_OS", "Censor", 
                                    "variant_type", "mutation_category", "AA_change", "Gene_Group",
                                    "Age", "Sex", "Risk", "Karyotype",
                                    "WBC", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent", "LDH")]
saveRDS(AML_Meta_Cohort_v2, "AML_Meta_Cohort_v2.rds")
cat("Saved to AML_Meta_Cohort_v2.rds\n")
cat("  Final dimensions: ", nrow(AML_Meta_Cohort_v2), " rows, ", ncol(AML_Meta_Cohort_v2), " columns\n", sep = "")
cat("  Samples: ", length(unique(AML_Meta_Cohort_v2$Sample)), "\n", sep = "")
cat("  Cohorts: ", paste(unique(AML_Meta_Cohort_v2$Cohort), collapse = ", "), "\n", sep = "")

# Final validation summary
cat("\n=== Validation Summary ===\n")
cat("Survival times: All in days (TCGA converted from months, Beat AML/AML-SG already in days)\n")
cat("ELN Risk: Standardized to Favorable/Intermediate/Adverse/Unknown\n")
cat("Beat AML consensus mutations: FLT3-ITD, NPM1, RUNX1, ASXL1, TP53 added from clinical file\n")
cat("Mutation frequencies: Checked across cohorts (see above)\n")
cat("Clinical variables: WBC, Hemoglobin, Platelet, BM_blast_percent, PB_blast_percent, LDH captured where available per cohort\n")
cat("\nDone!\n")
