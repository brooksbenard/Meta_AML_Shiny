# =============================================================================
# Build a single cleaned, MAF-style mutation table from all four cohorts.
# Sources: TCGA LAML, Beat AML2, AML-SG (Genetic + FLT3ITD), UK-NCRI.
# Run from Meta_AML_Shiny: Rscript build_aggregated_maf.R
# Output: Meta_AML_aggregated_mutations.maf.tsv (for review / future app download)
# =============================================================================

setwd(if (exists("ofile")) dirname(ofile) else getwd())
if (!dir.exists("laml_tcga") || !dir.exists("beataml2_data") || !dir.exists("amlsg_data") || !dir.exists("UK_NCRI_data")) {
  stop("Run from Meta_AML_Shiny root (laml_tcga, beataml2_data, amlsg_data, UK_NCRI_data must exist).")
}

# MAF-like output columns (standard where possible, plus Cohort/VAF/Variant_Caller)
MAF_COLS <- c(
  "Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Start_Position", "End_Position",
  "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
  "Tumor_Sample_Barcode", "HGVSp_Short", "AA_change", "VAF", "Cohort",
  "Variant_Caller", "cds_change", "Consequence", "source_file"
)

# Helper: normalize variant_classification -> Variant_Type (app logic)
variant_classification_to_type <- function(vc) {
  vc <- as.character(vc)
  out <- rep("Other", length(vc))
  out[vc %in% c("Sub", "SNV", "SNP", "missense_variant", "stop_gained", "Missense_Mutation", "Nonsense_Mutation", "missense", "nonsense")] <- "SNV"
  out[vc %in% c("Del", "DEL", "Frame_Shift_Del", "frame_shift_del", "Deletion")] <- "Deletion"
  out[vc %in% c("splice_acceptor_variant", "splice_donor_variant", "Splice_Site", "splice_site", "Splicing")] <- "Splicing"
  out[vc %in% c("Indel", "In_Frame_Del", "In_Frame_Ins", "in_frame_del", "in_frame_ins", "Frame_Shift_Ins", "frame_shift_ins", "Insertion", "frameshift_variant")] <- "Indel"
  out[vc %in% c("ITD", "Internal_Tandem_Duplication")] <- "ITD"
  out[vc %in% c("PTD", "Partial_Tandem_Duplication")] <- "PTD"
  out
}

# NPM1 frameshift in AA_change -> Indel
apply_npm1_frameshift <- function(df) {
  if (!"Hugo_Symbol" %in% colnames(df) && "Gene" %in% colnames(df)) df$Hugo_Symbol <- df$Gene
  if (!"AA_change" %in% colnames(df)) return(df)
  aa <- as.character(df$AA_change)
  aa[is.na(aa)] <- ""
  npm1 <- as.character(df$Hugo_Symbol) == "NPM1"
  frameshift <- grepl("frame_shift|frameshift|frame shift|fs\\*|fs\\d|fs$|fs[\\s,;]|ins.*fs|fs.*ins", aa, ignore.case = TRUE)
  if ("Variant_Type" %in% colnames(df)) df$Variant_Type[npm1 & frameshift] <- "Indel"
  df
}

# ---------- 1. TCGA LAML ----------
tcga_path <- "laml_tcga/data_mutations.txt"
tcga <- NULL
if (file.exists(tcga_path)) {
  tcga <- read.delim(tcga_path, stringsAsFactors = FALSE, check.names = FALSE)
  tcga$Cohort <- "TCGA"
  tcga$Tumor_Sample_Barcode <- as.character(tcga$Tumor_Sample_Barcode)
  tcga$Hugo_Symbol <- as.character(tcga$Hugo_Symbol)
  tcga$Variant_Caller <- if ("Sequencer" %in% colnames(tcga)) as.character(tcga$Sequencer) else NA_character_
  if ("HGVSp_Short" %in% colnames(tcga)) tcga$AA_change <- tcga$HGVSp_Short else tcga$AA_change <- NA_character_
  if ("t_alt_count" %in% colnames(tcga) && "t_ref_count" %in% colnames(tcga)) {
    tot <- as.numeric(tcga$t_alt_count) + as.numeric(tcga$t_ref_count)
    tcga$VAF <- ifelse(!is.na(tot) & tot > 0, 100 * as.numeric(tcga$t_alt_count) / tot, NA_real_)
  } else tcga$VAF <- NA_real_
  tcga$Variant_Type <- variant_classification_to_type(tcga$Variant_Classification)
  if ("Consequence" %in% colnames(tcga)) {
    vc2 <- variant_classification_to_type(tcga$Consequence)
    tcga$Variant_Type[vc2 != "Other"] <- vc2[vc2 != "Other"]
  }
  tcga$cds_change <- if ("HGVSc" %in% colnames(tcga)) tcga$HGVSc else NA_character_
  tcga$Consequence <- if ("Consequence" %in% colnames(tcga)) tcga$Consequence else tcga$Variant_Classification
  tcga$source_file <- "laml_tcga/data_mutations.txt"
  tcga <- apply_npm1_frameshift(tcga)
}

# ---------- 2. Beat AML2 ----------
beat_path <- "beataml2_data/mutations.txt"
beat <- NULL
if (file.exists(beat_path)) {
  beat <- read.delim(beat_path, stringsAsFactors = FALSE, check.names = FALSE)
  beat$Cohort <- "Beat AML"
  beat$Tumor_Sample_Barcode <- as.character(beat$dbgap_sample_id)
  beat$Hugo_Symbol <- as.character(beat$symbol)
  beat$Variant_Caller <- if ("genotyper" %in% colnames(beat)) as.character(beat$genotyper) else NA_character_
  vc <- tolower(as.character(beat$variant_classification))
  beat$Hugo_Symbol[beat$symbol == "FLT3" & vc == "inframe_insertion"] <- "FLT3-ITD"
  beat$Hugo_Symbol[beat$symbol == "FLT3" & vc %in% c("missense_variant", "inframe_deletion")] <- "FLT3-TKD"
  beat$Hugo_Symbol[beat$symbol == "FLT3" & beat$Hugo_Symbol == "FLT3"] <- "FLT3-TKD"
  beat$Chromosome <- if ("seqnames" %in% colnames(beat)) as.character(beat$seqnames) else NA_character_
  beat$Start_Position <- if ("pos_start" %in% colnames(beat)) as.integer(beat$pos_start) else NA_integer_
  beat$End_Position <- if ("pos_end" %in% colnames(beat)) as.integer(beat$pos_end) else NA_integer_
  beat$Reference_Allele <- if ("ref" %in% colnames(beat)) as.character(beat$ref) else NA_character_
  beat$Tumor_Seq_Allele1 <- if ("alt" %in% colnames(beat)) as.character(beat$alt) else NA_character_
  beat$Tumor_Seq_Allele2 <- NA_character_
  beat$HGVSp_Short <- if ("hgvsp_short" %in% colnames(beat)) as.character(beat$hgvsp_short) else NA_character_
  beat$AA_change <- beat$HGVSp_Short
  beat$VAF <- as.numeric(beat$t_vaf) * 100
  beat$Variant_Classification <- as.character(beat$variant_classification)
  beat$Variant_Type <- variant_classification_to_type(beat$Variant_Classification)
  beat$cds_change <- if ("hgvsc" %in% colnames(beat)) beat$hgvsc else NA_character_
  beat$Consequence <- beat$Variant_Classification
  beat$Entrez_Gene_Id <- NA_integer_
  beat$source_file <- "beataml2_data/mutations.txt"
  beat <- apply_npm1_frameshift(beat)
}

# ---------- 3. AML-SG Genetic ----------
amlsg_path <- "amlsg_data/AMLSG_Genetic.txt"
amlsg <- NULL
if (file.exists(amlsg_path)) {
  amlsg <- read.delim(amlsg_path, stringsAsFactors = FALSE, check.names = FALSE)
  amlsg$Cohort <- "AML-SG"
  amlsg$Tumor_Sample_Barcode <- as.character(amlsg$SAMPLE_NAME)
  amlsg$Hugo_Symbol <- as.character(amlsg$GENE)
  amlsg$Variant_Caller <- NA_character_
  amlsg$Chromosome <- if ("CHR" %in% colnames(amlsg)) as.character(amlsg$CHR) else NA_character_
  amlsg$Start_Position <- if ("POSITION" %in% colnames(amlsg)) suppressWarnings(as.integer(amlsg$POSITION)) else NA_integer_
  amlsg$End_Position <- if ("MAX_POSITION" %in% colnames(amlsg)) suppressWarnings(as.integer(amlsg$MAX_POSITION)) else amlsg$Start_Position
  amlsg$Reference_Allele <- if ("WT" %in% colnames(amlsg)) as.character(amlsg$WT) else NA_character_
  amlsg$Tumor_Seq_Allele1 <- if ("MT" %in% colnames(amlsg)) as.character(amlsg$MT) else NA_character_
  amlsg$Tumor_Seq_Allele2 <- NA_character_
  amlsg$AA_change <- if ("AA_CHANGE" %in% colnames(amlsg)) as.character(amlsg$AA_CHANGE) else NA_character_
  amlsg$HGVSp_Short <- amlsg$AA_change
  if ("VAF" %in% colnames(amlsg)) amlsg$VAF <- as.numeric(amlsg$VAF) else if ("%_MUT_IN_TUM" %in% colnames(amlsg)) amlsg$VAF <- as.numeric(amlsg$`%_MUT_IN_TUM`) else amlsg$VAF <- NA_real_
  amlsg$Variant_Classification <- if ("VARIANT_TYPE" %in% colnames(amlsg)) as.character(amlsg$VARIANT_TYPE) else "Unknown"
  amlsg$Variant_Type <- variant_classification_to_type(amlsg$Variant_Classification)
  amlsg$cds_change <- if ("CDNA_CHANGE" %in% colnames(amlsg)) as.character(amlsg$CDNA_CHANGE) else NA_character_
  amlsg$Consequence <- if ("CONSEQUENCE" %in% colnames(amlsg)) as.character(amlsg$CONSEQUENCE) else amlsg$Variant_Classification
  amlsg$Entrez_Gene_Id <- NA_integer_
  amlsg$source_file <- "amlsg_data/AMLSG_Genetic.txt"
  amlsg <- apply_npm1_frameshift(amlsg)
}

# ---------- 3b. AML-SG FLT3 ITD (add rows for samples with ITD not in Genetic) ----------
flt3itd_path <- "amlsg_data/AMLSG_FLT3ITD.txt"
if (file.exists(flt3itd_path) && !is.null(amlsg) && nrow(amlsg) > 0) {
  flt3itd <- read.delim(flt3itd_path, stringsAsFactors = FALSE)
  flt3itd <- flt3itd[tolower(as.character(flt3itd$FLT3_ITD_status)) == "itd", , drop = FALSE]
  if (nrow(flt3itd) > 0) {
    sid <- if ("ID" %in% colnames(flt3itd)) as.character(flt3itd$ID) else as.character(flt3itd[[grep("Sample|sample", colnames(flt3itd), value = TRUE)[1]]])
    in_genetic <- unique(amlsg$Tumor_Sample_Barcode[amlsg$Hugo_Symbol == "FLT3-ITD"])
    need_rows <- !sid %in% in_genetic
    if (any(need_rows)) {
      add <- amlsg[1, , drop = FALSE][rep(1, sum(need_rows)), , drop = FALSE]
      add[] <- NA
      add$Hugo_Symbol <- "FLT3-ITD"
      add$Chromosome <- "13"
      add$Variant_Classification <- "Internal_Tandem_Duplication"
      add$Variant_Type <- "ITD"
      add$Tumor_Seq_Allele1 <- "ITD"
      add$Tumor_Sample_Barcode <- sid[need_rows]
      add$Cohort <- "AML-SG"
      add$Variant_Caller <- NA_character_
      add$Consequence <- "ITD"
      add$source_file <- "amlsg_data/AMLSG_FLT3ITD.txt"
      if ("Read_count" %in% colnames(flt3itd) && "Coverage" %in% colnames(flt3itd))
        add$VAF <- 100 * as.numeric(flt3itd$Read_count[need_rows]) / pmax(as.numeric(flt3itd$Coverage[need_rows]), 1)
      amlsg <- rbind(amlsg, add)
    }
  }
}

# ---------- 4. UK-NCRI ----------
ncri_path <- "UK_NCRI_data/UK_NCRI_Mutations_data.csv"
ncri <- NULL
if (file.exists(ncri_path)) {
  raw <- readLines(ncri_path, n = 3)
  skip <- if (any(grepl("^S\\.Data|^#", raw))) 1 else 0
  ncri <- read.csv(ncri_path, skip = skip, stringsAsFactors = FALSE, check.names = FALSE)
  ncri$Cohort <- "UK-NCRI"
  sample_col <- if ("data_pd" %in% colnames(ncri)) "data_pd" else colnames(ncri)[1]
  ncri$Tumor_Sample_Barcode <- as.character(ncri[[sample_col]])
  ncri$Variant_Caller <- NA_character_
  ncri$Hugo_Symbol <- if ("gene" %in% colnames(ncri)) as.character(ncri$gene) else NA_character_
  ncri$Chromosome <- if ("chr" %in% colnames(ncri)) as.character(ncri$chr) else NA_character_
  ncri$Start_Position <- if ("start" %in% colnames(ncri)) suppressWarnings(as.integer(ncri$start)) else NA_integer_
  ncri$End_Position <- if ("end" %in% colnames(ncri)) suppressWarnings(as.integer(ncri$end)) else ncri$Start_Position
  ncri$Reference_Allele <- if ("wt" %in% colnames(ncri)) as.character(ncri$wt) else NA_character_
  ncri$Tumor_Seq_Allele1 <- if ("mt" %in% colnames(ncri)) as.character(ncri$mt) else NA_character_
  ncri$Tumor_Seq_Allele2 <- NA_character_
  ncri$AA_change <- if ("protein" %in% colnames(ncri)) as.character(ncri$protein) else NA_character_
  ncri$HGVSp_Short <- ncri$AA_change
  ncri$VAF <- if ("tum_vaf" %in% colnames(ncri)) as.numeric(ncri$tum_vaf) else NA_real_
  type_col <- if ("type" %in% colnames(ncri)) as.character(ncri$type) else "Unknown"
  ncri$Variant_Classification <- type_col
  ncri$Variant_Type <- variant_classification_to_type(type_col)
  ncri$cds_change <- if ("cds" %in% colnames(ncri)) as.character(ncri$cds) else NA_character_
  ncri$Consequence <- type_col
  ncri$Entrez_Gene_Id <- NA_integer_
  ncri$source_file <- "UK_NCRI_data/UK_NCRI_Mutations_data.csv"
  ncri <- apply_npm1_frameshift(ncri)
}

# ---------- Bind and select MAF columns ----------
out_list <- list()
if (!is.null(tcga)) out_list$tcga <- tcga
if (!is.null(beat)) out_list$beat <- beat
if (!is.null(amlsg)) out_list$amlsg <- amlsg
if (!is.null(ncri)) out_list$ncri <- ncri

if (length(out_list) == 0) stop("No mutation files found.")

# Ensure every table has the same columns (fill missing with NA) so rbind works
all_cols <- unique(c(MAF_COLS, unlist(lapply(out_list, colnames))))
for (i in seq_along(out_list)) {
  for (c in all_cols) if (!c %in% colnames(out_list[[i]])) out_list[[i]][[c]] <- NA
  out_list[[i]] <- out_list[[i]][, all_cols, drop = FALSE]
}
all_mut <- do.call(rbind, out_list)
rownames(all_mut) <- NULL

# Keep and order MAF columns only
present <- intersect(MAF_COLS, colnames(all_mut))
all_mut <- all_mut[, present, drop = FALSE]

out_file <- "Meta_AML_aggregated_mutations.maf.tsv"
write.table(all_mut, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Written ", nrow(all_mut), " rows to ", out_file)
message("Cohort counts: ")
print(table(all_mut$Cohort, useNA = "ifany"))
