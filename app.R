# =============================================================================
# Meta AML Explorer - Shiny App for Exploratory Analysis of AML Mutational Data
# Uses base R only to avoid data.table/dplyr "Argument 'x' is not a vector: list" errors
# =============================================================================

library(shiny)
library(ggplot2)
library(survival)
library(DT)
library(scales)

has_survminer <- requireNamespace("survminer", quietly = TRUE)
if (has_survminer) library(survminer)
has_gridExtra <- requireNamespace("gridExtra", quietly = TRUE)
if (has_gridExtra) library(gridExtra)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)
has_ComplexHeatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)
if (has_ComplexHeatmap) library(ComplexHeatmap)
has_maxstat <- requireNamespace("maxstat", quietly = TRUE)
if (has_maxstat) library(maxstat)
has_bestglm <- requireNamespace("bestglm", quietly = TRUE)
if (has_bestglm) library(bestglm)
has_g3viz <- requireNamespace("g3viz", quietly = TRUE)
if (has_g3viz) library(g3viz)
has_maftools <- requireNamespace("maftools", quietly = TRUE)
if (has_maftools) library(maftools)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggrepel) library(ggrepel)
has_rentrez <- requireNamespace("rentrez", quietly = TRUE)
if (has_rentrez) library(rentrez)
has_httr <- requireNamespace("httr", quietly = TRUE)
if (has_httr) library(httr)
has_UpSetR <- requireNamespace("UpSetR", quietly = TRUE)
if (has_UpSetR) library(UpSetR)

# Empty stub for legacy references (sidebar Cohort dropdown shows "All" only; Meta AML4 uses its own data)
final_data_matrix <- data.frame(Sample = character(), Gene = character(), Gene_for_analysis = character(), Cohort = character(), Subset = character(), mutation_category = character(), stringsAsFactors = FALSE)

# Normalize AML_Meta_Cohort column names for the app
normalize_AML_Meta_Cohort <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  # Map source columns -> app column names (select + rename)
  # Sample: map from various column names across four datasets
  # NCRI: data_pd, SAMPLE_NAME
  # AMLSG: SAMPLE_NAME
  # Beat AML: dbgap_sample_id
  # TCGA: Tumor_Sample_Barcode
  if ("SAMPLE_NAME" %in% colnames(d)) d$Sample <- as.character(d$SAMPLE_NAME)
  else if ("data_pd" %in% colnames(d)) d$Sample <- as.character(d$data_pd)
  else if ("dbgap_sample_id" %in% colnames(d)) d$Sample <- as.character(d$dbgap_sample_id)
  else if ("Tumor_Sample_Barcode" %in% colnames(d)) d$Sample <- as.character(d$Tumor_Sample_Barcode)
  else if ("patient_id" %in% colnames(d)) d$Sample <- as.character(d$patient_id)
  else if ("sample_id" %in% colnames(d)) d$Sample <- as.character(d$sample_id)
  # Gene, VAF, Cohort
  if (!"Sample" %in% colnames(d)) d$Sample <- paste0("S", seq_len(nrow(d)))
  if ("Gene" %in% colnames(d)) {
    d$Gene <- as.character(d$Gene)
    # Normalize gene names: extract base gene (remove variant suffixes like _ITD, _p.R140Q, etc.)
    # Keep original in Gene_for_analysis if needed, but Gene should be base name for grouping
    gene_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", d$Gene)
    has_suffix <- grepl("^([A-Z0-9]+)[_-]", d$Gene)
    d$Gene[has_suffix] <- gene_base[has_suffix]
  }
  if ("GENE" %in% colnames(d) && !"Gene" %in% colnames(d)) {
    d$Gene <- as.character(d$GENE)
    gene_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", d$Gene)
    has_suffix <- grepl("^([A-Z0-9]+)[_-]", d$Gene)
    d$Gene[has_suffix] <- gene_base[has_suffix]
  }
  if ("gene" %in% colnames(d) && !"Gene" %in% colnames(d)) {
    d$Gene <- as.character(d$gene)
    gene_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", d$Gene)
    has_suffix <- grepl("^([A-Z0-9]+)[_-]", d$Gene)
    d$Gene[has_suffix] <- gene_base[has_suffix]
  }
  if ("Hugo_Symbol" %in% colnames(d) && !"Gene" %in% colnames(d)) {
    d$Gene <- as.character(d$Hugo_Symbol)
    gene_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", d$Gene)
    has_suffix <- grepl("^([A-Z0-9]+)[_-]", d$Gene)
    d$Gene[has_suffix] <- gene_base[has_suffix]
  }
  if ("symbol" %in% colnames(d) && !"Gene" %in% colnames(d)) {
    d$Gene <- as.character(d$symbol)
    gene_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", d$Gene)
    has_suffix <- grepl("^([A-Z0-9]+)[_-]", d$Gene)
    d$Gene[has_suffix] <- gene_base[has_suffix]
  }
  if ("VAF" %in% colnames(d)) d$VAF <- as.numeric(d$VAF)
  if ("tum_vaf" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$tum_vaf)
  if ("TUM_VAF" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$TUM_VAF)
  if ("t_vaf" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$t_vaf)
  if ("%_MUT_IN_TUM" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$"%_MUT_IN_TUM")
  # TCGA VAF columns
  if ("Tumor_VAF" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$Tumor_VAF) * 100
  if ("tumor_vaf" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$tumor_vaf) * 100
  if ("TumorVAF" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$TumorVAF) * 100
  if ("Tumor_Alt_Allele_Freq" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$Tumor_Alt_Allele_Freq) * 100
  if ("tumor_alt_allele_freq" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$tumor_alt_allele_freq) * 100
  # TCGA GDC MAF: no precomputed VAF; derive from t_alt_count and t_ref_count (VAF = alt/(ref+alt) as %)
  if (!"VAF" %in% colnames(d) && "t_alt_count" %in% colnames(d) && "t_ref_count" %in% colnames(d)) {
    t_alt <- as.numeric(d$t_alt_count)
    t_ref <- as.numeric(d$t_ref_count)
    tot <- t_alt + t_ref
    d$VAF <- ifelse(!is.na(tot) & tot > 0, 100 * t_alt / tot, NA_real_)
  }
  # TCGA WashU / other: TumorVAF_WU often 0-1
  if ("TumorVAF_WU" %in% colnames(d) && !"VAF" %in% colnames(d)) d$VAF <- as.numeric(d$TumorVAF_WU) * 100
  if ("Cohort" %in% colnames(d)) {
    d$Cohort <- as.character(d$Cohort)
    d$Cohort[d$Cohort == "NCRI"] <- "UK-NCRI"
    d$Cohort[d$Cohort == "AMLSG"] <- "AML-SG"
    d$Cohort[d$Cohort == "Beat_AML"] <- "Beat AML"
  }
  # Subset: preserve if already present; normalize to capitalized display (De novo, Secondary, Relapse, Therapy, Other)
  if ("Subset" %in% colnames(d)) {
    d$Subset <- as.character(d$Subset)
    d$Subset[is.na(d$Subset) | d$Subset == ""] <- "Other"
    d$Subset[d$Subset %in% c("de_novo", "De_novo")] <- "De novo"
    d$Subset[d$Subset == "secondary"] <- "Secondary"
    d$Subset[d$Subset == "relapse"] <- "Relapse"
    d$Subset[d$Subset %in% c("therapy", "tAML")] <- "Therapy"
    d$Subset[d$Subset %in% c("other", "oAML", "unknown", "")] <- "Other"
  } else if ("type_AML" %in% colnames(d)) {
    t <- as.character(d$type_AML)
    d$Subset <- ifelse(t %in% c("De_novo", "de_novo"), "De novo",
      ifelse(t %in% c("Secondary", "secondary"), "Secondary",
      ifelse(t %in% c("Relapse", "relapse"), "Relapse",
      ifelse(t %in% c("tAML", "therapy"), "Therapy",
      ifelse(t %in% c("oAML", "other", "unknown", ""), "Other", "Other")))))
  } else if (all(c("De_novo", "Secondary") %in% colnames(d))) {
    d$Subset <- "Other"
    de_ok <- which(d$De_novo %in% c(TRUE, "TRUE", "True", 1))
    sec_ok <- which(d$Secondary %in% c(TRUE, "TRUE", "True", 1))
    rel_ok <- if ("Relapse" %in% colnames(d)) which(d$Relapse %in% c(TRUE, "TRUE", "True", 1)) else integer(0)
    if (length(de_ok)) d$Subset[de_ok] <- "De novo"
    if (length(sec_ok)) d$Subset[sec_ok] <- "Secondary"
    if (length(rel_ok)) d$Subset[rel_ok] <- "Relapse"
  } else d$Subset <- "Other"
  # Survival: os -> Time_to_OS, os_status -> Censor (1 = event, 0 = censored)
  if ("os" %in% colnames(d)) { d$Time_to_OS <- as.numeric(d$os); d$Time_to_OS[d$Time_to_OS < 0] <- NA }
  if ("os_status" %in% colnames(d)) d$Censor <- as.numeric(d$os_status)
  # variant_classification -> variant_type (Meta AML4 scheme: SNV, Deletion, Splicing, Indel, ITD, PTD, Other)
  # Handle VARIANT_TYPE (AMLSG), variant_classification (Beat AML, TCGA), Consequence (TCGA), type (NCRI)
  vc <- NULL
  if ("VARIANT_TYPE" %in% colnames(d)) vc <- as.character(d$VARIANT_TYPE)
  else if ("variant_classification" %in% colnames(d)) vc <- as.character(d$variant_classification)
  else if ("Consequence" %in% colnames(d)) vc <- as.character(d$Consequence)
  else if ("type" %in% colnames(d)) vc <- as.character(d$type)
  else if ("variant_type" %in% colnames(d)) {
    d$variant_type <- as.character(d$variant_type)
    d$variant_type[is.na(d$variant_type) | d$variant_type == ""] <- "Other"
    vc <- NULL
  }
  if (!is.null(vc)) {
    d$variant_type <- "Other"
    d$variant_type[vc %in% c("Sub", "SNV", "SNP", "missense_variant", "stop_gained", "Missense_Mutation", "Nonsense_Mutation", "missense", "nonsense")] <- "SNV"
    d$variant_type[vc %in% c("Del", "DEL", "Frame_Shift_Del", "frame_shift_del", "Deletion")] <- "Deletion"
    d$variant_type[vc %in% c("splice_acceptor_variant", "splice_donor_variant", "Splice_Site", "splice_site", "Splicing")] <- "Splicing"
    d$variant_type[vc %in% c("Indel", "In_Frame_Del", "In_Frame_Ins", "in_frame_del", "in_frame_ins")] <- "Indel"
    d$variant_type[vc %in% c("ITD", "Internal_Tandem_Duplication")] <- "ITD"
    d$variant_type[vc %in% c("PTD", "Partial_Tandem_Duplication")] <- "PTD"
    d$variant_type[is.na(d$variant_type) | d$variant_type == ""] <- "Other"
  } else if (!"variant_type" %in% colnames(d)) {
    d$variant_type <- "Other"
  }
  # Protein mutation columns for lollipop: unify from four datasets
  # NCRI: protein, AA_CHANGE
  # AMLSG: AA_CHANGE
  # Beat AML: hgvsp_short, hgvsp, protein_position, amino_acids
  # TCGA: HGVSp_Short, HGVSp, Protein_position
  if ("AA_CHANGE" %in% colnames(d) && !"AA_change" %in% colnames(d)) d$AA_change <- as.character(d$AA_CHANGE)
  if ("aa_change" %in% colnames(d) && !"AA_change" %in% colnames(d)) d$AA_change <- as.character(d$aa_change)
  if ("hgvsp_short" %in% colnames(d) && !"HGVSp_Short" %in% colnames(d)) d$HGVSp_Short <- as.character(d$hgvsp_short)
  if ("HGVSp" %in% colnames(d) && !"HGVSp_Short" %in% colnames(d)) d$HGVSp_Short <- as.character(d$HGVSp)
  if ("hgvsp" %in% colnames(d) && !"HGVSp_Short" %in% colnames(d)) d$HGVSp_Short <- as.character(d$hgvsp)
  if ("protein_change" %in% colnames(d) && !"Protein_Change" %in% colnames(d)) d$Protein_Change <- as.character(d$protein_change)
  if ("protein" %in% colnames(d)) {
    d$protein <- as.character(d$protein)
    if (!"AA_change" %in% colnames(d)) d$AA_change <- d$protein
    if (!"HGVSp_Short" %in% colnames(d)) {
      has_hgvsp <- !is.na(d$protein) & d$protein != "" & grepl("^p\\.", d$protein, ignore.case = TRUE)
      if (any(has_hgvsp)) {
        if (!"HGVSp_Short" %in% colnames(d)) d$HGVSp_Short <- rep(NA_character_, nrow(d))
        d$HGVSp_Short[has_hgvsp] <- d$protein[has_hgvsp]
      }
    }
  }
  # For Beat AML: if we have protein_position and amino_acids, construct HGVSp_Short
  # Format: amino_acids = "R/W" (ref/alt), protein_position = "882"
  if ("protein_position" %in% colnames(d) && "amino_acids" %in% colnames(d)) {
    pos <- as.character(d$protein_position)
    aa <- as.character(d$amino_acids)
    pos[is.na(pos) | pos == ""] <- ""
    aa[is.na(aa) | aa == ""] <- ""
    has_both <- pos != "" & aa != ""
    if (any(has_both)) {
      if (!"HGVSp_Short" %in% colnames(d)) d$HGVSp_Short <- rep(NA_character_, nrow(d))
      aa_split <- strsplit(aa[has_both], "/")
      for (i in which(has_both)) {
        if (length(aa_split[[i]]) >= 2 && aa_split[[i]][1] != "" && aa_split[[i]][2] != "") {
          ref_aa <- aa_split[[i]][1]
          alt_aa <- aa_split[[i]][2]
          aa_pos <- pos[i]
          if (aa_pos != "" && !is.na(aa_pos)) {
            d$HGVSp_Short[i] <- paste0("p.", ref_aa, aa_pos, alt_aa)
          }
        }
      }
    }
  }
  # For TCGA: Protein_position column (numeric)
  if ("Protein_position" %in% colnames(d) && !"HGVSp_Short" %in% colnames(d)) {
    if ("HGVSp_Short" %in% colnames(d)) {
      missing_hgvsp <- is.na(d$HGVSp_Short) | d$HGVSp_Short == ""
      if (any(missing_hgvsp)) {
        prot_pos <- d$Protein_position[missing_hgvsp]
        prot_pos[is.na(prot_pos)] <- ""
        if ("amino_acids" %in% colnames(d)) {
          aa <- as.character(d$amino_acids[missing_hgvsp])
          aa[is.na(aa)] <- ""
          has_both <- prot_pos != "" & aa != ""
          if (any(has_both)) {
            aa_split <- strsplit(aa[has_both], "/")
            for (i in seq_along(which(missing_hgvsp)[has_both])) {
              idx <- which(missing_hgvsp)[has_both][i]
              if (length(aa_split[[i]]) >= 2) {
                ref_aa <- aa_split[[i]][1]
                alt_aa <- aa_split[[i]][2]
                if (ref_aa != "" && alt_aa != "" && prot_pos[i] != "") {
                  d$HGVSp_Short[idx] <- paste0("p.", ref_aa, prot_pos[i], alt_aa)
                }
              }
            }
          }
        }
      }
    }
  }
  # mutation_category, Age, Sex, Risk (Sex: MALE->Male, FEMALE->Female; ELN2017: Adverse/Intermediate/Favorable else Unknown)
  if ("mutation_category" %in% colnames(d)) {
    d$mutation_category <- as.character(d$mutation_category)
    d$mutation_category[d$mutation_category %in% c("RTK/RAS Signaling", "RTK-RAS Signaling")] <- "RTK-RAS Signaling"
  }
  if ("Age" %in% colnames(d)) d$Age <- as.numeric(d$Age)
  if ("Sex" %in% colnames(d)) {
    sx <- as.character(d$Sex)
    d$Sex <- ifelse(sx == "MALE", "Male", ifelse(sx == "FEMALE", "Female", sx))
  } else if ("sex" %in% colnames(d)) {
    sx <- as.character(d$sex)
    d$Sex <- ifelse(sx == "MALE", "Male", ifelse(sx == "FEMALE", "Female", sx))
  }
  if ("Risk" %in% colnames(d)) {
    d$Risk <- as.character(d$Risk)
    d$Risk[is.na(d$Risk) | !d$Risk %in% c("Adverse", "Intermediate", "Favorable")] <- "Unknown"
  } else if ("ELN2017" %in% colnames(d)) {
    eln <- as.character(d$ELN2017)
    d$Risk <- ifelse(eln %in% c("Adverse", "Intermediate", "Favorable"), eln, "Unknown")
  } else if ("ELN_Risk" %in% colnames(d)) {
    eln <- as.character(d$ELN_Risk)
    d$Risk <- ifelse(eln %in% c("Adverse", "Intermediate", "Favorable"), eln, "Unknown")
  } else if ("ELN_risk" %in% colnames(d)) {
    eln <- as.character(d$ELN_risk)
    d$Risk <- ifelse(eln %in% c("Adverse", "Intermediate", "Favorable"), eln, "Unknown")
  } else if ("ELN" %in% colnames(d)) {
    eln <- as.character(d$ELN)
    d$Risk <- ifelse(eln %in% c("Adverse", "Intermediate", "Favorable"), eln, "Unknown")
  } else if ("ELN_2017" %in% colnames(d)) {
    eln <- as.character(d$ELN_2017)
    d$Risk <- ifelse(eln %in% c("Adverse", "Intermediate", "Favorable"), eln, "Unknown")
  } else if ("Karyotype" %in% colnames(d)) d$Risk <- as.character(d$Karyotype)
  else d$Risk <- NA_character_
  # Clinical vars: WBC, Hemoglobin, Platelet, LDH, BM_blast_percent, PB_blast_percent
  if ("WBC" %in% colnames(d)) d$WBC <- as.numeric(d$WBC)
  if ("hemoglobin" %in% colnames(d)) d$Hemoglobin <- as.numeric(d$hemoglobin)
  if ("Platelet" %in% colnames(d)) d$Platelet <- as.numeric(d$Platelet)
  if ("LDH" %in% colnames(d)) d$LDH <- as.numeric(d$LDH)
  if ("bm_blasts" %in% colnames(d)) d$BM_blast_percent <- as.numeric(d$bm_blasts)
  if ("pb_blasts" %in% colnames(d)) d$PB_blast_percent <- as.numeric(d$pb_blasts)
  # Karyotype: preserve if present (Complex / Normal / Other / Unknown)
  if ("Karyotype" %in% colnames(d)) {
    d$Karyotype <- as.character(d$Karyotype)
    d$Karyotype[is.na(d$Karyotype) | d$Karyotype == ""] <- "Unknown"
  } else d$Karyotype <- "Unknown"
  # Keep only columns the app uses (drop extras to avoid confusion); include protein mutation cols for lollipop
  if ("CCF" %in% colnames(d)) d$CCF <- as.numeric(d$CCF)
  if ("CN_at_locus" %in% colnames(d)) d$CN_at_locus <- as.integer(d$CN_at_locus)
  want <- c("Sample", "Gene", "VAF", "CCF", "CN_at_locus", "Cohort", "Subset", "Time_to_OS", "Censor", "variant_type",
    "mutation_category", "AA_change", "HGVSp_Short", "Protein_Change", "protein", "Gene_Group", "Age", "Sex", "Risk", "Karyotype", "WBC", "Hemoglobin", "Platelet", "LDH", "BM_blast_percent", "PB_blast_percent")
  keep <- want[want %in% colnames(d)]
  d <- d[, keep, drop = FALSE]
  # For associations and gene selection: use Gene_Group when set, otherwise Gene (so NA Gene_Group -> show Gene)
  if ("Gene_Group" %in% colnames(d)) {
    d$Gene_for_analysis <- ifelse(is.na(d$Gene_Group) | as.character(d$Gene_Group) == "", as.character(d$Gene), as.character(d$Gene_Group))
  } else {
    d$Gene_for_analysis <- as.character(d$Gene)
  }
  d
}

# Fallback mutation_category for genes not in reference (e.g. lower-frequency genes in Meta AML4); must be defined before AML_Meta_Cohort processing
GENE_MUT_CAT_FALLBACK <- c(
  SMC3 = "Chromatin/Cohesin", SMC1A = "Chromatin/Cohesin", STAG2 = "Chromatin/Cohesin", RAD21 = "Chromatin/Cohesin", NIPBL = "Chromatin/Cohesin",
  ASXL2 = "Chromatin/Cohesin", ASXL1 = "Chromatin/Cohesin", EZH2 = "Chromatin/Cohesin", BCOR = "Chromatin/Cohesin", BCORL1 = "Chromatin/Cohesin", KDM6A = "Chromatin/Cohesin",
  FBXW7 = "Tumor suppressors", TP53 = "Tumor suppressors", PHF6 = "Tumor suppressors", PIGA = "Other",
  CBL = "RTK-RAS Signaling", CBLB = "RTK-RAS Signaling", KIT = "RTK-RAS Signaling", BRAF = "RTK-RAS Signaling", NF1 = "RTK-RAS Signaling", JAK2 = "RTK-RAS Signaling", JAK3 = "RTK-RAS Signaling", MPL = "RTK-RAS Signaling", CSF3R = "RTK-RAS Signaling",
  SETBP1 = "Transcription", IKZF1 = "Transcription", GATA2 = "Transcription", ETV6 = "Transcription", WT1 = "Transcription",
  TET2 = "DNA Methylation", IDH1 = "DNA Methylation", IDH2 = "DNA Methylation", DNMT3A = "DNA Methylation",
  SRSF2 = "Splicing", U2AF1 = "Splicing", ZRSR2 = "Splicing", SF3B1 = "Splicing",
  KRAS = "RTK-RAS Signaling", NRAS = "RTK-RAS Signaling", PTPN11 = "RTK-RAS Signaling",
  RUNX1 = "Transcription", CEBPA = "Transcription", NPM1 = "NPM1"
)

# Load AML_Meta_Cohort_v2 (four cohorts merged) for Meta AML4 tab; fallback to v1
AML_Meta_Cohort <- NULL
meta4_uses_dedicated_file <- FALSE
# Try v2 first (built from raw data sources)
meta4_path_v2 <- file.path(getwd(), "AML_Meta_Cohort_v2.rds")
if (file.exists(meta4_path_v2)) {
  AML_Meta_Cohort <- as.data.frame(readRDS(meta4_path_v2), stringsAsFactors = FALSE)
  meta4_uses_dedicated_file <- TRUE
} else {
  # Fallback to v1 if v2 doesn't exist
  meta4_path_rdata <- file.path(getwd(), "AML_Meta_Cohort.RData")
  meta4_path_rds <- file.path(getwd(), "AML_Meta_Cohort.rds")
  if (file.exists(meta4_path_rdata)) {
    env_meta4 <- new.env()
    load(meta4_path_rdata, envir = env_meta4)
    nms <- ls(env_meta4)
    if ("AML_Meta_Cohort" %in% nms) AML_Meta_Cohort <- get("AML_Meta_Cohort", env_meta4)
    else if (length(nms) == 1L) AML_Meta_Cohort <- get(nms[1], env_meta4)
    if (!is.null(AML_Meta_Cohort)) meta4_uses_dedicated_file <- TRUE
  } else if (file.exists(meta4_path_rds)) {
    AML_Meta_Cohort <- as.data.frame(readRDS(meta4_path_rds), stringsAsFactors = FALSE)
    meta4_uses_dedicated_file <- TRUE
  }
}
if (!is.null(AML_Meta_Cohort) && is.data.frame(AML_Meta_Cohort)) {
  AML_Meta_Cohort <- normalize_AML_Meta_Cohort(AML_Meta_Cohort)
  if (!"mutation_category" %in% colnames(AML_Meta_Cohort)) AML_Meta_Cohort$mutation_category <- NA_character_
  # Assign mutation_category from GENE_MUT_CAT_FALLBACK lookup
  miss <- is.na(AML_Meta_Cohort$mutation_category) | as.character(AML_Meta_Cohort$mutation_category) == ""
  if (any(miss)) {
    g <- as.character(AML_Meta_Cohort$Gene[miss])
    AML_Meta_Cohort$mutation_category[miss] <- ifelse(g %in% names(GENE_MUT_CAT_FALLBACK), unname(GENE_MUT_CAT_FALLBACK[g]), "Other")
  }
} else {
  AML_Meta_Cohort <- final_data_matrix
}
tmp2 <- tempfile(fileext = ".csv")
utils::write.csv(AML_Meta_Cohort, tmp2, row.names = FALSE)
AML_Meta_Cohort <- utils::read.csv(tmp2, stringsAsFactors = FALSE)
unlink(tmp2)

# Gene -> chromosome arm mapping for CCF computation (used at runtime if CCF not in RDS)
GENE_CHR_MAP <- c(
  DNMT3A = "2p", IDH1 = "2q", SF3B1 = "2q",
  GATA2 = "3q",
  TET2 = "4q", KIT = "4q",
  NPM1 = "5q",
  ETV6 = "12p", KRAS = "12p",
  EZH2 = "7q",
  RAD21 = "8q",
  JAK2 = "9p", CDKN2A = "9p",
  SMC3 = "10q",
  WT1 = "11p", CBL = "11q", MLL = "11q", KMT2A = "11q", ATM = "11q",
  PTPN11 = "12q",
  FLT3 = "13q", "FLT3-ITD" = "13q", "FLT3-TKD" = "13q", "FLT3-Other" = "13q",
  IDH2 = "15q",
  CBFB = "16q",
  TP53 = "17p", NF1 = "17q", SRSF2 = "17q",
  SETBP1 = "18q",
  CEBPA = "19q",
  ASXL1 = "20q",
  RUNX1 = "21q", U2AF1 = "21q",
  NRAS = "1p", BCOR = "Xp", SMC1A = "Xp", KDM6A = "Xp", ZRSR2 = "Xp",
  STAG2 = "Xq", BCORL1 = "Xq", PHF6 = "Xq", ATRX = "Xq",
  ASXL2 = "2p", CSF3R = "1p", FBXW7 = "4q", PIGA = "Xp",
  IKZF1 = "7p", BRAF = "7q", JAK3 = "19p", MPL = "1p", CBLB = "3q"
)

# Runtime CCF computation for data that lacks precomputed CCF
# Assumes diploid CN=2 for autosomes, CN=1 for male X-linked genes
compute_ccf_runtime <- function(df) {
  if (!"VAF" %in% colnames(df)) return(df)
  n <- nrow(df)
  ccf <- rep(NA_real_, n)
  cn_used <- rep(NA_integer_, n)
  sex_col <- if ("Sex" %in% colnames(df)) df$Sex else rep(NA_character_, n)
  gene_col <- as.character(df$Gene)
  vaf_col <- df$VAF

  for (i in seq_len(n)) {
    vaf <- vaf_col[i]
    if (is.na(vaf)) next
    gene <- gene_col[i]
    sex <- as.character(sex_col[i])
    chr_arm <- GENE_CHR_MAP[gene]
    if (is.na(chr_arm)) {
      cn <- 2L
    } else {
      chr_num <- sub("[pq]$", "", chr_arm)
      if (chr_num == "X" && !is.na(sex) && sex == "Male") cn <- 1L
      else if (chr_num == "Y") cn <- ifelse(!is.na(sex) && sex == "Male", 1L, 0L)
      else cn <- 2L
    }
    cn_used[i] <- cn
    if (cn > 0) ccf[i] <- min(vaf / 100 * cn * 100, 100)
  }
  df$CCF <- round(ccf, 2)
  df$CN_at_locus <- cn_used
  df
}

# Compute CCF at runtime if not already present in the loaded data
if (!"CCF" %in% colnames(AML_Meta_Cohort) && "VAF" %in% colnames(AML_Meta_Cohort)) {
  AML_Meta_Cohort <- compute_ccf_runtime(AML_Meta_Cohort)
}

# Color palettes
PAL_COHORT <- c(
  Tyner = "#0073C2FF", TCGA = "#EFC000FF", Majeti = "#868686FF", Papaemmanuil = "#CD534CFF",
  Lindsley = "#7AA6DCFF", Wang = "#E64B35FF", Au = "#4DBBD5FF", Welch = "#00A087FF",
  Garg = "#3C5488FF", Greif = "#F39B7FFF", Hirsch = "#7E6148FF", Huet = "#B09C85FF"
)
PAL_SUBSET <- c(
  de_novo = "#C16622FF", secondary = "#767676FF", relapse = "#800000FF",
  other = "#FFA319FF", therapy = "#8F3931FF"
)
PAL_MUT_CAT <- c(
  "DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF",
  "RTK-RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF",
  "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF",
  "Other" = "#9E9E9EFF"
)
PAL_RISK <- c(Adverse = "#E64B35FF", Intermediate = "#8491B4FF", Favorable = "#00A087FF", Unknown = "#767676FF")
# Paper Figure 1B oncoprint colors (exact from manuscript code)
ONCO_COHORT_COL <- c(Tyner = "#0073C2FF", TCGA = "#EFC000FF", Majeti = "#868686FF", Papaemmanuil = "#CD534CFF",
  Lindsley = "#7AA6DCFF", Wang = "#E64B35FF", Au = "#4DBBD5FF", Welch = "#00A087FF", Garg = "#3C5488FF",
  Greif = "#F39B7FFF", Li = "#8491B4FF", Shlush = "#91D1C2FF", Parkin = "#DC0000FF", Hirsch = "#7E6148FF", Huet = "#B09C85FF")
ONCO_SUBSET_COL <- c(de_novo = "#C16622FF", secondary = "#767676FF", relapse = "#800000FF", other = "#FFA319FF",
  Remission = "#8A9045FF", Residual = "#155F83FF", therapy = "#8F3931FF", MDS = "#58593FFF", Unknown = "gray90")
ONCO_VAR_COL <- c(Deletion = "#374E55FF", INDEL = "#DF8F44FF", Insertion = "#00A1D5FF", ITD = "#79AF97FF",
  SNV = "#B24745FF", Splicing = "#6A6599FF", Unknown = "#80796BFF")
# Meta AML4 tab: uniform cohort/subset/variant colors (four cohorts: Beat AML, TCGA, AML-SG, UK-NCRI)
META4_COHORT_COL <- c("Beat AML" = "#155F83FF", "TCGA" = "#EFC000FF", "AML-SG" = "#CD534CFF", "UK-NCRI" = "#00A087FF")
META4_SUBSET_COL <- c("De novo" = "#C16622FF", "Secondary" = "#767676FF", "Relapse" = "#800000FF",
  "Therapy" = "#8F3931FF", "Other" = "#FFA319FF",
  De_novo = "#C16622FF", de_novo = "#C16622FF", secondary = "#767676FF", relapse = "#800000FF", therapy = "#8F3931FF", other = "#FFA319FF",
  oAML = "#FFA319FF", tAML = "#8F3931FF", unknown = "#8A9045FF")
META4_VAR_COL <- c(Deletion = "#374E55FF", Indel = "#DF8F44FF", Insertion = "#00A1D5FF", ITD = "#79AF97FF",
  SNV = "#B24745FF", Splicing = "#6A6599FF", PTD = "tan", Other = "#80796BFF", Unknown = "#80796BFF")

# Load BeatAML2 helpers (only if data exists)
beataml2_script <- file.path(getwd(), "load_beataml2.R")
if (file.exists(beataml2_script)) source(beataml2_script, local = FALSE)

# FLT3 annotation
prepare_data <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$Gene <- as.character(df$Gene)
  df$variant_type <- as.character(df$variant_type)
  flt3_itd <- df$Gene == "FLT3" & df$variant_type %in% c("ITD", "INDEL", "Indel")
  flt3_tkd <- df$Gene == "FLT3" & df$variant_type %in% c("SNV", "Deletion", "other", "Other", "Unknown")
  df$Gene[flt3_itd] <- "FLT3-ITD"
  df$Gene[flt3_tkd] <- "FLT3-TKD"
  df
}

# Normalize protein mutation strings from four datasets to HGVSp_Short and extract AA position for maftools lollipop
# Handles: p.Arg123Gly, p.R123G, Arg123Gly, R123G, 123, p.Val123_Ala125del, p.?, etc.
normalize_protein_change_for_maf <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "p.?"
  hgvsp <- character(length(x))
  position <- integer(length(x))
  for (i in seq_along(x)) {
    s <- trimws(x[i])
    if (s == "" || s == "p.?" || tolower(s) == "na" || s == ".") {
      hgvsp[i] <- "p.?"
      position[i] <- NA_integer_
      next
    }
    if (!grepl("^p\\.", s, ignore.case = TRUE)) s <- paste0("p.", gsub("^p\\.?", "", s))
    hgvsp[i] <- s
    pos_match <- regexpr("[0-9]+", s)
    if (pos_match > 0) {
      position[i] <- as.integer(substr(s, pos_match, pos_match + attr(pos_match, "match.length") - 1L))
    } else {
      position[i] <- NA_integer_
    }
  }
  list(hgvsp = hgvsp, position = position)
}

# UI
ui <- fluidPage(
  title = "Meta AML Explorer",
  tags$head(tags$style(HTML(
    ".nav-tabs{font-weight:600}.well{background:#f8f9fa;border-radius:8px}
     .app-banner{display:flex;align-items:center;justify-content:space-between;padding:12px 20px;background:#374e55;color:#fff;width:100vw;position:relative;left:50%;right:50%;margin-left:-50vw;margin-right:-50vw;margin-bottom:15px;}
     .app-banner .banner-title{font-size:22px;font-weight:700;margin:0;}
     .app-banner .banner-contact{display:flex;align-items:center;gap:12px;}
     .app-banner .banner-contact img{border-radius:50%;width:42px;height:42px;object-fit:cover;}
     .app-banner .banner-contact a{color:#fff;text-decoration:none;font-size:14px;}
     .app-banner .banner-contact a:hover{color:#df8f44;text-decoration:underline;}
     .welcome-page h2{color:#374e55;margin-top:1.2em;}
     .welcome-page h2:first-of-type{margin-top:0;}
     .welcome-page p{font-size:15px;line-height:1.6;}
     .welcome-page ul{font-size:15px;line-height:1.7;}
     #main_nav .nav-tabs{margin-bottom:0;border-bottom:2px solid #374e55;}
     #main_nav .nav-tabs .nav-link{font-weight:600;color:#374e55;}
     #main_nav .nav-tabs .nav-link.active{background:#374e55;color:#fff;border-color:#374e55;}
     #drug_subset_wrap .shiny-input-container{width:140px !important;min-width:140px;}
     /* Loading overlay spinner */
     #data-loading-overlay .loading-spinner{width:48px;height:48px;margin:0 auto 16px;border:4px solid #e2e8f0;border-top-color:#374e55;border-radius:50%;animation:data-load-spin 0.9s linear infinite;}
     @keyframes data-load-spin{to{transform:rotate(360deg);}}"
  ))),
  tags$head(
    tags$script(async = NA, src = "https://www.googletagmanager.com/gtag/js?id=G-BHRQ6LT7GG"),
    tags$script(HTML("
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-BHRQ6LT7GG');
    "))
  ),
  div(id = "data-loading-overlay", style = "display: none; position: fixed; top: 0; left: 0; right: 0; bottom: 0; background: rgba(0,0,0,0.15); z-index: 9999; align-items: center; justify-content: center;",
    tags$div(style = "background: #fff; border-radius: 12px; box-shadow: 0 4px 24px rgba(0,0,0,0.15); padding: 28px 36px; text-align: center; min-width: 260px; max-width: 320px;",
      div(class = "loading-spinner"),
      tags$div(style = "font-size: 18px; color: #374e55; margin-bottom: 8px; font-weight: 600;", "Loading data..."),
      tags$p(style = "color: #666; font-size: 14px; margin: 0;", "Preparing analyses. This may take a moment.")
    )
  ),
  tags$script(HTML("
    (function() {
      var overlayShownAt = 0;
      var minOverlayMs = 1500;
      var hideTimeout = null;
      var delayedShowTimeout = null;
      var maxOverlayTimeout = null;
      var overlayVisible = false;
      var currentMainNav = 'about';
      var LOADING_DELAY_MS = 1500;
      var MAX_OVERLAY_MS = 30000;

      function showOverlay() {
        if (overlayVisible) return;
        overlayVisible = true;
        overlayShownAt = Date.now();
        if (hideTimeout) { clearTimeout(hideTimeout); hideTimeout = null; }
        if (maxOverlayTimeout) { clearTimeout(maxOverlayTimeout); maxOverlayTimeout = null; }
        var el = document.getElementById('data-loading-overlay');
        if (el) el.style.display = 'flex';
        maxOverlayTimeout = setTimeout(function() {
          maxOverlayTimeout = null;
          if (overlayVisible) hideOverlay();
        }, MAX_OVERLAY_MS);
      }

      function hideOverlay() {
        if (!overlayVisible) return;
        overlayVisible = false;
        if (maxOverlayTimeout) { clearTimeout(maxOverlayTimeout); maxOverlayTimeout = null; }
        var elapsed = Date.now() - overlayShownAt;
        var wait = Math.max(0, minOverlayMs - elapsed);
        hideTimeout = setTimeout(function() {
          hideTimeout = null;
          $('#data-loading-overlay').hide();
        }, wait);
      }

      function cancelPending() {
        if (delayedShowTimeout) { clearTimeout(delayedShowTimeout); delayedShowTimeout = null; }
      }

      $(document).on('shiny:inputchanged', function(event) {
        if (event.name === 'main_nav') {
          currentMainNav = event.value;
          if (currentMainNav === 'about') {
            cancelPending();
            if (overlayVisible) hideOverlay();
          }
        }
      });

      $(document).on('shiny:busy', function() {
        if (currentMainNav === 'about') return;
        cancelPending();
        delayedShowTimeout = setTimeout(function() {
          delayedShowTimeout = null;
          showOverlay();
        }, LOADING_DELAY_MS);
      });

      $(document).on('shiny:idle', function() {
        cancelPending();
        if (overlayVisible) hideOverlay();
      });

      Shiny.addCustomMessageHandler('dataReady', function(tab) {});
    })();
  ")),
  div(class = "app-banner",
    div(class = "banner-title", "Meta AML Explorer")
  ),
  tabsetPanel(
    id = "main_nav",
    tabPanel("About", value = "about",
      div(class = "welcome-page", style = "max-width: 800px; margin: 0 auto; padding: 24px 15px;",
        h2("Welcome to Meta AML Explorer"),
        p("This site provides interactive exploration of acute myeloid leukemia (AML) mutational and clinical outcomes data. The ", strong("Meta AML4"), " tab offers the following analyses:"),
        tags$ul(
          tags$li("Cohort selection by AML type (e.g. de novo, secondary), dataset (e.g. UK-NCRI), karyotype (NK/Complex/Other), and minimum gene frequency (% of samples)."),
          tags$li("Single mutation associations (clinical variables, survival, hazard ratios)"),
          tags$li("Co-mutation (heatmap, odds-ratios, Kaplan-Meier by 2–3 genes)"),
          tags$li("Variant allele frequency (VAF) associations (clinical variables, survival, pairwise scatterplots)"),
          tags$li("Drug Sensitivity (Mutation and VAF associations for BeatAML2 data)")
        ),
        h2("Meta AML4 tab"),
        p("Meta AML4 merges ", strong("four of the largest molecularly profiled and clinically annotated AML datasets"), " into a single combined cohort of ", strong("~4,660 patients"), "."),
        div(style = "background: #f8f9fa; border-left: 4px solid #0066cc; padding: 14px 18px; margin: 12px 0; border-radius: 0 6px 6px 0;",
          p(strong("Datasets:"), style = "margin-top: 0; margin-bottom: 8px;"),
          p(style = "margin: 4px 0; line-height: 1.6;", strong("UK-NCRI"), " (2,113 patients) | ", tags$a("Tazi et al. (2022) (Nature Communications)", href = "https://www.nature.com/articles/s41467-022-32103-8", target = "_blank"), " | ", tags$a("Data", href = "https://github.com/papaemmelab/Tazi_NatureC_AML", target = "_blank")),
          p(style = "margin: 4px 0; line-height: 1.6;", strong("AML-SG"), " (1,540 patients) | ", tags$a("Papaemmanuil et al. (2016), NEJM", href = "https://www.nejm.org/doi/10.1056/NEJMoa1516192", target = "_blank"), " | ", tags$a("Data", href = "https://github.com/gerstung-lab/AML-multistage", target = "_blank")),
          p(style = "margin: 4px 0; line-height: 1.6;", strong("Beat AML"), " (805 patients) | ", tags$a("Tyner et al. (2022), Cancer Cell", href = "https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00312-9", target = "_blank"), " | ", tags$a("Data", href = "https://biodev.github.io/BeatAML2/", target = "_blank")),
          p(style = "margin: 4px 0; line-height: 1.6; margin-bottom: 0;", strong("TCGA LAML"), " (200 patients) | ", tags$a("NEJM", href = "https://www.nejm.org/doi/full/10.1056/NEJMoa1301689", target = "_blank"), " | ", tags$a("Data", href = "https://www.cbioportal.org/", target = "_blank"))
        ),
        h2("Caveats"),
        p("This site is intended for research purposes only and is in active development. Some initial data loadings may take a few seconds."),
        p("Estimated cancer cell fractions (CCF) are derived from VAF and cytogenetic copy number: CCF = VAF x CN (capped at 100%). Copy number is determined from UK-NCRI cytogenetics data, AML-SG clinical cytogenetic flags and ISCN karyotype strings, and Beat AML karyotype data. Sex-chromosome loci are corrected for patient sex (hemizygous CN=1 for males). Where cytogenetics are unavailable, diploid CN=2 is assumed. CCF values are estimates and do not account for tumor purity or subclonal copy number heterogeneity."),
        h2("About the author"),
        div(style = "display: flex; align-items: center; gap: 16px; margin-bottom: 10px;",
          tags$img(src = "linkedin_pic%20copy.jpeg", alt = "Brooks Benard", style = "border-radius: 50%; width: 64px; height: 64px; object-fit: cover;"),
          div(
            p(style = "margin: 0; font-size: 16px; font-weight: 600;", "Brooks Benard"),
            p(style = "margin: 2px 0 0; font-size: 14px; color: #555;", "Stanford University"),
            p(style = "margin: 4px 0 0;", tags$a(href = "mailto:bbenard@stanford.edu", "bbenard@stanford.edu"), " · ", tags$a("Website", href = "https://brooksbenard.github.io/", target = "_blank"), " · ", tags$a("GitHub", href = "https://github.com/brooksbenard", target = "_blank"))
          )
        ),
        hr(),
        p(style = "font-size: 0.9em; color: #666;",
          "© ", format(Sys.Date(), "%Y"), " Brooks Benard. Licensed under ",
          tags$a("MIT License", href = "https://opensource.org/licenses/MIT", target = "_blank"), "."
        )
      )
    ),
    tabPanel("Meta AML4", value = "meta_aml4",
      conditionalPanel(
        condition = "input.main_nav === 'analyses' || input.main_nav === 'meta_aml4'",
        sidebarLayout(
          sidebarPanel(
            width = 2,
            div(style = "background: #2c3e50; color: #fff; padding: 12px 14px; margin: -15px -15px 15px -15px; border-radius: 4px 4px 0 0;",
              h4(style = "margin: 0; font-weight: 700; font-size: 1.15em;", "Select cohort for analysis:")
            ),
            selectInput("subset", "AML Subset:", 
              choices = c("All", "De novo", "Secondary", "Relapse", "Therapy", "Other"), selected = "De novo"),
            selectInput("cohort", "Dataset:", 
              choices = local({
                coh <- c("All", sort(unique(as.character(na.omit(final_data_matrix$Cohort)))))
                coh[coh == "Beat_AML"] <- "Beat AML"
                coh
              }), selected = "All"),
            selectInput("karyotype", "Karyotype:", 
              choices = c("All", "Complex", "Normal", "Other", "Unknown"), selected = "All"),
            radioButtons("vaf_metric", "Allele metric:",
              choices = c("VAF" = "VAF", "CCF" = "CCF"), selected = "VAF", inline = TRUE),
            conditionalPanel(
              condition = "input.main_nav === 'meta_aml4'",
              hr(),
              div(style = "border-left: 3px solid #95a5a6; padding-left: 12px; margin-top: 8px; margin-bottom: 4px;",
                p(style = "margin: 0 0 6px 0; font-weight: 600; color: #34495e;", "Only interested in a specific gene?"),
                p(style = "margin: 0 0 8px 0; font-size: 13px; color: #5d6d7e;", "Summarize all associations for a single gene:"),
                selectInput("gene_summary", NULL, choices = c("Select gene..." = ""), selected = "")
              )
            )
          ),
          mainPanel(
            width = 10,
            uiOutput("meta4_fallback_banner"),
            uiOutput("main_tabs_ui")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {

  # Cache filtered data per tab to avoid recomputation when switching back
  cached_filtered_data <- reactiveValues(analyses = NULL, meta_aml4 = NULL)
  cached_data_params <- reactiveValues(analyses = NULL, meta_aml4 = NULL)
  # Preserve selected sub-tab when filters change (so we don't jump back to Overview)
  selected_meta4_tab <- reactiveVal(NULL)
  selected_analyses_tab <- reactiveVal(NULL)
  # Preserve inner tab inside "All Gene Associations" (Single, Co-mutation, VAF, Drug Sensitivity)
  selected_meta4_allgene_tab <- reactiveVal(NULL)
  # Preserve sub-tab inside Single Gene Associations (Clinical, Co-mutation, VAF, Drug Sensitivity)
  selected_gene_summary_sub_tab <- reactiveVal("clinical")

  # Precomputed Drug Sensitivity results (Mut vs WT, VAF vs AUC correlations, LOOCV) per subset to reduce memory on deploy
  precomputed_drug <- NULL
  tryCatch({
    path <- file.path(getwd(), "drug_sensitivity_precomputed.rds")
    if (file.exists(path)) precomputed_drug <- readRDS(path)
  }, error = function(e) precomputed_drug <- NULL)

  # Format p/q for display: exact numeric values (no < 0.01 style)
  format_pval_display <- function(p) {
    if (is.null(p) || length(p) == 0) return(p)
    out <- character(length(p))
    for (i in seq_along(p)) {
      if (is.na(p[i])) { out[i] <- NA_character_; next }
      out[i] <- format(p[i], scientific = TRUE, drop0trailing = TRUE)
    }
    out
  }

  # Map Gene_for_analysis (refined labels e.g. IDH2 R140, NRAS G12/13, CEBPA bi-allelic) to base gene name(s) used in Beat AML2 data
  gene_to_beataml_bases <- function(g) {
    if (is.null(g) || is.na(g) || g == "") return(character(0))
    g <- as.character(g)
    if (grepl("^IDH2", g)) return("IDH2")
    if (grepl("^NRAS", g)) return("NRAS")
    if (grepl("^CEBPA", g)) return("CEBPA")
    if (grepl("^U2AF1", g)) return("U2AF1")
    if (g %in% c("FLT3-ITD", "FLT3-TKD", "FLT3-Other")) return(g)
    return(g)
  }

  # Data source and Meta4 flag must be reactives on input$main_nav so outputs re-run
  # in the same flush as the tab change (no observer race).
  # Use isolate to prevent outputs from rendering until tab is properly set
  current_data_src <- reactive({
    nav <- req(input$main_nav)
    if (!nav %in% c("analyses", "meta_aml4")) return(NULL)
    if (identical(nav, "meta_aml4")) AML_Meta_Cohort else final_data_matrix
  })

  is_meta4 <- reactive(identical(input$main_nav, "meta_aml4"))

  # Render main tabs conditionally - only render the active tabsetPanel to avoid duplicate outputs and improve performance
  # Preserve selected tab via selected_meta4_tab/selected_analyses_tab so filter changes update the content without jumping back to Overview
  output$main_tabs_ui <- renderUI({
    nav <- input$main_nav
    if (identical(nav, "meta_aml4")) {
      df <- filtered_data()
      genes <- filtered_genes()
      gene_choices <- if (length(genes) == 0) {
        c("Select gene..." = "")
      } else {
        c("Select gene..." = "", setNames(genes, genes))
      }
      current_tab <- selected_meta4_tab()
      if (is.null(current_tab) || !current_tab %in% c("Cohort Overview", "All Gene Associations", "Single Gene Associations")) current_tab <- "Cohort Overview"
      inner_tab <- selected_meta4_allgene_tab()
      if (is.null(inner_tab) || !inner_tab %in% c("Single", "Co-mutation", "VAF", "Drug Sensitivity")) inner_tab <- "Single"
      tabsetPanel(
        id = "main_tabs_meta4",
        selected = current_tab,
        tabPanel("Cohort Overview",
          fluidRow(
            column(3, wellPanel(h4("Selected Cohort"), tableOutput("summary_table"))),
            column(9, wellPanel(h4("Samples by Cohort & Subset"), plotOutput("cohort_plot", height = 280)))
          ),
          fluidRow(
            column(6, wellPanel(h4("Survival by Dataset"), plotOutput("overview_surv_dataset_plot", height = 400))),
            column(6, wellPanel(h4("Age Distribution by Dataset"), plotOutput("overview_age_dataset_plot", height = 400)))
          ),
          fluidRow(
            column(6, wellPanel(h4("ELN Risk by Dataset"), plotOutput("overview_eln_dataset_plot", height = 400))),
            column(6, wellPanel(h4("Recurrently Mutated Genes"), p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Genes mutated in more than 1% of samples."), plotOutput("overview_gene_freq_plot", height = 400)))
          )
        ),
        tabPanel("All Gene Associations",
          tabsetPanel(
            id = "meta4_all_gene_tabs",
            selected = inner_tab,
            tabPanel("Single mutation", value = "Single",
              fluidRow(
                column(12, wellPanel(
                  fluidRow(column(8, h4("Clinical Variable by Mutation Status")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_clinical_plot", "Download plot")))),
                  selectInput("clin_var", "Variable:", choices = c("WBC", "Age", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent")),
                  p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", span(style = "color: #8B0000; font-weight: bold;", "Red"), " = mutated; ", span(style = "color: #808080; font-weight: bold;", "Gray"), " = WT. Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant (Wilcoxon)."),
                  plotOutput("clinical_plot", height = 400)))
              ),
              fluidRow(
                column(6, wellPanel(fluidRow(column(8, h4("Hazard Ratios")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_forest_plot", "Download plot")))), p(style = "font-size: 14px; color: #666; margin-bottom: 6px;", "Individual gene hazard ratios were calculated between mutated and wild-type patients using a univariate Cox proportional hazards model."), plotOutput("forest_plot", height = 500), p(style = "font-size: 13px; color: #666; margin-top: 6px;", span(style = "color: #762a83; font-weight: bold;", "Purple"), " = significant HR > 1 (higher risk); ", span(style = "color: #1b7837; font-weight: bold;", "Green"), " = significant HR < 1 (lower risk); ", span(style = "color: #9E9E9E; font-weight: bold;", "Gray"), " = not significant."))),
                column(6, wellPanel(fluidRow(column(8, h4("Kaplan-Meier")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_survival_plot", "Download plot")))), selectInput("surv_gene", NULL, choices = gene_choices), plotOutput("survival_plot", height = 500), p(style = "font-size: 13px; color: #666; margin-top: 6px;", span(style = "color: #8B0000; font-weight: bold;", "Red"), " = mutated; ", span(style = "color: #4D4D4D; font-weight: bold;", "Gray"), " = WT.")))
              ),
              fluidRow(column(12, wellPanel(fluidRow(column(8, h4("Survival Summary")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_survival_table", "Download table")))), p(style = "font-size: 13px; color: #666; margin-bottom: 8px;", "Cox proportional hazards (mutated vs WT) per gene; genes with ≥10 mutated and ≥10 WT. HR and 95% CI shown. ", em("p_adj"), " = false discovery rate (FDR, Benjamini–Hochberg)."), DTOutput("survival_table"))))
            ),
            tabPanel("Co-mutation", value = "Co-mutation",
              tabsetPanel(
                tabPanel("Oncoprint",
                  fluidRow(
                    column(12, wellPanel(
                      fluidRow(column(8, h4("Oncoprint")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_oncoprint_plot", "Download plot")))),
                      p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "All recurrent mutations in the selected cohort. Each column is a sample and genes are ordered by decreasing frequency."),
                      plotOutput("oncoprint_plot", height = 500)
                    ))
                  )
                ),
                tabPanel("Odds Ratio and Survival",
                  fluidRow(
                    column(6, wellPanel(fluidRow(column(8, h4("Co-mutation Odds Ratio")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_cooccurrence_plot", "Download plot")))), p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = co-occurring (OR > 1); ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = mutually exclusive (OR < 1); ", span(style = "color: #808080; font-weight: bold;", "Gray"), " = not significant."), div(style = "max-width: 100%; overflow: auto;", plotOutput("cooccurrence_plot", height = 500)))),
                    column(6, wellPanel(fluidRow(column(8, h4("Mutation Co-occurrence Statistics")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_cooccurrence_table", "Download table")))), p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", "Odds ratio (OR) and FDR q-value. OR > 1 = genes tend to co-occur; OR < 1 = mutually exclusive."), div(style = "height: 500px; overflow-y: auto;", DTOutput("cooccurrence_table"))))
                  ),
                  fluidRow(
                    column(12, wellPanel(
                      fluidRow(column(8, h4("Kaplan-Meier by Co-mutation Status")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_comut_survival_plot", "Download plot")))),
                      fluidRow(
                        column(4,
                          div(style = "background: #eaecef; border-radius: 8px; padding: 14px 16px; margin-top: 4px;",
                            h5(style = "margin-top: 0;", "Gene Selection"),
                            radioButtons("comut_n", "Number of genes:", choices = c("2 genes" = 2, "3 genes" = 3), selected = 2),
                            selectInput("comut_gene1", "Gene 1:", choices = gene_choices),
                            selectInput("comut_gene2", "Gene 2:", choices = gene_choices),
                            conditionalPanel(condition = "input.comut_n == 3",
                              selectInput("comut_gene3", "Gene 3:", choices = gene_choices))
                          ),
                          div(style = "margin-top: 12px;",
                            h5("Mutation Intersection (UpSet Plot)", style = "margin-top: 0;"),
                            p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", "Overlap of the selected genes across samples."),
                            plotOutput("upset_plot", height = 350),
                            div(style = "margin-top: 4px;", downloadButton("download_upset_plot", "Download plot"))
                          )
                        ),
                        column(8, uiOutput("comut_survival_plot_ui"))
                      )
                    ))
                  )
                )
              )
            ),
            tabPanel("Clonality", value = "VAF",
              tabsetPanel(
                tabPanel("Single Gene",
                  fluidRow(
                    column(6, wellPanel(
                      fluidRow(column(8, h4("VAF / CCF by Gene")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_gene_plot", "Download plot")))),
                      plotOutput("vaf_gene_plot", height = 500))),
                    column(6, wellPanel(
                      fluidRow(column(8, h4("VAF / CCF and Survival: Hazard Ratios (MaxStat)")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_survival_plot", "Download plot")))),
                      p(em("HR and 95% CI from Cox model using optimal MaxStat threshold per gene. Threshold computed from the selected allele metric (VAF or CCF). Genes with \u226520 mutated samples only.")),
                      plotOutput("vaf_survival_plot", height = 450)))
                  ),
                  fluidRow(
                    column(6, wellPanel(
                      fluidRow(column(8, h4("VAF / CCF by Mutation Category")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_category_plot", "Download plot")))),
                      plotOutput("vaf_category_plot", height = 450))),
                    column(6, wellPanel(
                      fluidRow(column(8, h4("Mutation HR vs MaxStat VAF/CCF HR")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_hr_vs_maxstat_plot", "Download plot")))),
                      p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Each point is a gene. Y-axis: mutation HR (mutated vs WT). X-axis: MaxStat VAF HR (high vs low VAF). Point size reflects optimal VAF cutpoint; color indicates significance (FDR q < 0.1)."),
                      plotOutput("hr_vs_maxstat_plot", height = 400)
                    ))
                  ),
                  fluidRow(
                    column(6, wellPanel(
                      fluidRow(column(8, h4("VAF / CCF Distribution by Dataset")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_cohort_plot", "Download plot")))),
                      plotOutput("vaf_cohort_plot", height = 450)))
                  )
                ),
                tabPanel("Pairwise Gene",
                  fluidRow(
                    column(6, wellPanel(h4("VAF / CCF Scatter (Clonal Ordering)"),
                      p("Compare VAF of two genes. Points above line = Gene 1 before Gene 2."),
                      selectInput("vaf_gene1", "Gene 1:", choices = gene_choices), selectInput("vaf_gene2", "Gene 2:", choices = gene_choices),
                      plotOutput("vaf_scatter_plot", height = 450))),
                    column(6, wellPanel(h4("Survival by Clonal Order"),
                      p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Kaplan-Meier comparing patients where Gene 1 has higher VAF (acquired first) vs Gene 2 has higher VAF (acquired first)."),
                      plotOutput("vaf_clonal_km_plot", height = 450)))
                  ),
                  fluidRow(
                    column(6, wellPanel(
                      fluidRow(column(8, h4("Clonal Ordering Survival: All Gene Pairs")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_clonal_hr_summary_plot", "Download plot")))),
                      p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Cox HR (95% CI) comparing survival when Gene A is clonally dominant vs Gene B, for all co-mutated gene pairs with \u226510 non-ambiguous orderings per group. Point size reflects total non-ambiguous cases."),
                      plotOutput("clonal_hr_summary_plot", height = 600)
                    )),
                    column(6, wellPanel(
                      fluidRow(column(8, h4("VAF vs CCF")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_vs_ccf_plot", "Download plot")))),
                      p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Each point is a mutation. Diagonal = CCF equals VAF (diploid, CN = 2). Points above the diagonal have CN > 2 (gain); points below have CN < 2 (loss)."),
                      plotOutput("vaf_vs_ccf_plot", height = 550)
                    ))
                  )
                )
              )
            ),
            tabPanel("Drug Sensitivity", value = "Drug Sensitivity",
              p("Data for mutation and VAF correlations with inhibitor sensitivity from BeatAML2 data.",
                a("biodev.github.io/BeatAML2", href = "https://biodev.github.io/BeatAML2/", target = "_blank")),
              fluidRow(
                column(12, wellPanel(
                  div(id = "drug_subset_wrap", selectInput("drug_subset", "Select Beat AML cohort:", choices = c("All", "de_novo", "secondary"), selected = "de_novo")),
                  uiOutput("drug_summary_ui")
                ))
              ),
              fluidRow(
                column(12, wellPanel(
                  fluidRow(column(8, h4("Mut vs WT AUC (≥5 mut samples)")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_mut_wt_heatmap", "Download plot")))),
                  uiOutput("drug_mut_wt_summary_ui"),
                  p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = mutated samples more resistant (higher AUC); ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = mutated samples more sensitive (lower AUC)."),
                  p(em("Delta AUC = mean AUC (mut) - mean AUC (WT). Only gene-inhibitor pairs with q < 0.1 (FDR) shown."), style = "margin-bottom: 4px;"),
                  div(style = "margin: 0; padding: 0; line-height: 0;", plotOutput("drug_mut_wt_heatmap", height = 400))
                ))
              ),
              fluidRow(
                column(12, wellPanel(
                  fluidRow(column(8, h4("VAF vs AUC correlations")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_summary_dotplot", "Download plot")))),
                  uiOutput("drug_vaf_auc_summary_ui"),
                  p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = negative slope (higher VAF \u2192 lower AUC); ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = positive slope (higher VAF \u2192 higher AUC). * = q < 0.1 (FDR)."),
                  plotOutput("drug_summary_dotplot", height = 500),
                  p(em("* q < 0.1 (FDR)"), style = "font-size: 13px; color: #666;")
                ))
              ),
              fluidRow(
                column(6, wellPanel(
                  fluidRow(column(8, h4("Individual Mutation and VAF Associations")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_scatter_plot", "Download plot")))),
                  p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", "Boxplot: ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = mut more resistant; ", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = mut more sensitive. Scatter: ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = positive VAF\u2013AUC; ", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = negative."),
                  p(em("Gene and inhibitor options match the subset used in Mut vs WT and VAF vs AUC analyses above.")),
                  uiOutput("drug_scatter_selects_ui"),
                  fluidRow(
                    column(6, plotOutput("drug_scatter_boxplot", width = "100%", height = 320)),
                    column(6, plotOutput("drug_scatter_plot", width = "100%", height = 320))
                  )
                )),
                column(6, wellPanel(
                  fluidRow(column(8, h4("Leave-One-Out Cross Validation for VAF vs AUC analysis")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_loo_heatmap", "Download plot")))),
                  p(em("RMSE = leave-one-out prediction error in AUC units. Lower RMSE = better VAF\u2013AUC fit."), style = "font-size: 13px; color: #666; margin-bottom: 4px;"),
                  plotOutput("drug_loo_heatmap", height = 350)
                ))
              ),
              fluidRow(
                column(6, wellPanel(fluidRow(column(8, h4("Mut vs WT Statistics")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_mut_wt_table", "Download table")))), p(em("Delta AUC = mean AUC (mut) \u2212 mean AUC (WT). q < 0.1 (FDR) shown."), style = "font-size: 13px; color: #666; margin-bottom: 4px;"), div(style = "height: 350px; overflow-y: auto;", DTOutput("drug_mut_wt_summary_table")))),
                column(6, wellPanel(fluidRow(column(8, h4("VAF vs AUC Statistics")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_correlation_table", "Download table")))), p(em("LOOCV RMSE = leave-one-out prediction error in AUC units (overfitting check)."), style = "font-size: 13px; color: #666; margin-bottom: 4px;"), div(style = "height: 350px; overflow-y: auto;", DTOutput("drug_correlation_table"))))
              )
            )
          )
        ),
        tabPanel("Single Gene Associations",
          uiOutput("gene_summary_ui")
        )
      )
    } else if (identical(nav, "analyses")) {
      df <- filtered_data()
      genes <- filtered_genes()
      gene_choices_analyses <- if (length(genes) == 0) {
        c("Select Gene..." = "")
      } else {
        c("Select Gene..." = "", setNames(genes, genes))
      }
      current_analyses_tab <- selected_analyses_tab()
      if (is.null(current_analyses_tab) || !current_analyses_tab %in% c("Cohort Overview", "Single mutation associations", "Co-mutation Associations", "VAF Associations", "Drug Sensitivity")) current_analyses_tab <- "Cohort Overview"
      tabsetPanel(
        id = "main_tabs_analyses",
        selected = current_analyses_tab,
        tabPanel("Cohort Overview",
          fluidRow(
            column(3, wellPanel(h4("Selected Cohort"), tableOutput("summary_table"))),
            column(9, wellPanel(h4("Samples by Cohort & Subset"), plotOutput("cohort_plot", height = 280)))
          ),
          fluidRow(
            column(6, wellPanel(h4("Survival by Dataset"), plotOutput("overview_surv_dataset_plot", height = 400))),
            column(6, wellPanel(h4("Age Distribution by Dataset"), plotOutput("overview_age_dataset_plot", height = 400)))
          ),
          fluidRow(
            column(6, wellPanel(h4("ELN Risk by Dataset"), plotOutput("overview_eln_dataset_plot", height = 400))),
            column(6, wellPanel(h4("Recurrently Mutated Genes"), p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Genes mutated in more than 1% of samples."), plotOutput("overview_gene_freq_plot", height = 400)))
          ),
          fluidRow(
            column(12, wellPanel(
              h4("Oncoprint"),
              p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "All recurrent mutations in the selected cohort. Each column is a sample and genes are ordered by decreasing frequency."),
              plotOutput("oncoprint_plot", height = 500)
            ))
          )
        ),
        tabPanel("Single mutation associations",
          fluidRow(
            column(12, wellPanel(
              fluidRow(column(8, h4("Clinical Variable by Mutation Status")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_clinical_plot", "Download plot")))),
              selectInput("clin_var", "Variable:", choices = c("WBC", "Age", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent")),
              p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", span(style = "color: #8B0000; font-weight: bold;", "Red"), " = mutated; ", span(style = "color: #808080; font-weight: bold;", "Gray"), " = WT. *** p<0.001, ** p<0.01, * p<0.05, ns = not significant (Wilcoxon)."),
              plotOutput("clinical_plot", height = 400)))
          ),
          fluidRow(
            column(6, wellPanel(fluidRow(column(8, h4("Hazard Ratios")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_forest_plot", "Download plot")))), p(style = "font-size: 14px; color: #666; margin-bottom: 6px;", "Individual gene hazard ratios were calculated between mutated and wild-type patients using a univariate Cox proportional hazards model."), plotOutput("forest_plot", height = 500), p(style = "font-size: 13px; color: #666; margin-top: 6px;", span(style = "color: #762a83; font-weight: bold;", "Purple"), " = significant HR > 1 (higher risk); ", span(style = "color: #1b7837; font-weight: bold;", "Green"), " = significant HR < 1 (lower risk); ", span(style = "color: #9E9E9E; font-weight: bold;", "Gray"), " = not significant."))),
            column(6, wellPanel(fluidRow(column(8, h4("Kaplan-Meier")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_survival_plot", "Download plot")))), selectInput("surv_gene", NULL, choices = gene_choices_analyses), plotOutput("survival_plot", height = 500), p(style = "font-size: 13px; color: #666; margin-top: 6px;", span(style = "color: #8B0000; font-weight: bold;", "Red"), " = mutated; ", span(style = "color: #4D4D4D; font-weight: bold;", "Gray"), " = WT.")))
          ),
          fluidRow(column(12, wellPanel(fluidRow(column(8, h4("Survival Summary")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_survival_table", "Download table")))), p(style = "font-size: 13px; color: #666; margin-bottom: 8px;", "Cox proportional hazards (mutated vs WT) per gene; genes with ≥10 mutated and ≥10 WT. HR and 95% CI shown. ", em("p_adj"), " = false discovery rate (FDR, Benjamini\u2013Hochberg)."), DTOutput("survival_table"))))
        ),
        tabPanel("Co-mutation Associations",
          tabsetPanel(
            tabPanel("Oncoprint",
              fluidRow(
                column(12, wellPanel(
                  fluidRow(column(8, h4("Oncoprint")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_oncoprint_plot", "Download plot")))),
                  p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "All recurrent mutations in the selected cohort. Each column is a sample and genes are ordered by decreasing frequency."),
                  plotOutput("oncoprint_plot", height = 500)
                ))
              )
            ),
            tabPanel("Odds Ratio and Survival",
              fluidRow(
                column(6, wellPanel(fluidRow(column(8, h4("Co-mutation Odds Ratio")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_cooccurrence_plot", "Download plot")))), p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = co-occurring (OR > 1); ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = mutually exclusive (OR < 1); ", span(style = "color: #808080; font-weight: bold;", "Gray"), " = not significant."), div(style = "max-width: 100%; overflow: auto;", plotOutput("cooccurrence_plot", height = 500)))),
                column(6, wellPanel(fluidRow(column(8, h4("Mutation Co-occurrence Statistics")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_cooccurrence_table", "Download table")))), p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", "Odds ratio (OR) and FDR q-value. OR > 1 = genes tend to co-occur; OR < 1 = mutually exclusive."), div(style = "height: 500px; overflow-y: auto;", DTOutput("cooccurrence_table"))))
              ),
              fluidRow(
                column(12, wellPanel(
                  fluidRow(column(8, h4("Kaplan-Meier by Co-mutation Status")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_comut_survival_plot", "Download plot")))),
                  fluidRow(
                    column(4,
                      div(style = "background: #eaecef; border-radius: 8px; padding: 14px 16px; margin-top: 4px;",
                        h5(style = "margin-top: 0;", "Gene Selection"),
                        radioButtons("comut_n", "Number of genes:", choices = c("2 genes" = 2, "3 genes" = 3), selected = 2),
                        selectInput("comut_gene1", "Gene 1:", choices = gene_choices_analyses),
                        selectInput("comut_gene2", "Gene 2:", choices = gene_choices_analyses),
                        conditionalPanel(condition = "input.comut_n == 3",
                          selectInput("comut_gene3", "Gene 3:", choices = gene_choices_analyses))
                      ),
                      div(style = "margin-top: 12px;",
                        h5("Mutation Intersection (UpSet Plot)", style = "margin-top: 0;"),
                        p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", "Overlap of the selected genes across samples."),
                        plotOutput("upset_plot", height = 350),
                        div(style = "margin-top: 4px;", downloadButton("download_upset_plot", "Download plot"))
                      )
                    ),
                    column(8, uiOutput("comut_survival_plot_ui"))
                  )
                ))
              )
            )
          )
        ),
        tabPanel("Clonality",
          tabsetPanel(
            tabPanel("Single Gene",
              fluidRow(
                column(6, wellPanel(
                  fluidRow(column(8, h4("VAF / CCF by Gene")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_gene_plot", "Download plot")))),
                  plotOutput("vaf_gene_plot", height = 500))),
                column(6, wellPanel(
                  fluidRow(column(8, h4("VAF / CCF and Survival: Hazard Ratios (MaxStat)")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_survival_plot", "Download plot")))),
                  p(em("HR and 95% CI from Cox model using optimal MaxStat threshold per gene. Threshold computed from the selected allele metric (VAF or CCF). Genes with \u226520 mutated samples only.")),
                  plotOutput("vaf_survival_plot", height = 450)))
              ),
              fluidRow(
                column(6, wellPanel(
                  fluidRow(column(8, h4("VAF / CCF by Mutation Category")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_category_plot", "Download plot")))),
                  plotOutput("vaf_category_plot", height = 450))),
                column(6, wellPanel(
                  fluidRow(column(8, h4("Mutation HR vs MaxStat VAF/CCF HR")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_hr_vs_maxstat_plot", "Download plot")))),
                  p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Each point is a gene. Y-axis: mutation HR (mutated vs WT). X-axis: MaxStat VAF HR (high vs low VAF). Point size reflects optimal VAF cutpoint; color indicates significance (FDR q < 0.1)."),
                  plotOutput("hr_vs_maxstat_plot", height = 400)
                ))
              ),
              fluidRow(
                column(6, wellPanel(
                  fluidRow(column(8, h4("VAF / CCF Distribution by Dataset")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_cohort_plot", "Download plot")))),
                  plotOutput("vaf_cohort_plot", height = 450)))
              )
            ),
            tabPanel("Pairwise Gene",
              fluidRow(
                column(6, wellPanel(h4("VAF / CCF Scatter (Clonal Ordering)"),
                  p("Compare VAF of two genes. Points above line = Gene 1 before Gene 2."),
                  selectInput("vaf_gene1", "Gene 1:", choices = gene_choices_analyses), selectInput("vaf_gene2", "Gene 2:", choices = gene_choices_analyses),
                  plotOutput("vaf_scatter_plot", height = 450))),
                column(6, wellPanel(h4("Survival by Clonal Order"),
                  p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Kaplan-Meier comparing patients where Gene 1 has higher VAF (acquired first) vs Gene 2 has higher VAF (acquired first)."),
                  plotOutput("vaf_clonal_km_plot", height = 450)))
              ),
              fluidRow(
                column(6, wellPanel(
                  fluidRow(column(8, h4("Clonal Ordering Survival: All Gene Pairs")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_clonal_hr_summary_plot", "Download plot")))),
                  p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Cox HR (95% CI) comparing survival when Gene A is clonally dominant vs Gene B, for all co-mutated gene pairs with \u226510 non-ambiguous orderings per group. Point size reflects total non-ambiguous cases."),
                  plotOutput("clonal_hr_summary_plot", height = 600)
                )),
                column(6, wellPanel(
                  fluidRow(column(8, h4("VAF vs CCF")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_vaf_vs_ccf_plot", "Download plot")))),
                  p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", "Each point is a mutation. Diagonal = CCF equals VAF (diploid, CN = 2). Points above the diagonal have CN > 2 (gain); points below have CN < 2 (loss)."),
                  plotOutput("vaf_vs_ccf_plot", height = 550)
                ))
              )
            )
          )
        ),
        tabPanel("Drug Sensitivity",
          p("Data for mutation and VAF correlations with inhibitor sensitivity from BeatAML2 data.",
            a("biodev.github.io/BeatAML2", href = "https://biodev.github.io/BeatAML2/", target = "_blank")),
          fluidRow(
            column(12, wellPanel(
              div(id = "drug_subset_wrap", selectInput("drug_subset", "Select Beat AML cohort:", choices = c("All", "de_novo", "secondary"), selected = "de_novo")),
              uiOutput("drug_summary_ui")
            ))
          ),
          fluidRow(
            column(12, wellPanel(
              fluidRow(column(8, h4("Mut vs WT AUC (\u22655 mut samples)")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_mut_wt_heatmap", "Download plot")))),
              uiOutput("drug_mut_wt_summary_ui"),
              p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = mutated samples more resistant (higher AUC); ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = mutated samples more sensitive (lower AUC)."),
              p(em("Delta AUC = mean AUC (mut) - mean AUC (WT). Only gene-inhibitor pairs with q < 0.1 (FDR) shown."), style = "margin-bottom: 4px;"),
              div(style = "margin: 0; padding: 0; line-height: 0;", plotOutput("drug_mut_wt_heatmap", height = 500))
            ))
          ),
          fluidRow(
            column(12, wellPanel(
              fluidRow(column(8, h4("VAF vs AUC correlations")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_summary_dotplot", "Download plot")))),
              uiOutput("drug_vaf_auc_summary_ui"),
              p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = negative slope (higher VAF \u2192 lower AUC); ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = positive slope (higher VAF \u2192 higher AUC). * = q < 0.1 (FDR)."),
              plotOutput("drug_summary_dotplot", height = 500),
              p(em("* q < 0.1 (FDR)"), style = "font-size: 13px; color: #666;")
            ))
          ),
          fluidRow(
            column(6, wellPanel(
              fluidRow(column(8, h4("Individual Mutation and VAF Associations")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_scatter_plot", "Download plot")))),
              p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", "Boxplot: ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = mut more resistant; ", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = mut more sensitive. Scatter: ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = positive VAF\u2013AUC; ", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = negative."),
              p(em("Gene and inhibitor options match the subset used in Mut vs WT and VAF vs AUC analyses above.")),
              uiOutput("drug_scatter_selects_ui"),
              fluidRow(
                column(6, plotOutput("drug_scatter_boxplot", width = "100%", height = 320)),
                column(6, plotOutput("drug_scatter_plot", width = "100%", height = 320))
              )
            )),
            column(6, wellPanel(
              fluidRow(column(8, h4("Leave-One-Out Cross Validation for VAF vs AUC analysis")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_loo_heatmap", "Download plot")))),
              p(em("RMSE = leave-one-out prediction error in AUC units. Lower RMSE = better VAF\u2013AUC fit."), style = "font-size: 13px; color: #666; margin-bottom: 4px;"),
              plotOutput("drug_loo_heatmap", height = 350)))
          ),
          fluidRow(
            column(6, wellPanel(fluidRow(column(8, h4("Mut vs WT Statistics")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_mut_wt_table", "Download table")))), p(em("Delta AUC = mean AUC (mut) \u2212 mean AUC (WT). q < 0.1 (FDR) shown."), style = "font-size: 13px; color: #666; margin-bottom: 4px;"), div(style = "height: 350px; overflow-y: auto;", DTOutput("drug_mut_wt_summary_table")))),
            column(6, wellPanel(fluidRow(column(8, h4("VAF vs AUC Statistics")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_drug_correlation_table", "Download table")))), p(em("LOOCV RMSE = leave-one-out prediction error in AUC units (overfitting check)."), style = "font-size: 13px; color: #666; margin-bottom: 4px;"), div(style = "height: 350px; overflow-y: auto;", DTOutput("drug_correlation_table"))))
          )
        )
      )
    } else {
      div()
    }
  })

  # Remember selected tab when user switches (so filter changes don't reset to Overview)
  observeEvent(input$main_tabs_meta4, {
    if (!is.null(input$main_tabs_meta4)) selected_meta4_tab(input$main_tabs_meta4)
  }, ignoreNULL = TRUE)
  observeEvent(input$main_tabs_analyses, {
    if (!is.null(input$main_tabs_analyses)) selected_analyses_tab(input$main_tabs_analyses)
  }, ignoreNULL = TRUE)
  observeEvent(input$meta4_all_gene_tabs, {
    if (!is.null(input$meta4_all_gene_tabs)) selected_meta4_allgene_tab(input$meta4_all_gene_tabs)
  })
  observeEvent(input$gene_summary_sub_tabs, {
    if (!is.null(input$gene_summary_sub_tabs)) selected_gene_summary_sub_tab(input$gene_summary_sub_tabs)
  }, ignoreNULL = TRUE)

  # When user selects a gene from "Summarize all associations for a single gene", switch to Single Gene Associations tab (Meta AML4)
  observeEvent(input$gene_summary, {
    if (!is.null(input$gene_summary) && input$gene_summary != "" &&
        identical(input$main_nav, "meta_aml4")) {
      updateTabsetPanel(session, "main_tabs_meta4", selected = "Single Gene Associations")
    }
  }, ignoreNULL = FALSE)

  output$meta4_fallback_banner <- renderUI({
    if (identical(input$main_nav, "meta_aml4") && !meta4_uses_dedicated_file) {
      div(class = "alert alert-warning", role = "alert", style = "margin-bottom: 1em;",
        strong("Meta AML4 data file not found."), " Add ", code("AML_Meta_Cohort_v2.rds"), " or ", code("AML_Meta_Cohort.rds"), " (or ", code("AML_Meta_Cohort.RData"), ") to the app directory and redeploy."
      )
    }
  })

  observeEvent(input$main_nav, {
    nav <- input$main_nav
    if (identical(nav, "meta_aml4")) {
      updateSelectInput(session, "cohort", choices = c("All", sort(unique(as.character(na.omit(AML_Meta_Cohort$Cohort))))), selected = "All")
    } else if (identical(nav, "analyses")) {
      updateSelectInput(session, "cohort", choices = c("All", sort(unique(as.character(na.omit(final_data_matrix$Cohort))))), selected = "All")
    }
  }, ignoreNULL = TRUE, ignoreInit = FALSE)
  
  # Invalidate cache when filters change (min_gene_pct does NOT affect sample filtering, so not included here)
  observeEvent(c(input$subset, input$cohort, input$karyotype), {
    nav <- input$main_nav
    if (nav %in% c("analyses", "meta_aml4")) {
      cached_filtered_data[[nav]] <- NULL
      cached_data_params[[nav]] <- NULL
    }
  }, ignoreInit = TRUE)

  # Update slider max based on maximum gene frequency in filtered dataset (removed: min_gene_pct fixed at 1%)
  # observeEvent for slider removed; min gene frequency is 1% everywhere

  # When switching main tab, clear the *other* tab's filtered-data cache so only one full subset is in memory (reduces OOM)
  observeEvent(input$main_nav, {
    nav <- input$main_nav
    if (nav == "meta_aml4") {
      cached_filtered_data$analyses <- NULL
      cached_data_params$analyses <- NULL
    } else if (nav == "analyses") {
      cached_filtered_data$meta_aml4 <- NULL
      cached_data_params$meta_aml4 <- NULL
    }
  }, ignoreInit = TRUE)

  # When switching to Single Gene Drug Sensitivity sub-tab, run gc() to free memory before loading Drug content (reduces OOM)
  observeEvent(input$gene_summary_sub_tabs, {
    if (identical(input$gene_summary_sub_tabs, "drug")) gc()
  }, ignoreInit = TRUE)

  # Use effective cohort so that when switching tabs the first render shows correct data
  # (stale cohort value from other tab is treated as "All" until dropdown updates)
  effective_cohort <- reactive({
    req(current_data_src())
    cohorts_in_data <- sort(unique(as.character(na.omit(current_data_src()$Cohort))))
    co <- input$cohort
    if (is.null(co) || co == "" || !co %in% c("All", cohorts_in_data)) "All" else co
  })

  filtered_data <- reactive({
    nav <- req(input$main_nav)
    if (!nav %in% c("analyses", "meta_aml4")) return(NULL)
    req(current_data_src())
    
    # Get effective cohort for cache key (compute it here to avoid reactive dependency issues)
    cohorts_in_data <- sort(unique(as.character(na.omit(current_data_src()$Cohort))))
    co <- input$cohort
    eff_cohort <- if (is.null(co) || co == "" || !co %in% c("All", cohorts_in_data)) "All" else co
    
    # Check cache: if same tab and same filter params, return cached data
    # Note: gene frequency filter does NOT affect which samples are included, only which genes are shown
    current_params <- list(subset = input$subset, cohort = eff_cohort, karyotype = input$karyotype)
    if (!is.null(cached_filtered_data[[nav]]) && 
        identical(cached_data_params[[nav]], current_params)) {
      return(cached_filtered_data[[nav]])
    }
    
    # Compute filtered data (only filter by subset/cohort/karyotype - keep ALL samples)
    df <- prepare_data(current_data_src())
    if (!"Gene_for_analysis" %in% colnames(df)) df$Gene_for_analysis <- as.character(df$Gene)
    if (input$subset != "All") df <- df[as.character(df$Subset) == input$subset, , drop = FALSE]
    if (eff_cohort != "All") df <- df[as.character(df$Cohort) == eff_cohort, , drop = FALSE]
    if (!is.null(input$karyotype) && input$karyotype != "All" && "Karyotype" %in% colnames(df)) {
      df <- df[as.character(df$Karyotype) == input$karyotype, , drop = FALSE]
    }
    # Do NOT filter genes here - gene frequency filter only affects which genes are shown in dropdowns/plots, not which samples are included
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    
    # Cache the result
    cached_filtered_data[[nav]] <- df
    cached_data_params[[nav]] <- current_params
    
    df
  })

  # Tell client that data is ready (but wait for outputs to render before hiding overlay)
  # Only send dataReady for the currently active tab to prevent race conditions
  observe({
    nav <- req(input$main_nav)
    if (!nav %in% c("analyses", "meta_aml4")) return()
    req(filtered_data())
    # Ensure we're still on the same tab before sending ready signal
    if (identical(input$main_nav, nav)) {
      # Send dataReady signal - JavaScript will wait for outputs to render before hiding
      session$sendCustomMessage("dataReady", nav)
    }
  })

  # Get genes filtered by frequency (for display in dropdowns/plots only, not for sample filtering)
  filtered_genes <- reactive({
    req(df <- filtered_data())
    min_pct <- 1
    if (min_pct > 0) {
      n_samples <- length(unique(df$Sample))
      if (n_samples > 0) {
        gene_counts <- table(as.character(df$Gene_for_analysis))
        gene_pct <- (gene_counts / n_samples) * 100
        keep_genes <- names(gene_pct)[gene_pct >= min_pct]
        sort(unique(as.character(keep_genes)))
      } else {
        character(0)
      }
    } else {
      sort(unique(as.character(df$Gene_for_analysis)))
    }
  })

  observe({
    req(genes <- filtered_genes())
    updateSelectInput(session, "surv_gene", choices = c("Select Gene..." = "", genes))
    updateSelectInput(session, "vaf_gene1", choices = c("Select..." = "", genes))
    updateSelectInput(session, "vaf_gene2", choices = c("Select..." = "", genes))
    updateSelectInput(session, "comut_gene1", choices = c("Select..." = "", genes))
    updateSelectInput(session, "comut_gene2", choices = c("Select..." = "", genes))
    updateSelectInput(session, "comut_gene3", choices = c("Select..." = "", genes))
    if (identical(input$main_nav, "meta_aml4")) {
      updateSelectInput(session, "gene_summary", choices = c("Select gene..." = "", genes))
    }
  })

  output$summary_table <- renderTable({
    req(input$main_nav %in% c("analyses", "meta_aml4"))
    req(df <- filtered_data())
    samples_surv <- length(unique(df$Sample[!is.na(df$Time_to_OS) & !is.na(df$Censor)]))
    vals <- c(length(unique(df$Sample)), nrow(df), length(unique(df$Gene_for_analysis)), samples_surv)
    data.frame(
      Metric = c("Samples", "Mutations", "Genes", "With survival data"),
      Value = as.integer(vals)
    )
  }, colnames = FALSE, width = "100%")

  output$cohort_plot <- renderPlot({
    req(input$main_nav %in% c("analyses", "meta_aml4"))
    req(df <- filtered_data())
    udf <- unique(df[c("Sample", "Cohort", "Subset")])
    tbl <- table(udf$Cohort, udf$Subset)
    tbl_long <- as.data.frame(tbl)
    names(tbl_long) <- c("Cohort", "Subset", "n")
    cohort_totals <- aggregate(n ~ Cohort, data = tbl_long, sum)
    cohort_order <- as.character(cohort_totals$Cohort[order(cohort_totals$n, decreasing = TRUE)])
    tbl_long$Cohort <- factor(as.character(tbl_long$Cohort), levels = cohort_order)
    cohort_totals$Cohort <- factor(as.character(cohort_totals$Cohort), levels = cohort_order)
    subset_pal <- if (is_meta4()) META4_SUBSET_COL else PAL_SUBSET
    # Filter palette to only include subsets present in data
    present_subsets <- unique(tbl_long$Subset[!is.na(tbl_long$Subset)])
    subset_pal <- subset_pal[names(subset_pal) %in% present_subsets]
    ggplot(tbl_long, aes(x = Cohort, y = n, fill = Subset)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(data = cohort_totals, aes(x = Cohort, y = n, label = n), vjust = -0.5, size = 5, inherit.aes = FALSE) +
      scale_fill_manual(values = subset_pal, na.value = "gray70") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(x = NULL, y = "Number of samples") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.text = element_text(size = 13), legend.title = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # ---- Overview: Survival by Dataset (KM colored by dataset) ----
  output$overview_surv_dataset_plot <- renderPlot({
    req(df <- filtered_data())
    surv_df <- df[!is.na(df$Time_to_OS) & !is.na(df$Censor) & !is.na(df$Cohort), , drop = FALSE]
    surv_df <- surv_df[!duplicated(surv_df$Sample), , drop = FALSE]
    surv_df$Time_to_OS <- as.numeric(surv_df$Time_to_OS) / 365
    if (nrow(surv_df) < 5 || length(unique(surv_df$Cohort)) < 1) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient survival data") + theme_void())
    surv_df$Cohort <- factor(as.character(surv_df$Cohort))  # drop unused levels
    n_cohorts <- length(unique(surv_df$Cohort))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Cohort, data = surv_df)
    if (has_survminer) {
      cohort_pal <- if (is_meta4()) META4_COHORT_COL else PAL_COHORT
      present <- levels(surv_df$Cohort)
      legend_labs <- present
      pal_vec <- unname(cohort_pal[present])
      pal_vec[is.na(pal_vec)] <- "gray50"
      if (n_cohorts == 1) {
        p <- survminer::ggsurvplot(fit, data = surv_df, risk.table = TRUE, pval = FALSE,
          title = "Overall Survival by Dataset", xlab = "Years",
          palette = pal_vec, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.05), risk.table.height = 0.3)
      } else {
        p <- survminer::ggsurvplot(fit, data = surv_df, risk.table = TRUE, pval = TRUE,
          title = "Overall Survival by Dataset", xlab = "Years",
          palette = pal_vec, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.05),
          risk.table.height = 0.3)
      }
      print(p)
    } else {
      plot(fit, col = seq_along(fit$strata), xlab = "Years", ylab = "Survival", main = "Overall Survival by Dataset")
    }
  })

  # ---- Overview: Age Distribution by Dataset ----
  output$overview_age_dataset_plot <- renderPlot({
    req(df <- filtered_data())
    if (!"Age" %in% colnames(df)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No age data") + theme_void())
    udf <- df[!duplicated(df$Sample), , drop = FALSE]
    udf <- udf[!is.na(udf$Age) & !is.na(udf$Cohort), , drop = FALSE]
    if (nrow(udf) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No age data") + theme_void())
    cohort_pal <- if (is_meta4()) META4_COHORT_COL else PAL_COHORT
    present <- unique(udf$Cohort)
    cohort_pal <- cohort_pal[names(cohort_pal) %in% present]
    ggplot(udf, aes(x = Cohort, y = Age, fill = Cohort)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_fill_manual(values = cohort_pal, na.value = "gray70") +
      labs(x = NULL, y = "Age (years)") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15),
            legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # ---- Overview: ELN Risk by Dataset ----
  output$overview_eln_dataset_plot <- renderPlot({
    req(df <- filtered_data())
    if (!"Risk" %in% colnames(df)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No ELN risk data") + theme_void())
    udf <- df[!duplicated(df$Sample), , drop = FALSE]
    udf <- udf[!is.na(udf$Risk) & !is.na(udf$Cohort), , drop = FALSE]
    udf <- udf[udf$Risk != "Unknown", , drop = FALSE]
    if (nrow(udf) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No ELN risk data") + theme_void())
    udf$Risk <- factor(udf$Risk, levels = c("Favorable", "Intermediate", "Adverse"))
    tbl <- as.data.frame(table(Cohort = udf$Cohort, Risk = udf$Risk), stringsAsFactors = FALSE)
    cohort_totals <- aggregate(Freq ~ Cohort, data = tbl, sum)
    tbl <- merge(tbl, cohort_totals, by = "Cohort", suffixes = c("", "_total"))
    tbl$Pct <- ifelse(tbl$Freq_total > 0, tbl$Freq / tbl$Freq_total * 100, 0)
    tbl$Risk <- factor(tbl$Risk, levels = c("Favorable", "Intermediate", "Adverse"))
    ggplot(tbl, aes(x = Cohort, y = Pct, fill = Risk)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = PAL_RISK, na.value = "gray70") +
      labs(x = NULL, y = "Percentage of patients (%)", fill = "ELN Risk") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15),
            legend.text = element_text(size = 13), legend.title = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # ---- Overview: Recurrently Mutated Genes (>1% of samples) ----
  output$overview_gene_freq_plot <- renderPlot({
    req(df <- filtered_data())
    n_samples <- length(unique(df$Sample))
    if (n_samples == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void())
    gene_sample <- df[!duplicated(paste(df$Sample, df$Gene_for_analysis)), , drop = FALSE]
    gene_counts <- as.data.frame(table(Gene = gene_sample$Gene_for_analysis), stringsAsFactors = FALSE)
    gene_counts$Pct <- gene_counts$Freq / n_samples * 100
    gene_counts <- gene_counts[gene_counts$Pct >= 1, , drop = FALSE]
    if (nrow(gene_counts) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No genes mutated in >1% of samples") + theme_void())
    gene_counts <- gene_counts[order(gene_counts$Pct, decreasing = TRUE), , drop = FALSE]
    gene_counts$Gene <- factor(gene_counts$Gene, levels = rev(gene_counts$Gene))
    ggplot(gene_counts, aes(x = Pct, y = Gene)) +
      geom_col(fill = "#4292c6", width = 0.7) +
      geom_text(aes(label = paste0(round(Pct, 1), "%")), hjust = -0.1, size = 3.5) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "% of samples mutated", y = NULL) +
      theme_minimal(base_size = 14) +
      theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 12), axis.title = element_text(size = 13))
  })

  oncoprint_data <- reactive({
    req(input$main_nav %in% c("analyses", "meta_aml4"))
    req(df <- filtered_data())
    final_data_matrix_3 <- df
    if (!is_meta4()) {
      for (i in seq_len(nrow(final_data_matrix_3))) {
        if (identical(as.character(final_data_matrix_3$variant_type[i]), "other")) {
          final_data_matrix_3$variant_type[i] <- "Unknown"
        }
      }
    }
    # Group FLT3-ITD, FLT3-TKD, FLT3-Other together as FLT3, with labeled variant types (oncoprint only)
    flt3_itd_idx <- which(final_data_matrix_3$Gene == "FLT3-ITD")
    flt3_tkd_idx <- which(final_data_matrix_3$Gene == "FLT3-TKD")
    flt3_other_idx <- which(final_data_matrix_3$Gene == "FLT3-Other")
    if (length(flt3_itd_idx) > 0) {
      final_data_matrix_3$Gene[flt3_itd_idx] <- "FLT3"
      final_data_matrix_3$variant_type[flt3_itd_idx] <- "ITD"
    }
    if (length(flt3_tkd_idx) > 0) {
      final_data_matrix_3$Gene[flt3_tkd_idx] <- "FLT3"
      final_data_matrix_3$variant_type[flt3_tkd_idx] <- "SNV"
    }
    if (length(flt3_other_idx) > 0) {
      final_data_matrix_3$Gene[flt3_other_idx] <- "FLT3"
      final_data_matrix_3$variant_type[flt3_other_idx] <- "Other"
    }
    gene_tbl <- table(as.character(final_data_matrix_3$Gene))
    if (nrow(final_data_matrix_3) == 0) return(NULL)
    # Apply gene frequency filter only to determine which genes to show (all samples still included)
    min_pct <- 1
    if (min_pct > 0) {
      n_samples <- length(unique(final_data_matrix_3$Sample))
      if (n_samples > 0) {
        gene_pct <- (gene_tbl / n_samples) * 100
        keep_genes <- names(gene_pct)[gene_pct >= min_pct]
        gene_tbl <- gene_tbl[names(gene_tbl) %in% keep_genes]
      }
    }
    # Cap at 30 genes for oncoprint (from genes that pass frequency filter)
    sorted <- names(sort(gene_tbl, decreasing = TRUE))
    top_genes <- head(sorted, 30L)
    # Filter final_data_matrix_3 to only show mutations for top genes (but all samples are still represented in the matrix)
    final_data_matrix_3 <- final_data_matrix_3[as.character(final_data_matrix_3$Gene) %in% top_genes, , drop = FALSE]
    # Get all samples from original filtered data (not just those with mutations in top genes)
    all_samples <- unique(as.character(df$Sample))
    samples <- all_samples
    if (nrow(final_data_matrix_3) == 0) return(NULL)
    c <- length(samples)
    r <- length(top_genes)
    temp_dat <- as.data.frame(matrix(NA_character_, nrow = r, ncol = c))
    colnames(temp_dat) <- samples
    rownames(temp_dat) <- top_genes
    final_data_matrix_3$variant_type <- as.character(final_data_matrix_3$variant_type)
    # Aggregate multiple mutations per gene-sample pair (e.g., FLT3 ITD + TKD in same sample)
    for (i in seq_len(nrow(final_data_matrix_3))) {
      pt <- final_data_matrix_3$Sample[i]
      gene <- final_data_matrix_3$Gene[i]
      var_type <- final_data_matrix_3$variant_type[i]
      if (is.na(var_type) || var_type == "") var_type <- "Unknown"
      k <- match(pt, colnames(temp_dat))
      l <- match(gene, rownames(temp_dat))
      if (!is.na(k) && !is.na(l)) {
        existing <- temp_dat[l, k]
        if (is.na(existing) || existing == "") {
          temp_dat[l, k] <- var_type
        } else {
          # Combine multiple variant types with semicolon (ComplexHeatmap oncoPrint format)
          existing_types <- strsplit(existing, ";")[[1]]
          if (!var_type %in% existing_types) {
            temp_dat[l, k] <- paste(c(existing_types, var_type), collapse = ";")
          }
        }
      }
    }
    temp_dat[is.na(temp_dat)] <- ""
    temp_dat <- as.matrix(temp_dat)
    # Enforce max 30 genes: keep only first 30 rows
    if (nrow(temp_dat) > 30L) {
      top_genes <- head(top_genes, 30L)
      temp_dat <- temp_dat[top_genes, , drop = FALSE]
    }
    clin_cols <- c("Sample", "Cohort", "Sex", "Risk", "Subset", "Time_to_OS")
    clin_cols <- clin_cols[clin_cols %in% colnames(final_data_matrix_3)]
    if (length(clin_cols) < 2) return(list(mut = NULL, temp_dat = temp_dat, anno_df = NULL, samples = samples, genes = top_genes))
    anno_df <- unique(final_data_matrix_3[clin_cols])
    anno_df <- anno_df[!duplicated(anno_df$Sample), , drop = FALSE]
    colnames(anno_df)[colnames(anno_df) == "Time_to_OS"] <- "Survival"
    if ("Survival" %in% colnames(anno_df)) {
      for (i in seq_len(nrow(anno_df))) {
        if (!is.na(anno_df$Survival[i]) && anno_df$Survival[i] != "") {
          anno_df$Survival[i] <- "Yes"
        }
      }
      anno_df$Survival[is.na(anno_df$Survival) | anno_df$Survival == ""] <- "No"
    }
    anno_df <- anno_df[match(samples, anno_df$Sample), , drop = FALSE]
    list(mut = final_data_matrix_3[, c("Sample", "Gene", "variant_type"), drop = FALSE], temp_dat = temp_dat, anno_df = anno_df, samples = samples, genes = top_genes)
  })

  output$oncoprint_plot <- renderPlot({
    req(input$main_nav %in% c("analyses", "meta_aml4"))
    od <- oncoprint_data()
    if (is.null(od) || !is.matrix(od$temp_dat) || nrow(od$temp_dat) == 0 || ncol(od$temp_dat) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    use_meta4 <- is_meta4()
    if (has_ComplexHeatmap) {
      temp_dat <- od$temp_dat
      anno_df <- od$anno_df
      samples <- od$samples
      genes <- od$genes
      col <- if (use_meta4) META4_VAR_COL else ONCO_VAR_COL
      vtypes <- unique(as.character(temp_dat[temp_dat != ""]))
      for (vt in setdiff(vtypes, names(col))) col[vt] <- "#80796BFF"
      # Vectorized alter_fun for performance with >100 rows/columns (avoids cell_fun warning)
      alter_fun_list <- list(
        background = function(x, y, w, h) NULL,
        Deletion = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Deletion"], col = "#374E55FF")) },
        INDEL = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["INDEL"], col = "#DF8F44FF")) },
        Indel = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Indel"], col = "#DF8F44FF")) },
        Insertion = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Insertion"], col = "#00A1D5FF")) },
        ITD = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["ITD"], col = "#79AF97FF")) },
        SNV = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["SNV"], col = "#B24745FF")) },
        Splicing = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Splicing"], col = "#6A6599FF")) },
        PTD = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["PTD"], col = "tan")) },
        Other = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Other"], col = "#80796BFF")) },
        Unknown = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Unknown"], col = "#80796BFF")) }
      )
      for (vt in setdiff(vtypes, names(alter_fun_list))) {
        alter_fun_list[[vt]] <- (function(fc) function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = fc, col = fc)) })(col[vt])
      }
      if (is.null(anno_df) || nrow(anno_df) == 0) {
        ha <- NULL
      } else {
        Sex <- anno_df[, "Sex"]
        Cohort <- anno_df[, "Cohort"]
        Risk <- anno_df[, "Risk"]
        Subset <- anno_df[, "Subset"]
        Survival <- if ("Survival" %in% colnames(anno_df)) anno_df[, "Survival"] else rep("No", nrow(anno_df))
        Sex[is.na(Sex)] <- "Unknown"; Cohort[is.na(Cohort)] <- "Unknown"; Risk[is.na(Risk)] <- "Unknown"
        Subset[is.na(Subset)] <- "Unknown"
        cohort_col <- if (use_meta4) META4_COHORT_COL else ONCO_COHORT_COL
        # Build color mapping for cohorts present in data
        unique_cohorts <- unique(Cohort[!is.na(Cohort)])
        cohort_col_final <- character(length(unique_cohorts))
        names(cohort_col_final) <- unique_cohorts
        for (co in unique_cohorts) {
          if (co %in% names(cohort_col)) {
            cohort_col_final[co] <- cohort_col[co]
          } else {
            cohort_col_final[co] <- "gray70"
          }
        }
        subset_col <- if (use_meta4) META4_SUBSET_COL else ONCO_SUBSET_COL
        # Build color mapping for subsets present in data
        unique_subsets <- unique(Subset[!is.na(Subset)])
        subset_col_final <- character(length(unique_subsets))
        names(subset_col_final) <- unique_subsets
        for (sub in unique_subsets) {
          if (sub %in% names(subset_col)) {
            subset_col_final[sub] <- subset_col[sub]
          } else {
            subset_col_final[sub] <- "gray70"
          }
        }
        ha <- ComplexHeatmap::HeatmapAnnotation(
          Sex = Sex, Cohort = Cohort, Risk = Risk, Subset = Subset, Survival = Survival,
          col = list(
            Sex = c("Male" = "#6a51a3", "Female" = "#43a2ca", "Unknown" = "gray90"),
            Survival = c("Yes" = "#252525", "No" = "#f0f0f0"),
            Risk = c("Adverse" = "#E64B35FF", "Intermediate" = "#8491B4FF", "Favorable" = "#00A087FF", "Unknown" = "#767676FF"),
            Subset = subset_col_final,
            Cohort = cohort_col_final
          ),
          annotation_height = grid::unit(c(3.5, 3.5, 3.5, 3.5, 3.5), "mm"),
          show_annotation_name = TRUE,
          annotation_legend_param = list(
            Sex = list(title = "Sex"), Survival = list(title = "Survival"),
            Risk = list(title = "Risk"), Subset = list(title = "Subset"), Cohort = list(title = "Cohort")
          )
        )
      }
      fig1B <- ComplexHeatmap::oncoPrint(temp_dat,
        col = col,
        row_names_side = "right",
        bottom_annotation = ha,
        alter_fun_is_vectorized = TRUE,
        alter_fun = alter_fun_list
      )
      ComplexHeatmap::draw(fig1B)
    } else {
      mut <- od$mut; clin <- od$anno_df; samples <- od$samples; genes <- od$genes
      mut$Gene <- factor(as.character(mut$Gene), levels = rev(genes))
      mut$Sample <- factor(as.character(mut$Sample), levels = samples)
      mut$Sample_num <- as.numeric(mut$Sample)
      mut$Gene_num <- as.numeric(mut$Gene)
      PAL_VAR <- if (use_meta4) c(Deletion = "#374E55", Indel = "#DF8F44", Insertion = "#00A1D5", ITD = "#79AF97",
        SNV = "#B24745", Splicing = "#6A6599", PTD = "tan", Other = "#80796B", Unknown = "#80796B") else
        c(Deletion = "#374E55", INDEL = "#DF8F44", Insertion = "#00A1D5", ITD = "#79AF97", SNV = "#B24745", Splicing = "#6A6599", Unknown = "#80796B")
      for (vt in setdiff(unique(mut$variant_type), names(PAL_VAR))) PAL_VAR[vt] <- "#80796B"
      bg <- expand.grid(Sample = factor(samples, levels = samples), Gene = factor(genes, levels = rev(genes)))
      bg$Sample_num <- as.numeric(bg$Sample)
      bg$Gene_num <- as.numeric(bg$Gene)
      p_mut <- ggplot() +
        geom_tile(data = bg, aes(x = Sample_num, y = Gene_num), fill = "gray97", width = 1, height = 1) +
        geom_rect(data = mut, aes(xmin = Sample_num - 0.25, xmax = Sample_num + 0.25, ymin = Gene_num - 0.45, ymax = Gene_num + 0.45, fill = variant_type)) +
        scale_fill_manual(values = PAL_VAR, name = "Variant") +
        scale_x_continuous(limits = c(0.5, length(samples) + 0.5), expand = c(0, 0), position = "top") +
        scale_y_continuous(breaks = seq_along(genes), labels = rev(genes), expand = c(0, 0)) +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "right",
          plot.margin = margin(5, 5, 2, 5))
      if (nrow(clin) > 0 && "Survival" %in% colnames(clin)) {
        annot_cols <- c("Sex", "Cohort", "Risk", "Subset", "Survival")
        annot_cols <- annot_cols[annot_cols %in% colnames(clin)]
        clin_long <- do.call(rbind, lapply(annot_cols, function(col) {
          data.frame(Sample = clin$Sample, annot = col, value = as.character(clin[[col]]), stringsAsFactors = FALSE)
        }))
        clin_long$Sample <- factor(as.character(clin_long$Sample), levels = samples)
        clin_long$Sample_num <- as.numeric(clin_long$Sample)
        clin_long$annot <- factor(clin_long$annot, levels = rev(annot_cols))
        clin_long$value[is.na(clin_long$value) | clin_long$value == ""] <- "NA"
        pal_clin <- if (use_meta4) c(META4_SUBSET_COL, META4_COHORT_COL, PAL_RISK, Male = "#6a51a3", Female = "#43a2ca", Yes = "#252525", No = "#f0f0f0") else c(PAL_SUBSET, PAL_COHORT, PAL_RISK, Male = "#6a51a3", Female = "#43a2ca", Yes = "#252525", No = "#f0f0f0")
        p_clin <- ggplot(clin_long, aes(x = Sample_num, y = annot, fill = value)) +
          geom_tile(color = "white", linewidth = 0.2) +
          scale_fill_manual(values = pal_clin, na.value = "gray90", name = "Clinical") +
          scale_x_continuous(limits = c(0.5, length(samples) + 0.5), expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +
          labs(x = NULL, y = NULL) +
          theme_minimal(base_size = 9) +
          theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "right",
            plot.margin = margin(2, 5, 5, 5))
        if (has_gridExtra) {
          gridExtra::grid.arrange(p_mut, p_clin, ncol = 1, heights = c(3, 1))
        } else {
          print(p_mut)
        }
      } else {
        print(p_mut)
      }
    }
  })

  # Reactive helpers for VAF/CCF metric toggle
  vaf_col <- reactive({ if (!is.null(input$vaf_metric) && input$vaf_metric == "CCF") "CCF" else "VAF" })
  vaf_label <- reactive({ if (!is.null(input$vaf_metric) && input$vaf_metric == "CCF") "Cancer Cell Fraction (%)" else "Variant Allele Frequency (%)" })

  output$vaf_gene_plot <- renderPlot({
    req(df <- filtered_data())
    vcol <- vaf_col()
    vlabel <- vaf_label()
    if (!vcol %in% colnames(df)) vcol <- "VAF"
    df <- df[!is.na(df[[vcol]]) & df[[vcol]] > 0, , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    gene_tbl <- table(as.character(df$Gene_for_analysis))
    min_pct <- 1
    if (min_pct > 0) {
      n_samples <- length(unique(df$Sample))
      if (n_samples > 0) {
        gene_pct <- (gene_tbl / n_samples) * 100
        keep_genes <- names(gene_pct)[gene_pct >= min_pct]
        gene_tbl <- gene_tbl[names(gene_tbl) %in% keep_genes]
      }
    }
    if (length(gene_tbl) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No genes pass the Min. gene frequency filter") + theme_void())
    plot_genes <- names(sort(gene_tbl, decreasing = TRUE))
    df <- df[as.character(df$Gene_for_analysis) %in% plot_genes, , drop = FALSE]
    df$.metric <- df[[vcol]]
    med <- aggregate(.metric ~ Gene_for_analysis, data = df, FUN = median)
    gene_order <- as.character(med$Gene_for_analysis[order(med$.metric)])
    df$Gene_for_analysis <- factor(as.character(df$Gene_for_analysis), levels = gene_order)
    n_cat <- if ("mutation_category" %in% colnames(df)) length(unique(df$mutation_category[!is.na(df$mutation_category) & df$mutation_category != ""])) else 0L
    fill_var <- if (n_cat >= 2L) "mutation_category" else "Gene_for_analysis"
    scale_y_shared <- scale_y_discrete(limits = gene_order, drop = FALSE)
    p_vaf <- ggplot(df, aes(x = .metric, y = Gene_for_analysis, fill = .data[[fill_var]])) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_x_continuous(limits = c(0, 100)) +
      scale_y_shared +
      (if (fill_var == "mutation_category") scale_fill_manual(name = "Mutation category", values = PAL_MUT_CAT, na.value = "gray70", guide = "legend") else scale_fill_discrete(guide = "none")) +
      labs(x = vlabel, y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.text = element_text(size = 13), legend.title = element_text(size = 14))
    # Right: N mutations per gene – horizontal bars, same y scale; colored by mutation_category when used
    if (fill_var == "mutation_category" && "mutation_category" %in% colnames(df)) {
      cat_counts <- as.data.frame(table(Gene_for_analysis = df$Gene_for_analysis, mutation_category = df$mutation_category), stringsAsFactors = FALSE)
      cat_counts$Gene_for_analysis <- factor(as.character(cat_counts$Gene_for_analysis), levels = gene_order)
      cat_counts$mutation_category[is.na(cat_counts$mutation_category) | cat_counts$mutation_category == ""] <- NA_character_
      p_count <- ggplot(cat_counts, aes(x = Freq, y = Gene_for_analysis, fill = mutation_category)) +
        geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
        geom_col(position = "stack", width = 0.7) +
        scale_y_shared +
        scale_x_continuous(breaks = function(x) { m <- max(x, na.rm = TRUE); if (m <= 0) 0 else round(c(0, m / 2, m)) }) +
        scale_fill_manual(values = PAL_MUT_CAT, na.value = "gray70", guide = "none") +
        labs(x = "# of mutations", y = NULL) +
        theme_minimal(base_size = 16) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_text(size = 15), axis.text = element_text(size = 14), plot.margin = margin(l = 6, r = 4, unit = "pt"))
    } else {
      count_df <- data.frame(Gene_for_analysis = factor(gene_order, levels = gene_order), n = as.integer(gene_tbl[gene_order]))
      count_df$n[is.na(count_df$n)] <- 0L
      p_count <- ggplot(count_df, aes(x = n, y = Gene_for_analysis)) +
        geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
        scale_y_shared +
        scale_x_continuous(breaks = function(x) { m <- max(x, na.rm = TRUE); if (m <= 0) 0 else round(c(0, m / 2, m)) }) +
        geom_col(width = 0.7, fill = "gray50") +
        labs(x = "# of mutations", y = NULL) +
        theme_minimal(base_size = 16) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_text(size = 15), axis.text = element_text(size = 14), plot.margin = margin(l = 6, r = 4, unit = "pt"))
    }
    if (has_patchwork) {
      # Patchwork: align plots by shared y-axis, legend collected and placed on top
      combined <- p_vaf + p_count +
        patchwork::plot_layout(ncol = 2, widths = c(2, 0.6), guides = "collect") +
        theme(legend.position = "top")
      print(combined)
    } else if (has_gridExtra) {
      p_vaf_noleg <- p_vaf + theme(legend.position = "none")
      leg_grob <- NULL
      if (fill_var == "mutation_category") {
        g <- ggplotGrob(p_vaf)
        idx <- which(g$layout$name == "guide-box")
        if (length(idx) > 0) leg_grob <- g$grobs[[idx[1]]]
      }
      if (!is.null(leg_grob)) {
        gridExtra::grid.arrange(leg_grob, gridExtra::arrangeGrob(p_vaf_noleg, p_count, ncol = 2, widths = c(2, 0.6)), nrow = 2, heights = c(0.12, 1))
      } else {
        gridExtra::grid.arrange(p_vaf_noleg, p_count, ncol = 2, widths = c(2.2, 0.6))
      }
    } else {
      print(p_vaf + theme(legend.position = "right"))
    }
  })

  output$vaf_cohort_plot <- renderPlot({
    req(df <- filtered_data())
    vcol <- vaf_col(); vlabel <- vaf_label()
    if (!vcol %in% colnames(df)) vcol <- "VAF"
    df <- df[!is.na(df[[vcol]]) & df[[vcol]] > 0 & !is.na(df$Cohort), , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    df$.metric <- df[[vcol]]
    cohort_pal <- if (is_meta4()) META4_COHORT_COL else PAL_COHORT
    present_cohorts <- unique(df$Cohort[!is.na(df$Cohort)])
    cohort_pal <- cohort_pal[names(cohort_pal) %in% present_cohorts]
    med <- aggregate(.metric ~ Cohort, data = df, FUN = median)
    cohort_order <- as.character(med$Cohort[order(med$.metric)])
    df$Cohort <- factor(as.character(df$Cohort), levels = cohort_order)
    scale_y_coh <- scale_y_discrete(limits = cohort_order, drop = FALSE)
    p_box <- ggplot(df, aes(x = .metric, y = Cohort, fill = Cohort)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_x_continuous(limits = c(0, 100)) +
      scale_y_coh +
      scale_fill_manual(values = cohort_pal, na.value = "gray70") +
      labs(x = vlabel, y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.position = "none")
    coh_counts <- as.data.frame(table(Cohort = df$Cohort), stringsAsFactors = FALSE)
    coh_counts$Cohort <- factor(as.character(coh_counts$Cohort), levels = cohort_order)
    p_count <- ggplot(coh_counts, aes(x = Freq, y = Cohort, fill = Cohort)) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      geom_col(width = 0.7) +
      scale_y_coh +
      scale_x_continuous(breaks = function(x) { m <- max(x, na.rm = TRUE); if (m <= 0) 0 else round(c(0, m / 2, m)) }) +
      scale_fill_manual(values = cohort_pal, na.value = "gray70", guide = "none") +
      labs(x = "# of mutations", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_text(size = 15), axis.text = element_text(size = 14), plot.margin = margin(l = 6, r = 4, unit = "pt"))
    if (has_patchwork) {
      combined <- p_box + p_count + patchwork::plot_layout(ncol = 2, widths = c(2, 0.6))
      print(combined)
    } else {
      print(p_box)
    }
  })

  output$vaf_category_plot <- renderPlot({
    req(df <- filtered_data())
    vcol <- vaf_col(); vlabel <- vaf_label()
    if (!vcol %in% colnames(df)) vcol <- "VAF"
    df <- df[!is.na(df[[vcol]]) & df[[vcol]] > 0 & !is.na(df$mutation_category) & df$mutation_category != "", , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    df$.metric <- df[[vcol]]
    med <- aggregate(.metric ~ mutation_category, data = df, FUN = median)
    cat_order <- as.character(med$mutation_category[order(med$.metric)])
    df$mutation_category <- factor(as.character(df$mutation_category), levels = cat_order)
    scale_y_cat <- scale_y_discrete(limits = cat_order, drop = FALSE)
    p_box <- ggplot(df, aes(x = .metric, y = mutation_category, fill = mutation_category)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_x_continuous(limits = c(0, 100)) +
      scale_y_cat +
      scale_fill_manual(values = PAL_MUT_CAT, na.value = "gray70") +
      labs(x = vlabel, y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.position = "none")
    cat_counts <- as.data.frame(table(mutation_category = df$mutation_category), stringsAsFactors = FALSE)
    cat_counts$mutation_category <- factor(as.character(cat_counts$mutation_category), levels = cat_order)
    p_count <- ggplot(cat_counts, aes(x = Freq, y = mutation_category, fill = mutation_category)) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
      geom_col(width = 0.7) +
      scale_y_cat +
      scale_x_continuous(breaks = function(x) { m <- max(x, na.rm = TRUE); if (m <= 0) 0 else round(c(0, m / 2, m)) }) +
      scale_fill_manual(values = PAL_MUT_CAT, na.value = "gray70", guide = "none") +
      labs(x = "# of mutations", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_text(size = 15), axis.text = element_text(size = 14), plot.margin = margin(l = 6, r = 4, unit = "pt"))
    if (has_patchwork) {
      combined <- p_box + p_count + patchwork::plot_layout(ncol = 2, widths = c(2, 0.6))
      print(combined)
    } else {
      print(p_box)
    }
  })

  vaf_survival_hr_summary <- reactive({
    req(df <- filtered_data())
    vcol <- vaf_col()
    if (!vcol %in% colnames(df)) vcol <- "VAF"
    if (!vcol %in% colnames(df)) return(NULL)
    surv_df <- survival_data()
    if (is.null(surv_df) || nrow(surv_df) == 0) return(NULL)
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    genes <- filtered_genes()
    if (length(genes) == 0) return(NULL)
    if (!has_maxstat) return(NULL)
    res_list <- list()
    for (gene in genes) {
      sub <- df[as.character(df$Gene_for_analysis) == gene & !is.na(df[[vcol]]) & df[[vcol]] > 0, c("Sample", vcol), drop = FALSE]
      if (nrow(sub) == 0) next
      colnames(sub)[colnames(sub) == vcol] <- ".metric"
      vaf_agg <- aggregate(.metric ~ Sample, data = sub, FUN = max)
      merge_df <- merge(vaf_agg, surv_uniq, by = "Sample", all = FALSE)
      if (nrow(merge_df) < 20) next
      mst <- tryCatch(maxstat::maxstat.test(Surv(Time_to_OS, as.numeric(Censor)) ~ .metric, data = merge_df, smethod = "LogRank", minprop = 0.2, maxprop = 0.8), error = function(e) NULL)
      if (is.null(mst)) next
      cutpoint <- mst$estimate
      merge_df$Group <- factor(ifelse(merge_df$.metric > cutpoint, "High", "Low"), levels = c("Low", "High"))
      if (length(unique(merge_df$Group)) < 2) next
      cx <- tryCatch(survival::coxph(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = merge_df), error = function(e) NULL)
      if (is.null(cx)) next
      ci <- tryCatch(exp(confint(cx)), error = function(e) NULL)
      if (is.null(ci)) next
      cx_sum <- summary(cx)
      hr <- exp(coef(cx))["GroupHigh"]
      ci_lo <- ci["GroupHigh", 1]
      ci_hi <- ci["GroupHigh", 2]
      pval <- cx_sum$coefficients["GroupHigh", "Pr(>|z|)"]
      n_high <- sum(merge_df$Group == "High")
      res_list[[length(res_list) + 1L]] <- data.frame(Gene = gene, HR = as.numeric(hr), CI_lower = as.numeric(ci_lo), CI_upper = as.numeric(ci_hi), p_value = as.numeric(pval), threshold = round(cutpoint, 1), n = nrow(merge_df), n_high = n_high, stringsAsFactors = FALSE)
    }
    if (length(res_list) == 0) return(NULL)
    out <- do.call(rbind, res_list)
    out <- out[order(out$HR), , drop = FALSE]
    out
  })

  output$vaf_survival_plot <- renderPlot({
    hr_summary <- vaf_survival_hr_summary()
    vcol <- vaf_col(); vlabel <- vaf_label()
    metric_short <- if (vcol == "CCF") "CCF" else "VAF"
    if (is.null(hr_summary) || nrow(hr_summary) == 0) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste0("No genes with \u226520 mutated samples (", metric_short, " + survival). Install maxstat for MaxStat threshold.")) + theme_void())
    }
    if (!has_maxstat) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste0("Install the maxstat package for optimal ", metric_short, " threshold.")) + theme_void())
    }
    hr_summary$Gene <- factor(hr_summary$Gene, levels = hr_summary$Gene)
    hr_summary$ColorGroup <- ifelse(hr_summary$p_value >= 0.05, "Not significant",
      ifelse(hr_summary$HR > 1, "Significant (HR > 1)", "Significant (HR < 1)"))
    ggplot(hr_summary, aes(x = HR, y = Gene, xmin = CI_lower, xmax = CI_upper, color = ColorGroup)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      geom_errorbarh(height = 0.2, linewidth = 0.6) +
      geom_point(aes(size = n_high), alpha = 0.8) +
      scale_color_manual(values = c("Significant (HR > 1)" = "#762a83", "Significant (HR < 1)" = "#1b7837", "Not significant" = "#9E9E9E"), name = NULL) +
      scale_size_continuous(name = paste0("# above\n", metric_short, "\nthreshold"), range = c(2, 7)) +
      scale_x_continuous(trans = "log10", name = paste0("Hazard ratio (High vs Low ", metric_short, ", 95% CI)")) +
      labs(y = NULL, title = paste0(metric_short, "\u2013Survival: HR per gene (MaxStat threshold)")) +
      theme_minimal(base_size = 14) +
      theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 13), plot.title = element_text(size = 14, face = "bold"))
  })

  survival_data <- reactive({
    df <- filtered_data()
    df <- df[!is.na(df$Time_to_OS) & !is.na(df$Censor), , drop = FALSE]
    df$Time_to_OS <- as.numeric(df$Time_to_OS) / 365
    df
  })

  # Minimal survival table for Single Gene KM only (reduces memory: one row per sample, 4 columns)
  gene_summary_survival_km_data <- reactive({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(NULL)
    surv_df <- survival_data()
    if (is.null(surv_df) || nrow(surv_df) == 0) return(NULL)
    mut_samples <- unique(surv_df$Sample[as.character(surv_df$Gene_for_analysis) == g])
    out <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor"), drop = FALSE]
    out$Mutation <- ifelse(out$Sample %in% mut_samples, paste0(g, " mut"), "WT")
    out
  })

  output$survival_plot <- renderPlot({
    gene <- input$surv_gene
    if (is.null(gene) || gene == "") return(
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene", size = 10) +
        theme_void() +
        theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA))
    )
    surv_df <- survival_data()
    mut_samples <- unique(surv_df$Sample[surv_df$Gene_for_analysis == gene])
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    surv_uniq$Mutation <- ifelse(surv_uniq$Sample %in% mut_samples, paste0(gene, " mut"), "WT")
    if (length(unique(surv_uniq$Mutation)) < 2) return(
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient groups", size = 10) +
        theme_void() +
        theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA))
    )
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Mutation, data = surv_uniq)
    strata_names <- names(fit$strata)
    is_mut <- grepl(" mut", strata_names, fixed = TRUE)
    surv_pal <- ifelse(is_mut, "#8B0000", "#4D4D4D")  # deep red for mutated, gray for WT
    if (has_survminer) {
      legend_labs <- gsub("^Mutation=", "", names(fit$strata))
      p <- survminer::ggsurvplot(fit, data = surv_uniq, risk.table = TRUE, pval = TRUE,
        title = paste("Survival by", gene, "mutation status"), xlab = "Years",
        palette = surv_pal, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.05))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install the survminer package for Kaplan-Meier plots with number at risk.", size = 4) + theme_void()
    }
  })

  # ---- Scatterplot: Mutation HR vs MaxStat VAF HR ----
  hr_vs_maxstat_data <- reactive({
    mut_hr <- forest_data()
    vaf_hr <- vaf_survival_hr_summary()
    if (is.null(mut_hr) || nrow(mut_hr) == 0 || is.null(vaf_hr) || nrow(vaf_hr) == 0) return(NULL)
    mut_hr <- mut_hr[, c("Gene", "HR", "p"), drop = FALSE]
    colnames(mut_hr) <- c("Gene", "Mutation_HR", "Mutation_p")
    vaf_hr <- vaf_hr[, c("Gene", "HR", "p_value", "threshold"), drop = FALSE]
    colnames(vaf_hr) <- c("Gene", "MaxStat_HR", "MaxStat_p", "VAF_threshold")
    merged <- merge(mut_hr, vaf_hr, by = "Gene")
    if (nrow(merged) == 0) return(NULL)
    merged$Mutation_q <- p.adjust(merged$Mutation_p, method = "fdr")
    merged$MaxStat_q <- p.adjust(merged$MaxStat_p, method = "fdr")
    merged$ColorGroup <- ifelse(merged$MaxStat_q >= 0.1, "Not significant",
      ifelse(merged$MaxStat_HR > 1, "Significant (HR > 1)", "Significant (HR < 1)"))
    merged$label_rank <- rank(merged$MaxStat_q, ties.method = "first")
    merged
  })

  output$hr_vs_maxstat_plot <- renderPlot({
    df <- hr_vs_maxstat_data()
    if (is.null(df) || nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") + theme_void())
    metric_short <- if (vaf_col() == "CCF") "CCF" else "VAF"
    sig_df <- df[df$ColorGroup != "Not significant", , drop = FALSE]
    p <- ggplot(df, aes(x = MaxStat_HR, y = Mutation_HR, size = VAF_threshold, color = ColorGroup)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("Significant (HR > 1)" = "#762a83", "Significant (HR < 1)" = "#1b7837", "Not significant" = "#9E9E9E"), name = NULL) +
      scale_size_continuous(name = paste(metric_short, "threshold"), range = c(2, 8)) +
      labs(x = paste0("MaxStat ", metric_short, " HR (High vs Low ", metric_short, ")"), y = "Mutation HR (Mutated vs WT)") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15),
            legend.text = element_text(size = 13), legend.title = element_text(size = 14))
    if (nrow(sig_df) > 0) {
      if (has_ggrepel) {
        p <- p + ggrepel::geom_label_repel(data = sig_df, aes(label = Gene), size = 3.5, max.overlaps = 20, show.legend = FALSE, fill = "white", alpha = 0.85, label.size = 0.2)
      } else {
        p <- p + geom_text(data = sig_df, aes(label = Gene), size = 3.5, vjust = -0.8, show.legend = FALSE)
      }
    }
    p
  })

  output$download_hr_vs_maxstat_plot <- downloadHandler(
    filename = function() paste0("mutation_hr_vs_maxstat_hr_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      df <- hr_vs_maxstat_data()
      if (is.null(df) || nrow(df) == 0) {
        ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 10, height = 7, dpi = 150)
        return(invisible(NULL))
      }
      metric_short <- if (vaf_col() == "CCF") "CCF" else "VAF"
      sig_df <- df[df$ColorGroup != "Not significant", , drop = FALSE]
      p <- ggplot(df, aes(x = MaxStat_HR, y = Mutation_HR, size = VAF_threshold, color = ColorGroup)) +
        geom_point(alpha = 0.8) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
        scale_color_manual(values = c("Significant (HR > 1)" = "#762a83", "Significant (HR < 1)" = "#1b7837", "Not significant" = "#9E9E9E"), name = NULL) +
        scale_size_continuous(name = paste(metric_short, "threshold"), range = c(2, 8)) +
        labs(x = paste0("MaxStat ", metric_short, " HR (High vs Low ", metric_short, ")"), y = "Mutation HR (Mutated vs WT)") +
        theme_minimal(base_size = 16) +
        theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15),
              legend.text = element_text(size = 13), legend.title = element_text(size = 14))
      if (nrow(sig_df) > 0) {
        if (has_ggrepel) {
          p <- p + ggrepel::geom_label_repel(data = sig_df, aes(label = Gene), size = 3.5, max.overlaps = 20, show.legend = FALSE, fill = "white", alpha = 0.85, label.size = 0.2)
        } else {
          p <- p + geom_text(data = sig_df, aes(label = Gene), size = 3.5, vjust = -0.8, show.legend = FALSE)
        }
      }
      ggsave(file, plot = p, width = 10, height = 7, dpi = 150)
    }
  )

  output$download_vaf_gene_plot <- downloadHandler(
    filename = function() paste0(tolower(vaf_col()), "_by_gene_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      df <- filtered_data(); vcol <- vaf_col(); vlabel <- vaf_label()
      if (is.null(df)) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 10, height = 7, dpi = 150); return(invisible(NULL)) }
      if (!vcol %in% colnames(df)) vcol <- "VAF"
      df <- df[!is.na(df[[vcol]]) & df[[vcol]] > 0, , drop = FALSE]
      if (nrow(df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 10, height = 7, dpi = 150); return(invisible(NULL)) }
      gene_tbl <- table(as.character(df$Gene_for_analysis))
      min_pct <- 1
      if (min_pct > 0) { n_samples <- length(unique(df$Sample)); if (n_samples > 0) { gene_pct <- (gene_tbl / n_samples) * 100; gene_tbl <- gene_tbl[names(gene_tbl) %in% names(gene_pct)[gene_pct >= min_pct]] } }
      if (length(gene_tbl) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No genes") + theme_void(), width = 10, height = 7, dpi = 150); return(invisible(NULL)) }
      plot_genes <- names(sort(gene_tbl, decreasing = TRUE))
      df <- df[as.character(df$Gene_for_analysis) %in% plot_genes, , drop = FALSE]
      df$.metric <- df[[vcol]]
      med <- aggregate(.metric ~ Gene_for_analysis, data = df, FUN = median)
      gene_order <- as.character(med$Gene_for_analysis[order(med$.metric)])
      df$Gene_for_analysis <- factor(as.character(df$Gene_for_analysis), levels = gene_order)
      n_cat <- if ("mutation_category" %in% colnames(df)) length(unique(df$mutation_category[!is.na(df$mutation_category) & df$mutation_category != ""])) else 0L
      fill_var <- if (n_cat >= 2L) "mutation_category" else "Gene_for_analysis"
      p <- ggplot(df, aes(x = .metric, y = Gene_for_analysis, fill = .data[[fill_var]])) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) + scale_x_continuous(limits = c(0, 100)) +
        (if (fill_var == "mutation_category") scale_fill_manual(name = "Mutation category", values = PAL_MUT_CAT, na.value = "gray70") else scale_fill_discrete(guide = "none")) +
        labs(x = vlabel, y = NULL) + theme_minimal(base_size = 14) +
        theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.position = "top")
      ggsave(file, plot = p, width = 10, height = max(5, 0.3 * length(plot_genes) + 2), dpi = 300)
    }
  )

  output$download_vaf_survival_plot <- downloadHandler(
    filename = function() paste0(tolower(vaf_col()), "_survival_hr_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      hr_summary <- vaf_survival_hr_summary()
      vcol <- vaf_col(); metric_short <- if (vcol == "CCF") "CCF" else "VAF"
      if (is.null(hr_summary) || nrow(hr_summary) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 6, dpi = 150); return(invisible(NULL)) }
      hr_summary$Gene <- factor(hr_summary$Gene, levels = hr_summary$Gene)
      hr_summary$ColorGroup <- ifelse(hr_summary$p_value >= 0.05, "Not significant", ifelse(hr_summary$HR > 1, "Significant (HR > 1)", "Significant (HR < 1)"))
      p <- ggplot(hr_summary, aes(x = HR, y = Gene, xmin = CI_lower, xmax = CI_upper, color = ColorGroup)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
        geom_errorbarh(height = 0.2, linewidth = 0.6) +
        geom_point(aes(size = n_high), alpha = 0.8) +
        scale_color_manual(values = c("Significant (HR > 1)" = "#762a83", "Significant (HR < 1)" = "#1b7837", "Not significant" = "#9E9E9E"), name = NULL) +
        scale_size_continuous(name = paste0("# above\n", metric_short, "\nthreshold"), range = c(2, 7)) +
        scale_x_continuous(trans = "log10", name = paste0("Hazard ratio (High vs Low ", metric_short, ", 95% CI)")) +
        labs(y = NULL, title = paste0(metric_short, "\u2013Survival: HR per gene (MaxStat threshold)")) +
        theme_minimal(base_size = 14) +
        theme(axis.text.y = element_text(size = 12), axis.title = element_text(size = 13), plot.title = element_text(size = 14, face = "bold"))
      ggsave(file, plot = p, width = 8, height = max(4, 0.3 * nrow(hr_summary) + 2), dpi = 300)
    }
  )

  output$download_vaf_category_plot <- downloadHandler(
    filename = function() paste0(tolower(vaf_col()), "_by_category_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      df <- filtered_data(); vcol <- vaf_col(); vlabel <- vaf_label()
      if (is.null(df)) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 6, dpi = 150); return(invisible(NULL)) }
      if (!vcol %in% colnames(df)) vcol <- "VAF"
      df <- df[!is.na(df[[vcol]]) & df[[vcol]] > 0 & !is.na(df$mutation_category) & df$mutation_category != "", , drop = FALSE]
      if (nrow(df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 6, dpi = 150); return(invisible(NULL)) }
      df$.metric <- df[[vcol]]
      med <- aggregate(.metric ~ mutation_category, data = df, FUN = median)
      cat_order <- as.character(med$mutation_category[order(med$.metric)])
      df$mutation_category <- factor(as.character(df$mutation_category), levels = cat_order)
      p <- ggplot(df, aes(x = .metric, y = mutation_category, fill = mutation_category)) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) + scale_x_continuous(limits = c(0, 100)) +
        scale_fill_manual(values = PAL_MUT_CAT, na.value = "gray70", guide = "none") +
        labs(x = vlabel, y = NULL) + theme_minimal(base_size = 14) +
        theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13))
      ggsave(file, plot = p, width = 8, height = max(4, 0.4 * length(cat_order) + 2), dpi = 300)
    }
  )

  output$download_vaf_cohort_plot <- downloadHandler(
    filename = function() paste0(tolower(vaf_col()), "_by_cohort_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      df <- filtered_data(); vcol <- vaf_col(); vlabel <- vaf_label()
      if (is.null(df)) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 6, dpi = 150); return(invisible(NULL)) }
      if (!vcol %in% colnames(df)) vcol <- "VAF"
      df <- df[!is.na(df[[vcol]]) & df[[vcol]] > 0 & !is.na(df$Cohort), , drop = FALSE]
      if (nrow(df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 6, dpi = 150); return(invisible(NULL)) }
      df$.metric <- df[[vcol]]
      cohort_pal <- if (is_meta4()) META4_COHORT_COL else PAL_COHORT
      present_cohorts <- unique(df$Cohort[!is.na(df$Cohort)]); cohort_pal <- cohort_pal[names(cohort_pal) %in% present_cohorts]
      med <- aggregate(.metric ~ Cohort, data = df, FUN = median)
      cohort_order <- as.character(med$Cohort[order(med$.metric)])
      df$Cohort <- factor(as.character(df$Cohort), levels = cohort_order)
      p <- ggplot(df, aes(x = .metric, y = Cohort, fill = Cohort)) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) + scale_x_continuous(limits = c(0, 100)) +
        scale_fill_manual(values = cohort_pal, na.value = "gray70", guide = "none") +
        labs(x = vlabel, y = NULL) + theme_minimal(base_size = 14) +
        theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13))
      ggsave(file, plot = p, width = 8, height = max(4, 0.5 * length(cohort_order) + 2), dpi = 300)
    }
  )

  # --- Clonal Ordering Survival: All Gene Pairs ---
  clonal_hr_all_pairs <- reactive({
    req(df <- filtered_data())
    vcol <- vaf_col()
    if (!vcol %in% colnames(df)) vcol <- "VAF"
    if (!vcol %in% colnames(df)) return(NULL)
    surv_df <- survival_data()
    if (is.null(surv_df) || nrow(surv_df) == 0) return(NULL)
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    genes <- filtered_genes()
    if (length(genes) < 2) return(NULL)
    pairs <- combn(genes, 2, simplify = FALSE)
    res_list <- list()
    for (pair in pairs) {
      g1 <- pair[1]; g2 <- pair[2]
      sub1 <- df[as.character(df$Gene_for_analysis) == g1 & !is.na(df[[vcol]]) & df[[vcol]] > 0, c("Sample", vcol), drop = FALSE]
      sub2 <- df[as.character(df$Gene_for_analysis) == g2 & !is.na(df[[vcol]]) & df[[vcol]] > 0, c("Sample", vcol), drop = FALSE]
      if (nrow(sub1) == 0 || nrow(sub2) == 0) next
      colnames(sub1)[2] <- "v1"; colnames(sub2)[2] <- "v2"
      sub1 <- aggregate(v1 ~ Sample, data = sub1, FUN = max)
      sub2 <- aggregate(v2 ~ Sample, data = sub2, FUN = max)
      both <- merge(sub1, sub2, by = "Sample")
      if (nrow(both) < 5) next
      diff_val <- both$v1 - both$v2
      both$Order <- ifelse(diff_val > 5, paste(g1, "first"), ifelse(diff_val < -5, paste(g2, "first"), "Ambiguous"))
      both <- both[both$Order != "Ambiguous", , drop = FALSE]
      n_g1 <- sum(both$Order == paste(g1, "first"))
      n_g2 <- sum(both$Order == paste(g2, "first"))
      if (n_g1 < 10 || n_g2 < 10) next
      km_df <- merge(both[, c("Sample", "Order"), drop = FALSE], surv_uniq, by = "Sample")
      if (nrow(km_df) < 20) next
      km_df$Order <- factor(km_df$Order, levels = c(paste(g2, "first"), paste(g1, "first")))
      cx <- tryCatch(survival::coxph(Surv(Time_to_OS, as.numeric(Censor)) ~ Order, data = km_df), error = function(e) NULL)
      if (is.null(cx)) next
      ci <- tryCatch(exp(confint(cx)), error = function(e) NULL)
      if (is.null(ci)) next
      cx_sum <- summary(cx)
      coef_name <- paste0("Order", paste(g1, "first"))
      hr <- exp(coef(cx))[coef_name]
      ci_lo <- ci[coef_name, 1]; ci_hi <- ci[coef_name, 2]
      pval <- cx_sum$coefficients[coef_name, "Pr(>|z|)"]
      res_list[[length(res_list) + 1L]] <- data.frame(
        Pair = paste(g1, "vs", g2), Gene_A = g1, Gene_B = g2,
        HR = as.numeric(hr), CI_lower = as.numeric(ci_lo), CI_upper = as.numeric(ci_hi),
        p_value = as.numeric(pval), n_A_first = n_g1, n_B_first = n_g2,
        n_total = n_g1 + n_g2, stringsAsFactors = FALSE)
    }
    if (length(res_list) == 0) return(NULL)
    out <- do.call(rbind, res_list)
    out$p_adj <- p.adjust(out$p_value, method = "fdr")
    out <- out[order(out$HR), , drop = FALSE]
    out
  })

  output$clonal_hr_summary_plot <- renderPlot({
    hr_df <- clonal_hr_all_pairs()
    vcol <- vaf_col()
    metric_short <- if (vcol == "CCF") "CCF" else "VAF"
    if (is.null(hr_df) || nrow(hr_df) == 0) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5,
        label = paste0("No gene pairs with \u226510 non-ambiguous orderings per group (", metric_short, " difference > 5%)"),
        size = 5) + theme_void())
    }
    hr_df$Pair <- factor(hr_df$Pair, levels = hr_df$Pair)
    hr_df$Significance <- ifelse(hr_df$p_adj < 0.05, "FDR < 0.05",
      ifelse(hr_df$p_value < 0.05, "p < 0.05", "Not significant"))
    sig_colors <- c("FDR < 0.05" = "#762a83", "p < 0.05" = "#e7298a", "Not significant" = "#9E9E9E")
    hr_df$label <- paste0("HR=", round(hr_df$HR, 2), " (", round(hr_df$CI_lower, 2), "\u2013", round(hr_df$CI_upper, 2), ")")
    ggplot(hr_df, aes(x = HR, y = Pair, xmin = CI_lower, xmax = CI_upper, color = Significance)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      geom_errorbarh(height = 0.25, linewidth = 0.6) +
      geom_point(aes(size = n_total), alpha = 0.8) +
      scale_color_manual(values = sig_colors, name = NULL) +
      scale_size_continuous(name = "Non-ambiguous\ncases", range = c(2, 8)) +
      scale_x_continuous(trans = "log10", name = paste0("Hazard Ratio (Gene A first vs Gene B first, 95% CI)")) +
      labs(y = NULL, title = paste0("Clonal ordering survival: all gene pairs (", metric_short, ")"),
        subtitle = paste0("HR > 1: worse survival when Gene A (left of 'vs') has higher ", metric_short)) +
      theme_minimal(base_size = 14) +
      theme(axis.text.y = element_text(size = 11), axis.title = element_text(size = 13),
        plot.title = element_text(size = 14, face = "bold"), plot.subtitle = element_text(size = 11, color = "gray40"))
  })

  output$download_clonal_hr_summary_plot <- downloadHandler(
    filename = function() paste0("clonal_hr_all_pairs_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      hr_df <- clonal_hr_all_pairs()
      if (is.null(hr_df) || nrow(hr_df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 10, height = 6, dpi = 150); return(invisible(NULL)) }
      vcol <- vaf_col(); metric_short <- if (vcol == "CCF") "CCF" else "VAF"
      hr_df$Pair <- factor(hr_df$Pair, levels = hr_df$Pair)
      hr_df$Significance <- ifelse(hr_df$p_adj < 0.05, "FDR < 0.05", ifelse(hr_df$p_value < 0.05, "p < 0.05", "Not significant"))
      sig_colors <- c("FDR < 0.05" = "#762a83", "p < 0.05" = "#e7298a", "Not significant" = "#9E9E9E")
      p <- ggplot(hr_df, aes(x = HR, y = Pair, xmin = CI_lower, xmax = CI_upper, color = Significance)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
        geom_errorbarh(height = 0.25, linewidth = 0.6) +
        geom_point(aes(size = n_total), alpha = 0.8) +
        scale_color_manual(values = sig_colors, name = NULL) +
        scale_size_continuous(name = "Non-ambiguous\ncases", range = c(2, 8)) +
        scale_x_continuous(trans = "log10", name = "Hazard Ratio (Gene A first vs Gene B first, 95% CI)") +
        labs(y = NULL, title = paste0("Clonal ordering survival: all gene pairs (", metric_short, ")"),
          subtitle = paste0("HR > 1: worse survival when Gene A (left of 'vs') has higher ", metric_short)) +
        theme_minimal(base_size = 14) +
        theme(axis.text.y = element_text(size = 11), axis.title = element_text(size = 13),
          plot.title = element_text(size = 14, face = "bold"), plot.subtitle = element_text(size = 11, color = "gray40"))
      ggsave(file, plot = p, width = 10, height = max(5, 0.35 * nrow(hr_df) + 2), dpi = 300)
    }
  )

  # --- VAF vs CCF scatterplot ---
  output$vaf_vs_ccf_plot <- renderPlot({
    req(df <- filtered_data())
    if (!"VAF" %in% colnames(df) || !"CCF" %in% colnames(df)) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "CCF column not available", size = 5) + theme_void())
    }
    df <- df[!is.na(df$VAF) & df$VAF > 0 & !is.na(df$CCF), , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data with both VAF and CCF") + theme_void())
    cn_col <- if ("CN_at_locus" %in% colnames(df)) "CN_at_locus" else NULL
    if (!is.null(cn_col)) {
      df$CN_label <- factor(ifelse(is.na(df[[cn_col]]), "Unknown",
        ifelse(df[[cn_col]] == 1, "CN = 1 (loss/hemizygous)",
        ifelse(df[[cn_col]] == 2, "CN = 2 (diploid)",
        ifelse(df[[cn_col]] == 3, "CN = 3 (gain)", paste0("CN = ", df[[cn_col]]))))),
        levels = c("CN = 1 (loss/hemizygous)", "CN = 2 (diploid)", "CN = 3 (gain)", "Unknown"))
      cn_pal <- c("CN = 1 (loss/hemizygous)" = "#d73027", "CN = 2 (diploid)" = "#4575b4", "CN = 3 (gain)" = "#1a9850", "Unknown" = "#999999")
      p <- ggplot(df, aes(x = VAF, y = CCF, color = CN_label)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
        geom_point(alpha = 0.4, size = 1.8) +
        scale_color_manual(values = cn_pal, name = "Copy number", drop = FALSE) +
        scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 100)) +
        labs(x = "VAF (%)", y = "CCF (%)", title = "VAF vs CCF") +
        theme_minimal(base_size = 14) +
        theme(legend.position = "bottom", legend.text = element_text(size = 11),
          axis.text = element_text(size = 12), axis.title = element_text(size = 13))
    } else {
      p <- ggplot(df, aes(x = VAF, y = CCF)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.5) +
        geom_point(alpha = 0.4, size = 1.8, color = "#4575b4") +
        scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 100)) +
        labs(x = "VAF (%)", y = "CCF (%)", title = "VAF vs CCF") +
        theme_minimal(base_size = 14) +
        theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13))
    }
    n_total <- nrow(df)
    n_above <- sum(df$CCF > df$VAF, na.rm = TRUE)
    n_on <- sum(abs(df$CCF - df$VAF) < 0.5, na.rm = TRUE)
    n_below <- sum(df$CCF < df$VAF, na.rm = TRUE)
    cor_val <- round(cor(df$VAF, df$CCF, use = "complete.obs"), 3)
    p + annotate("text", x = 5, y = 95, hjust = 0, vjust = 1, size = 3.8, color = "gray30",
      label = paste0("n = ", n_total, "  |  r = ", cor_val,
        "\nAbove diagonal (CN>2): ", n_above,
        "\nBelow diagonal (CN<2): ", n_below))
  })

  output$download_vaf_vs_ccf_plot <- downloadHandler(
    filename = function() paste0("vaf_vs_ccf_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      df <- filtered_data()
      if (is.null(df) || !"VAF" %in% colnames(df) || !"CCF" %in% colnames(df)) {
        ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 7, height = 7, dpi = 150)
        return(invisible(NULL))
      }
      df <- df[!is.na(df$VAF) & df$VAF > 0 & !is.na(df$CCF), , drop = FALSE]
      if (nrow(df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 7, height = 7, dpi = 150); return(invisible(NULL)) }
      cn_col <- if ("CN_at_locus" %in% colnames(df)) "CN_at_locus" else NULL
      if (!is.null(cn_col)) {
        df$CN_label <- factor(ifelse(is.na(df[[cn_col]]), "Unknown",
          ifelse(df[[cn_col]] == 1, "CN = 1 (loss/hemizygous)",
          ifelse(df[[cn_col]] == 2, "CN = 2 (diploid)",
          ifelse(df[[cn_col]] == 3, "CN = 3 (gain)", paste0("CN = ", df[[cn_col]]))))),
          levels = c("CN = 1 (loss/hemizygous)", "CN = 2 (diploid)", "CN = 3 (gain)", "Unknown"))
        cn_pal <- c("CN = 1 (loss/hemizygous)" = "#d73027", "CN = 2 (diploid)" = "#4575b4", "CN = 3 (gain)" = "#1a9850", "Unknown" = "#999999")
        p <- ggplot(df, aes(x = VAF, y = CCF, color = CN_label)) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
          geom_point(alpha = 0.4, size = 1.8) +
          scale_color_manual(values = cn_pal, name = "Copy number", drop = FALSE) +
          scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 100)) +
          labs(x = "VAF (%)", y = "CCF (%)", title = "VAF vs CCF") +
          theme_minimal(base_size = 14) + theme(legend.position = "bottom")
      } else {
        p <- ggplot(df, aes(x = VAF, y = CCF)) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
          geom_point(alpha = 0.4, size = 1.8, color = "#4575b4") +
          scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 100)) +
          labs(x = "VAF (%)", y = "CCF (%)", title = "VAF vs CCF") +
          theme_minimal(base_size = 14)
      }
      cor_val <- round(cor(df$VAF, df$CCF, use = "complete.obs"), 3)
      p <- p + annotate("text", x = 5, y = 95, hjust = 0, vjust = 1, size = 3.8, color = "gray30",
        label = paste0("n = ", nrow(df), "  |  r = ", cor_val))
      ggsave(file, plot = p, width = 7, height = 7, dpi = 300)
    }
  )

  forest_data <- reactive({
    surv_df <- survival_data()
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    # Apply gene frequency filter to determine which genes to show in forest plot (all samples still used for analysis)
    gene_tbl <- table(as.character(surv_df$Gene_for_analysis))
    min_pct <- 1
    if (min_pct > 0) {
      n_samples <- length(unique(surv_df$Sample))
      if (n_samples > 0) {
        gene_pct <- (gene_tbl / n_samples) * 100
        genes <- names(gene_pct)[gene_pct >= min_pct]
      } else {
        genes <- character(0)
      }
    } else {
      genes <- names(gene_tbl)
    }
    res <- list()
    for (g in genes) {
      mut_pts <- unique(surv_df$Sample[surv_df$Gene_for_analysis == g])
      surv_uniq$mut <- surv_uniq$Sample %in% mut_pts
      n_mut <- sum(surv_uniq$mut)
      n_wt <- sum(!surv_uniq$mut)
      if (n_mut >= 10 && n_wt >= 10) {
        m <- tryCatch(coxph(Surv(Time_to_OS, as.numeric(Censor)) ~ mut, data = surv_uniq), error = function(e) NULL)
        if (!is.null(m)) {
          ci <- confint(m)
          res[[g]] <- data.frame(Gene = g, n_mut = n_mut, n_wt = n_wt, HR = exp(coef(m))["mutTRUE"], lower = exp(ci[1]), upper = exp(ci[2]),
            p = summary(m)$coefficients["mutTRUE", "Pr(>|z|)"], stringsAsFactors = FALSE)
        }
      }
    }
    out <- do.call(rbind, res)
    if (!is.null(out)) out <- out[!is.na(out$HR) & out$lower < 10 & out$upper < 10, , drop = FALSE]
    else out <- data.frame(Gene = character(), n_mut = integer(), n_wt = integer(), HR = numeric(), lower = numeric(), upper = numeric(), p = numeric())
    out
  })

  # Co-mutation survival analysis
  comut_survival_data <- reactive({
    n_genes <- as.numeric(input$comut_n)
    g1 <- input$comut_gene1
    g2 <- input$comut_gene2
    g3 <- if (n_genes == 3) input$comut_gene3 else NULL
    if (is.null(g1) || g1 == "" || is.null(g2) || g2 == "") return(NULL)
    if (n_genes == 3 && (is.null(g3) || g3 == "")) return(NULL)
    genes <- c(g1, g2)
    if (n_genes == 3) genes <- c(g1, g2, g3)
    if (length(unique(genes)) != length(genes)) return(NULL)  # duplicate genes

    surv_df <- survival_data()
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    mut_list <- list()
    for (g in genes) {
      mut_list[[g]] <- unique(surv_df$Sample[surv_df$Gene_for_analysis == g])
    }

    if (n_genes == 2) {
      has_g1 <- surv_uniq$Sample %in% mut_list[[g1]]
      has_g2 <- surv_uniq$Sample %in% mut_list[[g2]]
      surv_uniq$Group <- ifelse(!has_g1 & !has_g2, "Neither",
        ifelse(has_g1 & !has_g2, paste0(g1, " only"),
          ifelse(!has_g1 & has_g2, paste0(g2, " only"), paste0(g1, " + ", g2))))
    } else {
      has_g1 <- surv_uniq$Sample %in% mut_list[[g1]]
      has_g2 <- surv_uniq$Sample %in% mut_list[[g2]]
      has_g3 <- surv_uniq$Sample %in% mut_list[[g3]]
      n_mut <- has_g1 + has_g2 + has_g3
      surv_uniq$Group <- character(nrow(surv_uniq))
      surv_uniq$Group[n_mut == 0] <- "None"
      surv_uniq$Group[n_mut == 1 & has_g1] <- paste0(g1, " only")
      surv_uniq$Group[n_mut == 1 & has_g2] <- paste0(g2, " only")
      surv_uniq$Group[n_mut == 1 & has_g3] <- paste0(g3, " only")
      surv_uniq$Group[n_mut == 2 & has_g1 & has_g2] <- paste0(g1, " + ", g2)
      surv_uniq$Group[n_mut == 2 & has_g1 & has_g3] <- paste0(g1, " + ", g3)
      surv_uniq$Group[n_mut == 2 & has_g2 & has_g3] <- paste0(g2, " + ", g3)
      surv_uniq$Group[n_mut == 3] <- paste0(g1, " + ", g2, " + ", g3)
    }
    surv_uniq
  })

  output$comut_survival_plot_ui <- renderUI({
    h <- if (as.character(input$comut_n) == "3") 750 else 500
    plotOutput("comut_survival_plot", height = h)
  })

  output$comut_survival_plot <- renderPlot({
    req(df <- comut_survival_data())
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    n_groups <- length(unique(df$Group))
    if (n_groups < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Need at least 2 groups for comparison"))
    if (any(table(df$Group) < 3)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Each group needs at least 3 patients"))

    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = df)
    sd <- survdiff(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = df)
    pval <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
    pval_txt <- if (pval < 0.001) "p < 0.001" else paste0("p = ", format.pval(pval, digits = 2))

    if (has_survminer) {
      n_groups <- length(unique(df$Group))
      pal_all <- unname(c(PAL_SUBSET, PAL_MUT_CAT, PAL_COHORT))
      pal_vec <- pal_all[seq_len(n_groups)]
      if (length(pal_vec) < n_groups) {
        pal_vec <- rep(pal_all, length.out = n_groups)
      }
      legend_labs <- gsub("^Group=", "", names(fit$strata))
      p <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, pval = pval_txt,
        title = "Survival by Co-mutation Status", xlab = "Years",
        palette = pal_vec, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.05))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install the survminer package for Kaplan-Meier plots with number at risk.", size = 4) + theme_void()
    }
  })

  output$download_oncoprint_plot <- downloadHandler(
    filename = function() paste0("oncoprint_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      png(file, width = 1400, height = 700, res = 150)
      od <- oncoprint_data()
      if (is.null(od) || !is.matrix(od$temp_dat) || nrow(od$temp_dat) == 0) {
        plot.new(); text(0.5, 0.5, "No data")
      } else if (has_ComplexHeatmap) {
        temp_dat <- od$temp_dat; anno_df <- od$anno_df; use_meta4 <- is_meta4()
        col <- if (use_meta4) META4_VAR_COL else ONCO_VAR_COL
        vtypes <- unique(as.character(temp_dat[temp_dat != ""]))
        for (vt in setdiff(vtypes, names(col))) col[vt] <- "#80796BFF"
        alter_fun_list <- list(
          background = function(x, y, w, h) NULL,
          Deletion = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Deletion"], col = "#374E55FF")) },
          Indel = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Indel"], col = "#DF8F44FF")) },
          ITD = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["ITD"], col = "#79AF97FF")) },
          SNV = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["SNV"], col = "#B24745FF")) },
          Splicing = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Splicing"], col = "#6A6599FF")) },
          PTD = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["PTD"], col = "tan")) },
          Other = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Other"], col = "#80796BFF")) },
          Unknown = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Unknown"], col = "#80796BFF")) }
        )
        for (vt in setdiff(vtypes, names(alter_fun_list))) alter_fun_list[[vt]] <- (function(fc) function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = fc, col = fc)) })(col[vt])
        fig <- ComplexHeatmap::oncoPrint(temp_dat, col = col, row_names_side = "right", alter_fun_is_vectorized = TRUE, alter_fun = alter_fun_list)
        ComplexHeatmap::draw(fig)
      } else {
        plot.new(); text(0.5, 0.5, "ComplexHeatmap not installed")
      }
      dev.off()
    }
  )

  output$download_cooccurrence_plot <- downloadHandler(
    filename = function() paste0("cooccurrence_odds_ratio_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      cooc <- cooccurrence_data()
      if (is.null(cooc$matrix) || nrow(cooc$matrix) < 2) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 8, dpi = 150); return(invisible(NULL)) }
      keep_genes <- filtered_genes()
      keep_genes <- intersect(keep_genes, rownames(cooc$matrix))
      if (length(keep_genes) < 2) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 8, dpi = 150); return(invisible(NULL)) }
      mat <- cooc$matrix[keep_genes, keep_genes, drop = FALSE]
      mat_log <- log2(mat + 0.01); mat_log[mat_log > 2] <- 2; mat_log[mat_log < -2] <- -2
      ord <- hclust(dist(mat_log))$order; mat_ord <- mat_log[ord, ord]
      df_plot <- data.frame(Gene1 = rep(rownames(mat_ord), ncol(mat_ord)), Gene2 = rep(colnames(mat_ord), each = nrow(mat_ord)), log2OR = as.vector(mat_ord))
      df_plot$Gene1 <- factor(df_plot$Gene1, levels = rownames(mat_ord)); df_plot$Gene2 <- factor(df_plot$Gene2, levels = colnames(mat_ord))
      p <- ggplot(df_plot, aes(x = Gene1, y = Gene2, fill = log2OR)) + geom_tile() + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) + coord_fixed() + labs(x = NULL, y = NULL, fill = "log2(OR)") + theme_minimal(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      ggsave(file, plot = p, width = 10, height = 10, dpi = 150)
    }
  )

  output$download_cooccurrence_table <- downloadHandler(
    filename = function() paste0("cooccurrence_pairs_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
    content = function(file) {
      cooc <- cooccurrence_data()
      if (is.null(cooc$pairs) || nrow(cooc$pairs) == 0) { write.csv(data.frame(), file, row.names = FALSE); return(invisible(NULL)) }
      keep_genes <- filtered_genes()
      df <- if (length(keep_genes) > 0) cooc$pairs[cooc$pairs$gene1 %in% keep_genes & cooc$pairs$gene2 %in% keep_genes, , drop = FALSE] else cooc$pairs
      df$pair <- paste(df$gene1, "+", df$gene2)
      df <- df[order(-df$n_cooccur), c("pair", "odds_ratio", "n_cooccur", "p", "q"), drop = FALSE]
      write.csv(df, file, row.names = FALSE)
    }
  )

  output$download_comut_survival_plot <- downloadHandler(
    filename = function() paste0("comut_survival_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      df <- comut_survival_data()
      if (is.null(df) || nrow(df) == 0 || length(unique(df$Group)) < 2) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 6, dpi = 150); return(invisible(NULL)) }
      fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = df)
      if (has_survminer) {
        n_groups <- length(unique(df$Group))
        pal_all <- unname(c(PAL_SUBSET, PAL_MUT_CAT, PAL_COHORT))
        pal_vec <- rep(pal_all, length.out = n_groups)
        legend_labs <- gsub("^Group=", "", names(fit$strata))
        p <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE, title = "Survival by Co-mutation Status", xlab = "Years", palette = pal_vec, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.05))
        png(file, width = 1200, height = 800, res = 150)
        print(p)
        dev.off()
      } else {
        ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "survminer not installed") + theme_void(), width = 8, height = 6, dpi = 150)
      }
    }
  )

  # ---- UpSet plot: uses only the genes selected in the co-mutation KM panel ----
  upset_matrix <- reactive({
    g1 <- input$comut_gene1
    g2 <- input$comut_gene2
    sel_genes <- unique(c(g1, g2))
    if (as.numeric(input$comut_n) == 3) sel_genes <- unique(c(sel_genes, input$comut_gene3))
    sel_genes <- sel_genes[!is.null(sel_genes) & sel_genes != ""]
    if (length(sel_genes) < 2) return(NULL)
    df <- filtered_data()[, c("Sample", "Gene_for_analysis"), drop = FALSE]
    df <- df[!duplicated(df), , drop = FALSE]
    df <- df[as.character(df$Gene_for_analysis) %in% sel_genes, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    wide <- as.data.frame.matrix(xtabs(~ Sample + Gene_for_analysis, data = df))
    wide[wide > 0] <- 1L
    if (ncol(wide) < 2) return(NULL)
    wide
  })

  # Build UpSetR queries to color intersection bars using the same KM palette
  upset_queries <- reactive({
    wide <- upset_matrix()
    if (is.null(wide) || ncol(wide) < 2) return(list())
    genes <- colnames(wide)
    n_genes <- length(genes)
    # Enumerate all non-empty subsets (same groups as the KM plot)
    combos <- list()
    for (k in seq_len(n_genes)) {
      cm <- combn(genes, k, simplify = FALSE)
      combos <- c(combos, cm)
    }
    # Build KM-style group labels to determine alphabetical color order
    km_labels <- character(length(combos))
    for (i in seq_along(combos)) {
      s <- combos[[i]]
      if (length(s) == n_genes) {
        km_labels[i] <- paste(s, collapse = " + ")
      } else {
        km_labels[i] <- paste0(paste(s, collapse = " + "), if (length(s) == 1) " only" else "")
      }
    }
    # KM groups include "Neither"/"None" as the first alphabetically if present;
    # survfit orders strata alphabetically on the Group factor
    all_labels <- sort(c(if (n_genes == 2) "Neither" else "None", km_labels))
    pal_all <- unname(c(PAL_SUBSET, PAL_MUT_CAT, PAL_COHORT))
    pal_vec <- rep(pal_all, length.out = length(all_labels))
    names(pal_vec) <- all_labels
    queries <- list()
    for (i in seq_along(combos)) {
      col <- pal_vec[km_labels[i]]
      queries[[i]] <- list(query = UpSetR::intersects, params = combos[[i]], color = unname(col), active = TRUE)
    }
    queries
  })

  upset_render <- function(wide, queries) {
    n_sets <- min(ncol(wide), 20)
    UpSetR::upset(wide, nsets = n_sets, order.by = "freq",
                  mainbar.y.label = "Intersection Size", sets.x.label = "Samples per Gene",
                  queries = queries, query.legend = "none",
                  point.size = 3.5, line.size = 1.2,
                  text.scale = c(1.8, 1.6, 1.5, 1.4, 1.8, 1.5))
  }

  output$upset_plot <- renderPlot({
    wide <- upset_matrix()
    if (is.null(wide) || ncol(wide) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for UpSet plot") + theme_void())
    if (!has_UpSetR) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "UpSetR package not installed") + theme_void())
    upset_render(wide, upset_queries())
  })

  output$download_upset_plot <- downloadHandler(
    filename = function() paste0("upset_mutations_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      wide <- upset_matrix()
      if (is.null(wide) || ncol(wide) < 2 || !has_UpSetR) {
        png(file, width = 1200, height = 800, res = 150)
        plot.new(); text(0.5, 0.5, "Insufficient data or UpSetR not installed")
        dev.off()
        return(invisible(NULL))
      }
      png(file, width = 1200, height = 800, res = 150)
      print(upset_render(wide, upset_queries()))
      dev.off()
    }
  )

  output$forest_plot <- renderPlot({
    df <- forest_data()
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    df <- df[order(df$HR), , drop = FALSE]
    df$Gene <- factor(df$Gene, levels = df$Gene)
    df$ColorGroup <- ifelse(df$p >= 0.05, "Not significant",
      ifelse(df$HR > 1, "Significant (HR > 1)", "Significant (HR < 1)"))
    n_genes <- nrow(df)
    y_axis_size <- max(7, min(14, round(420 / n_genes)))
    ggplot(df, aes(x = Gene, y = HR, ymin = lower, ymax = upper, color = ColorGroup)) +
      geom_pointrange(size = 0.5) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      coord_flip() +
      scale_color_manual(values = c("Significant (HR > 1)" = "#762a83", "Significant (HR < 1)" = "#1b7837", "Not significant" = "#9E9E9E"), name = NULL, guide = "none") +
      labs(x = NULL, y = "Hazard Ratio") +
      theme_minimal(base_size = 16) +
      theme(axis.text.y = element_text(size = y_axis_size), axis.text.x = element_text(size = 14), axis.title = element_text(size = 15))
  })

  output$survival_table <- renderDT({
    df <- forest_data()
    if (is.null(df) || nrow(df) == 0) return(datatable(data.frame(Gene = character(), n_mut = integer(), n_wt = integer(), HR_CI = character(), p_value = character(), p_adj = character()), options = list(pageLength = 10000, lengthChange = FALSE, scrollY = "350px", scrollCollapse = TRUE), rownames = FALSE))
    df$HR_CI <- sprintf("%.2f (%.2f-%.2f)", df$HR, df$lower, df$upper)
    df$p_value <- format_pval_display(df$p)
    df$p_adj <- format_pval_display(p.adjust(df$p, method = "fdr"))
    datatable(df[, c("Gene", "n_mut", "n_wt", "HR_CI", "p_value", "p_adj")], options = list(pageLength = 10000, lengthChange = FALSE, scrollY = "350px", scrollCollapse = TRUE), rownames = FALSE)
  })

  output$download_forest_plot <- downloadHandler(
    filename = function() paste0("survival_forest_plot_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      df <- forest_data()
      if (is.null(df) || nrow(df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 5, dpi = 150); return(invisible(NULL)) }
      df <- df[order(df$HR), , drop = FALSE]
      df$Gene <- factor(df$Gene, levels = df$Gene)
      df$ColorGroup <- ifelse(df$p >= 0.05, "Not significant", ifelse(df$HR > 1, "Significant (HR > 1)", "Significant (HR < 1)"))
      n_genes <- nrow(df)
      y_axis_size <- max(7, min(14, round(420 / n_genes)))
      p <- ggplot(df, aes(x = Gene, y = HR, ymin = lower, ymax = upper, color = ColorGroup)) +
        geom_pointrange(size = 0.5) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
        coord_flip() +
        scale_color_manual(values = c("Significant (HR > 1)" = "#762a83", "Significant (HR < 1)" = "#1b7837", "Not significant" = "#9E9E9E"), name = NULL, guide = "none") +
        labs(x = NULL, y = "Hazard Ratio") +
        theme_minimal(base_size = 16) +
        theme(axis.text.y = element_text(size = y_axis_size), axis.text.x = element_text(size = 14), axis.title = element_text(size = 15))
      ggsave(file, plot = p, width = max(7, n_genes * 0.2), height = max(5, n_genes * 0.15), dpi = 300)
    }
  )

  output$download_survival_plot <- downloadHandler(
    filename = function() paste0("survival_km_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      gene <- input$surv_gene
      if (is.null(gene) || gene == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene", size = 10) + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      surv_df <- survival_data()
      mut_samples <- unique(surv_df$Sample[surv_df$Gene_for_analysis == gene])
      surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
      surv_uniq$Mutation <- ifelse(surv_uniq$Sample %in% mut_samples, paste0(gene, " mut"), "WT")
      if (length(unique(surv_uniq$Mutation)) < 2) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient groups", size = 10) + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Mutation, data = surv_uniq)
      strata_names <- names(fit$strata)
      is_mut <- grepl(" mut", strata_names, fixed = TRUE)
      surv_pal <- ifelse(is_mut, "#8B0000", "#4D4D4D")
      if (has_survminer) {
        legend_labs <- gsub("^Mutation=", "", names(fit$strata))
        p <- survminer::ggsurvplot(fit, data = surv_uniq, risk.table = TRUE, pval = TRUE, title = paste("Survival by", gene, "mutation status"), xlab = "Years", palette = surv_pal, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.05))
        png(file, width = 7, height = 6, units = "in", res = 300)
        print(p)
        dev.off()
      } else {
        ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install survminer for KM plots") + theme_void(), width = 7, height = 5, dpi = 150)
      }
    }
  )

  output$download_survival_table <- downloadHandler(
    filename = function() paste0("survival_summary_table_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
    content = function(file) {
      df <- forest_data()
      if (is.null(df) || nrow(df) == 0) return(write("No data", file))
      df$HR_CI <- sprintf("%.2f (%.2f-%.2f)", df$HR, df$lower, df$upper)
      df$p_value <- df$p
      df$p_adj <- p.adjust(df$p, method = "fdr")
      out <- df[, c("Gene", "n_mut", "n_wt", "HR_CI", "p_value", "p_adj"), drop = FALSE]
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$clinical_plot <- renderPlot({
    var <- input$clin_var
    if (is.null(var) || var == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a variable"))
    df <- filtered_data()
    if (is.null(df) || !var %in% colnames(df)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Variable not in data"))
    df_one <- df[!duplicated(df$Sample), c("Sample", "Gene_for_analysis", var), drop = FALSE]
    df_one[[var]] <- as.numeric(df_one[[var]])
    df_one <- df_one[!is.na(df_one[[var]]), , drop = FALSE]
    # Apply same filters as in All Gene Associations tab
    if (var == "WBC") {
      df_one <- df_one[df_one[[var]] <= 200, , drop = FALSE]
    } else if (var == "Hemoglobin") {
      df_one <- df_one[df_one[[var]] <= 15, , drop = FALSE]
    } else if (var == "Platelet") {
      df_one <- df_one[df_one[[var]] <= 300, , drop = FALSE]
    }
    if (nrow(df_one) < 5) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data"))
    # Apply gene frequency filter to determine which genes to show
    gene_tbl <- table(df_one$Gene_for_analysis)
    min_pct <- 1
    if (min_pct > 0) {
      n_samples <- length(unique(df_one$Sample))
      if (n_samples > 0) {
        gene_pct <- (gene_tbl / n_samples) * 100
        keep_genes <- names(gene_pct)[gene_pct >= min_pct]
        gene_tbl <- gene_tbl[names(gene_tbl) %in% keep_genes]
      }
    }
    top_genes <- names(sort(gene_tbl, decreasing = TRUE))[1:min(20, length(gene_tbl))]
    df_one <- df_one[df_one$Gene_for_analysis %in% top_genes, , drop = FALSE]
    df_one$Gene <- factor(df_one$Gene_for_analysis, levels = top_genes)
    # Use same unit labels as in All Gene Associations tab
    unit_labels <- c("WBC" = "WBC (1e-9/l)", "Age" = "Age (years)", "Hemoglobin" = "Hemoglobin (g/dl)", 
                     "Platelet" = "Platelet (1e-9/l)", "BM_blast_percent" = "BM blasts (%)", "PB_blast_percent" = "PB blasts (%)")
    y_label <- if (var %in% names(unit_labels)) unit_labels[var] else var
    ggplot(df_one, aes(x = Gene, y = .data[[var]], fill = Gene)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_fill_discrete(guide = "none") +
      labs(x = NULL, y = y_label, title = paste("Distribution of", y_label, "by gene")) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
  })

  output$download_clinical_plot <- downloadHandler(
    filename = function() paste0("clinical_by_mutation_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      var <- input$clin_var
      if (is.null(var) || var == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a variable") + theme_void(), width = 8, height = 5, dpi = 150); return(invisible(NULL)) }
      df <- filtered_data()
      if (is.null(df) || !var %in% colnames(df)) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Variable not in data") + theme_void(), width = 8, height = 5, dpi = 150); return(invisible(NULL)) }
      df_one <- df[!duplicated(df$Sample), c("Sample", "Gene_for_analysis", var), drop = FALSE]
      df_one[[var]] <- as.numeric(df_one[[var]])
      df_one <- df_one[!is.na(df_one[[var]]), , drop = FALSE]
      if (var == "WBC") df_one <- df_one[df_one[[var]] <= 200, , drop = FALSE]
      else if (var == "Hemoglobin") df_one <- df_one[df_one[[var]] <= 15, , drop = FALSE]
      else if (var == "Platelet") df_one <- df_one[df_one[[var]] <= 300, , drop = FALSE]
      if (nrow(df_one) < 5) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") + theme_void(), width = 8, height = 5, dpi = 150); return(invisible(NULL)) }
      gene_tbl <- table(df_one$Gene_for_analysis)
      min_pct <- 1
      if (min_pct > 0) {
        n_samples <- length(unique(df_one$Sample))
        if (n_samples > 0) {
          gene_pct <- (gene_tbl / n_samples) * 100
          keep_genes <- names(gene_pct)[gene_pct >= min_pct]
          gene_tbl <- gene_tbl[names(gene_tbl) %in% keep_genes]
        }
      }
      top_genes <- names(sort(gene_tbl, decreasing = TRUE))[1:min(20, length(gene_tbl))]
      df_one <- df_one[df_one$Gene_for_analysis %in% top_genes, , drop = FALSE]
      df_one$Gene <- factor(df_one$Gene_for_analysis, levels = top_genes)
      unit_labels <- c("WBC" = "WBC (1e-9/l)", "Age" = "Age (years)", "Hemoglobin" = "Hemoglobin (g/dl)", "Platelet" = "Platelet (1e-9/l)", "BM_blast_percent" = "BM blasts (%)", "PB_blast_percent" = "PB blasts (%)")
      y_label <- if (var %in% names(unit_labels)) unit_labels[var] else var
      p <- ggplot(df_one, aes(x = Gene, y = .data[[var]], fill = Gene)) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
        scale_fill_discrete(guide = "none") +
        labs(x = NULL, y = y_label, title = paste("Distribution of", y_label, "by gene")) +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )

  # Co-occurrence computed on all genes (subset/cohort/karyotype only). Min. gene frequency filter is applied only for display.
  cooccurrence_data <- reactive({
    df <- filtered_data()[, c("Sample", "Gene_for_analysis"), drop = FALSE]
    df <- df[!duplicated(df), , drop = FALSE]
    df$Gene <- df$Gene_for_analysis
    df <- df[, c("Sample", "Gene"), drop = FALSE]
    # Use all genes with >= 15 samples for stable Fisher test (no % filter here)
    gene_tbl <- table(as.character(df$Gene))
    genes <- names(which(gene_tbl >= 15))
    if (length(genes) < 2) return(list(matrix = NULL, pairs = NULL))
    df <- df[as.character(df$Gene) %in% genes, , drop = FALSE]
    wide <- as.data.frame.matrix(xtabs(~ Sample + Gene, data = df))
    wide[wide > 0] <- 1
    gene_cols <- colnames(wide)
    if (length(gene_cols) < 2) return(list(matrix = NULL, pairs = NULL))
    pairs <- t(combn(gene_cols, 2))
    results <- list()
    for (i in seq_len(nrow(pairs))) {
      g1 <- pairs[i, 1]
      g2 <- pairs[i, 2]
      tbl <- table(wide[[g1]], wide[[g2]])
      if (any(dim(tbl) < 2)) next
      tryCatch({
        ft <- fisher.test(tbl)
        or <- (tbl[2,2] * tbl[1,1]) / (tbl[2,1] * tbl[1,2])
        if (!is.finite(or)) or <- 0.01
        results[[i]] <- data.frame(gene1 = g1, gene2 = g2, odds_ratio = or, n_cooccur = tbl[2,2], p = ft$p.value, stringsAsFactors = FALSE)
      }, error = function(e) NULL)
    }
    pairs_df <- do.call(rbind, results)
    if (is.null(pairs_df) || nrow(pairs_df) == 0) return(list(matrix = NULL, pairs = NULL))
    pairs_df$q <- p.adjust(pairs_df$p, method = "fdr")
    or_mat <- matrix(1, length(genes), length(genes))
    rownames(or_mat) <- colnames(or_mat) <- genes
    for (i in seq_len(nrow(pairs_df))) {
      g1 <- pairs_df$gene1[i]
      g2 <- pairs_df$gene2[i]
      or_mat[g1, g2] <- or_mat[g2, g1] <- pairs_df$odds_ratio[i]
    }
    list(matrix = or_mat, pairs = pairs_df)
  })

  output$cooccurrence_plot <- renderPlot({
    cooc <- cooccurrence_data()
    if (is.null(cooc$matrix) || nrow(cooc$matrix) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data"))
    # Filter to genes passing min. gene frequency for display only (underlying ORs unchanged)
    keep_genes <- filtered_genes()
    if (length(keep_genes) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No genes pass frequency filter"))
    keep_genes <- intersect(keep_genes, rownames(cooc$matrix))
    if (length(keep_genes) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data"))
    mat <- cooc$matrix[keep_genes, keep_genes, drop = FALSE]
    mat_log <- log2(mat + 0.01)
    mat_log[mat_log > 2] <- 2
    mat_log[mat_log < -2] <- -2
    ord <- hclust(dist(mat_log))$order
    mat_ord <- mat_log[ord, ord]
    n_genes <- nrow(mat_ord)
    # Dynamically size axis text based on number of genes: smaller for more genes, larger for fewer
    # Formula: start at 14 for 10 genes, scale down to ~8 for 50+ genes
    axis_text_size <- max(8, 14 - (n_genes - 10) * 0.15)
    df_plot <- data.frame(Gene1 = rep(rownames(mat_ord), ncol(mat_ord)),
      Gene2 = rep(colnames(mat_ord), each = nrow(mat_ord)),
      log2OR = as.vector(mat_ord))
    df_plot$Gene1 <- factor(df_plot$Gene1, levels = rownames(mat_ord))
    df_plot$Gene2 <- factor(df_plot$Gene2, levels = colnames(mat_ord))
    ggplot(df_plot, aes(x = Gene1, y = Gene2, fill = log2OR)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
      coord_fixed() +
      labs(x = NULL, y = NULL, fill = "log2(OR)", title = "Co-mutation Odds Ratio") +
      theme_minimal(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = axis_text_size), axis.text.y = element_text(size = axis_text_size), legend.title = element_text(size = 14), legend.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 16))
  })

  output$cooccurrence_table <- renderDT({
    cooc <- cooccurrence_data()
    if (is.null(cooc$pairs) || nrow(cooc$pairs) == 0) return(datatable(data.frame(pair = character(), odds_ratio = numeric(), n_cooccur = integer(), p = character(), q = character()), options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE))
    keep_genes <- filtered_genes()
    # Filter pairs to those where both genes pass min. gene frequency (for display only)
    if (length(keep_genes) > 0) {
      df <- cooc$pairs[cooc$pairs$gene1 %in% keep_genes & cooc$pairs$gene2 %in% keep_genes, , drop = FALSE]
    } else {
      df <- cooc$pairs
    }
    if (nrow(df) == 0) return(datatable(data.frame(pair = character(), odds_ratio = numeric(), n_cooccur = integer(), p = character(), q = character()), options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE))
    df$pair <- paste(df$gene1, "+", df$gene2)
    df$odds_ratio <- round(df$odds_ratio, 2)
    df$p <- format_pval_display(df$p)
    df$q <- format_pval_display(df$q)
    df <- df[order(-df$n_cooccur), , drop = FALSE]
    datatable(df[, c("pair", "odds_ratio", "n_cooccur", "p", "q")], options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE)
  })

  # Figure 2: Co-occurrence & Prognosis (de novo)
  fig2_data <- reactive({
    df <- filtered_data()
    df <- df[as.character(df$Subset) == "de_novo", , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df <- df[!is.na(df$Time_to_OS) & !is.na(df$Censor), , drop = FALSE]
    df$Time_to_OS <- as.numeric(df$Time_to_OS) / 365
    df$Gene <- df$Gene_for_analysis
    genes <- names(which(table(as.character(df$Gene)) >= 15))
    if (length(genes) < 2) return(NULL)
    df <- df[as.character(df$Gene) %in% genes, , drop = FALSE]
    wide <- as.data.frame.matrix(xtabs(~ Sample + Gene, data = unique(df[, c("Sample", "Gene")])))
    wide[wide > 0] <- 1
    surv_uniq <- df[!duplicated(df$Sample), c("Sample", "Time_to_OS", "Censor")]
    pairs <- t(combn(genes, 2))
    or_list <- list()
    hr_list <- list()
    for (i in seq_len(nrow(pairs))) {
      g1 <- pairs[i, 1]
      g2 <- pairs[i, 2]
      tbl <- table(wide[[g1]], wide[[g2]])
      if (any(dim(tbl) < 2)) next
      ft <- tryCatch(fisher.test(tbl), error = function(e) NULL)
      if (is.null(ft)) next
      or <- (tbl[2,2] * tbl[1,1]) / (tbl[2,1] * tbl[1,2])
      if (!is.finite(or) || or == 0) or <- 0.01
      n_co <- tbl[2, 2]
      if (n_co < 10) next
      or_list[[length(or_list) + 1]] <- data.frame(gene1 = g1, gene2 = g2, odds_ratio = or, p_or = ft$p.value, n_cooccur = n_co, stringsAsFactors = FALSE)
      has_both <- wide[[g1]] == 1 & wide[[g2]] == 1
      has_neither <- wide[[g1]] == 0 & wide[[g2]] == 0
      mut_samp <- rownames(wide)[has_both]
      wt_samp <- rownames(wide)[has_neither]
      surv_uniq$group <- ifelse(surv_uniq$Sample %in% mut_samp, "co_mut", ifelse(surv_uniq$Sample %in% wt_samp, "WT", NA))
      surv_sub <- surv_uniq[!is.na(surv_uniq$group), , drop = FALSE]
      if (sum(surv_sub$group == "co_mut") >= 10 && sum(surv_sub$group == "WT") >= 10) {
        m <- tryCatch(coxph(Surv(Time_to_OS, as.numeric(Censor)) ~ factor(group, levels = c("WT", "co_mut")), data = surv_sub), error = function(e) NULL)
        if (!is.null(m)) {
          ci <- confint(m)
          hr_list[[length(hr_list) + 1]] <- data.frame(gene1 = g1, gene2 = g2, HR = exp(coef(m))[1], lower = exp(ci[1]), upper = exp(ci[2]),
            p_hr = summary(m)$coefficients[1, 5], n_cooccur = n_co, stringsAsFactors = FALSE)
        }
      }
    }
    or_df <- if (length(or_list) > 0) { o <- do.call(rbind, or_list); o$q_or <- p.adjust(o$p_or, "fdr"); o } else NULL
    hr_df <- if (length(hr_list) > 0) { h <- do.call(rbind, hr_list); h$q_hr <- p.adjust(h$p_hr, "fdr"); h } else NULL
    merged <- NULL
    if (!is.null(or_df) && !is.null(hr_df)) {
      merged <- merge(or_df[, c("gene1", "gene2", "odds_ratio", "q_or", "n_cooccur")],
        hr_df[, c("gene1", "gene2", "HR", "q_hr", "lower", "upper")],
        by = c("gene1", "gene2"), all = TRUE)
    }
    list(or = or_df, hr = hr_df, merged = merged, wide = wide, surv = surv_uniq, df = df, genes = genes)
  })

  fig2_triple_data <- reactive({
    fd <- fig2_data()
    if (is.null(fd) || length(fd$genes) < 3) return(NULL)
    df <- fd$df
    wide <- fd$wide
    surv_uniq <- fd$surv
    genes <- fd$genes
    triples <- t(combn(genes, 3))
    res <- list()
    wide_cols <- colnames(wide)
    for (i in seq_len(nrow(triples))) {
      g1 <- triples[i, 1]; g2 <- triples[i, 2]; g3 <- triples[i, 3]
      if (!all(c(g1, g2, g3) %in% wide_cols)) next
      has_3 <- wide[[g1]] == 1 & wide[[g2]] == 1 & wide[[g3]] == 1
      has_2 <- (wide[[g1]] + wide[[g2]] + wide[[g3]]) == 2
      n3 <- sum(has_3)
      if (n3 < 10) next
      samp_3 <- rownames(wide)[has_3]
      samp_2 <- rownames(wide)[has_2]
      surv_uniq$grp <- ifelse(surv_uniq$Sample %in% samp_3, "triple", ifelse(surv_uniq$Sample %in% samp_2, "double", NA))
      surv_sub <- surv_uniq[surv_uniq$grp %in% c("triple", "double"), , drop = FALSE]
      if (nrow(surv_sub) < 20 || sum(surv_sub$grp == "triple") < 10) next
      m <- tryCatch(coxph(Surv(Time_to_OS, as.numeric(Censor)) ~ factor(grp, levels = c("double", "triple")), data = surv_sub), error = function(e) NULL)
      if (is.null(m)) next
      ci <- confint(m)
      pv <- summary(m)$coefficients[1, 5]
      genotype <- paste(g1, g2, g3, sep = " + ")
      res[[length(res) + 1]] <- data.frame(genotype = genotype, g1 = g1, g2 = g2, g3 = g3, n_3 = n3, HR = exp(coef(m))[1], lower = exp(ci[1]), upper = exp(ci[2]), p = pv, stringsAsFactors = FALSE)
    }
    if (length(res) == 0) return(NULL)
    out <- do.call(rbind, res)
    out$q <- p.adjust(out$p, "fdr")
    out
  })

  output$fig2a_plot <- renderPlot({
    fd <- fig2_data()
    if (is.null(fd$or) || nrow(fd$or) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    d <- fd$or
    d$sig <- "NS"
    d$sig[d$q_or < 0.05 & d$odds_ratio > 1] <- "Co-occurring"
    d$sig[d$q_or < 0.05 & d$odds_ratio < 1] <- "Mutually exclusive"
    d$neglog10q <- -log10(pmax(d$q_or, 1e-10))
    d$logOR <- log(pmax(d$odds_ratio, 0.01))
    ggplot(d, aes(x = logOR, y = neglog10q, color = sig, size = n_cooccur)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70") +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Co-occurring" = "#b2182b", "Mutually exclusive" = "#2166ac", "NS" = "gray70")) +
      scale_size_continuous(range = c(2, 8)) +
      labs(x = "log(Odds Ratio)", y = "-log10(q-value)", color = NULL, size = "n co-mut") +
      theme_minimal(base_size = 11)
  })

  output$fig2b_plot <- renderPlot({
    fd <- fig2_data()
    if (is.null(fd$hr) || nrow(fd$hr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    d <- fd$hr
    d$sig <- "NS"
    d$sig[d$q_hr < 0.05 & d$HR <= 1] <- "Favorable"
    d$sig[d$q_hr < 0.05 & d$HR > 1] <- "Unfavorable"
    d$neglog10q <- -log10(pmax(d$q_hr, 1e-10))
    ggplot(d, aes(x = HR, y = neglog10q, color = sig, size = n_cooccur)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray70") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70") +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Favorable" = "#1b7837", "Unfavorable" = "#762a83", "NS" = "gray70")) +
      scale_size_continuous(range = c(2, 8)) +
      labs(x = "Hazard Ratio", y = "-log10(q-value)", color = NULL, size = "n co-mut") +
      theme_minimal(base_size = 11)
  })

  output$fig2c_plot <- renderPlot({
    fd <- fig2_data()
    if (is.null(fd$merged) || nrow(fd$merged) < 5) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data"))
    d <- fd$merged
    d <- d[complete.cases(d[, c("odds_ratio", "HR")]), , drop = FALSE]
    if (nrow(d) < 5) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient overlap"))
    d$logOR <- log(pmax(d$odds_ratio, 0.01))
    d$logHR <- log(pmax(d$HR, 0.01))
    d$sig_hr <- "NS"
    d$sig_hr[!is.na(d$q_hr) & d$q_hr < 0.05 & d$HR <= 1] <- "Favorable"
    d$sig_hr[!is.na(d$q_hr) & d$q_hr < 0.05 & d$HR > 1] <- "Unfavorable"
    d$n_cooccur[is.na(d$n_cooccur)] <- 10
    p <- ggplot(d, aes(x = logOR, y = logHR, color = sig_hr, size = n_cooccur)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Favorable" = "#1b7837", "Unfavorable" = "#762a83", "NS" = "gray70")) +
      scale_size_continuous(range = c(2, 8)) +
      labs(x = "log(Odds Ratio)", y = "log(Hazard Ratio)", color = NULL, size = "n co-mut") +
      theme_minimal(base_size = 11)
    for (grp in c("Favorable", "Unfavorable", "NS")) {
      sub <- d[d$sig_hr == grp, , drop = FALSE]
      if (nrow(sub) >= 3) {
        fit <- lm(logHR ~ logOR, data = sub)
        xs <- seq(min(d$logOR), max(d$logOR), length.out = 50)
        pred <- predict(fit, newdata = data.frame(logOR = xs), se.fit = TRUE)
        df_band <- data.frame(x = xs, y = pred$fit, ymin = pred$fit - 1.96 * pred$se.fit, ymax = pred$fit + 1.96 * pred$se.fit)
        col <- c("Favorable" = "#1b7837", "Unfavorable" = "#762a83", "NS" = "gray70")[grp]
        p <- p + geom_ribbon(data = df_band, aes(x = x, y = y, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = col, alpha = 0.2) +
          geom_line(data = df_band, aes(x = x, y = y), inherit.aes = FALSE, color = col, linewidth = 0.8)
      }
    }
    p
  })

  output$fig2d_plot <- renderPlot({
    tr <- fig2_triple_data()
    if (is.null(tr) || nrow(tr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    tr <- tr[order(-tr$n_3), , drop = FALSE]
    tr <- head(tr, 25)
    tr$sig_surv <- factor(tr$p <= 0.05, levels = c(FALSE, TRUE), labels = c("NS", "p ≤ 0.05"))
    levs <- rev(unique(as.character(tr$genotype)))
    tr$genotype <- factor(as.character(tr$genotype), levels = levs)
    ggplot(tr, aes(x = genotype, y = n_3, fill = sig_surv)) +
      geom_col() +
      geom_text(aes(label = n_3), hjust = -0.2, size = 3, color = "black") +
      scale_fill_manual(values = c("NS" = "gray70", "p ≤ 0.05" = "#b2182b")) +
      coord_flip() +
      labs(x = NULL, y = "# patients with 3 mutations", fill = "Survival") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "right")
  })

  output$fig2e_plot <- renderPlot({
    tr <- fig2_triple_data()
    if (is.null(tr) || nrow(tr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    sig <- tr[tr$q <= 0.15 & tr$p <= 0.05, , drop = FALSE]
    if (nrow(sig) == 0) sig <- head(tr[tr$p <= 0.05, , drop = FALSE], 6)
    if (nrow(sig) == 0) sig <- head(tr[order(tr$p), ], 6)
    if (nrow(sig) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No genotypes to display"))
    sig$col <- ifelse(sig$HR > 1, "Unfavorable", "Favorable")
    sig$genotype <- factor(sig$genotype, levels = sig$genotype[order(sig$HR)])
    p1 <- ggplot(sig, aes(x = genotype, y = HR, ymin = lower, ymax = pmin(upper, 8), color = col)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      geom_pointrange(size = 0.5) +
      coord_flip() +
      scale_color_manual(values = c("Favorable" = "#1b7837", "Unfavorable" = "#762a83")) +
      labs(x = NULL, y = "Hazard Ratio (3 vs 2 muts)", title = "Forest plot") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none")
    fd <- fig2_data()
    if (is.null(fd)) return(p1)
    if (has_survminer) {
      splots <- list()
      for (i in seq_len(min(3, nrow(sig)))) {
        gs <- unlist(strsplit(as.character(sig$genotype[i]), " + ", fixed = TRUE))
        gs <- gs[gs %in% colnames(fd$wide)]
        if (length(gs) < 3) next
        wide <- fd$wide
        surv <- fd$surv
        has_3 <- rowSums(wide[, gs, drop = FALSE]) == 3
        has_2 <- rowSums(wide[, gs, drop = FALSE]) == 2
        surv$grp <- ifelse(surv$Sample %in% rownames(wide)[has_3], "3 mut", ifelse(surv$Sample %in% rownames(wide)[has_2], "2 mut", NA))
        surv_sub <- surv[!is.na(surv$grp), , drop = FALSE]
        if (nrow(surv_sub) >= 10) {
          fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ grp, data = surv_sub)
          sd <- survdiff(Surv(Time_to_OS, as.numeric(Censor)) ~ grp, data = surv_sub)
          pv <- 1 - pchisq(sd$chisq, 1)
          pv_txt <- if (pv < 0.001) "p < 0.001" else paste0("p = ", format.pval(pv, digits = 2))
          sn <- names(fit$strata)
          scale_vals <- ifelse(grepl("3 mut", sn), "#b2182b", "gray50")
          scale_vals <- setNames(scale_vals, sn)
          legend_labs <- gsub("^grp=", "", sn)
          sp <- survminer::ggsurvplot(fit, data = surv_sub, risk.table = TRUE, pval = pv_txt,
            title = as.character(sig$genotype[i]), xlab = "Years", palette = scale_vals, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.05))
          splots[[length(splots) + 1]] <- sp
        }
      }
      if (length(splots) == 0) return(p1)
      km_grob <- survminer::arrange_ggsurvplots(splots, print = FALSE, ncol = 2)
      if (has_gridExtra) {
        gridExtra::grid.arrange(p1, km_grob, ncol = 1, heights = c(0.4, 0.6))
      } else {
        print(p1)
      }
    } else {
      plots <- list(p1)
      for (i in seq_len(min(3, nrow(sig)))) {
        gs <- unlist(strsplit(as.character(sig$genotype[i]), " + ", fixed = TRUE))
        gs <- gs[gs %in% colnames(fd$wide)]
        if (length(gs) < 3) next
        wide <- fd$wide
        surv <- fd$surv
        has_3 <- rowSums(wide[, gs, drop = FALSE]) == 3
        has_2 <- rowSums(wide[, gs, drop = FALSE]) == 2
        surv$grp <- ifelse(surv$Sample %in% rownames(wide)[has_3], "3 mut", ifelse(surv$Sample %in% rownames(wide)[has_2], "2 mut", NA))
        surv_sub <- surv[!is.na(surv$grp), , drop = FALSE]
        if (nrow(surv_sub) >= 10) {
          fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ grp, data = surv_sub)
          sd <- survdiff(Surv(Time_to_OS, as.numeric(Censor)) ~ grp, data = surv_sub)
          pv <- 1 - pchisq(sd$chisq, 1)
          pv_txt <- if (pv < 0.001) "p < 0.001" else paste0("p = ", format.pval(pv, digits = 2))
          sf <- data.frame(time = fit$time, surv = fit$surv, strata = rep(names(fit$strata), fit$strata))
          sn <- names(fit$strata)
          scale_vals <- setNames(ifelse(grepl("3 mut", sn), "#b2182b", "gray50"), sn)
          pkm <- ggplot(sf, aes(x = time, y = surv, color = strata)) + geom_step(linewidth = 1) +
            scale_color_manual(values = scale_vals) +
            labs(title = as.character(sig$genotype[i]), subtitle = pv_txt, x = "Years", y = "Survival") +
            theme_minimal(base_size = 9)
          plots[[length(plots) + 1]] <- pkm
        }
      }
      if (length(plots) == 1) return(p1)
      if (has_gridExtra) {
        gridExtra::grid.arrange(grobs = plots, ncol = 2)
      } else {
        p1
      }
    }
  })

  output$clinical_plot <- renderPlot({
    var <- input$clin_var
    req(df <- filtered_data())
    cols <- c("Sample", "Gene_for_analysis", var)
    df <- df[, cols[cols %in% colnames(df)], drop = FALSE]
    df <- df[!is.na(df[[var]]), , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    sample_vals <- unique(df[, c("Sample", var), drop = FALSE])
    sample_vals$y_val <- suppressWarnings(as.numeric(sample_vals[[var]]))
    sample_vals <- sample_vals[!is.na(sample_vals$y_val), , drop = FALSE]
    if (nrow(sample_vals) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No numeric data"))
    # Remove points outside requested ranges (match Single Gene Associations); then restrict df to in-range samples for stats
    if (var == "WBC") sample_vals <- sample_vals[sample_vals$y_val <= 200, , drop = FALSE]
    else if (var == "Hemoglobin") sample_vals <- sample_vals[sample_vals$y_val <= 15, , drop = FALSE]
    else if (var == "Platelet") sample_vals <- sample_vals[sample_vals$y_val <= 300, , drop = FALSE]
    if (nrow(sample_vals) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data in range"))
    df <- df[df$Sample %in% sample_vals$Sample, , drop = FALSE]
    gene_tbl <- table(as.character(df$Gene_for_analysis))
    # Apply gene frequency filter to determine which genes to show
    min_pct <- 1
    if (min_pct > 0) {
      n_samples <- length(unique(df$Sample))
      if (n_samples > 0) {
        gene_pct <- (gene_tbl / n_samples) * 100
        keep_genes <- names(gene_pct)[gene_pct >= min_pct]
        gene_tbl <- gene_tbl[names(gene_tbl) %in% keep_genes]
      }
    }
    top_genes <- names(sort(gene_tbl, decreasing = TRUE))
    plot_list <- lapply(top_genes, function(g) {
      mut_samples <- unique(df$Sample[df$Gene_for_analysis == g])
      data.frame(
        Gene = g,
        Sample = sample_vals$Sample,
        Status = ifelse(sample_vals$Sample %in% mut_samples, "mut", "WT"),
        value = sample_vals$y_val,
        stringsAsFactors = FALSE
      )
    })
    plot_df <- do.call(rbind, plot_list)
    med_mut <- vapply(top_genes, function(g) median(plot_df$value[plot_df$Gene == g & plot_df$Status == "mut"], na.rm = TRUE), numeric(1))
    gene_levels <- names(sort(med_mut, decreasing = TRUE))
    plot_df$Gene <- factor(plot_df$Gene, levels = gene_levels)
    y_range <- diff(range(plot_df$value, na.rm = TRUE))
    y_range <- if (y_range == 0) 1 else y_range
    sig_df <- do.call(rbind, lapply(gene_levels, function(g) {
      sub <- plot_df[plot_df$Gene == g, ]
      n_mut <- sum(sub$Status == "mut", na.rm = TRUE)
      n_wt <- sum(sub$Status == "WT", na.rm = TRUE)
      if (n_mut < 2 || n_wt < 2) {
        pval <- NA
      } else {
        pval <- tryCatch(wilcox.test(value ~ Status, data = sub, exact = FALSE)$p.value, error = function(e) NA)
      }
      sym <- if (is.na(pval)) "" else if (pval < 0.001) "***" else if (pval < 0.01) "**" else if (pval < 0.05) "*" else "ns"
      y_pos <- max(sub$value, na.rm = TRUE) + 0.04 * y_range
      data.frame(Gene = factor(g, levels = gene_levels), y = y_pos, label = sym, stringsAsFactors = FALSE)
    }))
    sig_df <- sig_df[!is.na(sig_df$Gene), , drop = FALSE]
    p <- ggplot(plot_df, aes(x = Gene, y = value, fill = Status)) +
      geom_boxplot(alpha = 0.8, outlier.alpha = 0.3, position = position_dodge(0.85), width = 0.7) +
      scale_fill_manual(values = c("mut" = "#8B0000", "WT" = "lightgrey"), name = "Status") +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
      labs(x = NULL, y = var) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(size = 13), legend.title = element_text(size = 14))
    if (nrow(sig_df) > 0) {
      p <- p + geom_text(data = sig_df, aes(x = Gene, y = y, label = label), inherit.aes = FALSE, size = 5, vjust = -0.2)
    }
    p
  })

  vaf_clonal_wide <- reactive({
    g1 <- input$vaf_gene1
    g2 <- input$vaf_gene2
    if (is.null(g1) || g1 == "" || is.null(g2) || g2 == "" || g1 == g2) return(NULL)
    req(df <- filtered_data())
    vcol <- vaf_col()
    if (!vcol %in% colnames(df)) vcol <- "VAF"
    df <- df[as.character(df$Gene_for_analysis) %in% c(g1, g2), c("Sample", "Gene_for_analysis", vcol), drop = FALSE]
    df$Gene <- df$Gene_for_analysis
    df$.metric <- df[[vcol]]
    df <- df[, c("Sample", "Gene", ".metric"), drop = FALSE]
    df <- df[!is.na(df$.metric) & df$.metric > 0, , drop = FALSE]
    wide <- reshape(df, idvar = "Sample", timevar = "Gene", direction = "wide")
    colnames(wide) <- gsub("^\\.metric\\.", "", colnames(wide))
    if (!g1 %in% colnames(wide) || !g2 %in% colnames(wide)) return(NULL)
    wide <- wide[complete.cases(wide[, c(g1, g2), drop = FALSE]), , drop = FALSE]
    if (nrow(wide) < 5) return(NULL)
    diff <- wide[[g1]] - wide[[g2]]
    wide$Clonality <- ifelse(diff > 5, paste(g1, "first"), ifelse(diff < -5, paste(g2, "first"), "Ambiguous"))
    wide
  })

  output$vaf_scatter_plot <- renderPlot({
    g1 <- input$vaf_gene1
    g2 <- input$vaf_gene2
    vcol <- vaf_col(); vlabel <- vaf_label()
    metric_short <- if (vcol == "CCF") "CCF" else "VAF"
    if (is.null(g1) || g1 == "" || is.null(g2) || g2 == "" || g1 == g2) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select two different genes", size = 6) +
        theme_void() + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA)))
    }
    wide <- vaf_clonal_wide()
    if (is.null(wide)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient overlapping samples", size = 5) +
      theme_void() + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA)))
    color_vals <- c("Ambiguous" = "#999999", setNames(c("#01665e", "#8c510a"), c(paste(g1, "first"), paste(g2, "first"))))
    ggplot(wide, aes(x = .data[[g2]], y = .data[[g1]], color = Clonality)) +
      geom_point(size = 4, alpha = 0.8) +
      geom_abline(slope = 1, intercept = 5, linetype = "dashed", color = "gray50") +
      geom_abline(slope = 1, intercept = -5, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = color_vals) +
      scale_x_continuous(limits = c(0, 100)) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = paste(g2, metric_short, "(%)"), y = paste(g1, metric_short, "(%)"), color = "Clonal order") +
      theme_minimal(base_size = 14) +
      theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA),
            axis.text = element_text(size = 13), axis.title = element_text(size = 14), legend.text = element_text(size = 12))
  })

  output$vaf_clonal_km_plot <- renderPlot({
    g1 <- input$vaf_gene1
    g2 <- input$vaf_gene2
    if (is.null(g1) || g1 == "" || is.null(g2) || g2 == "" || g1 == g2) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select two different genes", size = 6) +
        theme_void() + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA)))
    }
    wide <- vaf_clonal_wide()
    if (is.null(wide)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient overlapping samples", size = 5) +
      theme_void() + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA)))
    surv_df <- survival_data()
    if (is.null(surv_df) || nrow(surv_df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No survival data") + theme_void())
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    km_df <- merge(wide[, c("Sample", "Clonality"), drop = FALSE], surv_uniq, by = "Sample")
    km_df <- km_df[km_df$Clonality != "Ambiguous", , drop = FALSE]
    if (nrow(km_df) < 5 || length(unique(km_df$Clonality)) < 2) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient non-ambiguous samples for KM", size = 5) +
        theme_void() + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA)))
    }
    km_df$Clonality <- factor(km_df$Clonality)
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Clonality, data = km_df)
    if (has_survminer) {
      color_vals <- c(setNames(c("#01665e", "#8c510a"), c(paste(g1, "first"), paste(g2, "first"))))
      legend_labs <- gsub("^Clonality=", "", names(fit$strata))
      pal_vec <- unname(color_vals[legend_labs])
      p <- survminer::ggsurvplot(fit, data = km_df, risk.table = TRUE, pval = TRUE,
        title = "Survival by Clonal Order", xlab = "Years",
        palette = pal_vec, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.05))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install survminer for KM plots", size = 5) + theme_void()
    }
  })

  # Drug Sensitivity tab (BeatAML2)
  # Two distinct datasets: 2026 = all samples in inhibitor_auc; 2021 = waves 1+2 only (from clinical)
  beataml_full <- reactive({
    if (!exists("load_beataml2")) return(NULL)
    load_beataml2()
  })

  beataml_wave12 <- reactive({
    b <- beataml_full()
    if (is.null(b) || !b$ok) return(NULL)
    if (exists("subset_beataml_to_wave12")) return(subset_beataml_to_wave12(b))
    b
  })

  beataml_for_drug <- reactive({
    if (identical(input$main_nav, "meta_aml4")) return(beataml_full())
    return(beataml_wave12())
  })

  drug_correlations <- reactive({
    input$main_nav
    input$drug_subset
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    nav <- input$main_nav
    if (!is.null(precomputed_drug) && identical(nav, "meta_aml4") && subset %in% names(precomputed_drug$meta_aml4)) {
      pc <- precomputed_drug$meta_aml4[[subset]]$correlations
      if (!is.null(pc) && nrow(pc) > 0) return(pc)
    }
    if (!is.null(precomputed_drug) && identical(nav, "analyses") && subset %in% names(precomputed_drug$analyses)) {
      pc <- precomputed_drug$analyses[[subset]]$correlations
      if (!is.null(pc) && nrow(pc) > 0) return(pc)
    }
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(NULL)
    if (!exists("compute_drug_vaf_correlations")) return(NULL)
    compute_drug_vaf_correlations(b, subset = subset)
  })

  drug_loo <- reactive({
    input$main_nav
    input$drug_subset
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    nav <- input$main_nav
    if (!is.null(precomputed_drug) && identical(nav, "meta_aml4") && subset %in% names(precomputed_drug$meta_aml4)) {
      pc <- precomputed_drug$meta_aml4[[subset]]$loo
      if (!is.null(pc) && nrow(pc) > 0) return(pc)
    }
    if (!is.null(precomputed_drug) && identical(nav, "analyses") && subset %in% names(precomputed_drug$analyses)) {
      pc <- precomputed_drug$analyses[[subset]]$loo
      if (!is.null(pc) && nrow(pc) > 0) return(pc)
    }
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(NULL)
    if (!exists("compute_drug_vaf_loo")) return(NULL)
    compute_drug_vaf_loo(b, subset = subset)
  })

  # Gene and inhibitor options for VAF vs AUC Scatter from same sources as Mut vs WT and VAF vs AUC analyses
  drug_scatter_choices <- reactive({
    input$drug_subset
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(list(genes = character(0), inhibitors = character(0)))
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    genes <- character(0)
    inhibitors <- character(0)
    corr <- drug_correlations()
    mut_wt <- drug_mut_wt_all()
    if (!is.null(corr) && nrow(corr) > 0) {
      genes <- c(genes, as.character(corr$Gene))
      inhibitors <- c(inhibitors, as.character(corr$Inhibitor))
    }
    if (!is.null(mut_wt) && nrow(mut_wt) > 0) {
      genes <- c(genes, as.character(mut_wt$Gene))
      inhibitors <- c(inhibitors, as.character(mut_wt$Inhibitor))
    }
    if (length(genes) == 0 && "mutations" %in% names(b) && nrow(b$mutations) > 0) {
      genes <- sort(unique(as.character(b$mutations$Gene[b$mutations$Sample %in% allowed])))
    }
    if (length(inhibitors) == 0 && "auc" %in% names(b) && nrow(b$auc) > 0) {
      inhibitors <- sort(unique(as.character(b$auc$inhibitor[b$auc$Sample %in% allowed])))
    }
    list(genes = sort(unique(genes)), inhibitors = sort(unique(inhibitors)))
  })

  output$drug_scatter_selects_ui <- renderUI({
    ch <- drug_scatter_choices()
    sel_g <- input$drug_gene
    sel_i <- input$drug_inhibitor
    if (is.null(sel_g)) sel_g <- ""
    if (is.null(sel_i)) sel_i <- ""
    if (length(ch$genes) > 0 && !sel_g %in% ch$genes) sel_g <- ""
    if (length(ch$inhibitors) > 0 && !sel_i %in% ch$inhibitors) sel_i <- ""
    fluidRow(
      column(6, selectInput("drug_gene", "Select Gene", choices = c("Select..." = "", ch$genes), selected = sel_g)),
      column(6, selectInput("drug_inhibitor", "Select Inhibitor", choices = c("Select..." = "", ch$inhibitors), selected = sel_i))
    )
  })

  # Shared data and y-axis range for Individual Mutation and VAF Associations (boxplot + scatter)
  drug_scatter_shared <- reactive({
    drug <- input$drug_inhibitor
    gene <- input$drug_gene
    if (is.null(drug) || drug == "" || is.null(gene) || gene == "") return(NULL)
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(NULL)
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    auc_sub <- b$auc[b$auc$inhibitor == drug & b$auc$Sample %in% allowed, c("Sample", "auc"), drop = FALSE]
    if (nrow(auc_sub) == 0) return(NULL)
    mut_samples <- unique(b$mutations$Sample[b$mutations$Gene == gene & b$mutations$Sample %in% allowed])
    auc_sub$Status <- ifelse(auc_sub$Sample %in% mut_samples, "Mut", "WT")
    auc_sub$Status <- factor(auc_sub$Status, levels = c("WT", "Mut"))
    mut_auc <- auc_sub$auc[auc_sub$Status == "Mut"]
    wt_auc <- auc_sub$auc[auc_sub$Status == "WT"]
    y_lim <- range(auc_sub$auc, na.rm = TRUE)
    y_lim <- y_lim + c(-0.05, 0.05) * diff(y_lim)
    if (diff(y_lim) < 1e-6) y_lim <- y_lim + c(-5, 5)
    tt <- if (length(mut_auc) >= 2 && length(wt_auc) >= 2) tryCatch(t.test(mut_auc, wt_auc), error = function(e) NULL) else NULL
    wilcox <- if (length(mut_auc) >= 2 && length(wt_auc) >= 2) tryCatch(wilcox.test(mut_auc, wt_auc), error = function(e) NULL) else NULL
    median_mut <- median(mut_auc, na.rm = TRUE)
    median_wt <- median(wt_auc, na.rm = TRUE)
    mut_sub <- b$mutations[b$mutations$Gene == gene & b$mutations$Sample %in% allowed, c("Sample", "VAF"), drop = FALSE]
    merged <- merge(auc_sub[, c("Sample", "auc")], mut_sub, by = "Sample")
    list(auc_sub = auc_sub, merged = merged, y_lim = y_lim, tt = tt, wilcox = wilcox, median_mut = median_mut, median_wt = median_wt, mean_mut = mean(mut_auc, na.rm = TRUE), mean_wt = mean(wt_auc, na.rm = TRUE), n_mut = length(mut_auc), n_wt = length(wt_auc))
  })

  # Mut vs WT t-test for all gene-inhibitor pairs (n_mut >= 5 with AUC data) for Drug Sensitivity tab heatmap
  drug_mut_wt_all <- reactive({
    input$drug_subset
    input$main_nav
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    nav <- input$main_nav
    if (!is.null(precomputed_drug) && identical(nav, "meta_aml4") && subset %in% names(precomputed_drug$meta_aml4)) {
      pc <- precomputed_drug$meta_aml4[[subset]]$mut_wt
      if (!is.null(pc) && nrow(pc) > 0) return(pc)
    }
    if (!is.null(precomputed_drug) && identical(nav, "analyses") && subset %in% names(precomputed_drug$analyses)) {
      pc <- precomputed_drug$analyses[[subset]]$mut_wt
      if (!is.null(pc) && nrow(pc) > 0) return(pc)
    }
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok || !"auc" %in% names(b) || !"mutations" %in% names(b)) return(NULL)
    if (exists("compute_mut_wt_all")) return(compute_mut_wt_all(b, subset))
    # Fallback inline if helper not loaded: require >= 5 mutated samples with inhibitor AUC data per gene-inhibitor
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    auc_wide <- b$auc[b$auc$Sample %in% allowed, , drop = FALSE]
    mut_wide <- b$mutations[b$mutations$Sample %in% allowed, , drop = FALSE]
    if (nrow(auc_wide) == 0 || nrow(mut_wide) == 0) return(NULL)
    genes <- names(which(table(mut_wide$Gene) >= 5))
    inhibitors <- unique(auc_wide$inhibitor)
    if (length(genes) == 0 || length(inhibitors) == 0) return(NULL)
    res_list <- list()
    for (g in genes) {
      mut_samples <- unique(mut_wide$Sample[mut_wide$Gene == g])
      if (length(mut_samples) < 5) next
      for (inh in inhibitors) {
        sub <- auc_wide[auc_wide$inhibitor == inh, c("Sample", "auc"), drop = FALSE]
        sub$mut <- sub$Sample %in% mut_samples
        mut_auc <- sub$auc[sub$mut]
        wt_auc <- sub$auc[!sub$mut]
        n_mut <- length(mut_auc)
        n_wt <- length(wt_auc)
        if (n_mut < 5) next
        delta <- mean(mut_auc, na.rm = TRUE) - mean(wt_auc, na.rm = TRUE)
        tt <- tryCatch(t.test(mut_auc, wt_auc), error = function(e) NULL)
        pval <- if (!is.null(tt)) tt$p.value else NA_real_
        res_list[[length(res_list) + 1L]] <- data.frame(Gene = g, Inhibitor = inh, delta_AUC_mut_wt = delta, p_value = pval, n_mut = n_mut, n_wt = n_wt, stringsAsFactors = FALSE)
      }
    }
    if (length(res_list) == 0) return(NULL)
    out <- do.call(rbind, res_list)
    out$q_value <- p.adjust(out$p_value, method = "BH")
    out
  })

  # Drug correlations and LOOCV for gene summary tab (uses main subset filter); use precomputed when subset is All/de_novo/secondary
  gene_summary_drug_correlations <- reactive({
    input$subset
    input$main_nav
    subset <- input$subset
    if (is.null(subset) || subset == "All") subset <- "de_novo"
    subset_map <- c("De novo" = "de_novo", "Secondary" = "secondary", "Relapse" = "relapse", "Therapy" = "therapy", "Other" = "other")
    subset <- if (subset %in% names(subset_map)) subset_map[subset] else tolower(subset)
    nav <- input$main_nav
    if (!is.null(precomputed_drug) && subset %in% c("All", "de_novo", "secondary")) {
      if (identical(nav, "meta_aml4") && subset %in% names(precomputed_drug$meta_aml4)) {
        pc <- precomputed_drug$meta_aml4[[subset]]$correlations
        if (!is.null(pc) && nrow(pc) > 0) return(pc)
      }
      if (identical(nav, "analyses") && subset %in% names(precomputed_drug$analyses)) {
        pc <- precomputed_drug$analyses[[subset]]$correlations
        if (!is.null(pc) && nrow(pc) > 0) return(pc)
      }
    }
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(NULL)
    if (!exists("compute_drug_vaf_correlations")) return(NULL)
    compute_drug_vaf_correlations(b, subset = subset)
  })

  # Single Gene Drug: only the correlation rows for the selected gene (uses gene_summary_drug_subset)
  gene_summary_drug_correlations_gene <- reactive({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(NULL)
    subset <- input$gene_summary_drug_subset
    if (is.null(subset) || subset == "") subset <- "de_novo"
    nav <- input$main_nav
    base_genes <- gene_to_beataml_bases(g)
    full <- NULL
    if (!is.null(precomputed_drug) && subset %in% c("All", "de_novo", "secondary")) {
      if (identical(nav, "meta_aml4") && subset %in% names(precomputed_drug$meta_aml4)) {
        full <- precomputed_drug$meta_aml4[[subset]]$correlations
      } else if (identical(nav, "analyses") && subset %in% names(precomputed_drug$analyses)) {
        full <- precomputed_drug$analyses[[subset]]$correlations
      }
    }
    if (is.null(full) || nrow(full) == 0) {
      b <- beataml_for_drug()
      if (is.null(b) || !b$ok || !exists("compute_drug_vaf_correlations")) return(NULL)
      full <- compute_drug_vaf_correlations(b, subset = subset)
    }
    if (is.null(full) || nrow(full) == 0) return(NULL)
    sub <- full[full$Gene %in% base_genes, , drop = FALSE]
    gc()
    sub
  })

  gene_summary_drug_loo <- reactive({
    input$subset
    input$main_nav
    subset <- input$subset
    if (is.null(subset) || subset == "All") subset <- "de_novo"
    subset_map <- c("De novo" = "de_novo", "Secondary" = "secondary", "Relapse" = "relapse", "Therapy" = "therapy", "Other" = "other")
    subset <- if (subset %in% names(subset_map)) subset_map[subset] else tolower(subset)
    nav <- input$main_nav
    if (!is.null(precomputed_drug) && subset %in% c("All", "de_novo", "secondary")) {
      if (identical(nav, "meta_aml4") && subset %in% names(precomputed_drug$meta_aml4)) {
        pc <- precomputed_drug$meta_aml4[[subset]]$loo
        if (!is.null(pc) && nrow(pc) > 0) return(pc)
      }
      if (identical(nav, "analyses") && subset %in% names(precomputed_drug$analyses)) {
        pc <- precomputed_drug$analyses[[subset]]$loo
        if (!is.null(pc) && nrow(pc) > 0) return(pc)
      }
    }
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(NULL)
    if (!exists("compute_drug_vaf_loo")) return(NULL)
    compute_drug_vaf_loo(b, subset = subset)
  })

  output$drug_summary_ui <- renderUI({
    b <- beataml_for_drug()
    if (is.null(b)) return(p("BeatAML2 data not loaded. Run setup_beataml2.R to download."))
    if (!b$ok) return(p(style = "color:red", b$msg))
    NULL
  })

  output$drug_vaf_auc_summary_ui <- renderUI({
    input$main_nav
    input$drug_subset
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(NULL)
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    corr <- drug_correlations()
    n_sig_p <- if (!is.null(corr) && nrow(corr) > 0) sum(corr$p_value < 0.05) else 0
    n_sig_q <- if (!is.null(corr) && nrow(corr) > 0) sum(corr$q_value < 0.1) else 0
    wave_note <- if (identical(input$main_nav, "analyses")) " (waves 1+2)" else " (all waves)"
    p(
      "Samples with mutation + drug data: ", length(allowed), wave_note, " | ",
      "Drug-gene pairs analyzed: ", if (!is.null(corr)) nrow(corr) else 0, " | ",
      "Significant (p < 0.05): ", n_sig_p, " | ",
      "Significant (FDR q < 0.1): ", n_sig_q
    )
  })

  output$drug_summary_dotplot <- renderPlot({
    corr <- drug_correlations()
    if (is.null(corr) || nrow(corr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No correlations"))
    # Paper-style: R² >= 0.25, VAF_range >= 25, AUC_range >= 75
    plot_df <- corr[corr$R_squared >= 0.25 & corr$VAF_range >= 25 & corr$AUC_range >= 75, , drop = FALSE]
    if (nrow(plot_df) == 0) {
      plot_df <- corr[corr$p_value < 0.05, , drop = FALSE]
      if (nrow(plot_df) == 0) plot_df <- head(corr[order(-corr$R_squared), ], 30)
    }
    if (nrow(plot_df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No correlations meeting criteria"))
    plot_df$star <- ifelse(plot_df$q_value < 0.1, "*", "")
    gene_meds <- aggregate(delta_AUC ~ Gene, data = plot_df, FUN = median)
    plot_df$Gene <- factor(plot_df$Gene, levels = rev(gene_meds$Gene[order(gene_meds$delta_AUC)]))
    inh_meds <- aggregate(delta_AUC ~ Inhibitor, data = plot_df, FUN = median)
    plot_df$Inhibitor <- factor(plot_df$Inhibitor, levels = inh_meds$Inhibitor[order(inh_meds$delta_AUC)])
    p <- ggplot(plot_df, aes(x = Inhibitor, y = Gene, fill = delta_AUC, size = VAF_range, label = star)) +
      geom_point(shape = 21, color = "black", stroke = 0.5) +
      geom_text(size = 5, color = "#525252", hjust = -0.2, vjust = 0.5, show.legend = FALSE) +
      scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0,
        name = expression(Delta~"AUC")) +
      scale_size_continuous(name = expression(Delta~"VAF"), range = c(2, 8)) +
      labs(x = "Inhibitor", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 15), legend.text = element_text(size = 13), legend.title = element_text(size = 14),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "#f0f0f0"),
        legend.position = "right")
    gc()
    p
  })

  output$drug_mut_wt_summary_ui <- renderUI({
    input$drug_subset
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(NULL)
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    mut_wt <- drug_mut_wt_all()
    n_pairs <- if (!is.null(mut_wt) && nrow(mut_wt) > 0) nrow(mut_wt) else 0
    n_sig_p <- if (!is.null(mut_wt) && nrow(mut_wt) > 0) sum(mut_wt$p_value < 0.05, na.rm = TRUE) else 0
    n_sig_q <- if (!is.null(mut_wt) && nrow(mut_wt) > 0) sum(mut_wt$q_value < 0.1, na.rm = TRUE) else 0
    wave_note <- if (identical(input$main_nav, "analyses")) " (waves 1+2)" else " (all waves)"
    p(
      "Samples with mutation + drug data: ", length(allowed), wave_note, " | ",
      "Drug-gene pairs analyzed: ", n_pairs, " | ",
      "Significant (p < 0.05): ", n_sig_p, " | ",
      "Significant (FDR q < 0.1): ", n_sig_q
    )
  })

  output$drug_mut_wt_summary_table <- renderDT({
    input$drug_subset
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(datatable(data.frame(Inhibitor = character(), Gene = character(), `n mut` = integer(), `n WT` = integer(), `Mean Mut AUC` = numeric(), `Mean WT AUC` = numeric(), `Delta AUC` = numeric(), `p value` = character(), `q value` = character(), `Resistance Trend` = character(), stringsAsFactors = FALSE), options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE))
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    mut_wt <- drug_mut_wt_all()
    if (is.null(mut_wt) || nrow(mut_wt) == 0) return(datatable(data.frame(Inhibitor = character(), Gene = character(), `n mut` = integer(), `n WT` = integer(), `Mean Mut AUC` = numeric(), `Mean WT AUC` = numeric(), `Delta AUC` = numeric(), `p value` = character(), `q value` = character(), `Resistance Trend` = character(), stringsAsFactors = FALSE), options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE))
    # Filter to q < 0.1 (same as heatmap)
    if ("q_value" %in% colnames(mut_wt)) {
      mut_wt <- mut_wt[!is.na(mut_wt$q_value) & mut_wt$q_value < 0.1, , drop = FALSE]
    } else {
      mut_wt <- mut_wt[!is.na(mut_wt$p_value) & mut_wt$p_value < 0.05, , drop = FALSE]
    }
    if (nrow(mut_wt) == 0) return(datatable(data.frame(Inhibitor = character(), Gene = character(), `n mut` = integer(), `n WT` = integer(), `Mean Mut AUC` = numeric(), `Mean WT AUC` = numeric(), `Delta AUC` = numeric(), `p value` = character(), `q value` = character(), `Resistance Trend` = character(), stringsAsFactors = FALSE), options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE))
    # Compute mean AUC for mut and WT groups
    auc_wide <- b$auc[b$auc$Sample %in% allowed, , drop = FALSE]
    mut_wide <- b$mutations[b$mutations$Sample %in% allowed, , drop = FALSE]
    mean_mut_auc <- numeric(nrow(mut_wt))
    mean_wt_auc <- numeric(nrow(mut_wt))
    for (i in seq_len(nrow(mut_wt))) {
      g <- mut_wt$Gene[i]
      inh <- mut_wt$Inhibitor[i]
      mut_samples <- unique(mut_wide$Sample[mut_wide$Gene == g])
      sub_auc <- auc_wide[auc_wide$inhibitor == inh, c("Sample", "auc"), drop = FALSE]
      mut_auc_vals <- sub_auc$auc[sub_auc$Sample %in% mut_samples]
      wt_auc_vals <- sub_auc$auc[!sub_auc$Sample %in% mut_samples]
      mean_mut_auc[i] <- mean(mut_auc_vals, na.rm = TRUE)
      mean_wt_auc[i] <- mean(wt_auc_vals, na.rm = TRUE)
    }
    # Determine resistance trend: Resistant if delta > 0 (mut higher AUC), Sensitive if delta < 0 (mut lower AUC)
    resistance_trend <- ifelse(mut_wt$delta_AUC_mut_wt > 0, "Resistant", ifelse(mut_wt$delta_AUC_mut_wt < 0, "Sensitive", "Neutral"))
    df <- data.frame(
      Inhibitor = mut_wt$Inhibitor,
      Gene = mut_wt$Gene,
      `n mut` = as.integer(mut_wt$n_mut),
      `n WT` = as.integer(mut_wt$n_wt),
      `Mean Mut AUC` = round(mean_mut_auc, 2),
      `Mean WT AUC` = round(mean_wt_auc, 2),
      `Delta AUC` = round(mut_wt$delta_AUC_mut_wt, 3),
      `p value` = format_pval_display(mut_wt$p_value),
      `q value` = if ("q_value" %in% colnames(mut_wt)) format_pval_display(mut_wt$q_value) else "",
      `Resistance Trend` = resistance_trend,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (ncol(df) == 10 && "q value" %in% names(df) && all(df$`q value` == "")) df$`q value` <- NULL
    datatable(df, filter = "top", options = list(pageLength = 10000, lengthChange = FALSE, scrollX = TRUE), rownames = FALSE)
  })

  output$drug_mut_wt_heatmap <- renderPlot({
    mut_wt <- drug_mut_wt_all()
    if (is.null(mut_wt) || nrow(mut_wt) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No Mut vs WT data (need ≥5 mut samples with AUC per gene-inhibitor)") + theme_void())
    # Only show q < 0.1 (FDR) results
    if ("q_value" %in% colnames(mut_wt)) {
      mut_wt <- mut_wt[!is.na(mut_wt$q_value) & mut_wt$q_value < 0.1, , drop = FALSE]
    } else {
      mut_wt <- mut_wt[!is.na(mut_wt$p_value) & mut_wt$p_value < 0.05, , drop = FALSE]
    }
    if (nrow(mut_wt) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No gene-inhibitor pairs with q < 0.1 (FDR)") + theme_void())
    gene_meds <- aggregate(delta_AUC_mut_wt ~ Gene, data = mut_wt, FUN = median, na.rm = TRUE)
    mut_wt$Gene <- factor(mut_wt$Gene, levels = rev(gene_meds$Gene[order(gene_meds$delta_AUC_mut_wt)]))
    inh_meds <- aggregate(delta_AUC_mut_wt ~ Inhibitor, data = mut_wt, FUN = median, na.rm = TRUE)
    mut_wt$Inhibitor <- factor(mut_wt$Inhibitor, levels = inh_meds$Inhibitor[order(inh_meds$delta_AUC_mut_wt)])
    n_genes <- length(unique(mut_wt$Gene))
    n_inhibitors <- length(unique(mut_wt$Inhibitor))
    aspect_ratio <- n_genes / n_inhibitors
    rng <- range(mut_wt$delta_AUC_mut_wt, na.rm = TRUE)
    leg_breaks <- c(rng[1], 0, rng[2])
    leg_labels <- list("Sensitive", expression(Delta~AUC), "Resistant")
    p <- ggplot(mut_wt, aes(x = Inhibitor, y = Gene, fill = delta_AUC_mut_wt)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0,
        name = "Mutated\nsamples\nare more:", breaks = leg_breaks, labels = leg_labels) +
      labs(x = "Inhibitor", y = NULL) +
      theme_minimal(base_size = 18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        aspect.ratio = aspect_ratio, plot.margin = margin(0, 0, 0, 0))
    gc()
    p
  })

  output$drug_scatter_boxplot <- renderPlot({
    d <- drug_scatter_shared()
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Select Gene and Inhibitor", size = 5) +
      theme_void() +
      theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), aspect.ratio = 1)
    if (is.null(d)) return(empty_plot)
    auc_sub <- d$auc_sub
    gene <- input$drug_gene
    drug <- input$drug_inhibitor
    mut_col <- if (d$median_mut > d$median_wt) "#2166ac" else "#b2182b"
    wilcox_ptxt <- if (!is.null(d$wilcox)) { pv <- d$wilcox$p.value; if (pv < 0.001) "Wilcoxon p < 0.001" else paste0("Wilcoxon p = ", format.pval(pv, digits = 2)) } else ""
    stats_label <- paste0("WT: n = ", d$n_wt, "  |  Mut: n = ", d$n_mut, if (nchar(wilcox_ptxt) > 0) paste0("  |  ", wilcox_ptxt) else "")
    p <- ggplot(auc_sub, aes(x = Status, y = auc, fill = Status)) +
      geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
      scale_fill_manual(values = c("WT" = "#4D4D4D", "Mut" = mut_col), guide = "none") +
      coord_cartesian(ylim = d$y_lim) +
      labs(title = paste0(gene, " + ", drug), subtitle = stats_label, x = paste0(gene, " Status"), y = "Drug AUC") +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), plot.title = element_text(size = 13, face = "bold"), plot.subtitle = element_text(size = 10, color = "gray40"))
    gc()
    p
  })

  output$drug_scatter_plot <- renderPlot({
    d <- drug_scatter_shared()
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Select Gene and Inhibitor", size = 5) +
      theme_void() +
      theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"), aspect.ratio = 1)
    if (is.null(d)) return(empty_plot)
    merged <- d$merged
    if (nrow(merged) < 5) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient samples") + theme_minimal(base_size = 14))
    gene <- input$drug_gene
    drug <- input$drug_inhibitor
    fit <- lm(auc ~ VAF, data = merged)
    slope <- coef(fit)[2]
    line_col <- if (slope > 0) "#2166ac" else "#b2182b"
    r2 <- round(summary(fit)$r.squared, 3)
    pv <- summary(fit)$coefficients[2, 4]
    ptxt <- if (pv < 0.001) "p < 0.001" else paste0("p = ", format.pval(pv, digits = 2))
    p <- ggplot(merged, aes(x = VAF, y = auc)) +
      geom_point(size = 3, alpha = 0.8, color = line_col) +
      geom_smooth(method = "lm", se = TRUE, color = line_col, fill = line_col, alpha = 0.2) +
      coord_cartesian(ylim = d$y_lim) +
      labs(title = paste0(drug, " vs ", gene, " VAF"), x = "VAF (%)", y = "Drug AUC",
        subtitle = paste0("R² = ", r2, ", ", ptxt)) +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), plot.title = element_text(size = 13, face = "bold"), plot.subtitle = element_text(size = 10, color = "gray40"))
    gc()
    p
  })

  output$drug_correlation_table <- renderDT({
    corr <- drug_correlations()
    if (is.null(corr) || nrow(corr) == 0) return(datatable(data.frame(Inhibitor = character(), Gene = character(), Direction = character(), R_squared = numeric(), p_value = character(), q_value = character(), LOOCV_RMSE = numeric(), n = integer()), options = list(pageLength = 10000, lengthChange = FALSE)))
    cols <- c("Inhibitor", "Gene", "Direction", "R_squared", "p_value", "q_value", "LOOCV_RMSE", "LOOCV_MSE", "LOOCV_MSE_sd", "n", "VAF_range", "AUC_range")
    cols <- cols[cols %in% colnames(corr)]
    df <- corr[order(corr$p_value), cols, drop = FALSE]
    for (nm in c("R_squared", "LOOCV_RMSE", "LOOCV_MSE", "LOOCV_MSE_sd", "VAF_range", "AUC_range")) {
      if (nm %in% colnames(df) && is.numeric(df[[nm]])) df[[nm]] <- round(df[[nm]], 2)
    }
    if ("p_value" %in% colnames(df)) df$p_value <- format_pval_display(df$p_value)
    if ("q_value" %in% colnames(df)) df$q_value <- format_pval_display(df$q_value)
    gc()
    datatable(df, filter = "top", options = list(pageLength = 10000, lengthChange = FALSE))
  })

  output$drug_loo_heatmap <- renderPlot({
    loo <- drug_loo()
    if (is.null(loo) || nrow(loo) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No LOOCV data") + theme_void())
    loo <- loo[!is.na(loo$RMSE), , drop = FALSE]
    if (nrow(loo) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No RMSE values") + theme_void())
    gene_meds <- aggregate(RMSE ~ Gene, data = loo, FUN = median)
    loo$Gene <- factor(loo$Gene, levels = rev(gene_meds$Gene[order(gene_meds$RMSE)]))
    inh_meds <- aggregate(RMSE ~ Inhibitor, data = loo, FUN = median)
    loo$Inhibitor <- factor(loo$Inhibitor, levels = inh_meds$Inhibitor[order(inh_meds$RMSE)])
    p <- ggplot(loo, aes(x = Inhibitor, y = Gene, fill = RMSE)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_viridis_c(option = "viridis", name = "RMSE") +
      labs(x = "Inhibitor", y = "Gene") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text = element_text(size = 11), legend.title = element_text(size = 12))
    gc()
    p
  })

  output$download_drug_mut_wt_heatmap <- downloadHandler(
    filename = function() paste0("drug_mut_wt_heatmap_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      mut_wt <- drug_mut_wt_all()
      if (is.null(mut_wt) || nrow(mut_wt) == 0) { p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(); ggsave(file, plot = p, width = 8, height = 5, dpi = 150); return(invisible(NULL)) }
      if ("q_value" %in% colnames(mut_wt)) mut_wt <- mut_wt[!is.na(mut_wt$q_value) & mut_wt$q_value < 0.1, , drop = FALSE]
      else mut_wt <- mut_wt[!is.na(mut_wt$p_value) & mut_wt$p_value < 0.05, , drop = FALSE]
      if (nrow(mut_wt) == 0) { p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No pairs with q < 0.1") + theme_void(); ggsave(file, plot = p, width = 8, height = 5, dpi = 150); return(invisible(NULL)) }
      gene_meds <- aggregate(delta_AUC_mut_wt ~ Gene, data = mut_wt, FUN = median, na.rm = TRUE)
      mut_wt$Gene <- factor(mut_wt$Gene, levels = rev(gene_meds$Gene[order(gene_meds$delta_AUC_mut_wt)]))
      inh_meds <- aggregate(delta_AUC_mut_wt ~ Inhibitor, data = mut_wt, FUN = median, na.rm = TRUE)
      mut_wt$Inhibitor <- factor(mut_wt$Inhibitor, levels = inh_meds$Inhibitor[order(inh_meds$delta_AUC_mut_wt)])
      n_genes <- length(unique(mut_wt$Gene)); n_inhibitors <- length(unique(mut_wt$Inhibitor))
      aspect_ratio <- n_genes / n_inhibitors
      rng <- range(mut_wt$delta_AUC_mut_wt, na.rm = TRUE)
      p <- ggplot(mut_wt, aes(x = Inhibitor, y = Gene, fill = delta_AUC_mut_wt)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0, name = "Delta AUC") +
        labs(x = "Inhibitor", y = NULL) + theme_minimal(base_size = 18) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12), axis.text.y = element_text(size = 14), aspect.ratio = aspect_ratio, plot.margin = margin(0, 0, 0, 0))
      ggsave(file, plot = p, width = max(8, n_inhibitors * 0.4), height = max(6, n_genes * 0.35), dpi = 300)
    }
  )

  output$download_drug_summary_dotplot <- downloadHandler(
    filename = function() paste0("drug_vaf_auc_dotplot_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      corr <- drug_correlations()
      if (is.null(corr) || nrow(corr) == 0) { p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(); ggsave(file, plot = p, width = 8, height = 5, dpi = 150); return(invisible(NULL)) }
      plot_df <- corr[corr$R_squared >= 0.25 & corr$VAF_range >= 25 & corr$AUC_range >= 75, , drop = FALSE]
      if (nrow(plot_df) == 0) { plot_df <- corr[corr$p_value < 0.05, , drop = FALSE]; if (nrow(plot_df) == 0) plot_df <- head(corr[order(-corr$R_squared), ], 30) }
      if (nrow(plot_df) == 0) { p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No correlations") + theme_void(); ggsave(file, plot = p, width = 8, height = 5, dpi = 150); return(invisible(NULL)) }
      plot_df$star <- ifelse(plot_df$q_value < 0.1, "*", "")
      gene_meds <- aggregate(delta_AUC ~ Gene, data = plot_df, FUN = median)
      plot_df$Gene <- factor(plot_df$Gene, levels = rev(gene_meds$Gene[order(gene_meds$delta_AUC)]))
      inh_meds <- aggregate(delta_AUC ~ Inhibitor, data = plot_df, FUN = median)
      plot_df$Inhibitor <- factor(plot_df$Inhibitor, levels = inh_meds$Inhibitor[order(inh_meds$delta_AUC)])
      p <- ggplot(plot_df, aes(x = Inhibitor, y = Gene, fill = delta_AUC, size = VAF_range, label = star)) +
        geom_point(shape = 21, color = "black", stroke = 0.5) + geom_text(size = 5, color = "#525252", hjust = -0.2, vjust = 0.5, show.legend = FALSE) +
        scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0, name = "Delta AUC") + scale_size_continuous(name = "Delta VAF", range = c(2, 8)) +
        labs(x = "Inhibitor", y = NULL) + theme_minimal(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text = element_text(size = 14), legend.position = "right")
      ggsave(file, plot = p, width = 12, height = max(6, length(unique(plot_df$Gene)) * 0.3), dpi = 300)
    }
  )

  output$download_drug_scatter_plot <- downloadHandler(
    filename = function() paste0("drug_scatter_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      d <- drug_scatter_shared()
      if (is.null(d)) { p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select Gene and Inhibitor") + theme_void(); ggsave(file, plot = p, width = 6, height = 4, dpi = 150); return(invisible(NULL)) }
      merged <- d$merged
      if (nrow(merged) < 5) { p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient samples") + theme_minimal(); ggsave(file, plot = p, width = 6, height = 4, dpi = 150); return(invisible(NULL)) }
      gene <- input$drug_gene; drug <- input$drug_inhibitor
      line_col <- if (cor(merged$VAF, merged$auc, use = "pairwise.complete.obs") > 0) "#2166ac" else "#b2182b"
      fit <- lm(auc ~ VAF, data = merged)
      r2 <- round(summary(fit)$r.squared, 3); pv <- summary(fit)$coefficients[2, 4]
      ptxt <- if (pv < 0.001) "p < 0.001" else paste0("p = ", format.pval(pv, digits = 2))
      p <- ggplot(merged, aes(x = VAF, y = auc)) + geom_point(size = 3, alpha = 0.8, color = line_col) + geom_smooth(method = "lm", se = TRUE, color = line_col, fill = line_col, alpha = 0.2) + coord_cartesian(ylim = d$y_lim) + labs(title = paste0(drug, " vs ", gene, " VAF"), x = "VAF (%)", y = "Drug AUC", subtitle = paste0("R² = ", r2, ", ", ptxt)) + theme_minimal(base_size = 14)
      ggsave(file, plot = p, width = 6, height = 5, dpi = 300)
    }
  )

  output$download_drug_loo_heatmap <- downloadHandler(
    filename = function() paste0("drug_loo_rmse_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      loo <- drug_loo()
      if (is.null(loo) || nrow(loo) == 0) { p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(); ggsave(file, plot = p, width = 6, height = 4, dpi = 150); return(invisible(NULL)) }
      loo <- loo[!is.na(loo$RMSE), , drop = FALSE]
      if (nrow(loo) == 0) { p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No RMSE values") + theme_void(); ggsave(file, plot = p, width = 6, height = 4, dpi = 150); return(invisible(NULL)) }
      gene_meds <- aggregate(RMSE ~ Gene, data = loo, FUN = median)
      loo$Gene <- factor(loo$Gene, levels = rev(gene_meds$Gene[order(gene_meds$RMSE)]))
      inh_meds <- aggregate(RMSE ~ Inhibitor, data = loo, FUN = median)
      loo$Inhibitor <- factor(loo$Inhibitor, levels = inh_meds$Inhibitor[order(inh_meds$RMSE)])
      p <- ggplot(loo, aes(x = Inhibitor, y = Gene, fill = RMSE)) + geom_tile(color = "white", linewidth = 0.3) + scale_fill_viridis_c(option = "viridis", name = "RMSE") + labs(x = "Inhibitor", y = "Gene") + theme_minimal(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      ggsave(file, plot = p, width = max(8, length(unique(loo$Inhibitor)) * 0.35), height = max(5, length(unique(loo$Gene)) * 0.3), dpi = 300)
    }
  )

  output$download_drug_mut_wt_table <- downloadHandler(
    filename = function() paste0("drug_mut_wt_table_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
    content = function(file) {
      b <- beataml_for_drug()
      if (is.null(b) || !b$ok) return(write("No Beat AML data", file))
      mut_wt <- drug_mut_wt_all()
      if (is.null(mut_wt) || nrow(mut_wt) == 0) return(write("No data", file))
      if ("q_value" %in% colnames(mut_wt)) mut_wt <- mut_wt[!is.na(mut_wt$q_value) & mut_wt$q_value < 0.1, , drop = FALSE]
      else mut_wt <- mut_wt[!is.na(mut_wt$p_value) & mut_wt$p_value < 0.05, , drop = FALSE]
      if (nrow(mut_wt) == 0) return(write("No pairs with q < 0.1", file))
      subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
      allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
      auc_wide <- b$auc[b$auc$Sample %in% allowed, , drop = FALSE]
      mut_wide <- b$mutations[b$mutations$Sample %in% allowed, , drop = FALSE]
      mean_mut_auc <- numeric(nrow(mut_wt)); mean_wt_auc <- numeric(nrow(mut_wt))
      for (i in seq_len(nrow(mut_wt))) {
        g <- mut_wt$Gene[i]; inh <- mut_wt$Inhibitor[i]
        mut_samples <- unique(mut_wide$Sample[mut_wide$Gene == g])
        sub_auc <- auc_wide[auc_wide$inhibitor == inh, c("Sample", "auc"), drop = FALSE]
        mean_mut_auc[i] <- mean(sub_auc$auc[sub_auc$Sample %in% mut_samples], na.rm = TRUE)
        mean_wt_auc[i] <- mean(sub_auc$auc[!sub_auc$Sample %in% mut_samples], na.rm = TRUE)
      }
      resistance_trend <- ifelse(mut_wt$delta_AUC_mut_wt > 0, "Resistant", ifelse(mut_wt$delta_AUC_mut_wt < 0, "Sensitive", "Neutral"))
      out <- data.frame(Inhibitor = mut_wt$Inhibitor, Gene = mut_wt$Gene, n_mut = as.integer(mut_wt$n_mut), n_WT = as.integer(mut_wt$n_wt), Mean_Mut_AUC = round(mean_mut_auc, 2), Mean_WT_AUC = round(mean_wt_auc, 2), Delta_AUC = round(mut_wt$delta_AUC_mut_wt, 3), p_value = mut_wt$p_value, q_value = if ("q_value" %in% colnames(mut_wt)) mut_wt$q_value else NA, Resistance_Trend = resistance_trend, stringsAsFactors = FALSE)
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_drug_correlation_table <- downloadHandler(
    filename = function() paste0("drug_vaf_auc_correlations_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
    content = function(file) {
      corr <- drug_correlations()
      if (is.null(corr) || nrow(corr) == 0) return(write("No data", file))
      cols <- c("Inhibitor", "Gene", "Direction", "R_squared", "p_value", "q_value", "LOOCV_RMSE", "n")
      cols <- cols[cols %in% colnames(corr)]
      write.csv(corr[, cols, drop = FALSE], file, row.names = FALSE)
    }
  )

  # ---- All Gene Associations tab: comprehensive summary for a single gene ----
  gene_summary_data <- reactive({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(NULL)
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    samples_with_gene <- unique(df$Sample[as.character(df$Gene_for_analysis) == g])
    df[df$Sample %in% samples_with_gene, , drop = FALSE]
  })

  # Fetch NCBI Gene Summary
  gene_ncbi_summary_text <- reactive({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(NULL)
    # Extract base gene name (remove suffixes like _ITD, _TKD, _mono, _bi)
    g_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", g)
    if (!grepl("^[A-Z0-9]+$", g_base)) g_base <- g
    
    # Try to fetch from NCBI using rentrez if available
    if (has_rentrez) {
      tryCatch({
        # Search for gene
        search_res <- rentrez::entrez_search(db = "gene", term = paste0(g_base, "[Gene] AND Homo sapiens[Organism]"), retmax = 1)
        if (length(search_res$ids) > 0) {
          # Get gene summary
          gene_summary <- rentrez::entrez_summary(db = "gene", id = search_res$ids[1])
          if (!is.null(gene_summary$summary)) {
            return(gene_summary$summary)
          }
        }
      }, error = function(e) {
        return(NULL)
      })
    }
    return(NULL)
  })
  
  output$gene_ncbi_summary <- renderUI({
    summary_text <- gene_ncbi_summary_text()
    if (is.null(summary_text) || summary_text == "") {
      return(NULL)
    }
    return(div(
      style = "margin-top: 10px; padding: 10px; background-color: #f0f0f0; border-left: 3px solid #374e55;",
      p(strong("NCBI Gene Summary:"), style = "margin-bottom: 5px;"),
      p(summary_text, style = "font-size: 13px; line-height: 1.5; margin: 0;")
    ))
  })

  output$gene_summary_ui <- renderUI({
    if (!input$main_nav %in% c("analyses", "meta_aml4")) return(NULL)
    g <- input$gene_summary
    if (is.null(g) || g == "") {
      return(wellPanel(
        p(strong("Select a gene"), " from ", strong("Summarize all associations for a single gene"), " in the sidebar to view a comprehensive summary of all associations for that gene.")
      ))
    }
    tagList(
      wellPanel(
        h3("All associations for: ", g),
        uiOutput("gene_ncbi_summary")
      ),
      tabsetPanel(
        id = "gene_summary_sub_tabs",
        type = "tabs",
        selected = if (is.null(selected_gene_summary_sub_tab()) || !selected_gene_summary_sub_tab() %in% c("clinical", "comut", "vaf", "drug")) "clinical" else selected_gene_summary_sub_tab(),
        tabPanel("Clinical", value = "clinical"),
        tabPanel("Co-mutation", value = "comut"),
        tabPanel("Clonality", value = "vaf"),
        tabPanel("Drug Sensitivity", value = "drug")
      ),
      uiOutput("gene_summary_sub_content")
    )
  })

  # Only render the active sub-tab's content so only those outputs exist (reduces memory)
  output$gene_summary_sub_content <- renderUI({
    if (!input$main_nav %in% c("analyses", "meta_aml4")) return(NULL)
    g <- input$gene_summary
    sub <- input$gene_summary_sub_tabs
    if (is.null(sub) || !sub %in% c("clinical", "comut", "vaf", "drug")) sub <- "clinical"
    if (is.null(g) || g == "") return(NULL)
    df <- filtered_data()
    genes_other <- setdiff(sort(unique(as.character(df$Gene_for_analysis))), g)
    switch(sub,
      clinical = tagList(
        wellPanel(
          fluidRow(column(8, h4("Clinical associations")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_clinical_plot", "Download plot")))),
          p(style = "font-size: 13px; color: #666; margin-bottom: 6px;", span(style = "color: #8B0000; font-weight: bold;", "Red"), " = mutated; ", span(style = "color: #808080; font-weight: bold;", "Gray"), " = WT. *** p<0.001, ** p<0.01, * p<0.05, ns = not significant (Wilcoxon)."),
          plotOutput("gene_summary_clinical_plot", height = 260)),
        fluidRow(
          column(6, wellPanel(style = "min-height: 560px; height: 560px; overflow: hidden; box-sizing: border-box;",
            fluidRow(column(8, h4("Mutation distribution")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_lollipop_plot", "Download plot")))),
            p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", "Protein positions with mutations; height = count of mutations at that position."),
            plotOutput("gene_summary_lollipop_plot", height = 460))),
          column(6, wellPanel(style = "min-height: 560px; height: 560px; overflow: hidden; box-sizing: border-box;",
            fluidRow(column(8, h4("Kaplan-Meier: Single gene")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_survival", "Download plot")))),
            p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", span(style = "color: #8B0000; font-weight: bold;", "Red"), " = mutated; ", span(style = "color: #4D4D4D; font-weight: bold;", "Gray"), " = WT."),
            plotOutput("gene_summary_survival", height = 460)))
        )
      ),
      comut = tagList(
        fluidRow(
          column(9, wellPanel(style = "min-height: 560px; height: 560px; overflow: hidden; box-sizing: border-box;",
            fluidRow(column(8, h4("Oncoprint")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_oncoprint", label = "", icon = icon("download"))))),
            p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", if (is.null(g) || g == "") "Co-occurring mutations in samples with this gene mutated." else paste0("Co-occurring mutations in samples with ", g, " mutated.")),
            plotOutput("gene_summary_oncoprint", height = 460))),
          column(3, wellPanel(style = "min-height: 560px; height: 560px; overflow: hidden; box-sizing: border-box;",
            fluidRow(column(8, h4("Odds Ratio")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_comut_heatmap", label = "", icon = icon("download"))))),
            p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = mutually exclusive; ", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = co-occurring."),
            p(em("Calculated from all samples in the filtered subset, not just those included in the oncoprint."), style = "font-size: 13px; color: #666; margin-bottom: 8px;"),
            plotOutput("gene_summary_comut_heatmap", height = 442)))
        ),
        fluidRow(
          column(6, div(style = "min-height: 720px;", wellPanel(style = "height: 100%; box-sizing: border-box; min-height: 720px;",
            fluidRow(column(8, h4("Pairwise co-mutation survival")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_comut2_plot", "Download plot")))),
            selectInput("gene_summary_gene2", "Second gene:", choices = c("Select..." = "", genes_other)),
            plotOutput("gene_summary_comut2_plot", height = 600)))),
          column(6, div(style = "min-height: 720px;", wellPanel(style = "height: 100%; box-sizing: border-box; min-height: 720px;",
            fluidRow(column(8, h4("Triple co-mutation survival")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_comut3_plot", "Download plot")))),
            fluidRow(column(6, selectInput("gene_summary_gene3a", "Second gene:", choices = c("Select..." = "", genes_other))), column(6, selectInput("gene_summary_gene3b", "Third gene:", choices = c("Select..." = "", genes_other))),
            plotOutput("gene_summary_comut3_plot", height = 600)))))
        )
      ),
      vaf = fluidRow(
        column(6, wellPanel(
          fluidRow(column(8, h4("VAF / CCF distribution")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_vaf_plot", "Download plot")))),
          p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", "Toggle VAF/CCF with the sidebar control"),
          plotOutput("gene_summary_vaf_plot", height = 350))),
        column(6, wellPanel(
          fluidRow(column(8, h4("VAF / CCF and survival (MaxStat)")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_vaf_survival_plot", "Download plot")))),
          p(style = "font-size: 13px; color: #666; margin-bottom: 4px;", "High vs Low split at MaxStat optimal cutpoint; toggle VAF/CCF with sidebar control."),
          plotOutput("gene_summary_vaf_survival_plot", height = 350)))
      ),
      drug = (tagList(
        wellPanel(
          style = "background-color: #e9ecef;",
          h4("Drug sensitivity (Beat AML)"),
          p(style = "font-size: 14px; color: #333; margin-bottom: 8px;", "Associations between mutations and drug sensitivity were analysed using the Beat AML database (", tags$a("here", href = "https://biodev.github.io/BeatAML2/", target = "_blank"), ")."),
          div(style = "margin-bottom: 10px;", selectInput("gene_summary_drug_subset", "Select Beat AML cohort for analysis:", choices = c("All" = "All", "De novo" = "de_novo", "Secondary" = "secondary"), selected = "de_novo")),
          p(style = "font-size: 12px; color: #555; margin-bottom: 8px;",
            strong("Mut vs WT:"), " Gene–inhibitor pairs included if there were ≥5 mutated samples (and ≥5 WT) with AUC data in the selected cohort. Statistical test: two-sample t-test (mutant vs WT mean AUC); p-values and FDR (Benjamini–Hochberg) for resistance/sensitivity trend.",
            br(),
            strong("VAF vs AUC:"), " Gene–inhibitor pairs included if there were ≥5 samples with both VAF and AUC, and VAF range ≥25% and AUC range ≥75. Linear regression (AUC ~ VAF) was used to test the association; p-value and R² from the slope. Leave-one-out cross-validation (LOOCV) RMSE is reported as a check for overfitting. FDR (Benjamini–Hochberg) applied for multiple testing."),
          fluidRow(
            column(6, wellPanel(
              fluidRow(column(8, h5("Mut vs WT Inhibitor Sensitivity")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_drug_volcano", "Download plot")))),
              p(style = "font-size: 10px; color: #666; margin-bottom: 4px;", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = sensitive (mut lower AUC); ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = resistant (mut higher AUC)."),
              plotOutput("gene_summary_drug_volcano", height = 400))),
            column(6, wellPanel(
              fluidRow(column(8, h5("Mutation VAF vs Inhibitor Sensitivity")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_drug_dotplot", "Download plot")))),
              p(style = "font-size: 10px; color: #666; margin-bottom: 4px;", span(style = "color: #b2182b; font-weight: bold;", "Red"), " = negative VAF–AUC slope (sensitive trend); ", span(style = "color: #2166ac; font-weight: bold;", "Blue"), " = positive slope (resistance trend)."),
              plotOutput("gene_summary_drug_dotplot", height = 400)))),
          fluidRow(
            column(6, wellPanel(style = "min-height: 460px;",
              fluidRow(column(8, h5("Mut vs WT Results Summary")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_drug_mut_wt_table", "Download table")))),
              div(style = "height: 350px; overflow-y: auto;", DTOutput("gene_summary_drug_mut_wt_table")))),
            column(6, wellPanel(style = "min-height: 460px;",
              fluidRow(column(8, h5("VAF–AUC correlations Summary")), column(4, div(style = "text-align: right; margin-top: 5px;", downloadButton("download_gene_summary_drug_table", "Download table")))),
              div(style = "height: 350px; overflow-y: auto;", DTOutput("gene_summary_drug_table")))))
        )
      )),
    NULL
    )
  })

  gene_summary_oncoprint_data <- reactive({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(NULL)
    df <- gene_summary_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    samples <- unique(df$Sample)
    df_full <- filtered_data()
    df_sub <- df_full[df_full$Sample %in% samples, , drop = FALSE]
    if (!is_meta4()) {
      for (i in seq_len(nrow(df_sub))) {
        if (identical(as.character(df_sub$variant_type[i]), "other")) df_sub$variant_type[i] <- "Unknown"
      }
    }
    df_sub$Gene <- as.character(df_sub$Gene_for_analysis)
    df_sub$Gene_base <- df_sub$Gene
    has_suffix <- grepl("^([A-Z0-9]+)[_-]", df_sub$Gene)
    df_sub$Gene_base[has_suffix] <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", df_sub$Gene[has_suffix])
    flt3_itd_idx <- which(df_sub$Gene %in% c("FLT3-ITD", "FLT3_ITD") | (df_sub$Gene_base == "FLT3" & df_sub$variant_type == "ITD"))
    flt3_tkd_idx <- which(df_sub$Gene %in% c("FLT3-TKD", "FLT3_TKD") | (df_sub$Gene_base == "FLT3" & df_sub$variant_type %in% c("SNV", "Deletion")))
    flt3_other_idx <- which(df_sub$Gene %in% c("FLT3-Other", "FLT3_Other"))
    if (length(flt3_itd_idx) > 0) { df_sub$Gene_base[flt3_itd_idx] <- "FLT3"; df_sub$variant_type[flt3_itd_idx] <- "ITD" }
    if (length(flt3_tkd_idx) > 0) { df_sub$Gene_base[flt3_tkd_idx] <- "FLT3"; df_sub$variant_type[flt3_tkd_idx] <- "SNV" }
    if (length(flt3_other_idx) > 0) { df_sub$Gene_base[flt3_other_idx] <- "FLT3"; df_sub$variant_type[flt3_other_idx] <- "Other" }
    df_sub$Gene <- df_sub$Gene_base
    gene_tbl <- table(as.character(df_sub$Gene))
    top_genes <- names(sort(gene_tbl, decreasing = TRUE))[1:min(30, length(gene_tbl))]
    df_sub <- df_sub[df_sub$Gene %in% top_genes, , drop = FALSE]
    if (nrow(df_sub) == 0 || !"variant_type" %in% colnames(df_sub)) return(NULL)
    samples <- unique(as.character(df_sub$Sample))
    r <- length(top_genes)
    c <- length(samples)
    temp_dat <- as.data.frame(matrix(NA_character_, nrow = r, ncol = c))
    colnames(temp_dat) <- samples
    rownames(temp_dat) <- top_genes
    df_sub$variant_type <- as.character(df_sub$variant_type)
    for (i in seq_len(nrow(df_sub))) {
      pt <- df_sub$Sample[i]; gene <- df_sub$Gene[i]; var_type <- df_sub$variant_type[i]
      if (is.na(var_type) || var_type == "") var_type <- "Unknown"
      k <- match(pt, colnames(temp_dat)); l <- match(gene, rownames(temp_dat))
      if (!is.na(k) && !is.na(l)) {
        existing <- temp_dat[l, k]
        if (is.na(existing) || existing == "") temp_dat[l, k] <- var_type
        else {
          existing_types <- strsplit(existing, ";")[[1]]
          if (!var_type %in% existing_types) temp_dat[l, k] <- paste(c(existing_types, var_type), collapse = ";")
        }
      }
    }
    temp_dat[is.na(temp_dat)] <- ""
    temp_dat <- as.matrix(temp_dat)
    clin_cols <- c("Sample", "Cohort", "Sex", "Risk", "Subset", "Time_to_OS")
    clin_cols <- clin_cols[clin_cols %in% colnames(df_sub)]
    anno_df <- if (length(clin_cols) < 2) NULL else {
      a <- unique(df_sub[clin_cols])
      a <- a[!duplicated(a$Sample), , drop = FALSE]
      colnames(a)[colnames(a) == "Time_to_OS"] <- "Survival"
      if ("Survival" %in% colnames(a)) {
        a$Survival[!is.na(a$Survival) & a$Survival != ""] <- "Yes"
        a$Survival[is.na(a$Survival) | a$Survival == ""] <- "No"
      }
      a[match(samples, a$Sample), , drop = FALSE]
    }
    gc()
    list(temp_dat = temp_dat, anno_df = anno_df, samples = samples, genes = top_genes)
  })

  output$gene_summary_oncoprint <- renderPlot({
    req(input$main_nav %in% c("analyses", "meta_aml4"))
    od <- gene_summary_oncoprint_data()
    if (is.null(od) || !is.matrix(od$temp_dat) || nrow(od$temp_dat) == 0 || ncol(od$temp_dat) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    use_meta4 <- is_meta4()
    if (has_ComplexHeatmap) {
      temp_dat <- od$temp_dat
      anno_df <- od$anno_df
      samples <- od$samples
      genes <- od$genes
      col <- if (use_meta4) META4_VAR_COL else ONCO_VAR_COL
      vtypes <- unique(as.character(temp_dat[temp_dat != ""]))
      for (vt in setdiff(vtypes, names(col))) col[vt] <- "#80796BFF"
      alter_fun_list <- list(
        background = function(x, y, w, h) NULL,
        Deletion = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Deletion"], col = "#374E55FF")) },
        INDEL = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["INDEL"], col = "#DF8F44FF")) },
        Indel = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Indel"], col = "#DF8F44FF")) },
        Insertion = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Insertion"], col = "#00A1D5FF")) },
        ITD = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["ITD"], col = "#79AF97FF")) },
        SNV = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["SNV"], col = "#B24745FF")) },
        Splicing = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Splicing"], col = "#6A6599FF")) },
        PTD = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["PTD"], col = "tan")) },
        Other = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Other"], col = "#80796BFF")) },
        Unknown = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Unknown"], col = "#80796BFF")) }
      )
      for (vt in setdiff(vtypes, names(alter_fun_list))) alter_fun_list[[vt]] <- (function(fc) function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = fc, col = fc)) })(col[vt])
      ha <- NULL
      if (!is.null(anno_df) && nrow(anno_df) > 0) {
        Sex <- anno_df[, "Sex"]; Cohort <- anno_df[, "Cohort"]; Risk <- anno_df[, "Risk"]; Subset <- anno_df[, "Subset"]
        Survival <- if ("Survival" %in% colnames(anno_df)) anno_df[, "Survival"] else rep("No", nrow(anno_df))
        Sex[is.na(Sex)] <- "Unknown"; Cohort[is.na(Cohort)] <- "Unknown"; Risk[is.na(Risk)] <- "Unknown"; Subset[is.na(Subset)] <- "Unknown"
        cohort_col <- if (use_meta4) META4_COHORT_COL else ONCO_COHORT_COL
        subset_col <- if (use_meta4) META4_SUBSET_COL else ONCO_SUBSET_COL
        unique_cohorts <- unique(Cohort[!is.na(Cohort)])
        cohort_col_final <- character(length(unique_cohorts))
        names(cohort_col_final) <- unique_cohorts
        for (co in unique_cohorts) cohort_col_final[co] <- if (co %in% names(cohort_col)) cohort_col[co] else "gray70"
        unique_subsets <- unique(Subset[!is.na(Subset)])
        subset_col_final <- character(length(unique_subsets))
        names(subset_col_final) <- unique_subsets
        for (sub in unique_subsets) subset_col_final[sub] <- if (sub %in% names(subset_col)) subset_col[sub] else "gray70"
        ha <- ComplexHeatmap::HeatmapAnnotation(
          Sex = Sex, Cohort = Cohort, Risk = Risk, Subset = Subset, Survival = Survival,
          col = list(Sex = c("Male" = "#6a51a3", "Female" = "#43a2ca", "Unknown" = "gray90"), Survival = c("Yes" = "#252525", "No" = "#f0f0f0"), Risk = c("Adverse" = "#E64B35FF", "Intermediate" = "#8491B4FF", "Favorable" = "#00A087FF", "Unknown" = "#767676FF"), Subset = subset_col_final, Cohort = cohort_col_final),
          annotation_height = grid::unit(c(3.5, 3.5, 3.5, 3.5, 3.5), "mm"), show_annotation_name = TRUE,
          annotation_legend_param = list(Sex = list(title = "Sex"), Survival = list(title = "Survival"), Risk = list(title = "Risk"), Subset = list(title = "Subset"), Cohort = list(title = "Cohort"))
        )
      }
      fig <- ComplexHeatmap::oncoPrint(temp_dat, col = col, row_names_side = "right", bottom_annotation = ha, alter_fun_is_vectorized = TRUE, alter_fun = alter_fun_list)
      ComplexHeatmap::draw(fig)
    } else {
      temp_dat <- od$temp_dat
      samples <- od$samples
      genes <- od$genes
      rows <- list()
      for (i in seq_len(nrow(temp_dat))) {
        for (j in seq_len(ncol(temp_dat))) {
          v <- temp_dat[i, j]
          if (is.na(v) || v == "") next
          vts <- strsplit(as.character(v), ";")[[1]]
          rows[[length(rows) + 1L]] <- data.frame(Sample = colnames(temp_dat)[j], Gene = rownames(temp_dat)[i], variant_type = vts[1], stringsAsFactors = FALSE)
        }
      }
      mut_long <- if (length(rows) > 0) do.call(rbind, rows) else NULL
      if (is.null(mut_long) || nrow(mut_long) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
      PAL_VAR <- c(Deletion = "#374E55", Indel = "#DF8F44", ITD = "#79AF97", SNV = "#B24745", Splicing = "#6A6599", PTD = "tan", Other = "#80796B", Unknown = "#80796B")
      mut_long$Gene <- factor(mut_long$Gene, levels = rev(genes))
      mut_long$Sample <- factor(mut_long$Sample, levels = samples)
      mut_long$Sample_num <- as.numeric(mut_long$Sample)
      mut_long$Gene_num <- as.numeric(mut_long$Gene)
      bg <- expand.grid(Sample = factor(samples, levels = samples), Gene = factor(genes, levels = rev(genes)))
      bg$Sample_num <- as.numeric(bg$Sample)
      bg$Gene_num <- as.numeric(bg$Gene)
      p <- ggplot() + geom_tile(data = bg, aes(x = Sample_num, y = Gene_num), fill = "gray97", width = 1, height = 1) +
        geom_rect(data = mut_long, aes(xmin = Sample_num - 0.25, xmax = Sample_num + 0.25, ymin = Gene_num - 0.45, ymax = Gene_num + 0.45, fill = variant_type)) +
        scale_fill_manual(values = PAL_VAR, name = "Variant") + scale_x_continuous(limits = c(0.5, length(samples) + 0.5), expand = c(0, 0)) +
        scale_y_continuous(breaks = seq_along(genes), labels = rev(genes), expand = c(0, 0)) + labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 10) + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "right")
      print(p)
    }
  })

  output$gene_summary_comut_heatmap <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    cooc <- cooccurrence_data()
    if (is.null(cooc$matrix) || !g %in% rownames(cooc$matrix)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Gene not in co-occurrence matrix"))
    mat <- cooc$matrix
    others <- setdiff(colnames(mat), g)
    if (length(others) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No other genes"))
    # Filter to genes passing min. gene frequency for display only
    keep_genes <- filtered_genes()
    if (length(keep_genes) > 0) others <- intersect(others, keep_genes)
    if (length(others) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No other genes pass frequency filter"))
    vec <- mat[others, g]
    vec_log <- log2(vec + 0.01)
    vec_log[vec_log > 2] <- 2
    vec_log[vec_log < -2] <- -2
    ord <- order(vec_log, decreasing = TRUE)
    genes_ord <- names(vec_log)[ord]
    df_plot <- data.frame(Gene = factor(genes_ord, levels = rev(genes_ord)), log2OR = vec_log[ord])
    ggplot(df_plot, aes(x = 1, y = Gene, fill = log2OR)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
      labs(x = NULL, y = NULL, fill = "log2(OR)", title = paste("Co-occurrence\nwith", g)) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  })

  output$gene_summary_lollipop_plot <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void())
    df <- gene_summary_data()
    if (is.null(df) || nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No mutations for", g)) + theme_void())
    aa_col <- if ("HGVSp_Short" %in% colnames(df)) "HGVSp_Short" else if ("AA_change" %in% colnames(df)) "AA_change" else if ("Protein_Change" %in% colnames(df)) "Protein_Change" else if ("protein" %in% colnames(df)) "protein" else NULL
    if (is.null(aa_col)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No protein-change column (need AA_change / HGVSp_Short)") + theme_void())
    g_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", g)
    if (!grepl("^[A-Z0-9]+$", g_base)) g_base <- g
    norm <- normalize_protein_change_for_maf(df[[aa_col]])
    has_position <- !is.na(norm$position) & norm$hgvsp != "p.?"
    if (sum(has_position) < 1L) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No valid AA positions for lollipop") + theme_void())
    pos <- norm$position[has_position]
    pos_range <- range(pos, na.rm = TRUE)
    if (diff(pos_range) < 1) pos_range <- pos_range + c(-5, 5)
    pos_counts <- as.data.frame(table(Position = pos), stringsAsFactors = FALSE)
    pos_counts$Position <- as.integer(pos_counts$Position)
    pos_counts$Freq <- as.numeric(pos_counts$Freq)
    if (nrow(pos_counts) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No positions to plot") + theme_void())
    ggplot(pos_counts, aes(x = Position, y = Freq)) +
      geom_segment(aes(x = Position, xend = Position, y = 0, yend = Freq), linewidth = 0.8, color = "#374e55") +
      geom_point(size = 3, color = "#8B0000", fill = "#8B0000", shape = 21) +
      labs(title = paste("Mutation distribution:", g_base), x = "Protein position", y = "Count") +
      theme_minimal(base_size = 12) +
      scale_x_continuous(limits = pos_range + c(-1, 1) * 0.02 * diff(pos_range))
  })

  output$gene_summary_survival <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    surv_uniq <- gene_summary_survival_km_data()
    if (is.null(surv_uniq) || nrow(surv_uniq) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No survival data"))
    if (length(unique(surv_uniq$Mutation)) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient groups"))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Mutation, data = surv_uniq)
    if (has_survminer) {
      p <- survminer::ggsurvplot(fit, data = surv_uniq, risk.table = TRUE, pval = TRUE, title = paste("Survival by", g), xlab = "Years", palette = c("#8B0000", "#4D4D4D"), legend.labs = gsub("^Mutation=", "", names(fit$strata)), pval.coord = c(0, 0.05))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install survminer for KM plot") + theme_void()
    }
  })

  gene_summary_comut2_data <- reactive({
    g1 <- input$gene_summary
    g2 <- input$gene_summary_gene2
    if (is.null(g1) || g1 == "" || is.null(g2) || g2 == "") return(NULL)
    surv_df <- survival_data()
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    m1 <- unique(surv_df$Sample[surv_df$Gene_for_analysis == g1])
    m2 <- unique(surv_df$Sample[surv_df$Gene_for_analysis == g2])
    surv_uniq$Group <- ifelse(!surv_uniq$Sample %in% m1 & !surv_uniq$Sample %in% m2, "Neither",
      ifelse(surv_uniq$Sample %in% m1 & !surv_uniq$Sample %in% m2, paste0(g1, " only"),
        ifelse(!surv_uniq$Sample %in% m1 & surv_uniq$Sample %in% m2, paste0(g2, " only"), paste0(g1, " + ", g2))))
    surv_uniq
  })

  output$gene_summary_comut2_plot <- renderPlot({
    df <- gene_summary_comut2_data()
    if (is.null(df) || nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select second gene", size = 8) + theme_void() + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA)))
    if (length(unique(df$Group)) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Need at least 2 groups", size = 8) + theme_void() + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA)))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = df)
    if (has_survminer) {
      n_groups <- length(unique(df$Group))
      pal_all <- unname(c(PAL_SUBSET, PAL_MUT_CAT, PAL_COHORT))
      pal_vec <- pal_all[seq_len(n_groups)]
      if (length(pal_vec) < n_groups) pal_vec <- rep(pal_all, length.out = n_groups)
      p <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE, title = "Pairwise co-mutation survival", xlab = "Years", palette = pal_vec, legend.labs = gsub("^Group=", "", names(fit$strata)), pval.coord = c(0, 0.05))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install survminer") + theme_void()
    }
  })

  gene_summary_comut3_data <- reactive({
    g1 <- input$gene_summary
    g2 <- input$gene_summary_gene3a
    g3 <- input$gene_summary_gene3b
    if (is.null(g1) || g1 == "" || is.null(g2) || g2 == "" || is.null(g3) || g3 == "") return(NULL)
    surv_df <- survival_data()
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    m1 <- surv_uniq$Sample %in% unique(surv_df$Sample[surv_df$Gene_for_analysis == g1])
    m2 <- surv_uniq$Sample %in% unique(surv_df$Sample[surv_df$Gene_for_analysis == g2])
    m3 <- surv_uniq$Sample %in% unique(surv_df$Sample[surv_df$Gene_for_analysis == g3])
    n_mut <- m1 + m2 + m3
    surv_uniq$Group <- character(nrow(surv_uniq))
    surv_uniq$Group[n_mut == 0] <- "None"
    surv_uniq$Group[n_mut == 1 & m1] <- paste0(g1, " only")
    surv_uniq$Group[n_mut == 1 & m2] <- paste0(g2, " only")
    surv_uniq$Group[n_mut == 1 & m3] <- paste0(g3, " only")
    surv_uniq$Group[n_mut == 2 & m1 & m2] <- paste0(g1, " + ", g2)
    surv_uniq$Group[n_mut == 2 & m1 & m3] <- paste0(g1, " + ", g3)
    surv_uniq$Group[n_mut == 2 & m2 & m3] <- paste0(g2, " + ", g3)
    surv_uniq$Group[n_mut == 3] <- paste0(g1, " + ", g2, " + ", g3)
    surv_uniq
  })

  output$gene_summary_comut3_plot <- renderPlot({
    df <- gene_summary_comut3_data()
    if (is.null(df) || length(unique(df$Group)) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select second and third genes", size = 8) + theme_void() + theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA)))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = df)
    if (has_survminer) {
      n_groups <- length(unique(df$Group))
      pal_all <- unname(c(PAL_SUBSET, PAL_MUT_CAT, PAL_COHORT))
      pal_vec <- pal_all[seq_len(n_groups)]
      if (length(pal_vec) < n_groups) pal_vec <- rep(pal_all, length.out = n_groups)
      p <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE, title = "Triple co-mutation survival", xlab = "Years", palette = pal_vec, legend.labs = gsub("^Group=", "", names(fit$strata)), risk.table.height = 0.4, pval.coord = c(0, 0.05))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install survminer") + theme_void()
    }
  })

  output$gene_summary_clinical_plot <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    df_full <- filtered_data()
    if (is.null(df_full)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    clin_vars <- c("Age", "BM_blast_percent", "PB_blast_percent", "WBC", "Hemoglobin", "Platelet")
    clin_vars <- clin_vars[clin_vars %in% colnames(df_full)]
    if (length(clin_vars) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No clinical variables in data"))
    cols_needed <- c("Sample", "Gene_for_analysis", clin_vars)
    df <- df_full[!duplicated(df_full$Sample), cols_needed[cols_needed %in% colnames(df_full)], drop = FALSE]
    samples_with_g <- unique(df_full$Sample[as.character(df_full$Gene_for_analysis) == g])
    df_full <- NULL
    df_one <- df[, c("Sample", clin_vars[clin_vars %in% colnames(df)]), drop = FALSE]
    mut_lab <- paste0(g, " mut")
    df_one$Mutation <- ifelse(df_one$Sample %in% samples_with_g, mut_lab, "WT")
    long_list <- list()
    pval_list <- list()
    unit_labels <- c("WBC" = "WBC (1e-9/l)", "Age" = "Age (years)", "Hemoglobin" = "Hemoglobin (g/dl)", 
                     "Platelet" = "Platelet (1e-9/l)", "BM_blast_percent" = "BM blasts (%)", "PB_blast_percent" = "PB blasts (%)")
    for (cv in clin_vars) {
      df_one[[cv]] <- as.numeric(df_one[[cv]])
      sub <- df_one[!is.na(df_one[[cv]]), c("Sample", "Mutation", cv), drop = FALSE]
      if (cv == "WBC") sub <- sub[sub[[cv]] <= 200, , drop = FALSE]
      else if (cv == "Hemoglobin") sub <- sub[sub[[cv]] <= 15, , drop = FALSE]
      else if (cv == "Platelet") sub <- sub[sub[[cv]] <= 300, , drop = FALSE]
      if (nrow(sub) >= 5) {
        names(sub)[3] <- "value"
        sub$variable <- cv
        vlbl <- if (cv %in% names(unit_labels)) unit_labels[cv] else cv
        sub$variable_label <- vlbl
        long_list[[length(long_list) + 1L]] <- sub[, c("Sample", "Mutation", "variable", "variable_label", "value")]
        # Wilcoxon rank-sum test (Mut vs WT)
        wt_vals <- sub$value[sub$Mutation == "WT"]
        mut_vals <- sub$value[sub$Mutation == mut_lab]
        if (length(wt_vals) >= 2 && length(mut_vals) >= 2) {
          pv <- tryCatch(stats::wilcox.test(value ~ Mutation, data = sub, exact = FALSE)$p.value, error = function(e) NA_real_)
        } else {
          pv <- NA_real_
        }
        pval_list[[vlbl]] <- pv
      }
    }
    if (length(long_list) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for any variable"))
    long_df <- do.call(rbind, long_list)
    desired_order <- c("Age (years)", "BM blasts (%)", "PB blasts (%)", "WBC (1e-9/l)", "Hemoglobin (g/dl)", "Platelet (1e-9/l)")
    present_order <- intersect(desired_order, unique(long_df$variable_label))
    long_df$variable_label <- factor(long_df$variable_label, levels = present_order)
    # Significance labels: *** p<0.001, ** p<0.01, * p<0.05, ns otherwise
    sig_label <- function(p) {
      if (is.na(p)) return("ns")
      if (p < 0.001) return("***")
      if (p < 0.01) return("**")
      if (p < 0.05) return("*")
      return("ns")
    }
    sig_labels <- vapply(present_order, function(vl) sig_label(pval_list[[vl]]), character(1))
    # Annotation df: one row per facet, with y = max value for that variable (for placement above boxes)
    ann_df <- data.frame(
      variable_label = factor(present_order, levels = present_order),
      sig = sig_labels,
      y = vapply(present_order, function(vl) max(long_df$value[long_df$variable_label == vl], na.rm = TRUE), numeric(1)),
      stringsAsFactors = FALSE
    )
    # Nudge y slightly above box for label
    y_range_per_var <- vapply(present_order, function(vl) {
      r <- range(long_df$value[long_df$variable_label == vl], na.rm = TRUE)
      r[2] - r[1]
    }, numeric(1))
    ann_df$y <- ann_df$y + 0.08 * pmax(y_range_per_var, 1)
    pal_clin <- setNames(c("gray70", "#8B0000"), c("WT", mut_lab))
    ggplot(long_df, aes(x = Mutation, y = value, fill = Mutation)) +
      geom_boxplot(alpha = 0.8) +
      geom_text(data = ann_df, aes(x = 1.5, y = y, label = sig), inherit.aes = FALSE, size = 5, fontface = "bold") +
      scale_fill_manual(values = pal_clin) +
      facet_wrap(~ variable_label, scales = "free_y", nrow = 1) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none", strip.text = element_text(size = 14, face = "bold"))
  })

  output$gene_summary_vaf_plot <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    vcol <- vaf_col(); vlabel <- vaf_label()
    metric_short <- if (vcol == "CCF") "CCF" else "VAF"
    df <- filtered_data()
    if (!vcol %in% colnames(df)) vcol <- "VAF"
    df <- df[as.character(df$Gene_for_analysis) == g & !is.na(df[[vcol]]) & df[[vcol]] > 0, , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No", metric_short, "data")))
    df$.metric <- df[[vcol]]
    ggplot(df, aes(x = .metric)) +
      geom_histogram(bins = 25, fill = "#374e55", alpha = 0.8) +
      labs(x = paste0(metric_short, " (%)"), y = "Count", title = paste(metric_short, "distribution for", g)) +
      theme_minimal(base_size = 14)
  })

  output$gene_summary_vaf_survival_plot <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    vcol <- vaf_col()
    metric_short <- if (vcol == "CCF") "CCF" else "VAF"
    df <- filtered_data()
    if (!vcol %in% colnames(df)) vcol <- "VAF"
    df <- df[as.character(df$Gene_for_analysis) == g & !is.na(df[[vcol]]) & df[[vcol]] > 0, c("Sample", vcol), drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No", metric_short, "data")))
    colnames(df)[colnames(df) == vcol] <- ".metric"
    vaf_agg <- aggregate(.metric ~ Sample, data = df, max)
    surv_df <- survival_data()
    surv_df <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    merge_df <- merge(vaf_agg, surv_df, by = "Sample", all = FALSE)
    if (nrow(merge_df) < 20 || !has_maxstat) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste0("Need \u226520 samples with ", metric_short, " and survival; maxstat required")) + theme_void())
    mst <- tryCatch(maxstat::maxstat.test(Surv(Time_to_OS, as.numeric(Censor)) ~ .metric, data = merge_df, smethod = "LogRank", minprop = 0.2, maxprop = 0.8), error = function(e) NULL)
    if (is.null(mst)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "MaxStat failed") + theme_void())
    high_label <- paste0("High ", metric_short)
    low_label <- paste0("Low ", metric_short)
    merge_df$Group <- factor(ifelse(merge_df$.metric > mst$estimate, high_label, low_label), levels = c(high_label, low_label))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = merge_df)
    if (has_survminer) {
      legend_labs <- gsub("^Group=", "", names(fit$strata))
      p <- survminer::ggsurvplot(fit, data = merge_df, risk.table = TRUE, pval = TRUE, title = paste0(g, ": High vs Low ", metric_short, " (", round(mst$estimate, 1), "%)"), xlab = "Years", palette = c("#8B0000", "#4D4D4D"), legend.labs = legend_labs, pval.coord = c(0, 0.05))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install survminer") + theme_void()
    }
  })

  gene_summary_drug_mut_wt <- reactive({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(NULL)
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok || !"auc" %in% names(b) || !"mutations" %in% names(b)) return(NULL)
    subset <- input$gene_summary_drug_subset
    if (is.null(subset) || subset == "") subset <- "de_novo"
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    base_genes <- gene_to_beataml_bases(g)
    mut_samples <- unique(b$mutations$Sample[b$mutations$Gene %in% base_genes & b$mutations$Sample %in% allowed])
    inhibitors <- unique(b$auc$inhibitor)
    res_list <- list()
    for (inh in inhibitors) {
      sub <- b$auc[b$auc$inhibitor == inh & b$auc$Sample %in% allowed, c("Sample", "auc"), drop = FALSE]
      if (nrow(sub) < 5) next
      sub$mut <- sub$Sample %in% mut_samples
      mut_auc <- sub$auc[sub$mut]
      wt_auc <- sub$auc[!sub$mut]
      n_mut <- length(mut_auc)
      n_wt <- length(wt_auc)
      if (n_mut < 5) next
      delta <- mean(mut_auc, na.rm = TRUE) - mean(wt_auc, na.rm = TRUE)
      tt <- tryCatch(t.test(mut_auc, wt_auc), error = function(e) NULL)
      pval <- if (!is.null(tt)) tt$p.value else NA_real_
      res_list[[length(res_list) + 1L]] <- data.frame(Inhibitor = inh, delta_AUC_mut_wt = delta, p_value = pval, n_mut = n_mut, n_wt = n_wt, stringsAsFactors = FALSE)
    }
    if (length(res_list) == 0) return(NULL)
    out <- do.call(rbind, res_list)
    gc()
    out
  })

  output$gene_summary_drug_volcano <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void())
    df <- gene_summary_drug_mut_wt()
    if (is.null(df) || nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No Beat AML drug data or insufficient mut/WT samples") + theme_void())
    df <- df[!is.na(df$p_value) & df$p_value > 0, , drop = FALSE]
    df$neglog10p <- -log10(df$p_value)
    df$Category <- "Not significant"
    df$Category[df$p_value < 0.05 & df$delta_AUC_mut_wt < 0] <- "Sensitive"
    df$Category[df$p_value < 0.05 & df$delta_AUC_mut_wt > 0] <- "Resistant"
    df_sig <- df[df$Category != "Not significant", , drop = FALSE]
    # Limit labels: max 5 Resistant and max 5 Sensitive (by most significant / largest effect)
    label_sensitive <- df_sig[df_sig$Category == "Sensitive", , drop = FALSE]
    label_resistant <- df_sig[df_sig$Category == "Resistant", , drop = FALSE]
    if (nrow(label_sensitive) > 5) {
      label_sensitive <- label_sensitive[order(label_sensitive$p_value)[seq_len(5)], , drop = FALSE]
    }
    if (nrow(label_resistant) > 5) {
      label_resistant <- label_resistant[order(label_resistant$p_value)[seq_len(5)], , drop = FALSE]
    }
    df_label <- rbind(label_sensitive, label_resistant)
    p <- ggplot(df, aes(x = delta_AUC_mut_wt, y = neglog10p, color = Category, label = Inhibitor, size = n_mut)) +
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("Sensitive" = "#b2182b", "Resistant" = "#2166ac", "Not significant" = "gray50"), name = NULL) +
      scale_size_continuous(name = "N (mut)", range = c(2, 10)) +
      labs(x = "Mean AUC (mut) - Mean AUC (WT)", y = "-log10(p)", title = paste(g, ": Mut vs WT drug AUC")) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right") +
      guides(
        color = guide_legend(ncol = 1, override.aes = list(size = 5)),
        size = guide_legend(ncol = 1)
      )
    if (nrow(df_label) > 0) {
      if (has_ggrepel) {
        p <- p + ggrepel::geom_label_repel(data = df_label, aes(x = delta_AUC_mut_wt, y = neglog10p, label = Inhibitor), size = 5, show.legend = FALSE, max.overlaps = Inf)
      } else {
        p <- p + geom_text(data = df_label, aes(x = delta_AUC_mut_wt, y = neglog10p, label = Inhibitor), hjust = -0.1, size = 5, show.legend = FALSE)
      }
    }
    gc()
    p
  })

  output$gene_summary_drug_mut_wt_table <- renderDT({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(datatable(data.frame(Gene = character(), Inhibitor = character(), `n mut` = integer(), `n WT` = integer(), `Mean Mut AUC` = numeric(), `Mean WT AUC` = numeric(), `Delta AUC` = numeric(), `p value` = character(), `q value` = character(), `Resistance Trend` = character(), stringsAsFactors = FALSE), options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE))
    mut_wt <- gene_summary_drug_mut_wt()
    if (is.null(mut_wt) || nrow(mut_wt) == 0) return(datatable(data.frame(Gene = character(), Inhibitor = character(), `n mut` = integer(), `n WT` = integer(), `Mean Mut AUC` = numeric(), `Mean WT AUC` = numeric(), `Delta AUC` = numeric(), `p value` = character(), `q value` = character(), `Resistance Trend` = character(), stringsAsFactors = FALSE), options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE))
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(datatable(data.frame(Gene = character(), Inhibitor = character(), `n mut` = integer(), `n WT` = integer(), `Mean Mut AUC` = numeric(), `Mean WT AUC` = numeric(), `Delta AUC` = numeric(), `p value` = character(), `q value` = character(), `Resistance Trend` = character(), stringsAsFactors = FALSE), options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE))
    subset <- input$gene_summary_drug_subset
    if (is.null(subset) || subset == "") subset <- "de_novo"
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    base_genes <- gene_to_beataml_bases(g)
    mut_samples <- unique(b$mutations$Sample[b$mutations$Gene %in% base_genes & b$mutations$Sample %in% allowed])
    auc_wide <- b$auc[b$auc$Sample %in% allowed, , drop = FALSE]
    mean_mut_auc <- numeric(nrow(mut_wt))
    mean_wt_auc <- numeric(nrow(mut_wt))
    for (i in seq_len(nrow(mut_wt))) {
      inh <- mut_wt$Inhibitor[i]
      sub_auc <- auc_wide[auc_wide$inhibitor == inh, c("Sample", "auc"), drop = FALSE]
      mut_auc_vals <- sub_auc$auc[sub_auc$Sample %in% mut_samples]
      wt_auc_vals <- sub_auc$auc[!sub_auc$Sample %in% mut_samples]
      mean_mut_auc[i] <- mean(mut_auc_vals, na.rm = TRUE)
      mean_wt_auc[i] <- mean(wt_auc_vals, na.rm = TRUE)
    }
    q_val <- p.adjust(mut_wt$p_value, method = "BH")
    resistance_trend <- ifelse(mut_wt$delta_AUC_mut_wt > 0, "Resistant", ifelse(mut_wt$delta_AUC_mut_wt < 0, "Sensitive", "Neutral"))
    df <- data.frame(
      Gene = g,
      Inhibitor = mut_wt$Inhibitor,
      `n mut` = as.integer(mut_wt$n_mut),
      `n WT` = as.integer(mut_wt$n_wt),
      `Mean Mut AUC` = round(mean_mut_auc, 2),
      `Mean WT AUC` = round(mean_wt_auc, 2),
      `Delta AUC` = round(mut_wt$delta_AUC_mut_wt, 3),
      `p value` = format_pval_display(mut_wt$p_value),
      `q value` = format_pval_display(q_val),
      `Resistance Trend` = resistance_trend,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(pageLength = 10000, lengthChange = FALSE, scrollX = TRUE), rownames = FALSE)
  })

  output$gene_summary_drug_dotplot <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void())
    corr <- gene_summary_drug_correlations_gene()
    if (is.null(corr) || nrow(corr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No drug data (Beat AML)") + theme_void())
    if ("n" %in% colnames(corr)) corr <- corr[!is.na(corr$n) & corr$n >= 10, , drop = FALSE]
    base_genes <- gene_to_beataml_bases(g)
    sub <- corr[corr$Gene %in% base_genes, , drop = FALSE]
    if (nrow(sub) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No VAF–AUC correlations for", g)) + theme_void())
    sub <- sub[!is.na(sub$p_value) & sub$p_value < 0.05, , drop = FALSE]
    if (nrow(sub) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No VAF–AUC correlations with p < 0.05 for", g)) + theme_void())
    sub <- sub[order(sub$p_value), , drop = FALSE]
    inhibitors <- unique(sub$Inhibitor)
    if (length(inhibitors) > 24) inhibitors <- head(inhibitors, 24)
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No Beat AML data") + theme_void())
    subset <- input$gene_summary_drug_subset
    if (is.null(subset) || subset == "") subset <- "de_novo"
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    scatter_list <- list()
    trend_list <- list()
    for (inh in inhibitors) {
      auc_sub <- b$auc[b$auc$inhibitor == inh & b$auc$Sample %in% allowed, c("Sample", "auc"), drop = FALSE]
      mut_sub <- b$mutations[b$mutations$Gene %in% base_genes & b$mutations$Sample %in% allowed, c("Sample", "VAF"), drop = FALSE]
      merged <- merge(auc_sub, mut_sub, by = "Sample")
      if (nrow(merged) >= 5) {
        corr_row <- sub[sub$Inhibitor == inh, , drop = FALSE]
        if (nrow(corr_row) > 0) {
          fit <- tryCatch(lm(auc ~ VAF, data = merged), error = function(e) NULL)
          slope <- if (!is.null(fit)) coef(fit)[2] else NA_real_
          trend <- if (is.na(slope)) "Unknown" else if (slope < 0) "Sensitive" else "Resistant"
          merged$Inhibitor <- inh
          merged$R_squared <- corr_row$R_squared[1]
          merged$p_value <- corr_row$p_value[1]
          merged$q_value <- corr_row$q_value[1]
          merged$Trend <- trend
          scatter_list[[length(scatter_list) + 1L]] <- merged
          trend_list[[length(trend_list) + 1L]] <- data.frame(Inhibitor = inh, Trend = trend, slope = slope, stringsAsFactors = FALSE)
        }
      }
    }
    if (length(scatter_list) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Insufficient data for significant VAF–AUC correlations for", g)) + theme_void())
    scatter_df <- do.call(rbind, scatter_list)
    trend_df <- do.call(rbind, trend_list)
    scatter_df$Inhibitor <- factor(scatter_df$Inhibitor, levels = inhibitors)
    scatter_df$Trend <- factor(scatter_df$Trend, levels = c("Sensitive", "Resistant", "Unknown"))
    label_df <- unique(scatter_df[, c("Inhibitor", "R_squared", "p_value", "q_value", "Trend"), drop = FALSE])
    label_df$star <- ifelse(label_df$q_value < 0.1, "*", "")
    pval_txt <- format_pval_display(label_df$p_value)
    label_df$label <- paste0("R² = ", round(label_df$R_squared, 3), ", p ", pval_txt, label_df$star)
    trend_colors <- c("Sensitive" = "#b2182b", "Resistant" = "#2166ac", "Unknown" = "gray50")
    p <- ggplot(scatter_df, aes(x = VAF, y = auc, color = Trend)) +
      geom_point(size = 2.5, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.8) +
      scale_color_manual(values = trend_colors, name = "Trend", guide = "none") +
      facet_wrap(~ Inhibitor, scales = "free", ncol = min(3, length(inhibitors))) +
      labs(x = "VAF (%)", y = "Drug AUC", title = paste("VAF vs AUC scatterplots:", g, "(p < 0.05)")) +
      theme_minimal(base_size = 16) +
      theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "#f0f0f0", color = NA),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(size = 17, face = "bold")
      ) +
      geom_text(aes(x = Inf, y = Inf, label = label), hjust = 1.1, vjust = 1.5, size = 4, inherit.aes = FALSE, data = label_df)
    gc()
    p
  })

  output$gene_summary_drug_table <- renderDT({
    empty_drug <- datatable(data.frame(Gene = character(), Inhibitor = character(), Direction = character(), R_squared = numeric(), p_value = character(), q_value = character(), LOOCV_RMSE = numeric(), n = integer()), options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE)
    g <- input$gene_summary
    if (is.null(g) || g == "") return(empty_drug)
    corr <- gene_summary_drug_correlations_gene()
    if (is.null(corr) || nrow(corr) == 0) return(empty_drug)
    if ("n" %in% colnames(corr)) corr <- corr[!is.na(corr$n) & corr$n >= 10, , drop = FALSE]
    base_genes <- gene_to_beataml_bases(g)
    sub <- corr[corr$Gene %in% base_genes, , drop = FALSE]
    if (nrow(sub) == 0) return(empty_drug)
    cols <- c("Gene", "Inhibitor", "Direction", "R_squared", "p_value", "q_value", "LOOCV_RMSE", "LOOCV_MSE", "n", "VAF_range", "AUC_range")
    cols <- cols[cols %in% colnames(sub)]
    df <- sub[order(sub$p_value), cols, drop = FALSE]
    for (nm in c("R_squared", "LOOCV_RMSE", "LOOCV_MSE", "VAF_range", "AUC_range")) {
      if (nm %in% colnames(df) && is.numeric(df[[nm]])) df[[nm]] <- round(df[[nm]], 2)
    }
    if ("p_value" %in% colnames(df)) df$p_value <- format_pval_display(df$p_value)
    if ("q_value" %in% colnames(df)) df$q_value <- format_pval_display(df$q_value)
    gc()
    datatable(df, options = list(pageLength = 10000, lengthChange = FALSE), rownames = FALSE)
  })

  output$gene_summary_drug_loo_heatmap <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void())
    corr <- gene_summary_drug_correlations_gene()
    if (is.null(corr) || nrow(corr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No drug data (Beat AML)") + theme_void())
    if ("n" %in% colnames(corr)) corr <- corr[!is.na(corr$n) & corr$n >= 10, , drop = FALSE]
    base_genes <- gene_to_beataml_bases(g)
    sub <- corr[corr$Gene %in% base_genes & !is.na(corr$LOOCV_RMSE), , drop = FALSE]
    if (nrow(sub) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No LOOCV RMSE data for", g)) + theme_void())
    inh_ord <- sub$Inhibitor[order(sub$LOOCV_RMSE)]
    sub$Inhibitor <- factor(sub$Inhibitor, levels = inh_ord)
    p <- ggplot(sub, aes(x = Inhibitor, y = Gene, fill = LOOCV_RMSE)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_viridis_c(option = "viridis", name = "RMSE") +
      labs(x = "Inhibitor", y = NULL, title = paste("LOOCV RMSE:", g)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    gc()
    p
  })

  output$download_gene_summary_clinical_plot <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_clinical_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      g <- input$gene_summary
      if (is.null(g) || g == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void(), width = 8, height = 4, dpi = 150); return(invisible(NULL)) }
      df_full <- filtered_data()
      if (is.null(df_full)) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void(), width = 8, height = 4, dpi = 150); return(invisible(NULL)) }
      clin_vars <- c("Age", "BM_blast_percent", "PB_blast_percent", "WBC", "Hemoglobin", "Platelet")
      clin_vars <- clin_vars[clin_vars %in% colnames(df_full)]
      if (length(clin_vars) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No clinical variables") + theme_void(), width = 8, height = 4, dpi = 150); return(invisible(NULL)) }
      cols_needed <- c("Sample", "Gene_for_analysis", clin_vars)
      df <- df_full[!duplicated(df_full$Sample), cols_needed[cols_needed %in% colnames(df_full)], drop = FALSE]
      samples_with_g <- unique(df_full$Sample[as.character(df_full$Gene_for_analysis) == g])
      df_one <- df[, c("Sample", clin_vars[clin_vars %in% colnames(df)]), drop = FALSE]
      mut_lab <- paste0(g, " mut")
      df_one$Mutation <- ifelse(df_one$Sample %in% samples_with_g, mut_lab, "WT")
      long_list <- list(); pval_list <- list()
      unit_labels <- c("WBC" = "WBC (1e-9/l)", "Age" = "Age (years)", "Hemoglobin" = "Hemoglobin (g/dl)", "Platelet" = "Platelet (1e-9/l)", "BM_blast_percent" = "BM blasts (%)", "PB_blast_percent" = "PB blasts (%)")
      for (cv in clin_vars) {
        df_one[[cv]] <- as.numeric(df_one[[cv]])
        sub <- df_one[!is.na(df_one[[cv]]), c("Sample", "Mutation", cv), drop = FALSE]
        if (cv == "WBC") sub <- sub[sub[[cv]] <= 200, , drop = FALSE]
        else if (cv == "Hemoglobin") sub <- sub[sub[[cv]] <= 15, , drop = FALSE]
        else if (cv == "Platelet") sub <- sub[sub[[cv]] <= 300, , drop = FALSE]
        if (nrow(sub) >= 5) {
          names(sub)[3] <- "value"; sub$variable <- cv
          vlbl <- if (cv %in% names(unit_labels)) unit_labels[cv] else cv
          sub$variable_label <- vlbl
          long_list[[length(long_list) + 1L]] <- sub[, c("Sample", "Mutation", "variable", "variable_label", "value")]
          pv <- if (length(sub$value[sub$Mutation == "WT"]) >= 2 && length(sub$value[sub$Mutation == mut_lab]) >= 2) tryCatch(stats::wilcox.test(value ~ Mutation, data = sub, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
          pval_list[[vlbl]] <- pv
        }
      }
      if (length(long_list) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") + theme_void(), width = 8, height = 4, dpi = 150); return(invisible(NULL)) }
      long_df <- do.call(rbind, long_list)
      desired_order <- c("Age (years)", "BM blasts (%)", "PB blasts (%)", "WBC (1e-9/l)", "Hemoglobin (g/dl)", "Platelet (1e-9/l)")
      present_order <- intersect(desired_order, unique(long_df$variable_label))
      long_df$variable_label <- factor(long_df$variable_label, levels = present_order)
      sig_label <- function(p) { if (is.na(p)) return("ns"); if (p < 0.001) return("***"); if (p < 0.01) return("**"); if (p < 0.05) return("*"); return("ns") }
      sig_labels <- vapply(present_order, function(vl) sig_label(pval_list[[vl]]), character(1))
      ann_df <- data.frame(variable_label = factor(present_order, levels = present_order), sig = sig_labels, y = vapply(present_order, function(vl) max(long_df$value[long_df$variable_label == vl], na.rm = TRUE), numeric(1)), stringsAsFactors = FALSE)
      y_range_per_var <- vapply(present_order, function(vl) { r <- range(long_df$value[long_df$variable_label == vl], na.rm = TRUE); r[2] - r[1] }, numeric(1))
      ann_df$y <- ann_df$y + 0.08 * pmax(y_range_per_var, 1)
      pal_clin <- setNames(c("gray70", "#8B0000"), c("WT", mut_lab))
      p <- ggplot(long_df, aes(x = Mutation, y = value, fill = Mutation)) + geom_boxplot(alpha = 0.8) + geom_text(data = ann_df, aes(x = 1.5, y = y, label = sig), inherit.aes = FALSE, size = 5, fontface = "bold") + scale_fill_manual(values = pal_clin) + facet_wrap(~ variable_label, scales = "free_y", nrow = 1) + labs(x = NULL, y = NULL) + theme_minimal(base_size = 12) + theme(legend.position = "none", strip.text = element_text(size = 14, face = "bold"))
      ggsave(file, plot = p, width = 10, height = 4, dpi = 300)
    }
  )

  output$download_gene_summary_lollipop_plot <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_lollipop_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      g <- input$gene_summary
      if (is.null(g) || g == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void(), width = 6, height = 5, dpi = 150); return(invisible(NULL)) }
      df <- gene_summary_data()
      if (is.null(df) || nrow(df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No mutations") + theme_void(), width = 6, height = 5, dpi = 150); return(invisible(NULL)) }
      aa_col <- if ("HGVSp_Short" %in% colnames(df)) "HGVSp_Short" else if ("AA_change" %in% colnames(df)) "AA_change" else if ("Protein_Change" %in% colnames(df)) "Protein_Change" else if ("protein" %in% colnames(df)) "protein" else NULL
      if (is.null(aa_col)) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No protein-change column") + theme_void(), width = 6, height = 5, dpi = 150); return(invisible(NULL)) }
      g_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", g); if (!grepl("^[A-Z0-9]+$", g_base)) g_base <- g
      norm <- normalize_protein_change_for_maf(df[[aa_col]])
      has_position <- !is.na(norm$position) & norm$hgvsp != "p.?"
      if (sum(has_position) < 1L) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No valid AA positions") + theme_void(), width = 6, height = 5, dpi = 150); return(invisible(NULL)) }
      pos <- norm$position[has_position]; pos_range <- range(pos, na.rm = TRUE); if (diff(pos_range) < 1) pos_range <- pos_range + c(-5, 5)
      pos_counts <- as.data.frame(table(Position = pos), stringsAsFactors = FALSE); pos_counts$Position <- as.integer(pos_counts$Position); pos_counts$Freq <- as.numeric(pos_counts$Freq)
      if (nrow(pos_counts) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No positions") + theme_void(), width = 6, height = 5, dpi = 150); return(invisible(NULL)) }
      p <- ggplot(pos_counts, aes(x = Position, y = Freq)) + geom_segment(aes(x = Position, xend = Position, y = 0, yend = Freq), linewidth = 0.8, color = "#374e55") + geom_point(size = 3, color = "#8B0000", fill = "#8B0000", shape = 21) + labs(title = paste("Mutation distribution:", g_base), x = "Protein position", y = "Count") + theme_minimal(base_size = 12) + scale_x_continuous(limits = pos_range + c(-1, 1) * 0.02 * diff(pos_range))
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )

  output$download_gene_summary_survival <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_survival_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      g <- input$gene_summary
      if (is.null(g) || g == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      surv_uniq <- gene_summary_survival_km_data()
      if (is.null(surv_uniq) || length(unique(surv_uniq$Mutation)) < 2) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient groups") + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Mutation, data = surv_uniq)
      if (has_survminer) { p <- survminer::ggsurvplot(fit, data = surv_uniq, risk.table = TRUE, pval = TRUE, title = paste("Survival by", g), xlab = "Years", palette = c("#8B0000", "#4D4D4D"), pval.coord = c(0, 0.05)); png(file, width = 7, height = 6, units = "in", res = 300); print(p); dev.off() } else { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "survminer required") + theme_void(), width = 7, height = 5, dpi = 150) }
    }
  )

  output$download_gene_summary_oncoprint <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_oncoprint_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      od <- gene_summary_oncoprint_data()
      if (is.null(od) || !is.matrix(od$temp_dat) || nrow(od$temp_dat) == 0 || ncol(od$temp_dat) == 0) { png(file, width = 6, height = 4, res = 150); plot.new(); text(0.5, 0.5, "No data"); dev.off(); return(invisible(NULL)) }
      temp_dat <- od$temp_dat; samples <- od$samples; genes <- od$genes
      rows <- list()
      for (i in seq_len(nrow(temp_dat))) {
        for (j in seq_len(ncol(temp_dat))) {
          v <- temp_dat[i, j]
          if (is.na(v) || v == "") next
          vts <- strsplit(as.character(v), ";")[[1]]
          rows[[length(rows) + 1L]] <- data.frame(Sample = colnames(temp_dat)[j], Gene = rownames(temp_dat)[i], variant_type = vts[1], stringsAsFactors = FALSE)
        }
      }
      mut_long <- if (length(rows) > 0) do.call(rbind, rows) else NULL
      if (is.null(mut_long) || nrow(mut_long) == 0) { png(file, width = 6, height = 4, res = 150); plot.new(); text(0.5, 0.5, "No data"); dev.off(); return(invisible(NULL)) }
      PAL_VAR <- c(Deletion = "#374E55", Indel = "#DF8F44", ITD = "#79AF97", SNV = "#B24745", Splicing = "#6A6599", PTD = "tan", Other = "#80796B", Unknown = "#80796B")
      mut_long$Gene <- factor(mut_long$Gene, levels = rev(genes)); mut_long$Sample <- factor(mut_long$Sample, levels = samples)
      mut_long$Sample_num <- as.numeric(mut_long$Sample); mut_long$Gene_num <- as.numeric(mut_long$Gene)
      bg <- expand.grid(Sample = factor(samples, levels = samples), Gene = factor(genes, levels = rev(genes)))
      bg$Sample_num <- as.numeric(bg$Sample); bg$Gene_num <- as.numeric(bg$Gene)
      p <- ggplot() + geom_tile(data = bg, aes(x = Sample_num, y = Gene_num), fill = "gray97", width = 1, height = 1) + geom_rect(data = mut_long, aes(xmin = Sample_num - 0.25, xmax = Sample_num + 0.25, ymin = Gene_num - 0.45, ymax = Gene_num + 0.45, fill = variant_type)) + scale_fill_manual(values = PAL_VAR, name = "Variant") + scale_x_continuous(limits = c(0.5, length(samples) + 0.5), expand = c(0, 0)) + scale_y_continuous(breaks = seq_along(genes), labels = rev(genes), expand = c(0, 0)) + labs(x = NULL, y = NULL) + theme_minimal(base_size = 10) + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "right")
      ggsave(file, plot = p, width = max(10, ncol(temp_dat) * 0.12), height = max(6, nrow(temp_dat) * 0.25), dpi = 200)
    }
  )

  output$download_gene_summary_comut_heatmap <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_comut_heatmap_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      g <- input$gene_summary
      if (is.null(g) || g == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void(), width = 5, height = 5, dpi = 150); return(invisible(NULL)) }
      cooc <- cooccurrence_data()
      if (is.null(cooc$matrix) || !g %in% rownames(cooc$matrix)) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Gene not in matrix") + theme_void(), width = 5, height = 5, dpi = 150); return(invisible(NULL)) }
      mat <- cooc$matrix; others <- setdiff(colnames(mat), g)
      keep_genes <- filtered_genes(); if (length(keep_genes) > 0) others <- intersect(others, keep_genes)
      if (length(others) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No other genes") + theme_void(), width = 5, height = 5, dpi = 150); return(invisible(NULL)) }
      vec <- mat[others, g]; vec_log <- log2(vec + 0.01); vec_log[vec_log > 2] <- 2; vec_log[vec_log < -2] <- -2
      ord <- order(vec_log, decreasing = TRUE); genes_ord <- names(vec_log)[ord]
      df_plot <- data.frame(Gene = factor(genes_ord, levels = rev(genes_ord)), log2OR = vec_log[ord])
      p <- ggplot(df_plot, aes(x = 1, y = Gene, fill = log2OR)) + geom_tile() + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) + labs(x = NULL, y = NULL, fill = "log2(OR)", title = paste("Co-occurrence with", g)) + theme_minimal(base_size = 14) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      ggsave(file, plot = p, width = max(5, length(others) * 0.15), height = 6, dpi = 300)
    }
  )

  output$download_gene_summary_comut2_plot <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_comut2_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      df <- gene_summary_comut2_data()
      if (is.null(df) || length(unique(df$Group)) < 2) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select second gene") + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = df)
      if (has_survminer) { p <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE, title = "Pairwise co-mutation survival", xlab = "Years", pval.coord = c(0, 0.05)); png(file, width = 7, height = 6, units = "in", res = 300); print(p); dev.off() } else { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "survminer required") + theme_void(), width = 7, height = 5, dpi = 150) }
    }
  )

  output$download_gene_summary_comut3_plot <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_comut3_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      df <- gene_summary_comut3_data()
      if (is.null(df) || length(unique(df$Group)) < 2) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select second and third genes") + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = df)
      if (has_survminer) { p <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE, title = "Triple co-mutation survival", xlab = "Years", pval.coord = c(0, 0.05)); png(file, width = 7, height = 6, units = "in", res = 300); print(p); dev.off() } else { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "survminer required") + theme_void(), width = 7, height = 5, dpi = 150) }
    }
  )

  output$download_gene_summary_vaf_plot <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_", tolower(vaf_col()), "_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      g <- input$gene_summary; vcol <- vaf_col(); metric_short <- if (vcol == "CCF") "CCF" else "VAF"
      if (is.null(g) || g == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void(), width = 5, height = 4, dpi = 150); return(invisible(NULL)) }
      df <- filtered_data(); if (!vcol %in% colnames(df)) vcol <- "VAF"
      df <- df[as.character(df$Gene_for_analysis) == g & !is.na(df[[vcol]]) & df[[vcol]] > 0, , drop = FALSE]
      if (nrow(df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No", metric_short, "data")) + theme_void(), width = 5, height = 4, dpi = 150); return(invisible(NULL)) }
      df$.metric <- df[[vcol]]
      p <- ggplot(df, aes(x = .metric)) + geom_histogram(bins = 25, fill = "#374e55", alpha = 0.8) + labs(x = paste0(metric_short, " (%)"), y = "Count", title = paste(metric_short, "distribution for", g)) + theme_minimal(base_size = 14)
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300)
    }
  )

  output$download_gene_summary_vaf_survival_plot <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_", tolower(vaf_col()), "_survival_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      g <- input$gene_summary; vcol <- vaf_col(); metric_short <- if (vcol == "CCF") "CCF" else "VAF"
      if (is.null(g) || g == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      df <- filtered_data(); if (!vcol %in% colnames(df)) vcol <- "VAF"
      df <- df[as.character(df$Gene_for_analysis) == g & !is.na(df[[vcol]]) & df[[vcol]] > 0, c("Sample", vcol), drop = FALSE]
      if (nrow(df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No", metric_short, "data")) + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      colnames(df)[colnames(df) == vcol] <- ".metric"
      vaf_agg <- aggregate(.metric ~ Sample, data = df, max)
      surv_df <- survival_data(); surv_df <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
      merge_df <- merge(vaf_agg, surv_df, by = "Sample", all = FALSE)
      if (nrow(merge_df) < 20 || !has_maxstat) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste0("Need \u226520 samples; maxstat required")) + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      mst <- tryCatch(maxstat::maxstat.test(Surv(Time_to_OS, as.numeric(Censor)) ~ .metric, data = merge_df, smethod = "LogRank", minprop = 0.2, maxprop = 0.8), error = function(e) NULL)
      if (is.null(mst)) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "MaxStat failed") + theme_void(), width = 7, height = 5, dpi = 150); return(invisible(NULL)) }
      high_label <- paste0("High ", metric_short); low_label <- paste0("Low ", metric_short)
      merge_df$Group <- factor(ifelse(merge_df$.metric > mst$estimate, high_label, low_label), levels = c(high_label, low_label))
      fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = merge_df)
      if (has_survminer) { legend_labs <- gsub("^Group=", "", names(fit$strata)); p <- survminer::ggsurvplot(fit, data = merge_df, risk.table = TRUE, pval = TRUE, title = paste0(g, ": High vs Low ", metric_short, " (", round(mst$estimate, 1), "%)"), xlab = "Years", palette = c("#8B0000", "#4D4D4D"), legend.labs = legend_labs, pval.coord = c(0, 0.05)); png(file, width = 7, height = 6, units = "in", res = 300); print(p); dev.off() } else { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "survminer required") + theme_void(), width = 7, height = 5, dpi = 150) }
    }
  )

  output$download_gene_summary_drug_volcano <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_drug_volcano_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      g <- input$gene_summary
      if (is.null(g) || g == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void(), width = 6, height = 5, dpi = 150); return(invisible(NULL)) }
      df <- gene_summary_drug_mut_wt()
      if (is.null(df) || nrow(df) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No drug data") + theme_void(), width = 6, height = 5, dpi = 150); return(invisible(NULL)) }
      df <- df[!is.na(df$p_value) & df$p_value > 0, , drop = FALSE]
      df$neglog10p <- -log10(df$p_value)
      df$Category <- "Not significant"; df$Category[df$p_value < 0.05 & df$delta_AUC_mut_wt < 0] <- "Sensitive"; df$Category[df$p_value < 0.05 & df$delta_AUC_mut_wt > 0] <- "Resistant"
      df_sig <- df[df$Category != "Not significant", , drop = FALSE]
      label_sensitive <- df_sig[df_sig$Category == "Sensitive", , drop = FALSE]
      label_resistant <- df_sig[df_sig$Category == "Resistant", , drop = FALSE]
      if (nrow(label_sensitive) > 5) label_sensitive <- label_sensitive[order(label_sensitive$p_value)[seq_len(5)], , drop = FALSE]
      if (nrow(label_resistant) > 5) label_resistant <- label_resistant[order(label_resistant$p_value)[seq_len(5)], , drop = FALSE]
      df_label <- rbind(label_sensitive, label_resistant)
      p <- ggplot(df, aes(x = delta_AUC_mut_wt, y = neglog10p, color = Category, size = n_mut)) + geom_point(alpha = 0.8) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") + geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") + scale_color_manual(values = c("Sensitive" = "#b2182b", "Resistant" = "#2166ac", "Not significant" = "gray50"), name = NULL) + scale_size_continuous(name = "N (mut)", range = c(2, 10)) + labs(x = "Mean AUC (mut) - Mean AUC (WT)", y = "-log10(p)", title = paste(g, ": Mut vs WT drug AUC")) + theme_minimal(base_size = 12) + theme(legend.position = "top")
      if (nrow(df_label) > 0) {
        if (has_ggrepel) {
          p <- p + ggrepel::geom_label_repel(data = df_label, aes(x = delta_AUC_mut_wt, y = neglog10p, label = Inhibitor), size = 5, show.legend = FALSE, max.overlaps = Inf)
        } else {
          p <- p + geom_text(data = df_label, aes(x = delta_AUC_mut_wt, y = neglog10p, label = Inhibitor), hjust = -0.1, size = 5, show.legend = FALSE)
        }
      }
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )

  output$download_gene_summary_drug_dotplot <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_drug_dotplot_", format(Sys.time(), "%Y%m%d_%H%M"), ".png"),
    content = function(file) {
      g <- input$gene_summary
      if (is.null(g) || g == "") { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void(), width = 10, height = 6, dpi = 150); return(invisible(NULL)) }
      corr <- gene_summary_drug_correlations_gene()
      if (is.null(corr) || nrow(corr) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No drug data") + theme_void(), width = 10, height = 6, dpi = 150); return(invisible(NULL)) }
      if ("n" %in% colnames(corr)) corr <- corr[!is.na(corr$n) & corr$n >= 10, , drop = FALSE]
      base_genes <- gene_to_beataml_bases(g); sub <- corr[corr$Gene %in% base_genes, , drop = FALSE]
      if (nrow(sub) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No correlations") + theme_void(), width = 10, height = 6, dpi = 150); return(invisible(NULL)) }
      sub <- sub[!is.na(sub$p_value) & sub$p_value < 0.05, , drop = FALSE]
      if (nrow(sub) == 0) sub <- corr[corr$Gene %in% base_genes, , drop = FALSE]
      if (nrow(sub) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No p < 0.05 correlations") + theme_void(), width = 10, height = 6, dpi = 150); return(invisible(NULL)) }
      sub <- sub[order(sub$p_value), , drop = FALSE]; inhibitors <- unique(sub$Inhibitor)
      if (length(inhibitors) > 24) inhibitors <- head(inhibitors, 24)
      b <- beataml_for_drug(); if (is.null(b) || !b$ok) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No Beat AML data") + theme_void(), width = 10, height = 6, dpi = 150); return(invisible(NULL)) }
      subset <- input$gene_summary_drug_subset; if (is.null(subset) || subset == "") subset <- "de_novo"
      allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
      scatter_list <- list()
      for (inh in inhibitors) {
        auc_sub <- b$auc[b$auc$inhibitor == inh & b$auc$Sample %in% allowed, c("Sample", "auc"), drop = FALSE]
        mut_sub <- b$mutations[b$mutations$Gene %in% base_genes & b$mutations$Sample %in% allowed, c("Sample", "VAF"), drop = FALSE]
        merged <- merge(auc_sub, mut_sub, by = "Sample")
        if (nrow(merged) >= 5) {
          corr_row <- sub[sub$Inhibitor == inh, , drop = FALSE]
          if (nrow(corr_row) > 0) { fit <- tryCatch(lm(auc ~ VAF, data = merged), error = function(e) NULL); slope <- if (!is.null(fit)) coef(fit)[2] else NA_real_; trend <- if (is.na(slope)) "Unknown" else if (slope < 0) "Sensitive" else "Resistant"; merged$Inhibitor <- inh; merged$R_squared <- corr_row$R_squared[1]; merged$p_value <- corr_row$p_value[1]; merged$q_value <- corr_row$q_value[1]; merged$Trend <- trend; scatter_list[[length(scatter_list) + 1L]] <- merged }
        }
      }
      if (length(scatter_list) == 0) { ggsave(file, plot = ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient scatter data") + theme_void(), width = 10, height = 6, dpi = 150); return(invisible(NULL)) }
      scatter_df <- do.call(rbind, scatter_list); scatter_df$Inhibitor <- factor(scatter_df$Inhibitor, levels = inhibitors); scatter_df$Trend <- factor(scatter_df$Trend, levels = c("Sensitive", "Resistant", "Unknown"))
      label_df <- unique(scatter_df[, c("Inhibitor", "R_squared", "p_value", "q_value", "Trend"), drop = FALSE]); label_df$star <- ifelse(label_df$q_value < 0.1, "*", ""); label_df$label <- paste0("R² = ", round(label_df$R_squared, 3), ", p ", format_pval_display(label_df$p_value), label_df$star)
      trend_colors <- c("Sensitive" = "#b2182b", "Resistant" = "#2166ac", "Unknown" = "gray50")
      p <- ggplot(scatter_df, aes(x = VAF, y = auc, color = Trend)) + geom_point(size = 2.5, alpha = 0.7) + geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.8) + scale_color_manual(values = trend_colors, name = "Trend", guide = "none") + facet_wrap(~ Inhibitor, scales = "free", ncol = min(3, length(inhibitors))) + labs(x = "VAF (%)", y = "Drug AUC", title = paste("VAF vs AUC:", g, "(p < 0.05)")) + theme_minimal(base_size = 16) + theme(axis.text = element_text(size = 14), strip.text = element_text(size = 14, face = "bold"), strip.background = element_rect(fill = "#f0f0f0", color = NA), panel.spacing = unit(1, "lines"), plot.title = element_text(size = 17, face = "bold")) + geom_text(aes(x = Inf, y = Inf, label = label), hjust = 1.1, vjust = 1.5, size = 4, inherit.aes = FALSE, data = label_df)
      ggsave(file, plot = p, width = 4 * min(3, length(inhibitors)), height = 4 * ceiling(length(inhibitors) / 3), dpi = 300)
    }
  )

  output$download_gene_summary_drug_mut_wt_table <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_drug_mut_wt_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
    content = function(file) {
      mut_wt <- gene_summary_drug_mut_wt()
      if (is.null(mut_wt) || nrow(mut_wt) == 0) return(write("No data", file))
      write.csv(mut_wt[, c("Inhibitor", "delta_AUC_mut_wt", "p_value", "n_mut", "n_wt")], file, row.names = FALSE)
    }
  )

  output$download_gene_summary_drug_table <- downloadHandler(
    filename = function() paste0("gene_", input$gene_summary, "_drug_correlations_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"),
    content = function(file) {
      corr <- gene_summary_drug_correlations_gene()
      if (is.null(corr) || nrow(corr) == 0) return(write("No data", file))
      g <- input$gene_summary; base_genes <- gene_to_beataml_bases(g); sub <- corr[corr$Gene %in% base_genes, , drop = FALSE]
      if (nrow(sub) == 0) return(write("No data", file))
      cols <- c("Inhibitor", "Gene", "Direction", "R_squared", "p_value", "q_value", "LOOCV_RMSE", "n"); cols <- cols[cols %in% colnames(sub)]
      write.csv(sub[, cols, drop = FALSE], file, row.names = FALSE)
    }
  )

  output$data_table <- renderDT({
    df <- req(filtered_data())
    cols <- c("Sample", "Gene_for_analysis", "Gene", "VAF", "CCF", "CN_at_locus", "variant_type", "AA_change", "Gene_Group", "Subset", "Cohort", "Age", "Sex", "Risk", "Time_to_OS", "Censor", "mutation_category")
    cols <- cols[cols %in% colnames(df)]
    df <- df[, cols, drop = FALSE]
    for (nm in c("VAF", "CCF", "Age", "Time_to_OS", "Censor")) {
      if (nm %in% colnames(df) && is.numeric(df[[nm]])) df[[nm]] <- round(df[[nm]], 2)
    }
    datatable(df, filter = "top", options = list(pageLength = 10000, lengthChange = FALSE, scrollX = TRUE))
  })
}

shinyApp(ui = ui, server = server)
