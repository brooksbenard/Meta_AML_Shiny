# =============================================================================
# Meta AML Explorer - Shiny App for Exploratory Analysis of AML Mutational Data
# Based on Benard et al. Nature Communications 2021
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

# Load data (app expects final_data_matrix.RData in the same directory as app.R)
data_path <- file.path(getwd(), "final_data_matrix.RData")
if (!file.exists(data_path)) {
  stop("Data file not found. Place final_data_matrix.RData in the app directory and run the app from that directory.\n  Looked for: ", normalizePath(data_path, mustWork = FALSE))
}
load(data_path)

# Convert to plain data.frame - CSV round-trip strips data.table/list columns
tmp <- tempfile(fileext = ".csv")
utils::write.csv(final_data_matrix, tmp, row.names = FALSE)
final_data_matrix <- utils::read.csv(tmp, stringsAsFactors = FALSE)
unlink(tmp)
# Assign RTK-RAS Signaling to all FLT3 mutations if they lack mutation_category; homogenize RTK/RAS -> RTK-RAS
if ("mutation_category" %in% colnames(final_data_matrix) && "Gene" %in% colnames(final_data_matrix)) {
  final_data_matrix$mutation_category[final_data_matrix$mutation_category %in% c("RTK/RAS Signaling", "RTK-RAS Signaling")] <- "RTK-RAS Signaling"
  flt3_idx <- which(final_data_matrix$Gene %in% c("FLT3", "FLT3-ITD", "FLT3-TKD", "FLT3-Other") &
                    (is.na(final_data_matrix$mutation_category) | final_data_matrix$mutation_category == ""))
  if (length(flt3_idx) > 0) {
    final_data_matrix$mutation_category[flt3_idx] <- "RTK-RAS Signaling"
  }
}
# Normalize Cohort and Subset for display across the app (Beat_AML -> Beat AML; de_novo -> De novo, etc.)
if ("Cohort" %in% colnames(final_data_matrix)) {
  final_data_matrix$Cohort[final_data_matrix$Cohort == "Beat_AML"] <- "Beat AML"
}
if ("Subset" %in% colnames(final_data_matrix)) {
  final_data_matrix$Subset[is.na(final_data_matrix$Subset) | final_data_matrix$Subset == ""] <- "Other"
  final_data_matrix$Subset[final_data_matrix$Subset %in% c("de_novo", "De_novo")] <- "De novo"
  final_data_matrix$Subset[final_data_matrix$Subset == "secondary"] <- "Secondary"
  final_data_matrix$Subset[final_data_matrix$Subset == "relapse"] <- "Relapse"
  final_data_matrix$Subset[final_data_matrix$Subset %in% c("therapy", "tAML")] <- "Therapy"
  final_data_matrix$Subset[final_data_matrix$Subset %in% c("other", "oAML", "unknown", "")] <- "Other"
}

# Normalize AML_Meta_Cohort to match final_data_matrix column names so existing analyses work
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
  want <- c("Sample", "Gene", "VAF", "Cohort", "Subset", "Time_to_OS", "Censor", "variant_type",
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

# Load AML_Meta_Cohort_v2 (four cohorts merged) for Meta AML4 tab; fallback to v1 or final_data_matrix if missing
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
  # Classify genes in Meta AML4 using same mutation_category as Benard et al. (2021) (final_data_matrix)
  ref <- final_data_matrix[!is.na(final_data_matrix$mutation_category) & as.character(final_data_matrix$mutation_category) != "", c("Gene", "mutation_category"), drop = FALSE]
  if (nrow(ref) > 0) {
    ref$Gene <- as.character(ref$Gene)
    ref$mutation_category <- as.character(ref$mutation_category)
    most_common <- aggregate(mutation_category ~ Gene, data = ref, FUN = function(x) names(sort(table(x), decreasing = TRUE))[1])
    gene_to_cat <- setNames(as.character(most_common$mutation_category), as.character(most_common$Gene))
    AML_Meta_Cohort$mutation_category <- unname(gene_to_cat[as.character(AML_Meta_Cohort$Gene)])
    # Homogenize RTK/RAS -> RTK-RAS Signaling and assign to all FLT3 mutations
    AML_Meta_Cohort$mutation_category[AML_Meta_Cohort$mutation_category %in% c("RTK/RAS Signaling", "RTK-RAS Signaling")] <- "RTK-RAS Signaling"
    flt3_idx <- which(AML_Meta_Cohort$Gene %in% c("FLT3", "FLT3-ITD", "FLT3-TKD", "FLT3-Other"))
    if (length(flt3_idx) > 0) {
      AML_Meta_Cohort$mutation_category[flt3_idx] <- "RTK-RAS Signaling"
    }
    mll_idx <- which(AML_Meta_Cohort$Gene == "MLL")
    if (length(mll_idx) > 0) {
      AML_Meta_Cohort$mutation_category[mll_idx] <- "Chromatin/Cohesin"
    }
  }
} else {
  AML_Meta_Cohort <- final_data_matrix
}
tmp2 <- tempfile(fileext = ".csv")
utils::write.csv(AML_Meta_Cohort, tmp2, row.names = FALSE)
AML_Meta_Cohort <- utils::read.csv(tmp2, stringsAsFactors = FALSE)
unlink(tmp2)

# Paper color palettes (Benard et al. Nat Commun 2021)
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
  "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF"
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
  div(id = "data-loading-overlay", style = "display: none; position: fixed; top: 0; left: 0; right: 0; bottom: 0; background: rgba(255,255,255,0.95); z-index: 9999; align-items: center; justify-content: center;",
    tags$div(style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); text-align: center; padding: 2em;",
      div(class = "loading-spinner"),
      tags$div(style = "font-size: 18px; color: #374e55; margin-bottom: 8px; font-weight: 600;", "Loading data..."),
      tags$p(style = "color: #666; font-size: 14px; margin: 0;", "Preparing analyses. This may take a moment.")
    )
  ),
  tags$script(HTML("
    (function() {
      var overlayShownAt = 0;
      var minOverlayMs = 2500;
      var currentTab = null;
      var tabLoaded = {analyses: false, meta_aml4: false};
      var hideTimeout = null;
      
      function showOverlay() {
        overlayShownAt = Date.now();
        if (hideTimeout) {
          clearTimeout(hideTimeout);
          hideTimeout = null;
        }
        $('#data-loading-overlay').css({'display': 'flex'});
      }
      
      function hideOverlay() {
        var elapsed = Date.now() - overlayShownAt;
        var wait = Math.max(0, minOverlayMs - elapsed);
        hideTimeout = setTimeout(function() { 
          $('#data-loading-overlay').hide(); 
          hideTimeout = null;
        }, wait);
      }
      
      // Show overlay immediately when clicking analysis tabs (only if not already loaded)
      $(document).on('click', 'a[data-value=\"analyses\"]', function() { 
        currentTab = 'analyses';
        if (!tabLoaded.analyses) {
          showOverlay();
        }
      });
      $(document).on('click', 'a[data-value=\"meta_aml4\"]', function() { 
        currentTab = 'meta_aml4';
        if (!tabLoaded.meta_aml4) {
          showOverlay();
        }
      });
      
      // Track which key outputs have rendered
      var renderedOutputs = {};
      
      // Hide overlay only when data is ready AND key outputs have rendered
      Shiny.addCustomMessageHandler('dataReady', function(tab) { 
        if (tab === currentTab) {
          tabLoaded[tab] = true;
          // Don't hide yet - wait for outputs to render
        }
      });
      
      // Track tab changes
      $(document).on('shiny:inputchanged', function(event) {
        if (event.name === 'main_nav') {
          currentTab = event.value;
          renderedOutputs = {}; // Reset when switching tabs
          if ((currentTab === 'analyses' || currentTab === 'meta_aml4') && !tabLoaded[currentTab]) {
            showOverlay();
          }
        }
      });
      
      // Check when outputs finish rendering - hide overlay only after key outputs are ready
      $(document).on('shiny:value', function(event) {
        if (currentTab && (currentTab === 'analyses' || currentTab === 'meta_aml4')) {
          // Track key outputs that indicate the tab is fully loaded
          var keyOutputs = ['summary_table', 'cohort_plot', 'oncoprint_plot'];
          if (keyOutputs.indexOf(event.name) !== -1) {
            renderedOutputs[event.name] = true;
            
            // Check if all key outputs have rendered and data is ready
            var allRendered = keyOutputs.every(function(output) {
              return renderedOutputs[output] === true;
            });
            
            if (allRendered && tabLoaded[currentTab]) {
              // All outputs rendered and data ready - wait a bit more then hide
              var tabToCheck = currentTab;
              setTimeout(function() {
                if (currentTab === tabToCheck) {
                  hideOverlay();
                }
              }, 500);
            }
          }
        }
      });
    })();
  ")),
  div(class = "app-banner",
    div(class = "banner-title", "Meta AML Explorer"),
    div(class = "banner-contact",
      tags$a(href = "mailto:bbenard@stanford.edu", "bbenard@stanford.edu"),
      tags$img(src = "linkedin_pic%20copy.jpeg", alt = "Contact", title = "Contact")
    )
  ),
  tabsetPanel(
    id = "main_nav",
    tabPanel("About", value = "about",
      div(class = "welcome-page", style = "max-width: 800px; margin: 0 auto; padding: 24px 15px;",
        h2("Welcome to Meta AML Explorer"),
        p("This site provides interactive exploration of acute myeloid leukemia (AML) mutational and clinical outcomes data. The ", strong("Meta AML4"), " and ", strong("Benard et al. (2021)"), " tabs offer slightly different aggregated data sources but the same core analyses:"),
        tags$ul(
          tags$li("Cohort selection by AML type (e.g. de novo, secondary), dataset (e.g. UK-NCRI), karyotype (NK/Complex/Other), and minimum mutation frequency."),
          tags$li("Single mutation associations (clinical variables, survival, hazard ratios)"),
          tags$li("Co-mutation (heatmap, odds-ratios, Kaplan-Meier by 2–3 genes)"),
          tags$li("Variant allele frequency (VAF) associations (clinical variables, survival, pairwise scatterplots)"),
          tags$li("Drug Sensitivity (Mutation and VAF associations for BeatAML2 data)")
        ),
        h2("Meta AML4 tab"),
        p("Meta AML4 merges ", strong("four of the largest molecularly profiled and clinically annotated AML datasets"), " into a single combined cohort of ", strong("~4,660 patients"), ". Meta AML4 also includes ", strong("refined mutation assignments"), " (e.g., IDH2 R140 vs R172, NRAS G12/13 vs Q61/62, CEBPA bi-allelic vs mono-allelic) for more granular associations."),
        div(style = "background: #f8f9fa; border-left: 4px solid #0066cc; padding: 14px 18px; margin: 12px 0; border-radius: 0 6px 6px 0;",
          p(strong("Datasets:"), style = "margin-top: 0; margin-bottom: 8px;"),
          p(style = "margin: 4px 0; line-height: 1.6;", strong("UK-NCRI"), " (2,113 patients) | ", tags$a("Tazi et al. (2022) (Nature Communications)", href = "https://www.nature.com/articles/s41467-022-32103-8", target = "_blank"), " | ", tags$a("Data", href = "https://github.com/papaemmelab/Tazi_NatureC_AML", target = "_blank")),
          p(style = "margin: 4px 0; line-height: 1.6;", strong("AML-SG"), " (1,540 patients) | ", tags$a("Papaemmanuil et al. (2016), NEJM", href = "https://www.nejm.org/doi/10.1056/NEJMoa1516192", target = "_blank"), " | ", tags$a("Data", href = "https://github.com/gerstung-lab/AML-multistage", target = "_blank")),
          p(style = "margin: 4px 0; line-height: 1.6;", strong("Beat AML"), " (805 patients) | ", tags$a("Tyner et al. (2022), Cancer Cell", href = "https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00312-9", target = "_blank"), " | ", tags$a("Data", href = "https://biodev.github.io/BeatAML2/", target = "_blank")),
          p(style = "margin: 4px 0; line-height: 1.6; margin-bottom: 0;", strong("TCGA LAML"), " (200 patients) | ", tags$a("NEJM", href = "https://www.nejm.org/doi/full/10.1056/NEJMoa1301689", target = "_blank"), " | ", tags$a("Data", href = "https://www.cbioportal.org/", target = "_blank"))
        ),
        h2("Benard et al. (2021) tab"),
        p("Data from ", strong("Benard et al. (2021)"), " ", em("Clonal architecture predicts clinical outcomes and drug sensitivity in acute myeloid leukemia."), " Nature Communications | ", tags$a("Paper", href = "https://www.nature.com/articles/s41467-021-27472-5", target = "_blank")),
        p("The paper aggregates mutation calls, variant allele frequencies (VAF), and clinical outcomes from ~2,800 patients across 12 cohorts (TCGA, Beat AML/Tyner, Papaemmanuil, Majeti, and others) to link clonal architecture to prognosis and therapy response."),
        h2("Caveats"),
        p("This site is intended for research purposes only and is in active development. Some initial data loadings may take a few seconds."),
        p("The current implementation does not account for copy number correction of variant allele frequencies (VAFs). A future implementation will incorporate copy number data to define and use cancer cell fractions of mutations instead of VAFs."),
        p("Available cytogenetic data will be incorporated in a future version."),
        h2("About the author"),
        p("Brooks Benard, Stanford University. ", tags$a("Website", href = "https://brooksbenard.github.io/", target = "_blank"), " · ", tags$a("GitHub", href = "https://github.com/brooksbenard", target = "_blank"), " · ", tags$a(href = "mailto:bbenard@stanford.edu", "bbenard@stanford.edu")),
        hr(),
        p(style = "font-size: 0.9em; color: #666;",
          "© ", format(Sys.Date(), "%Y"), " Brooks Benard. Licensed under ",
          tags$a("MIT License", href = "https://opensource.org/licenses/MIT", target = "_blank"), "."
        )
      )
    ),
    tabPanel("Meta AML4", value = "meta_aml4"),
    tabPanel("Benard et al. (2021)", value = "analyses"),
    conditionalPanel(
      condition = "input.main_nav === 'analyses' || input.main_nav === 'meta_aml4'",
      sidebarLayout(
        sidebarPanel(
          width = 2,
          h4("Filter analyses by the following:"),
          selectInput("subset", "AML Subset:", 
            choices = c("All", "De novo", "Secondary", "Relapse", "Therapy", "Other"), selected = "De novo"),
          selectInput("cohort", "Cohort:", 
            choices = local({
              coh <- c("All", sort(unique(as.character(na.omit(final_data_matrix$Cohort)))))
              coh[coh == "Beat_AML"] <- "Beat AML"
              coh
            }), selected = "All"),
          selectInput("karyotype", "Karyotype:", 
            choices = c("All", "Complex", "Normal", "Other", "Unknown"), selected = "All"),
          sliderInput("min_freq", "Min. mutation frequency:", min = 1, max = 100, value = 50, step = 1),
          conditionalPanel(
            condition = "input.main_nav === 'meta_aml4'",
            hr(),
            p(strong("Only interested in a specific gene?")),
            wellPanel(
              h5("Summarize all associations for a single gene:"),
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

  # Format p/q for display: < 0.001, < 0.01, < 0.05, or numeric
  format_pval_display <- function(p) {
    if (is.null(p) || length(p) == 0) return(p)
    out <- character(length(p))
    for (i in seq_along(p)) {
      if (is.na(p[i])) { out[i] <- NA_character_; next }
      if (p[i] < 0.001) out[i] <- "< 0.001"
      else if (p[i] < 0.01) out[i] <- "< 0.01"
      else if (p[i] < 0.05) out[i] <- "< 0.05"
      else out[i] <- as.character(round(p[i], 2))
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
      # Meta AML4: three main tabs — Overview, All Gene Associations (sub-tabs), Single Gene Associations
      df <- filtered_data()
      gene_choices <- if (is.null(df) || !"Gene_for_analysis" %in% colnames(df)) {
        c("Select gene..." = "")
      } else {
        genes <- sort(unique(as.character(df$Gene_for_analysis)))
        c("Select gene..." = "", setNames(genes, genes))
      }
      current_tab <- selected_meta4_tab()
      if (is.null(current_tab) || !current_tab %in% c("Overview", "All Gene Associations", "Single Gene Associations")) current_tab <- "Overview"
      inner_tab <- selected_meta4_allgene_tab()
      if (is.null(inner_tab) || !inner_tab %in% c("Single", "Co-mutation", "VAF", "Drug Sensitivity")) inner_tab <- "Single"
      tabsetPanel(
        id = "main_tabs_meta4",
        selected = current_tab,
        tabPanel("Overview",
          fluidRow(
            column(4, wellPanel(h4("Selected Cohort Summary"), tableOutput("summary_table"))),
            column(8, wellPanel(h4("Samples by Cohort & Subset"), plotOutput("cohort_plot", height = 280)))
          ),
          fluidRow(
            column(12, wellPanel(
              h4("OncoPrint: Mutations & Clinical Annotations"),
              plotOutput("oncoprint_plot", height = 500)
            ))
          )
        ),
        tabPanel("All Gene Associations",
          tabsetPanel(
            id = "meta4_all_gene_tabs",
            selected = inner_tab,
            tabPanel("Single",
              fluidRow(
                column(12, wellPanel(
                  h4("Clinical Variable by Mutation"),
                  selectInput("clin_var", "Variable:", choices = c("WBC", "Age", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent")),
                  plotOutput("clinical_plot", height = 400)))
              ),
              fluidRow(
                column(6, wellPanel(h4("Kaplan-Meier"), selectInput("surv_gene", "Gene:", choices = gene_choices), plotOutput("survival_plot", height = 500))),
                column(6, wellPanel(h4("Hazard Ratios"), plotOutput("forest_plot", height = 500)))
              ),
              fluidRow(column(12, wellPanel(h4("Survival Summary"), DTOutput("survival_table"))))
            ),
            tabPanel("Co-mutation",
              fluidRow(
                column(6, wellPanel(h4("Co-mutation Odds Ratio"), plotOutput("cooccurrence_plot", height = 500))),
                column(6, wellPanel(h4("Top Co-occurring Pairs"), div(style = "height: 500px; overflow-y: auto;", DTOutput("cooccurrence_table"))))
              ),
              fluidRow(
                column(4, wellPanel(
                  h4("Gene Selection"),
                  radioButtons("comut_n", "Number of genes:", choices = c("2 genes" = 2, "3 genes" = 3), selected = 2),
                  selectInput("comut_gene1", "Gene 1:", choices = gene_choices),
                  selectInput("comut_gene2", "Gene 2:", choices = gene_choices),
                  conditionalPanel(condition = "input.comut_n == 3",
                    selectInput("comut_gene3", "Gene 3:", choices = gene_choices))
                )),
                column(8, wellPanel(
                  h4("Kaplan-Meier by Co-mutation Status"),
                  uiOutput("comut_survival_plot_ui")
                ))
              )
            ),
            tabPanel("VAF",
              tabsetPanel(
                tabPanel("Single Gene",
                  fluidRow(
                    column(6, wellPanel(h4("VAF by Gene"), plotOutput("vaf_gene_plot", height = 450))),
                    column(6, wellPanel(
                      h4("VAF and Survival (MaxStat)"),
                      p(em("Optimal VAF threshold from maximally selected rank statistics (maxstat); log-rank. High vs Low VAF split among mutated patients.")),
                      selectInput("vaf_surv_gene", "Gene:", choices = gene_choices),
                      plotOutput("vaf_survival_plot", height = 450)))
                  ),
                  fluidRow(
                    column(6, wellPanel(h4("VAF by Mutation Category"), plotOutput("vaf_category_plot", height = 450))),
                    column(6, wellPanel(h4("VAF by Cohort"), plotOutput("vaf_cohort_plot", height = 450)))
                  )
                ),
                tabPanel("Pairwise Gene",
                  fluidRow(
                    column(12, wellPanel(h4("VAF Scatter (Clonal Ordering)"),
                      p("Compare VAF of two genes. Points above line = Gene 1 before Gene 2."),
                      selectInput("vaf_gene1", "Gene 1:", choices = gene_choices), selectInput("vaf_gene2", "Gene 2:", choices = gene_choices),
                      plotOutput("vaf_scatter_plot", height = 380)))
                  )
                )
              )
            ),
            tabPanel("Drug Sensitivity",
              p("VAF and mutation correlations with inhibitor sensitivity from BeatAML2 data.",
                a("biodev.github.io/BeatAML2", href = "https://biodev.github.io/BeatAML2/", target = "_blank")),
              fluidRow(
                column(12, wellPanel(
                  h4("Summary"),
                  div(id = "drug_subset_wrap", selectInput("drug_subset", "AML Subset:", choices = c("All", "de_novo", "secondary"), selected = "de_novo")),
                  uiOutput("drug_summary_ui")
                ))
              ),
              fluidRow(
                column(12, wellPanel(
                  h4("Mut vs WT AUC (≥10 mut samples)"),
                  uiOutput("drug_mut_wt_summary_ui"),
                  p(em("Delta AUC = mean AUC (mut) - mean AUC (WT). * q < 0.1 (FDR). Only gene-inhibitor pairs with p < 0.05 shown.")),
                  plotOutput("drug_mut_wt_heatmap", height = 500)
                ))
              ),
              fluidRow(
                column(12, wellPanel(
                  h4("VAF vs AUC correlations"),
                  uiOutput("drug_vaf_auc_summary_ui"),
                  plotOutput("drug_summary_dotplot", height = 500),
                  p(em("* q < 0.1 (FDR)"), style = "font-size: 11px; color: #666;")
                ))
              ),
              fluidRow(
                column(6, wellPanel(
                  h4("Individual Mutation and VAF Associations"),
                  p(em("Gene and inhibitor options match the subset used in Mut vs WT and VAF vs AUC analyses above.")),
                  uiOutput("drug_scatter_selects_ui"),
                  fluidRow(
                    column(6, plotOutput("drug_scatter_boxplot", height = 320)),
                    column(6, plotOutput("drug_scatter_plot", height = 320))
                  )
                )),
                column(6, wellPanel(
                  h4("Leave-One-Out Cross Validation for VAF vs AUC analysis"),
                  p(em("RMSE = leave-one-out prediction error in AUC units.")),
                  plotOutput("drug_loo_heatmap", height = 350)))
              ),
              fluidRow(
                column(12, wellPanel(h4("Top Correlations"), p(em("LOOCV RMSE = leave-one-out prediction error in AUC units (overfitting check)."), style = "font-size: 11px; color: #666;"), div(style = "height: 350px; overflow-y: auto;", DTOutput("drug_correlation_table"))))
              )
            )
          )
        ),
        tabPanel("Single Gene Associations",
          uiOutput("gene_summary_ui")
        )
      )
    } else if (identical(nav, "analyses")) {
      # Benard et al.: Overview, Single mutation associations, Co-mutation, VAF, Drug (Data Table hidden for now)
      df <- filtered_data()
      gene_choices_analyses <- if (is.null(df) || !"Gene_for_analysis" %in% colnames(df)) {
        c("Select..." = "")
      } else {
        genes <- sort(unique(as.character(df$Gene_for_analysis)))
        c("Select..." = "", setNames(genes, genes))
      }
      current_analyses_tab <- selected_analyses_tab()
      if (is.null(current_analyses_tab) || !current_analyses_tab %in% c("Overview", "Single mutation associations", "Co-mutation Associations", "VAF Associations", "Drug Sensitivity")) current_analyses_tab <- "Overview"
      tabsetPanel(
        id = "main_tabs_analyses",
        selected = current_analyses_tab,
        tabPanel("Overview",
          fluidRow(
            column(4, wellPanel(h4("Selected Cohort Summary"), tableOutput("summary_table"))),
            column(8, wellPanel(h4("Samples by Cohort & Subset"), plotOutput("cohort_plot", height = 280)))
          ),
          fluidRow(
            column(12, wellPanel(
              h4("OncoPrint: Mutations & Clinical Annotations"),
              plotOutput("oncoprint_plot", height = 500)
            ))
          )
        ),
        tabPanel("Single mutation associations",
          fluidRow(
            column(12, wellPanel(
              h4("Clinical Variable by Mutation"),
              selectInput("clin_var", "Variable:", choices = c("WBC", "Age", "Hemoglobin", "Platelet", "BM_blast_percent", "PB_blast_percent")),
              plotOutput("clinical_plot", height = 400)))
          ),
          fluidRow(
            column(6, wellPanel(h4("Kaplan-Meier"), selectInput("surv_gene", "Gene:", choices = gene_choices_analyses), plotOutput("survival_plot", height = 500))),
            column(6, wellPanel(h4("Hazard Ratios"), plotOutput("forest_plot", height = 500)))
          ),
          fluidRow(column(12, wellPanel(h4("Survival Summary"), DTOutput("survival_table"))))
        ),
        tabPanel("Co-mutation Associations",
          fluidRow(
            column(6, wellPanel(h4("Co-mutation Odds Ratio"), plotOutput("cooccurrence_plot", height = 500))),
            column(6, wellPanel(h4("Top Co-occurring Pairs"), div(style = "height: 500px; overflow-y: auto;", DTOutput("cooccurrence_table"))))
          ),
          fluidRow(
            column(4, wellPanel(
              h4("Gene Selection"),
              radioButtons("comut_n", "Number of genes:", choices = c("2 genes" = 2, "3 genes" = 3), selected = 2),
              selectInput("comut_gene1", "Gene 1:", choices = gene_choices_analyses),
              selectInput("comut_gene2", "Gene 2:", choices = gene_choices_analyses),
              conditionalPanel(condition = "input.comut_n == 3",
                selectInput("comut_gene3", "Gene 3:", choices = gene_choices_analyses))
            )),
            column(8, wellPanel(
              h4("Kaplan-Meier by Co-mutation Status"),
              uiOutput("comut_survival_plot_ui")
            ))
          )
        ),
        tabPanel("VAF Associations",
          tabsetPanel(
            tabPanel("Single Gene",
              fluidRow(
                column(6, wellPanel(h4("VAF by Gene"), plotOutput("vaf_gene_plot", height = 450))),
                column(6, wellPanel(
                  h4("VAF and Survival (MaxStat)"),
                  p(em("Optimal VAF threshold from maximally selected rank statistics (maxstat); log-rank. High vs Low VAF split among mutated patients.")),
                  selectInput("vaf_surv_gene", "Gene:", choices = gene_choices_analyses),
                  plotOutput("vaf_survival_plot", height = 450)))
              ),
              fluidRow(
                column(6, wellPanel(h4("VAF by Mutation Category"), plotOutput("vaf_category_plot", height = 450))),
                column(6, wellPanel(h4("VAF by Cohort"), plotOutput("vaf_cohort_plot", height = 450)))
              )
            ),
            tabPanel("Pairwise Gene",
              fluidRow(
                column(12, wellPanel(h4("VAF Scatter (Clonal Ordering)"),
                  p("Compare VAF of two genes. Points above line = Gene 1 before Gene 2."),
                  selectInput("vaf_gene1", "Gene 1:", choices = gene_choices_analyses), selectInput("vaf_gene2", "Gene 2:", choices = gene_choices_analyses),
                  plotOutput("vaf_scatter_plot", height = 380)))
              )
            )
          )
        ),
        tabPanel("Drug Sensitivity",
          p("VAF and mutation correlations with inhibitor sensitivity from BeatAML2 data.",
            a("biodev.github.io/BeatAML2", href = "https://biodev.github.io/BeatAML2/", target = "_blank")),
          fluidRow(
            column(12, wellPanel(
              h4("Summary"),
              div(id = "drug_subset_wrap", selectInput("drug_subset", "AML Subset:", choices = c("All", "de_novo", "secondary"), selected = "de_novo")),
              uiOutput("drug_summary_ui")
            ))
          ),
          fluidRow(
            column(12, wellPanel(
              h4("Mut vs WT AUC (≥10 mut samples)"),
              uiOutput("drug_mut_wt_summary_ui"),
              p(em("Delta AUC = mean AUC (mut) - mean AUC (WT). * q < 0.1 (FDR). Only gene-inhibitor pairs with p < 0.05 shown.")),
              plotOutput("drug_mut_wt_heatmap", height = 500)
            ))
          ),
          fluidRow(
            column(12, wellPanel(
              h4("VAF vs AUC correlations"),
              uiOutput("drug_vaf_auc_summary_ui"),
              plotOutput("drug_summary_dotplot", height = 500),
              p(em("* q < 0.1 (FDR)"), style = "font-size: 11px; color: #666;")
            ))
          ),
          fluidRow(
            column(6, wellPanel(
              h4("Individual Mutation and VAF Associations"),
              p(em("Gene and inhibitor options match the subset used in Mut vs WT and VAF vs AUC analyses above.")),
              uiOutput("drug_scatter_selects_ui"),
              fluidRow(
                column(6, plotOutput("drug_scatter_boxplot", height = 320)),
                column(6, plotOutput("drug_scatter_plot", height = 320))
              )
            )),
            column(6, wellPanel(
              h4("Leave-One-Out Cross Validation for VAF vs AUC analysis"),
              p(em("RMSE = leave-one-out prediction error in AUC units.")),
              plotOutput("drug_loo_heatmap", height = 350)))
          ),
          fluidRow(
            column(12, wellPanel(h4("Top Correlations"), p(em("LOOCV RMSE = leave-one-out prediction error in AUC units (overfitting check)."), style = "font-size: 11px; color: #666;"), div(style = "height: 350px; overflow-y: auto;", DTOutput("drug_correlation_table"))))
          )
        )
        # Data Table tab hidden for now; re-add tabPanel("Data Table", ...) here to restore
      )
    } else {
      # Default: empty div when neither tab is active
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
        strong("Meta AML4 data file not found."), " Add ", code("AML_Meta_Cohort_v2.rds"), " or ", code("AML_Meta_Cohort.rds"), " (or ", code("AML_Meta_Cohort.RData"), ") to the app directory and redeploy. This tab is currently showing the Benard et al. (2021) dataset."
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
  
  # Invalidate cache when filters change
  observeEvent(c(input$subset, input$cohort, input$karyotype, input$min_freq), {
    nav <- input$main_nav
    if (nav %in% c("analyses", "meta_aml4")) {
      cached_filtered_data[[nav]] <- NULL
      cached_data_params[[nav]] <- NULL
    }
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
    current_params <- list(subset = input$subset, cohort = eff_cohort, karyotype = input$karyotype, min_freq = input$min_freq)
    if (!is.null(cached_filtered_data[[nav]]) && 
        identical(cached_data_params[[nav]], current_params)) {
      return(cached_filtered_data[[nav]])
    }
    
    # Compute filtered data
    df <- prepare_data(current_data_src())
    if (!"Gene_for_analysis" %in% colnames(df)) df$Gene_for_analysis <- as.character(df$Gene)
    if (input$subset != "All") df <- df[as.character(df$Subset) == input$subset, , drop = FALSE]
    if (eff_cohort != "All") df <- df[as.character(df$Cohort) == eff_cohort, , drop = FALSE]
    if (!is.null(input$karyotype) && input$karyotype != "All" && "Karyotype" %in% colnames(df)) {
      df <- df[as.character(df$Karyotype) == input$karyotype, , drop = FALSE]
    }
    gene_counts <- table(as.character(df$Gene_for_analysis))
    keep_genes <- names(gene_counts)[gene_counts >= input$min_freq]
    df <- df[as.character(df$Gene_for_analysis) %in% keep_genes, , drop = FALSE]
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

  observe({
    req(df <- filtered_data())
    genes <- sort(unique(as.character(df$Gene_for_analysis)))
    updateSelectInput(session, "surv_gene", choices = c("Select..." = "", genes))
    updateSelectInput(session, "vaf_gene1", choices = c("Select..." = "", genes))
    updateSelectInput(session, "vaf_gene2", choices = c("Select..." = "", genes))
    updateSelectInput(session, "comut_gene1", choices = c("Select..." = "", genes))
    updateSelectInput(session, "comut_gene2", choices = c("Select..." = "", genes))
    updateSelectInput(session, "comut_gene3", choices = c("Select..." = "", genes))
    updateSelectInput(session, "vaf_surv_gene", choices = c("Select..." = "", genes))
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
      Value = round(vals, 2)
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
    keep_genes <- names(gene_tbl)[gene_tbl >= input$min_freq]
    final_data_matrix_3 <- final_data_matrix_3[as.character(final_data_matrix_3$Gene) %in% keep_genes, , drop = FALSE]
    if (nrow(final_data_matrix_3) == 0) return(NULL)
    # Cap at 30 genes for oncoprint regardless of min_freq
    sorted <- names(sort(gene_tbl[keep_genes], decreasing = TRUE))
    top_genes <- head(sorted, 30L)
    final_data_matrix_3 <- final_data_matrix_3[as.character(final_data_matrix_3$Gene) %in% top_genes, , drop = FALSE]
    samples <- unique(as.character(final_data_matrix_3$Sample))
    final_data_matrix_3 <- final_data_matrix_3[final_data_matrix_3$Sample %in% samples, , drop = FALSE]
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

  output$vaf_gene_plot <- renderPlot({
    req(df <- filtered_data())
    df <- df[!is.na(df$VAF) & df$VAF > 0, , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    gene_tbl <- table(as.character(df$Gene_for_analysis))
    top_genes <- names(sort(gene_tbl, decreasing = TRUE))[1:min(25, length(gene_tbl))]
    df <- df[as.character(df$Gene_for_analysis) %in% top_genes, , drop = FALSE]
    n_cat <- if ("mutation_category" %in% colnames(df)) length(unique(df$mutation_category[!is.na(df$mutation_category) & df$mutation_category != ""])) else 0L
    fill_var <- if (n_cat >= 2L) "mutation_category" else "Gene_for_analysis"
    ggplot(df, aes(x = VAF, y = reorder(Gene_for_analysis, VAF, median), fill = .data[[fill_var]])) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_x_continuous(limits = c(0, 100)) +
      (if (fill_var == "mutation_category") scale_fill_manual(values = PAL_MUT_CAT, na.value = "gray70", guide = "legend") else scale_fill_discrete(guide = "none")) +
      labs(x = "Variant Allele Frequency (%)", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.text = element_text(size = 13), legend.title = element_text(size = 14), legend.position = if (fill_var == "mutation_category") "right" else "none")
  })

  output$vaf_cohort_plot <- renderPlot({
    req(df <- filtered_data())
    df <- df[!is.na(df$VAF) & df$VAF > 0 & !is.na(df$Cohort), , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    cohort_pal <- if (is_meta4()) META4_COHORT_COL else PAL_COHORT
    # Filter palette to only include cohorts present in data
    present_cohorts <- unique(df$Cohort[!is.na(df$Cohort)])
    cohort_pal <- cohort_pal[names(cohort_pal) %in% present_cohorts]
    ggplot(df, aes(x = VAF, y = reorder(Cohort, VAF, median), fill = Cohort)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_x_continuous(limits = c(0, 100)) +
      scale_fill_manual(values = cohort_pal, na.value = "gray70") +
      labs(x = "Variant Allele Frequency (%)", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.position = "none")
  })

  output$vaf_category_plot <- renderPlot({
    req(df <- filtered_data())
    df <- df[!is.na(df$VAF) & df$VAF > 0 & !is.na(df$mutation_category) & df$mutation_category != "", , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    ggplot(df, aes(x = VAF, y = reorder(mutation_category, VAF, median), fill = mutation_category)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_x_continuous(limits = c(0, 100)) +
      scale_fill_manual(values = PAL_MUT_CAT, na.value = "gray70") +
      labs(x = "Variant Allele Frequency (%)", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.position = "none")
  })

  vaf_survival_data <- reactive({
    gene <- input$vaf_surv_gene
    if (is.null(gene) || gene == "") return(NULL)
    df <- filtered_data()
    if (!"VAF" %in% colnames(df)) return(NULL)
    df <- df[as.character(df$Gene_for_analysis) == gene & !is.na(df$VAF) & df$VAF > 0, c("Sample", "VAF"), drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    vaf_agg <- aggregate(VAF ~ Sample, data = df, max)
    surv_df <- survival_data()
    surv_df <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    merge_df <- merge(vaf_agg, surv_df, by = "Sample", all = FALSE)
    if (nrow(merge_df) < 20) return(NULL)
    merge_df
  })

  output$vaf_survival_plot <- renderPlot({
    gene <- input$vaf_surv_gene
    if (is.null(gene) || gene == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    vaf_surv_df <- vaf_survival_data()
    if (is.null(vaf_surv_df) || nrow(vaf_surv_df) < 20) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Need at least 20 mutated samples with VAF and survival data. Install the maxstat package for optimal threshold.") + theme_void())
    }
    if (!has_maxstat) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install the maxstat package for optimal VAF threshold (maximally selected rank statistics).") + theme_void())
    }
    mst <- tryCatch(maxstat::maxstat.test(Surv(Time_to_OS, as.numeric(Censor)) ~ VAF, data = vaf_surv_df, smethod = "LogRank", minprop = 0.2, maxprop = 0.8), error = function(e) NULL)
    if (is.null(mst)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "MaxStat failed for this gene.") + theme_void())
    cutpoint <- mst$estimate
    vaf_surv_df$Group <- factor(ifelse(vaf_surv_df$VAF > cutpoint, "High VAF", "Low VAF"), levels = c("High VAF", "Low VAF"))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = vaf_surv_df)
    sd <- survdiff(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = vaf_surv_df)
    pval <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
    pval_txt <- if (pval < 0.001) "p < 0.001" else paste0("p = ", format.pval(pval, digits = 2))
    if (has_survminer) {
      legend_labs <- gsub("^Group=", "", names(fit$strata))
      # Unnamed palette: first = High VAF (red), second = Low VAF (grey); Group factor ensures this order
      p <- survminer::ggsurvplot(fit, data = vaf_surv_df, risk.table = TRUE, pval = pval_txt,
        title = paste0(gene, ": High vs Low VAF (threshold ", round(cutpoint, 1), "%)"), xlab = "Years",
        palette = c("#8B0000", "#4D4D4D"), legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.1))
      print(p)
    } else {
      legend_labs <- gsub("^Group=", "", names(fit$strata))
      sf_dat <- data.frame(time = fit$time, surv = fit$surv, strata = rep(legend_labs, fit$strata))
      ggplot(sf_dat, aes(x = time, y = surv, color = strata)) +
        geom_step(linewidth = 1) +
        scale_color_manual(values = c("High VAF" = "#8B0000", "Low VAF" = "#4D4D4D")) +
        labs(title = paste0(gene, ": High vs Low VAF (", round(cutpoint, 1), "%)"), subtitle = pval_txt, x = "Years", y = "Survival probability") +
        theme_minimal() + theme(legend.position = "right")
    }
  })

  survival_data <- reactive({
    df <- filtered_data()
    df <- df[!is.na(df$Time_to_OS) & !is.na(df$Censor), , drop = FALSE]
    df$Time_to_OS <- as.numeric(df$Time_to_OS) / 365
    df
  })

  output$survival_plot <- renderPlot({
    gene <- input$surv_gene
    if (is.null(gene) || gene == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    surv_df <- survival_data()
    mut_samples <- unique(surv_df$Sample[surv_df$Gene_for_analysis == gene])
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    surv_uniq$Mutation <- ifelse(surv_uniq$Sample %in% mut_samples, paste0(gene, " mut"), "WT")
    if (length(unique(surv_uniq$Mutation)) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient groups"))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Mutation, data = surv_uniq)
    strata_names <- names(fit$strata)
    is_mut <- grepl(" mut", strata_names, fixed = TRUE)
    surv_pal <- ifelse(is_mut, "#8B0000", "#4D4D4D")  # deep red for mutated, gray for WT
    if (has_survminer) {
      legend_labs <- gsub("^Mutation=", "", names(fit$strata))
      p <- survminer::ggsurvplot(fit, data = surv_uniq, risk.table = TRUE, pval = TRUE,
        title = paste("Survival by", gene, "mutation status"), xlab = "Years",
        palette = surv_pal, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.1))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install the survminer package for Kaplan-Meier plots with number at risk.", size = 4) + theme_void()
    }
  })

  forest_data <- reactive({
    surv_df <- survival_data()
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    gene_tbl <- table(as.character(surv_df$Gene_for_analysis))
    genes <- names(gene_tbl)[gene_tbl >= input$min_freq]
    res <- list()
    for (g in genes) {
      mut_pts <- unique(surv_df$Sample[surv_df$Gene_for_analysis == g])
      surv_uniq$mut <- surv_uniq$Sample %in% mut_pts
      if (sum(surv_uniq$mut) >= 10 && sum(!surv_uniq$mut) >= 10) {
        m <- tryCatch(coxph(Surv(Time_to_OS, as.numeric(Censor)) ~ mut, data = surv_uniq), error = function(e) NULL)
        if (!is.null(m)) {
          ci <- confint(m)
          res[[g]] <- data.frame(Gene = g, HR = exp(coef(m))["mutTRUE"], lower = exp(ci[1]), upper = exp(ci[2]),
            p = summary(m)$coefficients["mutTRUE", "Pr(>|z|)"], stringsAsFactors = FALSE)
        }
      }
    }
    out <- do.call(rbind, res)
    if (!is.null(out)) out <- out[!is.na(out$HR) & out$lower < 10 & out$upper < 10, , drop = FALSE]
    else out <- data.frame(Gene = character(), HR = numeric(), lower = numeric(), upper = numeric(), p = numeric())
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
        palette = pal_vec, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.1))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install the survminer package for Kaplan-Meier plots with number at risk.", size = 4) + theme_void()
    }
  })

  output$forest_plot <- renderPlot({
    df <- forest_data()
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    df <- df[order(df$HR), , drop = FALSE]
    df$Gene <- factor(df$Gene, levels = df$Gene)
    ggplot(df, aes(x = Gene, y = HR, ymin = lower, ymax = upper, color = HR > 1)) +
      geom_pointrange(size = 0.5) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      coord_flip() +
      scale_color_manual(values = c("TRUE" = "#762a83", "FALSE" = "#1b7837"), guide = "none") +
      labs(x = NULL, y = "Hazard Ratio") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15))
  })

  output$survival_table <- renderDT({
    df <- forest_data()
    if (is.null(df) || nrow(df) == 0) return(datatable(data.frame(Gene = character(), HR_CI = character(), p = character(), p_adj = character()), options = list(pageLength = 15), rownames = FALSE))
    df$HR_CI <- sprintf("%.2f (%.2f-%.2f)", df$HR, df$lower, df$upper)
    df$p_adj <- format_pval_display(p.adjust(df$p, method = "fdr"))
    df$p <- format_pval_display(df$p)
    datatable(df[, c("Gene", "HR_CI", "p", "p_adj")], options = list(pageLength = 15), rownames = FALSE)
  })

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
    gene_tbl <- table(df_one$Gene_for_analysis)
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

  cooccurrence_data <- reactive({
    df <- filtered_data()[, c("Sample", "Gene_for_analysis"), drop = FALSE]
    df <- df[!duplicated(df), , drop = FALSE]
    df$Gene <- df$Gene_for_analysis
    df <- df[, c("Sample", "Gene"), drop = FALSE]
    genes <- names(which(table(as.character(df$Gene)) >= 15))
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
    if (is.null(cooc$matrix) || nrow(cooc$matrix) < 3) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data"))
    mat <- cooc$matrix
    mat_log <- log2(mat + 0.01)
    mat_log[mat_log > 2] <- 2
    mat_log[mat_log < -2] <- -2
    ord <- hclust(dist(mat_log))$order
    mat_ord <- mat_log[ord, ord]
    df_plot <- data.frame(Gene1 = rep(rownames(mat_ord), ncol(mat_ord)),
      Gene2 = rep(colnames(mat_ord), each = nrow(mat_ord)),
      log2OR = as.vector(mat_ord))
    df_plot$Gene1 <- factor(df_plot$Gene1, levels = rownames(mat_ord))
    df_plot$Gene2 <- factor(df_plot$Gene2, levels = colnames(mat_ord))
    ggplot(df_plot, aes(x = Gene1, y = Gene2, fill = log2OR)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
      labs(x = NULL, y = NULL, fill = "log2(OR)", title = "Co-mutation\nOdds Ratio") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_text(size = 14), legend.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 16))
  })

  output$cooccurrence_table <- renderDT({
    cooc <- cooccurrence_data()
    if (is.null(cooc$pairs) || nrow(cooc$pairs) == 0) return(datatable(data.frame(pair = character(), odds_ratio = numeric(), n_cooccur = integer(), p = character(), q = character()), options = list(pageLength = 10), rownames = FALSE))
    df <- cooc$pairs
    df$pair <- paste(df$gene1, "+", df$gene2)
    df$odds_ratio <- round(df$odds_ratio, 2)
    df$p <- format_pval_display(df$p)
    df$q <- format_pval_display(df$q)
    df <- df[order(-df$n_cooccur), , drop = FALSE]
    df <- head(df, 30)
    datatable(df[, c("pair", "odds_ratio", "n_cooccur", "p", "q")], options = list(pageLength = 10), rownames = FALSE)
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
            title = as.character(sig$genotype[i]), xlab = "Years", palette = scale_vals, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.1))
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

  output$vaf_scatter_plot <- renderPlot({
    g1 <- input$vaf_gene1
    g2 <- input$vaf_gene2
    if (is.null(g1) || g1 == "" || is.null(g2) || g2 == "" || g1 == g2) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select two different genes"))
    }
    req(df <- filtered_data())
    df <- df[as.character(df$Gene_for_analysis) %in% c(g1, g2), c("Sample", "Gene_for_analysis", "VAF"), drop = FALSE]
    df$Gene <- df$Gene_for_analysis
    df <- df[, c("Sample", "Gene", "VAF"), drop = FALSE]
    df <- df[!is.na(df$VAF) & df$VAF > 0, , drop = FALSE]
    wide <- reshape(df, idvar = "Sample", timevar = "Gene", direction = "wide")
    colnames(wide) <- gsub("^VAF\\.", "", colnames(wide))
    if (!g1 %in% colnames(wide) || !g2 %in% colnames(wide)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No overlap"))
    wide <- wide[complete.cases(wide[, c(g1, g2), drop = FALSE]), , drop = FALSE]
    if (nrow(wide) < 5) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient samples"))
    diff <- wide[[g1]] - wide[[g2]]
    wide$Clonality <- ifelse(diff > 5, paste(g1, "first"), ifelse(diff < -5, paste(g2, "first"), "Ambiguous"))
    color_vals <- c("Ambiguous" = "#999999", setNames(c("#01665e", "#8c510a"), c(paste(g1, "first"), paste(g2, "first"))))
    ggplot(wide, aes(x = .data[[g2]], y = .data[[g1]], color = Clonality)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_abline(slope = 1, intercept = 5, linetype = "dashed", color = "gray50") +
      geom_abline(slope = 1, intercept = -5, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = color_vals) +
      labs(x = paste(g2, "VAF (%)"), y = paste(g1, "VAF (%)"), color = "Clonal order") +
      theme_minimal(base_size = 12)
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

  # Mut vs WT t-test for all gene-inhibitor pairs (n_mut >= 10) for Drug Sensitivity tab heatmap
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
    # Fallback inline if helper not loaded
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    auc_wide <- b$auc[b$auc$Sample %in% allowed, , drop = FALSE]
    mut_wide <- b$mutations[b$mutations$Sample %in% allowed, , drop = FALSE]
    if (nrow(auc_wide) == 0 || nrow(mut_wide) == 0) return(NULL)
    genes <- names(which(table(mut_wide$Gene) >= 10))
    inhibitors <- unique(auc_wide$inhibitor)
    if (length(genes) == 0 || length(inhibitors) == 0) return(NULL)
    res_list <- list()
    for (g in genes) {
      mut_samples <- unique(mut_wide$Sample[mut_wide$Gene == g])
      if (length(mut_samples) < 10) next
      for (inh in inhibitors) {
        sub <- auc_wide[auc_wide$inhibitor == inh, c("Sample", "auc"), drop = FALSE]
        sub$mut <- sub$Sample %in% mut_samples
        mut_auc <- sub$auc[sub$mut]
        wt_auc <- sub$auc[!sub$mut]
        n_mut <- length(mut_auc)
        n_wt <- length(wt_auc)
        if (n_mut < 10) next
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
    p(em("Select AML subset. All panels below use this subset."))
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
    ggplot(plot_df, aes(x = Inhibitor, y = Gene, fill = delta_AUC, size = VAF_range, label = star)) +
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

  output$drug_mut_wt_heatmap <- renderPlot({
    mut_wt <- drug_mut_wt_all()
    if (is.null(mut_wt) || nrow(mut_wt) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No Mut vs WT data (need ≥10 mut samples per gene-inhibitor)") + theme_void())
    mut_wt <- mut_wt[!is.na(mut_wt$p_value) & mut_wt$p_value < 0.05, , drop = FALSE]
    if (nrow(mut_wt) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No gene-inhibitor pairs with p < 0.05") + theme_void())
    mut_wt$star <- ifelse(!is.na(mut_wt$q_value) & mut_wt$q_value < 0.1, "*", "")
    gene_meds <- aggregate(delta_AUC_mut_wt ~ Gene, data = mut_wt, FUN = median, na.rm = TRUE)
    mut_wt$Gene <- factor(mut_wt$Gene, levels = rev(gene_meds$Gene[order(gene_meds$delta_AUC_mut_wt)]))
    inh_meds <- aggregate(delta_AUC_mut_wt ~ Inhibitor, data = mut_wt, FUN = median, na.rm = TRUE)
    mut_wt$Inhibitor <- factor(mut_wt$Inhibitor, levels = inh_meds$Inhibitor[order(inh_meds$delta_AUC_mut_wt)])
    ggplot(mut_wt, aes(x = Inhibitor, y = Gene, fill = delta_AUC_mut_wt, label = star)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(size = 4, color = "black", fontface = "bold") +
      scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0, name = expression(Delta~"AUC (mut-WT)")) +
      labs(x = "Inhibitor", y = NULL, subtitle = "Delta AUC = mean AUC (mut) - mean AUC (WT). * q < 0.1 (FDR).") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text = element_text(size = 11),
        legend.title = element_text(size = 12), plot.subtitle = element_text(size = 10, color = "gray40"))
  })

  output$drug_scatter_boxplot <- renderPlot({
    d <- drug_scatter_shared()
    if (is.null(d)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select inhibitor and gene") + theme_void())
    auc_sub <- d$auc_sub
    gene <- input$drug_gene
    drug <- input$drug_inhibitor
    mut_col <- if (d$median_mut > d$median_wt) "#2166ac" else "#b2182b"
    wilcox_label <- if (!is.null(d$wilcox)) { pv <- d$wilcox$p.value; if (pv < 0.05) "*" else "ns" } else ""
    ptxt <- if (!is.null(d$tt)) { pv <- d$tt$p.value; if (pv < 0.001) "p < 0.001" else paste0("p = ", format.pval(pv, digits = 2)) } else ""
    stats_label <- paste0("WT: n = ", d$n_wt, ", median = ", round(d$median_wt, 1),
      "  |  Mut: n = ", d$n_mut, ", median = ", round(d$median_mut, 1),
      if (nchar(ptxt) > 0) paste0("  |  ", ptxt) else "",
      if (nchar(wilcox_label) > 0) paste0("  |  Wilcoxon: ", wilcox_label) else "")
    ggplot(auc_sub, aes(x = Status, y = auc, fill = Status)) +
      geom_boxplot(alpha = 0.8, outlier.alpha = 0.5) +
      scale_fill_manual(values = c("WT" = "#4D4D4D", "Mut" = mut_col), guide = "none") +
      coord_cartesian(ylim = d$y_lim) +
      labs(title = paste0(gene, " + ", drug), subtitle = stats_label, x = NULL, y = "Drug AUC") +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), plot.title = element_text(size = 13, face = "bold"), plot.subtitle = element_text(size = 10, color = "gray40"))
  })

  output$drug_scatter_plot <- renderPlot({
    d <- drug_scatter_shared()
    if (is.null(d)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select inhibitor and gene") + theme_minimal(base_size = 14))
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
    ggplot(merged, aes(x = VAF, y = auc)) +
      geom_point(size = 3, alpha = 0.8, color = line_col) +
      geom_smooth(method = "lm", se = TRUE, color = line_col, fill = line_col, alpha = 0.2) +
      coord_cartesian(ylim = d$y_lim) +
      labs(title = paste0(drug, " vs ", gene, " VAF"), x = "VAF (%)", y = "Drug AUC",
        subtitle = paste0("R² = ", r2, ", ", ptxt)) +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13), plot.title = element_text(size = 13, face = "bold"), plot.subtitle = element_text(size = 10, color = "gray40"))
  })

  output$drug_correlation_table <- renderDT({
    corr <- drug_correlations()
    if (is.null(corr) || nrow(corr) == 0) return(datatable(data.frame(Inhibitor = character(), Gene = character(), Direction = character(), R_squared = numeric(), p_value = character(), q_value = character(), LOOCV_RMSE = numeric(), n = integer()), options = list(pageLength = 15)))
    cols <- c("Inhibitor", "Gene", "Direction", "R_squared", "p_value", "q_value", "LOOCV_RMSE", "LOOCV_MSE", "LOOCV_MSE_sd", "n", "VAF_range", "AUC_range")
    cols <- cols[cols %in% colnames(corr)]
    df <- corr[order(corr$p_value), cols, drop = FALSE]
    for (nm in c("R_squared", "LOOCV_RMSE", "LOOCV_MSE", "LOOCV_MSE_sd", "VAF_range", "AUC_range")) {
      if (nm %in% colnames(df) && is.numeric(df[[nm]])) df[[nm]] <- round(df[[nm]], 2)
    }
    if ("p_value" %in% colnames(df)) df$p_value <- format_pval_display(df$p_value)
    if ("q_value" %in% colnames(df)) df$q_value <- format_pval_display(df$q_value)
    datatable(df, filter = "top", options = list(pageLength = 15))
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
    ggplot(loo, aes(x = Inhibitor, y = Gene, fill = RMSE)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_viridis_c(option = "viridis", name = "RMSE") +
      labs(x = "Inhibitor", y = "Gene") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text = element_text(size = 11), legend.title = element_text(size = 12))
  })

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

  # Fetch CIViC (Clinical Interpretation of Variants in Cancer) gene summary
  # API ref: https://griffithlab.github.io/civic-api-docs/
  gene_civic_summary_text <- reactive({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(NULL)
    g_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", g)
    if (!grepl("^[A-Z0-9]+$", g_base)) g_base <- g
    
    if (!has_httr) return(NULL)
    
    extract_description <- function(rec) {
      if (!is.list(rec)) return(NULL)
      errs <- rec$errors
      if (!is.null(errs) && (is.character(errs) && length(errs) > 0 || is.list(errs) && length(errs) > 0)) return(NULL)
      desc <- rec$description
      if (is.null(desc) || length(desc) == 0 || trimws(as.character(desc)) == "") return(NULL)
      as.character(desc)
    }
    
    tryCatch({
      # 1) Try detail endpoint: GET /api/genes/:symbol?identifier_type=entrez_symbol
      url <- paste0("https://civicdb.org/api/genes/", utils::URLencode(g_base, reserved = TRUE), "?identifier_type=entrez_symbol")
      resp <- httr::GET(url, httr::timeout(10), httr::user_agent("Meta-AML-Explorer/1.0"))
      if (httr::status_code(resp) == 200L) {
        rec <- httr::content(resp, as = "parsed", type = "application/json")
        # Response can be single object or array of one
        if (is.data.frame(rec)) rec <- as.list(rec)
        if (is.list(rec) && is.null(rec$records)) {
          desc <- extract_description(rec)
          if (!is.null(desc)) return(desc)
        }
        if (is.list(rec) && identical(names(rec), c("_meta", "records"))) {
          recs <- rec$records
          if (length(recs) > 0) {
            desc <- extract_description(recs[[1]])
            if (!is.null(desc)) return(desc)
          }
        }
        if (is.vector(rec) && length(rec) > 0 && is.list(rec[[1]])) {
          desc <- extract_description(rec[[1]])
          if (!is.null(desc)) return(desc)
        }
      }
      
      # 2) Fallback: search index (paginated) for gene by name or alias
      page <- 1
      count <- 100
      while (page <= 20) {
        idx_url <- sprintf("https://civicdb.org/api/genes?count=%d&page=%d", count, page)
        idx_resp <- httr::GET(idx_url, httr::timeout(10), httr::user_agent("Meta-AML-Explorer/1.0"))
        if (httr::status_code(idx_resp) != 200L) break
        idx <- httr::content(idx_resp, as = "parsed", type = "application/json")
        recs <- idx$records
        if (length(recs) == 0) break
        for (r in recs) {
          name_ok <- identical(toupper(as.character(r$name)), toupper(g_base))
          aliases <- r$aliases
          if (is.null(aliases)) aliases <- character(0)
          if (!is.character(aliases)) aliases <- as.character(aliases)
          alias_ok <- toupper(g_base) %in% toupper(aliases)
          if (name_ok || alias_ok) {
            desc <- extract_description(r)
            if (!is.null(desc)) return(desc)
            # Fetch full detail by CIViC id for description if index record has no description
            detail_url <- paste0("https://civicdb.org/api/genes/", r$id, "?identifier_type=civic_id")
            d_resp <- httr::GET(detail_url, httr::timeout(8), httr::user_agent("Meta-AML-Explorer/1.0"))
            if (httr::status_code(d_resp) == 200L) {
              det <- httr::content(d_resp, as = "parsed", type = "application/json")
              desc <- extract_description(det)
              if (!is.null(desc)) return(desc)
            }
            break
          }
        }
        meta <- idx$`_meta`
        next_link <- if (is.list(meta$links)) meta$links[["next"]] else NULL
        if (is.null(next_link) || identical(next_link, "")) break
        page <- page + 1
      }
      return(NULL)
    }, error = function(e) {
      return(NULL)
    })
  })

  output$gene_civic_summary <- renderUI({
    summary_text <- gene_civic_summary_text()
    if (is.null(summary_text) || summary_text == "") {
      return(NULL)
    }
    return(div(
      style = "margin-top: 10px; padding: 10px; background-color: #f5f5f5; border-left: 3px solid #0066cc;",
      p(strong("CIViC Summary:"), " ", a("(Clinical Interpretation of Variants in Cancer)", href = "https://civicdb.org/", target = "_blank", style = "font-size: 11px; color: #0066cc;"), style = "margin-bottom: 5px;"),
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
        uiOutput("gene_ncbi_summary"),
        uiOutput("gene_civic_summary")
      ),
      tabsetPanel(
        id = "gene_summary_sub_tabs",
        type = "tabs",
        selected = if (is.null(selected_gene_summary_sub_tab()) || !selected_gene_summary_sub_tab() %in% c("clinical", "comut", "vaf", "drug")) "clinical" else selected_gene_summary_sub_tab(),
        tabPanel("Clinical", value = "clinical"),
        tabPanel("Co-mutation", value = "comut"),
        tabPanel("VAF", value = "vaf"),
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
        wellPanel(h4("Clinical associations"), plotOutput("gene_summary_clinical_plot", height = 260)),
        fluidRow(
          column(6, wellPanel(h4("Mutation distribution"), plotOutput("gene_summary_lollipop_plot", height = 380))),
          column(6, wellPanel(h4("Kaplan-Meier: Single gene"), plotOutput("gene_summary_survival", height = 480)))
        )
      ),
      comut = tagList(
        fluidRow(
          column(9, wellPanel(h4("Oncoprint: Co-mutations in samples with this gene mutated"), plotOutput("gene_summary_oncoprint", height = 480))),
          column(3, wellPanel(h4("Co-mutation Odds Ratio"), plotOutput("gene_summary_comut_heatmap", height = 480)))
        ),
        fluidRow(
          column(6, div(style = "height: 700px;", wellPanel(style = "height: 100%;", h4("Pairwise co-mutation survival"), selectInput("gene_summary_gene2", "Second gene:", choices = c("Select..." = "", genes_other)), plotOutput("gene_summary_comut2_plot", height = 500)))),
          column(6, div(style = "height: 700px;", wellPanel(style = "height: 100%;", h4("Triple co-mutation survival"), fluidRow(column(6, selectInput("gene_summary_gene3a", "Second gene:", choices = c("Select..." = "", genes_other))), column(6, selectInput("gene_summary_gene3b", "Third gene:", choices = c("Select..." = "", genes_other))), div(style = "padding-left: 19px; padding-right: 19px;", plotOutput("gene_summary_comut3_plot", height = 520))))))
        )
      ),
      vaf = fluidRow(
        column(6, wellPanel(h4("VAF distribution"), plotOutput("gene_summary_vaf_plot", height = 350))),
        column(6, wellPanel(h4("VAF and survival (MaxStat)"), plotOutput("gene_summary_vaf_survival_plot", height = 350)))
      ),
      drug = wellPanel(
        style = "background-color: #e9ecef;",
        h4("Drug sensitivity (Beat AML)"),
        fluidRow(
          column(6, wellPanel(h5("Mut vs WT volcano (drug AUC)"), plotOutput("gene_summary_drug_volcano", height = 400))),
          column(6, wellPanel(h5("VAF vs AUC summary (this gene)"), plotOutput("gene_summary_drug_dotplot", height = 400)))
        ),
        fluidRow(
          column(6, wellPanel(h5("VAF–AUC correlations (this gene)"), div(style = "height: 300px; overflow-y: auto;", DTOutput("gene_summary_drug_table")))),
          column(6, wellPanel(h5("LOOCV RMSE heatmap (ordered by RMSE)"), plotOutput("gene_summary_drug_loo_heatmap", height = 350)))
        )
      ),
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
      labs(x = NULL, y = NULL, fill = "log2(OR)", title = paste("Co-occurrence with", g)) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  })

  output$gene_summary_lollipop_plot <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void())
    df <- filtered_data()
    if (is.null(df) || nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void())
    df <- df[as.character(df$Gene_for_analysis) == g, , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No mutations for", g)) + theme_void())
    aa_col <- if ("HGVSp_Short" %in% colnames(df)) "HGVSp_Short" else if ("AA_change" %in% colnames(df)) "AA_change" else if ("Protein_Change" %in% colnames(df)) "Protein_Change" else if ("protein" %in% colnames(df)) "protein" else NULL
    if (has_maftools && !is.null(aa_col)) {
      g_base <- gsub("^([A-Z0-9]+)[_-].*$", "\\1", g)
      if (!grepl("^[A-Z0-9]+$", g_base)) g_base <- g
      norm <- normalize_protein_change_for_maf(df[[aa_col]])
      hgvsp <- norm$hgvsp
      has_position <- !is.na(norm$position) & hgvsp != "p.?"
      if (sum(has_position) >= 1) {
        df_maf <- df[has_position & !is.na(df$Sample) & as.character(df$Sample) != "", , drop = FALSE]
        hgvsp_maf <- hgvsp[has_position & !is.na(df$Sample) & as.character(df$Sample) != ""]
        if (nrow(df_maf) < 1) {
          return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No valid samples with AA positions for lollipopPlot") + theme_void())
        }
        vt <- if ("variant_type" %in% colnames(df_maf)) as.character(df_maf$variant_type) else rep(NA_character_, nrow(df_maf))
        vc <- rep("Missense_Mutation", length(vt))
        vc[vt %in% c("Deletion", "Frame_Shift_Del", "Frame_Shift_Del")] <- "Frame_Shift_Del"
        vc[vt %in% c("Indel", "In_Frame_Del")] <- "In_Frame_Del"
        vc[vt %in% c("ITD", "Insertion", "In_Frame_Ins")] <- "In_Frame_Ins"
        vc[vt %in% c("Splicing", "Splice_Site")] <- "Splice_Site"
        vc[vt %in% c("PTD")] <- "In_Frame_Ins"
        vc[is.na(vt) | vt %in% c("Other", "Unknown", "")] <- "Unknown"
        vtype <- rep("SNP", length(vt))
        vtype[vt %in% c("ITD", "Insertion", "PTD")] <- "INS"
        vtype[vt %in% c("Deletion", "Indel")] <- "DEL"
        base <- 1000000L
        maf_df <- data.frame(
          Hugo_Symbol = rep(g_base, nrow(df_maf)),
          Chromosome = "1",
          Start_Position = base + seq_len(nrow(df_maf)),
          End_Position = base + seq_len(nrow(df_maf)),
          Reference_Allele = rep("N", nrow(df_maf)),
          Tumor_Seq_Allele1 = rep("N", nrow(df_maf)),
          Tumor_Seq_Allele2 = rep("N", nrow(df_maf)),
          Variant_Classification = vc,
          Variant_Type = vtype,
          Tumor_Sample_Barcode = as.character(df_maf$Sample),
          HGVSp_Short = hgvsp_maf,
          stringsAsFactors = FALSE
        )
        maf_obj <- tryCatch(maftools::read.maf(maf = maf_df, verbose = FALSE), error = function(e) {
          warning("read.maf failed: ", conditionMessage(e))
          NULL
        })
        if (!is.null(maf_obj) && inherits(maf_obj, "MAF")) {
          ok <- tryCatch({
            maftools::lollipopPlot(maf_obj, gene = g_base, AACol = "HGVSp_Short", showMutationRate = FALSE)
            TRUE
          }, error = function(e) {
            warning("lollipopPlot failed for ", g_base, ": ", conditionMessage(e))
            FALSE
          })
          if (isTRUE(ok)) return(invisible(NULL))
        }
      }
    }
    if (!has_maftools) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install maftools to render lollipop plots") + theme_void())
    }
    if (is.null(aa_col)) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No protein-change column available for lollipopPlot (need AA_change / HGVSp_Short)") + theme_void())
    }
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Unable to render lollipopPlot for this selection") + theme_void()
  })

  output$gene_summary_survival <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    surv_df <- survival_data()
    mut_samples <- unique(surv_df$Sample[surv_df$Gene_for_analysis == g])
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    surv_uniq$Mutation <- ifelse(surv_uniq$Sample %in% mut_samples, paste0(g, " mut"), "WT")
    if (length(unique(surv_uniq$Mutation)) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient groups"))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Mutation, data = surv_uniq)
    if (has_survminer) {
      p <- survminer::ggsurvplot(fit, data = surv_uniq, risk.table = TRUE, pval = TRUE, title = paste("Survival by", g), xlab = "Years", palette = c("#8B0000", "#4D4D4D"), legend.labs = gsub("^Mutation=", "", names(fit$strata)))
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
      p <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE, title = "Pairwise co-mutation survival", xlab = "Years", palette = pal_vec, legend.labs = gsub("^Group=", "", names(fit$strata)))
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
      p <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE, title = "Triple co-mutation survival", xlab = "Years", palette = pal_vec, legend.labs = gsub("^Group=", "", names(fit$strata)), risk.table.height = 0.4)
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install survminer") + theme_void()
    }
  })

  output$gene_summary_clinical_plot <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    df <- filtered_data()
    clin_vars <- c("Age", "BM_blast_percent", "PB_blast_percent", "WBC", "Hemoglobin", "Platelet")
    clin_vars <- clin_vars[clin_vars %in% colnames(df)]
    if (length(clin_vars) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No clinical variables in data"))
    samples_with_g <- unique(df$Sample[as.character(df$Gene_for_analysis) == g])
    df_one <- df[!duplicated(df$Sample), c("Sample", clin_vars), drop = FALSE]
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
    df <- filtered_data()
    df <- df[as.character(df$Gene_for_analysis) == g & !is.na(df$VAF) & df$VAF > 0, , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No VAF data"))
    ggplot(df, aes(x = VAF)) +
      geom_histogram(bins = 25, fill = "#374e55", alpha = 0.8) +
      labs(x = "VAF (%)", y = "Count", title = paste("VAF distribution for", g)) +
      theme_minimal(base_size = 14)
  })

  output$gene_summary_vaf_survival_plot <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    df <- filtered_data()
    df <- df[as.character(df$Gene_for_analysis) == g & !is.na(df$VAF) & df$VAF > 0, c("Sample", "VAF"), drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No VAF data"))
    vaf_agg <- aggregate(VAF ~ Sample, data = df, max)
    surv_df <- survival_data()
    surv_df <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    merge_df <- merge(vaf_agg, surv_df, by = "Sample", all = FALSE)
    if (nrow(merge_df) < 20 || !has_maxstat) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Need ≥20 samples with VAF and survival; maxstat required") + theme_void())
    mst <- tryCatch(maxstat::maxstat.test(Surv(Time_to_OS, as.numeric(Censor)) ~ VAF, data = merge_df, smethod = "LogRank", minprop = 0.2, maxprop = 0.8), error = function(e) NULL)
    if (is.null(mst)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "MaxStat failed") + theme_void())
    merge_df$Group <- factor(ifelse(merge_df$VAF > mst$estimate, "High VAF", "Low VAF"), levels = c("High VAF", "Low VAF"))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = merge_df)
    if (has_survminer) {
      legend_labs <- gsub("^Group=", "", names(fit$strata))
      # Unnamed palette: first = High VAF (red), second = Low VAF (grey); Group factor ensures this order
      p <- survminer::ggsurvplot(fit, data = merge_df, risk.table = TRUE, pval = TRUE, title = paste0(g, ": High vs Low VAF (", round(mst$estimate, 1), "%)"), xlab = "Years", palette = c("#8B0000", "#4D4D4D"), legend.labs = legend_labs)
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
    subset <- input$subset
    if (is.null(subset) || subset == "All") subset <- "de_novo"
    subset_map <- c("De novo" = "de_novo", "Secondary" = "secondary", "Relapse" = "relapse", "Therapy" = "therapy", "Other" = "other")
    subset <- if (subset %in% names(subset_map)) subset_map[subset] else tolower(subset)
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    base_genes <- gene_to_beataml_bases(g)
    mut_samples <- unique(b$mutations$Sample[b$mutations$Gene %in% base_genes & b$mutations$Sample %in% allowed])
    auc_wide <- b$auc[b$auc$Sample %in% allowed, , drop = FALSE]
    if (nrow(auc_wide) == 0) return(NULL)
    inhibitors <- unique(auc_wide$inhibitor)
    res_list <- list()
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
      res_list[[length(res_list) + 1L]] <- data.frame(Inhibitor = inh, delta_AUC_mut_wt = delta, p_value = pval, n_mut = n_mut, n_wt = n_wt, stringsAsFactors = FALSE)
    }
    if (length(res_list) == 0) return(NULL)
    do.call(rbind, res_list)
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
      theme(legend.position = "top")
    if (nrow(df_label) > 0) {
      if (has_ggrepel) {
        p <- p + ggrepel::geom_label_repel(data = df_label, aes(x = delta_AUC_mut_wt, y = neglog10p, label = Inhibitor), size = 5, show.legend = FALSE, max.overlaps = Inf)
      } else {
        p <- p + geom_text(data = df_label, aes(x = delta_AUC_mut_wt, y = neglog10p, label = Inhibitor), hjust = -0.1, size = 5, show.legend = FALSE)
      }
    }
    p
  })

  output$gene_summary_drug_dotplot <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void())
    corr <- gene_summary_drug_correlations()
    if (is.null(corr) || nrow(corr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No drug data (Beat AML)") + theme_void())
    if ("n" %in% colnames(corr)) corr <- corr[!is.na(corr$n) & corr$n >= 10, , drop = FALSE]
    base_genes <- gene_to_beataml_bases(g)
    sub <- corr[corr$Gene %in% base_genes, , drop = FALSE]
    if (nrow(sub) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No VAF–AUC correlations for", g)) + theme_void())
    sub <- sub[!is.na(sub$p_value) & sub$p_value < 0.05, , drop = FALSE]
    if (nrow(sub) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No VAF–AUC correlations with p < 0.05 for", g)) + theme_void())
    inhibitors <- unique(sub$Inhibitor)
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No Beat AML data") + theme_void())
    subset <- input$subset
    if (is.null(subset) || subset == "All") subset <- "de_novo"
    subset_map <- c("De novo" = "de_novo", "Secondary" = "secondary", "Relapse" = "relapse", "Therapy" = "therapy", "Other" = "other")
    subset <- if (subset %in% names(subset_map)) subset_map[subset] else tolower(subset)
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
    pval_txt <- ifelse(label_df$p_value < 0.001, "< 0.001", format.pval(label_df$p_value, digits = 2))
    label_df$label <- paste0("R² = ", round(label_df$R_squared, 3), ", p ", pval_txt, label_df$star)
    trend_colors <- c("Sensitive" = "#b2182b", "Resistant" = "#2166ac", "Unknown" = "gray50")
    ggplot(scatter_df, aes(x = VAF, y = auc, color = Trend)) +
      geom_point(size = 2.5, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.8) +
      scale_color_manual(values = trend_colors, name = "Trend", guide = "none") +
      facet_wrap(~ Inhibitor, scales = "free", ncol = min(3, length(inhibitors))) +
      labs(x = "VAF (%)", y = "Drug AUC", title = paste("VAF vs AUC scatterplots:", g, "(p < 0.05)")) +
      theme_minimal(base_size = 12) +
      theme(
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "#f0f0f0", color = NA),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(size = 13, face = "bold")
      ) +
      geom_text(aes(x = Inf, y = Inf, label = label), hjust = 1.1, vjust = 1.5, size = 3, inherit.aes = FALSE, data = label_df)
  })

  output$gene_summary_drug_table <- renderDT({
    empty_drug <- datatable(data.frame(Inhibitor = character(), Gene = character(), Direction = character(), R_squared = numeric(), p_value = character(), q_value = character(), LOOCV_RMSE = numeric(), n = integer()), options = list(pageLength = 10))
    g <- input$gene_summary
    if (is.null(g) || g == "") return(empty_drug)
    corr <- gene_summary_drug_correlations()
    if (is.null(corr) || nrow(corr) == 0) return(empty_drug)
    if ("n" %in% colnames(corr)) corr <- corr[!is.na(corr$n) & corr$n >= 10, , drop = FALSE]
    base_genes <- gene_to_beataml_bases(g)
    sub <- corr[corr$Gene %in% base_genes, , drop = FALSE]
    if (nrow(sub) == 0) return(empty_drug)
    cols <- c("Inhibitor", "Gene", "Direction", "R_squared", "p_value", "q_value", "LOOCV_RMSE", "LOOCV_MSE", "n", "VAF_range", "AUC_range")
    cols <- cols[cols %in% colnames(sub)]
    df <- sub[order(sub$p_value), cols, drop = FALSE]
    for (nm in c("R_squared", "LOOCV_RMSE", "LOOCV_MSE", "VAF_range", "AUC_range")) {
      if (nm %in% colnames(df) && is.numeric(df[[nm]])) df[[nm]] <- round(df[[nm]], 2)
    }
    if ("p_value" %in% colnames(df)) df$p_value <- format_pval_display(df$p_value)
    if ("q_value" %in% colnames(df)) df$q_value <- format_pval_display(df$q_value)
    datatable(df, options = list(pageLength = 10))
  })

  output$gene_summary_drug_loo_heatmap <- renderPlot({
    g <- input$gene_summary
    if (is.null(g) || g == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene") + theme_void())
    corr <- gene_summary_drug_correlations()
    if (is.null(corr) || nrow(corr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No drug data (Beat AML)") + theme_void())
    if ("n" %in% colnames(corr)) corr <- corr[!is.na(corr$n) & corr$n >= 10, , drop = FALSE]
    base_genes <- gene_to_beataml_bases(g)
    sub <- corr[corr$Gene %in% base_genes & !is.na(corr$LOOCV_RMSE), , drop = FALSE]
    if (nrow(sub) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("No LOOCV RMSE data for", g)) + theme_void())
    inh_ord <- sub$Inhibitor[order(sub$LOOCV_RMSE)]
    sub$Inhibitor <- factor(sub$Inhibitor, levels = inh_ord)
    ggplot(sub, aes(x = Inhibitor, y = Gene, fill = LOOCV_RMSE)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_viridis_c(option = "viridis", name = "RMSE") +
      labs(x = "Inhibitor", y = NULL, title = paste("LOOCV RMSE:", g)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  })

  output$data_table <- renderDT({
    df <- req(filtered_data())
    cols <- c("Sample", "Gene_for_analysis", "Gene", "VAF", "variant_type", "AA_change", "Gene_Group", "Subset", "Cohort", "Age", "Sex", "Risk", "Time_to_OS", "Censor", "mutation_category")
    cols <- cols[cols %in% colnames(df)]
    df <- df[, cols, drop = FALSE]
    for (nm in c("VAF", "Age", "Time_to_OS", "Censor")) {
      if (nm %in% colnames(df) && is.numeric(df[[nm]])) df[[nm]] <- round(df[[nm]], 2)
    }
    datatable(df, filter = "top", options = list(pageLength = 25, scrollX = TRUE))
  })
}

shinyApp(ui = ui, server = server)
