# =============================================================================
# Build AML_Meta_Cohort_v2.rds by merging aggregated MAF + clinical.
# Run from Meta_AML_Shiny after build_aggregated_clinical.R and build_aggregated_maf.R.
# Output: AML_Meta_Cohort_v2.rds (used by the app Meta AML4 tab)
# =============================================================================

setwd(if (exists("ofile")) dirname(ofile) else getwd())
maf_path <- "Meta_AML_aggregated_mutations.maf.tsv"
clin_path <- "Meta_AML_aggregated_clinical.tsv"
if (!file.exists(maf_path)) stop("Run build_aggregated_maf.R first.")
if (!file.exists(clin_path)) stop("Run build_aggregated_clinical.R first.")

maf <- read.delim(maf_path, stringsAsFactors = FALSE, check.names = FALSE)
clin <- read.delim(clin_path, stringsAsFactors = FALSE, check.names = FALSE)

# Merge: sample id in MAF is Tumor_Sample_Barcode, in clinical is Patient_ID
maf$Sample <- as.character(maf$Tumor_Sample_Barcode)
clin$Patient_ID <- as.character(clin$Patient_ID)
merged <- merge(maf, clin, by.x = "Sample", by.y = "Patient_ID", all.x = TRUE, suffixes = c("", "_clin"))

# Drop duplicate columns from clinical (Cohort_clin etc.) if any
dup_suffix <- "_clin"
dup_cols <- grep(paste0(dup_suffix, "$"), colnames(merged), value = TRUE)
if (length(dup_cols) > 0) merged <- merged[, !colnames(merged) %in% dup_cols, drop = FALSE]

out_path <- "AML_Meta_Cohort_v2_for_download.rds"
saveRDS(merged, out_path)
message("Written ", nrow(merged), " rows to ", out_path)
message("Cohort counts:")
print(table(merged$Cohort, useNA = "ifany"))
if ("ELN_2022_Risk" %in% colnames(merged))
  message("ELN_2022_Risk present: ", sum(!is.na(merged$ELN_2022_Risk)), " non-NA")
