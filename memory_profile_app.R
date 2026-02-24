# Memory profiling for Meta AML Shiny app
# Run from Meta_AML_Shiny directory: Rscript memory_profile_app.R
# Reports memory (MB) at each stage to find hotspots and accumulation.

options(stringsAsFactors = FALSE)

get_mem_mb <- function() {
  g <- gc(verbose = FALSE)
  # gc() matrix: col 1 = used (cells/bytes), col 2 = same in Mb
  if (is.matrix(g) && ncol(g) >= 2) sum(g[, 2]) else sum(g[, 1]) / 1024 / 1024
}

get_mem_max_mb <- function() {
  g <- gc(verbose = FALSE)
  if (is.matrix(g) && ncol(g) >= 7) sum(g[, 7], na.rm = TRUE) else 0
}

message("=== Meta AML Shiny memory profile ===\n")

# 1) Baseline after loading base R
gc(reset = TRUE)
mem_base <- get_mem_mb()
message(sprintf("1. After gc(reset=TRUE): %.1f MB", mem_base))

# 2) Load app libraries (same as app.R)
library(shiny)
library(ggplot2)
library(survival)
library(DT)
library(scales)
mem_libs <- get_mem_mb()
message(sprintf("2. After loading shiny/ggplot2/survival/DT/scales: %.1f MB (delta +%.1f)", mem_libs, mem_libs - mem_base))

# Optional packages (load if present)
for (pkg in c("survminer", "gridExtra", "patchwork", "ComplexHeatmap", "maxstat", "bestglm", "ggrepel", "UpSetR")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}
mem_opt <- get_mem_mb()
message(sprintf("3. After optional packages: %.1f MB (delta +%.1f)", mem_opt, mem_opt - mem_libs))

# 4) Load cohort RDS (same path as app)
if (!file.exists("AML_Meta_Cohort_v2.rds")) {
  message("AML_Meta_Cohort_v2.rds not found. Create it with build_AML_Meta_Cohort_v2.R or run from app directory.")
  quit(save = "no", status = 1)
}
AML_Meta_Cohort <- as.data.frame(readRDS("AML_Meta_Cohort_v2.rds"), stringsAsFactors = FALSE)
mem_rds <- get_mem_mb()
n_row <- nrow(AML_Meta_Cohort)
n_col <- ncol(AML_Meta_Cohort)
message(sprintf("4. After readRDS(AML_Meta_Cohort_v2.rds): %.1f MB (delta +%.1f) [%d rows x %d cols]",
                mem_rds, mem_rds - mem_opt, n_row, n_col))

# 5) Normalize (minimal replicate of app's normalize - full version trims columns etc.)
# For memory we care about: copy of data frame
normalize_AML_Meta_Cohort <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  if (!"Sample" %in% colnames(d)) d$Sample <- paste0("S", seq_len(nrow(d)))
  if ("Gene" %in% colnames(d)) d$Gene <- as.character(d$Gene)
  if ("Cohort" %in% colnames(d)) {
    d$Cohort <- as.character(d$Cohort)
    d$Cohort[d$Cohort == "Beat_AML"] <- "Beat AML"
  }
  if ("Subset" %in% colnames(d)) d$Subset <- as.character(d$Subset)
  if ("Gene_Group" %in% colnames(d)) {
    d$Gene_for_analysis <- ifelse(is.na(d$Gene_Group) | as.character(d$Gene_Group) == "", as.character(d$Gene), as.character(d$Gene_Group))
  } else {
    d$Gene_for_analysis <- as.character(d$Gene)
  }
  d
}
AML_Meta_Cohort <- normalize_AML_Meta_Cohort(AML_Meta_Cohort)
mem_norm <- get_mem_mb()
message(sprintf("5. After normalize_AML_Meta_Cohort: %.1f MB (delta +%.1f)", mem_norm, mem_norm - mem_rds))

# 6) Temp CSV write/read (app does this - duplicates memory briefly)
tmp2 <- tempfile(fileext = ".csv")
utils::write.csv(AML_Meta_Cohort, tmp2, row.names = FALSE)
AML_Meta_Cohort <- utils::read.csv(tmp2, stringsAsFactors = FALSE)
unlink(tmp2)
mem_csv <- get_mem_mb()
message(sprintf("6. After temp CSV write/read (app duplication): %.1f MB (delta +%.1f)", mem_csv, mem_csv - mem_norm))

# 7) Simulate filtered_data() for "All" (subset=De novo, cohort=All, karyotype=All)
prepare_data <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$Gene <- as.character(df$Gene)
  if ("variant_type" %in% colnames(df)) df$variant_type <- as.character(df$variant_type)
  df
}
df_all <- AML_Meta_Cohort
if ("Subset" %in% colnames(df_all)) df_all <- df_all[as.character(df_all$Subset) == "De novo", , drop = FALSE]
df_all <- prepare_data(df_all)
mem_filter_all <- get_mem_mb()
message(sprintf("7. Filtered data (De novo, All cohorts): %.1f MB (delta +%.1f) [%d rows]",
                mem_filter_all, mem_filter_all - mem_csv, nrow(df_all)))

# 8) Second filter (one cohort) - simulates switching cohort; app caches both tabs
cohorts <- sort(unique(as.character(na.omit(AML_Meta_Cohort$Cohort))))
one_cohort <- if (length(cohorts) > 0) cohorts[1] else "All"
df_cohort <- AML_Meta_Cohort
if ("Subset" %in% colnames(df_cohort)) df_cohort <- df_cohort[as.character(df_cohort$Subset) == "De novo", , drop = FALSE]
if (one_cohort != "All") df_cohort <- df_cohort[as.character(df_cohort$Cohort) == one_cohort, , drop = FALSE]
df_cohort <- prepare_data(df_cohort)
# Simulate keeping two caches (analyses + meta_aml4)
cached_analyses <- df_all
cached_meta4 <- df_cohort
mem_two_caches <- get_mem_mb()
message(sprintf("8. With two cached filter results (analyses + meta_aml4): %.1f MB (delta +%.1f)",
                mem_two_caches, mem_two_caches - mem_filter_all))

# 9) Co-occurrence matrix (genes x genes) - can be large
genes <- unique(as.character(df_all$Gene_for_analysis))
genes <- genes[!is.na(genes) & genes != ""]
n_genes <- length(genes)
# App builds OR matrix for pairs; approximate with a dense n_genes x n_genes matrix
if (n_genes > 1) {
  or_mat <- matrix(1, n_genes, n_genes)
  rownames(or_mat) <- colnames(or_mat) <- genes
}
mem_cooc <- get_mem_mb()
message(sprintf("9. After co-occurrence matrix (%d x %d): %.1f MB (delta +%.1f)",
                n_genes, n_genes, mem_cooc, mem_cooc - mem_two_caches))

# 10) Oncoprint-like: wide matrix samples x genes (can be huge)
samples <- unique(as.character(df_all$Sample))
n_samp <- length(samples)
# App builds a matrix for ComplexHeatmap/oncoprint - not always full dense, but worst case
# Approximate: one value per sample-gene for top genes
top_genes <- names(sort(table(df_all$Gene_for_analysis), decreasing = TRUE))[1:min(100, n_genes)]
if (length(top_genes) == 0) top_genes <- genes[1:min(50, n_genes)]
oncoprint_cells <- length(samples) * length(top_genes)
mem_onco <- get_mem_mb()
message(sprintf("10. After oncoprint prep (%d samples x %d genes): %.1f MB (delta +%.1f)",
                n_samp, length(top_genes), mem_onco, mem_onco - mem_cooc))

# 11) Precomputed drug (if exists)
precomputed_drug <- NULL
path_drug <- file.path(getwd(), "drug_sensitivity_precomputed.rds")
if (file.exists(path_drug)) {
  precomputed_drug <- readRDS(path_drug)
  mem_drug <- get_mem_mb()
  message(sprintf("11. After drug_sensitivity_precomputed.rds: %.1f MB (delta +%.1f)",
                  mem_drug, mem_drug - mem_onco))
} else {
  message("11. drug_sensitivity_precomputed.rds not found (skip)")
}

# 12) Beat AML2 load (if present)
beataml2_script <- file.path(getwd(), "load_beataml2.R")
if (file.exists(beataml2_script)) {
  mem_before_beat <- get_mem_mb()
  source(beataml2_script, local = FALSE)
  beataml <- tryCatch(load_beataml2(), error = function(e) list(ok = FALSE))
  if (is.list(beataml) && identical(beataml$ok, TRUE)) {
    mem_beat <- get_mem_mb()
    message(sprintf("12. After load_beataml2(): %.1f MB (delta +%.1f)", mem_beat, mem_beat - mem_before_beat))
  }
}

# Summary
message("\n=== Summary ===")
message(sprintf("Peak memory (max used): %.1f MB", get_mem_max_mb()))
message(sprintf("Current memory: %.1f MB", get_mem_mb()))
message("\nLargest contributors (deltas):")
message("  - Cohort RDS load + normalize + CSV round-trip: base data")
message("  - Caching two filtered datasets (Analyses + Meta AML4) doubles filtered data in memory")
message("  - Co-occurrence matrix and oncoprint builds add more; precomputed_drug and Beat AML if used.")
message("\nRecommendations:")
message("  1. Avoid keeping two full caches when possible (e.g. clear other tab cache when switching).")
message("  2. Consider removing temp CSV write/read in app.R if not required for compatibility.")
message("  3. Limit genes in co-occurrence/oncoprint for very large cohorts.")
message("  4. Use precomputed_drug to avoid recomputing drug sensitivity on each load.")
