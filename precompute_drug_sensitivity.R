# =============================================================================
# Precompute Drug Sensitivity results for All Gene Associations tab.
# Run once (e.g. before deploy) to generate drug_sensitivity_precomputed.rds.
# The Shiny app will load this file and use it instead of computing on the fly,
# reducing memory use and avoiding OOM on shinyapps.io free tier.
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

b_wave12 <- subset_beataml_to_wave12(b)
if (!b_wave12$ok) b_wave12 <- b  # fallback

subsets <- c("All", "de_novo", "secondary")
precomputed <- list(meta_aml4 = list(), analyses = list())

message("Precomputing for Meta AML4 (all waves)...")
for (sub in subsets) {
  message("  subset: ", sub)
  precomputed$meta_aml4[[sub]] <- list(
    mut_wt = compute_mut_wt_all(b, sub),
    correlations = compute_drug_vaf_correlations(b, subset = sub),
    loo = compute_drug_vaf_loo(b, subset = sub)
  )
}

message("Precomputing for Benard/analyses (waves 1+2)...")
for (sub in subsets) {
  message("  subset: ", sub)
  precomputed$analyses[[sub]] <- list(
    mut_wt = compute_mut_wt_all(b_wave12, sub),
    correlations = compute_drug_vaf_correlations(b_wave12, subset = sub),
    loo = compute_drug_vaf_loo(b_wave12, subset = sub)
  )
}

out_path <- "drug_sensitivity_precomputed.rds"
saveRDS(precomputed, out_path)
message("Saved to ", out_path)
