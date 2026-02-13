# =============================================================================
# Download BeatAML2 data for the Drug Sensitivity tab
# Run once from the Meta_AML_Shiny directory: source("setup_beataml2.R")
# Data from https://biodev.github.io/BeatAML2/
# =============================================================================

base_url <- "https://github.com/biodev/beataml2.0_data/raw/main"
data_dir <- file.path(getwd(), "beataml2_data")
dir.create(data_dir, showWarnings = FALSE)

files <- c(
  "beataml_wes_wv1to4_mutations_dbgap.txt" = "mutations.txt",
  "beataml_probit_curve_fits_v4_dbgap.txt" = "inhibitor_auc.txt",
  "beataml_wv1to4_clinical.xlsx" = "clinical.xlsx"
)

for (i in seq_along(files)) {
  dest <- file.path(data_dir, files[i])
  url <- paste(base_url, names(files)[i], sep = "/")
  if (!file.exists(dest) || file.info(dest)$size < 1000) {
    message("Downloading ", names(files)[i], " ...")
    tryCatch(
      download.file(url, dest, mode = "wb", quiet = TRUE),
      error = function(e) message("Failed: ", e$message)
    )
  } else {
    message("Already exists: ", dest)
  }
}

message("Done. BeatAML2 data in: ", data_dir)
