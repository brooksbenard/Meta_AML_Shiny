# =============================================================================
# Download AML-SG data from gerstung-lab/AML-multistage
# Run once from the Meta_AML_Shiny directory: source("setup_amlsg.R")
# Data from https://github.com/gerstung-lab/AML-multistage/tree/master/data
# Paper: Papaemmanuil et al. NEJM 2016; Precision oncology for AML (knowledge bank)
# =============================================================================

base_url <- "https://raw.githubusercontent.com/gerstung-lab/AML-multistage/master/data"
data_dir <- file.path(getwd(), "amlsg_data")
dir.create(data_dir, showWarnings = FALSE)

files <- c(
  "AMLSG_Classification.txt" = "AMLSG_Classification.txt",
  "AMLSG_Clinical_Anon.RData" = "AMLSG_Clinical_Anon.RData",
  "AMLSG_FLT3ITD.txt" = "AMLSG_FLT3ITD.txt",
  "AMLSG_Genetic.txt" = "AMLSG_Genetic.txt",
  "AMLSG_Karyotypes.txt" = "AMLSG_Karyotypes.txt",
  "README.md" = "README.md"
)

for (i in seq_along(files)) {
  dest <- file.path(data_dir, files[i])
  url <- paste(base_url, names(files)[i], sep = "/")
  if (!file.exists(dest) || file.info(dest)$size < 100) {
    message("Downloading ", names(files)[i], " ...")
    tryCatch(
      download.file(url, dest, mode = "wb", quiet = TRUE),
      error = function(e) message("Failed: ", e$message)
    )
  } else {
    message("Already exists: ", dest)
  }
}

message("Done. AML-SG data in: ", data_dir)
