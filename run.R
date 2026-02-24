# Meta AML Explorer - Launcher
# Run this script to start the Shiny app. The app runs only from Meta_AML_Shiny.
#
# From RStudio: Open app.R and click "Run App", or Source this file
# From terminal: Rscript run.R  (run from the Meta_AML_Shiny directory)

cohort_ok <- file.exists("AML_Meta_Cohort_v2.rds") || file.exists("AML_Meta_Cohort.rds") || file.exists("AML_Meta_Cohort.RData")
if (!cohort_ok) {
  stop("Cohort data not found. Run from the Meta_AML_Shiny directory (containing app.R and AML_Meta_Cohort_v2.rds or AML_Meta_Cohort.rds).")
}

message("Starting Meta AML Explorer (Meta_AML_Shiny)...")
message("The app will open in your browser. Close the browser tab to stop the app.")
message("")

shiny::runApp(".", launch.browser = TRUE)
