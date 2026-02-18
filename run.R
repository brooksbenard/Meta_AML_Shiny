# Meta AML Explorer - Launcher
# Run this script to start the Shiny app. The app runs only from Meta_AML_Shiny.
#
# From RStudio: Open app.R and click "Run App", or Source this file
# From terminal: Rscript run.R  (run from the Meta_AML_Shiny directory)

if (!file.exists("final_data_matrix.RData")) {
  stop("Data file not found. Run from the Meta_AML_Shiny directory (containing app.R and final_data_matrix.RData).")
}

message("Starting Meta AML Explorer (Meta_AML_Shiny)...")
message("The app will open in your browser. Close the browser tab to stop the app.")
message("")

shiny::runApp(".", launch.browser = TRUE)
