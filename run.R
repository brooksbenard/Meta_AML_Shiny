# Meta_AML Explorer - Launcher
# Run this script to start the Shiny app
#
# From RStudio: Open this file and click "Source", or use Run App on app.R
# From terminal: Rscript run.R  (run from the Meta_AML directory)

if (!file.exists("final_data_matrix.RData")) {
  stop("Data file not found. Please run from the Meta_AML directory (containing app.R and final_data_matrix.RData).")
}

message("Starting Meta_AML Explorer...")
message("The app will open in your browser. Close the browser tab to stop the app.")
message("")

shiny::runApp(".", launch.browser = TRUE)
