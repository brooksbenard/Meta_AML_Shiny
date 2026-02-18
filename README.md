# Meta AML Explorer (Shiny)

Interactive Shiny app for exploring acute myeloid leukemia (AML) mutational data. **This directory (Meta_AML_Shiny) is the only place to run the app**—do not run the app from the parent Meta_AML folder.

## Data sources in the app

- **Benard et al. (2021)** — Data from Benard et al. (2021) *Clonal architecture predicts clinical outcomes and drug sensitivity in acute myeloid leukemia*, Nature Communications (~2,800 patients, 12 cohorts).
- **Meta AML4** — Merged cohort of TCGA LAML (~200), BeatAML2 (805), AML-SG (1,540), and UK-NCRI (2,113) patients (~4,660 combined). Same analyses and filters.

Live app: [https://brooks-benard.shinyapps.io/MetaAML_Shiny/](https://brooks-benard.shinyapps.io/MetaAML_Shiny/).

## Running the app (only from this directory)

**From RStudio:**  
1. Open this folder (Meta_AML_Shiny) in RStudio.  
2. Open `app.R` and click **Run App**.

**From R:**  
```r
setwd("/path/to/Meta_AML_Shiny")
shiny::runApp(".", launch.browser = TRUE)
```

**From terminal:**  
```bash
cd /path/to/Meta_AML_Shiny
Rscript run.R
```

Required packages: `shiny`, `ggplot2`, `survival`, `DT`, `scales`. Optional: `survminer`, `gridExtra`, `ComplexHeatmap`, `maxstat`, `bestglm`, `readxl`.

## Setup (run once from this directory)

- **BeatAML2 (Drug Sensitivity tab):** `source("setup_beataml2.R")` — downloads into `beataml2_data/`.
- **AML-SG:** `source("setup_amlsg.R")` — downloads into `amlsg_data/`.
- **Meta AML4 cohort:** run `build_AML_Meta_Cohort_v2.R` to build `AML_Meta_Cohort_v2.rds` from TCGA, Beat AML, AML-SG, and UK-NCRI data.

Data folders used by the app: `beataml2_data/`, `amlsg_data/`, `laml_tcga_pan_can_atlas_2018/`, `UK_NCRI_data/`. The app expects `final_data_matrix.RData` and optionally `AML_Meta_Cohort_v2.rds` in this directory.

## Citation

Benard, B.A., Leak, L.B., Azizi, A. et al. Clonal architecture predicts clinical outcomes and drug sensitivity in acute myeloid leukemia. *Nat Commun* 12, 7244 (2021). https://doi.org/10.1038/s41467-021-27472-5
