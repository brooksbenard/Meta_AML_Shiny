# Meta AML Explorer (Shiny)

Interactive Shiny app for exploring acute myeloid leukemia (AML) mutational data. Two data sources:

- **Benard et al. (2021)** — Data from Benard et al. (2021) *Clonal architecture predicts clinical outcomes and drug sensitivity in acute myeloid leukemia*, Nature Communications (~2,800 patients, 12 cohorts).
- **Meta AML4** — Merged cohort of TCGA LAML (~200), BeatAML2 (805), AML-SG (1,540), and UK-NCRI (2,113) patients (~4,660 combined). Same analyses and filters.

## Quick start

### 1. Data files (required)

Place the main dataset in this directory:

- **`final_data_matrix.RData`** — Required. Contains the `final_data_matrix` object from the 2021 paper (mutation calls, VAF, clinical, survival). Not included in this repo due to size. Obtain from the [paper repository](https://github.com/brooksbenard/Meta_AML) or run the data aggregation pipeline there.

Optional (for the Meta AML4 tab):

- **`AML_Meta_Cohort.RData`** or **`AML_Meta_Cohort.rds`** — Merged four-cohort data. If missing, the Meta AML4 tab falls back to the Benard et al. (2021) dataset.

### 2. Run the app

**From RStudio:** Open `app.R` and click **Run App**.

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

### 3. R packages

```r
install.packages(c("shiny", "ggplot2", "survival", "DT", "scales"))
# Optional: survminer, gridExtra, ComplexHeatmap, maxstat, bestglm (for extra plots)
# For Drug Sensitivity subset filter: install.packages("readxl")
```

### 4. Drug Sensitivity tab (BeatAML2)

To enable the Drug Sensitivity tab, download BeatAML2 data once:

```r
setwd("/path/to/Meta_AML_Shiny")
source("setup_beataml2.R")
```

This populates `beataml2_data/`. If the folder is empty, the Drug Sensitivity tab will show a message to run the setup.

### 5. Banner image (optional)

For the top-banner profile image, place your image in the `www/` folder as:

`www/linkedin_pic%20copy.jpeg`

(or update the image filename in `app.R` to match your file). The app runs without it.

## Repository layout

```
Meta_AML_Shiny/
├── app.R              # Shiny app
├── load_beataml2.R    # BeatAML2 loading for Drug Sensitivity
├── setup_beataml2.R   # Download BeatAML2 data (run once)
├── run.R              # Launcher script
├── beataml2_data/     # BeatAML2 files (after setup_beataml2.R)
├── www/               # Static assets (e.g. banner image)
├── final_data_matrix.RData   # You add this
├── AML_Meta_Cohort.rds       # Optional, for Meta AML4 tab
└── README.md
```

## Citation

Benard, B.A., Leak, L.B., Azizi, A. et al. Clonal architecture predicts clinical outcomes and drug sensitivity in acute myeloid leukemia. *Nat Commun* 12, 7244 (2021). https://doi.org/10.1038/s41467-021-27472-5
