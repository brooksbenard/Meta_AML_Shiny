# Shinyapps.io deploy – required files

This directory is trimmed so only files needed for deployment are kept.

## Required for deploy (do not remove)

| File or folder | Purpose |
|----------------|--------|
| `app.R` | Main Shiny app |
| `load_beataml2.R` | Sourced by app; loads BeatAML2 mutations/AUC |
| `AML_Meta_Cohort_v2.rds` | Meta AML4 cohort (or `AML_Meta_Cohort.rds` / `AML_Meta_Cohort.RData`) |
| `beataml2_data/` | Must contain: `mutations.txt`, `inhibitor_auc.txt`; optional: `beataml_wv1to4_clinical.xlsx` or `clinical.xlsx` |
| `www/` | Static assets (e.g. `linkedin_pic copy.jpeg` used in UI) |
| `.rscignore` | Excludes `.git` etc. from upload |
| `rsconnect/` | Deploy target config (shinyapps.io) |

## Optional but recommended

- **`drug_sensitivity_precomputed.rds`** – Run `source("precompute_drug_sensitivity.R")` in this directory to create it. The app uses it for the Drug Sensitivity tab instead of computing on the fly (reduces memory and startup time). Include the generated file in your deploy bundle.
- **`precompute_drug_sensitivity.R`** – Script to regenerate the file above (can be run locally only).

## Not needed for deploy (removed)

- `build_AML_Meta_Cohort_v2.R` – Build script; only for rebuilding the cohort locally.
- `amlsg_data/`, `UK_NCRI_data/`, `laml_tcga/` – Raw inputs for the build script; not read by the app at runtime.
- `final_data_matrix.RData` – App uses in-memory stub and `AML_Meta_Cohort_v2.rds` for Meta AML4.
- `beataml2_data/sample_mapping.xlsx` – Not used by `load_beataml2.R`.

## Restoring build/data

To rebuild the cohort or run the build script again, restore from GitHub:

```bash
git clone https://github.com/brooksbenard/Meta_AML_Shiny .
# or: git checkout -- build_AML_Meta_Cohort_v2.R amlsg_data UK_NCRI_data laml_tcga
```

- **UK-NCRI survival and ELN risk:** The build script expects `aml_prognosis_updated.tsv` (with sample IDs, ELN 1/2/3, OS in years, and os_status) in `UK_NCRI_data/`, `UK_NCRI_data/data/`, or `~/Downloads/data/`. If this file is missing, UK-NCRI will have no survival or ELN risk in the cohort; the script will emit a warning.

## Run locally

From this directory: `Rscript run.R` or in RStudio open `app.R` and click Run App. `run.R` checks for `AML_Meta_Cohort_v2.rds` (or fallback cohort file) before starting.
