# Deploy-only copy

This directory is trimmed for shinyapps.io deployment (~24 MB without .git).

**Before deploy (reduces memory / avoids OOM):** From this directory run:
`source("precompute_drug_sensitivity.R")` to create `drug_sensitivity_precomputed.rds`. The app will use it for the Drug Sensitivity tab (Mut vs WT, VAF vs AUC correlations, LOOCV) instead of computing on the fly. Include the generated `.rds` in your deploy bundle.

**To restore build scripts and data** (e.g. to rebuild the app from scratch):

- From GitHub: clone or pull from https://github.com/brooksbenard/Meta_AML_Shiny  
- Locally: `git checkout -- .` or `git restore .` to restore deleted files from your repo.

**Required for deploy (kept):** `app.R`, `load_beataml2.R`, `final_data_matrix.RData`, `AML_Meta_Cohort_v2.rds` (or `AML_Meta_Cohort.rds`), `beataml2_data/`, `www/`. Optional but recommended: `drug_sensitivity_precomputed.rds` (run `precompute_drug_sensitivity.R` to generate).

**Excluded from upload:** `.git` is listed in `.rscignore` so it is not sent to shinyapps.io.
