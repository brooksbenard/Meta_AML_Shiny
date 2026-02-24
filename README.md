# Meta AML Explorer

A Shiny app for exploratory analysis of acute myeloid leukemia (AML) mutational data across multiple cohorts: **TCGA LAML**, **Beat AML**, **AML-SG**, and **UK-NCRI**. The app provides cohort overviews, survival and ELN risk by dataset, mutation frequency and co-occurrence, single-gene and co-mutation associations, VAF/CCF analyses, and drug sensitivity (Beat AML).

## Features

- **Cohort overview** – Survival by dataset, age and ELN risk distributions, recurrently mutated genes, oncoprint
- **Filter by cohort, subset (De novo / Secondary / Therapy), karyotype** – All analyses respect sidebar filters
- **Single-gene associations** – Clinical variables, co-mutation, VAF, survival (Kaplan–Meier, Cox), drug sensitivity
- **Co-mutation** – Odds ratios, Kaplan–Meier by 2–3 genes, UpSet-style overlap
- **VAF/CCF** – Gene-level distributions, survival by VAF/CCF, maxstat cutpoints
- **Drug sensitivity (Beat AML)** – Mut vs WT AUC, VAF vs AUC, LOOCV; optional precomputed results for faster load

## Requirements

- **R** (≥ 4.0 recommended)
- **CRAN:** `shiny`, `ggplot2`, `survival`, `DT`, `scales`
- **Optional (recommended):** `survminer`, `gridExtra`, `patchwork`, `ComplexHeatmap`, `maxstat`, `ggrepel`, `maftools`, `UpSetR` (for full plots and tables)

## Quick start

1. **Clone the repo**
   ```bash
   git clone https://github.com/brooksbenard/Meta_AML_Shiny.git
   cd Meta_AML_Shiny
   ```

2. **Get the cohort data**  
   The app expects a pre-built cohort file in the project root:
   - `AML_Meta_Cohort_v2.rds` (preferred), or  
   - `AML_Meta_Cohort.rds` / `AML_Meta_Cohort.RData`  
   These are not in the repo (see [Building the cohort](#building-the-cohort) below).

3. **Run the app**
   ```bash
   Rscript run.R
   ```
   Or in RStudio: open `app.R` and click **Run App**. The app will open in your browser.

## Building the cohort

To rebuild the meta cohort from source data, run the build script from the project root:

```bash
Rscript build_AML_Meta_Cohort_v2.R
```

This produces `AML_Meta_Cohort_v2.rds`. The script expects the following structure and files:

| Cohort    | Required inputs |
|----------|------------------|
| **TCGA** | `laml_tcga/` (cBioPortal LAML study: clinical patient/sample, mutations, CNA) |
| **Beat AML** | `beataml2_data/mutations.txt`, `beataml2_data/inhibitor_auc.txt`; optional clinical Excel |
| **AML-SG** | `amlsg_data/AMLSG_Genetic.txt`, `amlsg_data/AMLSG_Clinical_Anon.RData`, `amlsg_data/AMLSG_FLT3ITD.txt`, `amlsg_data/AMLSG_Karyotypes.txt` |
| **UK-NCRI** | `UK_NCRI_data/UK_NCRI_Mutations_data.csv`, `UK_NCRI_data/UK_NCRI_Clinical_data.csv`, `UK_NCRI_data/UK_NCRI_Cytogenetics_data.csv`; `aml_prognosis_updated.tsv` (survival + ELN) in `UK_NCRI_data/` or `UK_NCRI_data/data/`; optional `aml_molecular_bdp.tsv` for FLT3 ITD/TKD/Other |

- **UK-NCRI survival and ELN:** Place `aml_prognosis_updated.tsv` (columns: sample ID, `eln_2017` 1=Adverse/2=Intermediate/3=Favorable, `os`, `os_status`) in `UK_NCRI_data/`, `UK_NCRI_data/data/`, or `~/Downloads/data/`.
- **UK-NCRI FLT3:** If `aml_molecular_bdp.tsv` is not present, the build will try FLT3_ITD / FLT3_TKD / FLT3_other in `UK_NCRI_Clinical_data.csv` when available.

See `build_AML_Meta_Cohort_v2.R` and console output for exact paths and optional files.

## Deploying (Shinyapps.io)

See **[README_DEPLOY.md](README_DEPLOY.md)** for required files, optional precomputed drug sensitivity, and deploy steps.

## Memory and performance

- **Memory:** Typical peak ~400–500 MB per process. See **[MEMORY_PROFILE.md](MEMORY_PROFILE.md)** and run `Rscript memory_profile_app.R` to profile.
- **Drug sensitivity:** For faster startup, run `source("precompute_drug_sensitivity.R")` and include `drug_sensitivity_precomputed.rds` in your deploy or local run.

## License and citation

See repository license file. If you use this app or the meta cohort in publications, please cite the underlying data sources (TCGA, Beat AML, AML-SG, UK-NCRI) as appropriate.

## Links

- [Beat AML / BeatAML2](https://biodev.github.io/BeatAML2/)
- [TCGA LAML](https://www.cbioportal.org/study/summary?id=laml_tcga_pub)
