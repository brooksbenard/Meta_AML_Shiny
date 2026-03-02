# Meta AML Explorer

Interactive exploration of mutational, clinical outcomes, and drug sensitivity data in acute myeloid leukemia (AML). Meta AML4 merges four of the largest molecularly profiled and clinically annotated AML datasets into a single combined cohort of ~4,660 patients.

**[Launch the app on Shinyapps.io](https://brooks-benard.shinyapps.io/MetaAML_Shiny/)**

![Meta AML schematic](www/Meta_AML_schematic.png)

## Data aggregation

Clinical, mutation, and cytogenetic data from four cohorts are cleaned, harmonized, and re-annotated to common column names and units:

| Dataset | Patients | Source |
|---------|----------|--------|
| **UK-NCRI** | 2,113 | [Tazi et al. (2022), Nature Communications](https://www.nature.com/articles/s41467-022-32103-8) |
| **AML-SG** | 1,540 | [Papaemmanuil et al. (2016), NEJM](https://www.nejm.org/doi/10.1056/NEJMoa1516192) |
| **Beat AML** | 805 | [Tyner et al. (2022), Cancer Cell](https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00312-9) |
| **TCGA LAML** | 200 | [NEJM](https://www.nejm.org/doi/full/10.1056/NEJMoa1301689) |

- **Clinical data:** Cohort, demographics, AML type, ELN 2017 and ELN 2022 risk (derived where missing from cytogenetics and mutations), overall survival, key labs (WBC, hemoglobin, platelets, LDH, blasts), and simplified karyotype (Normal / Complex / Other).
- **Mutations:** MAF-style table with standard columns (Hugo_Symbol, Chromosome, positions, Variant_Classification, Variant_Type, HGVSp_Short, etc.), sample and cohort identifiers, and VAF. Variant classification and gene symbols are harmonized across sources.
- **Cytogenetics:** Arm-level gains and deletions from UK-NCRI; copy number used for CCF estimation (CCF = VAF × CN, capped at 100%).

## App functionality

- **Cohort selection** — Filter by AML type, dataset, karyotype, and ELN 2022 risk
- **Single mutation associations** — Clinical variables, survival, hazard ratios, Kaplan-Meier
- **Co-mutation analyses** — Oncoprint, odds ratios, Kaplan-Meier by 2–3 genes
- **VAF / CCF (clonality)** — VAF associations, MaxStat survival thresholds, pairwise scatterplots for clonal ordering
- **Drug sensitivity** — Mutation and VAF correlations with inhibitor AUC (Beat AML data)

## Run locally

- **RStudio:** Open `app.R` and click **Run App**.
- **Terminal:** `Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'`.

Requires `AML_Meta_Cohort_v2.rds` (or `AML_Meta_Cohort.rds` / `AML_Meta_Cohort.RData`) and `beataml2_data/` with `mutations.txt` and `inhibitor_auc.txt`.

## Deploy (Shinyapps.io)

Use the **rsconnect** configuration in this folder. The app reads: `app.R`, `load_beataml2.R`, the cohort RDS, `beataml2_data/`, and `www/`.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
