# UK-NCRI clinical enrichment (ELN and survival)

ELN 2017 risk and overall survival (OS) for the UK NCRI cohort are filled from **aml_prognosis_updated.tsv**.

- **Filename:** `aml_prognosis_updated.tsv`
- **Location (either):** `UK_NCRI_data/aml_prognosis_updated.tsv` or `~/Downloads/data/aml_prognosis_updated.tsv`
- **Columns used:** first column = sample ID (PD#####a/c); numeric ELN (1/2/3 → Favorable/Intermediate/Adverse); `os` (years → days); `os_status` (1 = event, 0 = censored).

The build script (`build_AML_Meta_Cohort_v2.R`) uses this file only; it covers all ~2,113 NCRI clinical samples.

**AML subset (de novo vs secondary):** From the clinical CSV or prognosis file, column `secondary`: **1 = de novo** (primary AML, no antecedent hematologic disorder; ~1,748), **2 = secondary AML** (~260), **3 = therapy-related** (~105). The `ahd` (antecedent hematologic disorder) column in `data/aml_prognosis_updated.tsv` matches this: ahd=0 ↔ secondary=1 (de novo), ahd=1 ↔ secondary=2 or 3.
