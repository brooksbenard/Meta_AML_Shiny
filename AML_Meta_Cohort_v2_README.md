# AML_Meta_Cohort_v2

Version 2 of the merged AML cohort data file, built directly from the four source datasets:
- **TCGA LAML** (~200 patients)
- **Beat AML** (805 patients, 942 specimens)
- **AML-SG** (1,540 patients)
- **UK-NCRI** (2,113 clinical; 1,950 with mutations in merged file)

**Total: 4,456 unique samples, 29,398 mutation records**

### Dataset makeup (after build)

| Cohort   | Unique samples | Mutation rows | Subset (main)        | Risk (main)                    | Survival |
|----------|----------------|---------------|-----------------------|--------------------------------|----------|
| TCGA     | 199            | 10,847        | de_novo               | Unknown (no ELN in source)     | Yes (from OS_MONTHS) |
| Beat_AML | 895            | 8,132         | de_novo, relapse, secondary, therapy, other | Favorable/Intermediate/Adverse/Unknown | Yes |
| AML-SG   | 1,412          | 3,908         | de_novo, secondary, therapy, other | Favorable/Intermediate/Adverse/Unknown | Yes |
| UK-NCRI  | 1,950          | 6,511         | de_novo, secondary   | Favorable/Intermediate/Adverse (from aml_prognosis_updated.tsv) | Yes |

- **Column names** in `AML_Meta_Cohort_v2.rds`: `Sample`, `Gene`, `VAF`, `Cohort`, `Subset`, `Time_to_OS`, `Censor`, `variant_type`, `mutation_category`, `Age`, `Sex`, `Risk` (all required by the Meta AML4 app).
- UK-NCRI ELN and survival come from **aml_prognosis_updated.tsv** only (full NCRI cohort).

## Build Script

Run `build_AML_Meta_Cohort_v2.R` to rebuild the dataset:

```r
source("build_AML_Meta_Cohort_v2.R")
```

## Data Sources

### TCGA (`laml_tcga_pan_can_atlas_2018/`)
- Mutations: `data_mutations.txt`
- Clinical: `data_clinical_patient.txt`, `data_clinical_sample.txt`
- Survival: Converted from months to days (OS_MONTHS × 30.44)
- Risk: Not available in TCGA clinical files (set to "Unknown")

### Beat AML (`beataml2_data/`)
- Mutations: `mutations.txt` + consensus mutations from clinical file
- Clinical: `clinical.xlsx` (or `beataml_wv1to4_clinical.xlsx`)
- **Consensus mutations added**: FLT3-ITD, NPM1, RUNX1, ASXL1, TP53 from clinical file (per BeatAML2 FAQ)
- Survival: Already in days
- Subset: de_novo, secondary, therapy, relapse (from isDenovo, isTransformed, specificDxAtAcquisition, diseaseStageAtSpecimenCollection)
- Risk: ELN2017 (Favorable/Intermediate/Adverse)

### AML-SG (`amlsg_data/`)
- Mutations: `AMLSG_Genetic.txt`
- Clinical: `AMLSG_Clinical_Anon.RData` (object: `clinicalData`)
- Survival: Already in days (OS column)
- Subset: de_novo (AML), secondary (sAML), therapy (tAML)
- Risk: M_Risk mapped from "Favorable", "Inter-1"/"Inter-2"→"Intermediate", "Adverse"

### UK-NCRI (`UK_NCRI_data/`)
- Mutations: `UK_NCRI_Mutations_data.csv`
- Clinical: `UK_NCRI_Clinical_data.csv` + **aml_prognosis_updated.tsv** (ELN and survival for all NCRI samples)
- Survival: From aml_prognosis_updated.tsv (`os` years → days, `os_status` → Censor)
- Subset: de_novo (secondary==0), secondary (secondary==1) from clinical CSV
- Risk: From aml_prognosis_updated.tsv (ELN 2017: 1/2/3 → Favorable/Intermediate/Adverse)

## Data Normalization

### Survival Times
- **Units**: All converted to **days**
- TCGA: Months × 30.44
- Beat AML: Already in days
- AML-SG: Already in days
- UK-NCRI: Not available

### ELN Risk Categories
- Standardized to: **Favorable**, **Intermediate**, **Adverse**, **Unknown**
- Beat AML: Direct from ELN2017 column
- AML-SG: Mapped from M_Risk (Inter-1/Inter-2 → Intermediate)
- TCGA: Not available (Unknown)
- UK-NCRI: Not available (Unknown)

### Subset Labels
- Standardized to: **de_novo**, **secondary**, **therapy**, **relapse**, **other**
- Beat AML: From isDenovo, isTransformed, specificDxAtAcquisition, diseaseStageAtSpecimenCollection
- AML-SG: From TypeAML (AML→de_novo, sAML→secondary, tAML→therapy)
- TCGA: All de_novo (primary AML)
- UK-NCRI: From secondary flag (0→de_novo, 1→secondary)

### Mutation Categories
- Assigned from `final_data_matrix.RData` reference (Benard et al. 2021)
- FLT3-ITD and FLT3-TKD: Assigned "RTK-RAS Signaling"

## Validation

### Mutation Frequencies
Top genes across cohorts (samples with mutation):
- NPM1: 1,410 total (TCGA: 54, Beat AML: 242, AML-SG: 440, UK-NCRI: 674)
- DNMT3A: 1,195 total
- FLT3-ITD: 219 (Beat AML only, from consensus calls)
- FLT3: 856 total (excluding ITD)

### Data Completeness
- Samples with survival data: 11,965 / 4,456 unique samples
- Samples with age data: 18,359 / 4,456 unique samples
- Samples with sex data: 10,349 / 4,456 unique samples

## File Format

Saved as `AML_Meta_Cohort_v2.rds` with columns:
- `Sample`: Sample/patient identifier
- `Gene`: Gene symbol
- `VAF`: Variant allele frequency (0-100%)
- `Cohort`: TCGA, Beat_AML, AML-SG, UK-NCRI
- `Subset`: de_novo, secondary, therapy, relapse, other
- `Time_to_OS`: Overall survival time (days)
- `Censor`: Censoring status (1=event, 0=censored)
- `variant_type`: SNV, Deletion, Splicing, Indel, ITD, PTD, Other
- `mutation_category`: DNA Methylation, Chromatin/Cohesin, RTK-RAS Signaling, Splicing, Transcription, NPM1, Tumor suppressors
- `Age`: Age at diagnosis
- `Sex`: Male, Female
- `Risk`: Favorable, Intermediate, Adverse, Unknown

## Usage in Shiny App

The app (`app.R`) automatically loads `AML_Meta_Cohort_v2.rds` for the **Meta AML4** tab if it exists, falling back to `AML_Meta_Cohort.rds` or `final_data_matrix.RData` if v2 is not found.

The **Benard et al. (2021)** tab continues to use `final_data_matrix.RData` as before.
