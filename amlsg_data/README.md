# AML-SG data (AMLSG)

Data from the [gerstung-lab/AML-multistage](https://github.com/gerstung-lab/AML-multistage) repository, accompanying *Precision oncology for acute myeloid leukemia using a knowledge bank approach*.

## Source

- **Repository:** https://github.com/gerstung-lab/AML-multistage/tree/master/data  
- **Classification:** Papaemmanuil et al. NEJM 2016

## Files

| File | Description |
|------|-------------|
| `AMLSG_Classification.txt` | Classification as in Papaemmanuil et al. NEJM 2016 |
| `AMLSG_Clinical_Anon.RData` | Anonymised clinical data table (R data dump) |
| `AMLSG_FLT3ITD.txt` | FLT3-ITD variant data for AMLSG cohort |
| `AMLSG_Genetic.txt` | Genetic data for AMLSG cohort |
| `AMLSG_Karyotypes.txt` | Karyotypic data for AMLSG cohort |

## Setup

From the `Meta_AML_Shiny` directory, run once:

```r
source("setup_amlsg.R")
```

This downloads the files above into `amlsg_data/`.
