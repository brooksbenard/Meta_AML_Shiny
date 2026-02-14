# Deploy Meta AML Explorer to shinyapps.io

## 1. Prerequisites

- **Data files in the app directory (same folder as `app.R`):**
  - **`final_data_matrix.RData`** — Required. Without it, the app will fail on startup.
  - **`AML_Meta_Cohort.rds`** (or **`AML_Meta_Cohort.RData`**) — Required for the **Meta AML4** tab to show the four-cohort merged data. If this file is missing when you deploy, the Meta AML4 tab will fall back to the Benard et al. (2021) dataset and the app will show a yellow warning banner. Include this file in the app directory before deploying so Meta AML4 uses the correct dataframe.
- R packages `rsconnect` (and everything the app needs: `shiny`, `ggplot2`, `survival`, `DT`, `scales`; optional: `survminer`, `gridExtra`, `ComplexHeatmap`, `maxstat`, `bestglm`).

## 2. Install rsconnect and authorize your account

In R:

```r
install.packages("rsconnect")
library(rsconnect)
```

**First-time setup:**

1. Create an account at **[shinyapps.io](https://www.shinyapps.io/)** (free tier available).
2. In the dashboard, go to **Account** → **Tokens** (or click your name → Tokens).
3. Click **Show** next to your token and copy the pre-generated `rsconnect::setAccountInfo(...)` command.
4. Paste and run that command in the R console once. It looks like:
   ```r
   rsconnect::setAccountInfo(
     name   = 'YOUR_ACCOUNT',
     token  = 'YOUR_TOKEN',
     secret = 'YOUR_SECRET'
   )
   ```
   After this, you don’t need to do it again on this machine.

**Alternative (RStudio):** **Tools** → **Global Options** → **Publishing** → connect your shinyapps.io account.

## 3. Deploy the app

From R, set the working directory to the folder that contains `app.R` and `final_data_matrix.RData`, then deploy:

```r
# Use the path to your app folder (e.g. Meta_AML_Shiny or Meta_AML)
setwd("/path/to/Meta_AML_Shiny")   # or wherever app.R and final_data_matrix.RData live

library(rsconnect)
deployApp(
  appDir = getwd(),
  appName = "meta-aml-explorer",   # URL will be youraccount.shinyapps.io/meta-aml-explorer
  appTitle = "Meta AML Explorer",
  launch.browser = TRUE
)
```

- **appName** must be unique in your account (letters, numbers, `_`, `-`). Change it if you get a conflict.
- **appTitle** is the label shown in your shinyapps.io dashboard.

Deployment will bundle everything in `appDir` (including `final_data_matrix.RData`, `AML_Meta_Cohort.rds` if present, `beataml2_data/`, `www/`, and the R scripts). The first deploy can take a few minutes if the data files are large.

## 4. Update the app after changes

From the same directory, run the same `deployApp(...)` call again. rsconnect will upload only changed files and redeploy.

## 5. Notes

- **Size limits:** Free tier has size and usage limits. If the bundle is too large, you may need a paid plan or to reduce data (e.g. a subset of `final_data_matrix`).
- **Startup time:** Loading `final_data_matrix.RData` and preparing data can make the first load slow; that’s normal.
- **Secrets:** Never commit your token/secret to git. The `setAccountInfo()` command is for your local R session only.
