# =============================================================================
# Meta AML Explorer - Shiny App for Exploratory Analysis of AML Mutational Data
# Based on Benard et al. Nature Communications 2021
# Uses base R only to avoid data.table/dplyr "Argument 'x' is not a vector: list" errors
# =============================================================================

library(shiny)
library(ggplot2)
library(survival)
library(DT)
library(scales)

has_survminer <- requireNamespace("survminer", quietly = TRUE)
if (has_survminer) library(survminer)
has_gridExtra <- requireNamespace("gridExtra", quietly = TRUE)
if (has_gridExtra) library(gridExtra)
has_ComplexHeatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)
if (has_ComplexHeatmap) library(ComplexHeatmap)
has_maxstat <- requireNamespace("maxstat", quietly = TRUE)
if (has_maxstat) library(maxstat)
has_bestglm <- requireNamespace("bestglm", quietly = TRUE)
if (has_bestglm) library(bestglm)

# Load data (app expects final_data_matrix.RData in the same directory as app.R)
data_path <- file.path(getwd(), "final_data_matrix.RData")
if (!file.exists(data_path)) {
  stop("Data file not found. Place final_data_matrix.RData in the app directory and run the app from that directory.\n  Looked for: ", normalizePath(data_path, mustWork = FALSE))
}
load(data_path)

# Convert to plain data.frame - CSV round-trip strips data.table/list columns
tmp <- tempfile(fileext = ".csv")
utils::write.csv(final_data_matrix, tmp, row.names = FALSE)
final_data_matrix <- utils::read.csv(tmp, stringsAsFactors = FALSE)
unlink(tmp)

# Normalize AML_Meta_Cohort to match final_data_matrix column names so existing analyses work
normalize_AML_Meta_Cohort <- function(d) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  # Map source columns -> app column names (select + rename)
  # Sample: patient_id or sample_id
  if ("patient_id" %in% colnames(d)) d$Sample <- as.character(d$patient_id)
  else if ("sample_id" %in% colnames(d)) d$Sample <- as.character(d$sample_id)
  # Gene, VAF, Cohort
  if (!"Sample" %in% colnames(d)) d$Sample <- paste0("S", seq_len(nrow(d)))
  if ("Gene" %in% colnames(d)) d$Gene <- as.character(d$Gene)
  if ("VAF" %in% colnames(d)) d$VAF <- as.numeric(d$VAF)
  if ("Cohort" %in% colnames(d)) {
    d$Cohort <- as.character(d$Cohort)
    d$Cohort[d$Cohort == "NCRI"] <- "UK-NCRI"
    d$Cohort[d$Cohort == "AMLSG"] <- "AML-SG"
  }
  # Subset: type_AML or De_novo / Secondary / Relapse
  if ("type_AML" %in% colnames(d)) {
    t <- as.character(d$type_AML)
    d$Subset <- ifelse(t %in% c("De_novo", "de_novo"), "de_novo",
      ifelse(t %in% c("Secondary", "secondary"), "secondary",
      ifelse(t %in% c("Relapse", "relapse"), "relapse",
      ifelse(t %in% c("tAML", "therapy"), "therapy",
      ifelse(t %in% c("oAML", "other", "unknown", ""), "other", "other")))))
  } else if (all(c("De_novo", "Secondary") %in% colnames(d))) {
    d$Subset <- "other"
    de_ok <- which(d$De_novo %in% c(TRUE, "TRUE", "True", 1))
    sec_ok <- which(d$Secondary %in% c(TRUE, "TRUE", "True", 1))
    rel_ok <- if ("Relapse" %in% colnames(d)) which(d$Relapse %in% c(TRUE, "TRUE", "True", 1)) else integer(0)
    if (length(de_ok)) d$Subset[de_ok] <- "de_novo"
    if (length(sec_ok)) d$Subset[sec_ok] <- "secondary"
    if (length(rel_ok)) d$Subset[rel_ok] <- "relapse"
  } else d$Subset <- "other"
  # Survival: os -> Time_to_OS, os_status -> Censor (1 = event, 0 = censored)
  if ("os" %in% colnames(d)) { d$Time_to_OS <- as.numeric(d$os); d$Time_to_OS[d$Time_to_OS < 0] <- NA }
  if ("os_status" %in% colnames(d)) d$Censor <- as.numeric(d$os_status)
  # variant_classification -> variant_type (Meta AML4 scheme: SNV, Deletion, Splicing, Indel, ITD, PTD, Other)
  if ("variant_classification" %in% colnames(d)) {
    vc <- as.character(d$variant_classification)
    d$variant_type <- "Other"
    d$variant_type[vc %in% c("Sub", "SNV", "SNP", "missense_variant", "stop_gained")] <- "SNV"
    d$variant_type[vc %in% c("Del", "DEL")] <- "Deletion"
    d$variant_type[vc %in% c("splice_acceptor_variant", "splice_donor_variant")] <- "Splicing"
    d$variant_type[vc %in% c("Indel")] <- "Indel"
    d$variant_type[vc %in% c("ITD")] <- "ITD"
    d$variant_type[vc %in% c("PTD")] <- "PTD"
  } else d$variant_type <- "Other"
  # mutation_category, Age, Sex, Risk (Sex: MALE->Male, FEMALE->Female; ELN2017: Adverse/Intermediate/Favorable else Unknown)
  if ("mutation_category" %in% colnames(d)) d$mutation_category <- as.character(d$mutation_category)
  if ("Age" %in% colnames(d)) d$Age <- as.numeric(d$Age)
  if ("Sex" %in% colnames(d)) {
    sx <- as.character(d$Sex)
    d$Sex <- ifelse(sx == "MALE", "Male", ifelse(sx == "FEMALE", "Female", sx))
  } else if ("sex" %in% colnames(d)) {
    sx <- as.character(d$sex)
    d$Sex <- ifelse(sx == "MALE", "Male", ifelse(sx == "FEMALE", "Female", sx))
  }
  if ("ELN2017" %in% colnames(d)) {
    eln <- as.character(d$ELN2017)
    d$Risk <- ifelse(eln %in% c("Adverse", "Intermediate", "Favorable"), eln, "Unknown")
  } else if ("Karyotype" %in% colnames(d)) d$Risk <- as.character(d$Karyotype)
  else d$Risk <- NA_character_
  # Clinical vars: WBC, Hemoglobin, Platelet, LDH, BM_blast_percent, PB_blast_percent
  if ("WBC" %in% colnames(d)) d$WBC <- as.numeric(d$WBC)
  if ("hemoglobin" %in% colnames(d)) d$Hemoglobin <- as.numeric(d$hemoglobin)
  if ("Platelet" %in% colnames(d)) d$Platelet <- as.numeric(d$Platelet)
  if ("LDH" %in% colnames(d)) d$LDH <- as.numeric(d$LDH)
  if ("bm_blasts" %in% colnames(d)) d$BM_blast_percent <- as.numeric(d$bm_blasts)
  if ("pb_blasts" %in% colnames(d)) d$PB_blast_percent <- as.numeric(d$pb_blasts)
  # Keep only columns the app uses (drop extras to avoid confusion)
  want <- c("Sample", "Gene", "VAF", "Cohort", "Subset", "Time_to_OS", "Censor", "variant_type",
    "mutation_category", "Age", "Sex", "Risk", "WBC", "Hemoglobin", "Platelet", "LDH", "BM_blast_percent", "PB_blast_percent")
  keep <- want[want %in% colnames(d)]
  d[, keep, drop = FALSE]
}

# Load AML_Meta_Cohort (four cohorts merged) for Meta AML4 tab; fallback to final_data_matrix if missing
AML_Meta_Cohort <- NULL
meta4_uses_dedicated_file <- FALSE
meta4_path_rdata <- file.path(getwd(), "AML_Meta_Cohort.RData")
meta4_path_rds <- file.path(getwd(), "AML_Meta_Cohort.rds")
if (file.exists(meta4_path_rdata)) {
  env_meta4 <- new.env()
  load(meta4_path_rdata, envir = env_meta4)
  nms <- ls(env_meta4)
  if ("AML_Meta_Cohort" %in% nms) AML_Meta_Cohort <- get("AML_Meta_Cohort", env_meta4)
  else if (length(nms) == 1L) AML_Meta_Cohort <- get(nms[1], env_meta4)
  if (!is.null(AML_Meta_Cohort)) meta4_uses_dedicated_file <- TRUE
} else if (file.exists(meta4_path_rds)) {
  AML_Meta_Cohort <- as.data.frame(readRDS(meta4_path_rds), stringsAsFactors = FALSE)
  meta4_uses_dedicated_file <- TRUE
}
if (!is.null(AML_Meta_Cohort) && is.data.frame(AML_Meta_Cohort)) {
  AML_Meta_Cohort <- normalize_AML_Meta_Cohort(AML_Meta_Cohort)
  # Classify genes in Meta AML4 using same mutation_category as Benard et al. (2021) (final_data_matrix)
  ref <- final_data_matrix[!is.na(final_data_matrix$mutation_category) & as.character(final_data_matrix$mutation_category) != "", c("Gene", "mutation_category"), drop = FALSE]
  if (nrow(ref) > 0) {
    ref$Gene <- as.character(ref$Gene)
    ref$mutation_category <- as.character(ref$mutation_category)
    most_common <- aggregate(mutation_category ~ Gene, data = ref, FUN = function(x) names(sort(table(x), decreasing = TRUE))[1])
    gene_to_cat <- setNames(as.character(most_common$mutation_category), as.character(most_common$Gene))
    AML_Meta_Cohort$mutation_category <- unname(gene_to_cat[as.character(AML_Meta_Cohort$Gene)])
  }
} else {
  AML_Meta_Cohort <- final_data_matrix
}
tmp2 <- tempfile(fileext = ".csv")
utils::write.csv(AML_Meta_Cohort, tmp2, row.names = FALSE)
AML_Meta_Cohort <- utils::read.csv(tmp2, stringsAsFactors = FALSE)
unlink(tmp2)

# Paper color palettes (Benard et al. Nat Commun 2021)
PAL_COHORT <- c(
  Tyner = "#0073C2FF", TCGA = "#EFC000FF", Majeti = "#868686FF", Papaemmanuil = "#CD534CFF",
  Lindsley = "#7AA6DCFF", Wang = "#E64B35FF", Au = "#4DBBD5FF", Welch = "#00A087FF",
  Garg = "#3C5488FF", Greif = "#F39B7FFF", Hirsch = "#7E6148FF", Huet = "#B09C85FF"
)
PAL_SUBSET <- c(
  de_novo = "#C16622FF", secondary = "#767676FF", relapse = "#800000FF",
  other = "#FFA319FF", therapy = "#8F3931FF"
)
PAL_MUT_CAT <- c(
  "DNA Methylation" = "#374E55FF", "Chromatin/Cohesin" = "#DF8F44FF",
  "RTK/RAS Signaling" = "#00A1D5FF", "Splicing" = "#B24745FF",
  "Transcription" = "#79AF97FF", "NPM1" = "#80796BFF", "Tumor suppressors" = "#6A6599FF"
)
PAL_RISK <- c(Adverse = "#E64B35FF", Intermediate = "#8491B4FF", Favorable = "#00A087FF", Unknown = "#767676FF")
# Paper Figure 1B oncoprint colors (exact from manuscript code)
ONCO_COHORT_COL <- c(Tyner = "#0073C2FF", TCGA = "#EFC000FF", Majeti = "#868686FF", Papaemmanuil = "#CD534CFF",
  Lindsley = "#7AA6DCFF", Wang = "#E64B35FF", Au = "#4DBBD5FF", Welch = "#00A087FF", Garg = "#3C5488FF",
  Greif = "#F39B7FFF", Li = "#8491B4FF", Shlush = "#91D1C2FF", Parkin = "#DC0000FF", Hirsch = "#7E6148FF", Huet = "#B09C85FF")
ONCO_SUBSET_COL <- c(de_novo = "#C16622FF", secondary = "#767676FF", relapse = "#800000FF", other = "#FFA319FF",
  Remission = "#8A9045FF", Residual = "#155F83FF", therapy = "#8F3931FF", MDS = "#58593FFF", Unknown = "gray90")
ONCO_VAR_COL <- c(Deletion = "#374E55FF", INDEL = "#DF8F44FF", Insertion = "#00A1D5FF", ITD = "#79AF97FF",
  SNV = "#B24745FF", Splicing = "#6A6599FF", Unknown = "#80796BFF")
# Meta AML4 tab: uniform cohort/subset/variant colors (four cohorts: Beat_AML, TCGA, AML-SG, UK-NCRI)
META4_COHORT_COL <- c(Beat_AML = "#155F83FF", TCGA = "#EFC000FF", `AML-SG` = "#CD534CFF", `UK-NCRI` = "#00A087FF")
META4_SUBSET_COL <- c(De_novo = "#C16622FF", Secondary = "#767676FF", Relapse = "#800000FF",
  oAML = "#FFA319FF", unknown = "#8A9045FF", tAML = "#8F3931FF",
  de_novo = "#C16622FF", secondary = "#767676FF", relapse = "#800000FF", therapy = "#8F3931FF", other = "#FFA319FF")
META4_VAR_COL <- c(Deletion = "#374E55FF", Indel = "#DF8F44FF", Insertion = "#00A1D5FF", ITD = "#79AF97FF",
  SNV = "#B24745FF", Splicing = "#6A6599FF", PTD = "tan", Other = "#80796BFF", Unknown = "#80796BFF")

# Load BeatAML2 helpers (only if data exists)
beataml2_script <- file.path(getwd(), "load_beataml2.R")
if (file.exists(beataml2_script)) source(beataml2_script, local = FALSE)

# FLT3 annotation
prepare_data <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$Gene <- as.character(df$Gene)
  df$variant_type <- as.character(df$variant_type)
  flt3_itd <- df$Gene == "FLT3" & df$variant_type %in% c("ITD", "INDEL", "Indel")
  flt3_tkd <- df$Gene == "FLT3" & df$variant_type %in% c("SNV", "Deletion", "other", "Other", "Unknown")
  df$Gene[flt3_itd] <- "FLT3-ITD"
  df$Gene[flt3_tkd] <- "FLT3-TKD"
  df
}

# UI
ui <- fluidPage(
  title = "Meta AML Explorer",
  tags$head(tags$style(HTML(
    ".nav-tabs{font-weight:600}.well{background:#f8f9fa;border-radius:8px}
     .app-banner{display:flex;align-items:center;justify-content:space-between;padding:12px 20px;background:#374e55;color:#fff;width:100vw;position:relative;left:50%;right:50%;margin-left:-50vw;margin-right:-50vw;margin-bottom:15px;}
     .app-banner .banner-title{font-size:22px;font-weight:700;margin:0;}
     .app-banner .banner-contact{display:flex;align-items:center;gap:12px;}
     .app-banner .banner-contact img{border-radius:50%;width:42px;height:42px;object-fit:cover;}
     .app-banner .banner-contact a{color:#fff;text-decoration:none;font-size:14px;}
     .app-banner .banner-contact a:hover{color:#df8f44;text-decoration:underline;}
     .welcome-page h2{color:#374e55;margin-top:1.2em;}
     .welcome-page h2:first-of-type{margin-top:0;}
     .welcome-page p{font-size:15px;line-height:1.6;}
     .welcome-page ul{font-size:15px;line-height:1.7;}
     #main_nav .nav-tabs{margin-bottom:0;border-bottom:2px solid #374e55;}
     #main_nav .nav-tabs .nav-link{font-weight:600;color:#374e55;}
     #main_nav .nav-tabs .nav-link.active{background:#374e55;color:#fff;border-color:#374e55;}
     #drug_subset_wrap .shiny-input-container{width:140px !important;min-width:140px;}"
  ))),
  div(id = "data-loading-overlay", style = "display: none; position: fixed; top: 0; left: 0; right: 0; bottom: 0; background: rgba(255,255,255,0.92); z-index: 9999; align-items: center; justify-content: center;",
    tags$div(style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); text-align: center; padding: 2em;",
      tags$div(style = "font-size: 18px; color: #374e55; margin-bottom: 12px;", strong("Loading data...")),
      tags$p(style = "color: #666; font-size: 14px; margin: 0;", "Preparing analyses. This may take a moment.")
    )
  ),
  tags$script(HTML("
    (function() {
      var analysesLoaded = false, meta4Loaded = false;
      function showOverlay() { $('#data-loading-overlay').css({'display': 'flex'}); }
      function hideOverlay() { $('#data-loading-overlay').hide(); }
      $(document).on('click', 'a[data-value=\"analyses\"]', function() {
        if (!analysesLoaded) showOverlay();
      });
      $(document).on('click', 'a[data-value=\"meta_aml4\"]', function() {
        if (!meta4Loaded) showOverlay();
      });
      Shiny.addCustomMessageHandler('dataReady', function(tab) {
        if (tab === 'analyses') analysesLoaded = true;
        if (tab === 'meta_aml4') meta4Loaded = true;
        hideOverlay();
      });
    })();
  ")),
  div(class = "app-banner",
    div(class = "banner-title", "Meta AML Explorer"),
    div(class = "banner-contact",
      tags$a(href = "mailto:bbenard@stanford.edu", "bbenard@stanford.edu"),
      tags$img(src = "linkedin_pic%20copy.jpeg", alt = "Contact", title = "Contact")
    )
  ),
  tabsetPanel(
    id = "main_nav",
    tabPanel("About", value = "about",
      div(class = "welcome-page", style = "max-width: 800px; margin: 0 auto; padding: 24px 15px;",
        h2("Welcome to Meta AML Explorer"),
        p("This site provides interactive exploration of acute myeloid leukemia (AML) mutational and clinical outcomes data. The ", strong("Benard et al. (2021)"), " and ", strong("Meta AML4"), " tabs offer different data sources and the same set of analyses (cohort summaries, survival, VAF, co-mutations, drug sensitivity)."),
        h2("Benard et al. (2021) tab"),
        p("Data from ", strong("Benard et al. (2021)"), " ", em("Clonal architecture predicts clinical outcomes and drug sensitivity in acute myeloid leukemia."), " Nature Communications | ", tags$a("Paper", href = "https://www.nature.com/articles/s41467-021-27472-5", target = "_blank")),
        p("The paper aggregates mutation calls, variant allele frequencies (VAF), and clinical outcomes from ~2,800 patients across 12 cohorts (TCGA, Beat AML/Tyner, Papaemmanuil, Majeti, and others) to link clonal architecture to prognosis and therapy response."),
        h2("Meta AML4 tab"),
        p("Meta AML4 merges ", strong("four of the largest molecularly profiled and clinically annotated AML datasets"), " into a single cohort: TCGA LAML (~200 patients), BeatAML2 (805 patients), AML-SG (1,540 patients), and UK-NCRI (2,113 patients)—", strong("~4,660 patients combined"), ". The same analyses (Overview, Single-Gene Associations, Co-mutation, VAF, Drug Sensitivity, Data Table) are available with cohort- and subset-specific filters."),
        p("Data sources: AML-SG and UK-NCRI from ", tags$a("Tazi et al. (2022)", href = "https://www.nature.com/articles/s41467-022-32103-8", target = "_blank"), " (Unified classification and risk-stratification in AML, Nature Communications); BeatAML2 from ", tags$a("biodev.github.io/BeatAML2", href = "https://biodev.github.io/BeatAML2/", target = "_blank"), "; TCGA LAML from ", tags$a("cBioPortal", href = "https://www.cbioportal.org/", target = "_blank"), "."),
        h2("What you can do in both tabs"),
        tags$ul(
          tags$li(strong("Overview"), "— Cohort summary, sample counts by cohort and AML subset, and an OncoPrint of mutations and clinical annotations."),
          tags$li(strong("Single-Gene Associations"), "— Clinical variables by mutation, Kaplan-Meier survival, hazard ratios, and survival summary tables."),
          tags$li(strong("Co-mutation Associations"), "— Mutation co-occurrence heatmap, top co-occurring pairs, and Kaplan-Meier by co-mutation status (2 or 3 genes)."),
          tags$li(strong("VAF Associations"), "— VAF by gene, VAF and survival (MaxStat), VAF by mutation category and cohort, and pairwise VAF scatter for clonal ordering."),
          tags$li(strong("Drug Sensitivity"), "— Correlations of VAF and mutation status with inhibitor sensitivity (AUC) from BeatAML2, including LOOCV and top correlations."),
          tags$li(strong("Data Table"), "— Filterable, searchable mutation table.")
        ),
        p("Use the ", strong("Filters"), " sidebar to restrict by AML subset, cohort, and minimum mutation frequency."),
        h2("Caveats"),
        p("This site is intended for research purposes only and is in active development. Some initial data loadings may take a few seconds."),
        p("The current implementation does not account for copy number correction of variant allele frequencies (VAFs). A future implementation will incorporate copy number data to define and use cancer cell fractions of mutations instead of VAFs."),
        h2("About the author"),
        p("Brooks Benard, Stanford University. ", tags$a("Website", href = "https://brooksbenard.github.io/", target = "_blank"), " · ", tags$a("GitHub", href = "https://github.com/brooksbenard", target = "_blank"), " · ", tags$a(href = "mailto:bbenard@stanford.edu", "bbenard@stanford.edu")),
        hr(),
        p(style = "font-size: 0.9em; color: #666;",
          "© ", format(Sys.Date(), "%Y"), " Brooks Benard. Licensed under ",
          tags$a("MIT License", href = "https://opensource.org/licenses/MIT", target = "_blank"), "."
        )
      )
    ),
    tabPanel("Benard et al. (2021)", value = "analyses"),
    tabPanel("Meta AML4", value = "meta_aml4"),
    conditionalPanel(
      condition = "input.main_nav === 'analyses' || input.main_nav === 'meta_aml4'",
      sidebarLayout(
        sidebarPanel(
          width = 2,
          h4("Filters"),
          selectInput("subset", "AML Subset:", 
            choices = c("All", "de_novo", "secondary", "relapse", "therapy", "other"), selected = "de_novo"),
          selectInput("cohort", "Cohort:", 
            choices = c("All", sort(unique(as.character(na.omit(final_data_matrix$Cohort))))), selected = "All"),
          sliderInput("min_freq", "Min. mutation frequency:", min = 1, max = 100, value = 50, step = 5),
          hr(),
          p(a("Paper", href = "https://www.nature.com/articles/s41467-021-27472-5", target = "_blank"), "|",
            a("GitHub", href = "https://github.com/brooksbenard/Meta_AML", target = "_blank"), class = "text-muted small")
        ),
        mainPanel(
          width = 10,
          uiOutput("meta4_fallback_banner"),
          tabsetPanel(
            tabPanel("Overview",
          fluidRow(
            column(4, wellPanel(h4("Selected Cohort Summary"), tableOutput("summary_table"))),
            column(8, wellPanel(h4("Samples by Cohort & Subset"), plotOutput("cohort_plot", height = 280)))
          ),
          fluidRow(
            column(12, wellPanel(
              h4("OncoPrint: Mutations & Clinical Annotations"),
              plotOutput("oncoprint_plot", height = 500)
            ))
          ),
        ),
        tabPanel("Single-Gene Associations",
          fluidRow(
            column(12, wellPanel(
              h4("Clinical Variable by Mutation"),
              selectInput("clin_var", "Variable:", choices = c("WBC", "Age", "Hemoglobin", "Platelet", "LDH", "BM_blast_percent", "PB_blast_percent")),
              plotOutput("clinical_plot", height = 400)))
          ),
          fluidRow(
            column(6, wellPanel(h4("Kaplan-Meier"), selectInput("surv_gene", "Gene:", choices = NULL), plotOutput("survival_plot", height = 500))),
            column(6, wellPanel(h4("Hazard Ratios"), plotOutput("forest_plot", height = 500)))
          ),
          fluidRow(column(12, wellPanel(h4("Survival Summary"), DTOutput("survival_table"))))
        ),
        tabPanel("Co-mutation Associations",
          fluidRow(
            column(6, wellPanel(h4("Mutation Co-occurrence"), plotOutput("cooccurrence_plot", height = 500))),
            column(6, wellPanel(h4("Top Co-occurring Pairs"), div(style = "height: 500px; overflow-y: auto;", DTOutput("cooccurrence_table"))))
          ),
          fluidRow(
            column(4, wellPanel(
              h4("Gene Selection"),
              radioButtons("comut_n", "Number of genes:", choices = c("2 genes" = 2, "3 genes" = 3), selected = 2),
              selectInput("comut_gene1", "Gene 1:", choices = NULL),
              selectInput("comut_gene2", "Gene 2:", choices = NULL),
              conditionalPanel(condition = "input.comut_n == 3",
                selectInput("comut_gene3", "Gene 3:", choices = NULL))
            )),
            column(8, wellPanel(
              h4("Kaplan-Meier by Co-mutation Status"),
              plotOutput("comut_survival_plot", height = 500)
            ))
          )
        ),
        tabPanel("VAF Associations",
          tabsetPanel(
            tabPanel("Single Gene",
              fluidRow(
                column(6, wellPanel(h4("VAF by Gene"), plotOutput("vaf_gene_plot", height = 450))),
                column(6, wellPanel(
                  h4("VAF and Survival (MaxStat)"),
                  p(em("Optimal VAF threshold from maximally selected rank statistics (maxstat); log-rank. High vs Low VAF split among mutated patients.")),
                  selectInput("vaf_surv_gene", "Gene:", choices = NULL),
                  plotOutput("vaf_survival_plot", height = 450)))
              ),
              fluidRow(
                column(6, wellPanel(h4("VAF by Mutation Category"), plotOutput("vaf_category_plot", height = 450))),
                column(6, wellPanel(h4("VAF by Cohort"), plotOutput("vaf_cohort_plot", height = 450)))
              )
            ),
            tabPanel("Pairwise Gene",
              fluidRow(
                column(12, wellPanel(h4("VAF Scatter (Clonal Ordering)"),
                  p("Compare VAF of two genes. Points above line = Gene 1 before Gene 2."),
                  selectInput("vaf_gene1", "Gene 1:", choices = NULL), selectInput("vaf_gene2", "Gene 2:", choices = NULL),
                  plotOutput("vaf_scatter_plot", height = 380)))
              )
            )
          )
        ),
        tabPanel("Drug Sensitivity",
          p("VAF and mutation correlations with inhibitor sensitivity from BeatAML2 data.",
            a("biodev.github.io/BeatAML2", href = "https://biodev.github.io/BeatAML2/", target = "_blank")),
          fluidRow(
            column(12, wellPanel(
              h4("Summary"),
              div(id = "drug_subset_wrap", selectInput("drug_subset", "AML Subset:", choices = c("All", "de_novo", "secondary"), selected = "de_novo")),
              uiOutput("drug_summary_ui"),
              plotOutput("drug_summary_dotplot", height = 500),
              p(em("* q < 0.1 (FDR)"), style = "font-size: 11px; color: #666;")
            ))
          ),
          fluidRow(
            column(6, wellPanel(
              h4("VAF vs AUC Scatter"),
              fluidRow(
                column(6, selectInput("drug_gene", "Gene for scatter:", choices = NULL)),
                column(6, selectInput("drug_inhibitor", "Inhibitor for scatter:", choices = NULL))
              ),
              plotOutput("drug_scatter_plot", height = 350)
            )),
            column(6, wellPanel(
              h4("Leave-One-Out Cross Validation"),
              p(em("RMSE = leave-one-out prediction error in AUC units.")),
              plotOutput("drug_loo_heatmap", height = 350)))
          ),
          fluidRow(
            column(12, wellPanel(h4("Top Correlations"), p(em("LOOCV RMSE = leave-one-out prediction error in AUC units (overfitting check)."), style = "font-size: 11px; color: #666;"), div(style = "height: 350px; overflow-y: auto;", DTOutput("drug_correlation_table"))))
          )
        ),
        tabPanel("Data Table",
          fluidRow(column(12, wellPanel(h4("Mutation Data"), DTOutput("data_table"))))
        )
      )
    )
    )
  )
  )
)

# Server
server <- function(input, output, session) {

  # Data source and Meta4 flag must be reactives on input$main_nav so outputs re-run
  # in the same flush as the tab change (no observer race).
  current_data_src <- reactive({
    req(input$main_nav)
    if (identical(input$main_nav, "meta_aml4")) AML_Meta_Cohort else final_data_matrix
  })

  is_meta4 <- reactive(identical(input$main_nav, "meta_aml4"))

  output$meta4_fallback_banner <- renderUI({
    if (identical(input$main_nav, "meta_aml4") && !meta4_uses_dedicated_file) {
      div(class = "alert alert-warning", role = "alert", style = "margin-bottom: 1em;",
        strong("Meta AML4 data file not found."), " Add ", code("AML_Meta_Cohort.rds"), " (or ", code("AML_Meta_Cohort.RData"), ") to the app directory and redeploy. This tab is currently showing the Benard et al. (2021) dataset."
      )
    }
  })

  observeEvent(input$main_nav, {
    nav <- input$main_nav
    if (identical(nav, "meta_aml4")) {
      updateSelectInput(session, "cohort", choices = c("All", sort(unique(as.character(na.omit(AML_Meta_Cohort$Cohort))))), selected = "All")
    } else if (identical(nav, "analyses")) {
      updateSelectInput(session, "cohort", choices = c("All", sort(unique(as.character(na.omit(final_data_matrix$Cohort))))), selected = "All")
    }
  }, ignoreNULL = TRUE, ignoreInit = FALSE)

  # Use effective cohort so that when switching tabs the first render shows correct data
  # (stale cohort value from other tab is treated as "All" until dropdown updates)
  effective_cohort <- reactive({
    req(current_data_src())
    cohorts_in_data <- sort(unique(as.character(na.omit(current_data_src()$Cohort))))
    co <- input$cohort
    if (is.null(co) || co == "" || !co %in% c("All", cohorts_in_data)) "All" else co
  })

  filtered_data <- reactive({
    req(current_data_src())
    df <- prepare_data(current_data_src())
    if (input$subset != "All") df <- df[as.character(df$Subset) == input$subset, , drop = FALSE]
    if (effective_cohort() != "All") df <- df[as.character(df$Cohort) == effective_cohort(), , drop = FALSE]
    gene_counts <- table(as.character(df$Gene))
    keep_genes <- names(gene_counts)[gene_counts >= input$min_freq]
    df <- df[as.character(df$Gene) %in% keep_genes, , drop = FALSE]
    as.data.frame(df, stringsAsFactors = FALSE)
  })

  # Tell client to hide loading overlay when data is ready (overlay shown on Benard et al. (2021) / 2026 tab click in JS)
  observe({
    req(filtered_data(), input$main_nav %in% c("analyses", "meta_aml4"))
    session$sendCustomMessage("dataReady", input$main_nav)
  })

  observe({
    req(df <- filtered_data())
    genes <- sort(unique(as.character(df$Gene)))
    updateSelectInput(session, "surv_gene", choices = c("Select..." = "", genes))
    updateSelectInput(session, "vaf_gene1", choices = c("Select..." = "", genes))
    updateSelectInput(session, "vaf_gene2", choices = c("Select..." = "", genes))
    updateSelectInput(session, "comut_gene1", choices = c("Select..." = "", genes))
    updateSelectInput(session, "comut_gene2", choices = c("Select..." = "", genes))
    updateSelectInput(session, "comut_gene3", choices = c("Select..." = "", genes))
    updateSelectInput(session, "vaf_surv_gene", choices = c("Select..." = "", genes))
  })

  output$summary_table <- renderTable({
    req(df <- filtered_data())
    samples_surv <- length(unique(df$Sample[!is.na(df$Time_to_OS) & !is.na(df$Censor)]))
    data.frame(
      Metric = c("Samples", "Mutations", "Genes", "With survival data"),
      Value = c(length(unique(df$Sample)), nrow(df), length(unique(df$Gene)), samples_surv)
    )
  }, colnames = FALSE, width = "100%")

  output$cohort_plot <- renderPlot({
    req(df <- filtered_data())
    udf <- unique(df[c("Sample", "Cohort", "Subset")])
    tbl <- table(udf$Cohort, udf$Subset)
    tbl_long <- as.data.frame(tbl)
    names(tbl_long) <- c("Cohort", "Subset", "n")
    cohort_totals <- aggregate(n ~ Cohort, data = tbl_long, sum)
    subset_pal <- if (is_meta4()) META4_SUBSET_COL else PAL_SUBSET
    ggplot(tbl_long, aes(x = Cohort, y = n, fill = Subset)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(data = cohort_totals, aes(x = Cohort, y = n, label = n), vjust = -0.5, size = 5, inherit.aes = FALSE) +
      scale_fill_manual(values = subset_pal, na.value = "gray70") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(x = NULL, y = "Number of samples") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.text = element_text(size = 13), legend.title = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
  })

  oncoprint_data <- reactive({
    req(df <- filtered_data())
    final_data_matrix_3 <- df
    if (!is_meta4()) {
      for (i in seq_len(nrow(final_data_matrix_3))) {
        if (identical(as.character(final_data_matrix_3$variant_type[i]), "other")) {
          final_data_matrix_3$variant_type[i] <- "Unknown"
        }
      }
    }
    gene_tbl <- table(as.character(final_data_matrix_3$Gene))
    keep_genes <- names(gene_tbl)[gene_tbl >= input$min_freq]
    final_data_matrix_3 <- final_data_matrix_3[as.character(final_data_matrix_3$Gene) %in% keep_genes, , drop = FALSE]
    if (nrow(final_data_matrix_3) == 0) return(NULL)
    top_genes <- names(sort(gene_tbl[keep_genes], decreasing = TRUE))
    final_data_matrix_3 <- final_data_matrix_3[as.character(final_data_matrix_3$Gene) %in% top_genes, , drop = FALSE]
    samples <- unique(as.character(final_data_matrix_3$Sample))
    final_data_matrix_3 <- final_data_matrix_3[final_data_matrix_3$Sample %in% samples, , drop = FALSE]
    if (nrow(final_data_matrix_3) == 0) return(NULL)
    c <- length(samples)
    r <- length(top_genes)
    temp_dat <- as.data.frame(matrix(NA_character_, nrow = r, ncol = c))
    colnames(temp_dat) <- samples
    rownames(temp_dat) <- top_genes
    final_data_matrix_3$variant_type <- as.character(final_data_matrix_3$variant_type)
    for (i in seq_len(nrow(final_data_matrix_3))) {
      pt <- final_data_matrix_3$Sample[i]
      gene <- final_data_matrix_3$Gene[i]
      var_type <- final_data_matrix_3$variant_type[i]
      if (is.na(var_type) || var_type == "") var_type <- "Unknown"
      k <- match(pt, colnames(temp_dat))
      l <- match(gene, rownames(temp_dat))
      if (!is.na(k) && !is.na(l)) temp_dat[l, k] <- var_type
    }
    temp_dat[is.na(temp_dat)] <- ""
    temp_dat <- as.matrix(temp_dat)
    clin_cols <- c("Sample", "Cohort", "Sex", "Risk", "Subset", "Time_to_OS")
    clin_cols <- clin_cols[clin_cols %in% colnames(final_data_matrix_3)]
    if (length(clin_cols) < 2) return(list(mut = NULL, temp_dat = temp_dat, anno_df = NULL, samples = samples, genes = top_genes))
    anno_df <- unique(final_data_matrix_3[clin_cols])
    anno_df <- anno_df[!duplicated(anno_df$Sample), , drop = FALSE]
    colnames(anno_df)[colnames(anno_df) == "Time_to_OS"] <- "Survival"
    if ("Survival" %in% colnames(anno_df)) {
      for (i in seq_len(nrow(anno_df))) {
        if (!is.na(anno_df$Survival[i]) && anno_df$Survival[i] != "") {
          anno_df$Survival[i] <- "Yes"
        }
      }
      anno_df$Survival[is.na(anno_df$Survival) | anno_df$Survival == ""] <- "No"
    }
    anno_df <- anno_df[match(samples, anno_df$Sample), , drop = FALSE]
    list(mut = final_data_matrix_3[, c("Sample", "Gene", "variant_type"), drop = FALSE], temp_dat = temp_dat, anno_df = anno_df, samples = samples, genes = top_genes)
  })

  output$oncoprint_plot <- renderPlot({
    od <- oncoprint_data()
    if (is.null(od) || !is.matrix(od$temp_dat) || nrow(od$temp_dat) == 0 || ncol(od$temp_dat) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    use_meta4 <- is_meta4()
    if (has_ComplexHeatmap) {
      temp_dat <- od$temp_dat
      anno_df <- od$anno_df
      samples <- od$samples
      genes <- od$genes
      col <- if (use_meta4) META4_VAR_COL else ONCO_VAR_COL
      vtypes <- unique(as.character(temp_dat[temp_dat != ""]))
      for (vt in setdiff(vtypes, names(col))) col[vt] <- "#80796BFF"
      # Vectorized alter_fun for performance with >100 rows/columns (avoids cell_fun warning)
      alter_fun_list <- list(
        background = function(x, y, w, h) NULL,
        Deletion = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Deletion"], col = "#374E55FF")) },
        INDEL = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["INDEL"], col = "#DF8F44FF")) },
        Indel = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Indel"], col = "#DF8F44FF")) },
        Insertion = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Insertion"], col = "#00A1D5FF")) },
        ITD = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["ITD"], col = "#79AF97FF")) },
        SNV = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["SNV"], col = "#B24745FF")) },
        Splicing = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Splicing"], col = "#6A6599FF")) },
        PTD = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["PTD"], col = "tan")) },
        Other = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Other"], col = "#80796BFF")) },
        Unknown = function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = col["Unknown"], col = "#80796BFF")) }
      )
      for (vt in setdiff(vtypes, names(alter_fun_list))) {
        alter_fun_list[[vt]] <- (function(fc) function(x, y, w, h) { if (length(x) > 0) grid.rect(x, y, w * 0.5, h * 0.9, gp = gpar(fill = fc, col = fc)) })(col[vt])
      }
      if (is.null(anno_df) || nrow(anno_df) == 0) {
        ha <- NULL
      } else {
        Sex <- anno_df[, "Sex"]
        Cohort <- anno_df[, "Cohort"]
        Risk <- anno_df[, "Risk"]
        Subset <- anno_df[, "Subset"]
        Survival <- if ("Survival" %in% colnames(anno_df)) anno_df[, "Survival"] else rep("No", nrow(anno_df))
        Sex[is.na(Sex)] <- "Unknown"; Cohort[is.na(Cohort)] <- "Unknown"; Risk[is.na(Risk)] <- "Unknown"
        Subset[is.na(Subset)] <- "Unknown"
        cohort_col <- if (use_meta4) c(META4_COHORT_COL, Unknown = "gray70") else c(ONCO_COHORT_COL, Unknown = "gray70")
        for (co in unique(Cohort)) if (!co %in% names(cohort_col)) cohort_col[co] <- "gray70"
        subset_col <- if (use_meta4) META4_SUBSET_COL else ONCO_SUBSET_COL
        ha <- ComplexHeatmap::HeatmapAnnotation(
          Sex = Sex, Cohort = Cohort, Risk = Risk, Subset = Subset, Survival = Survival,
          col = list(
            Sex = c("Male" = "#6a51a3", "Female" = "#43a2ca", "Unknown" = "gray90"),
            Survival = c("Yes" = "#252525", "No" = "#f0f0f0"),
            Risk = c("Adverse" = "#E64B35FF", "Intermediate" = "#8491B4FF", "Favorable" = "#00A087FF", "Unknown" = "#767676FF"),
            Subset = subset_col,
            Cohort = cohort_col
          ),
          annotation_height = grid::unit(c(3.5, 3.5, 3.5, 3.5, 3.5), "mm"),
          show_annotation_name = TRUE,
          annotation_legend_param = list(
            Sex = list(title = "Sex"), Survival = list(title = "Survival"),
            Risk = list(title = "Risk"), Subset = list(title = "Subset"), Cohort = list(title = "Cohort")
          )
        )
      }
      fig1B <- ComplexHeatmap::oncoPrint(temp_dat,
        col = col,
        row_names_side = "right",
        bottom_annotation = ha,
        alter_fun_is_vectorized = TRUE,
        alter_fun = alter_fun_list
      )
      ComplexHeatmap::draw(fig1B)
    } else {
      mut <- od$mut; clin <- od$anno_df; samples <- od$samples; genes <- od$genes
      mut$Gene <- factor(as.character(mut$Gene), levels = rev(genes))
      mut$Sample <- factor(as.character(mut$Sample), levels = samples)
      mut$Sample_num <- as.numeric(mut$Sample)
      mut$Gene_num <- as.numeric(mut$Gene)
      PAL_VAR <- if (use_meta4) c(Deletion = "#374E55", Indel = "#DF8F44", Insertion = "#00A1D5", ITD = "#79AF97",
        SNV = "#B24745", Splicing = "#6A6599", PTD = "tan", Other = "#80796B", Unknown = "#80796B") else
        c(Deletion = "#374E55", INDEL = "#DF8F44", Insertion = "#00A1D5", ITD = "#79AF97", SNV = "#B24745", Splicing = "#6A6599", Unknown = "#80796B")
      for (vt in setdiff(unique(mut$variant_type), names(PAL_VAR))) PAL_VAR[vt] <- "#80796B"
      bg <- expand.grid(Sample = factor(samples, levels = samples), Gene = factor(genes, levels = rev(genes)))
      bg$Sample_num <- as.numeric(bg$Sample)
      bg$Gene_num <- as.numeric(bg$Gene)
      p_mut <- ggplot() +
        geom_tile(data = bg, aes(x = Sample_num, y = Gene_num), fill = "gray97", width = 1, height = 1) +
        geom_rect(data = mut, aes(xmin = Sample_num - 0.25, xmax = Sample_num + 0.25, ymin = Gene_num - 0.45, ymax = Gene_num + 0.45, fill = variant_type)) +
        scale_fill_manual(values = PAL_VAR, name = "Variant") +
        scale_x_continuous(limits = c(0.5, length(samples) + 0.5), expand = c(0, 0), position = "top") +
        scale_y_continuous(breaks = seq_along(genes), labels = rev(genes), expand = c(0, 0)) +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "right",
          plot.margin = margin(5, 5, 2, 5))
      if (nrow(clin) > 0 && "Survival" %in% colnames(clin)) {
        annot_cols <- c("Sex", "Cohort", "Risk", "Subset", "Survival")
        annot_cols <- annot_cols[annot_cols %in% colnames(clin)]
        clin_long <- do.call(rbind, lapply(annot_cols, function(col) {
          data.frame(Sample = clin$Sample, annot = col, value = as.character(clin[[col]]), stringsAsFactors = FALSE)
        }))
        clin_long$Sample <- factor(as.character(clin_long$Sample), levels = samples)
        clin_long$Sample_num <- as.numeric(clin_long$Sample)
        clin_long$annot <- factor(clin_long$annot, levels = rev(annot_cols))
        clin_long$value[is.na(clin_long$value) | clin_long$value == ""] <- "NA"
        pal_clin <- if (use_meta4) c(META4_SUBSET_COL, META4_COHORT_COL, PAL_RISK, Male = "#6a51a3", Female = "#43a2ca", Yes = "#252525", No = "#f0f0f0") else c(PAL_SUBSET, PAL_COHORT, PAL_RISK, Male = "#6a51a3", Female = "#43a2ca", Yes = "#252525", No = "#f0f0f0")
        p_clin <- ggplot(clin_long, aes(x = Sample_num, y = annot, fill = value)) +
          geom_tile(color = "white", linewidth = 0.2) +
          scale_fill_manual(values = pal_clin, na.value = "gray90", name = "Clinical") +
          scale_x_continuous(limits = c(0.5, length(samples) + 0.5), expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +
          labs(x = NULL, y = NULL) +
          theme_minimal(base_size = 9) +
          theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "right",
            plot.margin = margin(2, 5, 5, 5))
        if (has_gridExtra) {
          gridExtra::grid.arrange(p_mut, p_clin, ncol = 1, heights = c(3, 1))
        } else {
          print(p_mut)
        }
      } else {
        print(p_mut)
      }
    }
  })

  output$vaf_gene_plot <- renderPlot({
    req(df <- filtered_data())
    df <- df[!is.na(df$VAF) & df$VAF > 0, , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    gene_tbl <- table(as.character(df$Gene))
    top_genes <- names(sort(gene_tbl, decreasing = TRUE))[1:min(25, length(gene_tbl))]
    df <- df[as.character(df$Gene) %in% top_genes, , drop = FALSE]
    n_cat <- if ("mutation_category" %in% colnames(df)) length(unique(df$mutation_category[!is.na(df$mutation_category) & df$mutation_category != ""])) else 0L
    fill_var <- if (n_cat >= 2L) "mutation_category" else "Gene"
    ggplot(df, aes(x = VAF, y = reorder(Gene, VAF, median), fill = .data[[fill_var]])) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_x_continuous(limits = c(0, 100)) +
      (if (fill_var == "mutation_category") scale_fill_manual(values = PAL_MUT_CAT, na.value = "gray70", guide = "legend") else scale_fill_discrete(guide = "none")) +
      labs(x = "Variant Allele Frequency (%)", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.text = element_text(size = 13), legend.title = element_text(size = 14), legend.position = if (fill_var == "mutation_category") "right" else "none")
  })

  output$vaf_cohort_plot <- renderPlot({
    req(df <- filtered_data())
    df <- df[!is.na(df$VAF) & df$VAF > 0 & !is.na(df$Cohort), , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    cohort_pal <- if (is_meta4()) META4_COHORT_COL else PAL_COHORT
    ggplot(df, aes(x = VAF, y = reorder(Cohort, VAF, median), fill = Cohort)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_x_continuous(limits = c(0, 100)) +
      scale_fill_manual(values = cohort_pal, na.value = "gray70") +
      labs(x = "Variant Allele Frequency (%)", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.position = "none")
  })

  output$vaf_category_plot <- renderPlot({
    req(df <- filtered_data())
    df <- df[!is.na(df$VAF) & df$VAF > 0 & !is.na(df$mutation_category) & df$mutation_category != "", , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    ggplot(df, aes(x = VAF, y = reorder(mutation_category, VAF, median), fill = mutation_category)) +
      geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
      scale_x_continuous(limits = c(0, 100)) +
      scale_fill_manual(values = PAL_MUT_CAT, na.value = "gray70") +
      labs(x = "Variant Allele Frequency (%)", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), legend.position = "none")
  })

  vaf_survival_data <- reactive({
    gene <- input$vaf_surv_gene
    if (is.null(gene) || gene == "") return(NULL)
    df <- filtered_data()
    if (!"VAF" %in% colnames(df)) return(NULL)
    df <- df[as.character(df$Gene) == gene & !is.na(df$VAF) & df$VAF > 0, c("Sample", "VAF"), drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    vaf_agg <- aggregate(VAF ~ Sample, data = df, max)
    surv_df <- survival_data()
    surv_df <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    merge_df <- merge(vaf_agg, surv_df, by = "Sample", all = FALSE)
    if (nrow(merge_df) < 20) return(NULL)
    merge_df
  })

  output$vaf_survival_plot <- renderPlot({
    gene <- input$vaf_surv_gene
    if (is.null(gene) || gene == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    vaf_surv_df <- vaf_survival_data()
    if (is.null(vaf_surv_df) || nrow(vaf_surv_df) < 20) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Need at least 20 mutated samples with VAF and survival data. Install the maxstat package for optimal threshold.") + theme_void())
    }
    if (!has_maxstat) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install the maxstat package for optimal VAF threshold (maximally selected rank statistics).") + theme_void())
    }
    mst <- tryCatch(maxstat::maxstat.test(Surv(Time_to_OS, as.numeric(Censor)) ~ VAF, data = vaf_surv_df, smethod = "LogRank", minprop = 0.2, maxprop = 0.8), error = function(e) NULL)
    if (is.null(mst)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "MaxStat failed for this gene.") + theme_void())
    cutpoint <- mst$estimate
    vaf_surv_df$Group <- ifelse(vaf_surv_df$VAF > cutpoint, "High VAF", "Low VAF")
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = vaf_surv_df)
    sd <- survdiff(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = vaf_surv_df)
    pval <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
    pval_txt <- if (pval < 0.001) "p < 0.001" else paste0("p = ", format.pval(pval, digits = 2))
    if (has_survminer) {
      legend_labs <- gsub("^Group=", "", names(fit$strata))
      p <- survminer::ggsurvplot(fit, data = vaf_surv_df, risk.table = TRUE, pval = pval_txt,
        title = paste0(gene, ": High vs Low VAF (threshold ", round(cutpoint, 1), "%)"), xlab = "Years",
        palette = c("High VAF" = "#8B0000", "Low VAF" = "#4D4D4D"), legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.1))
      print(p)
    } else {
      sf_dat <- data.frame(time = fit$time, surv = fit$surv, strata = rep(names(fit$strata), fit$strata))
      ggplot(sf_dat, aes(x = time, y = surv, color = strata)) +
        geom_step(linewidth = 1) +
        scale_color_manual(values = c("High VAF" = "#8B0000", "Low VAF" = "#4D4D4D")) +
        labs(title = paste0(gene, ": High vs Low VAF (", round(cutpoint, 1), "%)"), subtitle = pval_txt, x = "Years", y = "Survival probability") +
        theme_minimal() + theme(legend.position = "right")
    }
  })

  survival_data <- reactive({
    df <- filtered_data()
    df <- df[!is.na(df$Time_to_OS) & !is.na(df$Censor), , drop = FALSE]
    df$Time_to_OS <- as.numeric(df$Time_to_OS) / 365
    df
  })

  output$survival_plot <- renderPlot({
    gene <- input$surv_gene
    if (is.null(gene) || gene == "") return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select a gene"))
    surv_df <- survival_data()
    mut_samples <- unique(surv_df$Sample[surv_df$Gene == gene])
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    surv_uniq$Mutation <- ifelse(surv_uniq$Sample %in% mut_samples, paste0(gene, " mut"), "WT")
    if (length(unique(surv_uniq$Mutation)) < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient groups"))
    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Mutation, data = surv_uniq)
    strata_names <- names(fit$strata)
    is_mut <- grepl(" mut", strata_names, fixed = TRUE)
    surv_pal <- ifelse(is_mut, "#8B0000", "#4D4D4D")  # deep red for mutated, gray for WT
    if (has_survminer) {
      legend_labs <- gsub("^Mutation=", "", names(fit$strata))
      p <- survminer::ggsurvplot(fit, data = surv_uniq, risk.table = TRUE, pval = TRUE,
        title = paste("Survival by", gene, "mutation status"), xlab = "Years",
        palette = surv_pal, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.1))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install the survminer package for Kaplan-Meier plots with number at risk.", size = 4) + theme_void()
    }
  })

  forest_data <- reactive({
    surv_df <- survival_data()
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    gene_tbl <- table(as.character(surv_df$Gene))
    genes <- names(gene_tbl)[gene_tbl >= input$min_freq]
    res <- list()
    for (g in genes) {
      mut_pts <- unique(surv_df$Sample[surv_df$Gene == g])
      surv_uniq$mut <- surv_uniq$Sample %in% mut_pts
      if (sum(surv_uniq$mut) >= 10 && sum(!surv_uniq$mut) >= 10) {
        m <- tryCatch(coxph(Surv(Time_to_OS, as.numeric(Censor)) ~ mut, data = surv_uniq), error = function(e) NULL)
        if (!is.null(m)) {
          ci <- confint(m)
          res[[g]] <- data.frame(Gene = g, HR = exp(coef(m))["mutTRUE"], lower = exp(ci[1]), upper = exp(ci[2]),
            p = summary(m)$coefficients["mutTRUE", "Pr(>|z|)"], stringsAsFactors = FALSE)
        }
      }
    }
    out <- do.call(rbind, res)
    if (!is.null(out)) out <- out[!is.na(out$HR) & out$lower < 10 & out$upper < 10, , drop = FALSE]
    else out <- data.frame(Gene = character(), HR = numeric(), lower = numeric(), upper = numeric(), p = numeric())
    out
  })

  # Co-mutation survival analysis
  comut_survival_data <- reactive({
    n_genes <- as.numeric(input$comut_n)
    g1 <- input$comut_gene1
    g2 <- input$comut_gene2
    g3 <- if (n_genes == 3) input$comut_gene3 else NULL
    if (is.null(g1) || g1 == "" || is.null(g2) || g2 == "") return(NULL)
    if (n_genes == 3 && (is.null(g3) || g3 == "")) return(NULL)
    genes <- c(g1, g2)
    if (n_genes == 3) genes <- c(g1, g2, g3)
    if (length(unique(genes)) != length(genes)) return(NULL)  # duplicate genes

    surv_df <- survival_data()
    surv_uniq <- surv_df[!duplicated(surv_df$Sample), c("Sample", "Time_to_OS", "Censor")]
    mut_list <- list()
    for (g in genes) {
      mut_list[[g]] <- unique(surv_df$Sample[surv_df$Gene == g])
    }

    if (n_genes == 2) {
      has_g1 <- surv_uniq$Sample %in% mut_list[[g1]]
      has_g2 <- surv_uniq$Sample %in% mut_list[[g2]]
      surv_uniq$Group <- ifelse(!has_g1 & !has_g2, "Neither",
        ifelse(has_g1 & !has_g2, paste0(g1, " only"),
          ifelse(!has_g1 & has_g2, paste0(g2, " only"), paste0(g1, " + ", g2))))
    } else {
      has_g1 <- surv_uniq$Sample %in% mut_list[[g1]]
      has_g2 <- surv_uniq$Sample %in% mut_list[[g2]]
      has_g3 <- surv_uniq$Sample %in% mut_list[[g3]]
      n_mut <- has_g1 + has_g2 + has_g3
      surv_uniq$Group <- character(nrow(surv_uniq))
      surv_uniq$Group[n_mut == 0] <- "None"
      surv_uniq$Group[n_mut == 1 & has_g1] <- paste0(g1, " only")
      surv_uniq$Group[n_mut == 1 & has_g2] <- paste0(g2, " only")
      surv_uniq$Group[n_mut == 1 & has_g3] <- paste0(g3, " only")
      surv_uniq$Group[n_mut == 2 & has_g1 & has_g2] <- paste0(g1, " + ", g2)
      surv_uniq$Group[n_mut == 2 & has_g1 & has_g3] <- paste0(g1, " + ", g3)
      surv_uniq$Group[n_mut == 2 & has_g2 & has_g3] <- paste0(g2, " + ", g3)
      surv_uniq$Group[n_mut == 3] <- paste0(g1, " + ", g2, " + ", g3)
    }
    surv_uniq
  })

  output$comut_survival_plot <- renderPlot({
    req(df <- comut_survival_data())
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    n_groups <- length(unique(df$Group))
    if (n_groups < 2) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Need at least 2 groups for comparison"))
    if (any(table(df$Group) < 3)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Each group needs at least 3 patients"))

    fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = df)
    sd <- survdiff(Surv(Time_to_OS, as.numeric(Censor)) ~ Group, data = df)
    pval <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
    pval_txt <- if (pval < 0.001) "p < 0.001" else paste0("p = ", format.pval(pval, digits = 2))

    if (has_survminer) {
      n_groups <- length(unique(df$Group))
      pal_all <- unname(c(PAL_SUBSET, PAL_MUT_CAT, PAL_COHORT))
      pal_vec <- pal_all[seq_len(n_groups)]
      if (length(pal_vec) < n_groups) {
        pal_vec <- rep(pal_all, length.out = n_groups)
      }
      legend_labs <- gsub("^Group=", "", names(fit$strata))
      p <- survminer::ggsurvplot(fit, data = df, risk.table = TRUE, pval = pval_txt,
        title = "Survival by Co-mutation Status", xlab = "Years",
        palette = pal_vec, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.1))
      print(p)
    } else {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Install the survminer package for Kaplan-Meier plots with number at risk.", size = 4) + theme_void()
    }
  })

  output$forest_plot <- renderPlot({
    df <- forest_data()
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    df <- df[order(df$HR), , drop = FALSE]
    df$Gene <- factor(df$Gene, levels = df$Gene)
    ggplot(df, aes(x = Gene, y = HR, ymin = lower, ymax = upper, color = HR > 1)) +
      geom_pointrange(size = 0.5) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      coord_flip() +
      scale_color_manual(values = c("TRUE" = "#762a83", "FALSE" = "#1b7837"), guide = "none") +
      labs(x = NULL, y = "Hazard Ratio") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15))
  })

  output$survival_table <- renderDT({
    df <- forest_data()
    if (nrow(df) == 0) return(NULL)
    df$HR_CI <- sprintf("%.2f (%.2f-%.2f)", df$HR, df$lower, df$upper)
    df$p_adj <- p.adjust(df$p, method = "fdr")
    datatable(df[, c("Gene", "HR_CI", "p", "p_adj")], options = list(pageLength = 15), rownames = FALSE)
  })

  cooccurrence_data <- reactive({
    df <- filtered_data()[, c("Sample", "Gene")]
    df <- df[!duplicated(df), , drop = FALSE]
    genes <- names(which(table(as.character(df$Gene)) >= 15))
    if (length(genes) < 2) return(list(matrix = NULL, pairs = NULL))
    df <- df[as.character(df$Gene) %in% genes, , drop = FALSE]
    wide <- as.data.frame.matrix(xtabs(~ Sample + Gene, data = df))
    wide[wide > 0] <- 1
    gene_cols <- colnames(wide)
    if (length(gene_cols) < 2) return(list(matrix = NULL, pairs = NULL))
    pairs <- t(combn(gene_cols, 2))
    results <- list()
    for (i in seq_len(nrow(pairs))) {
      g1 <- pairs[i, 1]
      g2 <- pairs[i, 2]
      tbl <- table(wide[[g1]], wide[[g2]])
      if (any(dim(tbl) < 2)) next
      tryCatch({
        ft <- fisher.test(tbl)
        or <- (tbl[2,2] * tbl[1,1]) / (tbl[2,1] * tbl[1,2])
        if (!is.finite(or)) or <- 0.01
        results[[i]] <- data.frame(gene1 = g1, gene2 = g2, odds_ratio = or, n_cooccur = tbl[2,2], p = ft$p.value, stringsAsFactors = FALSE)
      }, error = function(e) NULL)
    }
    pairs_df <- do.call(rbind, results)
    if (is.null(pairs_df) || nrow(pairs_df) == 0) return(list(matrix = NULL, pairs = NULL))
    pairs_df$q <- p.adjust(pairs_df$p, method = "fdr")
    or_mat <- matrix(1, length(genes), length(genes))
    rownames(or_mat) <- colnames(or_mat) <- genes
    for (i in seq_len(nrow(pairs_df))) {
      g1 <- pairs_df$gene1[i]
      g2 <- pairs_df$gene2[i]
      or_mat[g1, g2] <- or_mat[g2, g1] <- pairs_df$odds_ratio[i]
    }
    list(matrix = or_mat, pairs = pairs_df)
  })

  output$cooccurrence_plot <- renderPlot({
    cooc <- cooccurrence_data()
    if (is.null(cooc$matrix) || nrow(cooc$matrix) < 3) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data"))
    mat <- cooc$matrix
    mat_log <- log2(mat + 0.01)
    mat_log[mat_log > 2] <- 2
    mat_log[mat_log < -2] <- -2
    ord <- hclust(dist(mat_log))$order
    mat_ord <- mat_log[ord, ord]
    df_plot <- data.frame(Gene1 = rep(rownames(mat_ord), ncol(mat_ord)),
      Gene2 = rep(colnames(mat_ord), each = nrow(mat_ord)),
      log2OR = as.vector(mat_ord))
    df_plot$Gene1 <- factor(df_plot$Gene1, levels = rownames(mat_ord))
    df_plot$Gene2 <- factor(df_plot$Gene2, levels = colnames(mat_ord))
    ggplot(df_plot, aes(x = Gene1, y = Gene2, fill = log2OR)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
      labs(x = NULL, y = NULL, fill = "log2(OR)") +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title = element_text(size = 14), legend.text = element_text(size = 12))
  })

  output$cooccurrence_table <- renderDT({
    cooc <- cooccurrence_data()
    if (is.null(cooc$pairs) || nrow(cooc$pairs) == 0) return(NULL)
    df <- cooc$pairs
    df$pair <- paste(df$gene1, "+", df$gene2)
    df <- df[order(-df$n_cooccur), , drop = FALSE]
    df <- head(df, 30)
    datatable(df[, c("pair", "odds_ratio", "n_cooccur", "p", "q")], options = list(pageLength = 10), rownames = FALSE)
  })

  # Figure 2: Co-occurrence & Prognosis (de novo)
  fig2_data <- reactive({
    df <- filtered_data()
    df <- df[as.character(df$Subset) == "de_novo", , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df <- df[!is.na(df$Time_to_OS) & !is.na(df$Censor), , drop = FALSE]
    df$Time_to_OS <- as.numeric(df$Time_to_OS) / 365
    genes <- names(which(table(as.character(df$Gene)) >= 15))
    if (length(genes) < 2) return(NULL)
    df <- df[as.character(df$Gene) %in% genes, , drop = FALSE]
    wide <- as.data.frame.matrix(xtabs(~ Sample + Gene, data = unique(df[, c("Sample", "Gene")])))
    wide[wide > 0] <- 1
    surv_uniq <- df[!duplicated(df$Sample), c("Sample", "Time_to_OS", "Censor")]
    pairs <- t(combn(genes, 2))
    or_list <- list()
    hr_list <- list()
    for (i in seq_len(nrow(pairs))) {
      g1 <- pairs[i, 1]
      g2 <- pairs[i, 2]
      tbl <- table(wide[[g1]], wide[[g2]])
      if (any(dim(tbl) < 2)) next
      ft <- tryCatch(fisher.test(tbl), error = function(e) NULL)
      if (is.null(ft)) next
      or <- (tbl[2,2] * tbl[1,1]) / (tbl[2,1] * tbl[1,2])
      if (!is.finite(or) || or == 0) or <- 0.01
      n_co <- tbl[2, 2]
      if (n_co < 10) next
      or_list[[length(or_list) + 1]] <- data.frame(gene1 = g1, gene2 = g2, odds_ratio = or, p_or = ft$p.value, n_cooccur = n_co, stringsAsFactors = FALSE)
      has_both <- wide[[g1]] == 1 & wide[[g2]] == 1
      has_neither <- wide[[g1]] == 0 & wide[[g2]] == 0
      mut_samp <- rownames(wide)[has_both]
      wt_samp <- rownames(wide)[has_neither]
      surv_uniq$group <- ifelse(surv_uniq$Sample %in% mut_samp, "co_mut", ifelse(surv_uniq$Sample %in% wt_samp, "WT", NA))
      surv_sub <- surv_uniq[!is.na(surv_uniq$group), , drop = FALSE]
      if (sum(surv_sub$group == "co_mut") >= 10 && sum(surv_sub$group == "WT") >= 10) {
        m <- tryCatch(coxph(Surv(Time_to_OS, as.numeric(Censor)) ~ factor(group, levels = c("WT", "co_mut")), data = surv_sub), error = function(e) NULL)
        if (!is.null(m)) {
          ci <- confint(m)
          hr_list[[length(hr_list) + 1]] <- data.frame(gene1 = g1, gene2 = g2, HR = exp(coef(m))[1], lower = exp(ci[1]), upper = exp(ci[2]),
            p_hr = summary(m)$coefficients[1, 5], n_cooccur = n_co, stringsAsFactors = FALSE)
        }
      }
    }
    or_df <- if (length(or_list) > 0) { o <- do.call(rbind, or_list); o$q_or <- p.adjust(o$p_or, "fdr"); o } else NULL
    hr_df <- if (length(hr_list) > 0) { h <- do.call(rbind, hr_list); h$q_hr <- p.adjust(h$p_hr, "fdr"); h } else NULL
    merged <- NULL
    if (!is.null(or_df) && !is.null(hr_df)) {
      merged <- merge(or_df[, c("gene1", "gene2", "odds_ratio", "q_or", "n_cooccur")],
        hr_df[, c("gene1", "gene2", "HR", "q_hr", "lower", "upper")],
        by = c("gene1", "gene2"), all = TRUE)
    }
    list(or = or_df, hr = hr_df, merged = merged, wide = wide, surv = surv_uniq, df = df, genes = genes)
  })

  fig2_triple_data <- reactive({
    fd <- fig2_data()
    if (is.null(fd) || length(fd$genes) < 3) return(NULL)
    df <- fd$df
    wide <- fd$wide
    surv_uniq <- fd$surv
    genes <- fd$genes
    triples <- t(combn(genes, 3))
    res <- list()
    wide_cols <- colnames(wide)
    for (i in seq_len(nrow(triples))) {
      g1 <- triples[i, 1]; g2 <- triples[i, 2]; g3 <- triples[i, 3]
      if (!all(c(g1, g2, g3) %in% wide_cols)) next
      has_3 <- wide[[g1]] == 1 & wide[[g2]] == 1 & wide[[g3]] == 1
      has_2 <- (wide[[g1]] + wide[[g2]] + wide[[g3]]) == 2
      n3 <- sum(has_3)
      if (n3 < 10) next
      samp_3 <- rownames(wide)[has_3]
      samp_2 <- rownames(wide)[has_2]
      surv_uniq$grp <- ifelse(surv_uniq$Sample %in% samp_3, "triple", ifelse(surv_uniq$Sample %in% samp_2, "double", NA))
      surv_sub <- surv_uniq[surv_uniq$grp %in% c("triple", "double"), , drop = FALSE]
      if (nrow(surv_sub) < 20 || sum(surv_sub$grp == "triple") < 10) next
      m <- tryCatch(coxph(Surv(Time_to_OS, as.numeric(Censor)) ~ factor(grp, levels = c("double", "triple")), data = surv_sub), error = function(e) NULL)
      if (is.null(m)) next
      ci <- confint(m)
      pv <- summary(m)$coefficients[1, 5]
      genotype <- paste(g1, g2, g3, sep = " + ")
      res[[length(res) + 1]] <- data.frame(genotype = genotype, g1 = g1, g2 = g2, g3 = g3, n_3 = n3, HR = exp(coef(m))[1], lower = exp(ci[1]), upper = exp(ci[2]), p = pv, stringsAsFactors = FALSE)
    }
    if (length(res) == 0) return(NULL)
    out <- do.call(rbind, res)
    out$q <- p.adjust(out$p, "fdr")
    out
  })

  output$fig2a_plot <- renderPlot({
    fd <- fig2_data()
    if (is.null(fd$or) || nrow(fd$or) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    d <- fd$or
    d$sig <- "NS"
    d$sig[d$q_or < 0.05 & d$odds_ratio > 1] <- "Co-occurring"
    d$sig[d$q_or < 0.05 & d$odds_ratio < 1] <- "Mutually exclusive"
    d$neglog10q <- -log10(pmax(d$q_or, 1e-10))
    d$logOR <- log(pmax(d$odds_ratio, 0.01))
    ggplot(d, aes(x = logOR, y = neglog10q, color = sig, size = n_cooccur)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70") +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Co-occurring" = "#b2182b", "Mutually exclusive" = "#2166ac", "NS" = "gray70")) +
      scale_size_continuous(range = c(2, 8)) +
      labs(x = "log(Odds Ratio)", y = "-log10(q-value)", color = NULL, size = "n co-mut") +
      theme_minimal(base_size = 11)
  })

  output$fig2b_plot <- renderPlot({
    fd <- fig2_data()
    if (is.null(fd$hr) || nrow(fd$hr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    d <- fd$hr
    d$sig <- "NS"
    d$sig[d$q_hr < 0.05 & d$HR <= 1] <- "Favorable"
    d$sig[d$q_hr < 0.05 & d$HR > 1] <- "Unfavorable"
    d$neglog10q <- -log10(pmax(d$q_hr, 1e-10))
    ggplot(d, aes(x = HR, y = neglog10q, color = sig, size = n_cooccur)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray70") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70") +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Favorable" = "#1b7837", "Unfavorable" = "#762a83", "NS" = "gray70")) +
      scale_size_continuous(range = c(2, 8)) +
      labs(x = "Hazard Ratio", y = "-log10(q-value)", color = NULL, size = "n co-mut") +
      theme_minimal(base_size = 11)
  })

  output$fig2c_plot <- renderPlot({
    fd <- fig2_data()
    if (is.null(fd$merged) || nrow(fd$merged) < 5) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data"))
    d <- fd$merged
    d <- d[complete.cases(d[, c("odds_ratio", "HR")]), , drop = FALSE]
    if (nrow(d) < 5) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient overlap"))
    d$logOR <- log(pmax(d$odds_ratio, 0.01))
    d$logHR <- log(pmax(d$HR, 0.01))
    d$sig_hr <- "NS"
    d$sig_hr[!is.na(d$q_hr) & d$q_hr < 0.05 & d$HR <= 1] <- "Favorable"
    d$sig_hr[!is.na(d$q_hr) & d$q_hr < 0.05 & d$HR > 1] <- "Unfavorable"
    d$n_cooccur[is.na(d$n_cooccur)] <- 10
    p <- ggplot(d, aes(x = logOR, y = logHR, color = sig_hr, size = n_cooccur)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c("Favorable" = "#1b7837", "Unfavorable" = "#762a83", "NS" = "gray70")) +
      scale_size_continuous(range = c(2, 8)) +
      labs(x = "log(Odds Ratio)", y = "log(Hazard Ratio)", color = NULL, size = "n co-mut") +
      theme_minimal(base_size = 11)
    for (grp in c("Favorable", "Unfavorable", "NS")) {
      sub <- d[d$sig_hr == grp, , drop = FALSE]
      if (nrow(sub) >= 3) {
        fit <- lm(logHR ~ logOR, data = sub)
        xs <- seq(min(d$logOR), max(d$logOR), length.out = 50)
        pred <- predict(fit, newdata = data.frame(logOR = xs), se.fit = TRUE)
        df_band <- data.frame(x = xs, y = pred$fit, ymin = pred$fit - 1.96 * pred$se.fit, ymax = pred$fit + 1.96 * pred$se.fit)
        col <- c("Favorable" = "#1b7837", "Unfavorable" = "#762a83", "NS" = "gray70")[grp]
        p <- p + geom_ribbon(data = df_band, aes(x = x, y = y, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = col, alpha = 0.2) +
          geom_line(data = df_band, aes(x = x, y = y), inherit.aes = FALSE, color = col, linewidth = 0.8)
      }
    }
    p
  })

  output$fig2d_plot <- renderPlot({
    tr <- fig2_triple_data()
    if (is.null(tr) || nrow(tr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    tr <- tr[order(-tr$n_3), , drop = FALSE]
    tr <- head(tr, 25)
    tr$sig_surv <- factor(tr$p <= 0.05, levels = c(FALSE, TRUE), labels = c("NS", "p ≤ 0.05"))
    levs <- rev(unique(as.character(tr$genotype)))
    tr$genotype <- factor(as.character(tr$genotype), levels = levs)
    ggplot(tr, aes(x = genotype, y = n_3, fill = sig_surv)) +
      geom_col() +
      geom_text(aes(label = n_3), hjust = -0.2, size = 3, color = "black") +
      scale_fill_manual(values = c("NS" = "gray70", "p ≤ 0.05" = "#b2182b")) +
      coord_flip() +
      labs(x = NULL, y = "# patients with 3 mutations", fill = "Survival") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "right")
  })

  output$fig2e_plot <- renderPlot({
    tr <- fig2_triple_data()
    if (is.null(tr) || nrow(tr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    sig <- tr[tr$q <= 0.15 & tr$p <= 0.05, , drop = FALSE]
    if (nrow(sig) == 0) sig <- head(tr[tr$p <= 0.05, , drop = FALSE], 6)
    if (nrow(sig) == 0) sig <- head(tr[order(tr$p), ], 6)
    if (nrow(sig) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No genotypes to display"))
    sig$col <- ifelse(sig$HR > 1, "Unfavorable", "Favorable")
    sig$genotype <- factor(sig$genotype, levels = sig$genotype[order(sig$HR)])
    p1 <- ggplot(sig, aes(x = genotype, y = HR, ymin = lower, ymax = pmin(upper, 8), color = col)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      geom_pointrange(size = 0.5) +
      coord_flip() +
      scale_color_manual(values = c("Favorable" = "#1b7837", "Unfavorable" = "#762a83")) +
      labs(x = NULL, y = "Hazard Ratio (3 vs 2 muts)", title = "Forest plot") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none")
    fd <- fig2_data()
    if (is.null(fd)) return(p1)
    if (has_survminer) {
      splots <- list()
      for (i in seq_len(min(3, nrow(sig)))) {
        gs <- unlist(strsplit(as.character(sig$genotype[i]), " + ", fixed = TRUE))
        gs <- gs[gs %in% colnames(fd$wide)]
        if (length(gs) < 3) next
        wide <- fd$wide
        surv <- fd$surv
        has_3 <- rowSums(wide[, gs, drop = FALSE]) == 3
        has_2 <- rowSums(wide[, gs, drop = FALSE]) == 2
        surv$grp <- ifelse(surv$Sample %in% rownames(wide)[has_3], "3 mut", ifelse(surv$Sample %in% rownames(wide)[has_2], "2 mut", NA))
        surv_sub <- surv[!is.na(surv$grp), , drop = FALSE]
        if (nrow(surv_sub) >= 10) {
          fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ grp, data = surv_sub)
          sd <- survdiff(Surv(Time_to_OS, as.numeric(Censor)) ~ grp, data = surv_sub)
          pv <- 1 - pchisq(sd$chisq, 1)
          pv_txt <- if (pv < 0.001) "p < 0.001" else paste0("p = ", format.pval(pv, digits = 2))
          sn <- names(fit$strata)
          scale_vals <- ifelse(grepl("3 mut", sn), "#b2182b", "gray50")
          scale_vals <- setNames(scale_vals, sn)
          legend_labs <- gsub("^grp=", "", sn)
          sp <- survminer::ggsurvplot(fit, data = surv_sub, risk.table = TRUE, pval = pv_txt,
            title = as.character(sig$genotype[i]), xlab = "Years", palette = scale_vals, legend = "right", legend.labs = legend_labs, pval.coord = c(0, 0.1))
          splots[[length(splots) + 1]] <- sp
        }
      }
      if (length(splots) == 0) return(p1)
      km_grob <- survminer::arrange_ggsurvplots(splots, print = FALSE, ncol = 2)
      if (has_gridExtra) {
        gridExtra::grid.arrange(p1, km_grob, ncol = 1, heights = c(0.4, 0.6))
      } else {
        print(p1)
      }
    } else {
      plots <- list(p1)
      for (i in seq_len(min(3, nrow(sig)))) {
        gs <- unlist(strsplit(as.character(sig$genotype[i]), " + ", fixed = TRUE))
        gs <- gs[gs %in% colnames(fd$wide)]
        if (length(gs) < 3) next
        wide <- fd$wide
        surv <- fd$surv
        has_3 <- rowSums(wide[, gs, drop = FALSE]) == 3
        has_2 <- rowSums(wide[, gs, drop = FALSE]) == 2
        surv$grp <- ifelse(surv$Sample %in% rownames(wide)[has_3], "3 mut", ifelse(surv$Sample %in% rownames(wide)[has_2], "2 mut", NA))
        surv_sub <- surv[!is.na(surv$grp), , drop = FALSE]
        if (nrow(surv_sub) >= 10) {
          fit <- survfit(Surv(Time_to_OS, as.numeric(Censor)) ~ grp, data = surv_sub)
          sd <- survdiff(Surv(Time_to_OS, as.numeric(Censor)) ~ grp, data = surv_sub)
          pv <- 1 - pchisq(sd$chisq, 1)
          pv_txt <- if (pv < 0.001) "p < 0.001" else paste0("p = ", format.pval(pv, digits = 2))
          sf <- data.frame(time = fit$time, surv = fit$surv, strata = rep(names(fit$strata), fit$strata))
          sn <- names(fit$strata)
          scale_vals <- setNames(ifelse(grepl("3 mut", sn), "#b2182b", "gray50"), sn)
          pkm <- ggplot(sf, aes(x = time, y = surv, color = strata)) + geom_step(linewidth = 1) +
            scale_color_manual(values = scale_vals) +
            labs(title = as.character(sig$genotype[i]), subtitle = pv_txt, x = "Years", y = "Survival") +
            theme_minimal(base_size = 9)
          plots[[length(plots) + 1]] <- pkm
        }
      }
      if (length(plots) == 1) return(p1)
      if (has_gridExtra) {
        gridExtra::grid.arrange(grobs = plots, ncol = 2)
      } else {
        p1
      }
    }
  })

  output$clinical_plot <- renderPlot({
    var <- input$clin_var
    req(df <- filtered_data())
    cols <- c("Sample", "Gene", var)
    df <- df[, cols[cols %in% colnames(df)], drop = FALSE]
    df <- df[!is.na(df[[var]]), , drop = FALSE]
    if (nrow(df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    sample_vals <- unique(df[, c("Sample", var), drop = FALSE])
    sample_vals$y_val <- suppressWarnings(as.numeric(sample_vals[[var]]))
    sample_vals <- sample_vals[!is.na(sample_vals$y_val), , drop = FALSE]
    if (nrow(sample_vals) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No numeric data"))
    gene_tbl <- table(as.character(df$Gene))
    top_genes <- names(sort(gene_tbl, decreasing = TRUE))
    plot_list <- lapply(top_genes, function(g) {
      mut_samples <- unique(df$Sample[df$Gene == g])
      data.frame(
        Gene = g,
        Sample = sample_vals$Sample,
        Status = ifelse(sample_vals$Sample %in% mut_samples, "mut", "WT"),
        value = sample_vals$y_val,
        stringsAsFactors = FALSE
      )
    })
    plot_df <- do.call(rbind, plot_list)
    med_mut <- vapply(top_genes, function(g) median(plot_df$value[plot_df$Gene == g & plot_df$Status == "mut"], na.rm = TRUE), numeric(1))
    gene_levels <- names(sort(med_mut, decreasing = TRUE))
    plot_df$Gene <- factor(plot_df$Gene, levels = gene_levels)
    y_range <- diff(range(plot_df$value, na.rm = TRUE))
    y_range <- if (y_range == 0) 1 else y_range
    sig_df <- do.call(rbind, lapply(gene_levels, function(g) {
      sub <- plot_df[plot_df$Gene == g, ]
      n_mut <- sum(sub$Status == "mut", na.rm = TRUE)
      n_wt <- sum(sub$Status == "WT", na.rm = TRUE)
      if (n_mut < 2 || n_wt < 2) {
        pval <- NA
      } else {
        pval <- tryCatch(wilcox.test(value ~ Status, data = sub, exact = FALSE)$p.value, error = function(e) NA)
      }
      sym <- if (is.na(pval)) "" else if (pval < 0.001) "***" else if (pval < 0.01) "**" else if (pval < 0.05) "*" else "ns"
      y_pos <- max(sub$value, na.rm = TRUE) + 0.04 * y_range
      data.frame(Gene = factor(g, levels = gene_levels), y = y_pos, label = sym, stringsAsFactors = FALSE)
    }))
    sig_df <- sig_df[!is.na(sig_df$Gene), , drop = FALSE]
    p <- ggplot(plot_df, aes(x = Gene, y = value, fill = Status)) +
      geom_boxplot(alpha = 0.8, outlier.alpha = 0.3, position = position_dodge(0.85), width = 0.7) +
      scale_fill_manual(values = c("mut" = "#8B0000", "WT" = "lightgrey"), name = "Status") +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
      labs(x = NULL, y = var) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(size = 13), legend.title = element_text(size = 14))
    if (nrow(sig_df) > 0) {
      p <- p + geom_text(data = sig_df, aes(x = Gene, y = y, label = label), inherit.aes = FALSE, size = 5, vjust = -0.2)
    }
    p
  })

  output$vaf_scatter_plot <- renderPlot({
    g1 <- input$vaf_gene1
    g2 <- input$vaf_gene2
    if (is.null(g1) || g1 == "" || is.null(g2) || g2 == "" || g1 == g2) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select two different genes"))
    }
    req(df <- filtered_data())
    df <- df[as.character(df$Gene) %in% c(g1, g2), c("Sample", "Gene", "VAF"), drop = FALSE]
    df <- df[!is.na(df$VAF) & df$VAF > 0, , drop = FALSE]
    wide <- reshape(df, idvar = "Sample", timevar = "Gene", direction = "wide")
    colnames(wide) <- gsub("^VAF\\.", "", colnames(wide))
    if (!g1 %in% colnames(wide) || !g2 %in% colnames(wide)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No overlap"))
    wide <- wide[complete.cases(wide[, c(g1, g2), drop = FALSE]), , drop = FALSE]
    if (nrow(wide) < 5) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient samples"))
    diff <- wide[[g1]] - wide[[g2]]
    wide$Clonality <- ifelse(diff > 5, paste(g1, "first"), ifelse(diff < -5, paste(g2, "first"), "Ambiguous"))
    color_vals <- c("Ambiguous" = "#999999", setNames(c("#01665e", "#8c510a"), c(paste(g1, "first"), paste(g2, "first"))))
    ggplot(wide, aes(x = .data[[g2]], y = .data[[g1]], color = Clonality)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_abline(slope = 1, intercept = 5, linetype = "dashed", color = "gray50") +
      geom_abline(slope = 1, intercept = -5, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = color_vals) +
      labs(x = paste(g2, "VAF (%)"), y = paste(g1, "VAF (%)"), color = "Clonal order") +
      theme_minimal(base_size = 12)
  })

  # Drug Sensitivity tab (BeatAML2)
  # Two distinct datasets: 2026 = all samples in inhibitor_auc; 2021 = waves 1+2 only (from clinical)
  beataml_full <- reactive({
    if (!exists("load_beataml2")) return(NULL)
    load_beataml2()
  })

  beataml_wave12 <- reactive({
    b <- beataml_full()
    if (is.null(b) || !b$ok) return(NULL)
    if (exists("subset_beataml_to_wave12")) return(subset_beataml_to_wave12(b))
    b
  })

  beataml_for_drug <- reactive({
    if (identical(input$main_nav, "meta_aml4")) return(beataml_full())
    return(beataml_wave12())
  })

  drug_correlations <- reactive({
    input$main_nav
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(NULL)
    if (!exists("compute_drug_vaf_correlations")) return(NULL)
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    compute_drug_vaf_correlations(b, subset = subset)
  })

  drug_loo <- reactive({
    input$main_nav
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(NULL)
    if (!exists("compute_drug_vaf_loo")) return(NULL)
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    compute_drug_vaf_loo(b, subset = subset)
  })

  observe({
    b <- beataml_for_drug()
    corr <- drug_correlations()
    if (!is.null(corr) && nrow(corr) > 0) {
      genes <- sort(unique(corr$Gene))
      drugs <- sort(unique(corr$Inhibitor))
      updateSelectInput(session, "drug_gene", choices = c("Select..." = "", genes))
      updateSelectInput(session, "drug_inhibitor", choices = c("Select..." = "", drugs))
    }
  })

  output$drug_summary_ui <- renderUI({
    input$main_nav
    b <- beataml_for_drug()
    if (is.null(b)) return(p("BeatAML2 data not loaded. Run setup_beataml2.R to download."))
    if (!b$ok) return(p(style = "color:red", b$msg))
    corr <- drug_correlations()
    n_sig_p <- if (!is.null(corr) && nrow(corr) > 0) sum(corr$p_value < 0.05) else 0
    n_sig_q <- if (!is.null(corr) && nrow(corr) > 0) sum(corr$q_value < 0.1) else 0
    wave_note <- if (identical(input$main_nav, "analyses")) " (waves 1+2)" else " (all waves)"
    p(
      "Samples with mutation + drug data: ", length(b$overlap_samples), wave_note, " | ",
      "Drug-gene pairs analyzed: ", if (!is.null(corr)) nrow(corr) else 0, " | ",
      "Significant (p < 0.05): ", n_sig_p, " | ",
      "Significant (FDR q < 0.1): ", n_sig_q
    )
  })

  output$drug_summary_dotplot <- renderPlot({
    corr <- drug_correlations()
    if (is.null(corr) || nrow(corr) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No correlations"))
    # Paper-style: R² >= 0.25, VAF_range >= 25, AUC_range >= 75
    plot_df <- corr[corr$R_squared >= 0.25 & corr$VAF_range >= 25 & corr$AUC_range >= 75, , drop = FALSE]
    if (nrow(plot_df) == 0) {
      plot_df <- corr[corr$p_value < 0.05, , drop = FALSE]
      if (nrow(plot_df) == 0) plot_df <- head(corr[order(-corr$R_squared), ], 30)
    }
    if (nrow(plot_df) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No correlations meeting criteria"))
    plot_df$star <- ifelse(plot_df$q_value < 0.1, "*", "")
    gene_meds <- aggregate(delta_AUC ~ Gene, data = plot_df, FUN = median)
    plot_df$Gene <- factor(plot_df$Gene, levels = rev(gene_meds$Gene[order(gene_meds$delta_AUC)]))
    inh_meds <- aggregate(delta_AUC ~ Inhibitor, data = plot_df, FUN = median)
    plot_df$Inhibitor <- factor(plot_df$Inhibitor, levels = inh_meds$Inhibitor[order(inh_meds$delta_AUC)])
    ggplot(plot_df, aes(x = Inhibitor, y = Gene, fill = delta_AUC, size = VAF_range, label = star)) +
      geom_point(shape = 21, color = "black", stroke = 0.5) +
      geom_text(size = 5, color = "#525252", hjust = -0.2, vjust = 0.5, show.legend = FALSE) +
      scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0,
        name = expression(Delta~"AUC")) +
      scale_size_continuous(name = expression(Delta~"VAF"), range = c(2, 8)) +
      labs(x = "Inhibitor", y = NULL) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 15), legend.text = element_text(size = 13), legend.title = element_text(size = 14),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "#f0f0f0"),
        legend.position = "right")
  })

  output$drug_scatter_plot <- renderPlot({
    drug <- input$drug_inhibitor
    gene <- input$drug_gene
    if (is.null(drug) || drug == "" || is.null(gene) || gene == "") {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Select inhibitor and gene"))
    }
    b <- beataml_for_drug()
    if (is.null(b) || !b$ok) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data"))
    subset <- if (is.null(input$drug_subset)) "de_novo" else input$drug_subset
    allowed <- if (exists("get_beataml_allowed_samples")) get_beataml_allowed_samples(b, subset) else b$overlap_samples
    auc_sub <- b$auc[b$auc$inhibitor == drug & b$auc$Sample %in% allowed, c("Sample", "auc"), drop = FALSE]
    mut_sub <- b$mutations[b$mutations$Gene == gene & b$mutations$Sample %in% allowed, c("Sample", "VAF"), drop = FALSE]
    merged <- merge(auc_sub, mut_sub, by = "Sample")
    if (nrow(merged) < 5) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient samples"))
    fit <- lm(auc ~ VAF, data = merged)
    r2 <- round(summary(fit)$r.squared, 3)
    pv <- summary(fit)$coefficients[2, 4]
    ptxt <- if (pv < 0.001) "p < 0.001" else paste0("p = ", format.pval(pv, digits = 2))
    ggplot(merged, aes(x = VAF, y = auc)) +
      geom_point(size = 3, alpha = 0.8, color = "#3C5488") +
      geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.2) +
      labs(title = paste0(drug, " vs ", gene, " VAF"), x = "VAF (%)", y = "Drug AUC",
        subtitle = paste0("R² = ", r2, ", ", ptxt)) +
      theme_minimal(base_size = 16) +
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 15), plot.title = element_text(size = 16), plot.subtitle = element_text(size = 13))
  })

  output$drug_correlation_table <- renderDT({
    corr <- drug_correlations()
    if (is.null(corr) || nrow(corr) == 0) return(NULL)
    cols <- c("Inhibitor", "Gene", "Direction", "R_squared", "p_value", "q_value", "LOOCV_RMSE", "LOOCV_MSE", "LOOCV_MSE_sd", "n", "VAF_range", "AUC_range")
    cols <- cols[cols %in% colnames(corr)]
    df <- corr[order(corr$p_value), cols, drop = FALSE]
    df$p_value <- format.pval(df$p_value, digits = 2)
    df$q_value <- format.pval(df$q_value, digits = 2)
    if ("LOOCV_RMSE" %in% colnames(df)) df$LOOCV_RMSE <- round(df$LOOCV_RMSE, 2)
    if ("LOOCV_MSE" %in% colnames(df)) df$LOOCV_MSE <- round(df$LOOCV_MSE, 4)
    if ("LOOCV_MSE_sd" %in% colnames(df)) df$LOOCV_MSE_sd <- round(df$LOOCV_MSE_sd, 4)
    datatable(df, filter = "top", options = list(pageLength = 15))
  })

  output$drug_loo_heatmap <- renderPlot({
    loo <- drug_loo()
    if (is.null(loo) || nrow(loo) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No LOOCV data") + theme_void())
    loo <- loo[!is.na(loo$RMSE), , drop = FALSE]
    if (nrow(loo) == 0) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No RMSE values") + theme_void())
    gene_meds <- aggregate(RMSE ~ Gene, data = loo, FUN = median)
    loo$Gene <- factor(loo$Gene, levels = rev(gene_meds$Gene[order(gene_meds$RMSE)]))
    inh_meds <- aggregate(RMSE ~ Inhibitor, data = loo, FUN = median)
    loo$Inhibitor <- factor(loo$Inhibitor, levels = inh_meds$Inhibitor[order(inh_meds$RMSE)])
    ggplot(loo, aes(x = Inhibitor, y = Gene, fill = RMSE)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_viridis_c(option = "viridis", name = "RMSE") +
      labs(x = "Inhibitor", y = "Gene") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text = element_text(size = 11), legend.title = element_text(size = 12))
  })

  output$data_table <- renderDT({
    df <- filtered_data()
    cols <- c("Sample", "Gene", "VAF", "variant_type", "Subset", "Cohort", "Age", "Sex", "Risk", "Time_to_OS", "Censor", "mutation_category")
    cols <- cols[cols %in% colnames(df)]
    datatable(df[, cols, drop = FALSE], filter = "top", options = list(pageLength = 25, scrollX = TRUE))
  })
}

shinyApp(ui = ui, server = server)
