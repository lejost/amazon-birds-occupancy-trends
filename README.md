# BirdNetOccDyn 🐦

This project is composed by R scripts, functions and R Markdown / Quarto documents for the analysis of Amazon bird population dynamics, from raw BirdNET vocalization detections all the way to publication-ready trend tables and figures.

As an undergraduate researcher (Statistics) at the [Ferraz Population Biology Lab (UFRGS)](http://ferrazlab.org/), I work on the data-preparation, analysis, interpretation and reporting stages of the project. The core Bayesian hierarchical occupancy models (scripts `4` to `7` in the `Scripts` folder) were written by Dr. Gonçalo Ferraz — **everything else in this repository is my own work, developed under Ferraz's guidance.**

The goal of the project is to describe how Amazon bird populations respond to the impacts of forest degradation over
time, comparing old-growth (OG) and secondary forest (SF) — where old-growth is forest that has always retained its
original structure, and secondary forest is one that has regrown after being cleared. This is done using the
outputs of a Bayesian hierarchical multi-species dynamic occupancy model. The model estimates, for each of the
surveyed bird species, the probabilities of initial presence (ψ₁), persistence (φ), colonization (γ) and detection
(pp), and lets us track how these change across five years (2010–2014). Two automated audio-processing engines are
compared throughout: BirdNET (deep-learning based) and PROTAX-Sound.

## Workflow

The scripts are numbered to enforce the execution order. The pipeline goes from raw detections to results in three stages:

**Stage 1 — Data preparation (Jost):** read and clean the raw BirdNET detections, validate and cluster the recording sites, standardize site/date/species names against the official South American bird taxonomy, and organize everything into the detection and effort arrays that the model needs.

**Stage 2 — Hierarchical model (Jost, Ferraz)** format the arrays as binary detection data, specify and run the JAGS multi-species dynamic occupancy model, and check its goodness of fit.

**Stage 3 — Analysis & reporting (Jost):** take the model's posterior chains and turn them into interpretable results — parameter estimates, temporal-trend tables, and species-richness figures comparing forest types and engines.

## Result

- A reproducible pipeline that starts from raw BirdNET output and produces the tables and figures used in the manuscript, with intermediate objects (`.rds` / `.RData`) saved at each step so any stage can be re-run on its own.
- Species-richness figures and temporal-trend tables comparing OG vs SF forests across two audio-processing engines (BirdNET and PROTAX-Sound).

## Future updates

- Finalize the Quarto table formatting to fully match the target journal template.

## Communication & dissemination

- Poster accepted at **SINAPE** (26th Brazilian Symposium on Probability and Statistics) and at **ISEC** (International Statistical Ecology Conference, 2027).
- Abstracts written for both statistics audiences (SINAPE, ISEC) and a general scientific-initiation audience (SIC), adapting the biology/statistics balance to each.

---

## Files by folder

### Scripts
The numbered analysis workflow. Scripts **4–7 are the Ferraz et al's** hierarchical model; all others are mine.

- **`0_ProjectSetupScript.R`**: sets up the project infrastructure — loads the `barracudar` utility functions, initializes the random seed and the logging system.
- **`1_ReadBirdNETOutput.R`**: reads every BirdNET detection `.txt` file and aggregates them into a single RDS file (audio file name, detection begin/end time, species, confidence score). The detections come from recordings made between 2010 and 2014, which were processed by BirdNET between April and July 2025.
- **`2_CheckSitesDates.R`**: validates the detection sites against the master site list, clusters sites by distance (110 m threshold, hierarchical clustering) and standardizes site and date names.
- **`3_DataPreparation.R`**: transforms the cleaned detections into standardized 3D detection-count arrays (site group × date × species) and effort arrays (recording minutes, number of sites), reconciling species names with the SACC taxonomy.
- **`4_FormatBinaryData.R`** *(Ferraz)*: converts detection counts into the 4D binary presence/absence array used by the model.
- **`5_WriteRunJAGSModel.R`** *(Ferraz)*: specifies and runs the JAGS multi-species dynamic occupancy model (forest type × year interactions) via MCMC.
- **`6_MakeGOFPlots.R`** / **`7_MakeGOFPlots.R`** *(Ferraz)*: posterior predictive checks and goodness-of-fit diagnostics.
- **`6_PrepSuppFig1.R`**: prepares the dynamic-parameter posteriors (φ, γ, pp, ψ) — community means and species-specific estimates with 95% credible intervals — for the supplementary figures.
- **`8_MakeTrendTables.R`**: fits one linear regression per MCMC iteration to estimate the temporal trend (slope) of each parameter, and marks significance with `*` (50% credible interval excludes zero) or `**` (95% excludes zero).
- **`9_Tab_trends_total.R`**: summarizes how many species show increasing vs. decreasing trends (with/without asterisks) per parameter and forest type, for both engines.
- **`10_PlotSumOfPsis.R`**: computes posterior species richness (sum of occupancy probabilities) per forest type per year with 95% credible intervals, and produces the comparison figure (BirdNET vs. PROTAX-Sound).
- **`11_OddsRatio.R`**: computes the ratio of decreasing to increasing species for each parameter, forest type and engine, in absolute and "extreme" (significant-only) versions.

### Functions
Helper functions used across the analysis scripts.  The `GetPsis` function and its variants were mostly written by Ferraz.

- **`CleanSiteName.R`**: standardizes site-name capitalization.
- **`NameFix.R`**: reconciles species names between the BirdNET detections and the official SACC taxonomy.
- **`GetPsiMat.R`**: extracts the initial-occupancy (ψ₁) posterior draws as a species × iteration matrix.
- **`GetPsis.R`** / **`GetPsisArray.R`**: extract time-varying equilibrium occupancy (ψ*) posteriors, as summary matrices or as full posterior arrays for summation.
- **`GetMeanPsis.R`**: computes time-averaged equilibrium occupancy pooled over years.
- **`Trend_sum.R`**: tabulates the frequency of increasing/decreasing trends by significance level and forest type.

### Markdown
R Markdown / Quarto documents that render the results into Word/PDF tables and figures:

- **`trend_tables_phi_gamma.Rmd`** / **`trend_tables_psi_pp.Rmd`**: generate the formatted trend tables (persistence & colonization; occupancy & detection).
- **`trend_plots.Rmd`**: generates the trend plots.
- **`OddsTables.qmd`**: Quarto document that builds the effect-measure and species-count tables for the manuscript.

### OriginalData
Raw, never-edited source data:

- **`BirdNETOutput/`**: BirdNET selection tables organized by year (2010–2014).
- **`ArquivosDeGravacao.csv`**: metadata for the audio recordings (file paths, volumes, sizes).
- **`locations_diurnal.txt`**: site metadata — name, forest type (old growth = 0, secondary forest = 1) and geographic coordinates.
- **`SACCList2025.csv`**: the South American Classification Committee bird taxonomy reference.
- **`CheckMissing.xlsx`**: flags for corrupted audio files.
- **`trend_tables.rds`**: comparable trend results from the PROTAX-Sound analysis.

### CleanedData
Processed and deduplicated data:

- **`BNetProcOutput.rds`**: all cleaned BirdNET detections.
- **`BNDetDatesSpotGroupFile.rds`**: detections joined with dates, sites and site groups.
- **`ValidAudioFiles.rds`**: validated audio files with their durations.
- **`groups.csv`**: site groupings used for clustering.

### DataObjects
Intermediate serialized objects:

- **`DataArrays.RData`**: detection-count and effort arrays (output of script 3).
- **`BinaryDataArray.RData`** / **`y.rds`**: the 4D binary detection array (site × survey × year × species) used by the model.
- **`TrendTables09.rds`** / **`TrendTables085.rds`**: trend-slope estimates at the 0.9 and 0.85 confidence thresholds.

### Outputs
Generated results, tables and figures:

- **`ManausMSOD_MCMCchains.rds`**: the complete JAGS model output (posterior chains).
- **`DynParEstimates.rds`**: posterior parameter estimates (mean and 95% credible intervals) per species and forest type.
- **`TabTrendsBirdNET.rds`** / **`TabTrendsProtaxSound.rds`**: compact trend-count summaries feeding the effect-measure tables.
- **`trend_table_phi_gamma.rds`** / **`trend_table_psi_pp.rds`** / **`trend_tables_montp.rds`**: rendered trend tables.
- **`species_richness_plot.png`**: species-richness figure comparing OG vs SF across both engines.

### Utils
Utility data objects:

- **`compdata.rds`**: species names in several formats, for compatibility with the JAGS output.
- **`spsnames.rds`**: vector of scientific species names.

### barracudar
An external utility library (project-infrastructure, plotting and logging helpers) used by the setup script.
