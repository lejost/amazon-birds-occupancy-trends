# Script to build trend tables for BirdNET (BirdNetOccDyn) 
# and PROTAX-Sound (MontpellierPaper) parameter outputs.
# Requires result tables from the 'trend_tables' Markdown.

# Load helper functions
source("Functions/trend_sum.R")

# -------
# BirdNET
# -------
# Load data frames
tab_df <- readRDS("Outputs/trend_table_phi_gamma.rds")
tab_df <- tab_df[-1,] # remove the intercept
tab_psi_pp <- readRDS("Outputs/trend_table_psi_pp.rds")
tab_psi_pp <- tab_psi_pp[-1,]

# rename columns 
names(tab_psi_pp) <- c("Species group and name", "psi OG", "pp OG", "psi SF", "pp SF")
tab_df <- cbind(tab_df, tab_psi_pp[,-1])

# calculate frequency of * and ** for each parameter
phi_tab   <- trend_sum(tab_df, "phi")
gamma_tab <- trend_sum(tab_df, "gamma")
psi_tab   <- trend_sum(tab_df, "psi")
pp_tab    <- trend_sum(tab_df, "pp")

# combine results into a single data frame
summary_tab <- rbind(phi_tab, gamma_tab, psi_tab, pp_tab)

# get totals
summary_tab$Total_pos <- summary_tab$pos_one_ast +
  summary_tab$pos_two_ast +
  summary_tab$pos_no_ast

summary_tab$Total_neg <- summary_tab$neg_one_ast +
  summary_tab$neg_two_ast +
  summary_tab$neg_no_ast

summary_tab$Total <- summary_tab$Total_pos + summary_tab$Total_neg

# -------
# Montpellier
# -------
#tab_montp <- readRDS("/Users/gferr/OneDrive/Projects/NSpsOccDynManaus/MontpellierPaper/Outputs/trend_tables.rds")
tab_montp <- readRDS("OriginalData/trend_tables.rds")
names(tab_montp) <- c("Species group and name", 
                       "phi OG", "gamma OG", "pp OG", "psi OG",
                       "phi SF", "gamma SF", "pp SF", "psi SF")
tab_montp <- tab_montp[-1,] # remove intercept

# apply the same calculation process as above
phi_tab   <- trend_sum(tab_montp, "phi")
gamma_tab <- trend_sum(tab_montp, "gamma")
psi_tab   <- trend_sum(tab_montp, "psi")
pp_tab    <- trend_sum(tab_montp, "pp")

summary_tab_montp <- rbind(phi_tab, gamma_tab, psi_tab, pp_tab)

summary_tab_montp$Total_pos <- summary_tab_montp$pos_one_ast +
  summary_tab_montp$pos_two_ast +
  summary_tab_montp$pos_no_ast

summary_tab_montp$Total_neg <- summary_tab_montp$neg_one_ast +
  summary_tab_montp$neg_two_ast +
  summary_tab_montp$neg_no_ast

summary_tab_montp$Total <- summary_tab_montp$Total_pos + summary_tab_montp$Total_neg

