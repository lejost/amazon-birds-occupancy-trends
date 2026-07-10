# Script to compute the ratio of species with decreasing over increasing linear
# trends for parameters phi, gamma, psi and pp in old growth (OG) and secondary
# forest (SF), according to BirdNET and PROTAX-Sound outputs.
#
# NOTE: this replaces an earlier version that called OddsR1() and returned
# (n_neg / n_pos)^2. Because the two categories (increasing vs. decreasing)
# are complementary within the same sample of species, that formula collapsed
# into the squared ratio and did not correspond to a real odds ratio. Here we
# report the simple, directly interpretable count ratio n_neg / n_pos ("there
# are X times more species declining than increasing").
#
# Two kinds of tables, following the format in 'BirdNetOccDyn_temp.docx':
#     - Absolute change: all species
#     - Extreme change: only species with one or two asterisks (* or **)

# Clean workspace
rm(list=ls())

# Simple count ratio: n_neg / n_pos, returning NA when n_pos == 0
count_ratio <- function(neg, pos){
  ifelse(pos == 0 | is.na(pos) | is.na(neg), NA_real_, neg / pos)
}

# Load summary tables built on script 9
BN_tb <- readRDS("Outputs/TabTrendsBirdNET.rds")  # BirdNET table
MP_tb <- readRDS("Outputs/TabTrendsProtaxSound.rds") # Montpellier table

# --- BirdNET ---
# Absolute changes
BN_a <- count_ratio(BN_tb$Total_neg, BN_tb$Total_pos)

# Extreme changes
BN_e <- count_ratio((BN_tb$neg_one_ast + BN_tb$neg_two_ast),
                    (BN_tb$pos_one_ast + BN_tb$pos_two_ast))

# --- ProtaxSound ---
# Absolute changes
MP_a <- count_ratio(MP_tb$Total_neg, MP_tb$Total_pos)

# Extreme changes
MP_e <- count_ratio((MP_tb$neg_one_ast + MP_tb$neg_two_ast),
                    (MP_tb$pos_one_ast + MP_tb$pos_two_ast))


# Count-ratio df
count_ratio_df <- data.frame(BN_tb$Parameter,
                             BN_tb$Forest_type,
                             BN_a, BN_e,
                             MP_a, MP_e)
colnames(count_ratio_df) <- c("Parameter", "Forest type", "Absolute change - BN",
                             "Extreme change - BN", "Absolute change - MP",
                             "Extreme change - MP")

##

tabela <- data.frame(
  Parameter = count_ratio_df$Parameter[c(7,3,1,5)],
  MP_OG_abs = count_ratio_df$`Absolute change - MP`[c(7,3,1,5)],
  MP_SF_abs = count_ratio_df$`Absolute change - MP`[c(8,4,2,6)],
  BN_OG_abs = count_ratio_df$`Absolute change - BN`[c(7,3,1,5)],
  BN_SF_abs = count_ratio_df$`Absolute change - BN`[c(8,4,2,6)],
  MP_OG_ext = count_ratio_df$`Extreme change - MP`[c(7,3,1,5)],
  MP_SF_ext = count_ratio_df$`Extreme change - MP`[c(8,4,2,6)],
  BN_OG_ext = count_ratio_df$`Extreme change - BN`[c(7,3,1,5)],
  BN_SF_ext = count_ratio_df$`Extreme change - BN`[c(8,4,2,6)])
