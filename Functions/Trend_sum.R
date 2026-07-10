# --------------------------------------
# FUNCTION trend_sum
# required packages: base R (sub, gsub, nchar)
# description: summarizes the frequency of positive and negative parameter 
#              estimates across Old Growth (OG) and Second Forest (SF) sites. 
#              Counts occurrences of significant results (marked by * or **) 
#              and non-significant results regarding florest type.
# inputs: trend_table (table from trend tables Markdown output), 
#         param (string for target parameter - psi, phi, pp or gamma)
# outputs: data frame containing the frequency of positive/negative estimates 
#          for three levels (0, 1, or 2 asterisks) per forest type.
########################################

trend_sum <- function(trend_table, param){
  
  OG_col <- paste(param, "OG")
  SF_col <- paste(param, "SF")
  
  OG_mean <- as.numeric(sub("/.*", "", trend_table[[OG_col]])) # [[]] bc we're not working with vectors
  SF_mean <- as.numeric(sub("/.*", "", trend_table[[SF_col]]))
  
  OG_ast <- nchar(trend_table[[OG_col]]) - nchar(gsub("\\*", "", trend_table[[OG_col]])) # get number of * and **
  SF_ast <- nchar(trend_table[[SF_col]]) - nchar(gsub("\\*", "", trend_table[[SF_col]]))
  
  data.frame(
    Forest_type = c("OG","SF"),
    Parameter = param,
    
    # One asterisk
    pos_one_ast = c(
      sum(OG_mean >= 0 & OG_ast == 1),
      sum(SF_mean >= 0 & SF_ast == 1)
    ),
    
    neg_one_ast = c(
      sum(OG_mean < 0 & OG_ast == 1),
      sum(SF_mean < 0 & SF_ast == 1)
    ),
    
    # Two asterisks
    pos_two_ast = c(
      sum(OG_mean >= 0 & OG_ast == 2),
      sum(SF_mean >= 0 & SF_ast == 2)
    ),
    
    neg_two_ast = c(
      sum(OG_mean < 0 & OG_ast == 2),
      sum(SF_mean < 0 & SF_ast == 2)
    ),
    
    # Without asterisks
    pos_no_ast = c(
      sum(OG_mean >= 0 & OG_ast == 0),
      sum(SF_mean >= 0 & SF_ast == 0)
    ),
    
    neg_no_ast = c(
      sum(OG_mean < 0 & OG_ast == 0),
      sum(SF_mean < 0 & SF_ast == 0)
    )
  )
}