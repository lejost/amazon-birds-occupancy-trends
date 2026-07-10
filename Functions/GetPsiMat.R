# --------------------------------------
# FUNCTION get_psis_mat
# required packages: 
# description: 
# inputs: fortyp (forest type)
# outputs: sxi matrix of posterior draws where s is the number of species and i 
#          is the number of MCMC iterations stored in the data object. 
########################################

## Function
get_psi1_mat <- function(output, fortyp) {
  
  # Safeguard
  if(!fortyp %in% c("OG","SF")) return("Invalid forest type")
  
  # Get number of iterations and species
  n.iter <- dim(output$sims.list$alpha.lpsi1)[1]
  n.spec <- dim(output$sims.list$alpha.lpsi1)[2]
  
  # Set up storing matrix (species x iterations)
  psi1_pos <- matrix(NA, n.spec, n.iter)
  
  # Obtain posterior information
  if(fortyp == "OG") { # Old growth
    
    lpsi1 <- output$sims.list$alpha.lpsi1   # [iter x species]
    
  } else { # Secondary forest
    
    lpsi1 <- output$sims.list$alpha.lpsi1 + 
      output$sims.list$beta.lpsi1.sf  # [iter x species]
  }
  
  # go to probability space
  psi.tmp <- plogis(lpsi1)  # [iter x species]
  
  # organize dimensions (species x iterations)
  psi1_pos <- t(psi.tmp)
  
  return(psi1_pos)
}

# 
# psis_og <- output$sims.list$alpha.lpsi1
# all_sps_og <- rowSums(psis_og)
# 
# psis_sf <- output$sims.list$alpha.lpsi1 + output$sims.list$beta.lpsi1.sf
# all_sps_sf <- rowSums(psis_sf)
