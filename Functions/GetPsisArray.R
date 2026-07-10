# --------------------------------------
# FUNCTION get_psis_arr
# required packages: 
# description: Extracts psi from the estimates of phi and gamma of all 
#              species in the JAGS output object. Works with occ dyn pars
#              from years t=1,2,3,4, therefore producing values of psi for
#              years t=2,3,4,5. The psi matrix for year 1 is harvested by
#              a different function get_psi1_mat. Data object from the
#              ManausMSOD_MCMCchains.rds must be available in the workspace.
# inputs: fortyp (forest type)
# outputs: sx4xi array where s is the number of species and i is the number of 
#          MCMC iterations stored in the data object. 
########################################
get_psis_arr <- function(output, fortyp) {
  
  # Safeguard
  if(!fortyp %in% c("OG","SF")) return("Invalid forest type")
  
  # dimensions
  n.iter <- dim(output$sims.list$alpha.lphi.year)[1]
  n.years <- 4
  n.spec <- dim(output$sims.list$alpha.lphi.year)[3]
  
  # Set up array (species x years x iterations)
  psis <- array(NA, dim = c(n.spec, n.years, n.iter))
  
  # Obtain posterior information
  if(fortyp == "OG") { # Old growth
    
    epsog <- 1 - plogis(output$sims.list$alpha.lphi.year)
    gammaog <- plogis(output$sims.list$alpha.lgamma.year)
    
    psi.tmp <- gammaog / (gammaog + epsog)
    
  } else { # Secondary forest
    
    epssf <- 1 - plogis(output$sims.list$alpha.lphi.year +
                          array(rep(output$sims.list$beta.lphi.sf, n.years),
                                dim = c(n.iter, n.years, n.spec)) +
                          output$sims.list$beta.lphi.sfXyear)
    
    gammasf <- plogis(output$sims.list$alpha.lgamma.year +
                        array(rep(output$sims.list$beta.lgamma.sf, n.years),
                              dim = c(n.iter, n.years, n.spec)) +
                        output$sims.list$beta.lgamma.sfXyear)
    
    psi.tmp <- gammasf / (gammasf + epssf)
  }
  
  # Species x years x iterations
  psis <- aperm(psi.tmp, c(3,2,1))
  
  return(psis)
}