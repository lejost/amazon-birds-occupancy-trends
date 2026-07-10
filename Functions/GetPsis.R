# --------------------------------------
# FUNCTION get_psis
# required packages: 
# description: extracts estimates of time varying predictions of equilibrium 
#              psi, or psi* (giving posterior mean and 95% bounds
#              for year t=1,2,3,4) of specific species OR random effect mean, 
#              in OG or SF sites, from a JAGS output of the script 
#              5_WriteRunJAGSModel.R. 
#              Objects output and compdata must be available in workspace.
# inputs: fortyp (forest type), spsn (species number or name, where 0 or "all" 
#         gives random effect means)
# outputs: 4x3 matrix with posterior mean, 95% lower bound and 95% upper bound
#          for the four years < 5
########################################
get_psis <- function(fortyp,sps=0) {
  
  # Safeguard
  if(!fortyp%in%c("OG","SF")) return("Invalid forest type")
  if(is.numeric(sps)) {
    if(sps<0 | sps>183) return("Species number out of bounds")
  } else {
    if(!sps%in%c("all",compdata$SACCspacenames)) return("Invalid species name")
  }
  
  # Turn species name into species number, if needed
  if(is.character(sps) & sps!="all"){
    spsn <- which(compdata$spacenames == sps)
  } else { # if not, keep the number
    spsn <- sps
  }
  
  # dimmensions
  n.iter <- length(output$sims.list$mu.alpha.lpsi1)
  n.years <- 4
  # Set up storing matrix
  psis <- matrix(NA,n.iter,n.years)
  
  # Obtain posterior information
  if(fortyp=="OG") { # in Old growth

    # get needed chains in probability space
    if(sps==0 | sps =="all") {
      epsog <- 1 - plogis(output$sims.list$mu.alpha.lphi.year)
      gammaog <- plogis(output$sims.list$mu.alpha.lgamma.year)
    } else {
      epsog <- 1 - plogis(output$sims.list$alpha.lphi.year[,,spsn])
      gammaog <- plogis(output$sims.list$alpha.lgamma.year[,,spsn])
    }
    # fill the four columns of psi* at once
    psis <- gammaog / (gammaog + epsog)
    # extract posterior summary
    mupsis <- apply(psis,2,mean)
    lbpsis <- apply(psis,2,quantile,probs=0.025) # lb = lower bound
    ubpsis <- apply(psis,2,quantile,probs=0.975) # ub = upper bound

  } else { # in Secondary forest

    # get needed chains in probability space
    if(sps==0 | sps =="all") { # for random effect means
      epssf <- 1 - plogis(output$sims.list$mu.alpha.lphi.year +
                          matrix(rep(output$sims.list$mu.beta.lphi.sf,4),ncol=4) +
                          matrix(c(rep(0,n.iter),output$sims.list$mu.alpha.lphi.sfXyear),ncol=4))
      gammasf <- plogis(output$sims.list$mu.alpha.lgamma.year +
                          matrix(rep(output$sims.list$mu.beta.lgamma.sf,4),ncol=4) +
                          matrix(c(rep(0,n.iter),output$sims.list$mu.alpha.lgamma.sfXyear),ncol=4))
    } else { # for individual species
      epssf <- 1 - plogis(output$sims.list$alpha.lphi.year[,,spsn] +
                          matrix(rep(output$sims.list$beta.lphi.sf[,spsn],4),ncol=4) +
                          output$sims.list$beta.lphi.sfXyear[,,spsn])
      gammasf <- plogis(output$sims.list$alpha.lgamma.year[,,spsn] +
                          matrix(rep(output$sims.list$beta.lgamma.sf[,spsn],4),ncol=4) +
                          output$sims.list$beta.lgamma.sfXyear[,,spsn])
    }
    # fill the four columns of psi* at once
    psis <- gammasf / (gammasf + epssf)
    # extract posterior summary
    mupsis <- apply(psis,2,mean)
    lbpsis <- apply(psis,2,quantile,probs=0.025) # lb = lower bound
    ubpsis <- apply(psis,2,quantile,probs=0.975) # ub = upper bound
    
  }
  
  # package for delivery
  output <- cbind(mupsis,lbpsis,ubpsis)
  colnames(output) <- c("mean","0.025%","0.975%")
  rownames(output) <- c("Psi*1","Psi*2","Psi*3","Psi*4")
  return(output)
  
}
