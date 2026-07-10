# --------------------------------------
# FUNCTION get_mean_psis
# required packages: 
# description: extracts estimate of time-averaged prediction of equilibrium 
#              psi, or psi* (giving one single posterior mean and 95% bounds) 
#              of either the species random effect mean or of specific species, 
#              in OG or SF sites, from a JAGS output of the script 
#              5_WriteRunJAGSModel.R. 
#              Objects output and compdata must be available in workspace.
# inputs: fortyp (forest type), spsn (species number or name, where 0 or "all" 
#         gives random effect means)
# outputs: 1x3 vector with posterior mean, 95% lower bound and 95% upper bound
#          of equilibrium psi predicted from the time averaged epsilon and 
#          gamma probabilities
########################################
get_mean_psis <- function(fortyp,sps=0) {
  
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
  
  # Dimmensions
  n.iter <- length(output$sims.list$mu.alpha.lpsi1)
  # Set up storing vector
  mpsis <- rep(NA,n.iter)
  
  # Obtain posterior information
  if(fortyp=="OG") { # in Old growth

    # get needed chains in probability space
    if(sps==0 | sps =="all") {
      lphiog <- output$sims.list$mu.alpha.lphi.year
      mulphiog <- rowSums(lphiog)/4
      muepsog <- 1 - plogis(mulphiog)
      lgammaog <- output$sims.list$mu.alpha.lgamma.year
      mulgammaog <- rowSums(lgammaog)/4
      mugammaog <- plogis(mulgammaog)
    } else {
      lphiog <- output$sims.list$alpha.lphi.year[,,spsn]
      mulphiog <- rowSums(lphiog)/4
      muepsog <- 1 - plogis(mulphiog)
      lgammaog <- output$sims.list$alpha.lgamma.year[,,spsn]
      mulgammaog <- rowSums(lgammaog)/4
      mugammaog <- plogis(mulgammaog)
    }

    # Compute time-averaged psi*
    mupsis <- mugammaog / (mugammaog + muepsog)
    # extract posterior summary
    mumupsis <- mean(mupsis)
    lbmupsis <- quantile(mupsis,probs=0.025) # lb = lower bound
    ubmupsis <- quantile(mupsis,probs=0.975) # ub = upper bound

  } else { # in Secondary forest

    # get needed chains in probability space
    if(sps==0 | sps =="all") {
      lphisf <- output$sims.list$mu.alpha.lphi.year +
                matrix(rep(output$sims.list$mu.beta.lphi.sf,4),ncol=4) +
                matrix(c(rep(0,n.iter),output$sims.list$mu.alpha.lphi.sfXyear),ncol=4)
      mulphisf <- rowSums(lphisf)/4
      muepssf <- 1 - plogis(mulphisf)
      lgammasf <- output$sims.list$mu.alpha.lgamma.year +
                  matrix(rep(output$sims.list$mu.beta.lgamma.sf,4),ncol=4) +
                  matrix(c(rep(0,n.iter),output$sims.list$mu.alpha.lgamma.sfXyear),ncol=4)
      mulgammasf <- rowSums(lgammasf)/4
      mugammasf <- plogis(mulgammasf)
    } else {
      lphisf <- output$sims.list$alpha.lphi.year[,,spsn] +
                matrix(rep(output$sims.list$beta.lphi.sf[,spsn],4),ncol=4) +
                output$sims.list$beta.lphi.sfXyear[,,spsn]
      mulphisf <- rowSums(lphisf)/4
      muepssf <- 1 - plogis(mulphisf)
      lgammasf <- output$sims.list$alpha.lgamma.year[,,spsn] +
                  matrix(rep(output$sims.list$beta.lgamma.sf[,spsn],4),ncol=4) +
                  output$sims.list$beta.lgamma.sfXyear[,,spsn]
      mulgammasf <- rowSums(lgammasf)/4
      mugammasf <- plogis(mulgammasf)
    }
    
    # Compute time-averaged psi*
    mupsis <- mugammasf / (mugammasf + muepssf)
    # extract posterior summary
    mumupsis <- mean(mupsis)
    lbmupsis <- quantile(mupsis,probs=0.025) # lb = lower bound
    ubmupsis <- quantile(mupsis,probs=0.975) # ub = upper bound
    
  }
  
  # package for delivery
  output <- cbind(mumupsis,lbmupsis,ubmupsis)
  colnames(output) <- c("mean","0.025%","0.975%")
  rownames(output) <- c("Time-averaged Psi*")
  return(output)
  
}
