# This script prepares a data object that will be used for 
# plotting n-species average and species specific values of
# persistence, colonization, detection, psi, equilibrium psi, 
# and turnover through time for old growth (OG) and secondary
# forest (SF).
#
# Adapted from Fábio's 'plots_dyn_param5.R'

# Load data and functions
output <- readRDS("Outputs/ManausMSOD_MCMCchains.rds")
compdata <- readRDS("Utils/compdata.rds")
source("Functions/GetDynPar.R")
source("Functions/GetPsi1.R")
source("Functions/GetPsit.R")
source("Functions/GetPsis.R")
source("Functions/GetMeanPsis.R")
source("Functions/GetTurnovert.R")
source("Functions/Getppt.R")

# Compute basic parameters from JAGS output
n.iter <- length(output$sims.list$beta.t)
n.year <- length(output$mean$mu.alpha.lpp.year)
n.sps <- length(output$mean$alpha.lpsi1)

# Prepare output objects ----
# list parameter names and stor number of parameters
parnames <- c("phi1","phi2","phi3","phi4",        # phi
              "gam1","gam2","gam3","gam4",        # gamma
              "det1","det2","det3","det4","det5", # pp
              "pss1","pss2","pss3","pss4",        # psi*
              "psi1","psi2","psi3","psi4","psi5", # mean psi
              "mpss", # year mean equilibrium psi
              "tur1","tur2","tur3","tur4") # turnover
# get number of estimates
nests <- length(parnames)
# store vector of names + Nsps mean
spsnames <- c("Nsps mean",compdata$SACCspacenames)

# Old Growth Estimates for plotting 
# mn = mean, lb = lower bound 95%, ub = upper bound 95%
# Extra row is for n-sps means
OGmn <- matrix(data=NA,nrow=n.sps+1,ncol=nests,dimnames = list(spsnames,parnames))
SFmn <- matrix(data=NA,nrow=n.sps+1,ncol=nests,dimnames = list(spsnames,parnames))
OGlb <- matrix(data=NA,nrow=n.sps+1,ncol=nests,dimnames = list(spsnames,parnames))
SFlb <- matrix(data=NA,nrow=n.sps+1,ncol=nests,dimnames = list(spsnames,parnames))
OGub <- matrix(data=NA,nrow=n.sps+1,ncol=nests,dimnames = list(spsnames,parnames))
SFub <- matrix(data=NA,nrow=n.sps+1,ncol=nests,dimnames = list(spsnames,parnames))

# Random effect means of all species ----
## Phi_t ----
### OG ----
# invoque get_dynpar to find posterior information
Nsps.pos <- get_dynpar("phi","OG")
OGmn[1,1:4] <- Nsps.pos[,1]
OGlb[1,1:4] <- Nsps.pos[,2]
OGub[1,1:4] <- Nsps.pos[,3]
### SF ----
# invoque get_dynpar to find posterior information
Nsps.pos <- get_dynpar("phi","SF")
# store time varying mean phi posterior summary (mean and bounds)
SFmn[1,1:4] <- Nsps.pos[,1]
SFlb[1,1:4] <- Nsps.pos[,2]
SFub[1,1:4] <- Nsps.pos[,3]
## Gamma_t ----
### OG ----
# invoque get_dynpar to find posterior information
Nsps.pos <- get_dynpar("gamma","OG")
# store time varying mean phi posterior summary (mean and bounds)
OGmn[1,5:8] <- Nsps.pos[,1]
OGlb[1,5:8] <- Nsps.pos[,2]
OGub[1,5:8] <- Nsps.pos[,3]
### SF ----
# invoque get_dynpar to find posterior information
Nsps.pos <- get_dynpar("gamma","SF")
# store time varying mean phi posterior summary (mean and bounds)
SFmn[1,5:8] <- Nsps.pos[,1]
SFlb[1,5:8] <- Nsps.pos[,2]
SFub[1,5:8] <- Nsps.pos[,3]
## Psi_1 ----
### OG ----
# invoque get_psi1 to obtain Psi1 posterior info
Nsps.pos <- get_psi1("OG")
# store result
OGmn[1,18] <- Nsps.pos[1]
OGlb[1,18] <- Nsps.pos[2]
OGub[1,18] <- Nsps.pos[3]
### SF ----
# invoque get_psi1 to obtain Psi1 posterior info
Nsps.pos <- get_psi1("SF")
# store result
SFmn[1,18] <- Nsps.pos[1]
SFlb[1,18] <- Nsps.pos[2]
SFub[1,18] <- Nsps.pos[3]
## Psi_t ----
### OG ----
# invoque get_psit to obtain Psi_t posterior info
Nsps.pos <- get_psit("OG")
OGmn[1,19:22] <- Nsps.pos[,1]
OGlb[1,19:22] <- Nsps.pos[,2]
OGub[1,19:22] <- Nsps.pos[,3]
### SF ----
# invoque get_psit to obtain Psi_t posterior info
Nsps.pos <- get_psit("SF")
SFmn[1,19:22] <- Nsps.pos[,1]
SFlb[1,19:22] <- Nsps.pos[,2]
SFub[1,19:22] <- Nsps.pos[,3]
## Psi*_t ----
### OG ----
# invoque get_psis to obtain Psi*_t posterior info
Nsps.pos <- get_psis("OG")
OGmn[1,14:17] <- Nsps.pos[,1]
OGlb[1,14:17] <- Nsps.pos[,2]
OGub[1,14:17] <- Nsps.pos[,3]
### SF ----
# invoque get_psis to obtain Psi*_t posterior info
Nsps.pos <- get_psis("SF")
SFmn[1,14:17] <- Nsps.pos[,1]
SFlb[1,14:17] <- Nsps.pos[,2]
SFub[1,14:17] <- Nsps.pos[,3]
## Psi*, time-averaged ----
### OG ----
# invoque get_mean_psis to obtain time averaged Psi* posterior info
Nsps.pos <- get_mean_psis("OG")
OGmn[1,23] <- Nsps.pos[1]
OGlb[1,23] <- Nsps.pos[2]
OGub[1,23] <- Nsps.pos[3]
### SF ----
# invoque get_mean_psis to obtain time averaged Psi* posterior info
Nsps.pos <- get_mean_psis("SF")
SFmn[1,23] <- Nsps.pos[1]
SFlb[1,23] <- Nsps.pos[2]
SFub[1,23] <- Nsps.pos[3]
## Turnover_t ----
### OG ----
# invoque get_turnovert to obtain time-varying Turnover posterior info
Nsps.pos <- get_turnovert("OG")
OGmn[1,24:27] <- Nsps.pos[,1]
OGlb[1,24:27] <- Nsps.pos[,2]
OGub[1,24:27] <- Nsps.pos[,3]
### SF ----
# invoque get_turnovert to obtain time-varying Turnover posterior info
Nsps.pos <- get_turnovert("SF")
SFmn[1,24:27] <- Nsps.pos[,1]
SFlb[1,24:27] <- Nsps.pos[,2]
SFub[1,24:27] <- Nsps.pos[,3]
## pp_t ----
### OG ----
# invoque get_ppt to obtain time-varying pp posterior info
Nsps.pos <- get_ppt("OG")
OGmn[1,9:13] <- Nsps.pos[,1]
OGlb[1,9:13] <- Nsps.pos[,2]
OGub[1,9:13] <- Nsps.pos[,3]
### SF ----
# invoque get_ppt to obtain time-varying pp posterior info
Nsps.pos <- get_ppt("SF")
SFmn[1,9:13] <- Nsps.pos[,1]
SFlb[1,9:13] <- Nsps.pos[,2]
SFub[1,9:13] <- Nsps.pos[,3]

# Species-specific results ----
for(s in 1:n.sps) {
  
  # find the position of the current species in the output objects that 
  # give separate species-specific results. csn = current species number
  csn <- which(compdata$spacenames == spsnames[-1][s]) 

  ## Phi_t ----
  ### OG ----
  # invoque get_dynpar to find posterior information
  Csps.pos <- get_dynpar("phi","OG",sps=csn)
  OGmn[(csn+1),1:4] <- Csps.pos[,1]
  OGlb[(csn+1),1:4] <- Csps.pos[,2]
  OGub[(csn+1),1:4] <- Csps.pos[,3]
  ### SF ----
  # invoque get_dynpar to find posterior information
  Csps.pos <- get_dynpar("phi","SF",sps=csn)
  # store time varying mean phi posterior summary (mean and bounds)
  SFmn[(csn+1),1:4] <- Csps.pos[,1]
  SFlb[(csn+1),1:4] <- Csps.pos[,2]
  SFub[(csn+1),1:4] <- Csps.pos[,3]

  ## Gamma_t ----
  ### OG ----
  # invoque get_dynpar to find posterior information
  Csps.pos <- get_dynpar("gamma","OG",sps=csn)
  # store time varying mean phi posterior summary (mean and bounds)
  OGmn[(csn+1),5:8] <- Csps.pos[,1]
  OGlb[(csn+1),5:8] <- Csps.pos[,2]
  OGub[(csn+1),5:8] <- Csps.pos[,3]
  ### SF ----
  # invoque get_dynpar to find posterior information
  Csps.pos <- get_dynpar("gamma","SF",sps=csn)
  # store time-varying mean phi posterior summary (mean and bounds)
  SFmn[(csn+1),5:8] <- Csps.pos[,1]
  SFlb[(csn+1),5:8] <- Csps.pos[,2]
  SFub[(csn+1),5:8] <- Csps.pos[,3]

  ## Psi_1 ----
  ### OG ----
  # invoque get_psi1 to obtain Psi1 posterior info
  Csps.pos <- get_psi1("OG",sps=csn)
  # store result
  OGmn[(csn+1),18] <- Csps.pos[1]
  OGlb[(csn+1),18] <- Csps.pos[2]
  OGub[(csn+1),18] <- Csps.pos[3]
  ### SF ----
  # invoque get_psi1 to obtain Psi1 posterior info
  Csps.pos <- get_psi1("SF",sps=csn)
  # store result
  SFmn[(csn+1),18] <- Csps.pos[1]
  SFlb[(csn+1),18] <- Csps.pos[2]
  SFub[(csn+1),18] <- Csps.pos[3]

  ## Psi_t ----
  ### OG ----
  # invoque get_psit to obtain Psi_t posterior info
  Csps.pos <- get_psit("OG",sps=csn)
  # store result
  OGmn[(csn+1),19:22] <- Csps.pos[,1]
  OGlb[(csn+1),19:22] <- Csps.pos[,2]
  OGub[(csn+1),19:22] <- Csps.pos[,3]
  ### SF ----
  # invoque get_psit to obtain Psi_t posterior info
  Csps.pos <- get_psit("SF",sps=csn)
  # store result
  SFmn[(csn+1),19:22] <- Csps.pos[,1]
  SFlb[(csn+1),19:22] <- Csps.pos[,2]
  SFub[(csn+1),19:22] <- Csps.pos[,3]

  ## Psi*_t ----
  ### OG ----
  # invoque get_psis to obtain Psi*_t posterior info
  Csps.pos <- get_psis("OG",sps=csn)
  OGmn[(csn+1),14:17] <- Csps.pos[,1]
  OGlb[(csn+1),14:17] <- Csps.pos[,2]
  OGub[(csn+1),14:17] <- Csps.pos[,3]
  ### SF ----
  # invoque get_psis to obtain Psi*_t posterior info
  Csps.pos <- get_psis("SF",sps=csn)
  SFmn[(csn+1),14:17] <- Csps.pos[,1]
  SFlb[(csn+1),14:17] <- Csps.pos[,2]
  SFub[(csn+1),14:17] <- Csps.pos[,3]
  
  ## Psi*, time-averaged ----
  ### OG ----
  # invoque get_mean_psis to obtain time averaged Psi* posterior info
  Csps.pos <- get_mean_psis("OG",sps=csn)
  OGmn[(csn+1),23] <- Csps.pos[1]
  OGlb[(csn+1),23] <- Csps.pos[2]
  OGub[(csn+1),23] <- Csps.pos[3]
  ### SF ----
  # invoque get_mean_psis to obtain time averaged Psi* posterior info
  Csps.pos <- get_mean_psis("SF",sps=csn)
  SFmn[(csn+1),23] <- Csps.pos[1]
  SFlb[(csn+1),23] <- Csps.pos[2]
  SFub[(csn+1),23] <- Csps.pos[3]
  
  ## Turnover_t ----
  ### OG ----
  # invoque get_turnovert to obtain time-varying Turnover posterior info
  Csps.pos <- get_turnovert("OG",sps=csn)
  OGmn[(csn+1),24:27] <- Csps.pos[,1]
  OGlb[(csn+1),24:27] <- Csps.pos[,2]
  OGub[(csn+1),24:27] <- Csps.pos[,3]
  ### SF ----
  # invoque get_turnovert to obtain time-varying Turnover posterior info
  Csps.pos <- get_turnovert("SF",sps=csn)
  SFmn[(csn+1),24:27] <- Csps.pos[,1]
  SFlb[(csn+1),24:27] <- Csps.pos[,2]
  SFub[(csn+1),24:27] <- Csps.pos[,3]
  
  ## pp_t ----
  ### OG ----
  # invoque get_ppt to obtain time-varying pp posterior info
  Csps.pos <- get_ppt("OG",sps=csn)
  OGmn[(csn+1),9:13] <- Csps.pos[,1]
  OGlb[(csn+1),9:13] <- Csps.pos[,2]
  OGub[(csn+1),9:13] <- Csps.pos[,3]
  ### SF ----
  # invoque get_ppt to obtain time-varying pp posterior info
  Csps.pos <- get_ppt("SF",sps=csn)
  SFmn[(csn+1),9:13] <- Csps.pos[,1]
  SFlb[(csn+1),9:13] <- Csps.pos[,2]
  SFub[(csn+1),9:13] <- Csps.pos[,3]
  
}

# Save data.frame with content of matrices to the Outputs folder
sps <- rep(spsnames,each=54)
yr <- rep(rep(c(rep(seq(2010,2013),2),seq(2010,2014),seq(2010,2013),seq(2010,2014),
        0,seq(2010,2013)),2),184)
ftype <- rep(c(rep("OG",27),rep("SF",27)),184)
varb <- rep(rep(colnames(OGmn),2),184)
mn <- lb <- ub <- rep(NA,length(varb))
for(i in 1:184) {
  # set dataframe row indices 
  sttOG <- (i-1)*54 + 1
  endOG <- (i-1)*54 + 27
  sttSF <- (i-1)*54 + 28
  endSF <- (i-1)*54 + 54
  # fill up
  mn[sttOG:endOG] <- OGmn[i,]
  lb[sttOG:endOG] <- OGlb[i,]
  ub[sttOG:endOG] <- OGub[i,]
  mn[sttSF:endSF] <- SFmn[i,]
  lb[sttSF:endSF] <- SFlb[i,]
  ub[sttSF:endSF] <- SFub[i,]
}
# assemble dataframe and save
ParEsts <- data.frame(sps,yr,ftype,varb,mn,lb,ub)
colnames(ParEsts) <- c("Species","year","ForType","parameter",
                       "mean","lb95","ub95")
saveRDS(object=ParEsts,
        file="Outputs/DynParEstimates.rds") 
