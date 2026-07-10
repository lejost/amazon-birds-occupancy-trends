## DCM_codigo 16
## Dynamic N-sps Model (DCM)
## Previous version of the file containing this scripe was named
## NSODMao_Code16_array - withPPC - v2 (clean).R
## Content renamed and saved with current file name '4_WriteRunJAGSModel.R'
## on December 14, 2023

## standard formulation where we have species-level random effects parameters 
## which are all drawn from community-level distributions with hyperparams that 
## characterize the community.

## model has the interaction between Forest type and Yearly interval 
## in the two dynamics rates

## in this model 16 (which includes interaction between forest type and time in
## pp, as does model 15) we build a quadratic effect of time on detection (pp).

# Model 16: with the traditional species.nested.within.community structure and 
#           interaction between Year and Forest Type
#
# ------------------------------------------------------------------------------
# This model now has the following linear model structures:
# - main effects of forest type (oldgrowth vs. secondgrowth) and time (years
# 1-5)
#        AND also the interaction terms between the two
# - these appear both as community effects (these are the hyperparameters)
#        and as species-level random effects
# - in this model 16 we try to apply interaction between forest type and time in
#   pp, like previously done for phi and gamma in model 14
# - quadratic effect of time on detection (pp)

# ==============================================================================
## Call libraries
library(jagsUI)
library(tictoc)

# ==============================================================================
# Load data input:

load("DataObjects/BinaryDataArray.RData")

# fill NA's in stjul with zeros - this will not affect the analysis because
# there are no samples corresponding to these days and the standardization was
# completed before replacing NA with zero
stjul0 <- stjul # standardized julian days of sampling
stjul0[which(is.na(stjul),arr.ind = TRUE)]<-0 # fill NAs with zero

# ==============================================================================
# Create "noccasion", a binary indicator for non-missingness
notNA <- y
notNA[notNA == 0] <- 1
notNA[is.na(notNA)] <- 0
sum(is.na(y))  ;  sum(is.na(notNA))   # Looks good
noccasion <- apply(notNA, c(1,3,4), sum, na.rm = TRUE) 
# Same under assumption that same non-missingness for all species
noccasion <- apply(notNA[,,,1], c(1,3), sum, na.rm = TRUE)
str(noccasion)

# ==============================================================================
# Prepare model input:

str(bdata <- list(y = y, nsite = dim(y)[1], nsurvey = dim(y)[2],
                  nyear = dim(y)[3], nspec = dim(y)[4], ss = sampsize,
                  sf = secfor, stjul = stjul0, noccasion = noccasion))

#List of 8
# $ y      : num [1:203301, 1:5] 0 0 0 0 0 0 0 0 0 1 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : chr [1:203301] "1" "2" "3" "4" ...
#  .. ..$ : chr [1:5] "y" "site" "visit" "year" ...
# $ nry    : int 203301
# $ nsite  : int 158
# $ nsurvey: int 13
# $ nyear  : int 5
# $ nspec  : int 63
# $ ss     : num [1:158, 1:13, 1:5] 1 1 1 1 1 1 1 1 1 1 ...
# $ sf     : int [1:158] 0 0 0 0 0 0 0 0 0 0 ...

## Specify model in BUGS language (as vector)
cat(file = "DCM16X.txt","
model {
  # Define models for parameters (link transformations, linear models, constraints)
  # -------------------------------------------------------------------------------
  # Initial occupancy probability psi1 (indexed by site and species)
  # Random species (intercept) effect 
  #     plus random, species-level effect of secondary forest
  for(i in 1:nsite){
    for(k in 1:nspec){
	    logit(psi1[i,k]) <- lpsi1[i,k]
	    lpsi1[i,k] <- alpha.lpsi1[k] + beta.lpsi1.sf[k] * sf[i]
	  } 
  }
  # Persistence and colonization probs. phi and gamma 
  #   (indexed by site, yearly interval and species)
  # This is now year + sf + year:sf for both phi and gamma
  for(i in 1:nsite){
    for(t in 1:(nyear-1)){
      for(k in 1:nspec){
        # Persistence
	      logit(phi[i,t,k]) <- lphi[i,t,k]
	      lphi[i,t,k] <- alpha.lphi.year[t,k] + 
	                     beta.lphi.sf[k] * sf[i] +
		                   beta.lphi.sfXyear[t,k] * sf[i] 
		                   # Last term is additional year effect only when sf == 1
		                   # and t > 1 
        # Colonization
	      logit(gamma[i,t,k]) <- lgamma[i,t,k]
	      lgamma[i,t,k] <- alpha.lgamma.year[t,k] + 
	                       beta.lgamma.sf[k] * sf[i] +
		                     beta.lgamma.sfXyear[t,k] * sf[i] 
		                     # Last term is additional year effect only when 
		                     # sf == 1 and t > 1
	    } 
    }
  }
  # Detection probability pp (indexed by site, year and species)
  # This is now year + sf + year:sf like for phi and gamma
  for(i in 1:nsite){
    for(t in 1:nyear){
      for(k in 1:nspec){
        logit(pp[i,t,k]) <- lpp[i,t,k]
	      lpp[i,t,k] <- alpha.lpp.year[t,k] + 
	                    beta.lpp.sf[k] * sf[i] +
		                  beta.lpp.sfXyear[t,k] * sf[i]  +
                      beta.t * stjul[i,t] +
                      beta.tt * (stjul[i,t]^2)
                      
		                  # 'beta.lpp.sfXyear' is additional year effect only when
		                  # sf == 1 and t > 1
		                  # 'beta.t' and 'beta.tt' are parameters of the quadratic 
		                  # effect of julian day on detection
	                    
        # Linear model for effort has year effect plus SF effect, 
        # both species-specific and without interaction
      }
    }
  }

  # Hyperpriors that define the relationship among species
  # ----------------------------------------------------------------------------------------

  # Regarding parameters in initial occupancy (psi1)
  for(k in 1:nspec){
    alpha.lpsi1[k] ~ dnorm(mu.alpha.lpsi1, tau.alpha.lpsi1)
    beta.lpsi1.sf[k] ~ dnorm(mu.beta.lpsi1.sf, tau.beta.lpsi1.sf)
  }
  mu.alpha.lpsi1 <- logit(mean.alpha.lpsi1)
  mean.alpha.lpsi1 ~ dunif(0, 1)
  tau.alpha.lpsi1 <- pow(sd.alpha.lpsi1, -2)
  sd.alpha.lpsi1 ~ dnt(0, 1/2.25^2, 1)I(0,)  # Half Cauchy with SD = 2.25 (as Broms et al. 2016)

  mu.beta.lpsi1.sf <- logit(mean.beta.lpsi1.sf)
  mean.beta.lpsi1.sf ~ dunif(0, 1)
  tau.beta.lpsi1.sf <- pow(sd.beta.lpsi1.sf, -2)
  sd.beta.lpsi1.sf ~ dnt(0, 1/2.25^2, 1)I(0,)# Half Cauchy with SD = 2.25 (as Broms et al. 2016)

  # Regarding parameters in the dynamics part and pp
  for(t in 1:(nyear-1)){
    for(k in 1:nspec){
      alpha.lphi.year[t,k] ~ dnorm(mu.alpha.lphi.year[t], tau.alpha.lphi.year[t])  
      alpha.lgamma.year[t,k] ~ dnorm(mu.alpha.lgamma.year[t], tau.alpha.lgamma.year[t])  
    }
    mu.alpha.lphi.year[t] <- logit(mean.alpha.lphi.year[t])
    mean.alpha.lphi.year[t] ~ dunif(0, 1)
    tau.alpha.lphi.year[t] <- pow(sd.alpha.lphi.year[t], -2)
	  sd.alpha.lphi.year[t] ~ dnt(0, 1/2.25^2, 1)I(0,)# Half Cauchy with SD = 2.25
   
    mu.alpha.lgamma.year[t] <- logit(mean.alpha.lgamma.year[t])
    mean.alpha.lgamma.year[t] ~ dunif(0, 1)
    tau.alpha.lgamma.year[t] <- pow(sd.alpha.lgamma.year[t], -2)
	  sd.alpha.lgamma.year[t] ~ dnt(0, 1/2.25^2, 1)I(0,) # Half Cauchy with SD = 2.25
  }
  for(k in 1:nspec){
    beta.lphi.sf[k] ~ dnorm(mu.beta.lphi.sf, tau.beta.lphi.sf)
    beta.lgamma.sf[k] ~ dnorm(mu.beta.lgamma.sf, tau.beta.lgamma.sf)
  }
  mu.beta.lphi.sf <- logit(mean.beta.lphi.sf)
  mean.beta.lphi.sf ~ dunif(0, 1)
  tau.beta.lphi.sf <- pow(sd.beta.lphi.sf, -2)
  sd.beta.lphi.sf ~ dnt(0, 1/2.25^2, 1)I(0,) # Half Cauchy with SD = 2.25
  mu.beta.lgamma.sf <- logit(mean.beta.lgamma.sf)
  mean.beta.lgamma.sf ~ dunif(0, 1)
  tau.beta.lgamma.sf <- pow(sd.beta.lgamma.sf, -2)
  sd.beta.lgamma.sf ~ dnt(0, 1/2.25^2, 1)I(0,)# Half Cauchy with SD = 2.25

  # Interaction terms in phi and gamma
  for(k in 1:nspec){           # Fix at zero first year to avoid overparameterization
	  beta.lphi.sfXyear[1,k] <- 0
	  beta.lgamma.sfXyear[1,k] <- 0
	  beta.lpp.sfXyear[1,k] <- 0
  }  
  for(t in 1:3){        # Need to start index at 1 or else get in trouble with hyperparams
    for(k in 1:nspec){
      beta.lphi.sfXyear[t+1,k] ~ dnorm(mu.alpha.lphi.sfXyear[t], tau.alpha.lphi.sfXyear[t])  
      beta.lgamma.sfXyear[t+1,k] ~ dnorm(mu.alpha.lgamma.sfXyear[t], tau.alpha.lgamma.sfXyear[t])  
    }
    mu.alpha.lphi.sfXyear[t] <- logit(mean.alpha.lphi.sfXyear[t])
    mean.alpha.lphi.sfXyear[t] ~ dunif(0, 1)
    tau.alpha.lphi.sfXyear[t] <- pow(sd.alpha.lphi.sfXyear[t], -2)
	  sd.alpha.lphi.sfXyear[t] ~ dnt(0, 1/2.25^2, 1)I(0,)# Half Cauchy with SD = 2.25
    mu.alpha.lgamma.sfXyear[t] <- logit(mean.alpha.lgamma.sfXyear[t])
    mean.alpha.lgamma.sfXyear[t] ~ dunif(0, 1)
    tau.alpha.lgamma.sfXyear[t] <- pow(sd.alpha.lgamma.sfXyear[t], -2)
	  sd.alpha.lgamma.sfXyear[t] ~ dnt(0, 1/2.25^2, 1)I(0,) # Half Cauchy with SD = 2.25
  }
  
  # Community models for parameters in the detection part
  # species-specific year effects plus species-specific forest effects 
  for(t in 1:nyear){
    for(k in 1:nspec){
      alpha.lpp.year[t,k] ~ dnorm(mu.alpha.lpp.year[t], tau.alpha.lpp.year[t])  
    }
    mu.alpha.lpp.year[t] <- logit(mean.alpha.lpp.year[t])
    mean.alpha.lpp.year[t] ~ dunif(0, 1)
    tau.alpha.lpp.year[t] <- pow(sd.alpha.lpp.year[t], -2)
	  sd.alpha.lpp.year[t] ~ dnt(0, 1/2.25^2, 1)I(0,) # Half Cauchy with SD = 2.25
  }
  for(k in 1:nspec){
    beta.lpp.sf[k] ~ dnorm(mu.beta.lpp.sf, tau.beta.lpp.sf)
  }
  mu.beta.lpp.sf <- logit(mean.beta.lpp.sf)
  mean.beta.lpp.sf ~ dunif(0, 1)
  tau.beta.lpp.sf <- pow(sd.beta.lpp.sf, -2)
  sd.beta.lpp.sf ~ dnt(0, 1/2.25^2, 1)I(0,)   # Half Cauchy with SD = 2.25
  
  for(t in 1:4){        # Need to start index at 1 or else get in trouble with hyperparams
    for(k in 1:nspec){
      beta.lpp.sfXyear[t+1,k] ~ dnorm(mu.alpha.lpp.sfXyear[t], tau.alpha.lpp.sfXyear[t])  
    }
    mu.alpha.lpp.sfXyear[t] <- logit(mean.alpha.lpp.sfXyear[t])
    mean.alpha.lpp.sfXyear[t] ~ dunif(0, 1)
    tau.alpha.lpp.sfXyear[t] <- pow(sd.alpha.lpp.sfXyear[t], -2)
	  sd.alpha.lpp.sfXyear[t] ~ dnt(0, 1/2.25^2, 1)I(0,)# Half Cauchy with SD = 2.25
  }
  
  # quadratic effect of julian day
  beta.t ~ dnorm(0,0.1)
  beta.tt ~ dnorm(0,0.1)
  
  # 'Likelihood' (or main structure of the model)
  # ---------------------------------------------
  # Ecological submodel: Define state conditional on parameters
  for (i in 1:nsite){             # Loop over sites
    for(k in 1:nspec){            # Loop over species
      # Initial conditions of system
      z[i,1,k] ~ dbern(psi1[i,k]) # Presence/absence at start of study
      # State transitions
      for (t in 2:nyear){         # Loop over yearly intervals
        z[i,t,k] ~ dbern(z[i,t-1,k] * phi[i,t-1,k] + (1-z[i,t-1,k]) * gamma[i,t-1,k])
      }
    }
  }

  # Observation model (now in array form)
  for (i in 1:nsite){ # Loop over sites
    for(k in 1:nspec){ # Loop over species
      for (j in 1:nsurvey){ # Loop over surveys
        for (t in 1:nyear){ # Loop over years
          pg[i,j,t,k] <- 1 - (1-pp[i,t,k])^ss[i,j,t] # ss e sample size dado
          # notar variacao temporal no pp (acima)
          y[i,j,t,k] ~ dbern(z[i,t,k] * pg[i,j,t,k])
        }
      }
    }
  }
  
  # GoF part of the model
  # --------------------
  # Aggregate for observed data and draw 'replicate data' under the same model
  # and then aggregate these in the same way. In this way we then can compare
  # the observed frequencies with the frequencies that we would expect for 
  # a fitting model. For storage and computing saving here we will compute only
  # 'dfrep1' and 'eval' objects, which have to be generated by JAGS. The
  # remaining objects can be computed based on them outside of JAGS.
  
  # (Note we change order of loops)
  for(k in 1:nspec){ # Loop over species
    for (t in 1:nyear){ # Loop over years
      for (i in 1:nsite){ # Loop over sites
        for (j in 1:nsurvey){ # Loop over surveys
          # Create a new data set at each iteration. This data set is 'perfect' since
		      # we simulate it under the exact structural assumptions of our model, plus
		      # we use the exact values of the parameters right at each iteration of the 
		      # MCMC algorithm
          yrep[i,j,t,k] ~ dbern(z[i,t,k] * pg[i,j,t,k]) # pg is zero when there is no effort
        }
        # Now we aggregate the replicate data, yrep
        # Notice that we will compute the observed data, y, outside of JAGS
        dfrep1[i,t,k] <- sum(yrep[i,,t,k]) # Same for 'ideal' data

        # Expected value at average p
		    eval[i,t,k] <- z[i,t,k] * noccasion[i,t] * mean(pg[i,,t,k])

      }
	  }
  }

}
")

# Initial values
zst <- apply(y,c(1,3,4),max,na.rm=TRUE) # Observed occurrence as inits for z
zst[zst == '-Inf'] <- 0
inits <- function(){ list(z = zst)}

# Parameters monitored + only essential PPC parameters ("dfrep1" and"eval")
params <- c(      # Random species-level parameters first ...
  "alpha.lpsi1", "beta.lpsi1.sf", "alpha.lphi.year", "beta.lphi.sf", 
  "alpha.lgamma.year", "beta.lgamma.sf", "alpha.lpp.year", "beta.lpp.sf",
  "beta.lphi.sfXyear", "beta.lgamma.sfXyear", "beta.lpp.sfXyear","beta.t","beta.tt",
  # ... then the hyperparameters
  "mu.alpha.lpsi1", "mean.alpha.lpsi1", "mu.beta.lpsi1.sf", "mean.beta.lpsi1.sf",
  "sd.alpha.lpsi1", "sd.beta.lpsi1.sf", "mu.alpha.lphi.year", "mean.alpha.lphi.year", 
  "mu.alpha.lgamma.year", "mean.alpha.lgamma.year", "sd.alpha.lphi.year",
  "sd.alpha.lgamma.year", "mu.beta.lphi.sf", "mean.beta.lphi.sf",
  "mu.beta.lgamma.sf", "mean.beta.lgamma.sf", "sd.beta.lphi.sf",
  "sd.beta.lgamma.sf", "mu.alpha.lpp.year", "mean.alpha.lpp.year",
  "mu.beta.lpp.sf", "mean.beta.lpp.sf", "sd.alpha.lpp.year", "sd.beta.lpp.sf",
  "mu.alpha.lphi.sfXyear", "mean.alpha.lphi.sfXyear", "mu.alpha.lgamma.sfXyear",
  "mean.alpha.lgamma.sfXyear", "sd.alpha.lphi.sfXyear", "sd.alpha.lgamma.sfXyear",
  "mu.alpha.lpp.sfXyear","mean.alpha.lpp.sfXyear","sd.alpha.lpp.sfXyear",
  "mu.beta.lt", "mean.beta.lt", "mu.beta.ltt", "mean.beta.ltt",
  # ... and finally only the essential posterior predictive check pars
  "dfrep1", "eval")

# MCMC settings
na <- 5000; ni <- 30000; nt <- 20; nb <- 10000; nc <- 3
# na <- 500; ni <- 200; nt <- 5; nb <- 50; nc <- 2


# Run model, check convergence and summarize posteriors
tic()
out16X <- jags(bdata, inits, params, "DCM16X.txt", n.adapt = na, n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)
toc()
traceplot(out16X)
which(out16X$summary[,8] > 1.1)                 # Check whether anybody > 1.1
out16X$summary[out16X$summary[,8] > 1.1,]       # Check how bad
print(out16X, 3)

# ------------------------------------------------------------------------------
## Save output

# clean output
rm(list = ls()[!ls() %in% c("out16X")])

# Save the object as a RDS file
saveRDS(object=out16X, file="Outputs/ManausMSOD_MCMCchains.rds")

# Save output as an "RData" file:
# save.image("Outputs/output_5_WriteRunJAGSModel.RData")

# Create compdata object similar to the compdata in Montpellier project
#load("DataObjects/DataArrays.RData")
sacc <- read.csv("OriginalData/SACCList2025.csv", sep = ",", h=T)
spsnames <- readRDS("Utils/spsnames.rds")

# Sort spsnames by SACC English names
spsnames <- spsnames[order(match(spsnames, sacc$Scientific.name))]

# Scientific with underscores
sps_sci_us <- gsub(" ", "_", spsnames)

# Psittacidae from SACC
psit_species <- sacc$Scientific.name[sacc$Family == "Psittacidae"]
psit_species <- gsub(" ", "_", psit_species)

# All species including psittacidae
spsnames_psit <- sps_sci_us

# Only species that are NOT psittacidae
spsnames_nonpsit <- sps_sci_us[!sps_sci_us %in% psit_species]

# Functions
space_name <- function(x) gsub("_", " ", x)
abr <- function(x) abbreviate(x, minlength = 3)

# Remove attributes (matching compdata original)
clean_vec <- function(x) { attributes(x) <- NULL; x }

# Build compdata with identical structure
compdata_bn <- list(
  timepsi = 1:5,
  timegam = 1:4,
  
  spsnames_psit   = clean_vec(spsnames_psit),
  abrnames_psit   = clean_vec(abr(space_name(spsnames_psit))),
  spacenames_psit = clean_vec(space_name(spsnames_psit)),
  
  spsnames       = clean_vec(spsnames_nonpsit),
  abrnames       = clean_vec(abr(space_name(spsnames_nonpsit))),
  spacenames     = clean_vec(space_name(spsnames_nonpsit)),
  
  SACCnames      = clean_vec(spsnames_nonpsit),
  SACCabrnames   = clean_vec(abr(space_name(spsnames_nonpsit))),
  SACCspacenames = clean_vec(space_name(spsnames_nonpsit))
)

str(compdata_bn) # seems to be working
saveRDS(object=compdata_bn, file="Utils/compdata.rds")
