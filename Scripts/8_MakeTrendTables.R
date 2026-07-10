# make matrix with coefficients of temporal variation in phi, gamma, pp and psi for 
# both OG and SF for cross-species mean (mu) and for each species:
# (tab has nsps+1 rows because the first row is for cross-sps values)

# find jags output and sps name list (for previous analysis)
#load("/Users/gferraz/OneDrive/Projects/NSpsOccDynManaus/Analysis/NSODManaus_Marc/Outputs from jags/output_NSODMao_Code16_vec.RData")
#load("/Users/gferraz/OneDrive/Projects/NSpsOccDynManaus/Analysis/NSODManaus_Marc/Outputs from jags/output_NSODMao_Code16X_array_deadSerious_trs_biggerequal85.RData")
#load("/Users/gferraz/OneDrive/Projects/NSpsOccDynManaus/Analysis/NSODManaus_Marc/Outputs from jags/output_NSODMao_Code16X_array_deadSerious_trs_biggerthan99.RData")


# For BirdNET Analysis:
compdata <- readRDS("Utils/compdata.rds")
output <- readRDS("Outputs/ManausMSOD_MCMCchains.rds")

nsps <- length(output$mean$alpha.lpsi1)              # get number of species
niter <- length(output$sims.list$mu.beta.lgamma.sf)  # get number of MCMC iters
years <- c(1:dim(output$sims.list$mean.alpha.lpp.year)[2])
nyears <- length(years)      # get number of years

below0 <- tabtrend <- tablow <- tabhigh <- matrix(NA, nsps+1, 8)
colnames(below0) <- c("phi OG", "gamma OG", "pp OG", "psi OG",
                      "phi SF", "gamma SF", "pp SF", "psi SF")
colnames(tabtrend) <- colnames(tablow) <- colnames(tabhigh) <- colnames(below0) 
sp_names <- c("All-species (intercept)", compdata$SACCspacenames)
# species names on their own column is convenient for Rmd knitting

# ------------------------------------------------------------------------------
# fill trends for cross-species phi OG
csps.lphiog <- output$sims.list$mu.alpha.lphi.year # logit space yearly mean phi
cy <- csps.lphiog # current values of y in the regression?
betas<-rep(NA,niter)  
for(i in 1:niter){
  betas[i]<-lm(cy[i,]~c(1:4))$coefficients[2]  
}
tabtrend[1,1] <- mean(betas)
quants <- quantile(betas,probs=c(0.025,0.975))
tablow[1,1] <- quants[1]
tabhigh[1,1] <- quants[2]
below0[1,1] <- sum(betas<0)/length(betas)
# ------------------------------------------------------------------------------
# cross-species gamma OG
csps.lgammaog <- output$sims.list$mu.alpha.lgamma.year
cy <- csps.lgammaog
betas<-rep(NA,niter)  
for(i in 1:niter){
  betas[i]<-lm(cy[i,]~c(1:4))$coefficients[2]  
}
tabtrend[1,2] <- mean(betas)
quants <- quantile(betas,probs=c(0.025,0.975))
tablow[1,2] <- quants[1]
tabhigh[1,2] <- quants[2]
below0[1,2] <- sum(betas<0)/length(betas)
# ------------------------------------------------------------------------------
# cross-species pp OG
csps.lppog <- output$sims.list$mu.alpha.lpp.year
cy <- csps.lppog
betas<-rep(NA,niter)  
for(i in 1:niter){
  betas[i]<-lm(cy[i,]~c(1:5))$coefficients[2]  
}
tabtrend[1,3] <- mean(betas)
quants <- quantile(betas,probs=c(0.025,0.975))
tablow[1,3] <- quants[1]
tabhigh[1,3] <- quants[2]
below0[1,3] <- sum(betas<0)/length(betas)
# ------------------------------------------------------------------------------
# cross-species psi OG
csps.psiog <- matrix(NA,niter,nyears)
csps.psiog[,1] <- plogis(output$sims.list$mu.alpha.lpsi1)
for (t in 2:nyears) {
  csps.psiog[,t] <- csps.psiog[,t-1] * plogis(csps.lphiog[,t-1]) +
                    (1 - csps.psiog[,t-1]) * plogis(csps.lgammaog[,t-1])
}
cy <- log( csps.psiog / (1-csps.psiog) )
betas<-rep(NA,niter)  
for(i in 1:niter){
  betas[i]<-lm(cy[i,]~c(1:5))$coefficients[2]  
}
tabtrend[1,4] <- mean(betas)
quants <- quantile(betas,probs=c(0.025,0.975))
tablow[1,4] <- quants[1]
tabhigh[1,4] <- quants[2]
below0[1,4] <- sum(betas<0)/length(betas)
# ------------------------------------------------------------------------------
# cross-species phi SF
csps.lphisf <- output$sims.list$mu.alpha.lphi.year + 
               matrix(rep(output$sims.list$mu.beta.lphi.sf,4),ncol=4) +
               cbind(rep(0,niter),output$sims.list$mu.alpha.lphi.sfXyear)
cy <- csps.lphisf
betas<-rep(NA,niter)  
for(i in 1:niter){
  betas[i]<-lm(cy[i,]~c(1:4))$coefficients[2]  
}
tabtrend[1,5] <- mean(betas)
quants <- quantile(betas,probs=c(0.025,0.975))
tablow[1,5] <- quants[1]
tabhigh[1,5] <- quants[2]
below0[1,5] <- sum(betas<0)/length(betas)
# ------------------------------------------------------------------------------
# cross-species gamma SF
csps.lgammasf <- output$sims.list$mu.alpha.lgamma.year + 
                 matrix(rep(output$sims.list$mu.beta.lgamma.sf,4),ncol=4) +
                 cbind(rep(0,niter),output$sims.list$mu.alpha.lgamma.sfXyear)
cy<- csps.lgammasf
betas<-rep(NA,niter)  
for(i in 1:niter){
  betas[i]<-lm(cy[i,]~c(1:4))$coefficients[2]  
}
tabtrend[1,6] <- mean(betas)
quants <- quantile(betas,probs=c(0.025,0.975))
tablow[1,6] <- quants[1]
tabhigh[1,6] <- quants[2]
below0[1,6] <- sum(betas<0)/length(betas)
# ------------------------------------------------------------------------------
# cross-species pp SF
csps.lppsf <- output$sims.list$mu.alpha.lpp.year + 
              matrix(rep(output$sims.list$mu.beta.lpp.sf,5),ncol=5) +
              cbind(rep(0,niter),output$sims.list$mu.alpha.lpp.sfXyear)
cy <- csps.lppsf
betas<-rep(NA,niter)  
for(i in 1:niter){
  betas[i]<-lm(cy[i,]~c(1:5))$coefficients[2]  
}
tabtrend[1,7] <- mean(betas)
quants <- quantile(betas,probs=c(0.025,0.975))
tablow[1,7] <- quants[1]
tabhigh[1,7] <- quants[2]
below0[1,7] <- sum(betas<0)/length(betas)
# ------------------------------------------------------------------------------
# cross-species psi SF
csps.psisf <- matrix(NA,niter,nyears)
csps.psisf[,1] <- plogis(output$sims.list$mu.alpha.lpsi1 + 
                         output$sims.list$mu.beta.lpsi1.sf)
for (t in 2:nyears) { # using values in probability space
  csps.psisf[,t] <- csps.psisf[,t-1] * plogis(csps.lphisf[,t-1]) +
                    (1 - csps.psisf[,t-1]) * plogis(csps.lgammasf[,t-1]) 
}
cy <- log(csps.psiog / (1-csps.psiog))
betas<-rep(NA,niter)  
for(i in 1:niter){
  betas[i]<-lm(cy[i,]~c(1:5))$coefficients[2]  
}
tabtrend[1,8] <- mean(betas)
quants <- quantile(betas,probs=c(0.025,0.975))
tablow[1,8] <- quants[1]
tabhigh[1,8] <- quants[2]
below0[1,8] <- sum(betas<0)/length(betas)

# ==============================================================================
## Species-specific phi, gamma, pp and psi OG and SF

for (k in 1:nsps) {
  print(k)
  # here we have added the "csn" index to generate a table in SACC species order
  csn <- which(compdata$spacenames == compdata$SACCspacenames[k]) # current species number
  # (adjust output sps order to SACC sps order)
  
  # notice that we're using the "k+1" row because the first row is being used for cross-species values
  
  # ----------------------------------------------------------------------------
  # species-specific phi OG
  lphiog <- output$sims.list$alpha.lphi.year[,,csn] # logit space yearly mean phi
  cy <- lphiog # current values of y in the regression?
  betas<-rep(NA,niter)  
  for(i in 1:niter){
    betas[i]<-lm(cy[i,]~c(1:4))$coefficients[2]  
  }
  tabtrend[k+1,1] <- mean(betas)
  quants <- quantile(betas,probs=c(0.025,0.975))
  tablow[k+1,1] <- quants[1]
  tabhigh[k+1,1] <- quants[2]
  below0[k+1,1] <- sum(betas<0)/length(betas)
  
  # ----------------------------------------------------------------------------
  # species-specific gamma OG
  lgammaog <- output$sims.list$alpha.lgamma.year[,,csn] # logit space yearly mean gamma
  cy <- lgammaog
  betas<-rep(NA,niter)  
  for(i in 1:niter){
    betas[i]<-lm(cy[i,]~c(1:4))$coefficients[2]  
  }
  tabtrend[k+1,2] <- mean(betas)
  quants <- quantile(betas,probs=c(0.025,0.975))
  tablow[k+1,2] <- quants[1]
  tabhigh[k+1,2] <- quants[2]
  below0[k+1,2] <- sum(betas<0)/length(betas)
  
  # ----------------------------------------------------------------------------
  # species-specific pp OG
  lppog <- output$sims.list$alpha.lpp.year[,,csn] # logit space yearly mean pp
  cy <- lppog
  betas<-rep(NA,niter)  
  for(i in 1:niter){
    betas[i]<-lm(cy[i,]~c(1:5))$coefficients[2]  
  }
  tabtrend[k+1,3] <- mean(betas)
  quants <- quantile(betas,probs=c(0.025,0.975))
  tablow[k+1,3] <- quants[1]
  tabhigh[k+1,3] <- quants[2]
  below0[k+1,3] <- sum(betas<0)/length(betas)
  
  # ----------------------------------------------------------------------------
  # species-specific psi OG
  psiog <- matrix(NA,niter,nyears)
  psiog[,1] <- plogis(output$sims.list$alpha.lpsi1[,csn])
  for (t in 2:nyears) {
    psiog[,t] <- psiog[,t-1] * plogis(lphiog[,t-1]) +
      (1 - psiog[,t-1]) * plogis(lgammaog[,t-1])
  }
  cy <- log( psiog / (1-psiog) )
  betas<-rep(NA,niter)  
  for(i in 1:niter){
    betas[i]<-lm(cy[i,]~c(1:5))$coefficients[2]  
  }
  tabtrend[k+1,4] <- mean(betas)
  quants <- quantile(betas,probs=c(0.025,0.975))
  tablow[k+1,4] <- quants[1]
  tabhigh[k+1,4] <- quants[2]
  below0[k+1,4] <- sum(betas<0)/length(betas)
  
  # ----------------------------------------------------------------------------
  # species-specific phi SF
  lphisf <- output$sims.list$alpha.lphi.year[,,csn] + 
    matrix(rep(output$sims.list$beta.lphi.sf[,csn],4),ncol=4) +
    output$sims.list$beta.lphi.sfXyear[,,csn]
  cy <- lphisf
  betas<-rep(NA,niter)  
  for(i in 1:niter){
    betas[i]<-lm(cy[i,]~c(1:4))$coefficients[2]  
  }
  tabtrend[k+1,5] <- mean(betas)
  quants <- quantile(betas,probs=c(0.025,0.975))
  tablow[k+1,5] <- quants[1]
  tabhigh[k+1,5] <- quants[2]
  below0[k+1,5] <- sum(betas<0)/length(betas)
  
  # ----------------------------------------------------------------------------
  # species-specific gamma SF
  lgammasf <- output$sims.list$alpha.lgamma.year[,,csn] + 
              matrix(rep(output$sims.list$beta.lgamma.sf[,csn],4),ncol=4) +
              output$sims.list$beta.lgamma.sfXyear[,,csn]
  cy <- lgammasf
  betas<-rep(NA,niter)  
  for(i in 1:niter){
    betas[i]<-lm(cy[i,]~c(1:4))$coefficients[2]  
  }
  tabtrend[k+1,6] <- mean(betas)
  quants <- quantile(betas,probs=c(0.025,0.975))
  tablow[k+1,6] <- quants[1]
  tabhigh[k+1,6] <- quants[2]
  below0[k+1,6] <- sum(betas<0)/length(betas)
  
  # ----------------------------------------------------------------------------
  # species-specific pp SF
  lppsf <- output$sims.list$alpha.lpp.year[,,csn] + 
    matrix(rep(output$sims.list$beta.lpp.sf[,csn],5),ncol=5) +
    output$sims.list$beta.lpp.sfXyear[,,csn]
  cy <- lppsf
  betas<-rep(NA,niter)  
  for(i in 1:niter){
    betas[i]<-lm(cy[i,]~c(1:5))$coefficients[2]  
  }
  tabtrend[k+1,7] <- mean(betas)
  quants <- quantile(betas,probs=c(0.025,0.975))
  tablow[k+1,7] <- quants[1]
  tabhigh[k+1,7] <- quants[2]
  below0[k+1,7] <- sum(betas<0)/length(betas)
  
  # ------------------------------------------------------------------------------
  # species-specific psi SF
  psisf <- matrix(NA,niter,nyears)
  psisf[,1] <- plogis(output$sims.list$alpha.lpsi1[,csn] + 
                        output$sims.list$beta.lpsi1.sf[,csn])
  for (t in 2:nyears) {
    psisf[,t] <- psisf[,t-1] * plogis(lphisf[,t-1]) +
      (1 - psisf[,t-1]) * plogis(lgammasf[,t-1])
  }
  cy <- log( psisf / (1-psisf) )
  betas<-rep(NA,niter)  
  for(i in 1:niter){
    betas[i]<-lm(cy[i,]~c(1:5))$coefficients[2]  
  }
  tabtrend[k+1,8] <- mean(betas)
  quants <- quantile(betas,probs=c(0.025,0.975))
  tablow[k+1,8] <- quants[1]
  tabhigh[k+1,8] <- quants[2]
  below0[k+1,8] <- sum(betas<0)/length(betas)
  
} # k

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# create list with all the tables and save in RData file
# make sure to get the right number in the names
# '09' means the criteria for detection was a probability of vocalization presence
# greater than 0.9 in Ulisses' Y matrix
trends09 <- list(mean09=tabtrend,low09=tablow,high09=tabhigh,under09=below0,sps=sp_names)
#save(trends09,file="TrendTables09.RData")
saveRDS(trends09, file="DataObjects/TrendTables09.rds")

trends085 <- list(mean085=tabtrend,low085=tablow,high085=tabhigh,under085=below0,sps=sp_names)
#save(trends085,file="TrendTables085.RData")
saveRDS(trends085, file="DataObjects/TrendTables085.rds")

#trends099 <- list(mean099=tabtrend,low099=tablow,high099=tabhigh,under099=below0,sps=sp_names)
#save(trends099,file="TrendTables099.RData")
