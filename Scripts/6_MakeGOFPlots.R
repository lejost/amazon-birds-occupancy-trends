## Code for extracting goodness of fits stats from minimum variables stored
## in jags output from model 16
## needs only "y", "dfrep1" and "eval"

################################################################################
## First get y array using code from the model script
Exc <- 1 # 1 for excluding and 0 for not excluding sites 
         # that were only sampled in year 2

## Read
load("/Users/gferraz/OneDrive/Projects/NSpsOccDynManaus/Analysis/NSODManaus_Marc/Data/NSODMaoUlisses.RData")
# load("C:/Users/fabio/OneDrive/NSpsOccDynManaus/Analysis/NSODManaus_Marc/Data/NSODMaoUlisses.RData")
## Create arguments and optionally clean sites that were only sampled in year 2:
maxdy <- c(6,8,13,6,5)
# get names of the first site of each group
gnames<-rep("x",length(unique(groups)))
for(i in 1:length(unique(groups))) {
  gnames[i]<-shortSITES$Spot[min(which(unique(groups)[i]==groups))]
}
## Exclude group of unidentified species from family Psittacidae
NDETSG <- NDETSG[,,-which(dimnames(NDETSG)[[3]]=="Fam_Psittacidae")]
## If Exc == 1, exclude sites that were only sampled during the second year
if(Exc==1) {
  # get indices of sites (grouped sites) that were only sampled on year 2
  check <- function(x) all(is.na(x))
  only2 <- which(apply(DATSG[c(1:6,15:38),],MARGIN = 2,FUN = check))
  keep <- c(1:214)[-only2]
  gnames <- gnames[keep]
  NDETSG <- NDETSG[keep,-c(13,14),]
  EFFGsits <- EFFGsits[keep,-c(13,14)]
  maxdy<-c(6,6,13,6,5)
}
nsites <- nrow(NDETSG)
nsurveys <- max(maxdy)
nyears <- length(maxdy)
nspecies <- dim(NDETSG)[[3]]
# preencher vetor secfor que informa quais grupos sao em
# floresta secundaria (0 = floresta primaria, 1 = floresta secundaria)
secfor <- rep(NA, nsites)
for(i in 1:length(secfor)) {
  secfor[i] <- shortSITES[which(shortSITES$Spot==gnames[i]),2]
}
## Create and fill "sampsize" (sample size):
sampsize <- array(NA, dim = c(nsites, nsurveys, nyears))
for(i in 1:nsites) {  # i
  for (t in 1:nyears) { # t
    if (t == 1) {
      sampsize[i,1:maxdy[t],t] <- EFFGsits[i, 1:maxdy[t]]
    } else {
      sampsize[i,1:maxdy[t],t] <- EFFGsits[i, (sum(maxdy[1:(t-1)])+1):(sum(maxdy[1:(t-1)])+maxdy[t])] }
  } # t
} # i
sampsize[which(is.na(sampsize))] <- 0
## Turn "NDETSG" into "NDETBIN" (binary NDETSG)
NDETBIN <- NDETSG
NDETBIN[which(NDETBIN>0)] <- 1
## Create and fill four-dimensional "NEWDATAY", which replaces "dat$y"
## and includes one dimension for species
NEWDATAY <- array(data = NA, dim = c(nsites, nsurveys, nyears, nspecies))
for (s in 1:nspecies) { # s
  for (t in 1:nyears) { # t
    if (t == 1) {
      NEWDATAY[, 1:maxdy[t], t, s] <- NDETBIN[, 1:maxdy[t], s]
    } else {
      NEWDATAY[, 1:maxdy[t], t, s] <- NDETBIN[,(sum(maxdy[1:(t-1)])+1):(sum(maxdy[1:(t-1)])+maxdy[t]), s] }
  } # t
} # s
y <- NEWDATAY

################################################################################
# get basic numbers from y
nsites <- dim(y)[1]
nsurveys <- dim(y)[2]
nyears <- dim(y)[3]
nsps <- dim(y)[4]
niter <- dim(output$sims.list$dfrep1)[1]
################################################################################
# get gof variables from MCMC output
eval <- output$sims.list$eval
dfrep1 <- output$sims.list$dfrep1
mean.dfrep1 <- output$mean$dfrep1
################################################################################
# derive  necessary gof objects
df1 <- apply(y,c(1,3,4),function(x)sum(x,na.rm=TRUE)) # missing the niter dim
diff1 <- sweep(-dfrep1,2:4,-df1) # sweep subtracts, so we need the minus signs
FTobs <- (sweep(-sqrt(eval),2:4,-df1))^2 # ERRADO! Falta a sqrt(df1)
FTrep <- (sqrt(dfrep1) - sqrt(eval))^2
FTratio <- FTobs / (FTrep + 0.0001) 

tmp <- ((df1-1)>=0)*1
tmprep <- ((dfrep1-1)>=0)*1

df2 <- apply(tmp,c(2,3),sum)
dfrep2 <- apply(tmprep,c(1,3,4),sum)
diff2 <- sweep(-dfrep2,2:3,-df2)

FTobs.tot <- apply(FTobs,1,sum)
FTrep.tot <- apply(FTrep,1,sum)




# test with objects from obtained directly from jags
load("/Users/gferraz/OneDrive/Projects/NSpsOccDynManaus/Analysis/NSODManaus_Marc/Outputs from jags/output_NSODMao_Code16X_array_deadSerious.RData")
load("C:/Users/fabio/OneDrive/NSpsOccDynManaus/Analysis/NSODManaus_Marc/Outputs from jags/output_NSODMao_Code16X_array_deadSerious.RData")

output<-out16X
FTobs.tot <- output$sims.list$FTobs.tot
FTrep.tot <- output$sims.list$FTrep.tot
mean.df1 <- output$mean$df1
mean.dfrep1 <- output$mean$dfrep1
mean.df2 <- output$mean$df2
mean.dfrep2 <- output$mean$dfrep2
mean.diff2 <- output$mean$diff2

# ------------------------------------------------------------------------------
# Make a PPC plot rightaway and compute a bayesian p-value
mrg<-0.1*(range(FTobs.tot,FTrep.tot)[2]-range(FTobs.tot,FTrep.tot)[1] )
minv <- range(FTobs.tot,FTrep.tot)[1] - mrg
maxv <- range(FTobs.tot,FTrep.tot)[2] + mrg
plot(FTobs.tot, FTrep.tot,
     xlab = 'FT statistic observed data', ylab = "FT statistic 'ideal' data",
     pch =16, col = rgb(0,0,0,0.3), xlim = c(minv,maxv), ylim = c(minv,maxv))
abline(0, 1)
mean(FTrep.tot > FTobs.tot)   # Bayesian p-value

# Same plot and color red those MCMC samples where rep > obs
TTT <- FTrep.tot > FTobs.tot
plot(FTobs.tot, FTrep.tot, 
     xlab = 'FT statistic observed data', ylab = "FT statistic 'ideal' data",
     pch = 16, col = rgb(0,0,0,0.3), xlim = c(6000, 7200), ylim = c(6000, 7200))
abline(0, 1)
points(FTobs.tot[TTT], FTrep.tot[TTT], 
       pch = 16, col = rgb(1,0,0,0.3))
# The Bayesian p-value is the proportion of red circles

# "Lack of fit ratio"
hist(FTobs.tot/FTrep.tot, breaks = 100, col = 'grey',
     main = "Ratio of FT statistic in actual and 'ideal' data ", xlim = c(0.95, 1.25))
abline(v = mean(FTobs.tot/FTrep.tot), lwd = 3, col = 'blue')
mean(FTobs.tot/FTrep.tot)
# [1] 1.068086, instead of previous "[1] 1.13673" using "out16X"


# Now we can plot the expected (under the 'ideal data') frequencies 
# against the observed ones: for df1
plot(mean.df1, mean.dfrep1, xlab = 'Actual data set', pch = 16, col = rgb(0,0,0,0.3),
     ylab = "'Ideal' data sets", frame = FALSE, 
     main = 'Number of occasions detected (out of 13)\nper site, year and species')
abline(0, 1)
# That looks pretty good

# Same for df2
plot(mean.df2, mean.dfrep2, xlab = 'Actual data set', pch = 16, col = rgb(0,0,0,0.3),
     ylab = "'Ideal' data sets", frame = FALSE, 
     main = 'Number of sites at which detected (out of 158)\nper species and year')
abline(0, 1)
# That looks almost too nice. Not sure it has any ability to assess goodness of fit ....

# Can summarize the differences by species and year to see when things don't fit so well
mean.diff2    # Diff actual minus expected under perfect data
# Overall tends to be negative, thus in our data we have a little fewer
# sites at which a species are actually observed than what we would expect
# under the model

# Can see which sites are different
# For this we summarize diff1 by site now
str(output$sims.list)
tmp <- apply(output$sims.list$diff1, 2:4, mean)
tmp <- apply(tmp, 1:2, sum)         # This is sum of mean difference of detection frequencies
# per site and year
tmp                 # Here can see which sites and years are 'striking' in the observed data
# compared to the kind of data one would expect under the model
