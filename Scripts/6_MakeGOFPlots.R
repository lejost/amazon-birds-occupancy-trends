## Code for extracting goodness of fits stats from minimum variables stored
## in jags output from model 16
## needs only "y", "dfrep1" and "eval"

################################################################################
## Read
y <- readRDS("DataObjects/y.rds")
output <- readRDS("Outputs/ManausMSOD_MCMCchains.rds")

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
FTobs <- (sweep(-sqrt(eval), 2:4, -sqrt(df1)))^2
FTrep <- (sqrt(dfrep1) - sqrt(eval))^2
FTratio <- FTobs / (FTrep + 0.0001) 

tmp <- ((df1-1)>=0)*1
tmprep <- ((dfrep1-1)>=0)*1

df2 <- apply(tmp,c(2,3),sum)
dfrep2 <- apply(tmprep,c(1,3,4),sum)
diff2 <- sweep(-dfrep2,2:3,-df2)

FTobs.tot <- apply(FTobs,1,sum)
FTrep.tot <- apply(FTrep,1,sum)

mean.dfrep1 <- output$mean$dfrep1
mean.dfrep2 <- apply(dfrep2,c(2,3),mean)
mean.diff2 <- apply(diff2,c(2,3),mean)

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
plot(df1, mean.dfrep1, xlab = 'Actual data set', pch = 16, col = rgb(0,0,0,0.3),
     ylab = "'Ideal' data sets", frame = FALSE, 
     main = 'Number of occasions detected (out of 13)\nper site, year and species')
abline(0, 1)
# That looks pretty good

# Same for df2
plot(df2, mean.dfrep2, xlab = 'Actual data set', pch = 16, col = rgb(0,0,0,0.3),
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
