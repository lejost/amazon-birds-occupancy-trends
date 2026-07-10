
## Code 4, for processing binary data

## Receives the output from code '3_DataPreparation.R' and turns the 3D array 
## "NDETSG" (which contains the number of detections by site, survey and 
## species) into a binary detection/non-detection version called "y", which will
## be input to both the occupancy model and the GOF plots. Also creates 
## "secfor", a vector which informs which groups are from secondary forest, and 
## "sampsize", a 3D array with the sample sizes
## -----------------------------------------------------------------------------

## Clean work space:
rm(list=ls())

## Load data input
load("DataObjects/DataArrays.RData")
sacc <- read.csv("OriginalData/SACCList2025.csv", sep = ",", h=T)

# 1 for excluding and 0 for not excluding sites that were only sampled in year 2
Exc <- 1

# Create arguments and optionally clean sites that were only sampled in year 2:
maxdy <- c(6,8,9,6,5)

# get names of the first site of each group
gnames <- names(groups)

# Exclude data on the family Psittacidae
NDETSG <- NDETSG[,,-c(87:101)]

# If Exc == 1, exclude sites that were only sampled during the second year
if(Exc==1) {
  # get indices of sites (grouped sites) that were only sampled on year 2
  check <- function(x) all(is.na(x))
  only2 <- which(apply(DATSG[c(1:6,15:34),],MARGIN = 2,FUN = check))
  keep <- c(1:209)[-only2]
  gnames <- gnames[keep]
  NDETSG <- NDETSG[keep,-c(13,14),] # after removing sites from only yr 2, cols 13 and 14 are empty
  EFFGsits <- EFFGsits[keep,-c(13,14)]
  EFFGmins <- EFFGmins[keep,-c(13,14)]
  DATSG <- DATSG[-c(13,14),keep]
  rownames(DATSG)<-seq(1,dim(DATSG)[1])
  maxdy<-c(6,6,9,6,5)
}

## Define array dimensions and species names
nsites <- nrow(NDETSG)
nsurveys <- max(maxdy)
nyears <- length(maxdy)
nspecies <- dim(NDETSG)[[3]]
spsnames <- dimnames(NDETSG)[[3]]
spsnames <- unname(setNames(sacc$Scientific.name, sacc$English.name)[spsnames])
dimnames(NDETSG)[[3]] <- spsnames
years <- c("2010","2011","2012","2013","2014")

# Fill 'secfor' vector which informs which groups are from secondary forest
# (0 = old growth, 1 = secondary forest)
secfor <- rep(NA, nsites)
for(i in 1:length(secfor)) {
  secfor[i] <- shortSITES[which(shortSITES$Spot==gnames[i]),2]
}

# Create a 3D array with BirdNET detection data:

## Create a 3D array ('det_array') with BirdNET detection data
# Dimensions:
#   [1] rows = years + total
#   [2] cols = metrics
#   [3] layers = species
rows <- c(years, "TOTAL")
cols <- c("n_det_OG", "n_det_SF", "n_sites_OG", "n_sites_SF")

det_array <- array(data=NA, dim = c(6, 4, nspecies),
                   dimnames = list(rows, cols, spsnames))

og <- which(secfor==0)
sf <- which(secfor==1)

# Fill det_array
# det_array dimensions: [year, metric, species]
for (i in 1:nspecies){ # loop through species
  cs <- NDETSG[ , , i] # cs = current species
  
  for (t in 1:nyears) {
    #print(t)
    fd <- sum(maxdy[1:t])-maxdy[t]+1    # first day for year t
    ld <- sum(maxdy[1:t])               # last day for year t
    cd <- cs[ , fd:ld ]                 # detections for year t
    
    det_array[t, 1, i] <- sum(cd[og, ], na.rm = TRUE)  # n_det_OG
    det_array[t, 2, i] <- sum(cd[sf, ], na.rm = TRUE)  # n_det_SF
    
    det_array[t, 3, i] <- sum(rowSums(cd[og,],na.rm=TRUE)>0) # n_sites_OG
    det_array[t, 4, i] <- sum(rowSums(cd[sf,],na.rm=TRUE)>0) # n_sites_SF
    
  }
  
  det_array[6, 1, i] <- sum(cs[og, ], na.rm = TRUE)  # n_det_OG
  det_array[6, 2, i] <- sum(cs[sf, ], na.rm = TRUE)  # n_det_SF
  
  det_array[6, 3, i] <- sum(rowSums(cs[og,],na.rm=TRUE)>0) # n_sites_OG
  det_array[6, 4, i] <- sum(rowSums(cs[sf,],na.rm=TRUE)>0) # n_sites_SF
  
}


## Select species based on data summarized in det_array
## 1. Identify species with no detections (just in case!)
nodets <- spsnames[which(colSums(det_array[6,1:4,])==0)]
## 2. Identify species that appear in only one type of forest 
##    and only one year
only1tpy <- spsnames[which(colSums(colSums(det_array[1:5,3:4,]>0))==1)]
## Supplemental stuff maybe not to use
## 3. Identify species that appear in only one site and year
#only1sit <- spsnames[which(colSums(det_array[6,3:4,])==1)]
## 4. Identify species that appear in only one type of forest
#only1typ <- spsnames[which(colSums(det_array[6,3:4,]==0)==1)]
#only1typ <- only1typ[which(!only1typ%in%only1sit)]

# Create a set with the union of only1tpy and nodets species names
set_union <- c(only1tpy, nodets)

# Remove excluded species
NDETSG <- NDETSG[,,-c(which(dimnames(NDETSG)[[3]] %in% set_union))]

# Update the vectors containing information related to the species
nspecies <- dim(NDETSG)[[3]]
spsnames <- dimnames(NDETSG)[[3]]

# Turn "NDETSG" into "NDETBIN" (binary NDETSG)
NDETBIN <- NDETSG
NDETBIN[which(NDETBIN>0)] <- 1

# Create and fill four-dimensional "NEWDATAY"
NEWDATAY <- array(data = NA, dim = c(nsites, nsurveys, nyears, nspecies))
# NEWDATAY: dim: 157 sites x 9 max survey days (max) x 5 years x 183 species

for (s in 1:nspecies) { # s species
  for (t in 1:nyears) { # t years
    
    if (t == 1) {
      NEWDATAY[, 1:maxdy[t], t, s] <- NDETBIN[, 1:maxdy[t], s]
    } else {
      NEWDATAY[, 1:maxdy[t], t, s] <- NDETBIN[,(sum(maxdy[1:(t-1)])+1):(sum(maxdy[1:(t-1)])+maxdy[t]), s] }
    
  } # t
} # s

y <- NEWDATAY # copy 4D array


# 3 - Construir matriz 2D com dia juliano mediano da amostragem de cada sitio
# (grupo) e ano.

# matrix site x year
stjul <- matrix(NA, nsites, nyears)

for (i in 1:nsites) {
  print(i)
  for (t in 1:nyears) {
    print(t)
    originday <- as.Date(paste(years[t],"-01-01", sep = ""))
    fd <- sum(maxdy[1:t])-maxdy[t]+1 # first day
    ld <- sum(maxdy[1:t]) # last day
    if(any(!is.na(DATSG[fd:ld,i]))) {
      firstdate <- DATSG[fd:ld,i][min(which(!is.na(DATSG[fd:ld,i])))]
      firstdate <- julian(firstdate, origin = originday)
      lastdate <- DATSG[fd:ld,i][max(which(!is.na(DATSG[fd:ld,i])))]
      lastdate <- julian(lastdate, origin = originday)
      stjul[i,t] <- median(seq(firstdate,lastdate))
    }
  }
}
# estandardizar os dias medianos em "stjul" pela subtracao da media e divisao
# pelo desvio padrao dos dias julianos de todos os anos (sem compartimentalizar
# a estandardizacao por ano, para nao perder a representacao de quais anos tem
# a amostragem mais cedo ou mais tarde).

# obs.: A "media" na frase anterior eh o dia juliano medio, ou seja, media de
# todos os dias julianos medianos ja obtidos em "stjul". O mesmo para o desvio.

# PROBLEMA: A media e o desvio padrao dos dias julianos devem ser calculados
# apenas com base naqueles grupos de sitios que entraram na analise -> na.rm = T

julmedia <- mean(stjul, na.rm = TRUE)
julsd <- sd(stjul, na.rm = TRUE)

stjul <- stjul - julmedia
stjul <- stjul/julsd

# ------------------------------------------------------------------------------
# Create and fill "sampsize" (sample size):
sampsize <- array(NA, dim = c(nsites, nsurveys, nyears))
# dim sampsize = 157 sites x 9 surveys x 5 years

for(i in 1:nsites) {  # i
  for (t in 1:nyears) { # t
    if (t == 1) {
      sampsize[i,1:maxdy[t],t] <- EFFGsits[i, 1:maxdy[t]]
    } else {
      sampsize[i,1:maxdy[t],t] <- EFFGsits[i, (sum(maxdy[1:(t-1)])+1):(sum(maxdy[1:(t-1)])+maxdy[t])] }
  } # t
} # i

sampsize[which(is.na(sampsize))] <- 0

# ------------------------------------------------------------------------------

## Save output

# clean output
rm(list = ls()[!ls() %in% c("y","secfor","sampsize","stjul")])

# Save output as an "RData" file:
save.image("DataObjects/BinaryDataArray.RData")
saveRDS()

# Save sps names object
saveRDS(object = spsnames, file = "Utils/spsnames.rds")
