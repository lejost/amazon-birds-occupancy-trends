## Processing 1:
#
# Script to organize BirdNET detection data so that the sites where BirdNET
# detected at least one species are the same sites as those in 'SITES',
# which was used as input in the previous occupancy dynamics model.
# After matching the sites, this script creates a vector that groups the sites
# by clustering them based on the maximum distance between them.
#-------------------------------------------------------------------------------
rm(list=ls())

# 1. Load Data and Functions
library(dplyr)
# library(readxl) for 2.1 Section
# library(tuneR) # for 'GetDuration' function

# Load data:
# "BNout" is an RDS file with the following dimensions:
# - Rows: bird species detections by BirdNet in May-July 2025
# - Columns: information related to the detection (file name, start time
#           of the audio segment, end time of the audio segment, species'
#           English name, and confidence)
BNout <- readRDS("CleanedData/BNetProcOutput.rds")

# "SITES" is a CSV file with the following dimensions:
# - Rows: sites
# - Columns: ID data of each site (name, forest type (primary = 0,
#           secondary = 1), geographic coordinates (X and Y))
SITES <- read.csv("OriginalData/locations_diurnal.txt", head = TRUE, sep = "\t")

# Load necessary functions
source("Functions/CleanSiteName.R")
# source("Functions/GetDuration.R")

# 2. Prepare and Clean Data
# 2.1 List Valid Audio Files ----
# No need to run this section if ValidAudioFiles.rds is already in the
# 'CleanedData" folder, in the intended form.
# That is, all audio files that are in HD "Thing 1" and are not corrupted
# List all audio files from the HD and prepare their base names
#hd <- "E:/GravOriginais" # get path to HD
#arq_hd_files_orig <- list.files(path = hd, #gets file names with path
#                                    pattern = "\\.(wav)$",
#                                    recursive = TRUE,
#                                    full.names = FALSE)

# Remove folder with corrupted files
#folder_torm <- "2013/corrompido/"
#folder_tokeep <- !grepl(folder_torm, arq_hd_files_orig, fixed = TRUE)
# 'TRUE' means the file path does NOT contain the 'folder_torm' string.
# Filter 'arq_hd_files_orig' to exclude files from the 'folder_torm'.
#arq_hd_files <- arq_hd_files_orig[folder_tokeep]
#audio_name <- basename(arq_hd_files)
#audio_name <- sub("\\.wav$", "", audio_name)

# Get corrupted files
# The spreadsheet "CheckMissing.xlsx" lists all audio files that are in 
# HD Thing 1 but are not mentioned in the BirdNet output. Some are corrupted
# and were identified with a 1, others are just without vocalizations identified
# with 0.
#corrupted_files <- read_excel("OriginalData/CheckMissing.xlsx") %>%
#  filter(corrupted == 1) %>%
#  pull(filename)
# valid files removes corrupted files from the list of all files in HD
#valid_files <- setdiff(audio_name, corrupted_files)

# Get full paths and durations for HD audio files
#audio_paths <- file.path(hd, arq_hd_files)
#durations <- sapply(audio_paths, GetDuration)

# Format duration to HH:MM:SS
#formatted_durations <- format(as.POSIXct(durations, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S")

#ValidAudioFiles <- data.frame(
#  name = audio_name,
#  duration = unname(formatted_durations)
#)

#ValidAudioFiles <- ValidAudioFiles[ValidAudioFiles$name %in% valid_files, ]
#saveRDS(ValidAudioFiles, file="CleanedData/ValidAudioFiles.rds")
# ----
# 2.2 Make corrections to the BirdNET spot names to ensure they match the names 
# used in SITES
# Load the 'ValidAudioFiles.rds' file which was created on section 2.1 of this script.
# It only needs to be created once.
ValidAudioFiles <- readRDS("CleanedData/ValidAudioFiles.rds")

# Extract only site and year information from detections made by BirdNET
BNout <- BNout[which(BNout$FileName %in% ValidAudioFiles$name == TRUE),]
BNoutSpot <- substr(BNout$FileName, 1, 8)

# Adjust site names and years where the character count does not
# match the expected pattern
nrowstofix13 <- BNout[which(substr(BNout$FileName, 1, 8) == "CLC00013"), 1]
nrowstofix14 <- BNout[which(substr(BNout$FileName, 1, 8) == "CLC00014"), 1]
BNoutSpot[which(substr(BNout$FileName, 1, 8) == "CLC00013")] <- substr(nrowstofix13, 1, 9)
BNoutSpot[which(substr(BNout$FileName, 1, 8) == "CLC00014")] <- substr(nrowstofix14, 1, 9)
rm(nrowstofix13, nrowstofix14)

# Match the site names in SITES$Spot and BNoutSpot
# Create a temporary 'Frankenstein' version of BNoutSpot with capitalization
# adjusted to match how it appears in SITES
FBSites <- sapply(BNoutSpot, CleanSiteName)

# Fix character count issues in the BNoutSpot vector
corrections <- data.frame(
  # Column with the original names
  original = c(
    # Fix character count issues
    "CLc00013x", "CLc00014x",
    # Correct site names from 2010 that were changed in 2011
    "QUm0H15x", "QUm0Q21x", "QUmCC19x", "QUmGG23x", "QUmKK13x", 
    "QUmKK17x", "QUmLH04x", "QUm0024x",
    # Fix typos in site names from Km41
    "QMu0008x", "QMu0010x",
    # Replace site names with synonyms
    "QUm0H06x", "QUm0M12x", "QUm0Q16x", "QUm0U16x", "QUm0X24x", 
    "QUmAB16x", "QUmAB20x", "QUmCC08x", "QUmCC24x", "QUmGG18x"
  ),
  # Column with the correct names
  corrected = c(
    # Fix character count issues
    "CLc0013x", "CLc0014x",
    # Correct site names from 2010 that were changed in 2011
    "QUm0004x", "QUm0006x", "QUm0012x", "QUm0015x", "QUm0016x",
    "QUm0017x", "QUm0018x", "QUmOO24x",
    # Fix typos in site names from Km41
    "QUm0008x", "QUm0010x",
    # Replace site names with synonyms
    "QUm0003x", "QUm0005x", "QUm0020x", "QUm0007x", "QUm0008x",
    "QUm0009x", "QUm0010x", "QUm0011x", "QUm0013x", "QUm0014x"
  ),
  stringsAsFactors = FALSE # Important to avoid names getting into factores
)

for (i in 1:nrow(corrections)) {
  wname <- corrections$original[i]
  rname <- corrections$corrected[i]
  FBSites[FBSites == wname] <- rname
}

# Convert FBSites to a factor and assign it to the final BNoutSpot object
BNoutSpot <- as.factor(FBSites)
names(BNoutSpot) <- seq(1, length(BNoutSpot))
# SECTION 2 in 'Exploratory_Analysis_CheckSitesDates.R' verifies whether all 
# FBSites are in SITES. Be sure to check it before cleaning the workspace.
# They are not - and the problem will be addressed below.
# Kill Frankenstein to clean workspace
rm(i, wname, rname, FBSites)

# 2.3 Deal with those sites that appear in the BirdNet processing but 
# are nowhere in the pre-existing list of Spots from object SITES
#
# a) Create an object equivalent to 'shortSITES', but with only the
# BirdNET-detected sites
sSITESbn <- SITES[which(SITES$Spot %in% levels(BNoutSpot)), ]
# b) Find out which SITES mentioned in the BirdNet output are not in sSITESbn
lostSITES <- as.character(unique(BNoutSpot[which(!BNoutSpot%in%SITES$Spot)]))
# lostSITES contains "CLc0013x" "CLc0014x" "CLc0015x"
# Deal with the lost sites 
# These site was only named in audio files from 2013. They all have related sites
# (named CLc0013O, CLc0014O, CLc0015N) named in 2011
# Solution: replace the lost 2013 names with the related 2011 names
levels(BNoutSpot)[which(levels(BNoutSpot)=="CLc0013x")] <- "CLc0013O"
levels(BNoutSpot)[which(levels(BNoutSpot)=="CLc0014x")] <- "CLc0014O"
levels(BNoutSpot)[which(levels(BNoutSpot)=="CLc0015x")] <- "CLc0015N"

# 2.4 Create a vector for corresponding detection dates
TempDate <- as.factor(substr(BNout$FileName, 16, 25))

# Standardize date format
levels(TempDate) <- gsub("^_+|_+$", "", levels(TempDate))  # Remove leading/trailing underscores
levels(TempDate) <- gsub("-", "_", levels(TempDate))       # Replace dashes with underscores
levels(TempDate) <- gsub("_0$", "", levels(TempDate))      # Remove "_0" suffix if exists

# Format "YYYYMMDD" to "YYYY_MM_DD"
levels(TempDate) <- ifelse(grepl("^\\d{8}$", levels(TempDate)), 
                         sub("(\\d{4})(\\d{2})(\\d{2})", "\\1_\\2_\\3", levels(TempDate)), 
                         levels(TempDate))

# Organize date vector and replace underscores with dashes
BNoutDates <- gsub("_", "-", as.character(TempDate))
rm(TempDate)

# ----
# 3. Site Clustering
# Define maximum distance in meters necessary to consider two sites as belonging
# to the same cluster (group).
dmax <- 110

# Obtain groups with the "single" method, which clusterizes based on the
# maximum distance between sites from each cluster and saves this groups in a
# dendrogram object:
dendrog <- hclust(dist(sSITESbn[, 3:4]), method = "single")

# Plot dendrogram:
plot(dendrog, hang = -1)

# Create vector with groups of sites based on the dendrogram:
groups <- cutree(dendrog, h = dmax)

# Plot representing groups with colors:
plot(sSITESbn[, 3:4], col = groups)

# Create matrix of groups' centroid coordinates:
CORCEN <- matrix(NA, nrow = max(groups), ncol = 2)
for (i in 1:max(groups)) {
  CORCEN[i, 1] <- mean(sSITESbn[which(groups == i), 3])
  CORCEN[i, 2] <- mean(sSITESbn[which(groups == i), 4])
}

# Add centroids to the plot:
points(CORCEN, pch = 4, cex = 0.5)


# ----
# 4. Review site names for 'ValidAudioFiles' and group membership
# Get the duration of audio files in minutes
fmints <- strptime(ValidAudioFiles$duration, format = "%H:%M:%S") # Turn it into a POSIX object
fmints <- floor(as.numeric(format(fmints, "%H")) * 60 + # hours to min
                               as.numeric(format(fmints, "%M"))) # min

# Clean site names from file names
# The goal here is to make sure that the same sites get the same name in 
# 'ValidAudioFiles' and in 'BNoutSpot'
fnames <- sub("_.*", "", ValidAudioFiles$name)
cleaned_sites <- sapply(fnames, CleanSiteName)

# Apply manual site name corrections
for (i in 1:nrow(corrections)) {
   cleaned_sites[cleaned_sites == corrections$original[i]] <- corrections$corrected[i]
}
# For the 2013 sites that originally had their names created in 2011
cleaned_sites[which(cleaned_sites=="CLc0013x")] <- "CLc0013O"
cleaned_sites[which(cleaned_sites=="CLc0014x")] <- "CLc0014O"
cleaned_sites[which(cleaned_sites=="CLc0015x")] <- "CLc0015N"

# Create a column for corresponding sites (related to groups)
vafi_sites <- cleaned_sites

# List unique site names from the list of valid audio files, regardless of 
# detections
# BNoutSpot <- as.character(BNoutSpot)
sitel <- sort(unique(cleaned_sites))

# Number of unique sites and groups
nsites <- length(sitel)
ngroups <- length(unique(groups))

# Create a list of site vectors for each group
sitegroup <- as.vector(rep(list(rep(NA, max(table(groups)))), ngroups))
for (i in 1:length(sitegroup)) { 
  sitegroup[[i]] <- sSITESbn[which(groups==i),1] 
}

# Group names: use the first site name of each group
groupl <- rep(NA,ngroups)
for(i in 1:ngroups) {
  groupl[i]<-sitegroup[[i]][1]
}
names(sitegroup) <- groupl

# Assign group names to the cleaned sites
file_group_name <- sapply(cleaned_sites, function(site) {
  if (site %in% unlist(sitegroup)) names(sitegroup)[grep(site,sitegroup)]
  else NA
})

# Create a column for corresponding detection dates
vafiDT <- as.factor(substr(ValidAudioFiles$name, 16, 25))

# Standardize date format
levels(vafiDT) <- gsub("^_+|_+$", "", levels(vafiDT))  # Remove leading/trailing underscores
levels(vafiDT) <- gsub("-", "_", levels(vafiDT))       # Replace dashes with underscores
levels(vafiDT) <- gsub("_0$", "", levels(vafiDT))      # Remove "_0" suffix if exists

# Format "YYYYMMDD" to "YYYY_MM_DD"
levels(vafiDT) <- ifelse(grepl("^\\d{8}$", levels(vafiDT)), 
                           sub("(\\d{4})(\\d{2})(\\d{2})", "\\1_\\2_\\3", levels(vafiDT)), 
                           levels(vafiDT))

# Organize date vector and replace underscores with dashes
vafi_dates <- gsub("_", "-", as.character(vafiDT))
rm(vafiDT)

# Final data frame
ValidAudioFilesInfo <- data.frame(
  file_name = ValidAudioFiles$name,
  group = file_group_name,
  site = vafi_sites,
  date = vafi_dates,
  duration_min = fmints,
  stringsAsFactors = FALSE
)

rm(ValidAudioFiles)

# ------------------------------------------------------------------------------
# 5. Workspace Cleanup and Save
## Save output

finalOut <- list(BNoutClean=BNout,Dates=BNoutDates,Spot=BNoutSpot,
                 groups=sitegroup,vafi=ValidAudioFilesInfo)
saveRDS(finalOut,file="CleanedData/BNDetDatesSpotGroupFile.rds")

# Clean work space:
rm(list=ls())


