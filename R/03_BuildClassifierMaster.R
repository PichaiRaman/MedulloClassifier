##################################################
# Purpose: Code to build medulloblastoma classifier
# Author: Pichai Raman
# Date: 5/22/2019
##################################################

# Build Models
source("R/01_BuildClassifierRNASeq.R")
source("R/02_BuildClassifierMicroarray.R")

# Get sets
medulloSetsRSQ <- readRDS("data/model/medulloSetsUpRNASeq.RDS")
medulloSetsMIA <- readRDS("data/model/medulloSetsUpGSE37418.RDS")

# Create Merged dataset
medulloGeneSetsUp <- list()
medulloGeneSetsUp$WNT <- union(medulloSetsRSQ$WNT, medulloSetsMIA$WNT)
medulloGeneSetsUp$SHH <- union(medulloSetsRSQ$SHH, medulloSetsMIA$SHH)
medulloGeneSetsUp$Group3 <- union(medulloSetsRSQ$Group3, medulloSetsMIA$Group3)
medulloGeneSetsUp$Group4 <- union(medulloSetsRSQ$Group4, medulloSetsMIA$Group4)
allGeneCombos <- c(medulloGeneSetsUp$WNT, medulloGeneSetsUp$SHH, medulloGeneSetsUp$Group3, medulloGeneSetsUp$Group4)

# Print length of each subtype gene ratio signature
print(paste("The number of gene ratios in the SHH signature is", length(medulloGeneSetsUp$SHH)))
print(paste("The number of gene ratios in the WNT signature is", length(medulloGeneSetsUp$WNT)))
print(paste("The number of gene ratios in the G3 signature is", length(medulloGeneSetsUp$Group3)))
print(paste("The number of gene ratios in the G4 signature is", length(medulloGeneSetsUp$Group4 )))
print(paste("The total number of gene ratios in the model is", length(allGeneCombos)))

# Save for use in other classes
saveRDS(allGeneCombos, "data/model/bestFeaturesNew.RDS")
saveRDS(medulloGeneSetsUp, "data/model/medulloSetsUp.RDS")

# Classify test datasets
# Source scripts
source("R/04_classifyGSE109401.R")
source("R/04_classifyGSE85217.R")
source("R/04_classifyGSE37418.R")
class109401 <- classifyGSE109401()
class37418 <- classifyGSE37418()
class85217 <- classifyGSE85217()






