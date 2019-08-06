# TEST SCRIPT
source("classifyGSE37418.R")
source("classifyGSE85217.R")
x <- classifyGSE85217(signatureProbesLoc="bestFeaturesNew.RDS",medulloGeneSetsUpLoc="medulloSetsUp.RDS")
source("classify.R")
classify(dataset = "loadedGSE_37418", signatureProbesLoc="bestFeaturesNew.RDS",medulloGeneSetsUpLoc="medulloSetsUp.RDS")
classify(dataset = "loadedGSE_85217", signatureProbesLoc="bestFeaturesNew.RDS",medulloGeneSetsUpLoc="medulloSetsUp.RDS")
classify(dataset = "loadedGSE_109401", signatureProbesLoc="bestFeaturesNew.RDS",medulloGeneSetsUpLoc="medulloSetsUp.RDS")
