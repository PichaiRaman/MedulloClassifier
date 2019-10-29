################################
# Purpose: Classify Medullo Samples from GSE85217
# Author: Pichai Raman
# Date: 6/7/2019
################################

# Let's get all the data
# library("GEOquery")
# require("preprocessCore")
# data_85217 <- getGEO('GSE85217',GSEMatrix=TRUE)
# annot_85217 <- pData(phenoData(data_85217[[1]]))
# exprs_85217 <- exprs(data_85217[[1]])
# save.image("data/loadedGSE_85217.RData")

library(pheatmap)
library(GSVA)
library(caret)
library(preprocessCore)
source(calcScore.R)
library(reshape2)
source('R/getGenes.R')
source('R/createRatio.R')

classifyGSE85217 <- function(signatureProbesLoc="data/model/bestFeaturesNew.RDS", medulloGeneSetsUpLoc="data/model/medulloSetsUp.RDS") {

  load("data/loadedGSE_85217.RData")
    
  # Convert expression matrix to gene symbol
  mapping <- read.delim("data/mappingRefseq.txt")
  mapping <- mapping[,c(1,3)]
  mapping <- unique(mapping)
  geneAnnot <- mapping[mapping[,2]!="",]
  # mapping <- mapping[!duplicated(mapping[,1]),]
  
  colnames(geneAnnot) <- c("ID", "HGNC.symbol")
  exprs_85217_Tmp <- normalize.quantiles(as.matrix(exprs_85217))
  rownames(exprs_85217_Tmp) <- gsub("_at", "", rownames(exprs_85217))
  colnames(exprs_85217_Tmp )<- colnames(exprs_85217)
  exprs_85217 <- exprs_85217_Tmp
  exprs_85217 <- data.frame(exprs_85217)
  exprs_85217[,"Max"] <- apply(exprs_85217, FUN=max, MARGIN=1)
  exprs_85217[,"Probe"] <- rownames(exprs_85217)
  exprs_85217 <- merge(exprs_85217, geneAnnot, by.x="Probe", by.y="ID")
  exprs_85217 <- exprs_85217[order(-exprs_85217[,"Max"]),]
  
  exprs_85217 <- exprs_85217[!is.na(exprs_85217[,"HGNC.symbol"]),]
  exprs_85217 <- exprs_85217[!duplicated(exprs_85217[,"HGNC.symbol"]),]
  rownames(exprs_85217) <- exprs_85217[,"HGNC.symbol"]
  exprs_85217 <- exprs_85217[-1]
  exprs_85217 <- exprs_85217[1:(ncol(exprs_85217)-2)]
  
  ################################
  # Now read in signature genes 
  # Filter matrix to signature genes
  # and Create Gene Ratios
  ################################
  signatureProbes <- readRDS(signatureProbesLoc)
  medulloGeneSetsUp <- readRDS(medulloGeneSetsUpLoc)
  
  signatureGenes <- sapply(signatureProbes, FUN=getGenes)
  signatureGenes <- as.character(signatureGenes)
  signatureGenes <- unique(signatureGenes)
  
  # Filter Matrix
  exprs_85217_SG <- exprs_85217[intersect(rownames(exprs_85217), signatureGenes), ]
  
  # Create Ratios
  corGenes <- cor(t(exprs_85217_SG))
  corGenes <- data.frame(melt(corGenes))
  corGenes <- corGenes[corGenes[,"value"]<.99,]
  print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))
  
  exprs_85217 <- as.matrix(exprs_85217)
  geneRatioOut <- apply(corGenes, FUN=createRatio, MARGIN=1)
  geneRatioOut <- data.frame(t(geneRatioOut))
  rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
  colnames(geneRatioOut) <- colnames(exprs_85217)
  print(paste("Gene Ratios created and processing ", nrow(geneRatioOut), "rows", sep=" "))
  
  ################################
  # Filter to signature ratios
  # Create Heatmap
  ################################
  
  geneRatioOut <- geneRatioOut[intersect(rownames(geneRatioOut), signatureProbes),]
  
  medulloGeneSetsUp$WNT <- intersect(medulloGeneSetsUp$WNT, rownames(geneRatioOut))
  medulloGeneSetsUp$SHH <- intersect(medulloGeneSetsUp$SHH, rownames(geneRatioOut))
  medulloGeneSetsUp$Group3 <- intersect(medulloGeneSetsUp$Group3, rownames(geneRatioOut))
  medulloGeneSetsUp$Group4 <- intersect(medulloGeneSetsUp$Group4, rownames(geneRatioOut))
  
  # Update class
  myClassActual <- as.character(annot_85217[,"subgroup:ch1"])
  myClassActual <- gsub("Group 3", "Group3", myClassActual)
  myClassActual <- gsub("Group 4", "Group4", myClassActual)
  print("Classified")
  
  myMat <- calcScoreMat(geneRatioOut, medulloGeneSetsUp)
  myClassPred <- colnames(myMat)[max.col(myMat,ties.method="first")]
  myScore <- sum(myClassPred==myClassActual)/(length(myClassActual))
  
  sampAnnot <- data.frame(myClassPred, myClassActual)
  colnames(sampAnnot) <- c("Pred", "Actual")
  rownames(sampAnnot) <- colnames(geneRatioOut)
  sampAnnot[,"Correct"] <- myClassPred==myClassActual
  print(confusionMatrix(sampAnnot[,1], sampAnnot[,2]))
  return(myScore)
}
