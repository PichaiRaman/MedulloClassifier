################################
#Purpose: Classify Medullo Samples from PDX Samples
#Author: Pichai Raman
#Date: 6/7/2019
################################

library("pheatmap")
library("GSVA");
library("caret")
library("preprocessCore");
source("calcScore.R")




classifyPDX_JLH <- function(signatureProbesLoc="../../results/V5_CombosUp/bestFeaturesNew.RDS", medulloGeneSetsUpLoc="../../results/V5_CombosUp/medulloSetsUp.RDS")
{

load("../../data/2019-07-03-mb-fpkm-mat.rda")

exprsPDX <- mb.mat
rm(mb.mat);
################################
#Now read in signature genes 
#Filter matrix to signature genes
#and Create Gene Ratios
################################
signatureProbes <- readRDS(signatureProbesLoc)
medulloGeneSetsUp <- readRDS(medulloGeneSetsUpLoc)

getGenes <- function(x)
{
	out <- strsplit(x, split="_")
	output <- c(out[[1]])
}

signatureGenes <- sapply(signatureProbes, FUN=getGenes)
signatureGenes <- as.character(signatureGenes);
signatureGenes <- unique(signatureGenes);

#Filter Matrix
exprsPDX <- exprsPDX[intersect(rownames(exprsPDX), signatureGenes), ]

#Create Ratios
createRatio <- function(x)
{
g1 <- x[1];
g2 <- x[2];
g1g2_ratio <- 2^(exprsPDX[g1,]-exprsPDX[g2,])
return(g1g2_ratio)
}
library(reshape2);
corGenes <- cor(t(exprsPDX));
corGenes <- data.frame(melt(corGenes));
corGenes <- corGenes[corGenes[,"value"]<.99,];
print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "));

exprsPDX <- log2(as.matrix(exprsPDX)+1);
geneRatioOut <- apply(corGenes, FUN=createRatio, MARGIN=1)
geneRatioOut <- data.frame(t(geneRatioOut));
rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_");
colnames(geneRatioOut) <- colnames(exprsPDX)
print(paste("Gene Ratios created and processing ", nrow(geneRatioOut), "rows", sep=" "));


################################
#Filter to signature ratios
#Create Heatmap
################################

geneRatioOut <- geneRatioOut[intersect(rownames(geneRatioOut), signatureProbes),]

medulloGeneSetsUp <- readRDS("../../results/V5_CombosUp/medulloSetsUp.RDS");

medulloGeneSetsUp$WNT <- intersect(medulloGeneSetsUp$WNT, rownames(geneRatioOut))
medulloGeneSetsUp$SHH <- intersect(medulloGeneSetsUp$SHH, rownames(geneRatioOut))
medulloGeneSetsUp$Group3 <- intersect(medulloGeneSetsUp$Group3, rownames(geneRatioOut))
medulloGeneSetsUp$Group4 <- intersect(medulloGeneSetsUp$Group4, rownames(geneRatioOut))


myMat <- calcScoreMat(geneRatioOut, medulloGeneSetsUp);
myClassPred <- colnames(myMat)[max.col(myMat,ties.method="first")]
output <- data.frame(myMat, myClassPred)
colnames(output)[5] <- "Prediction"
write.table(output, "PDX_Predictions.txt", sep="\t", row.names=F)
return(output)
}


















