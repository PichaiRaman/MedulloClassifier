###################################
#Purpose: Code to classify and compare Medullo samples
#Author: Pichai Raman
#Date: 6/11/19
###################################

#Call libraries
library("tidyverse")
library("pheatmap")
library("ggpubr")
source("calcScore.R")

######################################
#Part 1. Get Medullo Data
######################################

#Load data
res <- read.delim("../../data/PBTA_FPKM_mat_GeneSymbol.txt", sep="\t")
rownames(res) <- res[,1];
res <- res[-1];
#Read in Clinical Data and filter to medulloblastoma samples
clinData <- read.delim("../../data/cbttc_data_clinical_sample.txt", stringsAsFactors=F, skip=4);
clinData <- clinData[,c("SPECIMEN_ID", "CANCER_TYPE")]
clinDataMedullo <- clinData[clinData[,"CANCER_TYPE"]=="Medulloblastoma",]
clinDataMedulloSamps <- clinDataMedullo[,1];

getSamps <- function(x) { strsplit(x, split=";")[[1]]}
clinDataMedulloSamps <- sapply(clinDataMedulloSamps, FUN=getSamps)
clinDataMedulloSamps <- as.character(unlist(clinDataMedulloSamps));

resM <- res[,intersect(colnames(res), clinDataMedulloSamps)]

######################################
#Part 2. Classify Samples
######################################

#Now classify samples
signatureProbes <- readRDS("../../results/V5_CombosUp/bestFeaturesNew.RDS");

getGenes <- function(x)
{
	out <- strsplit(x, split="_")
	output <- c(out[[1]])
}

signatureGenes <- sapply(signatureProbes, FUN=getGenes)
signatureGenes <- as.character(signatureGenes);
signatureGenes <- unique(signatureGenes);

#Filter Matrix
resM_SG <- resM[intersect(rownames(resM), signatureGenes), ]

#Create Ratios
createRatio <- function(x)
{
g1 <- x[1];
g2 <- x[2];
g1g2_ratio <- 2^(log2(resM_SG[g1,]+1)-log2(resM_SG[g2,]+1))
return(g1g2_ratio)
}
library(reshape2);
corGenes <- cor(t(resM_SG));
corGenes <- data.frame(melt(corGenes));
corGenes <- corGenes[corGenes[,"value"]<.99,];
print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "));

resM_SG <- as.matrix(resM_SG);
geneRatioOut <- apply(corGenes, FUN=createRatio, MARGIN=1)
geneRatioOut <- data.frame(t(geneRatioOut));
rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_");
colnames(geneRatioOut) <- colnames(resM_SG)
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

#Classify
myMat <- calcScoreMat(geneRatioOut, medulloGeneSetsUp);
myClassPred <- colnames(myMat)[max.col(myMat,ties.method="first")]
myClassPred <- data.frame(myClassPred, "Medullo");
rownames(myClassPred) <- colnames(geneRatioOut)
ClassPred <- myClassPred[order(myClassPred[,1]),]
geneRatioOut <- geneRatioOut[,rownames(ClassPred)]
pheatmap(log2(geneRatioOut),cluster_cols=F, show_rownames=F, show_colnames=F, annotation_col=ClassPred[1],clustering_distance_cols="correlation")
write.table(myClassPred, "MedulloPrediction_V3.txt", sep="\t", row.names=T)
