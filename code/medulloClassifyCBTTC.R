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
load("../../data/cbttc_genes_fpkm_1110.RData")

#Read in Clinical Data and Mapping file
clinData <- read.delim("../../data/study_view_clinical_data.txt", stringsAsFactors=F);
mapping <- read.delim("../../data/mappingFile.txt", header=F, stringsAsFactors=F);

#Format clinical Data
clinData <- merge(clinData, mapping, by.x="Sample.ID", by.y="V1")
clinDataOrig <- clinData;
clinData <- clinData[!grepl("CL", clinData[,"CBTTC_PAIRED_IDS"]),]
clinData <- clinData[,c("V2", "Cancer.Type", "Cancer.Type.Detailed", "TUMOR_TISSUE_SITE")]
clinData <- unique(clinData);
clinData <- clinData[!duplicated(clinData[,"V2"]),]
rownames(clinData)<- clinData[,"V2"]

#Format expression data
rownames(res) <- res[,1];
res[,2] <- as.character(res[,2]);
res <- res[!grepl("-", res[,2]),]
res <- res[!grepl("\\.", res[,2]),]
res <- res[!grepl("_", res[,2]),]

#Now take gene with max value
res[,"max"] <- apply(res[3:ncol(res)], FUN=max, MARGIN=1)
res <- res[order(-res[,"max"]),]
res <- res[!duplicated(res[,2]),]
rownames(res) <- res[,2];
res <- res[-1:-2];
res <- res[-ncol(res)]; #Remove max

#Format so Clinical Data and Res have same ordering
sampOrder <- intersect(rownames(clinData), colnames(res))
clinData <- clinData[sampOrder,];
res <- res[,sampOrder];

#Now get it down to only Medullo's
clinDataM <- clinData[clinData[,"Cancer.Type"]=="Medulloblastoma",]
resM <- res[,rownames(clinDataM)]

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
write.table(myClassPred, "MedulloPrediction_V2.txt", sep="\t", row.names=T)
