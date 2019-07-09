################################
#Purpose: Classify Medullo Samples from GSE109401
#Author: Pichai Raman
#Date: 6/7/2019
################################



#Let's get all the data
#library("GEOquery");
#require("preprocessCore");
#data_109401 <- getGEO('GSE109401',GSEMatrix=TRUE);
#annot_109401 <- pData(phenoData(data_109401[[1]]));
#exprs_109401 <- exprs(data_109401[[1]]);
#annot_109401 <- annot_109401[annot_109401[,"subgroup:ch1"]!="Cerebellum",]
#exprs_109401 <- exprs_109401[,rownames(annot_109401)]
#save.image("loadedGSE_109401.RData")
library("pheatmap")
library("GSVA");
library("caret")
library("preprocessCore");
source("calcScore.R")




classifyGSE109401 <- function(signatureProbesLoc="../../results/V5_CombosUp/bestFeaturesNew.RDS", medulloGeneSetsUpLoc="../../results/V5_CombosUp/medulloSetsUp.RDS")
{

load("../../data/loadedGSE_109401.RData")


#Convert expression matrix to gene symbol
geneAnnot <- read.delim("../../data/GPL16686.txt", skip=8);
mapping <- read.delim("../../data/mappingRefseq.txt")
mapping <- mapping[3:4];
mapping <- unique(mapping);
mapping <- mapping[mapping[,2]!="",]
mapping <- mapping[!duplicated(mapping[,2]),]
geneAnnot <- merge(geneAnnot, mapping, by.x="GB_ACC", by.y="RefSeq.mRNA.ID", all.x=T)

geneAnnot <- geneAnnot[,c("ID", "HGNC.symbol")];
colnames(geneAnnot)[2] <- "Gene.Symbol";
exprs_109401_Tmp <- normalize.quantiles(as.matrix(exprs_109401));
rownames(exprs_109401_Tmp) <- rownames(exprs_109401);
colnames(exprs_109401_Tmp )<- colnames(exprs_109401);
exprs_109401 <- exprs_109401_Tmp
exprs_109401 <- data.frame(exprs_109401);
exprs_109401[,"Max"] <- apply(exprs_109401, FUN=max, MARGIN=1)
exprs_109401[,"Probe"] <- rownames(exprs_109401);
exprs_109401 <- merge(exprs_109401, geneAnnot, by.x="Probe", by.y="ID");
exprs_109401 <- exprs_109401[order(-exprs_109401[,"Max"]),]

exprs_109401 <- exprs_109401[!is.na(exprs_109401[,"Gene.Symbol"]),]
exprs_109401 <- exprs_109401[!duplicated(exprs_109401[,"Gene.Symbol"]),];
rownames(exprs_109401) <- exprs_109401[,"Gene.Symbol"];
exprs_109401 <- exprs_109401[-1];
exprs_109401 <- exprs_109401[1:(ncol(exprs_109401)-2)]

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
exprs_109401_SG <- exprs_109401[intersect(rownames(exprs_109401), signatureGenes), ]

#Create Ratios
createRatio <- function(x)
{
g1 <- x[1];
g2 <- x[2];
g1g2_ratio <- 2^(exprs_109401[g1,]-exprs_109401[g2,])
return(g1g2_ratio)
}
library(reshape2);
corGenes <- cor(t(exprs_109401_SG));
corGenes <- data.frame(melt(corGenes));
corGenes <- corGenes[corGenes[,"value"]<.99,];
print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "));

exprs_109401 <- as.matrix(exprs_109401);
geneRatioOut <- apply(corGenes, FUN=createRatio, MARGIN=1)
geneRatioOut <- data.frame(t(geneRatioOut));
rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_");
colnames(geneRatioOut) <- colnames(exprs_109401)
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

#Update class
myClassActual <- as.character(annot_109401[,"subgroup:ch1"]);
myClassActual <- gsub("Group 3", "Group3", myClassActual)
myClassActual <- gsub("Group 4", "Group4", myClassActual)
print("Classified");


myMat <- calcScoreMat(geneRatioOut, medulloGeneSetsUp);
myClassPred <- colnames(myMat)[max.col(myMat,ties.method="first")]
myScore <- sum(myClassPred==myClassActual)/(length(myClassActual))

sampAnnot <- data.frame(myClassPred, myClassActual);
colnames(sampAnnot) <- c("Pred", "Actual")
rownames(sampAnnot) <- colnames(geneRatioOut)
sampAnnot[,"Correct"] <- myClassPred==myClassActual
print(confusionMatrix(sampAnnot[,1], sampAnnot[,2]));
return(myScore)
}


















