################################
#Purpose: Classify Medullo Samples from GSE37418
#Author: Pichai Raman
#Date: 6/7/2019
################################



#Let's get all the data
#library("GEOquery");
#require("preprocessCore");
#data_37418 <- getGEO('GSE37418',GSEMatrix=TRUE);
#annot_37418 <- pData(phenoData(data_37418[[1]]));
#exprs_37418 <- exprs(data_37418[[1]]);
#save.image("loadedGSE_37418.RData")
library("pheatmap")
library("GSVA");
library("caret")
library("preprocessCore");
source("calcScore.R")



classifyGSE37418 <- function(signatureProbesLoc="../results/model/bestFeaturesNew.RDS", medulloGeneSetsUpLoc="../results/model/medulloSetsUp.RDS")
{
load("../data/loadedGSE_37418.RData")

#Convert expression matrix to gene symbol
geneAnnot <- read.delim("../data/GPL570-55999.txt", skip=16);
geneAnnot <- geneAnnot[,c("ID", "Gene.Symbol")];
exprs_37418_Tmp <- normalize.quantiles(as.matrix(exprs_37418));
rownames(exprs_37418_Tmp) <- rownames(exprs_37418);
colnames(exprs_37418_Tmp )<- colnames(exprs_37418);
exprs_37418 <- exprs_37418_Tmp
exprs_37418 <- data.frame(exprs_37418);
exprs_37418[,"Max"] <- apply(exprs_37418, FUN=max, MARGIN=1)
exprs_37418[,"Probe"] <- rownames(exprs_37418);
exprs_37418 <- merge(exprs_37418, geneAnnot, by.x="Probe", by.y="ID");
exprs_37418 <- exprs_37418[order(-exprs_37418[,"Max"]),]

exprs_37418 <- exprs_37418[exprs_37418[,"Gene.Symbol"]!="",]
exprs_37418 <- exprs_37418[!duplicated(exprs_37418[,"Gene.Symbol"]),];
rownames(exprs_37418) <- exprs_37418[,"Gene.Symbol"];
exprs_37418 <- exprs_37418[-1];
exprs_37418 <- exprs_37418[1:(ncol(exprs_37418)-2)]

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
exprs_37418_SG <- exprs_37418[intersect(rownames(exprs_37418), signatureGenes), ]

#Create Ratios
createRatio <- function(x)
{
g1 <- x[1];
g2 <- x[2];
g1g2_ratio <- 2^(exprs_37418[g1,]-exprs_37418[g2,])
return(g1g2_ratio)
}
library(reshape2);
corGenes <- cor(t(exprs_37418_SG));
corGenes <- data.frame(melt(corGenes));
corGenes <- corGenes[corGenes[,"value"]<.99,];
print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "));

exprs_37418 <- as.matrix(exprs_37418);
geneRatioOut <- apply(corGenes, FUN=createRatio, MARGIN=1)
geneRatioOut <- data.frame(t(geneRatioOut));
rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_");
colnames(geneRatioOut) <- colnames(exprs_37418)
print(paste("Gene Ratios created and processing", nrow(geneRatioOut), "rows", sep=" "));


################################
#Filter to signature ratios
#Create Heatmap
################################

geneRatioOut <- geneRatioOut[intersect(rownames(geneRatioOut), signatureProbes),]

medulloGeneSetsUp$WNT <- intersect(medulloGeneSetsUp$WNT, rownames(geneRatioOut))
medulloGeneSetsUp$SHH <- intersect(medulloGeneSetsUp$SHH, rownames(geneRatioOut))
medulloGeneSetsUp$Group3 <- intersect(medulloGeneSetsUp$Group3, rownames(geneRatioOut))
medulloGeneSetsUp$Group4 <- intersect(medulloGeneSetsUp$Group4, rownames(geneRatioOut))

#Update class
myClassActual <- as.character(annot_37418[,43]);
myClassActual <- gsub("G3", "Group3", myClassActual)
myClassActual <- gsub("G4", "Group4", myClassActual)
myClassActual <- gsub("SHH OUTLIER", "SHH", myClassActual)
print("Classified");


myMat <- calcScoreMat(geneRatioOut, medulloGeneSetsUp);
myClassPred <- colnames(myMat)[max.col(myMat,ties.method="first")]
myScore <- sum(myClassPred==myClassActual)/(length(myClassActual)-2)

sampAnnot <- data.frame(myClassPred, myClassActual);
colnames(sampAnnot) <- c("Pred", "Actual")
rownames(sampAnnot) <- colnames(geneRatioOut)
sampAnnot[,"Correct"] <- myClassPred==myClassActual
#pheatmap(geneRatioOut ,show_rownames=F, show_colnames=F, annotation_col=sampAnnot,clustering_distance_cols="correlation")

sampAnnot <- sampAnnot[sampAnnot[,2]!="U",]
sampAnnot[,2] <- factor(sampAnnot[,2], levels=c("Group3", "Group4", "WNT", "SHH"))
sampAnnot[,1] <- factor(sampAnnot[,1], levels=c("Group3", "Group4", "WNT", "SHH"))

print(confusionMatrix(sampAnnot[,1], sampAnnot[,2]));
return(myScore)
}

















