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
source("runLimma.R")



PullFeaturesGSE37418 <- function()
{

load("../../data/loadedGSE_37418.RData")

#Convert expression matrix to gene symbol
geneAnnot <- read.delim("../../data/GPL570-55999.txt", skip=16);
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

signatureProbes <- readRDS("../../results/V5_CombosUp/bestFeaturesNew.RDS")
medulloGeneSetsUp <- readRDS("../../results/V5_CombosUp/medulloSetsUp.RDS")

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

#geneRatioOut <- geneRatioOut[intersect(rownames(geneRatioOut), signatureProbes),]

sampAnnot <- annot_37418[,c("Sex:ch1", "subgroup:ch1")]
sampAnnot[,2] <- gsub("G3", "Group3", sampAnnot[,2])
sampAnnot[,2] <- gsub("G4", "Group4", sampAnnot[,2])
sampAnnot[,2] <- gsub("SHH OUTLIER", "SHH", sampAnnot[,2])

outputGR <- runLimma(sampAnnot[2], 
	cont=c("SHH-Group4", "SHH-Group3", "SHH-WNT",
		"WNT-Group4", "WNT-Group3", "WNT-SHH",
		"Group4-WNT", "Group4-SHH", "Group4-Group3",
		"Group3-WNT", "Group3-SHH", "Group3-Group4"),log2(geneRatioOut))


getGenesComboUp <- function(x, myNumGenes=250, pvalcutoff=0.01, lfccutoff=2)
{
	tmpTable <- topTable(outputGR[[1]], x, 50000)
	tmpTable <- tmpTable[tmpTable[,"adj.P.Val"]<pvalcutoff,]
	tmpTable <- tmpTable[tmpTable[,"logFC"]>lfccutoff,]
	tmpTable <- tmpTable[order(tmpTable[,"adj.P.Val"]),][1:myNumGenes,]
	return(rownames(tmpTable));
}

getGenesComboDown <- function(x, myNumGenes=250, pvalcutoff=0.01, lfccutoff=2)
{
	tmpTable<- topTable(outputGR[[1]], x, 50000)
	tmpTable <- tmpTable[tmpTable[,"adj.P.Val"]<pvalcutoff,]
	tmpTable <- tmpTable[tmpTable[,"logFC"]<((-1)*lfccutoff),]
	tmpTable <- tmpTable[order(tmpTable[,"adj.P.Val"]),][1:myNumGenes,]
	return(rownames(tmpTable));
}
#Had 500 and 5000 before

shhGeneCombosUp <- Reduce(intersect, list(getGenesComboUp(1, 200),getGenesComboUp(2, 200),getGenesComboUp(3, 200)))
WNTGeneCombosUp <- Reduce(intersect, list(getGenesComboUp(4, 200),getGenesComboUp(5, 200),getGenesComboUp(6, 200)))
g4GeneCombosUp <- Reduce(intersect, list(getGenesComboUp(7, 3000),getGenesComboUp(8, 3000),getGenesComboUp(9, 3000)))
g3GeneCombosUp <- Reduce(intersect, list(getGenesComboUp(10, 3000),getGenesComboUp(11, 3000),getGenesComboUp(12, 3000)))

medulloGeneSetsUp$WNT <- intersect(WNTGeneCombosUp, rownames(geneRatioOut))
medulloGeneSetsUp$SHH <- intersect(shhGeneCombosUp, rownames(geneRatioOut))
medulloGeneSetsUp$Group3 <- intersect(g3GeneCombosUp, rownames(geneRatioOut))
medulloGeneSetsUp$Group4 <- intersect(g4GeneCombosUp, rownames(geneRatioOut))

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

saveRDS(medulloGeneSetsUp, "../../results/V5_CombosUp/medulloSetsUpGSE37418.RDS")

}















