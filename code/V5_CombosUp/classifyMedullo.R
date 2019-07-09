##################################################
#Purpose: Code to build medulloblastoma classifier
#Author: Pichai Raman
#Date: 5/22/2019
##################################################


#Libraries
library("tidyverse");
library("preprocessCore");
library("cluster")
library("Rtsne")
library("pheatmap")
library("GSVA");
library("genalg");
library("ggplot2");
library("caret")


source("runLimma.R");
source("calcScore.R");
source("classifyGSE109401.R");
source("classifyGSE37418.R");
source("classifyGSE85217.R");
source("PullFeaturesGSE37418.R");

#Read data
mtData <- readRDS("../../data/Medullo_MikeTaylor_Genes.RDS");
sampAnnot <- unique(mtData[,c("Sample", "Subgroup")]);
rownames(sampAnnot) <- sampAnnot[,1];

#Read and Transform 
keepCols <- c("gene_id", "gene_symbol", "FPKM", "Sample")
expDataF <- mtData[,keepCols]
expDataFts <- spread(expDataF, key="Sample", value="FPKM")
geneAnnot <- expDataFts[1:2]
rownames(geneAnnot) <- geneAnnot[,1]
rownames(expDataFts) <- expDataFts[,1]
expDataFts <- data.frame(expDataFts[-1:-2])
dim(expDataFts) # 60498 rows

####################
#Filter Matrix
######################
#Remove genes that have less than 20 FPKM max
expDataFts <- expDataFts[apply(expDataFts, FUN=max, MARGIN=1)>20,] # 24069 genes

#Remove low CV
myCV <- function(x) { mean(x)/sd(x); }
allCVs <- log2(apply(expDataFts, FUN=myCV, MARGIN=1))
allCVs <- (allCVs-mean(allCVs))/sd(allCVs)
expDataFts <- expDataFts[allCVs>(-1),]  # 20153 genes

#Take max value and use to have one gene per
expDataFts[,"Max"] <- apply(expDataFts, FUN=max, MARGIN=1)
expDataFts <- expDataFts[order(-expDataFts[,"Max"]),]
expDataFts[,"GeneName"] <- geneAnnot[rownames(expDataFts),2]
expDataFts <- expDataFts[!duplicated(expDataFts[,"GeneName"]),]
rownames(expDataFts)<- expDataFts[,"GeneName"]
expDataFts <- expDataFts[1:(ncol(expDataFts)-2)] # 19871 rows

#Remove MT- and RPL and RPS genes
expDataFts <- expDataFts[!grepl("^MT-", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("Metazoa_SRP", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("^RPS", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("^RPL", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("^SNORD", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("-", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("\\.", rownames(expDataFts)),]

###################################
#Filter based on DE
###################################
expDataFts_QN <- normalize.quantiles(as.matrix(log2(expDataFts+1)))
rownames(expDataFts_QN) <- rownames(expDataFts);
colnames(expDataFts_QN) <- colnames(expDataFts);
expDataFts_QN <- data.frame(expDataFts_QN)
colnames(expDataFts_QN) <- rownames(sampAnnot)
output <- runLimma(sampAnnot[2], 
	cont=c("SHH-Group4", "SHH-Group3", "SHH-WNT",
		"WNT-Group4", "WNT-Group3", "WNT-SHH",
		"Group4-WNT", "Group4-SHH", "Group4-Group3",
		"Group3-WNT", "Group3-SHH", "Group3-Group4"),expDataFts_QN)

getGenesUp <- function(x, pvalcutoff=0.05, lfccutoff=log2(2))
{
	tmpTable <- topTable(output[[1]], x, 20000)
	tmpTable <- tmpTable[tmpTable[,"adj.P.Val"]<pvalcutoff,]
	tmpTable <- tmpTable[tmpTable[,"logFC"]>lfccutoff,]
	return(rownames(tmpTable));
}

getGenesDown <- function(x, pvalcutoff=0.05, lfccutoff=log2(2))
{
	tmpTable<- topTable(output[[1]], x, 20000)
	tmpTable <- tmpTable[tmpTable[,"adj.P.Val"]<pvalcutoff,]
	tmpTable <- tmpTable[tmpTable[,"logFC"]<((-1)*lfccutoff),]
	return(rownames(tmpTable));
}

getGenesUpTop <- function(x, pvalcutoff=0.05, top=250)
{
	tmpTable <- topTable(output[[1]], x, 20000)
	tmpTable <- tmpTable[tmpTable[,"logFC"]>0,]
	tmpTable <- tmpTable[order(tmpTable[,"adj.P.Val"]),][1:top,]
	return(rownames(tmpTable));
}

getGenesDownTop <- function(x, pvalcutoff=0.05, top=250)
{
	tmpTable<- topTable(output[[1]], x, 20000)
	tmpTable <- tmpTable[tmpTable[,"logFC"]<0,]
	tmpTable <- tmpTable[order(tmpTable[,"adj.P.Val"]),][1:top,]
	return(rownames(tmpTable));
}

###################################
#Get all Differentially expressed genes
###################################

#Genes that are DE in one vs all i.e. intersection
shhGenesUpInt <- Reduce(intersect, list(getGenesUp(1),getGenesUp(2),getGenesUp(3)))
WNTGenesUpInt <- Reduce(intersect, list(getGenesUp(4),getGenesUp(5),getGenesUp(6)))
g4GenesUpInt <- Reduce(intersect, list(getGenesUp(7, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesUp(8, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesUp(9, pvalcutoff=0.05, lfccutoff=log2(1.5))))
g3GenesUpInt <- Reduce(intersect, list(getGenesUp(10, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesUp(11, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesUp(12, pvalcutoff=0.05, lfccutoff=log2(1.5))))

shhGenesDownInt <- Reduce(intersect, list(getGenesDown(1),getGenesDown(2),getGenesDown(3)))
WNTGenesDownInt <- Reduce(intersect, list(getGenesDown(4),getGenesDown(5),getGenesDown(6)))
g4GenesDownInt <- Reduce(intersect, list(getGenesDown(7, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesDown(8, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesDown(9, pvalcutoff=0.05, lfccutoff=log2(1.5))))
g3GenesDownInt <- Reduce(intersect, list(getGenesDown(10, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesDown(11, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesDown(12, pvalcutoff=0.05, lfccutoff=log2(1.5))))

#Genes that are DE in one vs each i.e. union
shhGenesUpUnion <- Reduce(union, list(getGenesUpTop(1),getGenesUpTop(2),getGenesUpTop(3)))
WNTGenesUpUnion <- Reduce(union, list(getGenesUpTop(4),getGenesUpTop(5),getGenesUpTop(6)))
g4GenesUpUnion <- Reduce(union, list(getGenesUpTop(7),getGenesUpTop(8),getGenesUpTop(9)))
g3GenesUpUnion <- Reduce(union, list(getGenesUpTop(10),getGenesUpTop(11),getGenesUpTop(12)))

shhGenesDownUnion <- Reduce(union, list(getGenesDownTop(1),getGenesDownTop(2),getGenesDownTop(3)))
WNTGenesDownUnion <- Reduce(union, list(getGenesDownTop(4),getGenesDownTop(5),getGenesDownTop(6)))
g4GenesDownUnion <- Reduce(union, list(getGenesDownTop(7),getGenesDownTop(8),getGenesDownTop(9)))
g3GenesDownUnion <- Reduce(union, list(getGenesDownTop(10),getGenesDownTop(11),getGenesDownTop(12)))

#Merge lists i.e. union & intersection
shhGenesUp <- c(shhGenesUpInt, shhGenesUpUnion)
shhGenesDown <- c(shhGenesDownInt, shhGenesDownUnion)
WNTGenesUp <- c(WNTGenesUpInt, WNTGenesUpUnion)
WNTGenesDown <- c(WNTGenesDownInt, WNTGenesDownUnion)
g4GenesUp <- c(g4GenesUpInt, g4GenesUpUnion)
g4GenesDown <- c(g4GenesDownInt, g4GenesDownUnion)
g3GenesUp <- c(g3GenesUpInt, g4GenesUpUnion)
g3GenesDown <- c(g3GenesDownInt, g3GenesDownUnion)

#Subtype specific Genes
shhGenes <- c(shhGenesUp, shhGenesDown)
wntGenes <- c(WNTGenesUp, WNTGenesDown)
g3Genes <- c(g3GenesUp, g3GenesDown)
g4Genes <- c(g4GenesUp, g4GenesDown)

#All UpGenes
upGenes <- c(shhGenesUp, WNTGenesUp, g4GenesUp, g3GenesUp)

#All DownGenes
downGenes <- c(shhGenesDown, WNTGenesDown, g4GenesDown, g3GenesDown)

#Best genes are up in some subtypes and down in others
bestGenes <- sort(intersect(upGenes, downGenes)); #193 genes

#Heatmap
png("../../results/V5_CombosUp/heatmapBestGenes.png", width=800, height=800, res=150)
pheatmap(expDataFts_QN[bestGenes,], scale="row", show_rownames=F, show_colnames=F, annotation_col=sampAnnot[2],clustering_distance_cols="correlation")
dev.off();

#TSNE with all data
tsneOut <- Rtsne(t(expDataFts_QN), initial_dims=100, perplexity=20, max_iter=1000)
tsneOut <- data.frame(tsneOut$Y, sampAnnot)
png("../../results/V5_CombosUp/scatterPlotTSNE_AllGenes.png", width=800, height=800, res=150)
ggplot(tsneOut, aes(X1, X2, shape=Subgroup, color=Subgroup))+geom_point()+theme_bw()+ggtitle("T-SNE Medulloblastoma All Genes");
dev.off()

#TSNE with all data
tsneOut <- Rtsne(t(expDataFts_QN[bestGenes,]), initial_dims=100, perplexity=20, max_iter=1000)
tsneOut <- data.frame(tsneOut$Y, sampAnnot)
png("../../results/V5_CombosUp/scatterPlotTSNE_NMFGenes.png", width=800, height=800, res=150)
ggplot(tsneOut, aes(X1, X2, shape=Subgroup, color=Subgroup))+geom_point()+theme_bw()+ggtitle("T-SNE Medulloblastoma NMF Genes");
dev.off()


#################################
#Create Gene Ratio
#################################
library(reshape2);
createRatio <- function(x)
{
g1 <- x[1];
g2 <- x[2];
g1g2_ratio <- 2^(expDataFts_QNMat[g1,]-expDataFts_QNMat[g2,])
return(g1g2_ratio)
}

corGenes <- cor(t(expDataFts_QN[bestGenes,]));
corGenes[lower.tri(corGenes)] <- 1
corGenes <- data.frame(melt(corGenes));
corGenes <- corGenes[corGenes[,"value"]<.99,];
print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "));

expDataFts_QNMat <- as.matrix(expDataFts_QN)
geneRatioOut <- apply(corGenes, FUN=createRatio, MARGIN=1)
geneRatioOut <- data.frame(t(geneRatioOut));
rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_");
colnames(geneRatioOut) <- colnames(expDataFts_QN)
print(paste("Gene Ratios created and processing ", nrow(geneRatioOut), "rows", sep=" "));

#Filter to only high and low ratios first
geneRatioOutM <- melt(geneRatioOut);
hist(log2(geneRatioOutM[,2]+1))

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

shhGeneCombosUp <- Reduce(intersect, list(getGenesComboUp(1, 200),getGenesComboUp(2, 200),getGenesComboUp(3, 200)))
WNTGeneCombosUp <- Reduce(intersect, list(getGenesComboUp(4, 200),getGenesComboUp(5, 200),getGenesComboUp(6, 200)))
g4GeneCombosUp <- Reduce(intersect, list(getGenesComboUp(7, 10000),getGenesComboUp(8, 10000),getGenesComboUp(9, 10000)))
g3GeneCombosUp <- Reduce(intersect, list(getGenesComboUp(10, 10000),getGenesComboUp(11, 10000),getGenesComboUp(12, 1000)))

###########################
#G3 and G4 specific Genes
#To get better separation
###########################
g4GeneCombosUpvsG3<- getGenesComboUp(9, 200)
g3GeneCombosUpvsG4<- getGenesComboUp(12, 200)
g4GeneCombosUp <- union(g4GeneCombosUp, g4GeneCombosUpvsG3);
g3GeneCombosUp <- union(g3GeneCombosUp, g3GeneCombosUpvsG4);

#All UpGenes
allGeneCombos <- c(shhGeneCombosUp, WNTGeneCombosUp, g4GeneCombosUp, g3GeneCombosUp)

#####################
#Plots
#####################
#Heatmap
png("../../results/V5_CombosUp/geneRatios_Heatmap.png", width=800, height=800, res=150)
pheatmap(log2(geneRatioOut[allGeneCombos,]+1), scale="row", show_rownames=F, show_colnames=F, annotation_col=sampAnnot[2],clustering_distance_cols="correlation")
dev.off()

#TSNE 
tsneOut <- Rtsne(t(log2(geneRatioOut[allGeneCombos,]+1)), initial_dims=100, perplexity=20, max_iter=500)
tsneOut <- data.frame(tsneOut$Y, sampAnnot)
png("../../results/V5_CombosUp/scatterPlotTSNE_geneRatios.png", width=800, height=800, res=150)
ggplot(tsneOut, aes(X1, X2, shape=Subgroup, color=Subgroup))+geom_point()+theme_bw()+ggtitle("T-SNE Medulloblastoma Genes");
dev.off();

medulloGeneSetsUp <- list("WNT"=WNTGeneCombosUp, "SHH"=shhGeneCombosUp, "Group3"=g3GeneCombosUp, "Group4"=g4GeneCombosUp)

#Add in variables from GSE37418
saveRDS(allGeneCombos, "../../results/V5_CombosUp/bestFeaturesNew.RDS")
saveRDS(medulloGeneSetsUp, "../../results/V5_CombosUp/medulloSetsUp.RDS")

PullFeaturesGSE37418();



medulloGeneSets37418 <- readRDS("../../results/V5_CombosUp/medulloSetsUpGSE37418.RDS")


################################################
#Special Section to interrogate / Examine Genes
################################################


IntWNT <- intersect(medulloGeneSetsUp$WNT, medulloGeneSets37418$WNT);
IntSHH <- intersect(medulloGeneSetsUp$SHH, medulloGeneSets37418$SHH);
IntGroup3 <- intersect(medulloGeneSetsUp$Group3, medulloGeneSets37418$Group3);
IntGroup4 <- intersect(medulloGeneSetsUp$Group4, medulloGeneSets37418$Group4);


getPlot <- function(x)
{
	tmp <- data.frame(t(geneRatioOut[x,]), sampAnnot[,2]);
	colnames(tmp)[ncol(tmp)] <- "Class"
	tmp[,"Sample"] <- rownames(tmp)
	tmp <- gather(tmp, key="Feature", value="Ratio", -Class, -Sample)
	ggplot(tmp, aes(Class, log2(Ratio), fill=Class))+geom_boxplot()+geom_point()+facet_wrap(~Feature)+theme_bw();

}

getPlot(IntSHH) #Keep ATOH1_NEUROG1, NDP_NEUROG1 (ATOH1, NDP, NEUROG1)
getPlot(IntWNT) #Keep AXIN2_DCX, AXIN2_SAYSD1, AXIN2_CITED2 (ATOH1, NDP, NEUROG1)
getPlot(IntGroup3) #Nothing that great at discriminating
getPlot(IntGroup4) #Nothing that great at discriminating

################################################
#Special Section to interrogate / Examine Genes
################################################

medulloGeneSetsUp$WNT <- union(medulloGeneSetsUp$WNT, intersect(rownames(geneRatioOut), medulloGeneSets37418$WNT));
medulloGeneSetsUp$SHH <- union(medulloGeneSetsUp$SHH, intersect(rownames(geneRatioOut), medulloGeneSets37418$SHH));
medulloGeneSetsUp$Group3 <- union(medulloGeneSetsUp$Group3, intersect(rownames(geneRatioOut), medulloGeneSets37418$Group3));
medulloGeneSetsUp$Group4 <- union(medulloGeneSetsUp$Group4, intersect(rownames(geneRatioOut), medulloGeneSets37418$Group4));

myMat <- calcScoreMat(geneRatioOut, medulloGeneSetsUp);
myClass <- colnames(myMat)[max.col(myMat,ties.method="first")]
myScore <- sum(myClass==as.character(sampAnnot[,2]))/97
confusionMatrix(factor(myClass), sampAnnot[,2])
outputMat <- data.frame(myMat, myClass, as.character(sampAnnot[,2]), myClass==as.character(sampAnnot[,2]))

saveRDS(allGeneCombos, "../../results/V5_CombosUp/bestFeaturesNew.RDS")
saveRDS(medulloGeneSetsUp, "../../results/V5_CombosUp/medulloSetsUp.RDS")

class109401 <- classifyGSE109401();
class37418 <- classifyGSE37418();
class85217 <- classifyGSE85217();

print(paste("WO filtering Accuracy is", round(class109401,4), "and", round(class37418,4), "and", round(class85217,4)));
