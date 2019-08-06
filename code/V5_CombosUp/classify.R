################################
#Purpose: Classify Medullo Samples
#Author: Sherjeel Arif
#Date: 8/5/2019
################################

#######################
# Required Libraries
#######################

library("pheatmap")
library("GSVA");
library("caret")
library("preprocessCore");
source("calcScore.R")


classify <- function(dataset = NULL, signatureProbesLoc="bestFeaturesNew.RDS",medulloGeneSetsUpLoc="medulloSetsUp.RDS")
{
############# 
# Initialize
#############
load(paste("../../data/", dataset,".RData", sep=""))

###############
# Process Data
###############
# processData(dataset)

if (dataset == "loadedGSE_109401") {

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
  
	exprs <- exprs_109401
	annot <- annot_109401
	exprs_Tmp <- normalize.quantiles(as.matrix(exprs));
	rownames(exprs_Tmp) <- rownames(exprs);
	colnames(exprs_Tmp )<- colnames(exprs);
	
	symbol = "Gene.Symbol"

	} else if (dataset == "loadedGSE_85217") {
		symbol = "HGNC.Symbol"
		#Convert expression matrix to gene symbol
		mapping <- read.delim("../../data/mappingRefseq.txt")
		mapping <- mapping[,c(1,3)];
		mapping <- unique(mapping);
		geneAnnot <- mapping[mapping[,2]!="",]
		#mapping <- mapping[!duplicated(mapping[,1]),]

		colnames(geneAnnot) <- c("ID", symbol);
		exprs <- exprs_85217
		annot <- annot_85217
		exprs_Tmp <- normalize.quantiles(as.matrix(exprs));
		rownames(exprs_Tmp) <- gsub("_at", "", rownames(exprs));
		colnames(exprs_Tmp )<- colnames(exprs);
		
		


	} else if (dataset == "loadedGSE_37418") {
		symbol = "Gene.Symbol"
		geneAnnot <- read.delim("../../data/GPL570-55999.txt", skip=16);
		geneAnnot <- geneAnnot[,c("ID", symbol)];
		exprs <- exprs_37418
		annot <- annot_37418
		exprs_Tmp <- normalize.quantiles(as.matrix(exprs));
		rownames(exprs_Tmp) <- rownames(exprs);
		colnames(exprs_Tmp )<- colnames(exprs);
		
		# Update Class
		# myClassActual <- as.character(annot[,43]);
		# myClassActual <- gsub("G3", "Group3", myClassActual)
		# myClassActual <- gsub("G4", "Group4", myClassActual)
		# myClassActual <- gsub("SHH OUTLIER", "SHH", myClassActual)
		# print("Classified");

	} else {
		print("Invalid Dataset Entered")
	}


exprs <- exprs_Tmp
exprs <- data.frame(exprs);
exprs[,"Max"] <- apply(exprs, FUN=max, MARGIN=1)
exprs[,"Probe"] <- rownames(exprs);
exprs <- merge(exprs, geneAnnot, by.x="Probe", by.y="ID");
exprs <- exprs[order(-exprs[,"Max"]),]

exprs <- exprs[!is.na(exprs[,symbol]),]
exprs <- exprs[!duplicated(exprs[,symbol]),];
rownames(exprs) <- exprs[,symbol];
exprs <- exprs[-1];
exprs <- exprs[1:(ncol(exprs)-2)]



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
exprs_SG <- exprs[intersect(rownames(exprs), signatureGenes), ]

#Create Ratios
createRatio <- function(x)
{
g1 <- x[1];
g2 <- x[2];
g1g2_ratio <- 2^(exprs[g1,]-exprs[g2,])
return(g1g2_ratio)
}
library(reshape2);
corGenes <- cor(t(exprs_SG));
corGenes <- data.frame(melt(corGenes));
corGenes <- corGenes[corGenes[,"value"]<.99,];
print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "));

exprs <- as.matrix(exprs);
geneRatioOut <- apply(corGenes, FUN=createRatio, MARGIN=1)
geneRatioOut <- data.frame(t(geneRatioOut));
rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_");
colnames(geneRatioOut) <- colnames(exprs)
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


if (dataset == "loadedGSE_85217" | dataset == "loadedGSE_109401") {
  #Update class
  myClassActual <- as.character(annot[,"subgroup:ch1"]);
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
  
} else {
  #Update class
  myClassActual <- as.character(annot[,43]);
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





}


