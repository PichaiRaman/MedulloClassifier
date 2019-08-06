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



processData <- function(dataset = NULL)
{

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
		exprs_Tmp <- normalize.quantiles(as.matrix(exprs));
		rownames(exprs_Tmp) <- gsub("_at", "", rownames(exprs));
		colnames(exprs_Tmp )<- colnames(exprs);
		
		


	} else if (dataset == "loadedGSE_37418") {
		symbol = "Gene.Symbol"
		geneAnnot <- read.delim("../../data/GPL570-55999.txt", skip=16);
		geneAnnot <- geneAnnot[,c("ID", symbol)];
		exprs_Tmp <- normalize.quantiles(as.matrix(exprs));
		rownames(exprs_Tmp) <- rownames(exprs);
		colnames(exprs_Tmp )<- colnames(exprs);
		


	} else {
		print("Invalid Dataset Entered")
	}