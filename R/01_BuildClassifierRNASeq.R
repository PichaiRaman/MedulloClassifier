##################################################
# Purpose: Code to build medulloblastoma classifier
# Author: Pichai Raman
# Date: 5/22/2019
##################################################

# Libraries
library(tidyverse)
library(preprocessCore)
library(pheatmap)
library(Rtsne)
library(reshape2)

# Source scripts 
source("R/utils/runLimma.R")
source("R/utils/pubTheme.R")
source("R/utils/createRatio.R")

######################
# 1. Read in RNA-Seq Data and Format Matrix
######################

# Read RNA-Seq data
mtData <- readRDS("data/Medullo_MikeTaylor_Genes.RDS")

# Get the unique samples
sampAnnot <- unique(mtData[,c("Sample", "Subgroup")])
rownames(sampAnnot) <- sampAnnot[,1]

# mtData is Tall-skinny so going to transform to short-wide, keeping only certain columns i.e. removing expected_count, study, and subgroup
keepCols <- c("gene_id", "gene_symbol", "FPKM", "Sample")
expDataF <- mtData[,keepCols]
expDataFts <- tidyr::spread(expDataF, key="Sample", value="FPKM")

# Let's remove the gene_id (Ensembl) and gene_symbol (entrez) because we will want gene symbols as our rows (unique gene symbols)
geneAnnot <- expDataFts[1:2]
rownames(geneAnnot) <- geneAnnot[,1]

# Set row names as gene_id, this will be updated in the filter section to gene symbol
rownames(expDataFts) <- expDataFts[,1]
expDataFts <- data.frame(expDataFts[-1:-2])
print(paste("The total number of genes is", nrow(expDataFts)))

####################
# 2. Inital Filter of RNA-Seq Matrix
######################

# Remove genes that have less than 20 FPKM max
maxFPKMperGene <- apply(log2(expDataFts), FUN=max, MARGIN=1)
png("results/plots/SuppFig1A.png", width=960, height=960, res=150)
hist(maxFPKMperGene, breaks=1000, xlab="Log2 (FPKM)", main="Histogram of Maximum FPKM Per Gene")
abline(v=log(20), col="red", lwd=3, lty=2)
dev.off()
expDataFts <- expDataFts[apply(expDataFts, FUN=max, MARGIN=1)>20,] # 12334 genes
print(paste("After filtering by FPKM, the total number of genes is", nrow(expDataFts)))

# Remove genes with low coefficient of variance
myCV <- function(x) { mean(x)/sd(x) }
allCVs <- log2(apply(expDataFts, FUN=myCV, MARGIN=1))
allCVs <- (allCVs-mean(allCVs))/sd(allCVs)
png("results/plots/SuppFig1B.png", width=960, height=960, res=150)
hist(allCVs, breaks=1000, xlab="Z-score of CVs (Log2 FPKM)", main="Histogram of Standardized CVs per Gene")
abline(v=(-1), col="red", lwd=3, lty=2)
dev.off()
expDataFts <- expDataFts[allCVs>(-1),]  # 10335 genes
print(paste("After filtering by CV, the total number of genes is", nrow(expDataFts)))

# This is the step where we will get one Gene symbol per row and set gene symbols as our rownames
expDataFts[,"Max"] <- apply(expDataFts, FUN=max, MARGIN=1)
expDataFts <- expDataFts[order(-expDataFts[,"Max"]),]
expDataFts[,"GeneName"] <- geneAnnot[rownames(expDataFts),2]
expDataFts <- expDataFts[!duplicated(expDataFts[,"GeneName"]),]
rownames(expDataFts)<- expDataFts[,"GeneName"]
expDataFts <- expDataFts[1:(ncol(expDataFts)-2)] # 19871 rows
print(paste("After getting one gene symbol per row", nrow(expDataFts)))

# Remove MT- and RPL and RPS genes
expDataFts <- expDataFts[!grepl("^MT-", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("Metazoa_SRP", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("^RPS", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("^RPL", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("^SNORD", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("-", rownames(expDataFts)),]
expDataFts <- expDataFts[!grepl("\\.", rownames(expDataFts)),]
print(paste("After removing MT and RPL genes, the total number of genes is", nrow(expDataFts)))


###################################
# 3. Filter of RNASeq based on differential expression of genes
###################################
expDataFts_QN <- preprocessCore::normalize.quantiles(as.matrix(log2(expDataFts+1))) #Quantile normalization first and then set row/col names
rownames(expDataFts_QN) <- rownames(expDataFts)
colnames(expDataFts_QN) <- colnames(expDataFts)
expDataFts_QN <- data.frame(expDataFts_QN)
colnames(expDataFts_QN) <- rownames(sampAnnot)

# Run differential expression analysis using Limma 
output <- runLimma(sampAnnot[2], 
	cont=c("SHH-Group4", "SHH-Group3", "SHH-WNT",
		"WNT-Group4", "WNT-Group3", "WNT-SHH",
		"Group4-WNT", "Group4-SHH", "Group4-Group3",
		"Group3-WNT", "Group3-SHH", "Group3-Group4"), expDataFts_QN)

# Returns the upregulated genes based on FC cutoff
getGenesUp <- function(x, pvalcutoff=0.05, lfccutoff=log2(2)) {
	tmpTable <- topTable(output[[1]], x, 20000)
	tmpTable <- tmpTable[tmpTable[,"adj.P.Val"]<pvalcutoff,]
	tmpTable <- tmpTable[tmpTable[,"logFC"]>lfccutoff,]
	return(rownames(tmpTable))
}

# Returns the down-regulated genes based on FC cutoff
getGenesDown <- function(x, pvalcutoff=0.05, lfccutoff=log2(2)) {
	tmpTable<- topTable(output[[1]], x, 20000)
	tmpTable <- tmpTable[tmpTable[,"adj.P.Val"]<pvalcutoff,]
	tmpTable <- tmpTable[tmpTable[,"logFC"]<((-1)*lfccutoff),]
	return(rownames(tmpTable))
}

# Returns the most 250 (or user specified)  up-regulated genes
getGenesUpTop <- function(x, pvalcutoff=0.05, top=250) {
	tmpTable <- topTable(output[[1]], x, 20000)
	tmpTable <- tmpTable[tmpTable[,"logFC"]>0,]
	tmpTable <- tmpTable[order(tmpTable[,"adj.P.Val"]),][1:top,]
	return(rownames(tmpTable))
}

# Returns the most 250 (or user specified) down-regulated genes
getGenesDownTop <- function(x, pvalcutoff=0.05, top=250) {
	tmpTable<- topTable(output[[1]], x, 20000)
	tmpTable <- tmpTable[tmpTable[,"logFC"]<0,]
	tmpTable <- tmpTable[order(tmpTable[,"adj.P.Val"]),][1:top,]
	return(rownames(tmpTable))
}

###################################
# 4. Get all Differentially expressed genes associated with each subtype
# We are filtering to genes that are specific for each subtype
# Filters differ between subtypes so that we can get at least a few genes for each
###################################

# Genes that are DE in one subtype vs rest - Upregulated & Downregulated & take intersection
shhGenesUpInt <- Reduce(intersect, list(getGenesUp(1),getGenesUp(2),getGenesUp(3)))
WNTGenesUpInt <- Reduce(intersect, list(getGenesUp(4),getGenesUp(5),getGenesUp(6)))
g4GenesUpInt <- Reduce(intersect, list(getGenesUp(7, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesUp(8, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesUp(9, pvalcutoff=0.05, lfccutoff=log2(1.5))))
g3GenesUpInt <- Reduce(intersect, list(getGenesUp(10, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesUp(11, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesUp(12, pvalcutoff=0.05, lfccutoff=log2(1.5))))

shhGenesDownInt <- Reduce(intersect, list(getGenesDown(1),getGenesDown(2),getGenesDown(3)))
WNTGenesDownInt <- Reduce(intersect, list(getGenesDown(4),getGenesDown(5),getGenesDown(6)))
g4GenesDownInt <- Reduce(intersect, list(getGenesDown(7, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesDown(8, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesDown(9, pvalcutoff=0.05, lfccutoff=log2(1.5))))
g3GenesDownInt <- Reduce(intersect, list(getGenesDown(10, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesDown(11, pvalcutoff=0.05, lfccutoff=log2(1.5)),getGenesDown(12, pvalcutoff=0.05, lfccutoff=log2(1.5))))

# Genes the top 250 DE genes in one subtype vs rest and take union
shhGenesUpUnion <- Reduce(union, list(getGenesUpTop(1),getGenesUpTop(2),getGenesUpTop(3)))
WNTGenesUpUnion <- Reduce(union, list(getGenesUpTop(4),getGenesUpTop(5),getGenesUpTop(6)))
g4GenesUpUnion <- Reduce(union, list(getGenesUpTop(7),getGenesUpTop(8),getGenesUpTop(9)))
g3GenesUpUnion <- Reduce(union, list(getGenesUpTop(10),getGenesUpTop(11),getGenesUpTop(12)))

shhGenesDownUnion <- Reduce(union, list(getGenesDownTop(1),getGenesDownTop(2),getGenesDownTop(3)))
WNTGenesDownUnion <- Reduce(union, list(getGenesDownTop(4),getGenesDownTop(5),getGenesDownTop(6)))
g4GenesDownUnion <- Reduce(union, list(getGenesDownTop(7),getGenesDownTop(8),getGenesDownTop(9)))
g3GenesDownUnion <- Reduce(union, list(getGenesDownTop(10),getGenesDownTop(11),getGenesDownTop(12)))

# Merge lists i.e. union & intersection
shhGenesUp <- c(shhGenesUpInt, shhGenesUpUnion)
shhGenesDown <- c(shhGenesDownInt, shhGenesDownUnion)
WNTGenesUp <- c(WNTGenesUpInt, WNTGenesUpUnion)
WNTGenesDown <- c(WNTGenesDownInt, WNTGenesDownUnion)
g4GenesUp <- c(g4GenesUpInt, g4GenesUpUnion)
g4GenesDown <- c(g4GenesDownInt, g4GenesDownUnion)
g3GenesUp <- c(g3GenesUpInt, g4GenesUpUnion)
g3GenesDown <- c(g3GenesDownInt, g3GenesDownUnion)

# Subtype specific Genes
shhGenes <- c(shhGenesUp, shhGenesDown)
wntGenes <- c(WNTGenesUp, WNTGenesDown)
g3Genes <- c(g3GenesUp, g3GenesDown)
g4Genes <- c(g4GenesUp, g4GenesDown)

# All UpGenes
upGenes <- c(shhGenesUp, WNTGenesUp, g4GenesUp, g3GenesUp)

# All DownGenes
downGenes <- c(shhGenesDown, WNTGenesDown, g4GenesDown, g3GenesDown)

# Best genes are up in some subtypes and down in others
bestGenes <- sort(intersect(upGenes, downGenes)) #1399 genes

###################################
# 5. Now print some plots / figures for paper to show how genes can discriminate subtypes
###################################
# Start with a heatmap
# heatmap top genes
png("results/plots/Figure1B.png", width=1000, height=800, res=150)
sampAnnot$Subgroup <- factor(sampAnnot$Subgroup, levels = c("Group3", "Group4", "SHH", "WNT"))
ann_colors = list(
  Subgroup = c(Group3 = "#F8766D", Group4 = "#7CAE00", SHH = "#00BFC4", WNT = "#C77CFF")
)
pheatmap::pheatmap(expDataFts_QN[bestGenes,], 
                   show_rownames=F, show_colnames=F, 
                   annotation_col= sampAnnot[2],
                   annotation_colors =  ann_colors,
                   clustering_distance_cols="correlation")
dev.off()

# TSNE with all data
tsneOut <- Rtsne::Rtsne(t(expDataFts_QN), initial_dims=100, perplexity=20, max_iter=1000)
tsneOut <- data.frame(tsneOut$Y, sampAnnot)
p <- ggplot2::ggplot(tsneOut, aes(X1, X2, shape = Subgroup, color = Subgroup))+
	ggplot2::geom_point(size = 5, alpha = 0.6)+
	ggplot2::theme_bw()+
	ggplot2::ggtitle("T-SNE Medulloblastoma All Genes")+
	theme_Publication()
ggsave(plot = p, filename = "results/plots/scatterPlotTSNE_AllGenes.png", width = 7, height = 6)

# TSNE only using Best Genes data
tsneOut <- Rtsne(t(log2(expDataFts+1)[bestGenes,]), initial_dims=100, perplexity=20, max_iter=1000)
tsneOut <- data.frame(tsneOut$Y, sampAnnot)
ggplot2::ggplot(tsneOut, aes(X1, X2, shape = Subgroup, color = Subgroup))+
	ggplot2::geom_point(size = 5, alpha = 0.6)+
	ggplot2::theme_bw()+
	ggplot2::ggtitle("T-SNE Medulloblastoma Best Genes")+
	theme_Publication()
ggsave(filename = "results/plots/scatterPlotTSNE_TopGenes.png", width = 7, height = 6)

#################################
# 6. Create Gene Ratios and filter 
#################################

# This correlation is just performed to get all pairs of gene, there may be an easier way of course
corGenes <- cor(t(expDataFts_QN[bestGenes,]))
corGenes[lower.tri(corGenes)] <- 1
corGenes <- data.frame(melt(corGenes))
corGenes <- corGenes[corGenes[,"value"]<.99,] #remove when both the same gene or highly correlated
print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))

# Now create all gene ratios, this part will take a while
expDataFts_QNMat <- as.matrix(expDataFts_QN)
geneRatioOut <- apply(corGenes, FUN = function(x) createRatio(exprs = expDataFts_QNMat, x = x), MARGIN=1)
geneRatioOut <- data.frame(t(geneRatioOut))
rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
colnames(geneRatioOut) <- colnames(expDataFts_QN)
print(paste("Gene Ratios created and processing ", nrow(geneRatioOut), "rows", sep=" "))

# Filter to only high and low ratios first
geneRatioOutM <- reshape2::melt(geneRatioOut)

# Plot of gene ratios in dataset
png("results/plots/SuppFig2B.png", width=800, height=800, res=150)
hist(log2(geneRatioOutM[,2]), breaks=1000, main="Histogram of Gene Ratios (Log2)", xlab="Log2 Gene Ratio")
dev.off()

# Now run Limma on all Gene Ratios
outputGR <- runLimma(sampAnnot[2], 
	cont=c("SHH-Group4", "SHH-Group3", "SHH-WNT",
		"WNT-Group4", "WNT-Group3", "WNT-SHH",
		"Group4-WNT", "Group4-SHH", "Group4-Group3",
		"Group3-WNT", "Group3-SHH", "Group3-Group4"),log2(geneRatioOut))

# Pick out salient features from this
getGenesComboUp <- function(x, percTop=0.01, pvalcutoff=0.01, lfccutoff=4) {
	tmpTable <- topTable(outputGR[[1]], x, 5000000)
	myNumGenes <- round(nrow(tmpTable)*percTop)
	tmpTable <- tmpTable[order(tmpTable[,"adj.P.Val"]),][1:myNumGenes,]
	tmpTable <- tmpTable[tmpTable[,"adj.P.Val"]<pvalcutoff,]
	tmpTable <- tmpTable[tmpTable[,"logFC"]>lfccutoff,]
	return(rownames(tmpTable))
}

shhGeneCombosUp <- Reduce(intersect, list(getGenesComboUp(1, 0.005),getGenesComboUp(2, 0.005),getGenesComboUp(3, 0.005)))
WNTGeneCombosUp <- Reduce(intersect, list(getGenesComboUp(4, 0.005),getGenesComboUp(5, 0.005),getGenesComboUp(6, 0.005)))
g4GeneCombosUp <- Reduce(intersect, list(getGenesComboUp(7, 0.01),getGenesComboUp(8, 0.01)))
g3GeneCombosUp <- Reduce(intersect, list(getGenesComboUp(10, 0.005),getGenesComboUp(11, 0.005)))

# G3 and G4 specific Genes To get better separation
g4GeneCombosUpvsG3<- getGenesComboUp(9, 0.001, lfccutoff=5)
g3GeneCombosUpvsG4<- getGenesComboUp(12, 0.001, lfccutoff=5)
g4GeneCombosUp <- union(g4GeneCombosUp, g4GeneCombosUpvsG3)
g3GeneCombosUp <- union(g3GeneCombosUp, g3GeneCombosUpvsG4)

# All the gene combos
allGeneCombos <- c(shhGeneCombosUp, WNTGeneCombosUp, g4GeneCombosUp, g3GeneCombosUp)

# Print length of each subtype gene ratio signature
print(paste("The number of gene ratios in the SHH signature is", length(shhGeneCombosUp)))
print(paste("The number of gene ratios in the WNT signature is", length(WNTGeneCombosUp)))
print(paste("The number of gene ratios in the G3 signature is", length(g3GeneCombosUp)))
print(paste("The number of gene ratios in the G4 signature is", length(g4GeneCombosUp)))
print(paste("The total number of gene ratios in the model is", length(allGeneCombos)))


#####################
# 7. Plots of gene ratios across data
#####################
#Heatmap
png("results/plots/geneRatios_Heatmap.png", width=800, height=800, res=150)
sampAnnot$Subgroup <- factor(sampAnnot$Subgroup, levels = c("Group3", "Group4", "SHH", "WNT"))
ann_colors = list(
  Subgroup = c(Group3 = "#F8766D", Group4 = "#7CAE00", SHH = "#00BFC4", WNT = "#C77CFF")
)
pheatmap::pheatmap(log2(geneRatioOut[allGeneCombos,]), show_rownames=F, 
                   show_colnames=F, annotation_col=sampAnnot[2],
                   annotation_colors = ann_colors,
                   clustering_distance_cols="correlation", 
                   main="Heatmap Gene Ratios")
dev.off()

# TSNE 
tsneOut <- Rtsne(t(log2(geneRatioOut[allGeneCombos,])), initial_dims=200, perplexity=10, max_iter=500)
tsneOut <- data.frame(tsneOut$Y, sampAnnot)
p <- ggplot(tsneOut, aes(X1, X2, shape=Subgroup, color=Subgroup))+
	geom_point(size = 5, alpha = 0.6)+
 	theme_Publication()+
 	ggtitle("T-SNE Medulloblastoma Gene Ratios")
ggsave(plot = p, filename = "results/plots/scatterPlotTSNE_geneRatios.png", width = 7, height = 6)

#####################
# 8. Save model
#####################
medulloGeneSetsUp <- list("WNT"=WNTGeneCombosUp, "SHH"=shhGeneCombosUp, "Group3"=g3GeneCombosUp, "Group4"=g4GeneCombosUp)
saveRDS(allGeneCombos, "data/model/bestFeaturesNewRNASeq.RDS")
saveRDS(medulloGeneSetsUp, "data/model/medulloSetsUpRNASeq.RDS")
saveRDS(list(geneRatioOut, sampAnnot), "data/RNASeqDataForPlot.RDS")
