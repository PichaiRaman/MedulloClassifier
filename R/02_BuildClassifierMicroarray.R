################################
# Purpose: Build Classifier from microarray dataset GSE37418
# Author: Pichai Raman
# Date: 6/7/2019
################################

# Load libraries
library(reshape2)
library(preprocessCore)
library(GEOquery)

# Source Functions
source("R/utils/runLimma.R")
source("R/utils/pubTheme.R")
source("R/utils/getGenes.R")
source("R/utils/createRatio.R")

################################
# 1. Load and format data to one gene 
# symbol per row
################################

# Load data
if(file.exists('data/loadedGSE_37418.RData')){
  load("data/loadedGSE_37418.RData")
} else {
  data_37418 <- getGEO('GSE37418',GSEMatrix=TRUE)
  annot_37418 <- pData(phenoData(data_37418[[1]]))
  exprs_37418 <- exprs(data_37418[[1]])
  save.image("data/loadedGSE_37418.RData")
}

# Convert expression matrix to gene symbol
geneAnnot <- read.delim("data/GPL570-55999.txt", skip=16)
geneAnnot <- geneAnnot[,c("ID", "Gene.Symbol")]
exprs_37418_Tmp <- normalize.quantiles(as.matrix(exprs_37418))
rownames(exprs_37418_Tmp) <- rownames(exprs_37418)
colnames(exprs_37418_Tmp )<- colnames(exprs_37418)
exprs_37418 <- exprs_37418_Tmp
exprs_37418 <- data.frame(exprs_37418)
exprs_37418[,"Max"] <- apply(exprs_37418, FUN=max, MARGIN=1)
exprs_37418[,"Probe"] <- rownames(exprs_37418)
exprs_37418 <- merge(exprs_37418, geneAnnot, by.x="Probe", by.y="ID")
exprs_37418 <- exprs_37418[order(-exprs_37418[,"Max"]),]

# Remove rows without gene symbols & make sure one symbol per row
exprs_37418 <- exprs_37418[exprs_37418[,"Gene.Symbol"]!="",]
exprs_37418 <- exprs_37418[!duplicated(exprs_37418[,"Gene.Symbol"]),]
rownames(exprs_37418) <- exprs_37418[,"Gene.Symbol"]
exprs_37418 <- exprs_37418[-1]
exprs_37418 <- exprs_37418[1:(ncol(exprs_37418)-2)]

################################
# 2. Now read in signature genes 
# Filter matrix to signature genes
# and Create Gene Ratios
################################

# Read in signature probes
signatureProbes <- readRDS("data/model/bestFeaturesNewRNASeq.RDS")
medulloGeneSetsUp <- readRDS("data/model/medulloSetsUpRNASeq.RDS")

# Pull out signature genes
signatureGenes <- sapply(signatureProbes, FUN=getGenes)
signatureGenes <- as.character(signatureGenes)
signatureGenes <- unique(signatureGenes)

# Filter Matrix
exprs_37418_SG <- exprs_37418[intersect(rownames(exprs_37418), signatureGenes), ]

# Create Ratios
# This correlation is just performed to get all pairs of gene, there may be an easier way of course
corGenes <- cor(t(exprs_37418_SG))
corGenes <- data.frame(reshape2::melt(corGenes))
corGenes <- corGenes[corGenes[,"value"]<.99,]
print(paste("Cor Matrix Created and processing", nrow(corGenes), "rows", sep=" "))

# Now create all gene ratios, this part will take a while
exprs_37418 <- as.matrix(exprs_37418)
geneRatioOut <- apply(corGenes, FUN = function(x) createRatio(exprs = exprs_37418, x = x), MARGIN=1)
geneRatioOut <- data.frame(t(geneRatioOut))
rownames(geneRatioOut) <- paste(corGenes[,1], corGenes[,2], sep="_")
colnames(geneRatioOut) <- colnames(exprs_37418)
print(paste("Gene Ratios created and processing", nrow(geneRatioOut), "rows", sep=" "))

# Filter to only high and low ratios first
geneRatioOutM <- reshape2::melt(geneRatioOut)

# Plot of gene ratios in dataset
png("results/plots/SuppFig2A.png", width=800, height=800, res=150)
hist(log2(geneRatioOutM[,2]), breaks=1000, main="Histogram of Gene Ratios (Log2)", xlab="Log2 Gene Ratio")
dev.off()

################################
# 3. Filter to signature ratios
################################

# Format samples
sampAnnot <- annot_37418[,c("Sex:ch1", "subgroup:ch1")]
sampAnnot[,2] <- gsub("G3", "Group3", sampAnnot[,2])
sampAnnot[,2] <- gsub("G4", "Group4", sampAnnot[,2])
sampAnnot[,2] <- gsub("SHH OUTLIER", "SHH", sampAnnot[,2])

# Run Limma
outputGR <- runLimma(sampAnnot[2], 
	cont=c("SHH-Group4", "SHH-Group3", "SHH-WNT",
		"WNT-Group4", "WNT-Group3", "WNT-SHH",
		"Group4-WNT", "Group4-SHH", "Group4-Group3",
		"Group3-WNT", "Group3-SHH", "Group3-Group4"),log2(geneRatioOut))

# Function to get gene combos
getGenesComboUp <- function(x, percTop=0.01, pvalcutoff=0.01, lfccutoff=3) {
	tmpTable <- topTable(outputGR[[1]], x, 5000000)
	myNumGenes <- round(nrow(tmpTable)*percTop)
	tmpTable <- tmpTable[order(tmpTable[,"adj.P.Val"]),][1:myNumGenes,]
	tmpTable <- tmpTable[tmpTable[,"adj.P.Val"]<pvalcutoff,]
	tmpTable <- tmpTable[tmpTable[,"logFC"]>lfccutoff,]
	return(rownames(tmpTable))
}

# Create list of features
shhGeneCombosUp <- Reduce(intersect, list(getGenesComboUp(1, 0.005),getGenesComboUp(2, 0.005),getGenesComboUp(3, 0.005)))
WNTGeneCombosUp <- Reduce(intersect, list(getGenesComboUp(4, 0.005),getGenesComboUp(5, 0.005),getGenesComboUp(6, 0.005)))
g4GeneCombosUp <- Reduce(intersect, list(getGenesComboUp(7, 0.01),getGenesComboUp(8, 0.01)))
g3GeneCombosUp <- Reduce(intersect, list(getGenesComboUp(10, 0.01),getGenesComboUp(11, 0.01)))

# G3 and G4 specific Genes To get better separation
g4GeneCombosUpvsG3<- getGenesComboUp(9, 0.005, lfccutoff=3)
g3GeneCombosUpvsG4<- getGenesComboUp(12, 0.005, lfccutoff=3)
g4GeneCombosUp <- union(g4GeneCombosUp, g4GeneCombosUpvsG3)
g3GeneCombosUp <- union(g3GeneCombosUp, g3GeneCombosUpvsG4)

# Put all gene ratios in model
medulloGeneSetsUp$WNT <- intersect(WNTGeneCombosUp, rownames(geneRatioOut))
medulloGeneSetsUp$SHH <- intersect(shhGeneCombosUp, rownames(geneRatioOut))
medulloGeneSetsUp$Group3 <- intersect(g3GeneCombosUp, rownames(geneRatioOut))
medulloGeneSetsUp$Group4 <- intersect(g4GeneCombosUp, rownames(geneRatioOut))

# All the gene combos
allGeneCombos <- c(shhGeneCombosUp, WNTGeneCombosUp, g4GeneCombosUp, g3GeneCombosUp)

# Print length of each subtype gene ratio signature
print(paste("The number of gene ratios in the SHH signature is", length(shhGeneCombosUp)))
print(paste("The number of gene ratios in the WNT signature is", length(WNTGeneCombosUp)))
print(paste("The number of gene ratios in the G3 signature is", length(g3GeneCombosUp)))
print(paste("The number of gene ratios in the G4 signature is", length(g4GeneCombosUp)))
print(paste("The total number of gene ratios in the model is", length(allGeneCombos)))

#####################
# 4. Save model
#####################
saveRDS(medulloGeneSetsUp, "data/model/medulloSetsUpGSE37418.RDS")
