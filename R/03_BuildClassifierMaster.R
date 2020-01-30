##################################################
# Purpose: Code to build medulloblastoma classifier
# Author: Pichai Raman
# Date: 5/22/2019
##################################################

# load libraries
library(ggplot2)
library(Rtsne)
library(tidyverse)
source('R/utils/pubTheme.R')

# build models
source('R/utils/pubTheme.R')
if(!file.exists('data/model/medulloSetsUpRNASeq.RDS')){
  source('R/01_BuildClassifierRNASeq.R')
}
if(!file.exists('data/model/medulloSetsUpGSE37418.RDS')){
  source('R/02_BuildClassifierMicroarray.R')
}

# Get sets
medulloSetsRSQ <- readRDS('data/model/medulloSetsUpRNASeq.RDS')
medulloSetsMIA <- readRDS('data/model/medulloSetsUpGSE37418.RDS')

# Create Merged dataset
medulloGeneSetsUp <- list()
medulloGeneSetsUp$WNT <- union(medulloSetsRSQ$WNT, medulloSetsMIA$WNT)
medulloGeneSetsUp$SHH <- union(medulloSetsRSQ$SHH, medulloSetsMIA$SHH)
medulloGeneSetsUp$Group3 <- union(medulloSetsRSQ$Group3, medulloSetsMIA$Group3)
medulloGeneSetsUp$Group4 <- union(medulloSetsRSQ$Group4, medulloSetsMIA$Group4)
allGeneCombos <- c(medulloGeneSetsUp$WNT, medulloGeneSetsUp$SHH, medulloGeneSetsUp$Group3, medulloGeneSetsUp$Group4)

# Print length of each subtype gene ratio signature
print(paste("The number of gene ratios in the SHH signature is", length(medulloGeneSetsUp$SHH)))
print(paste("The number of gene ratios in the WNT signature is", length(medulloGeneSetsUp$WNT)))
print(paste("The number of gene ratios in the G3 signature is", length(medulloGeneSetsUp$Group3)))
print(paste("The number of gene ratios in the G4 signature is", length(medulloGeneSetsUp$Group4 )))
print(paste("The total number of gene ratios in the model is", length(allGeneCombos)))

# Figure 1C: barplot of GERs per subtype
medulloSigTS <- stack(medulloGeneSetsUp)
colnames(medulloSigTS) <- c("GeneRatio", "Subtype")
write.csv(medulloSigTS, file = 'results/tables/SuppTable2.csv', quote = F, row.names = F)
medulloSigTSCount <- data.frame(table(medulloSigTS[,"Subtype"]))
medulloSigTSCount$Var1 <- factor(medulloSigTSCount$Var1, levels = c("Group3", "Group4", "SHH", "WNT", "Unknown"))
colnames(medulloSigTSCount) <- c("subtype","gers")
fig2c  <- ggplot(medulloSigTSCount, aes(x = subtype, y = gers, fill = subtype)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = gers, hjust = 0.5, vjust = 2), color = 'white') +
  xlab("Molecular Subtype") + ylab("Count of GERs") +  
  theme_Publication(base_size = 10) + guides(fill = FALSE)
ggsave(plot = fig2c, filename = "results/plots/Figure2C.png", width = 6, height = 4)
save(fig2c, file = 'results/Fig_2C.RData')

# Save for use in other classes
saveRDS(allGeneCombos, "data/model/bestFeaturesNew.RDS")
saveRDS(medulloGeneSetsUp, "data/model/medulloSetsUp.RDS")

# Figure 2B: T-SNE for DS1
ds1GER <- readRDS("data/RNASeqDataForPlotDS1.RDS")
ds1geneRatioOut <- ds1GER[[1]]
ds1sampAnnot <- ds1GER[[2]]

set.seed(42)
tsneOut <- Rtsne(t(log2(ds1geneRatioOut[intersect(allGeneCombos, rownames(ds1geneRatioOut)),])), initial_dims=200, perplexity=10, max_iter=500)
tsneOut <- data.frame(tsneOut$Y, ds1sampAnnot)
tsneOut$Subgroup <- factor(tsneOut$Subgroup, levels = c("Group3", "Group4", "SHH", "WNT", "Unknown"))
fig3a <- ggplot(tsneOut, aes(X1, X2, shape = Subgroup, color = Subgroup))+
  geom_point(size = 5, alpha = 0.6)+
  theme_bw()+
  # ggtitle("T-SNE Medulloblastoma Gene Ratios - DS1") +
  theme_Publication(base_size = 10) + xlab("PC1") + ylab("PC2")  +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_color_manual(values = c("Group3" = "#F8766D", 
                                "Group4" = "#7CAE00", 
                                "SHH" = "#00BFC4", 
                                "WNT" = "#C77CFF", 
                                "Unknown" = "#000000"))
fig3a
ggsave(plot = fig3a, filename = "results/plots/Figure3A.png", width = 7, height = 6)

# Figure 2C: T-SNE for DS2
ds2GER <- readRDS("data/RNASeqDataForPlotDS2.RDS")
ds2geneRatioOut <- ds2GER[[1]]
ds2sampAnnot <- ds2GER[[2]]

set.seed(150)
tsneOut <- Rtsne(t(log2(ds2geneRatioOut[intersect(allGeneCombos, rownames(ds2geneRatioOut)),])), initial_dims=200, perplexity=10, max_iter=500)
tsneOut <- data.frame(tsneOut$Y, ds2sampAnnot)
colnames(tsneOut)[4] <- "Subgroup"
tsneOut$Subgroup[tsneOut$Subgroup == "U"] <- "Unknown"
tsneOut$Subgroup <- factor(tsneOut$Subgroup, levels = c("Group3", "Group4", "SHH", "WNT", "Unknown"))
fig3b <- ggplot(tsneOut, aes(X1, X2, shape = Subgroup, color = Subgroup))+
  geom_point(size = 5, alpha = 0.6)+
  theme_bw()+
  # ggtitle("T-SNE Medulloblastoma Gene Ratios - DS2") +
  theme_Publication(base_size = 10) + xlab("PC1") + ylab("PC2") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_color_manual(values = c("Group3" = "#F8766D", 
                                "Group4" = "#7CAE00", 
                                "SHH" = "#00BFC4", 
                                "WNT" = "#C77CFF", 
                                "Unknown" = "#000000"))
fig3b
ggsave(plot = q, filename = "results/plots/Figure3B.png", width = 7, height = 6)

# Figure 2C: Frequency Plot
medulloGeneSetsUpTS <- stack(medulloGeneSetsUp)
getSets <- function(x) {
  out <- strsplit(x, split="_")[[1]];
  upG <- out[1]
  downG <- out[2];
  return(c(upG, downG));
}
medulloGeneSetsUpTS <- data.frame(medulloGeneSetsUpTS, t(sapply(medulloGeneSetsUpTS[,1], FUN=getSets)));
colnames(medulloGeneSetsUpTS)[3:4] <- c("upGene", "downGene")
upGeneMat <- data.frame(table(medulloGeneSetsUpTS[,c("upGene", "ind")]));
downGeneMat <- data.frame(table(medulloGeneSetsUpTS[,c("downGene", "ind")]));

getTopXBar <- function(myMat=NULL, topx=5) {
  myTabTmp <- myMat %>% group_by(ind) %>% top_n(topx, Freq)
  myTabTmp <- data.frame(myTabTmp);
  myTabTmp <- myTabTmp[order(-myTabTmp[,3]),];
  myTabTmp[,1] <- factor(myTabTmp[,1], levels=unique(myTabTmp[,1]))
  colnames(myTabTmp)[1] <- "Gene";
  p <- ggplot(myTabTmp, aes(Gene, Freq)) + 
    geom_bar(stat="identity") + facet_grid(~ind, scales="free")
  p <- p + theme_Publication(base_size = 10) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  p <- p + ylab("Frequency")
  return(list(myTabTmp, p))
}

# Top 5
upGeneMat$ind <- factor(upGeneMat$ind, levels = c("Group3", "Group4", "SHH", "WNT"))
fig3c <- getTopXBar(upGeneMat)[[2]]
ggsave(plot = fig3c, filename = "results/plots/Figure3C.png", width = 7, height = 6)

# Bottom 5
downGeneMat$ind <- factor(downGeneMat$ind, levels = c("Group3", "Group4", "SHH", "WNT"))
fig3d <- getTopXBar(downGeneMat)[[2]]
ggsave(plot = s, filename = "results/plots/Figure3D.png", width = 7, height = 6)

# Combine all Figure 3 plots
library(ggpubr)
library(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend <- g_legend(fig3b)
ggexport(ggarrange(ggarrange(fig3a + theme(legend.position="none"), 
                             fig3b + theme(legend.position="none"), nrow = 1, labels = c("A", "B")),
                   mylegend,
                   ggarrange(fig3c, fig3d, nrow = 1, labels = c("C", "D")), nrow = 3, heights = c(5, 1, 6)),
         filename = "results/plots/Figure3.pdf", width = 10, height = 8)

# Figure4: Frequency Plot - Figure 4
ds1geneRatioOut <- ds1GER[[1]]
ds1sampAnnot <- ds1GER[[2]]

IntWNT <- intersect(medulloSetsRSQ$WNT, medulloSetsMIA$WNT)
IntSHH <- intersect(medulloSetsRSQ$SHH, medulloSetsMIA$SHH)
IntGroup3 <- intersect(medulloSetsRSQ$Group3, medulloSetsMIA$Group3)
IntGroup4 <- intersect(medulloSetsRSQ$Group4, medulloSetsMIA$Group4)

getBoxPlot <- function(x, types) {
  tmp <- data.frame(t(ds1geneRatioOut[x,]), ds1sampAnnot[,2]);
  colnames(tmp)[ncol(tmp)] <- "Class"
  tmp[,"Sample"] <- rownames(tmp)
  tmp <- gather(tmp, key="Feature", value="Ratio", -Class, -Sample)
  tmp[,"Dataset"] <- "DS1"
  
  tmp2 <- data.frame(t(ds2geneRatioOut[x,]), ds2sampAnnot[,2]);
  colnames(tmp2)[ncol(tmp2)] <- "Class"
  tmp2[,"Sample"] <- rownames(tmp2)
  tmp2 <- gather(tmp2, key="Feature", value="Ratio", -Class, -Sample)
  tmp2[,"Dataset"] <- "DS2"
  tmp <- rbind(tmp, tmp2);
  
  # Now need to replace features with class represented
  myRep <- paste(types, x, sep=": ")
  for(i in 1:length(myRep)) {
    tmp[,"Feature"] <- gsub(x[i], myRep[i], tmp[,"Feature"])
  }
  ggplot(tmp, aes(x = Class, y = log2(Ratio), fill = Class)) + 
    geom_boxplot(lwd = 0.5, fatten = 0.5, outlier.shape = 1, width = 0.8, outlier.size = 1, aes(fill = Class)) +
    geom_jitter(width = 0.1, pch = 21, stroke = 0.2, aes(fill = Class)) +
    #  theme_bw() +
    theme_Publication2(base_size = 10) +
    theme(axis.text.x = element_blank()) +
    facet_grid(Dataset~Feature) 
}
myRatios <- c("MAK_RGL1", "PTPN5_ROBO1", "ATOH1_OTX2", "AXIN2_DCX")
types <- c("Group3", "Group4", "SHH", "WNT")
getBoxPlot(myRatios, types)
ggsave("results/plots/Figure4.pdf", width = 10, height = 5)

# Classify test datasets
# Source scripts
source("R/04_classifyGSE109401.R")
source("R/04_classifyGSE85217.R")
source("R/04_classifyGSE37418.R")
class109401 <- classifyGSE109401()
class37418 <- classifyGSE37418()
class85217 <- classifyGSE85217()
