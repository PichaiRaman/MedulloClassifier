##################################################
# Purpose: Code to create manuscript figures
# Author: Pichai Raman, Komal S. Rathi
# Date: 11/07/2019
##################################################

ds1 <- readRDS('data/RNASeqDataForPlotDS1.RDS')
ds2 <- readRDS('data/RNASeqDataForPlotDS2.RDS')
source('R/utils/pubTheme.R')

library(ggplot2)
library(Rtsne)

# Figure 2A
expDataFts_QN <-  ds1[[1]]
sampAnnot <- ds1[[2]]
set.seed(42)
tsneOut <- Rtsne::Rtsne(t(expDataFts_QN), initial_dims=100, perplexity=20, max_iter=1000)
tsneOut <- data.frame(tsneOut$Y, sampAnnot)
tsneOut$Subgroup <- factor(tsneOut$Subgroup, levels = c("Group3", "Group4", "SHH", "WNT", "U"))
p <- ggplot2::ggplot(tsneOut, aes(X1, X2, shape = Subgroup, color = Subgroup))+
  ggplot2::geom_point(size = 5, alpha = 0.6)+
  ggplot2::theme_bw()+
  ggplot2::ggtitle("T-SNE Medulloblastoma Gene Ratios - DS1")+
  theme_Publication() +
  scale_color_manual(values = c("Group3" = "#F8766D", "Group4" = "#7CAE00", "SHH" = "#00BFC4", "WNT" = "#C77CFF", "U" = "#000000"))
p
ggsave(plot = p, filename = "results/plots/Figure2A.png", width = 7, height = 6)

# Figure 2B
expDataFts_QN <-  ds2[[1]]
sampAnnot <- ds2[[2]]
set.seed(150)
tsneOut <- Rtsne::Rtsne(t(expDataFts_QN), initial_dims=100, perplexity=20, max_iter=1000)
tsneOut <- data.frame(tsneOut$Y, sampAnnot)
tsneOut$Subgroup <- factor(tsneOut$Subgroup, levels = c("Group3", "Group4", "SHH", "WNT", "U"))
p <- ggplot2::ggplot(tsneOut, aes(X1, X2, shape = Subgroup, color = Subgroup))+
  ggplot2::geom_point(size = 5, alpha = 0.6)+
  ggplot2::theme_bw()+
  ggplot2::ggtitle("T-SNE Medulloblastoma Gene Ratios - DS2")+
  theme_Publication() +
  scale_color_manual(values = c("Group3" = "#F8766D", "Group4" = "#7CAE00", "SHH" = "#00BFC4", "WNT" = "#C77CFF", "U" = "#000000"))

p
ggsave(plot = p, filename = "results/plots/Figure2B.png", width = 7, height = 6)

#
