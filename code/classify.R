################################
#Purpose: Classify Medullo Samples
#Author: Sherjeel Arif
#Date: 8/5/2019
################################


classify <- function(exprs = NULL, signatureProbesLoc="../../results/V5_CombosUp/bestFeaturesNew.RDS",medulloGeneSetsUpLoc="../../results/V5_CombosUp/medulloSetsUp.RDS")
{
  
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
  
  myMat <- calcScoreMat(geneRatioOut, medulloGeneSetsUp);
  myClassPred <- colnames(myMat)[max.col(myMat,ties.method="first")]
  return(myClassPred)
  
}


