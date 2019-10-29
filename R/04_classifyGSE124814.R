# Author: Komal S. Rathi
# Date: 10/19/2019
# Function: Classify medulloblastoma subtypes for 1286 samples across 15 datasets (GSE124814)

training.sets <- c("GSE37418","GSE109401")
sets.without.subtype.info <- c("GSE22569","GSE25219","GSE3526","GSE35974","GSE4036","GSE44971","GSE60862")

# install if not present
if (!require("medulloPackage")) devtools::install_github("d3b-center/medullo-classifier-package")
library(medulloPackage)
library(xlsx)
library(dplyr)
library(data.table)

# download test data (quantile normalized full dataset n = 1641) and metadata
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124814/suppl/GSE124814_HW_expr_matrix.tsv.gz
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124814/suppl/GSE124814_sample_descriptions.xlsx

# format merged expression data
dat <- fread('data/GSE124814_HW_expr_matrix.tsv')
dat <- as.data.frame(dat)
rownames(dat) <- dat$Gene_Symbol
dat$Gene_Symbol <- NULL

# format meta data
meta <- read.xlsx('data/GSE124814_sample_descriptions.xlsx', sheetIndex = 1, startRow = 2)
meta <- meta[,c('Sample.name','title','characteristics..subgroup.supplied.original','characteristics..subgroup.supplied.renamed','characteristics..subgroup.relabeled')]
meta$title <- gsub('-.*','',meta$title)
colnames(meta)[3:5] <- c("original","renamed","relabeled")
meta <- meta %>%
  dplyr::mutate_all(as.character) %>%
  as.data.frame()

######################################## QC ######################################## 
# Format metadata to include only 5 types: WNT, SHH, Group3, Group4 and U (for Unknown)
# original subtypes
meta$original <- ifelse(meta$original %in% c("A", "WNT", "MB_WNT", "WNT_beta", "WNT_alpha"), "WNT",
                        ifelse(meta$original %in% c("B", "SHH", "MB_SHH", "SHH_beta", "SHH_alpha", "SHH_gamma", "SHH_delta"), "SHH",
                               ifelse(meta$original %in% c("E", "G3", "Group 3", "MB_GRP3", "Group3_alpha", "Group3_beta", "Group3_gamma"), "Group3",
                                      ifelse(meta$original %in% c("C", "D", "G4", "Group 4", "MB_GRP4", "Group4_gamma", "Group4_beta", "Group4_alpha"), "Group4",
                                             ifelse(meta$original %in% c("Unknown","NA", "SHH OUTLIER"), "U", meta$original)))))

# renamed subtypes
meta$renamed <- ifelse(meta$renamed == "G3", "Group3", 
                       ifelse(meta$renamed == "G4", "Group4",
                              ifelse(meta$renamed %in% c("Unknown","NA"), "U", meta$renamed)))

# original and renamed subtypes should be identical
identical(meta$original, meta$renamed) # TRUE

# relabeled subtypes
meta$relabeled <- ifelse(meta$relabeled == "G3", "Group3", 
                         ifelse(meta$relabeled == "G4", "Group4",
                                ifelse(meta$relabeled %in% c("Unknown","NA"), "U", meta$relabeled)))
######################################## QC ######################################## 

# remove training sets from meta and dat
rownames(meta) <-  meta$Sample.name
meta <- meta[-which(meta$title %in% c(training.sets, sets.without.subtype.info)),] # 1286/1641
dat <- dat[,rownames(meta)]
actual <- data.frame(meta[,c('Sample.name','title','relabeled')])
colnames(actual) <- c("id","study","subtype")

# save formatted and clean data for classification and summary plots
if(identical(colnames(dat), rownames(actual))){
  print("rownames and colnames match!")
  save(dat, file = 'data/expr.RData')
  save(actual, file = 'data/meta.RData')
} else {
  print("match rownames and colnames")
  return(NULL)
}

######################################## QC ######################################## 
# how different are renamed and relabeled subtypes
dim(meta[which(meta$renamed != meta$relabeled),])  # 178
######################################## QC ######################################## 

# classify subtypes
pred124814 <- classify(exprs = dat) 
meta$predicted <- pred124814$best.fit # add to metadata

######################################## QC ######################################## 
# how different are predicted vs renamed and relabled subtypes
dim(meta[which(meta$predicted != meta$relabeled),]) # 96
dim(meta[which(meta$predicted != meta$renamed),]) # 188
######################################## QC ######################################## 

# calculate summary statistics
# stats
acc1 <- calcStats(myClassActual = meta$original, pred124814$best.fit) # original subtypes 96%
acc2 <- calcStats(myClassActual = meta$renamed, pred124814$best.fit) # renamed subtype 96%
acc3 <- calcStats(myClassActual = meta$relabeled, pred124814$best.fit) # relabeled subtypes 98%
confusion.matrix <- acc3[[1]]
overall.stats <- acc3[[2]]
class.stats <- acc3[[3]]
save(confusion.matrix, overall.stats, class.stats, file = 'results/GSE124814_results.RData')
