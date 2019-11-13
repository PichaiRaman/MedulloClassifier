# Author: Komal S. Rathi
# Date: 10/19/2019
# Function: Classify medulloblastoma subtypes (GSE124814)
# Split by dataset (22 datasets minus 7 which have only Us)
# 1286 samples across 15 datasets
# Create summary plots

library(medulloPackage)
library(plyr)
library(metafor)
library(meta)
library(ggplot2)
library(reshape2)
source('R/utils/pubTheme.R')

# function to classify all 15 datasets and save the output
if(file.exists('results/GSE124814_split_results.RData')){
  print("file exists..")
  load('results/GSE124814_split_results.RData')
}  else {
  # load data (generated using 05_analysis.R)
  load('data/expr.RData')
  load('data/meta.RData')
  study.ct <- plyr::count(actual$study)

  classify.and.summarize <- function(x, dat){
    print(paste0('Study: ', unique(x$study)))
    rownames(x) <- x$id
    dat.sub <- dat[,rownames(x)]
    actual <- x$subtype
    pred <- classify(exprs = dat.sub)
    res <- calcStats(myClassActual = actual,  pred$best.fit)
    accuracy <- res[[2]]$stats[[1]]
    sens_spec <- res[[3]][1:2]
    rownames(sens_spec) <- gsub("Class: ","", rownames(sens_spec))
    sens_spec <- melt(as.matrix(sens_spec))
    df <- data.frame(var = paste0(sens_spec$Var1,"_", sens_spec$Var2), val = sens_spec$value)
    df <- rbind(df, data.frame(var = "Accuracy", val = accuracy))
    return(df)
  }
  
  big.res <- ddply(actual, .variables = "study", .fun = function(x) classify.and.summarize(x, dat))
  big.res <- dcast(big.res, study~var, value.var = 'val')
  save(big.res, study.ct, file = 'results/GSE124814_split_results.RData')
}

# remove studies that have all Us (15 studies left)
big.res <- big.res[which(big.res$Accuracy != "NaN%"),]
accuracy <- big.res[,c('study','Accuracy')]
cor <- as.numeric(gsub('%','',accuracy$Accuracy))
studies <- accuracy$study
n <- study.ct[study.ct$x %in% studies,'freq']

# format to add labels
for.meta <- data.frame(studies, n, accuracy = cor)
rownames(for.meta) <- for.meta$studies
for.meta$label = paste0(for.meta$studies, ' (n = ',for.meta$n,')')
for.meta <- for.meta[order(for.meta$accuracy, for.meta$n),]

# Accuracy: waterfall plot
for.meta$label = paste0(for.meta$studies, '\n(n = ',for.meta$n,')')
for.meta$label <- factor(for.meta$label, levels = for.meta$label)
p <- ggplot(for.meta, aes(x = label, y = accuracy, fill = accuracy)) + 
  geom_bar(stat = "identity")  + 
  geom_text(aes(label = paste0(accuracy,"%"), vjust = 2), size = 2.5, color = "white") +
  theme_Publication2(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle(paste0("Accuracy across 15 datasets (median = ",median(sort(for.meta$accuracy)),"%)")) +
  ylab("Accuracy (%)") + xlab("") +
  geom_abline(slope = 0, intercept = median(sort(for.meta$accuracy)),  col = "red", lty=2) + 
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") + guides(fill = FALSE)
p
ggsave(p, filename = 'results/plots/Figure4A.png', device = "png", width = 7, height = 3)

# Sensitivity and Specificity: Line plot
sens_spec <- melt(big.res, id.vars = "study")
sens_spec <- sens_spec[which(sens_spec$variable != "Accuracy"),]
sens_spec <- cbind(sens_spec, colsplit(sens_spec$variable, "_", names = c("var","type")))
sens_spec$value <- as.numeric(gsub("%","",sens_spec$value))
sens_spec <- merge(sens_spec, study.ct, by.x = 'study', by.y = 'x')
sens_spec$study <- paste0(sens_spec$study,"\n(n = ",sens_spec$freq,")")
sens_spec$type <- factor(sens_spec$type, levels = c("Specificity","Sensitivity"))
q <- ggplot(sens_spec[!is.na(sens_spec$value),], aes(x = study, value, group = type)) + 
  geom_line(aes(color = type), size = 0.5, alpha = 0.8, position=position_dodge(width=0.2)) +
  geom_point(aes(color = type, shape = type), size = 2, alpha = 0.6, position=position_dodge(width=0.2)) + 
  ylab("Metric (in %)") +
  facet_wrap(~var, nrow = 4) + theme_Publication(base_size = 8) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(shape = 'Type', fill = 'Type', color = 'Type') +
  scale_color_manual(values=c("darkorchid3", "darkorange3")) +
  theme(legend.text=element_text(size = 8), legend.title = element_text(size = 10)) +
  ylim(c(1, 100))
q
ggsave(q, filename = 'results/plots/Figure4B.png', device = "png", width = 7, height = 5)
