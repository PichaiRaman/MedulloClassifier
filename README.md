# MedulloClassifier

## Objective

The goal of this project is to develop a model/function that can accurately predict amongst 4 molecular subtypes of Medulloblastoma, Sonic Hedgehog (SHH), WNT, Group 3, and Group 4. These subtypes were first identified using Non-negative matrix factorization on microarray data  (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4874239/). Since then, these subtypes are widely used in both research and clinical practice (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4105823/). It is important to note there are other studies that classify medulloblastoma's into many more subgroups. This classifier is not capable of that currently as we don't have adequete training data (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6163053/). 

## Method

The classifier was built using mulitiple microarray and RNA-Seq datasets. As there is tremendous variability in the dynamic range and values of microarray (between different platforms) and RNA-Seq, gene ratio's instead of gene expression were used as features for the model. Discrimanating features (gene ratios) were identified using the limma package (https://academic.oup.com/nar/article/43/7/e47/2414268) and samples are classified using an un-weighted sum of normalized scores.  The data was trained using 2 datasets and tested in 2 separate datasets. The classifier is still being refined but can currently differentiate between the 4 subtypes at a greater than 95% accuracy. 



## Results





