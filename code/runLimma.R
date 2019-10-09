##############################################
#Author: Pichai Raman
#Purpose: Run limma
#Date: 9/28/2017
##############################################

#Call libraries
library(limma);
library(scales);

#Code to run limma
runLimma <- function(targ, cont, myData)
{    
    targDF <- targ;
    targ <- targ[,1];
  
    #Main code
    design <- model.matrix(~0+targ);
    colnames(design) <- gsub("targ", "", colnames(design));
    dataExpTmp <- myData[, rownames(targDF)];
    fit <- lmFit(dataExpTmp, design);
    
    contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list
    (cont),levels=list(design))))
    
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    outputAll <- topTable(fit2, number=60000)
        
    #Sort according to p-value
    outputAll <- outputAll[order(outputAll[,"P.Value"]),];
    outputAll <- unique(outputAll);
    return(list(fit2, outputAll));       
}
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

#volcano plot, takes in limma analysis
plotVolcano <- function(result, title="Volcano Plot", colorCol=NULL, otherCol="blue", lfcCol=NULL, pvalCol=NULL)
{
p <- ggplot(result, aes_string(x=lfcCol, y=pvalCol, color=colorCol))+geom_point();
p <- p+scale_y_continuous(trans=reverselog_trans(10))+ggtitle(title)+scale_colour_manual(values =c("gray", otherCol));
p <- p+theme_Publication()+labs(color="DEG", size="")+guides(size=F);
p <- p+xlab("Log FC")+ylab("-Log10 P-Value")
return(p);
}

#volcano plot text, takes in limma analysis
plotVolcanoText <- function(result, title="Volcano Plot", colorCol=NULL, otherCol="blue", lfcCol=NULL, pvalCol=NULL, myLabel=NULL)
{
p <- ggplot(result, aes_string(x=lfcCol, y=pvalCol, label=myLabel,color=colorCol))+geom_text();
p <- p+scale_y_continuous(trans=reverselog_trans(10))+ggtitle(title)+scale_colour_manual(values =c("gray", otherCol));
p <- p+theme_Publication()+labs(color="DEG", size="")+guides(size=F);
p <- p+xlab("Log FC")+ylab("-Log10 P-Value")
return(p);

}
