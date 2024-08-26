### EWAS analysis
# INPUT:
# Methylation M-value matrix (rows=CpGs, columns=sampleID)
# Metadata file including the sampleID, covariables and
# estimated blood cell types (cd4t, cd8t, nk, mono, neu, bcell)

#load packages
library(readxl)
library(tableone)
library(dplyr)
library(matrixStats)
library(limma)
library(data.table)

#create functions
#Extreme outlier exclusion according to Tukey
IQR.removal <- function(meth)
  {
  rowIQR <- rowIQRs(meth, na.rm = T)
  row2575 <- rowQuantiles(meth, probs = c(0.25, 0.75), na.rm = T)
  maskL <- meth < row2575[,1] - 3 * rowIQR 
  maskU <- meth > row2575[,2] + 3 * rowIQR 
  meth[maskL] <- NA
  meth[maskU] <- NA
  meth
}

#Ewas function
ewas.function <-  function(meth, pheno, variable.of.interest)
  {
  model.covariates <- colnames(pheno)[-which(colnames(pheno) %in% c(variable.of.interest,"sample.id"))]
  des = model.matrix(reformulate(paste0("pheno$",c(variable.of.interest,model.covariates))))
  fit = lmFit(meth, des, method="robust")
  fit.ebayes = eBayes(fit)
  n = rowSums(!is.na(meth))
  se = (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$stdev.unscaled))])
  res = data.frame(n=n,
                   coef=fit.ebayes$coefficient[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$coefficient))],
                   se=se,
                   p=fit.ebayes$p.value[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$p.value))])
  res
}

#data path
setwd() #set path for stored data

#load methylation Mvalues
load("Mval.RData") #name according to your methylation matrix

# load and prepare metadata
metadata <- read_excel("Pheno.xlsx") %>%
  mutate_at(vars(age, PC1, PC2, PC3, cd8t:neu), as.numeric) %>%
  mutate(sex = ifelse(sex == 0, "Male", "Female"))
metadata$maternal_t1d<-as.factor(metadata$maternal_t1d)

#check matching and ordering of methylation matrix and metadata
Mval<-Mval[,as.character(colnames(Mval)) %in% as.character(metadata$sample_ID)]
Mval<-Mval[,order(match(as.character(colnames(Mval)), as.character(metadata$sample_ID)))]
Mval<-as.matrix(Mval)

#run EWAS
cell.names <- c("cd8t", "cd4t", "nk", "bcell", "mono", "neu")
covs <- c("maternal_t1d","age", "sex", "PC1", "PC2", "PC3", cell.names)

Mval <- IQR.removal(Mval)
res <- ewas.function(
    meth=Mval, 
    pheno=metadata[,colnames(metadata) %in% covs], 
    variable.of.interest = "maternal_t1d")

#add P FDR according to Benjamini-Hochberg
res$P.fdr<-p.adjust(res$p, method = "BH")
res<-res[order(res$P.fdr),]