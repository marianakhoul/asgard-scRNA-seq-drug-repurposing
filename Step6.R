#!/usr/bin/env Rscript

# Load libraries
library('Seurat')
library('Asgard')
library('Hmisc')

load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/Step5.RData")
#Select drug for individual clusters
Final.drugs<-TopDrug(SC.integrated=SC.integrated,
                   Drug.data=Drug.ident.res,
                   Drug.FDR=0.1,
                   FDA.drug.only=TRUE,
                   Case=Case
)
save(Final.drugs, file = "/home/sas1782/asgard-scRNA-seq-drug-repurposing/Step6Final.RData")
