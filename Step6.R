#!/usr/bin/env Rscript

# Load libraries
library('Seurat')
library('Asgard')
library('Hmisc')

#load data
load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/ToStep2.RData")
load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/Step5.RData")

Case=c("PDX-110","PDX-332")
Control=c("Normal1","Normal2","Normal3")

#Select drug for individual clusters
Final.drugs<-TopDrug(SC.integrated=SC.integrated,
                   Drug.data=Drug.ident.res,
                   Drug.FDR=0.1,
                   FDA.drug.only=TRUE,
                   Case=Case
)
save(Final.drugs, file = "/home/sas1782/asgard-scRNA-seq-drug-repurposing/Step6Final.RData")
