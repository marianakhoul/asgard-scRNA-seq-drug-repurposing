#!/usr/bin/env Rscript

# Load libraries
library('Seurat')
library('Asgard')
library('Hmisc')

## Step 1: Load single-cell RNA-seq data
load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/data/FinalQC_ch1newT_nodoublets.RData")
load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/data/FinalQC_ch1n_nodoublets.RData")

## Step 2: Single-cell alignment
SC.list<-list(TNBC.PDX2=TNBC.PDX2,TNBC.PDX3=TNBC.PDX3,Epithelial2=Epithelial2,Epithelial3=Epithelial3,Epithelial4=Epithelial4)
CellCycle=FALSE #Set it TRUE if you want to do Cell Cycle Regression
anchor.features=2000
