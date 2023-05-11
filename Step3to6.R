#!/usr/bin/env Rscript

# Load libraries
library('Seurat')
library('Asgard')
library('Hmisc')

load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/ToStep2.RData")

## Step 3: Single-cell comparison
#Case sample names
Case=c("PDX-110","PDX-332")

#Control sample names
Control=c("Normal1","Normal2","Normal3")

#Get differential genes from Seurat (Wilcoxon Rank Sum test); we will only need this for the samples we have since 1 tumor and 1 normal
DefaultAssay(SC.integrated) <- "RNA"
set.seed(123456)
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.integrated@meta.data$celltype)){
  Idents(SC.integrated) <- "celltype"
  c_cells <- subset(SC.integrated, celltype == i)
  Idents(c_cells) <- "type"
  C_data <- FindMarkers(c_cells, ident.1 = "TNBC.PDX", ident.2 = "Normal")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val)
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names

##Step 4: Mono-drug repurposing for every cell type

#Load tissue specific drug reference produced by PrepareReference function as mentioned above. Please select proper tissue accroding to the disease.
my_gene_info<-read.table(file="/home/sas1782/asgard-scRNA-seq-drug-repurposing/DrugReference/breast_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="/home/sas1782/asgard-scRNA-seq-drug-repurposing/DrugReference/breast_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = '/home/sas1782/asgard-scRNA-seq-drug-repurposing/DrugReference/breast_rankMatrix.txt',
                               probe.to.genes = my_gene_info, 
                               drug.info = my_drug_info)

#Repurpose mono-drugs for every cell type                               
Drug.ident.res = GetDrug(gene.data = Gene.list, 
                        drug.ref.profiles = drug.ref.profiles, 
                        repurposing.unit = "drug", 
                        connectivity = "negative", 
                        drug.type = "FDA")
                    
## Step 5: Estimation of Drug Score

GSE92742.gctx.path="/home/sas1782/asgard-scRNA-seq-drug-repurposing/GEOProfiles/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
GSE70138.gctx.path="/home/sas1782/asgard-scRNA-seq-drug-repurposing/GEOProfiles/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
Tissue="breast"
Drug.score<-DrugScore(SC.integrated=SC.integrated,
                     Gene.data=Gene.list,
                     Cell.type=NULL, 
                     Drug.data=Drug.ident.res,
                     FDA.drug.only=TRUE,
                     Case=Case, 
                     Tissue="breast",
                     GSE92742.gctx=GSE92742.gctx.path,
                     GSE70138.gctx=GSE70138.gctx.path)
#Cell.type: select cell types/clusters to be used for drug score estimation
#Case: select samples to be used for drug score estimation

## Step 6: Select mono-drug therapies
#Select drug using drug socre

Final.drugs<-subset(Drug.score,Drug.therapeutic.score>quantile(Drug.score$Drug.therapeutic.score, 0.99,na.rm=T) & FDR <0.05)


#Select drug for individual clusters
Final.drugs<-TopDrug(SC.integrated=SC.integrated,
                   Drug.data=Drug.ident.res,
                   Drug.FDR=0.1,
                   FDA.drug.only=TRUE,
                   Case=Case.samples,
                   DrugScore=FALSE
)
save(Final.drugs,Drug.score,Drug.ident.res,my_drug_info,my_gene_info,Gene.list,drug.ref.profiles, file = "/home/sas1782/asgard-scRNA-seq-drug-repurposing/Step6Final.RData")

