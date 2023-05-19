#!/usr/bin/env Rscript

# Load libraries
library('Seurat')
library('Asgard')
library('Hmisc')

## Step 1: Load single-cell RNA-seq data
load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/data/CH103_T_Merged_Nodoublets_Integrated.RData")
load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/data/FinalQC_ch103n_nodoublets.RData")
ch103n_nodoublets_updated <- UpdateSeuratObject(ch103n_nodoublets)
ch103t_nodoublets_updated <- UpdateSeuratObject(ch103t_merged_nodoublets)

filtered_ch103t_nodoublets_updated <- subset(x = ch103t_nodoublets_updated, 
                         subset= (nUMI >= 500) & 
                           (nGene > 250) & (nGene<2500) &
                           (log10GenesPerUMI > 0.80) & 
                           (percent.mt < 20))

filtered_ch103n_nodoublets_updated <- subset(x = ch103n_nodoublets_updated, 
                         subset= (nUMI >= 500) & 
                           (nGene > 250) & (nGene<2500) &
                           (log10GenesPerUMI > 0.80) & 
                           (percent.mt < 20))

## Step 2: Single-cell alignment
#list of single cell objects
SC.list<-list(ch103t_nodoublets_updated=filtered_ch103t_nodoublets_updated,ch103n_nodoublets_updated=filtered_ch103n_nodoublets_updated)
anchor.features=2000

for (i in 1:length(SC.list)) {
     SC.list[[i]] <- NormalizeData(SC.list[[i]], verbose = FALSE)
     SC.list[[i]] <- FindVariableFeatures(SC.list[[i]], selection.method = "vst",
                           nfeatures = anchor.features, verbose = FALSE)
    }
SC.anchors <- FindIntegrationAnchors(object.list = SC.list,anchor.features = anchor.features, dims = 1:15)
SC.integrated <- IntegrateData(anchorset = SC.anchors, dims = 1:15)
DefaultAssay(SC.integrated) <- "integrated"
SC.integrated <- ScaleData(SC.integrated, verbose = FALSE)
SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
SC.integrated <- RunUMAP(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindNeighbors(SC.integrated, reduction = "pca", dims = 1:15)
SC.integrated <- FindClusters(SC.integrated, algorithm = 1, resolution = 0.4)

#new cell types
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Kidney" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = SC.integrated[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(SC.integrated@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(SC.integrated@meta.data[SC.integrated@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(SC.integrated@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

SC.integrated@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  SC.integrated@meta.data$customclassif[SC.integrated@meta.data$seurat_clusters == j] = as.character(paste0(cl_type$cluster[1],"_",cl_type$type[1]))
}

SC.integrated@meta.data$celltype <- SC.integrated@meta.data$customclassif

sample<-SC.integrated@meta.data$sample
sample[which(sample=="CH1_NewT")]<-"Tumor"
sample[which(sample=="CH1_N")]<-"Normal"
SC.integrated@meta.data$sample<-sample

## Step 3: Single-cell comparison
#Case sample names
Case=c("Tumor")
#Control sample names
Control=c("Normal")

DefaultAssay(SC.integrated) <- "RNA"
set.seed(123456)
Gene.list <- list()
C_names <- NULL
for(i in unique(SC.integrated@meta.data$celltype)){
  Idents(SC.integrated) <- "celltype"
  c_cells <- subset(SC.integrated, celltype == i)
  Idents(c_cells) <- "type"
  C_data <- FindMarkers(c_cells, ident.1 = "Tumor", ident.2 = "Normal")
  C_data_for_drug <- data.frame(row.names=row.names(C_data),score=C_data$avg_log2FC,adj.P.Val=C_data$p_val_adj,P.Value=C_data$p_val)
  Gene.list[[i]] <- C_data_for_drug
  C_names <- c(C_names,i)
}
names(Gene.list) <- C_names

##Step 4: Mono-drug repurposing for every cell type

#Load tissue specific drug reference produced by PrepareReference function as mentioned above. Please select proper tissue accroding to the disease.
my_gene_info<-read.table(file="/home/sas1782/asgard-scRNA-seq-drug-repurposing/DrugReference/kidney_gene_info.txt",sep="\t",header = T,quote = "")
my_drug_info<-read.table(file="/home/sas1782/asgard-scRNA-seq-drug-repurposing/DrugReference/kidney_drug_info.txt",sep="\t",header = T,quote = "")
drug.ref.profiles = GetDrugRef(drug.response.path = '/home/sas1782/asgard-scRNA-seq-drug-repurposing/DrugReference/kidney_rankMatrix.txt',
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
Tissue="kidney"
Drug.score<-DrugScore(SC.integrated=SC.integrated,
                     Gene.data=Gene.list,
                     Cell.type=NULL, 
                     Drug.data=Drug.ident.res,
                     FDA.drug.only=TRUE,
                     Case=Case, 
                     Tissue="kidney",
                     GSE92742.gctx=GSE92742.gctx.path,
                     GSE70138.gctx=GSE70138.gctx.path)
#Cell.type: select cell types/clusters to be used for drug score estimation
#Case: select samples to be used for drug score estimation
#Please use " " instead of "-" in tissue name, for example, while haematopoietic-and-lymphoid-tissue is the prefix of the drug reference files, the corresponding tissue name is "haematopoietic and lymphoid tissue". 


## Step 6: Select mono-drug therapies
#Select drug using drug socre

Final.drugs<-subset(Drug.score,Drug.therapeutic.score>quantile(Drug.score$Drug.therapeutic.score, 0.99,na.rm=T) & FDR <0.05)

#Select drug for individual clusters
Final.drugs<-TopDrug(SC.integrated=SC.integrated,
                   Drug.data=Drug.ident.res,
                   Drug.FDR=0.1,
                   FDA.drug.only=TRUE,
                   Case=Case
)


