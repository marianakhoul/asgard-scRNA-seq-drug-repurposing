#!/usr/bin/env Rscript

# Load libraries
library('Seurat')
library('Asgard')
library('Hmisc')

## Step 1: Load single-cell RNA-seq data
load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/data/FinalQC_ch1newT_nodoublets.RData")
load("/home/sas1782/asgard-scRNA-seq-drug-repurposing/data/FinalQC_ch1n_nodoublets.RData")
ch1n_nodoublets_updated <- UpdateSeuratObject(ch1n_nodoublets)
ch1t_nodoublets_updated <- UpdateSeuratObject(ch1newT_nodoublets)

## Step 2: Single-cell alignment
#list of single cell objects
SC.list<-list(ch1t_nodoublets_updated=ch1t_nodoublets_updated,ch1n_nodoublets_updated=ch1n_nodoublets_updated)
CellCycle=FALSE #Set it TRUE if you want to do Cell Cycle Regression
anchor.features=2000

for (i in 1:length(SC.list)) {
     SC.list[[i]] <- NormalizeData(SC.list[[i]], verbose = FALSE)
     SC.list[[i]] <- FindVariableFeatures(SC.list[[i]], selection.method = "vst",
                           nfeatures = anchor.features, verbose = FALSE)
    }
    SC.anchors <- FindIntegrationAnchors(object.list = SC.list,anchor.features = anchor.features, dims = 1:15)
    SC.integrated <- IntegrateData(anchorset = SC.anchors, dims = 1:15)
    DefaultAssay(SC.integrated) <- "integrated"
    if(CellCycle){
    ##Cell Cycle Regression
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    SC.integrated <- CellCycleScoring(SC.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    SC.integrated <- ScaleData(SC.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SC.integrated))
    SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
    }else{
     ##Run the standard workflow for visualization and clustering
     SC.integrated <- ScaleData(SC.integrated, verbose = FALSE)
     SC.integrated <- RunPCA(SC.integrated, npcs = 15, verbose = FALSE)
    }
    ##t-SNE and Clustering
    SC.integrated <- RunUMAP(SC.integrated, reduction = "pca", dims = 1:15)
    SC.integrated <- FindNeighbors(SC.integrated, reduction = "pca", dims = 1:15)
    SC.integrated <- FindClusters(SC.integrated, algorithm = 1, resolution = 0.4)

    ##Cell Type Annotation, set by.CellType=TRUE if you want to annotate cell  type.
    by.CellType=FALSE
    if(by.CellType == TRUE){
     data <- as.matrix(SC.integrated@assays$RNA@data)
     hpca.se <- HumanPrimaryCellAtlasData()
     pred.hpca <- SingleR(test = data, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
     cell.label <- data.frame(row.names = row.names(pred.hpca),celltype=pred.hpca$labels)
     if(length(SC.integrated@meta.data$celltype)>0){
      SC.integrated@meta.data$celltype <- cell.label$celltype
     }else{
       SC.integrated@meta.data <- cbind(SC.integrated@meta.data,cell.label)
     }
     new.cells <- data.frame()
     for(i in unique(SC.integrated$seurat_clusters)){
      sub.data <- subset(SC.integrated,seurat_clusters==i)
      temp <- table(sub.data@meta.data$celltype)
      best.cell <- names(which(temp==temp[which.max(temp)]))
      cells.temp <- data.frame(cell.id=row.names(sub.data@meta.data),celltype=best.cell)
      new.cells <- rbind(new.cells,cells.temp)
     }
     cell.meta <- SC.integrated@meta.data
     cell.id <- rownames(cell.meta)
     row.names(new.cells) <- new.cells[,1]
     new.cells <- new.cells[cell.id,]
     SC.integrated@meta.data$celltype <- new.cells$celltype
    }else{
     SC.integrated@meta.data$celltype <- paste0("C",as.numeric(SC.integrated@meta.data$seurat_clusters))
    }
    
#Change sample names
sample<-SC.integrated@meta.data$sample
sample[which(sample=="CH1_NewT")]<-"Tumor"
sample[which(sample=="CH1_N")]<-"Normal"
SC.integrated@meta.data$sample<-sample

pdf(file = "/home/sas1782/asgard-scRNA-seq-drug-repurposing/umapch1_data.pdf")
DimPlot(SC.integrated, reduction = "umap", split.by = "sample",group.by = "celltype")
dev.off()

save(anchor.features, celltype,SC.integrated,SC.anchors,sample,i,SC.list,data,common,ch1n_nodoublets_updated,ch1t_nodoublets_updated, file = "/home/sas1782/asgard-scRNA-seq-drug-repurposing/ToStep2_ch1_data.RData")

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

save(Drug.ident.res,my_drug_info,my_gene_info,Gene.list,drug.ref.profiles, file = "/home/sas1782/asgard-scRNA-seq-drug-repurposing/Step4_ch1_data.RData")

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

save(Final.drugs,Drug.score,Drug.ident.res,my_drug_info,my_gene_info,Gene.list,drug.ref.profiles, file = "/home/sas1782/asgard-scRNA-seq-drug-repurposing/Step5_ch1_data.RData")

#Select drug for individual clusters
Final.drugs<-TopDrug(SC.integrated=SC.integrated,
                   Drug.data=Drug.ident.res,
                   Drug.FDR=0.1,
                   FDA.drug.only=TRUE,
                   Case=Case
)
save(Final.drugs, file = "/home/sas1782/asgard-scRNA-seq-drug-repurposing/Step6_ch1_data_Final.RData")
