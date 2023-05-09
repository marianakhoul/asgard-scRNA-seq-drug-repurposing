

library('Seurat')

#Load normal sample Ind5 from GSE113197 dataset
celltype<-read.table(file="https://raw.githubusercontent.com/lanagarmire/Single-cell-drug-repositioning/master/Normal_celltype.txt",header = T,check.names=FALSE)
data<-read.table(file="GSM3099847_Ind5_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype2<-subset(celltype,sample=="Ind5" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype2))
data<-data[,common]
Epithelial2 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype2,cell=colnames(data),type="Normal"))

#Load normal sample Ind6 from GSE113197 dataset
data<-read.table(file="GSM3099848_Ind6_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype3<-subset(celltype,sample=="Ind6" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype3))
data<-data[,common]
Epithelial3 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype3,cell=colnames(data),type="Normal"))

#Load normal sample Ind7 from GSE113197 dataset
data<-read.table(file="GSM3099849_Ind7_Expression_Matrix.txt",header = T,check.names=FALSE)
row.names(data)<-data[,1]
data<-data[,-1]
celltype4<-subset(celltype,sample=="Ind7" & celltype %in% c("Luminal_L2_epithelial_cells","Luminal_L1.1_epithelial_cells", "Luminal_L1.2_epithelial_cells", "Basal_epithelial_cells"))
common <- intersect(colnames(data), rownames(celltype4))
data<-data[,common]
Epithelial4 <- CreateSeuratObject(counts = data, project = "Epithelial", min.cells = 3, min.features = 200,meta.data=data.frame(celltype4,cell=colnames(data),type="Normal"))

#Load cancer sample PDX110 from GSE123926 dataset
TNBC_PDX.data<- Read10X(data.dir = "GSM3516947_PDX110")
TNBC.PDX2 <- CreateSeuratObject(counts = TNBC_PDX.data, project = "TNBC", min.cells = 3, min.features = 200, meta.data=data.frame(row.names=colnames(TNBC_PDX.data), cell=colnames(TNBC_PDX.data), sample="PDX-110",type="TNBC.PDX"))

#Load cancer sample PDX322 from GSE123926 dataset
TNBC_PDX.data<- Read10X(data.dir = "GSM3516948_PDX322")
TNBC.PDX3 <- CreateSeuratObject(counts = TNBC_PDX.data, project = "TNBC", min.cells = 3, min.features = 200, meta.data=data.frame(row.names=colnames(TNBC_PDX.data), cell=colnames(TNBC_PDX.data), sample="PDX-332",type="TNBC.PDX"))

