# Asgard scRNA-seq Drug Repurposing Pipeline 
This pipeline is built to run [Asgard](https://github.com/lanagarmire/Asgard) on SLURM job scheduling high performance cluster.
It follows the steps highlighted in their github README file.

## Steps:
1. **Prepare Drug Referecne Library**
* Generate tissue specific drug references from GSE70138 and GSE92742
2. **Drug Repurposing**
* Load single-cell RNA-seq data
* Single-cell alignment
* Single-cell comparison
* Mono-drug repurposing for every cell type
* Estimation of drug score
* Select mono-drug therapies
* (Optional) Drug combination analysis

Download Example Data:
```
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123926/suppl/GSE123926_RAW.tar'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113197/suppl/GSE113197_RAW.tar'

mkdir GSE123926
mkdir GSE113197

tar -xf GSE123926_RAW.tar -C ./GSE123926
tar -xf GSE113197_RAW.tar -C ./GSE113197

cd ./GSE113197
gunzip GSM3099847_Ind5_Expression_Matrix.txt.gz
gunzip GSM3099848_Ind6_Expression_Matrix.txt.gz
gunzip GSM3099849_Ind7_Expression_Matrix.txt.gz
```
