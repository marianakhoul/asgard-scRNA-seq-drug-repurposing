# Asgard scRNA-seq Drug Repurposing Pipeline 
This pipeline is built to run [Asgard](https://github.com/lanagarmire/Asgard) on SLURM job scheduling high performance cluster.
It follows the steps highlighted in their github README file.

Download Example Data:
```
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123926/suppl/GSE123926_RAW.tar'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113197/suppl/GSE113197_RAW.tar'

mkdir GSE123926
mkdir GSE113197

tar -xf GSE123926_RAW.tar -C ./GSE123926
tar -xf GSE113197_RAW.tar -C ./GSE113197
```
