#!/bin/bash

#make directory and cd into it for the reference data download
mkdir GEOProfiles
cd GEOProfiles

#download the files for the reference
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz

#unzip them using gunzip since .gz files
gunzip GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz
gunzip GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz
gunzip GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
gunzip GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz
gunzip GSE92742_Broad_LINCS_cell_info.txt.gz
gunzip  GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz 
gunzip GSE92742_Broad_LINCS_sig_info.txt.gz
