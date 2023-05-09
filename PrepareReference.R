#!/usr/bin/env Rscript

library('Asgard')

#Please replace Your_local_path with your real local folder
#Also, the names of these files may be different after unzipping, change them when unzipped.
PrepareReference(cell.info="GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "DrugReference/")# change this directory and the above file locations

#Please note that it takes more than one hour to produce drug references in a standard computer with RAM>64GB.
