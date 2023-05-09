#!/bin/bash

#SBATCH -p short
#SBATCH –n 4
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH –o %j.out
#SBATCH –e %j.err
#SBATCH --mail-type=ALL
#SBATCH –-mail-user=ssalem5@bwh.harvard.edu

module load gcc/9.2.0 R/4.2.1

./PrepareReference.R
