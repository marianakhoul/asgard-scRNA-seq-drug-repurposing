#! /bin/bash

#SBATCH -p short
#SBATCH –n 4
#SBATCH --mem=8G
#SBATCH –o %j.out
#SBATCH –e %j.err
#SBATCH -J bowtie2_run1
#SBATCH --mail-type=ALL
#SBATCH –-mail-user=mfk8@med.harvard.edu
