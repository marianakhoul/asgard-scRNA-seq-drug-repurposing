#!/bin/bash

#SBATCH -p short
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --mail-type=ALL

module load gcc/9.2.0 R/4.2.1
./AnalyzingTestData.R
