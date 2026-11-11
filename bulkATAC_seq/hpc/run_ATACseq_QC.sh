#!/bin/bash

#SBATCH --job-name="QC_atacseq"
#SBATCH --cpus-per-task=8
#SBATCH --mem=128000
#SBATCH --partition=long                    
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=zuleima.pava@burnet.edu.au
#SBATCH --time=00-96:00:00


#Load necesary module
module load R/4.4

#Execute the script
Rscript /working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/ATACseq_QC2.R