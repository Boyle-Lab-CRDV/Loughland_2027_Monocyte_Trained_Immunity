#!/bin/bash

#SBATCH --job-name="QC_atac_sabah1"
#SBATCH --cpus-per-task=8
#SBATCH --mem=100000
#SBATCH --partition=long                    
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=zuleima.pava@burnet.edu.au
#SBATCH --time=00-48:00:00


#Load necesary module
module load R/4.4

#Execute the script
Rscript /working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/ATACseq_QC_sabah1_m.R