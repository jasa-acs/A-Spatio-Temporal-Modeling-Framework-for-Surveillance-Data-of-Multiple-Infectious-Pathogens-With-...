#!/bin/bash

#SBATCH --job-name=balanced_combo1
#SBATCH --account=epi
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xytang@ufl.edu
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=90:00:00
#SBATCH --qos=epi
#SBATCH --array=1-500
#SBATCH --output=mcmc_%a.out

module load gcc
module load R

index_chain=$(((${SLURM_ARRAY_TASK_ID} - 1) % 5 + 1))
index_epi=$(((${SLURM_ARRAY_TASK_ID} - 1) / 5 + 1))

./run_mcmc balanced_combo1 ${index_epi} ${index_chain} 2
Rscript simulation/code/txt2RData.R balanced_combo1 ${index_epi} ${index_chain} 2


