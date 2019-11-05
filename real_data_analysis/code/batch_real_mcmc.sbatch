#!/bin/bash

#SBATCH --job-name=case_study
#SBATCH --account=epi
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xytang@ufl.edu
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=90:00:00
#SBATCH —-qos=epi
#SBATCH --array=1-10
#SBATCH —-output=real_mcmc_%a.out

module load gcc
module load R

./run_real_mcmc ${SLURM_ARRAY_TASK_ID}
Rscript real_data_analysis/code/txt2RData.R ${SLURM_ARRAY_TASK_ID}


