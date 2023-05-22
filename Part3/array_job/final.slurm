#!/bin/bash

#SBATCH --job-name=ModelOXBS
#SBATCH --account=gaydojac0
#SBATCH --partition=standard
#SBATCH --time=24:00:00 

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000m

#SBATCH --output=/dev/null

#SBATCH --mail-type=NONE
#SBATCH --mail-user=petrofrl@umich.edu


cd /home/petrofrl/Documents/PFAS_allmodels/Part3/array_job

module load R
Rscript --vanilla array_part3FINAL.r
