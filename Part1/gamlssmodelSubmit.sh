#!/bin/bash

#SBATCH --job-name=FitModels
#SBATCH --account=gaydojac0
#SBATCH --partition=standard
#SBATCH --time=24:00:00 

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=48000m

#SBATCH --output=/home/petrofrl/Documents/PFAS_allmodels/output/log_ModelFitSubmit.log
#SBATCH --error=/home/petrofrl/Documents/PFAS_allmodels/output/error_ModelFitSubmit.log

#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=petrofrl@umich.edu


cd /home/petrofrl/Documents/PFAS_allmodels/

module load R
Rscript gamlss.model.r
