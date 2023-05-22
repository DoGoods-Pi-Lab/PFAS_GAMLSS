#!/bin/bash

#SBATCH --job-name=Interaction-TermsOXBS
#SBATCH --account=gaydojac0
#SBATCH --partition=standard
#SBATCH --time=12:30:00 

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000m

#SBATCH --output=/home/petrofrl/Documents/PFAS_allmodels/Part4/output/output.log
#SBATCH --error=/home/petrofrl/Documents/PFAS_allmodels/Part4/output/error.log

#SBATCH --mail-type=end,FAIL
#SBATCH --mail-user=petrofrl@umich.edu


cd /home/petrofrl/Documents/PFAS_allmodels/Part4

module load R
Rscript --vanilla oxbs_independent.r
