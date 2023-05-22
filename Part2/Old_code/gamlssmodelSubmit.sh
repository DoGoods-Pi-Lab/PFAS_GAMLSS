#!/bin/bash

#SBATCH --job-name=ModelBSonly
#SBATCH --account=gaydojac0
#SBATCH --partition=standard
#SBATCH --time=18:30:00 

#SBATCH --nodes=1
#SBATCH --cpus-per-task=11
#SBATCH --mem-per-cpu=5000m

#SBATCH --output=/home/petrofrl/Documents/PFAS_allmodels/Part2/output/log1_ModelFitSubmit.log
#SBATCH --error=/home/petrofrl/Documents/PFAS_allmodels/Part2/output/error1_ModelFitSubmit.log

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=petrofrl@umich.edu


cd /home/petrofrl/Documents/PFAS_allmodels/Part2

module load R
Rscript --vanilla part2.r
