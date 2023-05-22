#!/bin/bash

#SBATCH --job-name=ModelBSfin
#SBATCH --account=gaydojac0
#SBATCH --partition=standard
#SBATCH --time=12:00:00 

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000m

#SBATCH --output=/dev/null

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=petrofrl@umich.edu


cd /home/petrofrl/Documents/PFAS_allmodels/Part2/array_output

module load R
Rscript --vanilla array_part2FINAL.r
