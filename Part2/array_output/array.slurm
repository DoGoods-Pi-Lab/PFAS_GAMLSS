#!/bin/bash

#SBATCH --job-name=ModelBSonly
#SBATCH --account=gaydojac0
#SBATCH --partition=standard
#SBATCH --time=18:30:00 

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000m

#SBATCH --output=/dev/null

#SBATCH --mail-type=NONE
#SBATCH --mail-user=petrofrl@umich.edu
#SBATCH --array=1-62


cd /home/petrofrl/Documents/PFAS_allmodels/Part2

module load R
Rscript --vanilla array_part2.r
