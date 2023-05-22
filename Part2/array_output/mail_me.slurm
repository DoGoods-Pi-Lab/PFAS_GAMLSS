#!/bin/bash

#SBATCH --job-name=ModelBSonly
#SBATCH --account=gaydojac0
#SBATCH --partition=standard
#SBATCH --time=1:00

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000m

#SBATCH --output=/dev/null

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=petrofrl@umich.edu

