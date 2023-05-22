#!/bin/bash
for i in $(seq 1 7); do
  export exposure=$i
  sbatch array.slurm
done
