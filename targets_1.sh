#!/bin/bash

#SBATCH --job-name=test1
#SBATCH --mail-user=paul.kruse@nih.gov
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm/t1_%j.err
#SBATCH --output=slurm/t1_%j.out

Rscript -e "targets::tar_make()"