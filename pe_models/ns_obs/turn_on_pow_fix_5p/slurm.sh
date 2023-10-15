#!/bin/bash
#SBATCH --job-name=hyper
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=500
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate bilby

# The directory logs/ must exist before you submit this file

# Note npool must match the ntasks-per-node given above

/home/yzq/miniconda3/envs/bilby/bin/python hyper.py
