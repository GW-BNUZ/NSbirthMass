#!/bin/bash
#SBATCH --job-name=hyper
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=500
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

#source ~/anaconda3/etc/profile.d/conda.sh
#conda activate bilby

# The directory logs/ must exist before you submit this file

# Note npool must match the ntasks-per-node given above

/data/users/yzq/anaconda3/envs/bilby/bin/python hyper.py
