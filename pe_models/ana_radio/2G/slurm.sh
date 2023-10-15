#!/bin/bash
#SBATCH --job-name=ns_radio
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20 
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=800
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

#source ~/miniconda3/etc/profile.d/conda.sh
#conda activate bilby

# The directory logs/ must exist before you submit this file

# Note ntasks-per-node indicates the threads per node for parallel calculation
# Note npool must match the ntasks-per-node given above
# Before submitting a parallel task, make sure to switch the Python environment, which can be determined using the "which python" command.

/home/yzq/miniconda3/envs/ns_mass_test/bin/python hyper.py
