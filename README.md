# NSbirthMass
NSbirthMass is a Python source code containing tools for 1. computing the neutron star birth mass distribution under various accretion modes. 2. performing rapid parameter estimation and model selection for different neutron star mass models.

# Overview
NSbirthMass is designed to offer a wider range of neutron star mass distribution models, faster parameter estimation and model selection. It includes commonly used Gaussian series models, and utilizes the Nested sampling method for parameter estimation. This source code aims to provide a user-friendly and efficient way to conduct calculations. NSbirthMass can be directly downloaded from GitHub and is compatible with the GNU/Linux operating system.

# System Requirements

## Hardware requirements

NSbirthMass source code requires only a standard computer with enough RAM and cores to support parallel computing.

## Software requirements

### OS Requirements
This package is supported for Linux. The source code has been tested on the following systems
#### Linux: Ubuntu 20.04

### Python Dependencies
#### python 3.18.16
#### Slurm
#### bilby==1.1.1
#### dynesty==1.0.1
#### numpy=1.19.5
#### astropy==5.0
#### pandas==1.4.4
#### GalDynPsr
#### galpy==1.8.2
#### scipy==1.10.1
#### matplotlib==3.6.3
#### random
#### corner==2.2.1
#### openpyxl==3.1.2
#### math
#### SciencePlots


# Setting up the development environment:

## To set up Slurm for parallel parameters estimation：
### 1.install Slurm
### 2.switch Python environment to your own in slurm.sh files. For example: cd NSbirthMass/pe_models/ns_obs/2G/ modefiy the last line command in slurm.sh,"python hyper.py", to you own Python environment, "/home/xxx/miniconda/envs/xxx/bin/python hyper.py"

## To run parallel parameters estimation based on task management system:
### cd /NSbirthMass/pemodel/xxx/xxx ; sbatch slurm.sh, or cd /NSbirthMass/pemodel/xxx; bash sbatch 

### Note
#### We employed the Slurm task management system. In this computation, we used single-node parallel processing. When performing parallel calculations, ensure that the "logs/" directory is created before submitting the "slurm.sh" file. The "ntasks-per-node" specifies the number of threads per node for parallel computation. Please note that "npool" should match the "ntasks-per-node" value specified in the "slurm.sh" file. Before submitting a parallel task, make sure to switch the Python environment in "slurm.sh" file, the last line of command, which can be determined using the "which python" command.

## To run notebooks
### cd demos
### run jupyter notebook, xxx.ipynb
