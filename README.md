# NSbirthMass

**NSbirthMass** is a Python-based toolkit providing source data and tools for:
1. Computing neutron star birth mass distributions using various accretion models.
2. Performing rapid parameter estimation and model selection for neutron star mass models.

---

## Overview

NSbirthMass is designed to provide a diverse range of neutron star mass distribution models and efficient parameter estimation. It includes Gaussian series models and the Nested Sampling method for parameter estimation. This toolkit is user-friendly and efficient for the GNU/Linux operating system.

### Folder Contents
- **`demos`**: Contains partial source data and Python code for generating main datasets and figures.
- **`figures`**: Includes figures and source data from the publication ([DOI: 10.48550/arXiv.2412.05524](https://doi.org/10.48550/arXiv.2412.05524)).
- **`NSmassData`**: Provides neutron star mass data, flux values, AP4 equations of state, and posterior data for models.
- **`pe_models`**: Contains hyperparameter estimation code for 2G, TOP, and other models using observed, analytical, and phenomenological neutron-star mass data.

---

## System Requirements

### Hardware Requirements
- A standard computer with sufficient RAM and CPU cores for parallel computing.

### Software Requirements
#### Supported Operating Systems
- **Linux**: Tested on **Ubuntu 20.04**.

#### Python Dependencies
- All dependencies are listed in the `requirements.txt` file.

---

## Setting Up the Development Environment

### To Set Up Slurm for Parallel Parameter Estimation:

1. **Install Slurm**: Ensure that Slurm is installed on your system.
2. **Configure Python Environment**:
   - Navigate to `NSbirthMass/pe_models/ns_obs/2G/` and locate the `slurm.sh` file.
   - Modify the last line of the file (e.g., `python hyper.py`) to point to your Python environment:
     ```bash
     /home/xxx/miniconda/envs/xxx/bin/python hyper.py
     ```
3. **Generate Neutron Star Mass Data**:
   - Use the mass generation code(Main_obtain_ANA_xxx.ipynb) in the `NSbirthMass/demos/` directory to generate neutron star mass data under different accretion models.
4. **Copy Generated Data**:
   - Move the generated neutron star mass data from `NSbirthMass/demos/` to the parameter estimation folder: `NSbirthMass/pe_models/`. There is already some data (like ana-, obs- and phe- NS mass data) in this folder.

### To Run Parallel Parameter Estimation:

1. Navigate to the desired folder:
   ```bash
   cd /NSbirthMass/pe_models/xxx/yyy
   ```
   - `xxx`: Accretion model (e.g., analytical approach `ana_dns`).
   - `yyy`: Mass model (e.g., Gaussian distribution `G`).

2. Submit the job using Slurm:
   ```bash
   sbatch slurm.sh
   ```
   Alternatively, run the script directly:
   ```bash
   bash sbatch
   ```

**Note**:
- Ensure the `logs/` directory in Mass model folder (/NSbirthMass/pe_models/xxx/yyy) exists before submitting the `slurm.sh` file.
- Match `ntasks-per-node` in the `slurm.sh` file with the `npool` parameter, where npool is the cpu cores pool in cluster.
- Update the Python environment path (/home/xxx/miniconda3/envs/xxx/bin/python) in `slurm.sh` using the `which python` command.

### To Run Jupyter Notebooks:

1. Navigate to the `demos` directory:
   ```bash
   cd demos
   ```
2. Launch the Jupyter Notebook:
   ```bash
   jupyter notebook xxx.ipynb
   ```

---

## License
This project is covered under the **Apache 2.0 License**.

---

## Citation Guide
To cite this work, please refer to [DOI: 10.48550/arXiv.2412.05524](https://doi.org/10.48550/arXiv.2412.05524).
