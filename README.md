# SupportingCRE
## Install SupportingCRE conda environment
This will install an environment with R4.3.3 and Snakemake, as well as all the dependencies necessary to run the pipeline. It will take some time to install.
```bash
￼
cd SupportingCRE
conda env create -f conda_env/SupportingCRE_env.yml
```
## Activate the conda environment
```bash
￼
conda activate SupportingCRE
```
## Run the pipeline 
### Install all R packages
The first time the pipeline is launched, need R to check and install all the necessary packages to run the analysis. This will take some time to install.
```bash
￼
snakemake --cores 1 -f install_packages
```
### Run the analysis
Once all the packages are installed, the pipeline can be run using the following command:
```bash
￼
snakemake --cores 8
```