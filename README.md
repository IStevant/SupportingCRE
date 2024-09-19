# Mapping mouse gonadal supporting cell cis-regulatory elements

## Introduction

The following pipeline was created to analyse paired time-series ATAC and RNA-seq data of purified pre-granulosa and Sertoli cells from fetal gonads.
The analysis consist first in analysing the RNA-seq data by performing sex-differentially expressed gene analysis, embryonic stage differentially expressed genes, and GO term enrichment analysis of the obtained group of genes.
Second, the ATAC-seq analysis consists in performing sex-differentially accessible open chromatin regions, embryonic stage differentially accessible open chromatin regions, and TF motif enrichment analysis of the obtained group of chromatin regions.
Finally, we combine the RNA and the ATAC data to predict cis-regulatory regions and their target genes by correlation analysis.

## How to run the pipeline

### Install Conda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
# Specify your installation directory
```

### Clone the SupportingCRE git project

The following command downloads the last version of the pipeline in the SupportingCRE directory.

```bash
git clone git@github.com:IStevant/SupportingCRE.git
```

### Install the Conda SupportingCRE environment

This installs a conda environment with R4.3.3, Snakemake, as well as all the dependencies necessary to run the pipeline. It can some time to install.

```bash
cd SupportingCRE
conda env create -f conda_env/SupportingCRE_env.yml
```

### Install the R packages

#### Install missing dependencies

```bash
snakemake --core 1 -f install_packages
```

### Configure the analysis parameters

#### The genome version

#### The path to the needed files

#### The different thresholds

### Run the analysis on a HPC using Slurm (recommended)

#### Run the full pipeline

```bash
snakemake --slurm --profile=smk_env/ 
```

#### Run only the RNA-seq analysis

```bash
snakemake --slurm --profile=smk_env/ -f RNA_analysis
```

#### Run only the ATAC-seq analysis

```bash
snakemake --slurm --profile=smk_env/ -f ATAC_analysis
```
