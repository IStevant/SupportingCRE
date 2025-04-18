# Mapping mouse gonadal supporting cell cis-regulatory elements

## Introduction

The following pipeline was created to analyse paired time-series ATAC and RNA-seq data of purified pre-granulosa and Sertoli cells from fetal gonads.

The analysis consists first in analysing the RNA-seq data by performing sex-differentially expressed gene analysis, embryonic stage differentially expressed genes using DESeq2, and GO term enrichment analysis of the obtained group of genes with ClusterProfiler.

Second, the ATAC-seq analysis consists in performing sex-differentially accessible open chromatin regions, embryonic stage differentially accessible open chromatin regions using DESeq2, and TF motif enrichment analysis of the obtained group of chromatin regions.
Finally, we combine the RNA and the ATAC data to predict cis-regulatory regions and their target genes by correlation analysis with monaLISA.

We also performed ATAC footprint analysis with TOBIAS to predict TF bonding sites and recontruct gene regulatory networks.

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

#### Run the full pipeline

```bash
snakemake --profile=slurm/ 
```

#### Run only the RNA-seq analysis

```bash
snakemake --profile=slurm/ -f RNA_analysis
```

#### Run only the ATAC-seq analysis

```bash
snakemake --profile=slurm/ -f ATAC_analysis
```
