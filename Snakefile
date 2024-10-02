'''
Author: Isabelle Stévant
Affiliation: University of Bar Ilan
Date: 24/09/2024
Licence: MIT

Snakemake rules to run the full pipeline.

The parameters of the analsysis are defined in the analysis_parameters.yaml configuration file.

Pipeline created to analyse the paired time series RNA and ATAC-seq of the supporting cells from mouse fetal gonads.
Stévant et al. 2024
doi:

'''

# Generate the report
# report: "report/workflow.rst"

# Import the different pipeline modules
include: "workflow/rules/RNA_analysis.smk"
include: "workflow/rules/ATAC_analysis.smk"
include: "workflow/rules/MULTI_analysis.smk"
include: "workflow/rules/TOBIAS_analysis.smk"

###########################################
#                                         #
#                  Rules                  #
#                                         #
###########################################

# Run the whole pipeline
rule all:
	input:
		rule_RNA_input_list + rule_ATAC_input_list + rule_MULTI_input_list + rule_TOBIAS_input_list

# Install the necessary R packages using Renv
rule install_packages:
	script:
		"renv/restore.R"

# Run the RNA-seq analysis only
rule RNA_analysis:
	input:
		rule_RNA_input_list

# Run the ATAC-seq analysis only
rule ATAC_analysis:
	input:
		rule_ATAC_input_list

# Run the Multiomics analysis only
rule MULTI_analysis:
	input:
		rule_MULTI_input_list
