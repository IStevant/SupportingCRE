'''
Author: Isabelle Stévant
Affiliation: University of Bar Ilan
Date: 18/09/2024
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


# List of output files
# rule_all_input_list = [
# 	f"{output_pdf}/MULTI_gene2peak_plots.pdf",
# 	f"{output_png}/MULTI_TFBS_motifs_peak_XX_genes.png",
# 	f"{output_png}/MULTI_TFBS_motifs_peak_XY_genes.png",
# 	# f"{processed_data}/plot_example_1.log",
# 	# f"{processed_data}/plot_example_2.log",
# 	f"{processed_data}/plot_example_3.log"
# ]

# Install the necessary R packages using Renv
rule install_packages:
	script:
		"renv/restore.R"

# Run the whole pipeline
rule all:
	input:
		rule_RNA_input_list + rule_ATAC_input_list

# Run the RNA-seq analysis only
rule RNA_analysis:
	input:
		rule_RNA_input_list

# Run the ATAC-seq analysis only
rule ATAC_analysis:
	input:
		rule_ATAC_input_list



################################################################################################
rule MULTI_Get_all_gene_peak_correlation:
	input:
		RNA_samplesheet=f"{processed_data}/RNA_samplesheet.csv",
		ATAC_samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
		RNA_norm_counts=f"{output_tables}/RNA_norm_counts.csv",
		ATAC_norm_counts=f"{processed_data}/ATAC_norm_counts.csv",
		chrom_size=f"{input_data}/chrom.size",
		genes=f"{input_data}/gene_standard.bed",
		gtf=f"{genome}"
	params:
		sex="all",
		distance=config["MULTI_peak_gene_distance"],
		min_cor=config["MULTI_peak_gene_min_cor"],
		FDR=config["MULTI_peak_gene_FDR"]
	output:
		linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/MULTI_gene2peak.R"

rule MULTI_Plot_gene_peak_correlation:
	input:
		genome=f"{genome}",
		gene_bed=f"{input_data}/gene_standard.bed",
		peaks=f"{processed_data}/ATAC_norm_counts.csv",
		linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
		gene_list=config["gene_peak_examples"],
	params:
		bw_folder="results/processed_data/mm10/ATAC_bigwig"
	output:
		pdf=f"{output_pdf}/MULTI_gene2peak_plots.pdf"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/MULTI_plot_gene2peak.R"

rule MULTI_TFBS_motifs_peak_XX_genes:
	input:
		linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
		sex_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj",
		TPM=f"{output_tables}/RNA_TPM.csv"
	params:
		background="conditions",
		genome=config["genome_version"],
		save_folder=f"{output_tables}",
		sex="XX"
	output:
		pdf=f"{output_pdf}/MULTI_TFBS_motifs_peak_XX_genes.pdf",
		png=f"{output_png}/MULTI_TFBS_motifs_peak_XX_genes.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/MULTI_linked_OCR_motif_enrichment.R"

rule MULTI_TFBS_motifs_peak_XY_genes:
	input:
		linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
		sex_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj",
		TPM=f"{output_tables}/RNA_TPM.csv"
	params:
		background="conditions",
		genome=config["genome_version"],
		save_folder=f"{output_tables}",
		sex="XY"
	output:
		pdf=f"{output_pdf}/MULTI_TFBS_motifs_peak_XY_genes.pdf",
		png=f"{output_png}/MULTI_TFBS_motifs_peak_XY_genes.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/MULTI_linked_OCR_motif_enrichment.R"
