'''
Author: Isabelle Stévant
Affiliation: University of Bar Ilan
Date: 18/09/2024
Licence: MIT


Snakemake rules concerning the ATAC-seq analysis

The pipeline support mm10 and mm39 mouse genome versions
The genome version is specified in the smk_env/workflow_config.yaml configuration file.

Pipeline created to analyse the paired time series RNA and ATAC-seq of the supporting cells from mouse fetal gonads.
Stévant et al. 2024
doi:

'''

configfile: "smk_env/workflow_config.yaml"


# Get file path according to the genome version
input_data = f'{config["path_to_data"]}{config["genome_version"]}'
processed_data = f'{config["path_to_process"]}{config["genome_version"]}'
output_png = f'{config["path_to_graphs"]}{config["genome_version"]}/PNG'
output_pdf = f'{config["path_to_graphs"]}{config["genome_version"]}/PDF'
output_tables = f'{config["path_to_tables"]}{config["genome_version"]}'


# List of output files
rule_ATAC_input_list = [
	f"{output_png}/ATAC_corr_pca_all_samples.png",
	f"{output_png}/ATAC_consensus_peak_distribution.png",
	f"{output_png}/ATAC_all_consensus_peak_annotation.png",
	f"{output_png}/ATAC_sex_DAR_histograms.png",
	f"{output_png}/ATAC_sig_sex_DARs_annotation.png",
	f"{output_png}/ATAC_sex_DAR_upset.png",
]


if config["genome_version"] == "mm10":
	genome = f'{input_data}/gencode.vM25.annotation.gtf.gz'
	ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm10"]
else :
	genome = f'{input_data}/gencode.vM34.annotation.gtf.gz'
	ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm39"]

## Uncomment if you want to run only this pipeline
# rule all:
# 	input:
# 		rule_ATAC_input_list

rule ATAC_Get_matrices:
	input:
		counts=f'{input_data}/{config["ATAC_counts"]}',
	params:
		minReads=config["ATAC_minReads"],
	output:
		counts=f"{processed_data}/ATAC_raw_counts.csv",
		norm_counts=f"{processed_data}/ATAC_norm_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
		size_factors=f"{processed_data}/ATAC_size_factors.csv"
	resources:
		cpus_per_task=4,
		mem_mb=16000
	script:
		"../scripts/ATAC_clean_matrices.R"

# Force the execution if you need to re-normalize the bigwig files
rule ATAC_Normalize_bigwig:
	input:
		size_factors=f"{processed_data}/ATAC_size_factors.csv"
	params:
		bigwig_folder=f"{ATAC_bigwig_folder}",
		new_bigwig_folder=f"{processed_data}/ATAC_bigwig"
	output:
		output_file=f"{processed_data}/ATAC_bigwig/scale/ATAC_size_factors.csv"
	resources:
		cpus_per_task=12,
		mem_mb=16000
	script:
		"../scripts/MULTI_norm_bigwig.R"

rule ATAC_Plot_nb_consensus_peak:
	input:
		peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
		norm_data=f"{processed_data}/ATAC_norm_counts.csv"
	output:
		pdf=f"{output_pdf}/ATAC_consensus_peak_distribution.pdf",
		png=f"{output_png}/ATAC_consensus_peak_distribution.png"
	resources:
		cpus_per_task=4,
		mem_mb=16000
	script:
		"../scripts/ATAC_peak_distribution_per_sex.R"

rule ATAC_Plot_consensus_peak_annotation:
	input:
		genome=f"{genome}",
		peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj"
	params:
		promoter=config["ATAC_promoter_distance"]
	output:
		anno_list=f"{processed_data}/ATAC_all_consensus_peak_annotation.Robj",
		pdf=f"{output_pdf}/ATAC_all_consensus_peak_annotation.pdf",
		png=f"{output_png}/ATAC_all_consensus_peak_annotation.png"
	resources:
		cpus_per_task=4,
		mem_mb=16000
	script:
		"../scripts/ATAC_peak_annotation_per_sex.R"

rule ATAC_corr_PCA:
	input:
		norm_data=f"{processed_data}/ATAC_norm_counts.csv"
	params:
		corr_method=config["ATAC_corr_met"]
	output:
		pdf=f"{output_pdf}/ATAC_corr_pca_all_samples.pdf",
		png=f"{output_png}/ATAC_corr_pca_all_samples.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/Corr_pca.R"

# rule ATAC_Plot_peak_examples:
# 	input:
# 		genome=f"{genome}",
# 		gene_bed=f"{input_data}/gene_standard.bed",
# 		peaks=f"{processed_data}/ATAC_norm_counts.csv",
# 		linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
# 		gene_list=config["peak_examples"],
# 	params:
# 		bw_folder="results/processed_data/mm10/ATAC_bigwig",
# 		save_folder=f"{output_png}"
# 	output: 
# 		log= f"{processed_data}/plot_example_1.log"
# 	resources:
# 		cpus_per_task=12,
# 		mem_mb=64000
# 	script:
# 		"../scripts/MULTI_plot_genomic_tracks_examples.R"


rule ATAC_Get_sex_DARs:
	input:
		counts=f"{processed_data}/ATAC_raw_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
		gtf=f"{genome}"
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		promoter=config["ATAC_promoter_distance"],
		save_folder=f"{output_tables}"
	output:
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		sig_DARs_GR=f"{processed_data}/ATAC_sig_SexDARs_GR.Robj"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/ATAC_sex_DAR.R"

rule ATAC_Plot_sex_DAR_histogram:
	input:
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv"
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_histograms.pdf",
		png=f"{output_png}/ATAC_sex_DAR_histograms.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/ATAC_plot_sex_DAR_hist.R"

rule ATAC_Plot_sex_DAR_peak_annotation:
	input:
		genome=f"{genome}",
		peak_list=f"{processed_data}/ATAC_sig_SexDARs_GR.Robj"
	params:
		promoter=config["ATAC_promoter_distance"]
	output:
		anno_list=f"{processed_data}/ATAC_sig_SexDARs_annotation.Robj",
		pdf=f"{output_pdf}/ATAC_sig_sex_DARs_annotation.pdf",
		png=f"{output_png}/ATAC_sig_sex_DARs_annotation.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/ATAC_peak_annotation_per_sex.R"

rule ATAC_Plot_sex_DAR_upset:
	input:
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj"
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_upset.pdf",
		png=f"{output_png}/ATAC_sex_DAR_upset.png"
	resources:
		cpus_per_task=24,
		mem_mb=64000
	script:
		"../scripts/ATAC_sex_DAR_upset.R"
