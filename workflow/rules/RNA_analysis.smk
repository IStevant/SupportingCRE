'''
Author: Isabelle Stévant
Affiliation: University of Bar Ilan
Date: 18/09/2024
Licence: MIT


Snakemake rules concerning the RNA-seq analysis

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
rule_RNA_input_list = [
	f"{output_png}/RNA_corr_pca_all_samples.png",
	f"{output_png}/RNA_corr_pca.png",
	f"{output_png}/RNA_marker_genes.png",
	f"{output_png}/RNA_marker_genes_enrichment.png",
	f"{output_png}/RNA_sex_DEG_histograms.png",
	f"{output_png}/RNA_sex_DEG_volcano.png",
	f"{output_png}/RNA_sex_DEG_upset.png",
	f"{output_png}/RNA_XX_DEG_stage_heatmap.png",
	f"{output_png}/RNA_XY_DEG_stage_heatmap.png",
	f"{output_png}/RNA_sex_stage_common_DEGs.png",
]

# If there is no outliers, do not run the analysis that discard them
if len(config["RNA_outliers"])<1:
	rule_RNA_input_list.remove(f"{output_png}/RNA_corr_pca_all_samples.png")


## Uncomment if you want to run only this pipeline
# rule all:
# 	input:
# 		rule_RNA_input_list

rule RNA_Get_matrices:
	input:
		counts=f'{input_data}/{config["RNA_counts"]}',
		tpm=f'{input_data}/{config["RNA_TPM"]}',
		protein_genes=config["protein_genes"]
	params:
		minReads=config["RNA_minReads"],
		minTPM=config["RNA_minTPM"],
		RNA_outliers=config["RNA_outliers"]
	output:
		tpm_all=f"{processed_data}/RNA_TPM_all_samples.csv",
		tpm=f"{output_tables}/RNA_TPM.csv",
		counts=f"{output_tables}/RNA_raw_counts.csv",
		norm_counts=f"{output_tables}/RNA_norm_counts.csv",
		norm_counts_all=f"{processed_data}/RNA_norm_counts_all.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_clean_matrices.R"

rule RNA_corr_PCA_with_outliers:
	input:
		norm_data=f"{processed_data}/RNA_norm_counts_all.csv"
	params:
		corr_method=config["RNA_corr_met"]
	output:
		pdf=f"{output_pdf}/RNA_corr_pca_all_samples.pdf",
		png=f"{output_png}/RNA_corr_pca_all_samples.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/Corr_pca.R"

rule RNA_corr_PCA:
	input:
		norm_data=f"{output_tables}/RNA_norm_counts.csv"
	params:
		corr_method="spearman"
	output:
		pdf=f"{output_pdf}/RNA_corr_pca.pdf",
		png=f"{output_png}/RNA_corr_pca.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/Corr_pca.R"

rule RNA_Plot_marker_genes:
	input:
		marker_genes=config["marker_genes"],
		whole_gonad=f'{input_data}/{config["whole_gonad_RNAseq"]}',
		tpm=f"{output_tables}/RNA_TPM.csv"
	output:
		pdf1=f"{output_pdf}/RNA_marker_genes.pdf",
		png1=f"{output_png}/RNA_marker_genes.png",
		pdf2=f"{output_pdf}/RNA_marker_genes_enrichment.pdf",
		png2=f"{output_png}/RNA_marker_genes_enrichment.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_plot_marker_genes.R"

rule RNA_Get_sex_DEGs:
	input:
		counts=f"{output_tables}/RNA_raw_counts.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv",
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"]
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		save_folder=f"{output_tables}"
	output:
		all_DEGs=f"{processed_data}/RNA_all_SexDEGs.Robj",
		sig_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_sex_DEG.R"

rule RNA_Plot_sex_DEG_histogram:
	input:
		sig_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	output:
		pdf=f"{output_pdf}/RNA_sex_DEG_histograms.pdf",
		png=f"{output_png}/RNA_sex_DEG_histograms.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_plot_sex_DEG_hist.R"

rule RNA_Plot_sex_DEG_volcano_GO:
	input:
		all_DEGs=f"{processed_data}/RNA_all_SexDEGs.Robj",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		path=f"{output_tables}"
	output:
		pdf=f"{output_pdf}/RNA_sex_DEG_volcano.pdf",
		png=f"{output_png}/RNA_sex_DEG_volcano.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_plot_sex_volcano_GO.R"

rule RNA_Plot_sex_DEG_upset:
	input:
		sig_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj"
	params:
		output_folder=f"{output_tables}/"
	output:
		pdf=f"{output_pdf}/RNA_sex_DEG_upset.pdf",
		png=f"{output_png}/RNA_sex_DEG_upset.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_sex_DEG_upset.R"

rule RNA_Get_XX_dynamic_DEGs:
	input:
		counts=f"{output_tables}/RNA_raw_counts.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv",
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"]
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		sex="XX"
	output:
		tsv=f"{output_tables}/RNA_XX_DEG_stage.tsv",
		sig_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XX.Robj"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_stage_DEG.R"

rule RNA_Get_XY_dynamic_DEGs:
	input:
		counts=f"{output_tables}/RNA_raw_counts.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv",
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"]
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		sex="XY"
	output:
		tsv=f"{output_tables}/RNA_XY_DEG_stage.tsv",
		sig_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XY.Robj"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_stage_DEG.R"

rule RNA_Plot_heatmap_GO_XX:
	input:
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"],
		sig_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XX.Robj",
		norm_counts=f"{output_tables}/RNA_norm_counts.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	params:
		sex="XX",
		clusters=config["RNA_XX_stage_DEG_clusters"]
	output:
		GO=f"{output_tables}/RNA_XX_GO_DEG_stage.csv",
		cluster_file=f"{output_tables}/RNA_XX_DEG_stage_heatmap_clusters.csv",
		pdf=f"{output_pdf}/RNA_XX_DEG_stage_heatmap.pdf",
		png=f"{output_png}/RNA_XX_DEG_stage_heatmap.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_stage_DEG_heatmap.R"

rule RNA_Plot_heatmap_GO_XY:
	input:
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"],
		sig_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XY.Robj",
		norm_counts=f"{output_tables}/RNA_norm_counts.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	params:
		sex="XY",
		clusters=config["RNA_XY_stage_DEG_clusters"]
	output:
		GO=f"{output_tables}/RNA_XY_GO_DEG_stage.csv",
		cluster_file=f"{output_tables}/RNA_XY_DEG_stage_heatmap_clusters.csv",
		pdf=f"{output_pdf}/RNA_XY_DEG_stage_heatmap.pdf",
		png=f"{output_png}/RNA_XY_DEG_stage_heatmap.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_stage_DEG_heatmap.R"

rule RNA_Plot_sex_stage_common_DEGs:
	input:
		sex_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj",
		XY_stage_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XY.Robj",
		XX_stage_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XX.Robj",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	output:
		pdf=f"{output_pdf}/RNA_sex_stage_common_DEGs.pdf",
		png=f"{output_png}/RNA_sex_stage_common_DEGs.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/RNA_overlap_sex_stage_DEG.R"
