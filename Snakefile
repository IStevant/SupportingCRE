'''
Author: Isabelle Stévant
Affiliation: University of Bar Ilan
Date: 18/09/2024
Licence: MIT

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
# 	f"{output_png}/RNA_corr_pca_all_samples.png",
# 	f"{output_png}/RNA_corr_pca.png",
# 	f"{output_png}/RNA_marker_genes.png",
# 	f"{output_png}/RNA_marker_genes_enrichment.png",
# 	f"{output_png}/RNA_sex_DEG_histograms.png",
# 	f"{output_png}/RNA_sex_DEG_volcano.png",
# 	f"{output_png}/RNA_sex_DEG_upset.png",
# 	f"{output_png}/RNA_XX_DEG_stage_heatmap.png",
# 	f"{output_png}/RNA_XY_DEG_stage_heatmap.png",
# 	f"{output_png}/RNA_sex_stage_common_DEGs.png",
# 	f"{output_png}/ATAC_corr_pca_all_samples.png",
# 	f"{output_png}/ATAC_consensus_peak_distribution.png",
# 	f"{output_png}/ATAC_all_consensus_peak_annotation.png",
# 	f"{output_png}/ATAC_sex_DAR_histograms.png",
# 	f"{output_png}/ATAC_sig_sex_DARs_annotation.png",
# 	f"{output_png}/ATAC_sex_DAR_upset.png",
# 	f"{output_pdf}/ATAC_sex_DAR_TF_motifs_rdm_bg.pdf",
# 	f"{output_pdf}/ATAC_sex_DAR_TF_motifs_sex_bg.pdf",
# 	f"{output_png}/ATAC_sex_DAR_TF_motifs_merged.png",
# 	f"{output_png}/ATAC_XX_DAR_stage_heatmap.png",
# 	f"{output_png}/ATAC_XY_DAR_stage_heatmap.png",
# 	f"{output_png}/ATAC_XX_stage_DAR_TF_motifs_sex_bg.png",
# 	f"{output_png}/ATAC_XY_stage_DAR_TF_motifs_sex_bg.png",
# 	f"{output_png}/ATAC_XX_stage_DAR_TF_motifs_random_bg.png",
# 	f"{output_png}/ATAC_XY_stage_DAR_TF_motifs_random_bg.png",
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











rule ATAC_TFBS_motifs_sex_rdm_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		TPM=f"{output_tables}/RNA_TPM.csv"
	params:
		minTPM=config["RNA_minTPM"],
		background="genome",
		genome=config["genome_version"],
		save_folder=f"{output_tables}"
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_rdm_bg.pdf",
		png=f"{output_png}/ATAC_sex_DAR_TF_motifs_rdm_bg.png"
	resources:
		cpus_per_task=24,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_sex_DAR_motif_enrich.R"

rule ATAC_TFBS_motifs_sex_cond_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		TPM=f"{output_tables}/RNA_TPM.csv"
	params:
		minTPM=config["RNA_minTPM"],
		background="conditions",
		genome=config["genome_version"],
		save_folder=f"{output_tables}"
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_sex_bg.pdf",
		png=f"{output_png}/ATAC_sex_DAR_TF_motifs_sex_bg.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_sex_DAR_motif_enrich.R"

rule ATAC_TFBS_motifs_sex_DAR_merged_stages:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		TPM=f"{output_tables}/RNA_TPM.csv"
	params:
		minTPM=config["RNA_minTPM"],
		background="conditions",
		genome=config["genome_version"],
		save_folder=f"{output_tables}"
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_merged.pdf",
		png=f"{output_png}/ATAC_sex_DAR_TF_motifs_merged.png"
	resources:
		cpus_per_task=24,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_sex_DAR_motif_enrich_merge_stages.R"


rule ATAC_Get_XX_dynamic_DARs:
	input:
		counts=f"{processed_data}/ATAC_raw_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
		peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj",
		gtf=f"{genome}"
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		promoter=config["ATAC_promoter_distance"],
		sex="XX"
	output:
		tsv=f"{output_tables}/ATAC_XX_DAR_stage.tsv",
		sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XX.Robj"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_stage_DAR.R"

# rule ATAC_Plot_XX_stage_DAR_peak_examples:
# 	input:
# 		genome=f"{genome}",
# 		peaks=f"{processed_data}/ATAC_norm_counts.csv",
# 		peak_list=config["DAR_peak_examples_XX"],
# 	params:
# 		bw_folder="results/processed_data/mm10/ATAC_bigwig",
# 		save_folder=f"{output_png}",
# 		sex="XX"
# 	output: 
# 		log= f"{processed_data}/plot_example_2.log"
# 	resources:
# 		cpus_per_task=12,
# 		mem_mb=64000
# 	script:
# 		"workflow/scripts/ATAC_plot_DAR_peak_examples.R"


rule ATAC_Get_XY_dynamic_DARs:
	input:
		counts=f"{processed_data}/ATAC_raw_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
		peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj",
		gtf=f"{genome}"
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		promoter=config["ATAC_promoter_distance"],
		sex="XY"
	output:
		tsv=f"{output_tables}/ATAC_XY_DAR_stage.tsv",
		sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XY.Robj"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_stage_DAR.R"

# rule ATAC_Plot_XY_stage_DAR_peak_examples:
# 	input:
# 		genome=f"{genome}",
# 		peaks=f"{processed_data}/ATAC_norm_counts.csv",
# 		peak_list=config["DAR_peak_examples_XY"],
# 	params:
# 		bw_folder="results/processed_data/mm10/ATAC_bigwig",
# 		save_folder=f"{output_png}",
# 		sex="XY"
# 	output: 
# 		log= f"{processed_data}/plot_example_3.log"
# 	resources:
# 		cpus_per_task=12,
# 		mem_mb=64000
# 	script:
# 		"workflow/scripts/ATAC_plot_DAR_peak_examples.R"


rule ATAC_Plot_heatmap_dyn_DARs_XX:
	input:
		sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XX.Robj",
		norm_counts=f"{processed_data}/ATAC_norm_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv"
	params:
		sex="XX",
		clusters=config["ATAC_XX_stage_DAR_clusters"]
	output:
		clusters=f"{output_tables}/ATAC_XX_DAR_stage_heatmap_clusters.csv",
		pdf=f"{output_pdf}/ATAC_XX_DAR_stage_heatmap.pdf",
		png=f"{output_png}/ATAC_XX_DAR_stage_heatmap.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_stage_DAR_heatmap.R"

rule ATAC_Plot_heatmap_dyn_DARs_XY:
	input:
		sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XY.Robj",
		norm_counts=f"{processed_data}/ATAC_norm_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv"
	params:
		sex="XY",
		clusters=config["ATAC_XY_stage_DAR_clusters"]
	output:
		clusters=f"{output_tables}/ATAC_XY_DAR_stage_heatmap_clusters.csv",
		pdf=f"{output_pdf}/ATAC_XY_DAR_stage_heatmap.pdf",
		png=f"{output_png}/ATAC_XY_DAR_stage_heatmap.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_stage_DAR_heatmap.R"


rule ATAC_XX_TFBS_motifs_stage_cond_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{output_tables}/ATAC_XX_DAR_stage_heatmap_clusters.csv",
		TPM=f"{output_tables}/RNA_TPM.csv"
	params:
		minTPM=config["RNA_minTPM"],
		background="conditions",
		logos="TRUE",
		nbTFs=8,
		genome=config["genome_version"],
		save_folder=f"{output_tables}",
		sex="XX"
	output:
		pdf=f"{output_pdf}/ATAC_XX_stage_DAR_TF_motifs_sex_bg.pdf",
		png=f"{output_png}/ATAC_XX_stage_DAR_TF_motifs_sex_bg.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_stage_DAR_motif_enrich.R"

rule ATAC_XY_TFBS_motifs_stage_cond_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{output_tables}/ATAC_XY_DAR_stage_heatmap_clusters.csv",
		TPM=f"{output_tables}/RNA_TPM.csv"
	params:
		minTPM=config["RNA_minTPM"],
		background="conditions",
		logos="TRUE",
		nbTFs=8,
		genome=config["genome_version"],
		save_folder=f"{output_tables}",
		sex="XY"
	output:
		pdf=f"{output_pdf}/ATAC_XY_stage_DAR_TF_motifs_sex_bg.pdf",
		png=f"{output_png}/ATAC_XY_stage_DAR_TF_motifs_sex_bg.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_stage_DAR_motif_enrich.R"

rule ATAC_XX_TFBS_motifs_stage_rand_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{output_tables}/ATAC_XX_DAR_stage_heatmap_clusters.csv",
		TPM=f"{output_tables}/RNA_TPM.csv"
	params:
		minTPM=config["RNA_minTPM"],
		background="genome",
		logos="TRUE",
		nbTFs=8,
		genome=config["genome_version"],
		save_folder=f"{output_tables}",
		sex="XX"
	output:
		pdf=f"{output_pdf}/ATAC_XX_stage_DAR_TF_motifs_random_bg.pdf",
		png=f"{output_png}/ATAC_XX_stage_DAR_TF_motifs_random_bg.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_stage_DAR_motif_enrich.R"

rule ATAC_XY_TFBS_motifs_stage_rand_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{output_tables}/ATAC_XY_DAR_stage_heatmap_clusters.csv",
		TPM=f"{output_tables}/RNA_TPM.csv"
	params:
		minTPM=config["RNA_minTPM"],
		background="genome",
		logos="TRUE",
		nbTFs=8,
		genome=config["genome_version"],
		save_folder=f"{output_tables}",
		sex="XY"
	output:
		pdf=f"{output_pdf}/ATAC_XY_stage_DAR_TF_motifs_random_bg.pdf",
		png=f"{output_png}/ATAC_XY_stage_DAR_TF_motifs_random_bg.png"
	resources:
		cpus_per_task=12,
		mem_mb=80000
	script:
		"workflow/scripts/ATAC_stage_DAR_motif_enrich.R"

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
