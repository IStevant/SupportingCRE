'''
Author: Isabelle Stévant
Affiliation: Mammalian sex determination lab, University of Bar Ilan
Date: 19/09/2024
Licence: MIT


Snakemake rules concerning the ATAC-seq analysis

The pipeline support mm10 and mm39 mouse genome versions
The genome version is specified in the smk_env/workflow_config.yaml configuration file.

Pipeline created to analyse the paired time series RNA and ATAC-seq of the supporting cells from mouse fetal gonads.
Stévant et al. 2024
doi:

'''

configfile: "smk_env/workflow_config.yaml"

sexes = config["sexes"]
stages = config["stages"]
TFBS_background = ["genome", "conditions"]

# Get file path according to the genome version
input_data = f'{config["path_to_data"]}{config["genome_version"]}'
processed_data = f'{config["path_to_process"]}{config["genome_version"]}/ATAC'
output_png = f'{config["path_to_graphs"]}{config["genome_version"]}/ATAC'
output_pdf = f'{config["path_to_graphs"]}{config["genome_version"]}/ATAC/PDF'
output_tables = f'{config["path_to_tables"]}{config["genome_version"]}/ATAC'

if config["genome_version"] == "mm10":
	genome = f'{input_data}/gencode.vM25.annotation.gtf.gz'
	ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm10"]
else :
	genome = f'{input_data}/gencode.vM34.annotation.gtf.gz'
	ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm39"]

# List of output files
rule_ATAC_input_list = [
	f"{output_png}/ATAC_corr_pca_all_samples.png",
	f"{output_png}/ATAC_consensus_peak_distribution.png",
	f"{output_png}/ATAC_all_consensus_peak_annotation.png",
	f"{output_png}/ATAC_sex_DAR_histograms.png",
	f"{output_png}/ATAC_sig_sex_DARs_annotation.png",
	f"{output_png}/ATAC_sex_DAR_upset.png",
	f"{output_pdf}/ATAC_sex_DAR_TF_motifs_rdm_bg.pdf",
	f"{output_pdf}/ATAC_sex_DAR_TF_motifs_sex_bg.pdf",
	f"{output_png}/ATAC_sex_DAR_TF_motifs_merged.png",
	f"{output_png}/ATAC_XX_DAR_stage_heatmap.png",
	f"{output_png}/ATAC_XY_DAR_stage_heatmap.png",
	f"{output_png}/ATAC_XX_stage_DAR_TF_motifs_sex_bg.png",
	f"{output_png}/ATAC_XY_stage_DAR_TF_motifs_sex_bg.png",
	f"{output_png}/ATAC_XX_stage_DAR_TF_motifs_random_bg.png",
	f"{output_png}/ATAC_XY_stage_DAR_TF_motifs_random_bg.png"
]

## Uncomment if you want to run only this pipeline
# rule all:
# 	input:
# 		rule_ATAC_input_list


###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

# Generate the expression matrices that will be used for the rest of the analysis
rule ATAC_Get_matrices:
	input:
		counts=f'{input_data}/{config["ATAC_counts"]}', # Raw read count per peak comming from the nf-core/atacseq pipeline
	params:
		minReads=config["ATAC_minReads"],  # Minimal number of raw reads from which we consider a peak accessible
	output:
		counts=f"{output_tables}/ATAC_raw_counts.csv",          # Filtered raw read counts
		norm_counts=f"{output_tables}/ATAC_norm_counts.csv",    # Filtered normalized read counts (normalization by the library size)
		samplesheet=f"{output_tables}/ATAC_samplesheet.csv",    # Description of the samples for downstream analysis
		size_factors=f"{output_tables}/ATAC_size_factors.csv"   # Size factors used to normalize the data (was used to normalize the BigWig files)
	resources:
		cpus_per_task=1,
		mem_mb=4000
	script:
		"../scripts/ATAC_clean_matrices.R"

# BigWig files from the nf-core/atacseq pipeline are normalized by read per million.
# This normalization is not appropriate when samples have different FRiP scores (i.e. different level of background noise).
# We decided to normalize the BigWig using the same normalization applied for the peaks in the differential accessibility analysis.
# The rational is that most of the peaks are unchanged between conditions, so the signal can be normalized using the peak signals only.
# We provide the normalized bigwig files on GEO.
# Force the execution if you need to re-normalize the bigwig files (not recommanded, apply only on output nf-core/atacseq files normalized by read per million).
rule ATAC_Normalize_bigwig:
	input:
		size_factors=f"{output_tables}/ATAC_size_factors.csv"   # Size factors used to normalize the data
	params:
		bigwig_folder=f"{ATAC_bigwig_folder}",             # Path to the original nf-core/atacseq BigWig files
		new_bigwig_folder=f"{processed_data}/ATAC_bigwig"  # Path to the normalized files
	output:
		output_file=f"{processed_data}/ATAC_bigwig/scale/ATAC_size_factors.csv"  # Save the size factors with the new BigWig (to signify snakemake that the rule executed well)
	resources:
		cpus_per_task=12,
		mem_mb=16000
	script:
		"../scripts/MULTI_norm_bigwig.R"

# Get differentially accessible regions between XX and XY at each embryonic stage using DESeq2 Wald test
rule ATAC_Get_sex_DARs:
	input:
		counts=f"{output_tables}/ATAC_raw_counts.csv",        # Filtered raw read counts
		samplesheet=f"{output_tables}/ATAC_samplesheet.csv",  # Description of the samples
		gtf=f"{genome}"                                       # GTF file
	params:
		adjpval=config["ATAC_adjpval"],             # Adjusted p-value threshold to condider a region differentially accessible
		log2FC=config["ATAC_log2FC"],               # Log2 fold change threshold to condider a region differentially accessible
		promoter=config["ATAC_promoter_distance"],  # distance to the TSS defining the promoter region
		save_folder=f"{output_tables}"              # Location where the result tables are saved
	output:
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",      # R object containing the filtered DESeq2 results
		sig_DARs_GR=f"{processed_data}/ATAC_sig_SexDARs_GR.Robj" # R object containing DAR regions as a GRanges list
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"../scripts/ATAC_sex_DAR.R"

# Get TFBS enrichment in sex-biased OCRs against random genomic background and against the opposite sex as background.
# Retirns the plots and the enrichment scores as tables.
rule ATAC_TFBS_motifs_sex_DAR:
	input:
		TF_genes=config["TF_genes"],                         # List of known mouse transcripotion factors
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",  # R object containing the filtered DESeq2 results
		TPM=f"{output_tables}/TPM.csv"                       # Gene expression matrix
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_{{background}}_bg.pdf",  # Figure as PDF
		png=f"{output_png}/ATAC_sex_DAR_TF_motifs_{{background}}_bg.png"   # Figure as PNG
	params:
		minTPM=config["RNA_minTPM"],                        # Minimal number of TPM from which we consider a gene expressed
		background=lambda wildcards: wildcards.background,  # Background to use to do the enrichment analysis (either "genome" or "conditions")
		genome=config["genome_version"],                    # Version of the genome ("mm10" or "mm39")
		save_folder=f"{output_tables}"                      # Location where the result tables are saved
	resources:
		cpus_per_task=24,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_sex_DAR_motif_enrich.R"

# rule ATAC_TFBS_motifs_sex_rdm_bg_DAR:
# 	input:
# 		TF_genes=config["TF_genes"],
# 		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
# 		TPM=f"{output_tables}/TPM.csv"
# 	params:
# 		minTPM=config["RNA_minTPM"],
# 		background="genome",
# 		genome=config["genome_version"],
# 		save_folder=f"{output_tables}"
# 	output:
# 		pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_rdm_bg.pdf",
# 		png=f"{output_png}/ATAC_sex_DAR_TF_motifs_rdm_bg.png"
# 	resources:
# 		cpus_per_task=24,
# 		mem_mb=64000
# 	script:
# 		"workflow/scripts/ATAC_sex_DAR_motif_enrich.R"

# rule ATAC_TFBS_motifs_sex_cond_bg_DAR:
# 	input:
# 		TF_genes=config["TF_genes"],
# 		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
# 		TPM=f"{output_tables}/TPM.csv"
# 	params:
# 		minTPM=config["RNA_minTPM"],
# 		background="conditions",
# 		genome=config["genome_version"],
# 		save_folder=f"{output_tables}"
# 	output:
# 		pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_sex_bg.pdf",
# 		png=f"{output_png}/ATAC_sex_DAR_TF_motifs_sex_bg.png"
# 	resources:
# 		cpus_per_task=12,
# 		mem_mb=64000
# 	script:
# 		"workflow/scripts/ATAC_sex_DAR_motif_enrich.R"

rule ATAC_TFBS_motifs_sex_DAR_merged_stages:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		TPM=f"{output_tables}/TPM.csv"
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
		counts=f"{output_tables}/ATAC_raw_counts.csv",
		samplesheet=f"{output_tables}/ATAC_samplesheet.csv",
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

rule ATAC_Get_XY_dynamic_DARs:
	input:
		counts=f"{output_tables}/ATAC_raw_counts.csv",
		samplesheet=f"{output_tables}/ATAC_samplesheet.csv",
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

rule ATAC_XX_TFBS_motifs_stage_cond_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{output_tables}/ATAC_XX_DAR_stage_heatmap_clusters.csv",
		TPM=f"{output_tables}/TPM.csv"
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
		TPM=f"{output_tables}/TPM.csv"
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
		TPM=f"{output_tables}/TPM.csv"
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
		TPM=f"{output_tables}/TPM.csv"
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

###########################################
#                                         #
#                 Plots                   #
#                                         #
###########################################

rule ATAC_Plot_nb_consensus_peak:
	input:
		peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj",
		samplesheet=f"{output_tables}/ATAC_samplesheet.csv",
		norm_data=f"{output_tables}/ATAC_norm_counts.csv"
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
		norm_data=f"{output_tables}/ATAC_norm_counts.csv"
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
# 		peaks=f"{output_tables}/ATAC_norm_counts.csv",
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




rule ATAC_Plot_sex_DAR_histogram:
	input:
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		samplesheet=f"{output_tables}/ATAC_samplesheet.csv"
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


# rule ATAC_Plot_XX_stage_DAR_peak_examples:
# 	input:
# 		genome=f"{genome}",
# 		peaks=f"{output_tables}/ATAC_norm_counts.csv",
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



# rule ATAC_Plot_XY_stage_DAR_peak_examples:
# 	input:
# 		genome=f"{genome}",
# 		peaks=f"{output_tables}/ATAC_norm_counts.csv",
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
		norm_counts=f"{output_tables}/ATAC_norm_counts.csv",
		samplesheet=f"{output_tables}/ATAC_samplesheet.csv"
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
		norm_counts=f"{output_tables}/ATAC_norm_counts.csv",
		samplesheet=f"{output_tables}/ATAC_samplesheet.csv"
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


