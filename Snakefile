configfile: "smk_env/workflow_config.yaml"

# Generate the report
# report: "report/workflow.rst"

# Get file path according to the genome version
input_data = f'{config["path_to_data"]}{config["genome_version"]}'
processed_data = f'{config["path_to_process"]}{config["genome_version"]}'
output_png = f'{config["path_to_graphs"]}{config["genome_version"]}/PNG'
output_pdf = f'{config["path_to_graphs"]}{config["genome_version"]}/PDF'
output_tables = f'{config["path_to_tables"]}{config["genome_version"]}'

if config["genome_version"] == "mm10":
	genome = f'{input_data}/gencode.vM25.annotation.gtf.gz'
	RNA_bigwig_folder = config["RNA_bigwig_folder_mm10"]
	ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm10"]
else :
	genome = f'{input_data}/gencode.vM34.annotation.gtf.gz'
	RNA_bigwig_folder = config["RNA_bigwig_folder_mm39"]
	ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm39"]

# List of output files
rule_all_input_list = [
	f"{output_png}/RNA_corr_pca_all_samples.png",
	f"{output_png}/RNA_corr_pca.png",
	f"{output_png}/RNA_marker_genes.png",
	f"{output_png}/RNA_marker_genes_enrichment.png",
	f"{output_png}/RNA_sex_DEG_histograms.png",
	f"{output_png}/RNA_sex_DEG_volcano.png",
	f"{output_png}/RNA_sex_DEG_double_heatmap.png",
	f"{output_png}/RNA_sex_DEG_upset.png",
	f"{output_png}/RNA_XX_DEG_stage_heatmap.png",
	f"{output_png}/RNA_XY_DEG_stage_heatmap.png",
	f"{output_png}/RNA_sex_stage_common_DEGs.png",
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
	f"{output_pdf}/gene2peak_plots.pdf",
	f"{output_png}/MULTI_TFBS_motifs_peak_XX_genes.png",
	f"{output_png}/MULTI_TFBS_motifs_peak_XY_genes.png",
	f"{processed_data}/plot_example_1.log",
	f"{processed_data}/plot_example_2.log",
	f"{processed_data}/plot_example_3.log",
	f"{processed_data}/plot_example_4.log"
]

# If there is no outliers, do not run the analysis that discard them
if len(config["RNA_outliers"])<1:
	rule_all_input_list.remove(f"{output_png}/RNA_corr_pca_all_samples.png")

rule all:
	input:
		rule_all_input_list

rule install_packages:
	script:
		"renv/restore.R"

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
		tpm=f"{processed_data}/RNA_TPM.csv",
		counts=f"{processed_data}/RNA_raw_counts.csv",
		norm_counts=f"{processed_data}/RNA_norm_counts.csv",
		norm_counts_all=f"{processed_data}/RNA_norm_counts_all.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv",
		size_factors=f"{processed_data}/RNA_size_factors.csv"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/RNA_clean_matrices.R"

# Force the execution if you need to re-normalize the bigwig files
rule RNA_Normalize_bigwig:
	input:
		size_factors=f"{processed_data}/RNA_size_factors.csv"
	params:
		bigwig_folder=f"{RNA_bigwig_folder}",
		new_bigwig_folder=f"{processed_data}/RNA_bigwig"
	resources:
		cpus_per_task=12,
		mem_mb=16000
	script:
		"workflow/scripts/RNA_norm_bigwig.R"

rule RNA_corr_PCA_with outliers:
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
		"workflow/scripts/Corr_pca.R"

rule RNA_corr_PCA:
	input:
		norm_data=f"{processed_data}/RNA_norm_counts.csv"
	params:
		corr_method="spearman"
	output:
		pdf=f"{output_pdf}/RNA_corr_pca.pdf",
		png=f"{output_png}/RNA_corr_pca.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/Corr_pca.R"

rule RNA_Plot_marker_genes:
	input:
		marker_genes=config["marker_genes"],
		whole_gonad=f'{input_data}/{config["whole_gonad_RNAseq"]}',
		tpm=f"{processed_data}/RNA_TPM.csv"
	output:
		pdf1=f"{output_pdf}/RNA_marker_genes.pdf",
		png1=f"{output_png}/RNA_marker_genes.png",
		pdf2=f"{output_pdf}/RNA_marker_genes_enrichment.pdf",
		png2=f"{output_png}/RNA_marker_genes_enrichment.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/RNA_plot_marker_genes.R"

rule RNA_Plot_peak_examples:
	input:
		genome=f"{genome}",
		gene_bed=f"{input_data}/gene_standard.bed",
		gene_list=config["peak_examples"],
	params:
		bw_folder="results/processed_data/mm10/RNA_bigwig",
		save_folder=f"{output_png}"
	output: 
		log= f"{processed_data}/plot_example_4.log"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/RNA_plot_genomic_tracks_examples.R"

rule RNA_Get_sex_DEGs:
	input:
		counts=f"{processed_data}/RNA_raw_counts.csv",
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
		"workflow/scripts/RNA_sex_DEG.R"

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
		"workflow/scripts/RNA_plot_sex_DEG_hist.R"

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
		"workflow/scripts/RNA_plot_sex_volcano_GO.R"

rule RNA_Plot_sex_DEG_double_heatmap:
	input:
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"],
		sig_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj",
		norm_counts=f"{processed_data}/RNA_norm_counts.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	params:
		clusters=config["RNA_sex_double_heatmap_clusters"]
	output:
		pdf=f"{output_pdf}/RNA_sex_DEG_double_heatmap.pdf",
		png=f"{output_png}/RNA_sex_DEG_double_heatmap.png",
		clusters=f"{output_tables}/RNA_sex_DEG_double_heatmap_clustering.csv"
		# TFs=f"{output_tables}/RNA_sex_DEG_double_heatmap_TFs.csv"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/RNA_sex_DEG_double_heatmap.R"

rule RNA_Plot_sex_DEG_upset:
	input:
		sig_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj"
	output:
		pdf=f"{output_pdf}/RNA_sex_DEG_upset.pdf",
		png=f"{output_png}/RNA_sex_DEG_upset.png"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/RNA_sex_DEG_upset.R"

rule RNA_Get_XX_dynamic_DEGs:
	input:
		counts=f"{processed_data}/RNA_raw_counts.csv",
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
		"workflow/scripts/RNA_stage_DEG.R"

rule RNA_Get_XY_dynamic_DEGs:
	input:
		counts=f"{processed_data}/RNA_raw_counts.csv",
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
		"workflow/scripts/RNA_stage_DEG.R"

rule RNA_Plot_heatmap_GO_XX:
	input:
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"],
		sig_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XX.Robj",
		norm_counts=f"{processed_data}/RNA_norm_counts.csv",
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
		"workflow/scripts/RNA_stage_DEG_heatmap.R"

rule RNA_Plot_heatmap_GO_XY:
	input:
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"],
		sig_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XY.Robj",
		norm_counts=f"{processed_data}/RNA_norm_counts.csv",
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
		"workflow/scripts/RNA_stage_DEG_heatmap.R"

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
		"workflow/scripts/RNA_overlap_sex_stage_DEG.R"

################################################################################################
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
		"workflow/scripts/ATAC_clean_matrices.R"

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
		"workflow/scripts/MULTI_norm_bigwig.R"

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
		"workflow/scripts/ATAC_peak_distribution_per_sex.R"

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
		"workflow/scripts/ATAC_peak_annotation_per_sex.R"

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
		"workflow/scripts/Corr_pca.R"

rule ATAC_Plot_peak_examples:
	input:
		genome=f"{genome}",
		gene_bed=f"{input_data}/gene_standard.bed",
		peaks=f"{processed_data}/ATAC_norm_counts.csv",
		linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
		gene_list=config["peak_examples"],
	params:
		bw_folder="results/processed_data/mm10/ATAC_bigwig",
		save_folder=f"{output_png}"
	output: 
		log= f"{processed_data}/plot_example_1.log"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/MULTI_plot_genomic_tracks_examples.R"


rule ATAC_Get_sex_DARs:
	input:
		counts=f"{processed_data}/ATAC_raw_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		save_folder=f"{output_tables}"
	output:
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		sig_DARs_GR=f"{processed_data}/ATAC_sig_SexDARs_GR.Robj"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_sex_DAR.R"

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
		"workflow/scripts/ATAC_plot_sex_DAR_hist.R"

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
		"workflow/scripts/ATAC_peak_annotation_per_sex.R"

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
		"workflow/scripts/ATAC_sex_DAR_upset.R"

rule ATAC_TFBS_motifs_sex_rdm_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		TPM=f"{processed_data}/RNA_TPM.csv"
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
		TPM=f"{processed_data}/RNA_TPM.csv"
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
		TPM=f"{processed_data}/RNA_TPM.csv"
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
		peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj"
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		sex="XX"
	output:
		csv=f"{output_tables}/ATAC_XX_DEG_stage.csv",
		sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XX.Robj"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_stage_DAR.R"

rule ATAC_Plot_XX_stage_DAR_peak_examples:
	input:
		genome=f"{genome}",
		peaks=f"{processed_data}/ATAC_norm_counts.csv",
		peak_list=config["DAR_peak_examples_XX"],
	params:
		bw_folder="results/processed_data/mm10/ATAC_bigwig",
		save_folder=f"{output_png}",
		sex="XX"
	output: 
		log= f"{processed_data}/plot_example_2.log"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_plot_DAR_peak_examples.R"


rule ATAC_Get_XY_dynamic_DARs:
	input:
		counts=f"{processed_data}/ATAC_raw_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
		peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj"
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		sex="XY"
	output:
		csv=f"{output_tables}/ATAC_XY_DEG_stage.csv",
		sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XY.Robj"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_stage_DAR.R"

rule ATAC_Plot_XY_stage_DAR_peak_examples:
	input:
		genome=f"{genome}",
		peaks=f"{processed_data}/ATAC_norm_counts.csv",
		peak_list=config["DAR_peak_examples_XY"],
	params:
		bw_folder="results/processed_data/mm10/ATAC_bigwig",
		save_folder=f"{output_png}",
		sex="XY"
	output: 
		log= f"{processed_data}/plot_example_3.log"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/ATAC_plot_DAR_peak_examples.R"


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
		TPM=f"{processed_data}/RNA_TPM.csv"
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
		TPM=f"{processed_data}/RNA_TPM.csv"
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
################################################################################################
rule MULTI_Get_all_gene_peak_correlation:
	input:
		RNA_samplesheet=f"{processed_data}/RNA_samplesheet.csv",
		ATAC_samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
		RNA_norm_counts=f"{processed_data}/RNA_norm_counts.csv",
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
		pdf=f"{output_pdf}/gene2peak_plots.pdf"
	resources:
		cpus_per_task=12,
		mem_mb=64000
	script:
		"workflow/scripts/MULTI_plot_gene2peak.R"

rule MULTI_TFBS_motifs_peak_XX_genes:
	input:
		linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
		sex_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj",
		TPM=f"{processed_data}/RNA_TPM.csv"
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
		TPM=f"{processed_data}/RNA_TPM.csv"
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
