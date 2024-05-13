configfile: "smk_env/config.yaml"

# Get file path according to the genome version
input_data = f'{config["path_to_data"]}{config["genome_version"]}'
processed_data = f'{config["path_to_process"]}{config["genome_version"]}'
output_png = f'{config["path_to_graphs"]}{config["genome_version"]}/PNG'
output_pdf = f'{config["path_to_graphs"]}{config["genome_version"]}/PDF'
output_tables = f'{config["path_to_tables"]}{config["genome_version"]}'

if config["genome_version"] == "mm10":
	genome = f'{input_data}/iGenome_mm10_ucsc_genes.gtf.gz'
else :
	genome = f'{input_data}/gencode.vM34.annotation.gtf.gz'

# List of output files
rule_all_input_list = [
	f"{output_png}/RNA_corr_pca_all_samples.png",
	f"{output_png}/RNA_corr_pca.png",
	f"{output_png}/RNA_marker_genes.png",
	f"{output_png}/RNA_sex_DEG_histograms.png",
	f"{output_png}/RNA_sex_DEG_volcano.png",
	f"{output_png}/RNA_sex_DEG_double_heatmap.png",
	f"{output_png}/RNA_sex_DEG_upset.png",
	f"{output_png}/RNA_XX_DEG_stage_heatmap.png",
	f"{output_png}/RNA_XY_DEG_stage_heatmap.png",
	f"{output_png}/RNA_sex_stage_common_DEGs.png",
	f"{output_png}/ATAC_corr_pca_all_samples.png",
	f"{output_png}/ATAC_all_consensus_peak_annotation.png",
	f"{output_png}/ATAC_sex_DAR_histograms.png",
	f"{output_png}/ATAC_sig_sex_DARs_annotation.png",
	f"{output_png}/ATAC_sex_DAR_upset.png",
	f"{output_pdf}/ATAC_sex_DAR_TF_motifs_rdm_bg.pdf",
	f"{output_pdf}/ATAC_sex_DAR_TF_motifs_sex_bg.pdf",
	f"{output_png}/ATAC_XX_DAR_stage_heatmap.png",
	f"{output_png}/ATAC_XY_DAR_stage_heatmap.png",
	# f"{output_tables}/XX_sig_gene2peak_linkage.csv",
	# f"{output_tables}/XY_sig_gene2peak_linkage.csv",
	f"{output_tables}/all_sig_gene2peak_linkage.csv"
]

# If there is no outliers, do not run the analysis that discard the
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
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	script:
		"workflow/scripts/RNA_clean_matrices.R"

rule RNA_corr_PCA_with outliers:
	input:
		norm_data=f"{processed_data}/RNA_norm_counts_all.csv"
	params:
		corr_method=config["RNA_corr_met"]
	output:
		pdf=f"{output_pdf}/RNA_corr_pca_all_samples.pdf",
		png=f"{output_png}/RNA_corr_pca_all_samples.png"
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
	script:
		"workflow/scripts/Corr_pca.R"

rule RNA_Plot_marker_genes:
	input:
		marker_genes=config["marker_genes"],
		whole_gonad=f'{input_data}/{config["whole_gonad_RNAseq"]}',
		tpm=f"{processed_data}/RNA_TPM.csv"
	output:
		pdf=f"{output_pdf}/RNA_marker_genes.pdf",
		png=f"{output_png}/RNA_marker_genes.png"
	script:
		"workflow/scripts/RNA_plot_marker_genes.R"

rule RNA_Get_sex_DEGs:
	input:
		counts=f"{processed_data}/RNA_raw_counts.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv",
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		save_folder=f"{output_tables}"
	output:
		all_DEGs=f"{processed_data}/RNA_all_SexDEGs.Robj",
		sig_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj"
	script:
		"workflow/scripts/RNA_sex_DEG.R"

rule RNA_Plot_sex_DEG_histogram:
	input:
		sig_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	output:
		pdf=f"{output_pdf}/RNA_sex_DEG_histograms.pdf",
		png=f"{output_png}/RNA_sex_DEG_histograms.png"
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
	script:
		"workflow/scripts/RNA_sex_DEG_double_heatmap.R"

rule RNA_Plot_sex_DEG_upset:
	input:
		sig_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj"
	output:
		pdf=f"{output_pdf}/RNA_sex_DEG_upset.pdf",
		png=f"{output_png}/RNA_sex_DEG_upset.png"
	script:
		"workflow/scripts/RNA_sex_DEG_upset.R"

rule RNA_Get_XX_dynamic_DEGs:
	input:
		counts=f"{processed_data}/RNA_raw_counts.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		sex="XX"
	output:
		csv=f"{output_tables}/RNA_XX_DEG_stage.csv",
		sig_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XX.Robj"
	script:
		"workflow/scripts/RNA_stage_DEG.R"

rule RNA_Get_XY_dynamic_DEGs:
	input:
		counts=f"{processed_data}/RNA_raw_counts.csv",
		samplesheet=f"{processed_data}/RNA_samplesheet.csv"
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		sex="XY"
	output:
		csv=f"{output_tables}/RNA_XY_DEG_stage.csv",
		sig_DEGs=f"{processed_data}/RNA_sig_stage_DEGs_XY.Robj"
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
		clusters=f"{output_tables}/RNA_XX_DEG_stage_heatmap_clusters.csv",
		pdf=f"{output_pdf}/RNA_XX_DEG_stage_heatmap.pdf",
		png=f"{output_png}/RNA_XX_DEG_stage_heatmap.png"
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
		clusters=f"{output_tables}/RNA_XY_DEG_stage_heatmap_clusters.csv",
		pdf=f"{output_pdf}/RNA_XY_DEG_stage_heatmap.pdf",
		png=f"{output_png}/RNA_XY_DEG_stage_heatmap.png"
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
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv"
	script:
		"workflow/scripts/ATAC_clean_matrices.R"

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
	script:
		"workflow/scripts/Corr_pca.R"

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
	script:
		"workflow/scripts/ATAC_sex_DAR.R"

rule ATAC_Plot_sex_DAR_histogram:
	input:
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv"
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_histograms.pdf",
		png=f"{output_png}/ATAC_sex_DAR_histograms.png"
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
	script:
		"workflow/scripts/ATAC_peak_annotation_per_sex.R"

rule ATAC_Plot_sex_DAR_upset:
	input:
		sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj"
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_upset.pdf",
		png=f"{output_png}/ATAC_sex_DAR_upset.png"
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
		logos="TRUE",
		nbTFs=40,
		genome=config["genome_version"],
		save_folder=f"{output_tables}"
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_rdm_bg.pdf",
		png=f"{output_png}/ATAC_sex_DAR_TF_motifs_rdm_bg.png"
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
		logos="TRUE",
		nbTFs=40,
		genome=config["genome_version"],
		save_folder=f"{output_tables}"
	output:
		pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_sex_bg.pdf",
		png=f"{output_png}/ATAC_sex_DAR_TF_motifs_sex_bg.png"
	script:
		"workflow/scripts/ATAC_sex_DAR_motif_enrich.R"

rule ATAC_Get_XX_dynamic_DARs:
	input:
		counts=f"{processed_data}/ATAC_raw_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		sex="XX"
	output:
		csv=f"{output_tables}/ATAC_XX_DEG_stage.csv",
		sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XX.Robj"
	script:
		"workflow/scripts/ATAC_stage_DAR.R"

rule ATAC_Get_XY_dynamic_DARs:
	input:
		counts=f"{processed_data}/ATAC_raw_counts.csv",
		samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		sex="XY"
	output:
		csv=f"{output_tables}/ATAC_XY_DEG_stage.csv",
		sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XY.Robj"
	script:
		"workflow/scripts/ATAC_stage_DAR.R"

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
	script:
		"workflow/scripts/ATAC_stage_DAR_heatmap.R"

################################################################################################
# rule MULTI_Get_XX_gene_peak_correlation:
# 	input:
# 		RNA_samplesheet=f"{processed_data}/RNA_samplesheet.csv",
# 		ATAC_samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
# 		RNA_norm_counts=f"{processed_data}/RNA_norm_counts.csv",
# 		ATAC_norm_counts=f"{processed_data}/ATAC_norm_counts.csv",
# 		chrom_size=f"{input_data}/chrom.size",
# 		genes=f"{input_data}/gene_standard.bed",
# 		gtf=f"{genome}"
# 	params:
# 		sex="XX",
# 		distance=config["MULTI_peak_gene_distance"],
# 		min_cor=config["MULTI_peak_gene_min_cor"],
# 		FDR=config["MULTI_peak_gene_FDR"]
# 	output:
# 		linkage=f"{output_tables}/XX_sig_gene2peak_linkage.csv",
# 	script:
# 		"workflow/scripts/MULTI_gene2peak.R"

# rule MULTI_Get_XY_gene_peak_correlation:
# 	input:
# 		RNA_samplesheet=f"{processed_data}/RNA_samplesheet.csv",
# 		ATAC_samplesheet=f"{processed_data}/ATAC_samplesheet.csv",
# 		RNA_norm_counts=f"{processed_data}/RNA_norm_counts.csv",
# 		ATAC_norm_counts=f"{processed_data}/ATAC_norm_counts.csv",
# 		chrom_size=f"{input_data}/chrom.size",
# 		genes=f"{input_data}/gene_standard.bed",
# 		gtf=f"{genome}"
# 	params:
# 		sex="XY",
# 		distance=config["MULTI_peak_gene_distance"],
# 		min_cor=config["MULTI_peak_gene_min_cor"],
# 		FDR=config["MULTI_peak_gene_FDR"]
# 	output:
# 		linkage=f"{output_tables}/XY_sig_gene2peak_linkage.csv",
# 	script:
# 		"workflow/scripts/MULTI_gene2peak.R"

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
	script:
		"workflow/scripts/MULTI_gene2peak.R"
