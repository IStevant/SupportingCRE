configfile: "smk_env/config.yaml"

rule_all_input_list = [
	"results/graphs/PNG/RNA_corr_pca_all_samples.png",
	"results/graphs/PNG/RNA_corr_pca.png",
	"results/graphs/PNG/RNA_marker_genes.png",
	"results/graphs/PNG/RNA_sex_DEG_histograms.png",
	"results/graphs/PNG/RNA_sex_DEG_volcano.png",
	"results/graphs/PNG/RNA_sex_DEG_double_heatmap.png",
	"results/graphs/PNG/RNA_sex_DEG_upset.png",
	"results/graphs/PNG/RNA_XX_DEG_stage_heatmap.png",
	"results/graphs/PNG/RNA_XY_DEG_stage_heatmap.png",
	"results/graphs/PNG/RNA_sex_stage_common_DEGs.png",
	"results/graphs/PNG/ATAC_corr_pca_all_samples.png",
	"results/graphs/PNG/ATAC_all_consensus_peak_annotation.png",
	"results/graphs/PNG/ATAC_sex_DAR_histograms.png",
	"results/graphs/PNG/ATAC_sig_sex_DARs_annotation.png",
	"results/graphs/PNG/ATAC_sex_DAR_upset.png",
	"results/graphs/PDF/ATAC_sex_DAR_TF_motifs_rdm_bg.pdf",
	"results/graphs/PDF/ATAC_sex_DAR_TF_motifs_sex_bg.pdf",
	"results/graphs/PNG/ATAC_XX_DAR_stage_heatmap.png",
	"results/graphs/PNG/ATAC_XY_DAR_stage_heatmap.png",
	"results/tables/XX_sig_gene2peak_linkage.csv",
	"results/tables/XY_sig_gene2peak_linkage.csv"
]

if len(config["RNA_outliers"])<1:
	rule_all_input_list.remove("results/graphs/PNG/RNA_corr_pca_all_samples.png")

rule all:
	input:
		rule_all_input_list

rule install_packages:
	script:
		"renv/restore.R"


rule RNA_Get_matrices:
	input:
		counts=config["RNA_counts"],
		tpm=config["RNA_TPM"],
		protein_genes=config["protein_genes"]
	params:
		minReads=config["RNA_minReads"],
		minTPM=config["RNA_minTPM"],
		RNA_outliers=config["RNA_outliers"]
	output:
		tpm_all="results/processed_data/RNA_TMP_all_samples.csv",
		tpm="results/processed_data/RNA_TMP.csv",
		counts="results/processed_data/RNA_raw_counts.csv",
		norm_counts="results/processed_data/RNA_norm_counts.csv",
		samplesheet="results/processed_data/RNA_samplesheet.csv"
	script:
		"workflow/scripts/RNA_clean_matrices.R"

rule RNA_corr_PCA_with outliers:
	input:
		norm_data="results/processed_data/RNA_TMP_all_samples.csv"
	params:
		corr_method=config["RNA_corr_met"]
	output:
		pdf="results/graphs/PDF/RNA_corr_pca_all_samples.pdf",
		png="results/graphs/PNG/RNA_corr_pca_all_samples.png"
	script:
		"workflow/scripts/Corr_pca.R"

rule RNA_corr_PCA:
	input:
		norm_data="results/processed_data/RNA_TMP.csv"
	params:
		corr_method="spearman"
	output:
		pdf="results/graphs/PDF/RNA_corr_pca.pdf",
		png="results/graphs/PNG/RNA_corr_pca.png"
	script:
		"workflow/scripts/Corr_pca.R"

rule RNA_Plot_marker_genes:
	input:
		marker_genes=config["marker_genes"],
		whole_gonad=config["whole_gonad_RNAseq"],
		tpm="results/processed_data/RNA_TMP.csv"

	output:
		pdf="results/graphs/PDF/RNA_marker_genes.pdf",
		png="results/graphs/PNG/RNA_marker_genes.png"
	script:
		"workflow/scripts/RNA_plot_marker_genes.R"

rule RNA_Get_sex_DEGs:
	input:
		counts="results/processed_data/RNA_raw_counts.csv",
		samplesheet="results/processed_data/RNA_samplesheet.csv",
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"]
	output:
		all_DEGs="results/processed_data/RNA_all_SexDEGs.Robj",
		sig_DEGs="results/processed_data/RNA_sig_SexDEGs.Robj"
	script:
		"workflow/scripts/RNA_sex_DEG.R"

rule RNA_Plot_sex_DEG_histogram:
	input:
		sig_DEGs="results/processed_data/RNA_sig_SexDEGs.Robj",
		samplesheet="results/processed_data/RNA_samplesheet.csv"
	output:
		pdf="results/graphs/PDF/RNA_sex_DEG_histograms.pdf",
		png="results/graphs/PNG/RNA_sex_DEG_histograms.png"
	script:
		"workflow/scripts/RNA_plot_sex_DEG_hist.R"

rule RNA_Plot_sex_DEG_volcano_GO:
	input:
		all_DEGs="results/processed_data/RNA_all_SexDEGs.Robj",
		samplesheet="results/processed_data/RNA_samplesheet.csv"
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"]
	output:
		pdf="results/graphs/PDF/RNA_sex_DEG_volcano.pdf",
		png="results/graphs/PNG/RNA_sex_DEG_volcano.png"
	script:
		"workflow/scripts/RNA_plot_sex_volcano_GO.R"

rule RNA_Plot_sex_DEG_double_heatmap:
	input:
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"],
		sig_DEGs="results/processed_data/RNA_sig_SexDEGs.Robj",
		norm_counts="results/processed_data/RNA_norm_counts.csv",
		samplesheet="results/processed_data/RNA_samplesheet.csv"
	params:
		clusters=config["RNA_sex_double_heatmap_clusters"]
	output:
		pdf="results/graphs/PDF/RNA_sex_DEG_double_heatmap.pdf",
		png="results/graphs/PNG/RNA_sex_DEG_double_heatmap.png",
		clusters="results/tables/RNA_sex_DEG_double_heatmap_clustering.csv"
		# TFs="results/tables/RNA_sex_DEG_double_heatmap_TFs.csv"
	script:
		"workflow/scripts/RNA_sex_DEG_double_heatmap.R"

rule RNA_Plot_sex_DEG_upset:
	input:
		sig_DEGs="results/processed_data/RNA_sig_SexDEGs.Robj"
	output:
		pdf="results/graphs/PDF/RNA_sex_DEG_upset.pdf",
		png="results/graphs/PNG/RNA_sex_DEG_upset.png"
	script:
		"workflow/scripts/RNA_sex_DEG_upset.R"

rule RNA_Get_XX_dynamic_DEGs:
	input:
		counts="results/processed_data/RNA_raw_counts.csv",
		samplesheet="results/processed_data/RNA_samplesheet.csv"
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		sex="XX"
	output:
		csv="results/tables/RNA_XX_DEG_stage.csv",
		sig_DEGs="results/processed_data/RNA_sig_stage_DEGs_XX.Robj"
	script:
		"workflow/scripts/RNA_stage_DEG.R"

rule RNA_Get_XY_dynamic_DEGs:
	input:
		counts="results/processed_data/RNA_raw_counts.csv",
		samplesheet="results/processed_data/RNA_samplesheet.csv"
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		sex="XY"
	output:
		csv="results/tables/RNA_XY_DEG_stage.csv",
		sig_DEGs="results/processed_data/RNA_sig_stage_DEGs_XY.Robj"
	script:
		"workflow/scripts/RNA_stage_DEG.R"

rule RNA_Plot_heatmap_GO_XX:
	input:
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"],
		sig_DEGs="results/processed_data/RNA_sig_stage_DEGs_XX.Robj",
		norm_counts="results/processed_data/RNA_norm_counts.csv",
		samplesheet="results/processed_data/RNA_samplesheet.csv"
	params:
		sex="XX",
		clusters=config["RNA_XX_stage_DEG_clusters"]
	output:
		GO="results/tables/RNA_XX_GO_DEG_stage.csv",
		clusters="results/tables/RNA_XX_DEG_stage_heatmap_clusters.csv",
		pdf="results/graphs/PDF/RNA_XX_DEG_stage_heatmap.pdf",
		png="results/graphs/PNG/RNA_XX_DEG_stage_heatmap.png"
	script:
		"workflow/scripts/RNA_stage_DEG_heatmap.R"

rule RNA_Plot_heatmap_GO_XY:
	input:
		TF_genes=config["TF_genes"],
		TF_pheno=config["TF_pheno"],
		sig_DEGs="results/processed_data/RNA_sig_stage_DEGs_XY.Robj",
		norm_counts="results/processed_data/RNA_norm_counts.csv",
		samplesheet="results/processed_data/RNA_samplesheet.csv"
	params:
		sex="XY",
		clusters=config["RNA_XY_stage_DEG_clusters"]
	output:
		GO="results/tables/RNA_XY_GO_DEG_stage.csv",
		clusters="results/tables/RNA_XY_DEG_stage_heatmap_clusters.csv",
		pdf="results/graphs/PDF/RNA_XY_DEG_stage_heatmap.pdf",
		png="results/graphs/PNG/RNA_XY_DEG_stage_heatmap.png"
	script:
		"workflow/scripts/RNA_stage_DEG_heatmap.R"

rule RNA_Plot_sex_stage_common_DEGs:
	input:
		sex_DEGs="results/processed_data/RNA_sig_SexDEGs.Robj",
		XY_stage_DEGs="results/processed_data/RNA_sig_stage_DEGs_XY.Robj",
		XX_stage_DEGs="results/processed_data/RNA_sig_stage_DEGs_XX.Robj",
		samplesheet="results/processed_data/RNA_samplesheet.csv"
	output:
		pdf="results/graphs/PDF/RNA_sex_stage_common_DEGs.pdf",
		png="results/graphs/PNG/RNA_sex_stage_common_DEGs.png"
	script:
		"workflow/scripts/RNA_overlap_sex_stage_DEG.R"

################################################################################################
rule ATAC_Get_matrices:
	input:
		counts=config["ATAC_counts"],
	params:
		minReads=config["ATAC_minReads"],
	output:
		counts="results/processed_data/ATAC_raw_counts.csv",
		norm_counts="results/processed_data/ATAC_norm_counts.csv",
		samplesheet="results/processed_data/ATAC_samplesheet.csv"
	script:
		"workflow/scripts/ATAC_clean_matrices.R"

rule ATAC_Plot_consensus_peak_annotation:
	input:
		genome=config["genome"],
		peak_list="workflow/data/ATAC_all_consensus_peaks_2rep_list.Robj"
	params:
		promoter=config["ATAC_promoter_distance"]
	output:
		anno_list="results/processed_data/ATAC_all_consensus_peak_annotation.Robj",
		pdf="results/graphs/PDF/ATAC_all_consensus_peak_annotation.pdf",
		png="results/graphs/PNG/ATAC_all_consensus_peak_annotation.png"
	script:
		"workflow/scripts/ATAC_peak_annotation_per_sex.R"

rule ATAC_corr_PCA:
	input:
		norm_data="results/processed_data/ATAC_norm_counts.csv"
	params:
		corr_method=config["ATAC_corr_met"]
	output:
		pdf="results/graphs/PDF/ATAC_corr_pca_all_samples.pdf",
		png="results/graphs/PNG/ATAC_corr_pca_all_samples.png"
	script:
		"workflow/scripts/Corr_pca.R"

rule ATAC_Get_sex_DARs:
	input:
		counts="results/processed_data/ATAC_raw_counts.csv",
		samplesheet="results/processed_data/ATAC_samplesheet.csv",
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"]
	output:
		sig_DARs="results/processed_data/ATAC_sig_SexDARs.Robj",
		sig_DARs_GR="results/processed_data/ATAC_sig_SexDARs_GR.Robj"
	script:
		"workflow/scripts/ATAC_sex_DAR.R"

rule ATAC_Plot_sex_DAR_histogram:
	input:
		sig_DARs="results/processed_data/ATAC_sig_SexDARs.Robj",
		samplesheet="results/processed_data/ATAC_samplesheet.csv"
	output:
		pdf="results/graphs/PDF/ATAC_sex_DAR_histograms.pdf",
		png="results/graphs/PNG/ATAC_sex_DAR_histograms.png"
	script:
		"workflow/scripts/ATAC_plot_sex_DAR_hist.R"

rule ATAC_Plot_sex_DAR_peak_annotation:
	input:
		genome=config["genome"],
		peak_list="results/processed_data/ATAC_sig_SexDARs_GR.Robj"
	params:
		promoter=config["ATAC_promoter_distance"]
	output:
		anno_list="results/processed_data/ATAC_sig_SexDARs_annotation.Robj",
		pdf="results/graphs/PDF/ATAC_sig_sex_DARs_annotation.pdf",
		png="results/graphs/PNG/ATAC_sig_sex_DARs_annotation.png"
	script:
		"workflow/scripts/ATAC_peak_annotation_per_sex.R"

rule ATAC_Plot_sex_DAR_upset:
	input:
		sig_DARs="results/processed_data/ATAC_sig_SexDARs.Robj"
	output:
		pdf="results/graphs/PDF/ATAC_sex_DAR_upset.pdf",
		png="results/graphs/PNG/ATAC_sex_DAR_upset.png"
	script:
		"workflow/scripts/ATAC_sex_DAR_upset.R"

rule ATAC_TFBS_motifs_sex_rdm_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs="results/processed_data/ATAC_sig_SexDARs.Robj",
		TPM="results/processed_data/RNA_TMP.csv"
	params:
		minTPM=config["RNA_minTPM"],
		background="genome",
		logos="FALSE",
		nbTFs=20
	output:
		pdf="results/graphs/PDF/ATAC_sex_DAR_TF_motifs_rdm_bg.pdf",
		png="results/graphs/PNG/ATAC_sex_DAR_TF_motifs_rdm_bg.png"
	script:
		"workflow/scripts/ATAC_sex_DAR_motif_enrich.R"

rule ATAC_TFBS_motifs_sex_cond_bg_DAR:
	input:
		TF_genes=config["TF_genes"],
		sig_DARs="results/processed_data/ATAC_sig_SexDARs.Robj",
		TPM="results/processed_data/RNA_TMP.csv"
	params:
		minTPM=config["RNA_minTPM"],
		background="conditions",
		logos="TRUE",
		nbTFs=10
	output:
		pdf="results/graphs/PDF/ATAC_sex_DAR_TF_motifs_sex_bg.pdf",
		png="results/graphs/PNG/ATAC_sex_DAR_TF_motifs_sex_bg.png"
	script:
		"workflow/scripts/ATAC_sex_DAR_motif_enrich.R"

rule ATAC_Get_XX_dynamic_DARs:
	input:
		counts="results/processed_data/ATAC_raw_counts.csv",
		samplesheet="results/processed_data/ATAC_samplesheet.csv",
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		sex="XX"
	output:
		csv="results/tables/ATAC_XX_DEG_stage.csv",
		sig_DARs="results/processed_data/ATAC_sig_stage_DARs_XX.Robj"
	script:
		"workflow/scripts/ATAC_stage_DAR.R"

rule ATAC_Get_XY_dynamic_DARs:
	input:
		counts="results/processed_data/ATAC_raw_counts.csv",
		samplesheet="results/processed_data/ATAC_samplesheet.csv",
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"],
		sex="XY"
	output:
		csv="results/tables/ATAC_XY_DEG_stage.csv",
		sig_DARs="results/processed_data/ATAC_sig_stage_DARs_XY.Robj"
	script:
		"workflow/scripts/ATAC_stage_DAR.R"

rule ATAC_Plot_heatmap_dyn_DARs_XX:
	input:
		sig_DARs="results/processed_data/ATAC_sig_stage_DARs_XX.Robj",
		norm_counts="results/processed_data/ATAC_norm_counts.csv",
		samplesheet="results/processed_data/ATAC_samplesheet.csv"
	params:
		sex="XX",
		clusters=config["ATAC_XX_stage_DAR_clusters"]
	output:
		clusters="results/tables/ATAC_XX_DAR_stage_heatmap_clusters.csv",
		pdf="results/graphs/PDF/ATAC_XX_DAR_stage_heatmap.pdf",
		png="results/graphs/PNG/ATAC_XX_DAR_stage_heatmap.png"
	script:
		"workflow/scripts/ATAC_stage_DAR_heatmap.R"

rule ATAC_Plot_heatmap_dyn_DARs_XY:
	input:
		sig_DARs="results/processed_data/ATAC_sig_stage_DARs_XY.Robj",
		norm_counts="results/processed_data/ATAC_norm_counts.csv",
		samplesheet="results/processed_data/ATAC_samplesheet.csv"
	params:
		sex="XY",
		clusters=config["ATAC_XY_stage_DAR_clusters"]
	output:
		clusters="results/tables/ATAC_XY_DAR_stage_heatmap_clusters.csv",
		pdf="results/graphs/PDF/ATAC_XY_DAR_stage_heatmap.pdf",
		png="results/graphs/PNG/ATAC_XY_DAR_stage_heatmap.png"
	script:
		"workflow/scripts/ATAC_stage_DAR_heatmap.R"

################################################################################################
rule MULTI_Get_XX_gene_peak_correlation:
	input:
		RNA_samplesheet="results/processed_data/RNA_samplesheet.csv",
		ATAC_samplesheet="results/processed_data/ATAC_samplesheet.csv",
		RNA_norm_counts="results/processed_data/RNA_norm_counts.csv",
		ATAC_norm_counts="results/processed_data/ATAC_norm_counts.csv"
	params:
		sex="XX",
		distance=config["MULTI_peak_gene_distance"],
		min_cor=config["MULTI_peak_gene_min_cor"],
		FDR=config["MULTI_peak_gene_FDR"]
	output:
		linkage="results/tables/XX_sig_gene2peak_linkage.csv",
	script:
		"workflow/scripts/MULTI_gene2peak.R"

rule MULTI_Get_XY_gene_peak_correlation:
	input:
		RNA_samplesheet="results/processed_data/RNA_samplesheet.csv",
		ATAC_samplesheet="results/processed_data/ATAC_samplesheet.csv",
		RNA_norm_counts="results/processed_data/RNA_norm_counts.csv",
		ATAC_norm_counts="results/processed_data/ATAC_norm_counts.csv"
	params:
		sex="XY",
		distance=config["MULTI_peak_gene_distance"],
		min_cor=config["MULTI_peak_gene_min_cor"],
		FDR=config["MULTI_peak_gene_FDR"]
	output:
		linkage="results/tables/XY_sig_gene2peak_linkage.csv",
	script:
		"workflow/scripts/MULTI_gene2peak.R"
