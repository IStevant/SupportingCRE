configfile: "config.yaml"

rule_all_input_list = [
	"graphs/PNG/RNA_corr_pca_all_samples.png",
	"graphs/PNG/RNA_corr_pca.png",
	"graphs/PNG/RNA_marker_genes.png",
	"graphs/PNG/RNA_sex_DEG_histograms.png",
	"graphs/PNG/RNA_sex_DEG_volcano.png",
	"graphs/PNG/RNA_sex_DEG_double_heatmap.png",
	"graphs/PNG/RNA_sex_DEG_upset.png",
	"graphs/PNG/RNA_XX_DEG_stage_heatmap.png",
	"graphs/PNG/RNA_XY_DEG_stage_heatmap.png",
	"graphs/PNG/RNA_sex_stage_common_DEGs.png",
	"graphs/PNG/ATAC_corr_pca_all_samples.png",
	"graphs/PNG/ATAC_all_consensus_peak_annotation.png",
	"graphs/PNG/ATAC_sex_DAR_histograms.png"
]

if len(config["RNA_outliers"])<1:
    rule_all_input_list.remove("graphs/PNG/RNA_corr_pca_all_samples.png")

rule all:
	input:
		rule_all_input_list

rule install_packages:
	script:
		"renv/restore.R"

rule RNA_Get_matrices:
	input:
		counts=config["RNA_counts"],
		tpm=config["RNA_TPM"]
	params:
		minReads=config["RNA_minReads"],
		minTPM=config["RNA_minTPM"],
		RNA_outliers=config["RNA_outliers"]
	output:
		tpm_all="processed_data/RNA_TMP_all_samples.csv",
		tpm="processed_data/RNA_TMP.csv",
		counts="processed_data/RNA_raw_counts.csv",
		norm_counts="processed_data/RNA_norm_counts.csv",
		samplesheet="processed_data/RNA_samplesheet.csv"
	script:
		"scripts/RNA_clean_matrices.R"


rule RNA_corr_PCA_all:
	input:
		norm_data="processed_data/RNA_TMP_all_samples.csv"
	params:
		corr_method=config["RNA_corr_met"]
	output:
		pdf="graphs/PDF/RNA_corr_pca_all_samples.pdf",
		png="graphs/PNG/RNA_corr_pca_all_samples.png"
	script:
		"scripts/Corr_pca.R"


rule RNA_corr_PCA:
	input:
		norm_data="processed_data/RNA_TMP.csv"
	params:
		corr_method="spearman"
	output:
		pdf="graphs/PDF/RNA_corr_pca.pdf",
		png="graphs/PNG/RNA_corr_pca.png"
	script:
		"scripts/Corr_pca.R"


rule RNA_Plot_marker_genes:
	input:
		tpm="processed_data/RNA_TMP.csv"

	output:
		pdf="graphs/PDF/RNA_marker_genes.pdf",
		png="graphs/PNG/RNA_marker_genes.png"
	script:
		"scripts/RNA_plot_marker_genes.R"


rule RNA_Get_sex_DEGs:
	input:
		counts="processed_data/RNA_raw_counts.csv",
		samplesheet="processed_data/RNA_samplesheet.csv",
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"]
	output:
		all_DEGs="processed_data/RNA_all_SexDEGs.Robj",
		sig_DEGs="processed_data/RNA_sig_SexDEGs.Robj"
	script:
		"scripts/RNA_sex_DEG.R"


rule RNA_Plot_sex_DEG_histogram:
	input:
		sig_DEGs="processed_data/RNA_sig_SexDEGs.Robj",
		samplesheet="processed_data/RNA_samplesheet.csv"
	output:
		pdf="graphs/PDF/RNA_sex_DEG_histograms.pdf",
		png="graphs/PNG/RNA_sex_DEG_histograms.png"
	script:
		"scripts/RNA_plot_sex_DEG_hist.R"


rule RNA_Plot_sex_DEG_volcano_GO:
	input:
		all_DEGs="processed_data/RNA_all_SexDEGs.Robj",
		samplesheet="processed_data/RNA_samplesheet.csv"
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"]
	output:
		pdf="graphs/PDF/RNA_sex_DEG_volcano.pdf",
		png="graphs/PNG/RNA_sex_DEG_volcano.png"
	script:
		"scripts/RNA_plot_sex_volcano_GO.R"


rule RNA_Plot_sex_DEG_double_heatmap:
	input:
		sig_DEGs="processed_data/RNA_sig_SexDEGs.Robj",
		norm_counts="processed_data/RNA_norm_counts.csv",
		samplesheet="processed_data/RNA_samplesheet.csv"
	params:
		clusters=config["RNA_sex_double_heatmap_clusters"]
	output:
		pdf="graphs/PDF/RNA_sex_DEG_double_heatmap.pdf",
		png="graphs/PNG/RNA_sex_DEG_double_heatmap.png",
		clusters="results/RNA_sex_DEG_double_heatmap_clustering.csv"
		# TFs="results/RNA_sex_DEG_double_heatmap_TFs.csv"
	script:
		"scripts/RNA_sex_DEG_double_heatmap.R"


rule RNA_Plot_sex_DEG_upset:
	input:
		sig_DEGs="processed_data/RNA_sig_SexDEGs.Robj"
	output:
		pdf="graphs/PDF/RNA_sex_DEG_upset.pdf",
		png="graphs/PNG/RNA_sex_DEG_upset.png"
	script:
		"scripts/RNA_sex_DEG_upset.R"

rule RNA_Get_XX_dynamic_DEGs:
	input:
		counts="processed_data/RNA_raw_counts.csv",
		samplesheet="processed_data/RNA_samplesheet.csv"
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		sex="XX"
	output:
		csv="results/RNA_XX_DEG_stage.csv",
		sig_DEGs="processed_data/RNA_sig_stage_DEGs_XX.Robj"
	script:
		"scripts/RNA_stage_DEG.R"


rule RNA_Get_XY_dynamic_DEGs:
	input:
		counts="processed_data/RNA_raw_counts.csv",
		samplesheet="processed_data/RNA_samplesheet.csv"
	params:
		adjpval=config["RNA_adjpval"],
		log2FC=config["RNA_log2FC"],
		sex="XY"
	output:
		csv="results/RNA_XY_DEG_stage.csv",
		sig_DEGs="processed_data/RNA_sig_stage_DEGs_XY.Robj"
	script:
		"scripts/RNA_stage_DEG.R"


rule RNA_Plot_heatmap_GO_XX:
	input:
		sig_DEGs="processed_data/RNA_sig_stage_DEGs_XX.Robj",
		norm_counts="processed_data/RNA_norm_counts.csv",
		samplesheet="processed_data/RNA_samplesheet.csv"
	params:
		sex="XX",
		clusters=config["RNA_XX_stage_DEG_clusters"]
	output:
		GO="results/RNA_XX_GO_DEG_stage.csv",
		clusters="results/RNA_XX_DEG_stage_heatmap_clusters.csv",
		pdf="graphs/PDF/RNA_XX_DEG_stage_heatmap.pdf",
		png="graphs/PNG/RNA_XX_DEG_stage_heatmap.png"
	script:
		"scripts/RNA_stage_DEG_heatmap.R"


rule RNA_Plot_heatmap_GO_XY:
	input:
		sig_DEGs="processed_data/RNA_sig_stage_DEGs_XY.Robj",
		norm_counts="processed_data/RNA_norm_counts.csv",
		samplesheet="processed_data/RNA_samplesheet.csv"
	params:
		sex="XY",
		clusters=config["RNA_XY_stage_DEG_clusters"]
	output:
		GO="results/RNA_XY_GO_DEG_stage.csv",
		clusters="results/RNA_XY_DEG_stage_heatmap_clusters.csv",
		pdf="graphs/PDF/RNA_XY_DEG_stage_heatmap.pdf",
		png="graphs/PNG/RNA_XY_DEG_stage_heatmap.png"
	script:
		"scripts/RNA_stage_DEG_heatmap.R"


rule RNA_Plot_sex_stage_common_DEGs:
	input:
		sex_DEGs="processed_data/RNA_sig_SexDEGs.Robj",
		XY_stage_DEGs="processed_data/RNA_sig_stage_DEGs_XY.Robj",
		XX_stage_DEGs="processed_data/RNA_sig_stage_DEGs_XX.Robj",
		samplesheet="processed_data/RNA_samplesheet.csv"
	output:
		pdf="graphs/PDF/RNA_sex_stage_common_DEGs.pdf",
		png="graphs/PNG/RNA_sex_stage_common_DEGs.png"
	script:
		"scripts/RNA_overlap_sex_stage_DEG.R"

###################################################################

rule ATAC_Get_matrices:
	input:
		counts=config["ATAC_counts"],
	params:
		minReads=config["ATAC_minReads"],
	output:
		counts="processed_data/ATAC_raw_counts.csv",
		norm_counts="processed_data/ATAC_norm_counts.csv",
		samplesheet="processed_data/ATAC_samplesheet.csv"
	script:
		"scripts/ATAC_clean_matrices.R"


rule ATAC_Plot_consensus_peak_annotation:
	input:
		peak_list="data/ATAC_all_consensus_peaks_2rep_list.Robj"
	params:
		promoter=config["ATAC_promoter_distance"]
	output:
		anno_list="processed_data/ATAC_all_consensus_peak_annotation.Robj",
		pdf="graphs/PDF/ATAC_all_consensus_peak_annotation.pdf",
		png="graphs/PNG/ATAC_all_consensus_peak_annotation.png"
	script:
		"scripts/ATAC_all_peak_annotation.R"


rule ATAC_corr_PCA:
	input:
		norm_data="processed_data/ATAC_norm_counts.csv"
	params:
		corr_method=config["ATAC_corr_met"]
	output:
		pdf="graphs/PDF/ATAC_corr_pca_all_samples.pdf",
		png="graphs/PNG/ATAC_corr_pca_all_samples.png"
	script:
		"scripts/Corr_pca.R"


rule ATAC_Get_sex_DARs:
	input:
		counts="processed_data/ATAC_raw_counts.csv",
		samplesheet="processed_data/ATAC_samplesheet.csv",
	params:
		adjpval=config["ATAC_adjpval"],
		log2FC=config["ATAC_log2FC"]
	output:
		sig_DARs="processed_data/ATAC_sig_SexDARs.Robj"
	script:
		"scripts/ATAC_sex_DAR.R"


rule ATAC_Plot_sex_DAR_histogram:
	input:
		sig_DARs="processed_data/ATAC_sig_SexDARs.Robj",
		samplesheet="processed_data/ATAC_samplesheet.csv"
	output:
		pdf="graphs/PDF/ATAC_sex_DAR_histograms.pdf",
		png="graphs/PNG/ATAC_sex_DAR_histograms.png"
	script:
		"scripts/ATAC_plot_sex_DAR_hist.R"
