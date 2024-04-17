configfile: "config.yaml"

rule_all_input_list = [
	"processed_data/RNA_raw_counts.csv",
	"graphs/PNG/RNA_corr_pca_all_samples.png",
	"graphs/PNG/RNA_corr_pca.png",
	"graphs/PNG/RNA_marker_genes.png",
	"processed_data/RNA_all_SexDEGs.Robj",
	"graphs/PNG/RNA_sex_DEG_histograms.png",
	"graphs/PNG/RNA_sex_DEG_volcano.png"
]

if len(config["RNA_outliers"])<1:
    rule_all_input_list.remove("graphs/PNG/RNA_corr_pca_all_samples.png")

rule all:
	input:
		rule_all_input_list


rule RNA-Get_RNA_matrices:
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
		samplesheet="processed_data/RNA_samplesheet.csv"
	script:
		"scripts/01.RNA_clean_matrices.R"


rule RNA-corr_PCA_all:
	input:
		tpm="processed_data/RNA_TMP_all_samples.csv"
	params:
		corr_method=config["RNA_corr_met"]
	output:
		pdf="graphs/PDF/RNA_corr_pca_all_samples.pdf",
		png="graphs/PNG/RNA_corr_pca_all_samples.png"
	script:
		"scripts/02.RNA_corr_pca.R"


rule RNA-corr_PCA:
	input:
		tpm="processed_data/RNA_TMP.csv"
	params:
		corr_method="spearman"
	output:
		pdf="graphs/PDF/RNA_corr_pca.pdf",
		png="graphs/PNG/RNA_corr_pca.png"
	script:
		"scripts/02.RNA_corr_pca.R"


rule RNA-Plot_marker_genes:
	input:
		tpm="processed_data/RNA_TMP.csv"

	output:
		pdf="graphs/PDF/RNA_marker_genes.pdf",
		png="graphs/PNG/RNA_marker_genes.png"
	script:
		"scripts/03.RNA_plot_marker_genes.R"


rule RNA-Get_sex_DEGs:
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
		"scripts/04.RNA_sex_DEG.R"


rule RNA-Plot_sex_DEG_histogram:
	input:
		sig_DEGs="processed_data/RNA_sig_SexDEGs.Robj",
		samplesheet="processed_data/RNA_samplesheet.csv"
	output:
		pdf="graphs/PDF/RNA_sex_DEG_histograms.pdf",
		png="graphs/PNG/RNA_sex_DEG_histograms.png"
	script:
		"scripts/05.RNA_plot_sex_DEG_hist.R"


rule RNA-Plot_sex_DEG_volcano_GO:
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
		"scripts/06.RNA_plot_sex_volcano_GO.R"


rule RNA-Plot_sex_DEG_double_heatmap:
	input:
		sig_DEGs="processed_data/RNA_sig_SexDEGs.Robj",
		samplesheet="processed_data/RNA_samplesheet.csv"
	params:
		clusters=config["RNA_sex_double_heatmap_clusters"]
	output:
		pdf="graphs/PDF/RNA_sex_DEG_double_heatmap.pdf",
		png="graphs/PNG/RNA_sex_DEG_double_heatmap.png",
		clusters="results/RNA_sex_DEG_double_heatmap_clustering.csv"
		# TFs="results/RNA_sex_DEG_double_heatmap_TFs.csv"
	script:
		"scripts/07.RNA_plot_sex_DEG_double_heatmap.R"
