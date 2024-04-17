rule all:
	input:
		"graphs/PNG/RNA_corr_pca_all_samples.png",
		"graphs/PNG/RNA_corr_pca.png",
		"graphs/PNG/RNA_marker_genes.png",
		"processed_data/RNA_all_SexDEGs.Robj",
		"graphs/PNG/RNA_sex_DEG_histograms.png"
		# "graphs/PNG/RNA_sex_DEG_volcano.png"


rule get_RNA_matrices:
	input:
		counts="data/XX_Enh8-mCherry_XY_SOX9-IRES-GFP_read_count.csv",
		tpm="data/XX_Enh8-mCherry_XY_SOX9-IRES-GFP_TPM.csv"
	output:
		tpm_all="processed_data/RNA_TMP_all_samples.csv",
		tpm="processed_data/RNA_TMP.csv",
		counts="processed_data/RNA_raw_counts.csv",
		samplesheet="processed_data/RNA_samplesheet.csv"
	script:
		"scripts/01.RNA_clean_matrices.R"


rule RNA_corr_PCA_all:
	input:
		tpm="processed_data/RNA_TMP_all_samples.csv"
	params:
		corr_method="spearman"
	output:
		pdf="graphs/PDF/RNA_corr_pca_all_samples.pdf",
		png="graphs/PNG/RNA_corr_pca_all_samples.png"
	script:
		"scripts/02.RNA_corr_pca.R"


rule RNA_corr_PCA:
	input:
		tpm="processed_data/RNA_TMP.csv"
	params:
		corr_method="spearman"
	output:
		pdf="graphs/PDF/RNA_corr_pca.pdf",
		png="graphs/PNG/RNA_corr_pca.png"
	script:
		"scripts/02.RNA_corr_pca.R"


rule Plot_marker_genes:
	input:
		tpm="processed_data/RNA_TMP.csv"

	output:
		pdf="graphs/PDF/RNA_marker_genes.pdf",
		png="graphs/PNG/RNA_marker_genes.png"
	script:
		"scripts/03.RNA_plot_marker_genes.R"


rule Get_sex_DEGs:
	input:
		counts="processed_data/RNA_raw_counts.csv",
		samplesheet="processed_data/RNA_samplesheet.csv",
	params:
		adjpval=[0.01],
		log2FC=[0.5]
	output:
		all_DEGs="processed_data/RNA_all_SexDEGs.Robj",
		sig_DEGs="processed_data/RNA_sig_SexDEGs.Robj"
	script:
		"scripts/04.RNA_sex_DEG.R"


rule Plot_sex_DEG_histogram:
	input:
		sig_DEGs="processed_data/RNA_sig_SexDEGs.Robj",
		samplesheet="processed_data/RNA_samplesheet.csv"
	output:
		pdf="graphs/PDF/RNA_sex_DEG_histograms.pdf",
		png="graphs/PNG/RNA_sex_DEG_histograms.png"
	script:
		"scripts/05.RNA_plot_sex_DEG_hist.R"


# rule Plot_sex_DEG_volcano:
# 	input:
# 		all_DEGs="processed_data/RNA_all_SexDEGs.Robj",
# 		samplesheet="processed_data/RNA_samplesheet.csv"
# 	params:
# 		adjpval=[0.01],
# 		log2FC=[0.5]
# 	output:
# 		pdf="graphs/PDF/RNA_sex_DEG_volcano.pdf",
# 		png="graphs/PNG/RNA_sex_DEG_volcano.png"
# 	script:
# 		"scripts/06.RNA_plot_sex_volcano_GO.R"
