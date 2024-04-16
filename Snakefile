rule all:
	input:
		"graphs/RNA_corr_pca_all_samples.png",
		"graphs/RNA_corr_pca.png"


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
		pdf="graphs/RNA_corr_pca_all_samples.pdf",
		png="graphs/RNA_corr_pca_all_samples.png"
	script:
		"scripts/02.RNA_corr_pca.R"


rule RNA_corr_PCA:
	input:
		tpm="processed_data/RNA_TMP.csv"
	params:
		corr_method="spearman"
	output:
		pdf="graphs/RNA_corr_pca.pdf",
		png="graphs/RNA_corr_pca.png"
	script:
		"scripts/02.RNA_corr_pca.R"