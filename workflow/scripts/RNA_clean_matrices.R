source(".Rprofile")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

count_file <- snakemake@input[['counts']]
TPM_file <- snakemake@input[['tpm']]
protein_genes <- snakemake@input[['protein_genes']]
minReads <- snakemake@params[['minReads']]
minTPM <- snakemake@params[['minTPM']]
outlierSamples <- snakemake@params[['RNA_outliers']]

#################################################################################################################################

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Generate the read count matrix
#' @param csv_file Path to the read count matrix.
#' @param genes Keep all the genes or only the protein coding genes. Values can be "all" or "protein".
#' @return Return a dataframe.
get_gene_matrix <- function(csv_file, genes){
	# load file
	raw_counts <- read.csv(file=csv_file, row.names=1)
	# Transform values as integers for DESeq2 that does not support floats
	raw_counts <- round(raw_counts, digits = 0)
	# Select protein coding genes
	if (genes=="protein"){
		protein_coding_genes <- read.csv(file=protein_genes)
		raw_counts <- raw_counts[rownames(raw_counts) %in% protein_coding_genes$external_gene_name,]
		rm(protein_coding_genes)
	}
	raw_counts <- raw_counts[, order(names(raw_counts))]
	return(raw_counts)
}

#' Generate the tpm matrix
#' @param csv_file Path to the read count matrix.
#' @param genes Keep all the genes or only the protein coding genes. Values can be "all" or "protein".
#' @return Return a dataframe.
get_TPM_counts <- function(csv_file, genes){
	# load file
	tpm <- read.csv(file=csv_file, row.names=1)
	# Select protein coding genes
	if (genes=="protein"){
		protein_coding_genes <- read.csv(file=protein_genes)
		tpm <- tpm[rownames(tpm) %in% protein_coding_genes$external_gene_name,]
		rm(protein_coding_genes)
	}
	return(tpm)
}

#' Normalize the read count using the size factor normalization from DESeq2
#' @param raw_counts Read count matrix.
#' @param samplesheet Samplesheet for DESeq2.
#' @return Return a dataframe.
get_normalized_counts <- function(raw_counts, samplesheet){
	dds <- DESeq2::DESeqDataSetFromMatrix(
		countData = raw_counts,
		colData = samplesheet,
		design = ~conditions
	)
	dds <- DESeq2::estimateSizeFactors(dds)
	norm_counts <- DESeq2::counts(dds, normalized=TRUE)
	# norm_counts <- assay(vst(dds, blind=FALSE))

	return(norm_counts)
}

#' Extract the size factor normalization from DESeq2 to later normalize the bigwig files
#' @param raw_counts Read count matrix.
#' @param samplesheet Samplesheet for DESeq2.
#' @return Return a dataframe.
get_size_factors <- function(raw_counts, samplesheet){
	dds <- DESeq2::DESeqDataSetFromMatrix(
		countData = raw_counts,
		colData = samplesheet,
		design = ~conditions
	)
	dds <- DESeq2::estimateSizeFactors(dds)
	size_factors <- DESeq2::sizeFactors(dds)
	# norm_counts <- DESeq2::counts(dds, normalized=TRUE)
	# norm_counts <- SummarizedExperiment::assay(DESeq2::vst(dds, blind=FALSE))
	size_factors <- data.frame(
		sample = names(size_factors),
		SizeFactor = size_factors
	)
	return(size_factors)
}

#' When the maximum expression value (TPM or read count) of a gene between samples is under a certain threshold, we considere it is not properly expressed and the expression values are set to 0.
#' @param row Current gene expression data.
#' @param col_names Save sample names to keep the colnames in the final dataframe.
#' @param minExp Minimum expression value. Default is 10.
#' @return Return a dataframe.
filter_low_counts <- function(row, col_names, minExp) {
	if (max(row) < minExp) {
		return(setNames(rep(0, length(row)), col_names))
	} else {
		return(setNames(row, col_names))
	}
}

#' Apply the expression filtration for lowly expressed genes.
#' @param data Read count or TPM matrix.
#' @param minExp Minimum expression value. Default is 10.
#' @return Return a dataframe.
run_filter_low_counts <- function(data, minExp) {
	col_names <- colnames(data)
	data <- t(apply(data, 1, filter_low_counts, col_names = col_names, minExp = minExp))
	data <- as.data.frame(data)
	data <- data[rowSums(data[])>0,]
	return(data)
}

#################################################################################################################################

###########################################
#                                         #
#              Get matrices               #
#                                         #
###########################################

raw_counts <- get_gene_matrix(
	count_file, 
	"protein"
)

TPM <- get_gene_matrix(
	TPM_file, 
	"protein"
)

###########################################
#                                         #
#       Filter low expressed genes        #
#                                         #
###########################################

# Remove gene expression if max value < x TPM
TPM <- run_filter_low_counts(TPM, minTPM)
kept_genes <- rownames(TPM[rowSums(TPM)>0,])

# Remove gene expression if max value < x reads
raw_counts <- run_filter_low_counts(raw_counts, minReads)

# Remove genes with low TPM from the count matrix
raw_counts <- raw_counts[rownames(raw_counts) %in% kept_genes,]

###########################################
#                                         #
#          Remove outlier samples         #
#                                         #
###########################################

# If there is an outlier sample, remove it from the matrices but also keep the full version to plot them for QCs
if(length(outlierSamples)>0){
	# Keep the full TPM matrix
	TPM_all <- TPM
	# Exclude the outliers
	raw_counts_all <- raw_counts
	raw_counts <- raw_counts[, !colnames(raw_counts) %in% outlierSamples]
	TPM <- TPM[, !colnames(TPM) %in% outlierSamples]
}

###########################################
#                                         #
#             Get samplesheet             #
#                                         #
###########################################

# Generate the reduced samplesheet whith all samples
stage <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 1)
sex <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 2)
conditions <- paste(sex, stage, sep=" ")
conditions_all <- paste(
	sapply(strsplit(colnames(TPM_all), "_"), `[`, 2), 
	sapply(strsplit(colnames(TPM_all), "_"), `[`, 1), 
	sep=" "
)
replicate <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 4)
transgene <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 3)

samplesheet <- data.frame(
	sample = colnames(raw_counts),
	stages = stage,
	sex = sex,
	conditions = conditions,
	replicate = replicate,
	transgene = transgene
)

# Generate the reduced samplesheet whithout outliers
if(length(outlierSamples)>0){
	stage <- sapply(strsplit(colnames(raw_counts_all), "_"), `[`, 1)
	sex <- sapply(strsplit(colnames(raw_counts_all), "_"), `[`, 2)
	conditions <- paste(sex, stage, sep=" ")
	conditions_all <- paste(
		sapply(strsplit(colnames(TPM_all), "_"), `[`, 2), 
		sapply(strsplit(colnames(TPM_all), "_"), `[`, 1), 
		sep=" "
	)
	replicate <- sapply(strsplit(colnames(raw_counts_all), "_"), `[`, 4)
	transgene <- sapply(strsplit(colnames(raw_counts_all), "_"), `[`, 3)

	samplesheet_all <- data.frame(
		sample = colnames(raw_counts_all),
		stages = stage,
		sex = sex,
		conditions = conditions,
		replicate = replicate,
		transgene = transgene
	)
}

###########################################
#                                         #
#           Normalize read counts         #
#                                         #
###########################################

# Generate the normalized read count matrice with the final samples for the analysis (without outliers)
if(length(outlierSamples)>0){
	norm_counts <- get_normalized_counts(
		raw_counts, 
		samplesheet
	)

	norm_counts_all <- get_normalized_counts(
		raw_counts_all, 
		samplesheet_all
	)
} else {
	norm_counts <- get_normalized_counts(
		raw_counts, 
		samplesheet
	)
}

# Extract the size factors used for the normalization
size_factors <- get_size_factors(
	raw_counts, 
	samplesheet
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

# If there are outliers, white the TPM and the normalized counts matrices including them
if(length(outlierSamples)>0){
	write.csv(TPM_all, snakemake@output[['tpm_all']])
	write.csv(norm_counts_all, snakemake@output[['norm_counts_all']])
}

# Save the matrices without ourliers fort the analysis
write.csv(TPM, snakemake@output[['tpm']])
write.csv(raw_counts, snakemake@output[['counts']])
write.csv(norm_counts, snakemake@output[['norm_counts']])
write.csv(samplesheet, snakemake@output[['samplesheet']])
write.csv(size_factors, snakemake@output[['size_factors']])