source("scripts/00.functions.R")

###########################################
#                                         #
#               data files                #
#                                         #
###########################################

count_file <- snakemake@input[['counts']]
TPM_file <- snakemake@input[['tpm']]

###########################################
#                                         #
#              Get matrices               #
#                                         #
###########################################

raw_counts <- get_read_counts(
	count_file, 
	"protein"
)

TPM <- get_read_counts(
	TPM_file, 
	"protein"
)

###########################################
#                                         #
#       Filter low expressed genes        #
#                                         #
###########################################

# Remove gene expression if max value < 5 TPM
TPM <- run_filter_low_counts(TPM, 2)
kept_genes <- rownames(TPM[rowSums(TPM)>0,])

# Remove gene expression if max value < 10 reads
raw_counts <- run_filter_low_counts(raw_counts, 15)

# Remove genes with low TPM from the count matrix
raw_counts <- raw_counts[rownames(raw_counts) %in% kept_genes,]

###########################################
#                                         #
#        Remove XX E11.5 replicate        #
#                                         #
###########################################

TPM_all <- TPM

# Exclude XX E11.5 rep2 sample
raw_counts <- raw_counts[, !colnames(raw_counts) %in% "E11.5_XX_enh8.mcherry_rep2"]
TPM <- TPM[, !colnames(TPM) %in% "E11.5_XX_enh8.mcherry_rep2"]

###########################################
#                                         #
#             Get samplesheet             #
#                                         #
###########################################

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

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

write.csv(TPM_all, snakemake@output[['tpm_all']])
write.csv(TPM, snakemake@output[['tpm']])
write.csv(raw_counts, snakemake@output[['counts']])
write.csv(samplesheet, snakemake@output[['samplesheet']])
