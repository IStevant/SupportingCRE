source(".Rprofile")
# source("scripts/00.color_palettes.R")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

raw_counts <- read.csv(file=snakemake@input[['counts']], row.names=1)
samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)
adj.pval <- snakemake@params[['adjpval']]
log2FC <- snakemake@params[['log2FC']]
save_folder <- snakemake@params[['save_folder']]

#################################################################################################################################

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Extract the significantly differentially expressed genes and annotate the result table to indicate in which sex the gene is up.
#' @param dds DESeq2 result object.
#' @param stage Embryonic stage.
#' @param p.adj Maximal adjusted p-value threshold.
#' @param p.adj Minimal log2FoldChange threshold.
#' @return Return a datatable.
get_sex_DEG_per_stage <- function(dds, stage, p.adj, log2FC){
	res <- DESeq2::results(dds, contrast=c("conditions", paste("XX", stage), paste("XY", stage)))

	de_res <- as.data.frame(res)
	res <- dplyr::mutate(
		de_res, 
		Diff.Exp. = dplyr::case_when(
			log2FoldChange >= log2FC & padj <= p.adj ~ "Up in XX",
			log2FoldChange <= (-log2FC) & padj <= p.adj ~ "Up in XY",
			TRUE ~ "non sig."
		)
	)
	sig.DE <- subset(res, padj < p.adj)
	sig.DE <- subset(sig.DE, abs(log2FoldChange) > log2FC)

	return(sig.DE)
}

#################################################################################################################################

###########################################
#                                         #
#     DESeq2 analysis sex per stages      #
#                                         #
###########################################

SexDEGs <- DESeq2::DESeqDataSetFromMatrix(
		countData = raw_counts,
		colData = samplesheet,
		design = ~conditions
)

SexDEGs <- DESeq2::DESeq(SexDEGs)

# Get embryonic stages
stages <- unique(samplesheet$stages)

# For each stages, extract significant DEGs
filtered_SexDEGs <- lapply(stages, function(stg) get_sex_DEG_per_stage(SexDEGs, stg, adj.pval, log2FC))

##########################################
#                                        #
#               Save files               #
#                                        #
##########################################

# For each stages, write DEG results into separated files
export <- lapply(seq_along(stages), function(stg) write.csv(filtered_SexDEGs[stg], paste0(save_folder, "/RNA_DEG_sex_", stages[stg], ".csv")))

# Save the Robj of the results for reuse
save(SexDEGs, file=snakemake@output[['all_DEGs']])
names(filtered_SexDEGs) <- stages
save(filtered_SexDEGs, file=snakemake@output[['sig_DEGs']])
