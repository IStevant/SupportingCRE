source(".Rprofile")
# source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Extract the significantly differentially accessible peaks and annotate the result table to indicate in which sex the peak is more accessible.
#' @param dds DESeq2 result object.
#' @param stage Embryonic stage.
#' @param p.adj Maximal adjusted p-value threshold.
#' @param log2FC Minimal log2FoldChange threshold.
#' @return Return a datatable.
get_sex_DAR_per_stage <- function(dds, stage, p.adj, log2FC){
	res <- DESeq2::results(dds, contrast=c("conditions", paste("XX", stage), paste("XY", stage)))

	da_res <- as.data.frame(res)
	res <- dplyr::mutate(
		da_res, 
		Diff.Acc. = dplyr::case_when(
			log2FoldChange >= log2FC & padj <= p.adj ~ "More in XX",
			log2FoldChange <= (-log2FC) & padj <= p.adj ~ "More in XY",
			TRUE ~ "non sig."
		)
	)
	sig.DA <- subset(res, padj < p.adj)
	sig.DA <- subset(sig.DA, abs(log2FoldChange) > log2FC)

	return(sig.DA)
}

get_GR <- function(dataframe){
	# Test is the dataframe is not empty
	if (nrow(dataframe)>0){
		# Get DAR names
		DAR <- rownames(dataframe)
		# Expression is null
		expr <- rep(0, length(DAR))
	# If the dataframe is emtpy
	} else {
		# Generate an empty GR object
		gr <- GenomicRanges::GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL, name=NULL))
		return(gr)
	}
	# generate the GR object by extracting the DAR locus coordinates from its name
	gr <- GenomicRanges::GRanges(
		DAR,
		name = DAR,
		Diff.Acc.=dataframe$Diff.Acc.
	)
	return(gr)
}



#################################################################################################################################

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

###########################################
#                                         #
#     DESeq2 analysis sex per stages      #
#                                         #
###########################################

SexDARs <- DESeq2::DESeqDataSetFromMatrix(
		countData = raw_counts,
		colData = samplesheet,
		design = ~conditions
)

SexDARs <- DESeq2::DESeq(SexDARs)

# Get embryonic stages
stages <- unique(samplesheet$stages)

# For each stages, extract significant DARs
filtered_SexDARs <- lapply(stages, function(stg) get_sex_DAR_per_stage(SexDARs, stg, adj.pval, log2FC))

# For each stages, write DAR results into separated files
export <- lapply(seq_along(stages), function(stg) write.csv(filtered_SexDARs[stg], paste0(save_folder, "/ATAC_DAR_sex_", stages[stg], ".csv")))

names(filtered_SexDARs) <- stages
save(filtered_SexDARs, file=snakemake@output[['sig_DARs']])

sexDAR_GR_list <- lapply(filtered_SexDARs, get_GR)
XX_peaks <- lapply(sexDAR_GR_list, function(stgGR) stgGR[(GenomicRanges::elementMetadata(stgGR)[, "Diff.Acc."] == "More in XX")])
XY_peaks <- lapply(sexDAR_GR_list, function(stgGR) stgGR[(GenomicRanges::elementMetadata(stgGR)[, "Diff.Acc."] == "More in XY")])

peak_list <- list(
	XX=XX_peaks,
	XY=XY_peaks
)
save(peak_list, file=snakemake@output[['sig_DARs_GR']])
