source(".Rprofile")

suppressPackageStartupMessages({
  library("GenomicRanges")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

raw_counts <- read.csv(file = snakemake@input[["counts"]], row.names = 1)
samplesheet <- read.csv(file = snakemake@input[["samplesheet"]], row.names = 1)
load(file = snakemake@input[["peak_list"]])
adj.pval <- snakemake@params[["adjpval"]]
log2FC <- snakemake@params[["log2FC"]]
sex <- snakemake@params[["sex"]]


# raw_counts <- read.csv(file="results/processed_data/mm39/ATAC_raw_counts.csv", row.names=1)
# samplesheet <- read.csv(file="results/processed_data/mm39/ATAC_samplesheet.csv", row.names=1)
# load(file="workflow/data/mm39/ATAC_all_consensus_peaks_2rep_list.Robj")
# adj.pval <- 0.01
# log2FC <- 1
# sex <- "XY"

###########################################
#                                         #
#            DESeq analysis               #
#                                         #
###########################################

# The consensus chromatin regions includes peaks found in both sexes.
# Because we are doing our analysis in only one of the two sexes, we recalculate the consensus peaks present oin the analysed sex, and discard the other regions from the quantification matrix

sex_peaks <- unlist(peak_list[grep(sex, names(peak_list))])
merged_sex_peaks <- reduce(unlist(as(sex_peaks, "GRangesList")))

all_peaks <- GRanges(rownames(raw_counts))

common_peaks <- subsetByOverlaps(all_peaks, merged_sex_peaks)

sex_peaks <- paste0(seqnames(common_peaks), ":", start(common_peaks), "-", end(common_peaks))

sex_count <- raw_counts[, grepl(sex, colnames(raw_counts))]
sex_count <- sex_count[rownames(sex_count) %in% sex_peaks, ]


sex_samplesheet <- samplesheet[grepl(sex, samplesheet$sample), ]


dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = sex_count,
  colData = sex_samplesheet,
  design = ~stages
)

ddsTC <- DESeq2::DESeq(dds, test = "LRT", reduced = ~1)
resTC <- DESeq2::results(ddsTC)
resTC$symbol <- GenomicRanges::mcols(ddsTC)$symbol


filtered_StageDARs <- subset(resTC, padj < adj.pval)
filtered_StageDARs <- subset(filtered_StageDARs, abs(log2FoldChange) > log2FC)

write.csv(filtered_StageDARs[order(filtered_StageDARs$padj), ], file = snakemake@output[["csv"]], quote = FALSE)

filtered_StageDARs <- rownames(filtered_StageDARs)

save(filtered_StageDARs, file = snakemake@output[["sig_DARs"]])
