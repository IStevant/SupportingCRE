source(".Rprofile")
###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

raw_counts <- read.csv(file=snakemake@input[['counts']], row.names=1)
samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)

adj.pval <- snakemake@params[['adjpval']]
log2FC <- snakemake@params[['log2FC']]
sex <- snakemake@params[['sex']]

###########################################
#                                         #
#            DESeq analysis               #
#                                         #
###########################################


sex_count <- raw_counts[, grepl(sex , colnames(raw_counts))]
sex_samplesheet <- samplesheet[grepl(sex , samplesheet$sample),]


dds <- DESeq2::DESeqDataSetFromMatrix(
	countData = sex_count,
	colData = sex_samplesheet,
	design = ~stages
)

ddsTC <- DESeq2::DESeq(dds, test="LRT", reduced = ~1)
resTC <- DESeq2::results(ddsTC)
resTC$symbol <- GenomicRanges::mcols(ddsTC)$symbol


filtered_StageDARs <- subset(resTC, padj < adj.pval)
filtered_StageDARs <- subset(filtered_StageDARs, abs(log2FoldChange) > log2FC)

write.csv(filtered_StageDARs[order(filtered_StageDARs$padj),], file=snakemake@output[['csv']], quote = FALSE)

filtered_StageDARs <- rownames(filtered_StageDARs)

save(filtered_StageDARs, file=snakemake@output[['sig_DARs']])
