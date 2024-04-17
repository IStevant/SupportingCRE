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


sig.DE <- subset(resTC, padj < adj.pval)
sig.DE <- subset(sig.DE, abs(log2FoldChange) > log2FC)

write.csv(sig.DE[order(sig.DE$padj),], file=snakemake@output[['csv']], quote = FALSE)

sig.DE <- rownames(sig.DE)

save(sig.DE, file=snakemake@output[['sig_DEGs']])
