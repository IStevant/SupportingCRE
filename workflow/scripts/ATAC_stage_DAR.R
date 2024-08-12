source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

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
gtf <- snakemake@input[["gtf"]]

# Promoter region, i.e. distance to TSS
promoter <- snakemake@params[["promoter"]]

genome_gtf <- rtracklayer::import(gtf)
gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id
###########################################
#                                         #
#            DESeq analysis               #
#                                         #
###########################################

# The consensus chromatin regions includes peaks found in both sexes.
# Because we are doing our analysis in only one of the two sexes, we recalculate the consensus peaks present in the analysed sex, and discard the other regions from the quantification matrix

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


sig.DA <- subset(resTC, padj < adj.pval)
sig.DA <- subset(sig.DA, abs(log2FoldChange) > log2FC)

DAR_GR <- GenomicRanges::GRanges(rownames(sig.DA))

TxDb <- GenomicFeatures::makeTxDbFromGFF(gtf)

DA_anno <- as.data.frame(
  ChIPseeker::annotatePeak(
    DAR_GR,
    genomicAnnotationPriority = c("Promoter", "5UTR", "Exon", "Intron", "3UTR", "Downstream", "Intergenic"),
    tssRegion = c(-promoter, 0),
    TxDb = TxDb,
    level = "gene",
    overlap = "all"
  )
)

DA_anno$geneId <- gene2symbol[DA_anno$geneId, "gene_name"]

filtered_StageDARs <- data.frame(
  sig.DA[, -c(3:4)],
  annotation = DA_anno$annotation,
  nearest.gene = DA_anno$geneId,
  distanceToTSS = DA_anno$distanceToTSS
)

write.table(filtered_StageDARs[order(filtered_StageDARs$padj), ], file = snakemake@output[["tsv"]], quote = FALSE, sep="\t")

filtered_StageDARs <- rownames(filtered_StageDARs)

save(filtered_StageDARs, file = snakemake@output[["sig_DARs"]])
