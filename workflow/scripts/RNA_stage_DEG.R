source(".Rprofile")
###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

raw_counts <- read.csv(file=snakemake@input[['counts']], row.names=1)
samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)
TF_list <- read.csv(file=snakemake@input[['TF_genes']], header=FALSE)$V1
pheno_TFs <- read.table(file=snakemake@input[['TF_pheno']], header=FALSE, sep="\t")
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


filtered_StageDEGs <- subset(resTC, padj < adj.pval)
filtered_StageDEGs <- subset(filtered_StageDEGs, abs(log2FoldChange) > log2FC)

filtered_StageDEGs <- filtered_StageDEGs[,!colnames(filtered_StageDEGs) %in% c("lfcSE","stat")]


filtered <- data.frame(
	filtered_StageDEGs,
	is.TF = ifelse(rownames(filtered_StageDEGs) %in% TF_list, "Yes", "-")
)

rownames(pheno_TFs) <- pheno_TFs[,1]
colnames(pheno_TFs) <- c("genes", "Phenotype")
merged_filtered <- merge(filtered, pheno_TFs[,2, drop=FALSE], by = 0, all.x=TRUE)
merged_filtered[is.na(merged_filtered)] <- "-"
names(merged_filtered)[names(merged_filtered) == 'Row.names'] <- 'Genes'
rownames(merged_filtered) <- merged_filtered$Genes

filtered_StageDEGs <- merged_filtered

write.table(filtered_StageDEGs[order(filtered_StageDEGs$padj),], file=snakemake@output[['tsv']], quote=FALSE, row.names=FALSE, sep="\t")


filtered_StageDEGs <- rownames(filtered_StageDEGs)

save(filtered_StageDEGs, file=snakemake@output[['sig_DEGs']])
