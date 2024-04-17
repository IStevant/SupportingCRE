source("scripts/00.functions.R")

raw_counts <- read.csv(file=snakemake@input[['counts']], row.names=1)
samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)

adj.pval <- snakemake@params[['adjpval']]
log2FC <- snakemake@params[['log2FC']]

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

# Save the Robj of the results for reuse
save(SexDEGs, file=snakemake@output[['all_DEGs']])

# Get embryonic stages
stages <- unique(samplesheet$stages)

# For each stages, extract significant DEGs
filtered_SexDEGs <- lapply(stages, function(stg) get_sex_DEG_per_stage(SexDEGs, stg, adj.pval, log2FC))

# For each stages, write DEG results into separated files
export <- lapply(seq_along(stages), function(stg) write.csv(filtered_SexDEGs[stg], paste0("results/RNA_DEG_sex_", stages[stg], ".csv")))

names(filtered_SexDEGs) <- stages
save(filtered_SexDEGs, file=snakemake@output[['sig_DEGs']])