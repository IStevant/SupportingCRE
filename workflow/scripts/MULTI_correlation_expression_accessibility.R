# For each stage
# 1- get promoters +/- 100bp TSS
# 2- get accessibility + mean triplicates
# 3- get gene expression + mean triplicates
# 4- annotate if sex-DE gene and annotate
# 5- plot log2 acc vs log2 expr + color if gene sex DE

source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  if (snakemake@params[["genome"]] == "mm10") {
    library("BSgenome.Mmusculus.UCSC.mm10")
  } else {
    library("BSgenome.Mmusculus.UCSC.mm39")
  }
  library("dplyr")
  library("cowplot")
  library("grid")
  library("ggplot2")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

TPM <- read.csv(file = snakemake@input[["TPM"]], header = TRUE, row.names = 1)
genome_version <- snakemake@params[["genome"]]
# filtered_SexDEGs
load(snakemake@input[["sex_DEGs"]])
gtf <- snakemake@input[["gtf"]]

atac <- read.csv(file = snakemake@input[["ATAC_raw_counts"]], header = TRUE, row.names = 1)
minReads <- snakemake@params[["minReads"]]

# # filtered_SexDARs
# load(snakemake@input[["atac_sex_DAR"]])

background="genome"

JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db

# Interactive

TPM <- read.csv(file = "results/processed_data/mm10/RNA_TPM.csv", header = TRUE, row.names = 1)
genome_version <- "mm10"
# filtered_SexDEGs
load("results/processed_data/mm10/RNA_sig_SexDEGs.Robj")


atac <- read.csv(file = "results/processed_data/mm10/ATAC_norm_counts.csv", header = TRUE, row.names = 1)
# # filtered_SexDARs
# load("results/processed_data/mm10/ATAC_sig_SexDARs.Robj")


gtf <- "workflow/data/mm10/gencode.vM25.annotation.gtf.gz"
# Promoter region, i.e. distance to TSS
promoter <- 100
genome_gtf <- rtracklayer::import(gtf)
# gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("transcript_id", "gene_name")]
# gene2symbol <- na.omit(as.data.frame(unique(gene2symbol)))
# rownames(gene2symbol) <- gene2symbol$transcript_id


###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' When the maximum value (read count) of a peak between samples is under a certain threshold, we considere it is not relevant and the values are set to 0.
#' @param data Read count matrix.
#' @param minExp Minimum value. Default is 5.
#' @return Return a dataframe.
run_filter_low_counts <- function(data, minExp = 5) {
  col_names <- colnames(data)
  data <- t(apply(data, 1, filter_low_counts, col_names = col_names, minExp = minExp))
  data <- as.data.frame(data)
  data <- data[rowSums(data[]) > 0, ]
  return(data)
}
filter_low_counts <- function(row, col_names, minExp) {
  if (max(row) < minExp) {
    return(setNames(rep(0, length(row)), col_names))
  } else {
    return(setNames(row, col_names))
  }
}

# Taken from https://rdrr.io/github/frankRuehle/systemsbio/src/R/HLP_geneRanges.R
subsetByOverlaps.keepAllMeta <- function(gr1, gr2, write.ranges.tofile = NULL, addStart=0) {
  
  ranges <- subsetByOverlaps(gr1, gr2) # query, subject
  
  hits <- findOverlaps(ranges, gr2) 
  idx <- unique(subjectHits(hits)) 
  
  names <- CharacterList(split(names(gr2)[subjectHits(hits)], queryHits(hits))) # row names of gr2 object added to meta data
  names <- sapply(names, function(x) {paste(unique(x), collapse="; ") })
  mcols(ranges) <- DataFrame(mcols(ranges), names)
  
  for(m in names(mcols(gr2))) { # meta columns of gr2 summarized and added to meta data of gr1
    meta <- mcols(gr2)[subjectHits(hits),m]
    if (is.factor(meta)) {meta <- as.character(meta)}
    meta <- CharacterList(split(meta, queryHits(hits)))
    mcols(ranges) <- DataFrame(mcols(ranges), metaname=meta)
    names(mcols(ranges))[names(mcols(ranges))=="metaname"] <- m
    if(class(mcols(ranges)[,m])=="CompressedCharacterList") {
      mcols(ranges)[,m] <-  sapply(mcols(ranges)[,m], function(x) {paste(unique(x), collapse="; ") })
    } else {
      mcols(ranges)[,m] <-  mcols(ranges)[,m]
    }
  }
  
  if(!is.null(write.ranges.tofile)) {
    df <- granges2df(ranges, addStart=addStart) 
    write.table(df, write.ranges.tofile, sep="\t", quote = F, row.names = F)
  }
  
  return(ranges)
} 


#' @describeIn subsetByOverlaps.keepAllMeta Convert GRanges object to dataframe
granges2df <- function(gr1, addStart=0) {
  
  if(is.null(names(gr1))) {names(gr1) <- 1:length(gr1)}
  df <- data.frame(names=names(gr1),
                   seqnames=seqnames(gr1),
                   start=start(gr1) + addStart,  # -1: BED uses 0-based coordinates
                   end=end(gr1),
                   strand=strand(gr1))
  
  if(addStart!=0) {cat(addStart, "bp added to start coordinate.\n")}
  
  for(m in names(mcols(gr1))) {
    if(class(mcols(gr1)[,m])=="CompressedCharacterList") {
      df[,m] <-  sapply(mcols(gr1)[,m], function(x) {paste(unique(x), collapse="; ") })
    } else {
      df[,m] <-  mcols(gr1)[,m]
    }
  }
  
  return(df)   
}

###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

stages <- unique(sapply(strsplit(colnames(TPM), "_"), `[`, 1))
sex <- unique(sapply(strsplit(colnames(TPM), "_"), `[`, 2))

DEGs <- lapply(filtered_SexDEGs, function(DEG){
  DEGs <- DEG[,c("Genes", "Diff.Exp.")]
  nonDEGs <- rownames(TPM[!rownames(TPM) %in% rownames(DEGs),])
  nonDEGs <- data.frame(
    genes = nonDEGs,
    exp = rep("n.s.", length(nonDEGs))
  )
  colnames(nonDEGs) <- colnames(DEGs)
  all <- rbind(DEGs, nonDEGs)
  rownames(all) <- all$Genes
  return(all)
})

gene_exp <- lapply(stages, function(stage){
  Genes <- DEGs[[stage]]
  expr_XX <- as.data.frame(log2(rowMeans(TPM[,grep(paste0(stage, "_XX"), colnames(TPM))])))
  expr_XX$Genes <- rownames(expr_XX)
  expr_XY <- as.data.frame(log2(rowMeans(TPM[,grep(paste0(stage, "_XY"), colnames(TPM))])))
  expr_XY$Genes <- rownames(expr_XY)

  Gene_exp <- merge(x = Genes, y = expr_XX, by = "Genes")
  Gene_exp <- merge(x = Gene_exp, y = expr_XY, by = "Genes")
  colnames(Gene_exp) <- c("Genes", "Diff.Exp.", "XX_Expr.", "XY_Expr.")
  return(Gene_exp)
})

names(gene_exp) <- stages

TxDb <- GenomicFeatures::makeTxDbFromGFF(gtf)

genes <- GenomicFeatures::genes(TxDb)
gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- na.omit(as.data.frame(unique(gene2symbol)))
rownames(gene2symbol) <- gene2symbol$gene_id

all_TSS <- unique(GenomicRanges::resize(genes, width=2, fix='start'))
GenomicRanges::mcols(all_TSS)$gene_name <- gene2symbol[all_TSS$gene_id, "gene_name"]
all_promoters <- GenomicRanges::promoters(unique(all_TSS), upstream=100, downstream=100)

expr_gene_prom <- unique(all_promoters[all_promoters$gene_name %in% rownames(TPM),])
expr_gene_prom$interval <- as.character(expr_gene_prom, ignore.strand=T)

gr_atac <- GenomicRanges::GRanges(rownames(atac))
gr_atac$interval <- rownames(atac)

gr_atac_prom <- IRanges::subsetByOverlaps(gr_atac, expr_gene_prom)

gr_atac_prom2 <- IRanges::subsetByOverlaps(expr_gene_prom, gr_atac)



atac_prom <- atac[rownames(atac) %in% gr_atac_prom$interval,]

gene_exp_acc <- lapply(stages, function(stage){
  expr_table <- gene_exp[[stage]]
  acc_prom <- data.frame(
    interval = rownames(atac_prom),
    XX_Acc. = as.data.frame(log2(rowMeans(atac_prom[, grep(paste0(stage, "_XX"), colnames(atac_prom))]))),
    XY_Acc. = as.data.frame(log2(rowMeans(atac_prom[, grep(paste0(stage, "_XY"), colnames(atac_prom))])))
  )
  colnames(acc_prom) <- c("interval", "XX_Acc.", "XY_Acc.")

  acc_prom_gr <- GenomicRanges::GRanges(acc_prom$interval)
  acc_prom_gr$XX_Acc. <- acc_prom$XX_Acc.
  acc_prom_gr$XY_Acc. <- acc_prom$XY_Acc.

  acc_gene_prom_gr <- subsetByOverlaps.keepAllMeta(acc_prom_gr, expr_gene_prom)

  acc_gene_prom <- granges2df(acc_gene_prom_gr)
  acc_gene_prom <- acc_gene_prom[,c("gene_name", "XX_Acc.", "XY_Acc.", "interval")]
  colnames(acc_gene_prom) <- c("Genes", "XX_Acc.", "XY_Acc.", "interval")

  expr_acc_table <- merge(expr_table, acc_gene_prom, by="Genes")

  return(expr_acc_table)
})

names(gene_exp_acc) <- stages

plots <- lapply(stages, function(stage){
  expr_table <- gene_exp_acc[[stage]]
  expr_table <- expr_table[expr_table$Diff.Exp.=="Up in XY",]
  p <- ggplot(expr_table, aes(x=XY_Expr., y=XY_Acc., color=Diff.Exp.))+
        geom_point()+
        theme(
          aspect.ratio = 1
        )
})

figure <- plot_grid(
  plotlist = plots,
  labels = "AUTO",
  ncol = 2
)


##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################


save_plot(
  "test.pdf",
  figure,
  base_width = 28,
  base_height = 28,
  units = c("cm"),
  dpi = 300
)


# save_plot(
#   snakemake@output[["png"]],
#   heatmap_list,
#   base_width = 30,
#   base_height = 40,
#   units = c("cm"),
#   dpi = 300,
#   bg = "white"
# )

# save_plot(
#   snakemake@output[["pdf"]],
#   heatmap_list,
#   base_width = 30,
#   base_height = 40,
#   units = c("cm"),
#   dpi = 300
# )
