# 2- get promoters
# 3- check if open in the corresponding sex or both (binary or quantitative?)
# get TF motif enrichment in sex-spe and non-sex spe accessible

source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################


suppressPackageStartupMessages({
  library("monaLisa")
  if (snakemake@params[["genome"]] == "mm10") {
    library("BSgenome.Mmusculus.UCSC.mm10")
  } else {
    library("BSgenome.Mmusculus.UCSC.mm39")
  }
  library("dplyr")
  library("cowplot")
  library("grid")
  library("ggplot2")
  library("gtable")
  library("SummarizedExperiment")
  library("ComplexHeatmap")
  library("circlize")
  library("motifStack")
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

# TPM <- read.csv(file = "results/processed_data/mm10/RNA_TPM.csv", header = TRUE, row.names = 1)
# genome_version <- "mm10"
# sex <- "XY"
# # filtered_SexDEGs
# load("results/processed_data/mm10/RNA_sig_SexDEGs.Robj")
# gtf <- "workflow/data/mm10/gencode.vM25.annotation.gtf.gz"
# # Promoter region, i.e. distance to TSS
# promoter <- 100
# genome_gtf <- rtracklayer::import(gtf)
# # gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("transcript_id", "gene_name")]
# # gene2symbol <- na.omit(as.data.frame(unique(gene2symbol)))
# # rownames(gene2symbol) <- gene2symbol$transcript_id

# atac <- read.csv(file = "results/processed_data/mm10/ATAC_raw_counts.csv", header = TRUE, row.names = 1)
# minReads <- 50
# # filtered_SexDARs
# load("results/processed_data/mm10/ATAC_sig_SexDARs.Robj")

# background="genome"

# JASPAR <- JASPAR2020::JASPAR2020
# JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db


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

#' Get TFBS motifs enrichment using the monaLisa package.
#' @param TPM Read count matrix.
#' @param peaks Open chromatin regions linked to genes.
#' @return Return a monaLisa enrichment object.
get_enriched_TFs <- function(TPM, peaks) {
  gene_exp <- TPM
  gene_exp <- run_filter_low_counts(gene_exp, 5)
  gene_exp <- gene_exp[rowSums(gene_exp) > 0, ]

  female <- peaks[[1]]
  male <- peaks[[2]]
  all <- c(female, male)

  # Get the peak sequences
  if (genome_version == "mm10") {
    sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm10, all)
  } else {
    sequences <- getSeq(BSgenome.Mmusculus.UCSC.mm39, all)
  }

  # Define which sequences are male or female specific
  bins <- rep(c("Pre.gran.", "Sertoli"), c(length(female), length(male)))
  bins <- factor(bins)

  # Get all vertebrate TF matrices
  pwms <- TFBSTools::getMatrixSet(
    JASPAR,
    opts = list(
      matrixtype = "PWM",
      tax_group = "vertebrates"
    )
  )

  # Calculate motif enrichments
  # If background is "genome", run the analysis against radom genomic regions
  # Else, compare the two sexes
  if (background == "genome") {
    if (genome_version == "mm10") {
      se2 <- calcBinnedMotifEnrR(
        seqs = sequences,
        bins = bins,
        pwmL = pwms,
        background = "genome",
        genome = BSgenome.Mmusculus.UCSC.mm10,
        genome.regions = NULL, # sample from full genome
        genome.oversample = 2,
        BPPARAM = BiocParallel::MulticoreParam(12)
      )
    } else {
      se2 <- calcBinnedMotifEnrR(
        seqs = sequences,
        bins = bins,
        pwmL = pwms,
        background = "genome",
        genome = BSgenome.Mmusculus.UCSC.mm39,
        genome.regions = NULL, # sample from full genome
        genome.oversample = 2,
        BPPARAM = BiocParallel::MulticoreParam(12)
      )
    }
  } else {
    se2 <- calcBinnedMotifEnrR(
      seqs = sequences,
      bins = bins,
      pwmL = pwms,
      BPPARAM = BiocParallel::MulticoreParam(12)
    )
  }

  # Get gene expression of the TFs
  TF_names <- tolower(as.vector(rowData(se2)$motif.name))
  names(TF_names) <- rowData(se2)$motif.id
  genes <- gene_exp
  rownames(genes) <- tolower(rownames(genes))

  TF_exp <- matrix(0, nrow = length(TF_names), ncol = ncol(genes))
  rownames(TF_exp) <- TF_names
  colnames(TF_exp) <- colnames(genes)

  for (TF in TF_names) {
    if (TF %in% rownames(genes)) {
      TF_exp[TF, ] <- unlist(genes[TF, ])
    }
  }

  rownames(TF_exp) <- names(TF_names)

  TF_expression <- data.frame(
    Pre.gran. = rowMeans(TF_exp),
    Sertoli = rowMeans(TF_exp)
  )

  assays(se2)$expr <- TF_expression

  # Select the motifs enriched with a -Log10Padj > 10
  sel2 <- apply(
    SummarizedExperiment::assay(se2, "negLog10Padj"), 1,
    function(x) max(abs(x), 0, na.rm = TRUE)
  ) > 5
  seSel <- se2[sel2, ]

  # Select the motif of the TFs expressed
  sel2 <- apply(
    SummarizedExperiment::assay(seSel, "expr"), 1,
    function(x) max(x, 0, na.rm = TRUE)
  ) > 0
  seSel <- seSel[sel2, ]

  # Make TF names upper case
  seSel@elementMetadata$motif.name <- toupper(seSel@elementMetadata$motif.name)

  # Transform expression into log10
  assays(seSel)$expr <- log10(assays(seSel)$expr)

  TF_summary <- data.frame(
    SummarizedExperiment::assay(seSel, "log2enr"),
    TF.name = seSel@elementMetadata$motif.name
  )

  return(seSel)
}

#' Merge the enrichment result by TFBS motif similarity and plot the results as heatmap.
#' @param seSel monaLisa enrichment object.
#' @return Return a grid object.
merge_TF_motifs <- function(seSel) {

  # seSel <- enrichments

  # Cluster motifs by enrichment
  TF_enrichment <- SummarizedExperiment::assay(seSel, "log2enr")
  # rownames(TF_enrichment) <- seSel@elementMetadata$motif.name
  if (background == "genome") {
    nbCluster <- 2
  } else {
    nbCluster <- 2
  }

  hcl <- hclust(dist(TF_enrichment))
  clustering <- cutree(hcl, k = nbCluster)

  motif_sig <- lapply(1:nbCluster, function(cl) {
    cl <- 2
    TFs <- names(clustering[clustering == cl])
    matrices <- seSel@elementMetadata$motif.pfm[TFs]
    pfms <- universalmotif::convert_motifs(matrices, class = "motifStack-pcm")
    hc <- motifStack::clusterMotifs(pfms)
    phylog <- ade4::hclust2phylog(hc)
    # extract the motif signatures
    motifSig <- motifSignature(pfms, phylog, cutoffPval = 0.05, min.freq = 1)
    ## get the signatures from object of motifSignature
    sig <- signatures(motifSig)
    return(sig)
  })

  motif_sig_enr <- lapply(motif_sig, function(cl) {
    lapply(cl, function(motif) {
      names <- motif@name
      lapply(names, function(TF) {
        tf_vector <- toupper(unlist(strsplit(TF, ";")))
        enrichment <- TF_enrichment
        rownames(enrichment) <- toupper(seSel@elementMetadata$motif.name)
        tf_enrichment <- enrichment[tf_vector, , drop = FALSE]
        if (nrow(tf_enrichment) > 1) {
          extract_prefix <- function(gene) {
            prefix <- sub("([a-zA-Z]+).*", "\\1", gene)
            if (prefix == "NR") {
              return(gene)
            } else if (grepl("^HOX", prefix)) {
              prefix <- stringr::str_sub(prefix, end = -2)
              return(prefix)
            } else {
              return(prefix)
            }
          }

          grouped_genes <- split(tf_vector, sapply(tf_vector, extract_prefix))

          new_tf_vector <- sapply(grouped_genes, function(gene_group) {
            common_prefix <- extract_prefix(gene_group[1])
            if (common_prefix != "HOX") {
              suffixes <- gsub(paste0("^", common_prefix), "", gene_group)
              numeric_suffixes <- sort(as.numeric(suffixes[grepl("^\\d+$", suffixes)]), na.last = TRUE)
              non_numeric_suffixes <- sort(suffixes[!grepl("^\\d+$", suffixes)])
              all_suffixes <- c(non_numeric_suffixes, numeric_suffixes)
              concatenated_suffixes <- paste(all_suffixes, collapse = "/")
              paste0(common_prefix, concatenated_suffixes)
            } else {
              paste0(common_prefix, "s")
            }
          })

          new_tf_vector <- new_tf_vector[order(new_tf_vector)]
          paste(new_tf_vector, collapse = ";")
          tf_names <- paste(new_tf_vector, collapse = ";")

          tf_enrichment <- as.data.frame(colMedians(tf_enrichment))
          colnames(tf_enrichment) <- tf_names
          tf_enrichment <- as.data.frame(t(tf_enrichment))
          tf_enrichment$mat_names <- TF
        } else {
          tf_enrichment <- as.data.frame(tf_enrichment)
          tf_enrichment$mat_names <- TF
        }

        return(tf_enrichment)
      })
    })
  })

  enrichment <- do.call("rbind", do.call("rbind", do.call("rbind", motif_sig_enr)))
  enrichment <- enrichment[!duplicated(enrichment), ]

  merged_pfms <- do.call("c", do.call("c", do.call("c", motif_sig)))
  names(merged_pfms) <- unlist(lapply(merged_pfms, function(mat) mat@name))


  motifs <- merged_pfms[enrichment$mat_names]
  motifs_pfms <- universalmotif::convert_motifs(motifs, class = "TFBSTools-PFMatrix")

  maxwidth <- max(unlist(lapply(motifs_pfms, function(x) ncol(x@profileMatrix))))

  grobL <- lapply(motifs_pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")

  names(grobL) <- rownames(enrichment)

  matrix <- enrichment[, -3]
  matrix[matrix > 1] <- 1
  matrix[matrix < (-1)] <- (-1)

  hmSeqlogo <- HeatmapAnnotation(
    logo = annoSeqlogo(
      grobL = grobL, which = "row",
      space = unit(0.5, "mm"),
      width = unit(1.5, "inch")
    ),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "row"
  )

  bincols <- c(
    Pos = "#EF6351",
    Neg = "#2191FB"
  )

  names(bincols) <- colnames(matrix)
  conditions <- colnames(matrix)

  stage_anno <- HeatmapAnnotation(
    Stages = anno_block(
      gp = gpar(fill = bincols, col = 0),
      labels = names(bincols),
      labels_gp = gpar(col = "white", fontsize = 14, fontface = "bold"),
      height = unit(7, "mm")
    )
  )

  cold <- colorRampPalette(c("#04bbc6", "#52d0cf", "#8addd8", "#c1f0e0", "#fffee8"))
  warm <- colorRampPalette(c("#fffee8", "#ffd9cb", "#ffb1ad", "#f5808d", "#eb2d62"))
  TYP <- c(cold(12), warm(12))

  mypalette <- colorRamp2(
    breaks = seq(-1, 1, length.out = 24),
    colors = TYP
  )

  print(head(matrix))

  ht_list <- Heatmap(
    matrix,
    clustering_method_rows = "ward.D2",
    name = "Log2 enrichment",
    right_annotation = hmSeqlogo,
    top_annotation = stage_anno,
    column_split = conditions,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    column_title = NULL,
    row_title = NULL,
    col = mypalette,
    width = ncol(matrix) * unit(22, "mm"),
    height = nrow(matrix) * unit(8, "mm"),
    row_names_max_width = unit(17, "cm"),
    heatmap_legend_param = list(direction = "horizontal")
  )

  ht_list_2 <- grid.grabExpr(draw(ht_list, heatmap_legend_side = "bottom"), wrap.grobs = TRUE)

  return(ht_list_2)
}

###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

xx_spe_genes <- unique(unlist(lapply(filtered_SexDEGs, function(x) rownames(x[x$Diff.Exp. == "Up in XX", ]))))
xy_spe_genes <- unique(unlist(lapply(filtered_SexDEGs, function(x) rownames(x[x$Diff.Exp. == "Up in XY", ]))))


TxDb <- GenomicFeatures::makeTxDbFromGFF(gtf)

genome_gtf <- rtracklayer::import(gtf)
genes <- GenomicFeatures::genes(TxDb)
gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- na.omit(as.data.frame(unique(gene2symbol)))
rownames(gene2symbol) <- gene2symbol$gene_id

all_TSS <- unique(GenomicRanges::resize(genes, width=2, fix='start'))

GenomicRanges::mcols(all_TSS)$gene_name <- gene2symbol[all_TSS$gene_id, "gene_name"]

all_promoters <- GenomicRanges::promoters(unique(all_TSS), upstream=100, downstream=0)


xx_DE_gene_prom <- unique(all_promoters[all_promoters$gene_name %in% xx_spe_genes,])
xy_DE_gene_prom <- unique(all_promoters[all_promoters$gene_name %in% xy_spe_genes,])

gr_atac <- GenomicRanges::GRanges(rownames(atac))

xx_OCR_promoters <- IRanges::subsetByOverlaps(gr_atac, xx_DE_gene_prom)
xy_OCR_promoters <- IRanges::subsetByOverlaps(gr_atac, xy_DE_gene_prom)

OCR_promoters <- list(
  Pre.Gran.=xx_OCR_promoters,
  Sertoli=xy_OCR_promoters
)

enrichments <- get_enriched_TFs(TPM, OCR_promoters)

heatmap_list <- merge_TF_motifs(enrichments)


##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################

save_plot(
  snakemake@output[["png"]],
  heatmap_list,
  base_width = 30,
  base_height = 40,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)

save_plot(
  snakemake@output[["pdf"]],
  heatmap_list,
  base_width = 30,
  base_height = 40,
  units = c("cm"),
  dpi = 300
)
