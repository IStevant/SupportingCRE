source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

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
  library("doParallel")
  library("foreach")
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

doParallel::registerDoParallel(cores = 12)

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# filtered_SexDARs
load(snakemake@input[["sig_DARs"]])
# Load RNA-seq TPM matrix to filter the TFs that are expressed in the gonads
TPM <- read.csv(file = snakemake@input[["TPM"]], header = TRUE, row.names = 1)
# Load the mouse TFs list
TFs <- as.vector(read.csv(snakemake@input[["TF_genes"]], header = FALSE)[, 1])
# Minimum TPM value to considere a gene expressed
minTPM <- snakemake@params[["minTPM"]]
# Run analysis using the genome bakground or calculating the enrichment compared to the conditions
background <- snakemake@params[["background"]]

save_folder <- snakemake@params[["save_folder"]]

genome_version <- snakemake@params[["genome"]]

stage <- sapply(strsplit(colnames(TPM), "_"), `[`, 1)
sex <- sapply(strsplit(colnames(TPM), "_"), `[`, 2)
conditions <- unique(paste(sex, stage, sep = " "))
names(conditions_color) <- conditions[order(conditions)]
XX_colors <- conditions_color[grepl("XX", names(conditions_color))]
XY_colors <- conditions_color[grepl("XY", names(conditions_color))]

# Load Jaspar 2024 database
JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db

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
  return(as.data.frame(data))
}
filter_low_counts <- function(row, col_names, minExp) {
  if (max(row) < minExp) {
    return(setNames(rep(0, length(row)), col_names))
  } else {
    return(setNames(row, col_names))
  }
}

#' Get TFBS motifs enrichment using the monaLisa package.
#' @param DARs Table of the differentially accessible regions, with the peak coordinates as rownames.
#' @param TPM Read count matrix.
#' @param minTPM Minimal expression value.
#' @param save_folder Minimum value. Default is 5.
#' @return Return a monaLisa enrichment object.
get_enriched_TFs <- function(DARs, TPM, minTPM, save_folder) {
  # Select the genes expressed at a specific stage for both sexes
  genes <- run_filter_low_counts(TPM, minTPM)
  # Select only the TFs
  stg_TFs <- genes[which(genes %in% TFs)]

  # Get all vertebrate TF matrices
  pwms <- TFBSTools::getMatrixSet(
    JASPAR,
    opts = list(
      matrixtype = "PWM",
      tax_group = "vertebrates"
    )
  )

  SexDARs <- lapply(DARs, function(sexDARs) {
    sexDARs <- data.frame(
      sexDARs,
      coord = rownames(sexDARs)
    )
  })

  sex_peaks <- do.call(rbind, SexDARs)

  # generate GRanges objects
  female <- unique(GenomicRanges::GRanges(sex_peaks[sex_peaks$Diff.Acc. == "More in XX", "coord"]))
  male <- unique(GenomicRanges::GRanges(sex_peaks[sex_peaks$Diff.Acc. == "More in XY", "coord"]))
  all <- c(female, male)

  # Get the peak sequences
  if (genome_version == "mm10") {
    sequences <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, all)
  } else {
    sequences <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm39, all)
  }

  # Define which sequences are male or female specific
  bins <- rep(c("XX", "XY"), c(length(female), length(male)))
  bins <- factor(bins)

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
  genes <- TPM
  rownames(genes) <- tolower(rownames(genes))

  TF_exp <- matrix(0, nrow = length(TF_names), ncol = ncol(genes))
  rownames(TF_exp) <- TF_names
  colnames(TF_exp) <- colnames(genes)

  for (TF in TF_names) {
    if (TF %in% rownames(genes)) {
      TF_exp[TF, ] <- unlist(genes[TF, ])
    }
  }

  TF_expression <- data.frame(
    XX = rowMeans(TF_exp[, grep("XX", colnames(TF_exp))]),
    XY = rowMeans(TF_exp[, grep("XY", colnames(TF_exp))])
  )

  rownames(TF_expression) <- rownames(SummarizedExperiment::assay(se2, "negLog10Padj"))

  assays(se2)$expr <- TF_expression

  # Select the motifs enriched with a -Log10Padj > 10
  sel2 <- apply(
    SummarizedExperiment::assay(se2, "negLog10Padj"), 1,
    function(x) max(abs(x), 0, na.rm = TRUE)
  ) > 10.0
  seSel <- se2[sel2, ]

  # Select the motif of the TFs expressed with at least 5 TPM
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

  write.csv(TF_summary, file = paste0(save_folder, "/ATAC_sex_DAR_TF_", background, "_bg.csv"))
  return(seSel)
}

#' Merge the enrichment result by TFBS motif similarity and plot the results as heatmap.
#' @param seSel monaLisa enrichment object.
#' @return Return a grid object.
merge_TF_motifs <- function(seSel) {
  # Cluster motifs by enrichment
  TF_enrichment <- SummarizedExperiment::assay(seSel, "log2enr")
  # rownames(TF_enrichment) <- seSel@elementMetadata$motif.name
  if (background == "genome") {
    nbCluster <- 4
  } else {
    nbCluster <- 2
  }

  hcl <- hclust(dist(TF_enrichment))
  clustering <- cutree(hcl, k = nbCluster)

  motif_sig <- lapply(1:nbCluster, function(cl) {
    TFs <- names(clustering[clustering == cl])
    # TFs <- names(clustering[clustering==cl])
    # print(TFs)
    matrices <- seSel@elementMetadata$motif.pfm[TFs]
    # print(length(matrices))
    pfms <- universalmotif::convert_motifs(matrices, class = "motifStack-pcm")
    hc <- motifStack::clusterMotifs(pfms)
    phylog <- ade4::hclust2phylog(hc)
    # extract the motif signatures
    motifSig <- motifSignature(pfms, phylog, cutoffPval = 0.001, min.freq = 1)
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
            } else if (grepl("^HOX|^FOX", prefix)) {
              prefix <- stringr::str_sub(prefix, end = -2)
              return(prefix)
            } else {
              return(prefix)
            }
          }

          grouped_genes <- split(tf_vector, sapply(tf_vector, extract_prefix))

          new_tf_vector <- sapply(grouped_genes, function(gene_group) {
            common_prefix <- extract_prefix(gene_group[1])
            if (grepl("^HOX|^FOX", common_prefix)) {
              paste0(common_prefix, "s")
            } else {
              suffixes <- gsub(paste0("^", common_prefix), "", gene_group)
              numeric_suffixes <- sort(as.numeric(suffixes[grepl("^\\d+$", suffixes)]), na.last = TRUE)
              non_numeric_suffixes <- sort(suffixes[!grepl("^\\d+$", suffixes)])
              all_suffixes <- c(non_numeric_suffixes, numeric_suffixes)
              concatenated_suffixes <- paste(all_suffixes, collapse = "/")
              paste0(common_prefix, concatenated_suffixes)            }
          })

          new_tf_vector <- new_tf_vector[order(new_tf_vector)]
          paste(new_tf_vector, collapse = ";")
          tf_names <- paste(new_tf_vector, collapse = ";")

          tf_enrichment <- as.data.frame(colMedians(tf_enrichment))
          colnames(tf_enrichment) <- tf_names
          tf_enrichment <- as.data.frame(t(tf_enrichment))
          tf_enrichment$mat_names <- TF
          # print(tf_enrichment)
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

  fileConn <- file(snakemake@output[["matrices_txt"]], "w")

  for (i in seq_along(motifs_pfms)) {
    pfm <- motifs_pfms[[i]]
    # tf_name <- pfm@name
    tf_name <- rownames(enrichment)[i]
    writeLines(paste0(">", tf_name), fileConn)
    pfm_matrix <- pfm@profileMatrix
    bases <- rownames(pfm_matrix)
    for (j in 1:nrow(pfm_matrix)) {
      pfm_values <- paste(pfm_matrix[j, ], collapse = "   ")
      writeLines(paste(bases[j], "[", pfm_values, "]"), fileConn)
    }
  }

  close(fileConn)

  maxwidth <- max(unlist(lapply(motifs_pfms, function(x) ncol(x@profileMatrix))))

  grobL <- lapply(motifs_pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")

  names(grobL) <- rownames(enrichment)

  matrix <- enrichment[, -3]
  matrix[matrix > 1] <- 1
  matrix[matrix < (-1)] <- (-1)

  matrix <- matrix[order(matrix$XX, decreasing = TRUE),]

  save(matrix, file=snakemake@output[["heatmap_matrice"]])

  grobL <- grobL[order(match(names(grobL),rownames(matrix)))]

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
    XX = "#FFB100",
    XY = "#339989"
  )
  names(bincols) <- c("Pre-gr.", "Sert.")
  conditions <- c("Pre-gr.", "Sert.")

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
    cluster_rows = FALSE,
    column_title = NULL,
    row_title = NULL,
    col = mypalette,
    width = ncol(matrix) * unit(20, "mm"),
    height = nrow(matrix) * unit(5, "mm"),
    row_names_max_width = unit(17, "cm"),
    heatmap_legend_param = list(direction = "horizontal")
  )

  ht_list_2 <- grid.grabExpr(draw(ht_list, heatmap_legend_side = "bottom"), wrap.grobs = TRUE)
  return(ht_list_2)
}

#################################################################################################################################

###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

stages <- unique(sapply(strsplit(colnames(TPM), "_"), `[`, 1))

enrichments <- get_enriched_TFs(filtered_SexDARs, TPM, minTPM, save_folder)

heatmap_list <- merge_TF_motifs(enrichments)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  heatmap_list,
  base_width = 30,
  base_height = 22,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  heatmap_list,
  base_width = 30,
  base_height = 22,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
