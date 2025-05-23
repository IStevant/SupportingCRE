source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

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
  library("monaLisa")
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
  library("seriation")
})

doParallel::registerDoParallel(cores = 48)

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# filtered_SexDARs
load(snakemake@input[["sig_DARs"]])
# Load RNA-seq TPM matrix to filter the TFs that are expressed in the gonads
TPM <- read.csv(file = snakemake@input[["TPM"]], header = TRUE, row.names = 1)
# # filtered_SexDEGs
# load(snakemake@input[["DEG"]])

# DEGs <- lapply(filtered_SexDEGs, function(stg){
#   df <- data.frame(
#     gene = rownames(stg),
#     sex = stg$Diff.Exp.
#   )
#   return(df)
# })

# DEGs <- do.call(rbind, DEGs)
# rownames(DEGs) <- NULL
# DEGs <- unique(DEGs)

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
#' @param stg String, current embryonic stage.
#' @param TPM Read count matrix.
#' @param minTPM Minimal expression value.
#' @param save_folder Minimum value. Default is 5.
#' @return Return a monaLisa enrichment object.
get_enriched_TFs <- function(DARs, sex, TPM, minTPM, save_folder) {
  # Select the genes expressed at a specific stage for both sexes
  genes <- TPM[, grep(sex, colnames(TPM))]
  # Discard lowly expressed genes
  genes <- run_filter_low_counts(genes, minTPM)
  genes <- rownames(genes[rowSums(genes) > 0, ])
  # Select only the TFs
  sex_TFs <- genes[which(genes %in% TFs)]

  # sex_DEGs <- unique(DEGs[grep(sex,DEGs$sex),"gene"])

  # sex_TFs <-sex_TFs[sex_TFs %in% sex_DEGs]

  # Get all vertebrate TF matrices
  pwms <- TFBSTools::getMatrixSet(
    JASPAR,
    opts = list(
      matrixtype = "PWM",
      tax_group = "vertebrates"
    )
  )

  stages <- names(DARs)
  stage_sex_peaks <- lapply(DARs, function(DAR){
    sex_DAR <- GenomicRanges::GRanges(rownames(DAR[DAR$Diff.Acc. == paste0("More in ",sex), ]))
  })
  names(stage_sex_peaks) <- stages

  # Create one dataframe per embryonic stage
  for (name in names(stage_sex_peaks)) {
    assign(name, stage_sex_peaks[[name]])
  }

  all <- c(E11.5, E12.5, E13.5, E15.5)
  # # generate GRanges objects
  # female <- GenomicRanges::GRanges(rownames(sex_peaks[sex_peaks$Diff.Acc. == "More in XX", ]))
  # male <- GenomicRanges::GRanges(rownames(sex_peaks[sex_peaks$Diff.Acc. == "More in XY", ]))
  # all <- c(female, male)

  # Get the peak sequences
  if (genome_version == "mm10") {
    sequences <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, all)
  } else {
    sequences <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm39, all)
  }

  # Define which sequences are male or female specific
  bins <- rep(stages, c(length(E11.5), length(E12.5), length(E13.5), length(E15.5)))
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
  genes <- TPM[, grep(sex, colnames(TPM))]
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
    E11.5 = rowMeans(TF_exp[, grep(paste0("E11.5_",sex), colnames(TF_exp))]),
    E12.5 = rowMeans(TF_exp[, grep(paste0("E12.5_",sex), colnames(TF_exp))]),
    E13.5 = rowMeans(TF_exp[, grep(paste0("E13.5_",sex), colnames(TF_exp))]),
    E15.5 = rowMeans(TF_exp[, grep(paste0("E15.5_",sex), colnames(TF_exp))])
  )

  rownames(TF_expression) <- rownames(SummarizedExperiment::assay(se2, "negLog10Padj"))

  assays(se2)$expr <- TF_expression

  # Select the motifs enriched with a -Log10Padj > 10
  sel2 <- apply(
    SummarizedExperiment::assay(se2, "negLog10Padj"), 1,
    function(x) max(abs(x), 0, na.rm = TRUE)
  ) > 5.0
  seSel <- se2[sel2, ]

  # Select the motif of the TFs expressed with at least 5 TPM
  sel2 <- apply(
    SummarizedExperiment::assay(seSel, "expr"), 1,
    function(x) max(x, 0, na.rm = TRUE)
  ) > minTPM
  seSel <- seSel[sel2, ]

  # Make TF names upper case
  seSel@elementMetadata$motif.name <- toupper(seSel@elementMetadata$motif.name)

  # Transform expression into log10
  assays(seSel)$expr <- log10(assays(seSel)$expr)

  TF_summary <- data.frame(
    TF.name = seSel@elementMetadata$motif.name,
    TF.matrix = rownames(SummarizedExperiment::assay(seSel, "log2enr")),
    SummarizedExperiment::assay(seSel, "log2enr"),
    10^(-SummarizedExperiment::assay(seSel, "negLog10Padj"))
  )

  write.table(TF_summary, file = paste0(save_folder, "/ATAC_sex_DAR_TF_", sex, "_", background, "_bg.csv"), row.names=FALSE, quote=FALSE, sep="\t")
  return(seSel)
}

#' Merge the enrichment result by TFBS motif similarity and plot the results as heatmap.
#' @param seSel monaLisa enrichment object.
#' @return Return a grid object.
merge_TF_motifs <- function(seSel, stg, sex, save_folder) {
  # Cluster motifs by enrichment
  TF_enrichment <- SummarizedExperiment::assay(seSel, "log2enr")
  TF_pVal <- 10^(-SummarizedExperiment::assay(seSel, "negLog10Padj"))
  colnames(TF_pVal) <- paste0("p-val_", colnames(TF_pVal))
  TF_enrichment <- cbind(TF_enrichment, TF_pVal)
  if (background == "genome") {
    nbCluster <- 1
  } else {
    nbCluster <- 2
  }

  hcl <- hclust(dist(TF_enrichment), method = "ward.D")
  clustering <- cutree(hcl, k = nbCluster)

  motif_sig <- lapply(1:nbCluster, function(cl) {
    TFs <- names(clustering[clustering == cl])
    matrices <- seSel@elementMetadata$motif.pfm[TFs]
    pfms <- universalmotif::convert_motifs(matrices, class = "motifStack-pcm")
    hc <- motifStack::clusterMotifs(pfms)
    phylog <- ade4::hclust2phylog(hc)
    # extract the motif signatures
    motifSig <- motifSignature(pfms, phylog, cutoffPval = 0.001, min.freq = 1)
    # get the signatures from object of motifSignature
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

  # print("Prepare to plot: 4")
  grobL <- lapply(motifs_pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")

  names(grobL) <- rownames(enrichment)

  write.table(enrichment, file = paste0(save_folder, "/ATAC_sex_DAR_TF_merged-motifs_", sex, "_", background, "_bg.csv"), row.names=FALSE, quote=FALSE, sep="\t")

  matrix <- enrichment[,1:4]
  matrix_noNA <- matrix
  row_dend <- hclust(dist(matrix_noNA), method = "ward.D2")
  row_dend <- reorder(row_dend, dist(matrix_noNA), method = "OLO")

  pVal <- enrichment[,5:8]

  matrix[pVal>0.05] <- NA

  # matrix[matrix > 1] <- 1
  # matrix[matrix < (-1)] <- (-1)
  # matrix <- t(scale(t(matrix)))

  hmSeqlogo <- HeatmapAnnotation(
    logo = annoSeqlogo(
      grobL = grobL, which = "row",
      space = unit(0.5, "mm"),
      width = unit(2, "inch")
    ),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "row"
  )

  # bincols <- c(
  #   XX = XX_colors[grep(stg, names(XX_colors))],
  #   XY = XY_colors[grep(stg, names(XY_colors))]
  # )

  bincols <- conditions_color[grepl(sex, names(conditions_color))]

  names(bincols) <- paste(sex, stg)
  conditions <- paste(sex, stg)

  stage_anno <- HeatmapAnnotation(
    Stages = anno_block(
      gp = gpar(fill = bincols, col = 0),
      labels = names(bincols),
      labels_gp = gpar(col = "white", fontsize = 14, fontface = "bold"),
      height = unit(7, "mm")
    )
  )

  cold <- colorRampPalette(c("#04bbc6", "#52d0cf", "#8addd8", "#c1f0e0"))
  warm <- colorRampPalette(c("#fffee8", "#ffd9cb", "#ffb1ad", "#f5808d", "#eb2d62"))
  TYP <- c(cold(2), warm(12))

  mypalette <- colorRamp2(
    breaks = seq(0, 1.5, length.out = 12),
    colors = warm(12)
  )

  # mypalette <- warm(12)

  # Order the rows by decreasing values
  # matrix <- matrix[order(-matrix[[1]], -matrix[[2]], -matrix[[3]], -matrix[[4]]), ]

  ht_list <- Heatmap(
    matrix,
    # clustering_method_rows = "ward.D2",
    name = "Log2 enrichment",
    right_annotation = hmSeqlogo,
    top_annotation = stage_anno,
    column_split = conditions,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = row_dend,
    # row_split = 6,
    row_names_gp = grid::gpar(fontsize = 14),
    column_title = NULL,
    row_title = NULL,
    col = mypalette,
    na_col = "#04bbc6",
    width = ncol(matrix) * unit(25, "mm"),
    row_names_max_width = unit(14, "cm"),
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

stages <- unique(sapply(strsplit(colnames(TPM), "_"), `[`, 1))
sexes <- unique(sapply(strsplit(colnames(TPM), "_"), `[`, 2))

# For each sex, get TFBS motif enrichments
# enrichments <- lapply(sexes, function(sex) {
#   get_enriched_TFs(filtered_SexDARs, sex, TPM, minTPM, save_folder)
# })

# enrichments <- get_enriched_TFs(filtered_SexDARs, "XX", TPM, minTPM, save_folder)

# enrichments <- get_enriched_TFs(filtered_SexDARs, "XY", TPM, minTPM, save_folder)

enrichments <- lapply(sexes, function(sex) {
  get_enriched_TFs(filtered_SexDARs, sex, TPM, minTPM, save_folder)
})


heatmap_list <- lapply(seq_along(sexes), function(i) {
  htm <- merge_TF_motifs(enrichments[[i]], stages, sexes[i], save_folder)
  return(htm)
})


figure <- plot_grid(
  plotlist = heatmap_list,
  labels = "AUTO",
  ncol = 2,
  greedy = FALSE
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 70,
  base_height = 42,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 70,
  base_height = 42,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
