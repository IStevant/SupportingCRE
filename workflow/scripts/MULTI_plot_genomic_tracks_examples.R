source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################
suppressPackageStartupMessages({
  library("Gviz")
  library("GenomicInteractions")
  library("InteractionSet")
  library("dplyr")
  library("doParallel")
  library("foreach")
  library("cowplot")
  library("grid")
  library("ggplotify")
})

doParallel::registerDoParallel(cores = 12)

options(ucscChromosomeNames = FALSE)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

import_bw_files <- function(folder, files, locus, window) {
  locus <- gene_coord[gene_coord$name %in% locus, ]
  full_locus <- locus
  start(full_locus) <- start(full_locus) - window[1]
  end(full_locus) <- end(full_locus) + window[2]

  imports <- foreach(file = files, .final = function(file) setNames(file, files)) %dopar% {
    sex <- sapply(strsplit(file, "_"), `[`, 2)
    stage <- sapply(strsplit(file, "_"), `[`, 1)
    gr <- rtracklayer::import(
      paste(folder, file, sep = "/"),
      which = full_locus
    )
    if (length(gr) < 1) {
      gr <- locus
      gr$score <- 0.1
    }
    gr$cond <- rep(paste("           ", sex, stage), length(gr))
    return(gr)
  }
  sex <- unique(sapply(strsplit(files, "_"), `[`, 2))
  stages <- unique(sapply(strsplit(files, "_"), `[`, 1))
  gr_list <- lapply(sex, function(sx) {
    sex_gr <- imports[grep(sx, names(imports))]
    lapply(stages, function(stg) {
      stg_gr <- sex_gr[grep(stg, names(sex_gr))]
    })
  })
  names(gr_list) <- sex
  names(gr_list[[1]]) <- stages
  names(gr_list[[2]]) <- stages
  return(gr_list)
}

generate_genomic_tracks <- function(gr_list, locus, window, max_score, colors) {
  locus <- locus <- gene_coord[gene_coord$name %in% locus, ]

  cond <- stringr::str_extract(names(gr_list), "^.{8}")
  color <- colors[cond]
  res <- 2000
  gTrack_list <- lapply(
    gr_list,
    function(gr) {
      DataTrack(
        range = gr,
        ucscChromosomeNames = FALSE,
        type = "hist",
        baseline = 0,
        lwd.baseline = 1,
        window = res,
        chromosome = as.character(seqnames(locus)),
        name = unique(gr$cond),
        col.baseline = color,
        col.histogram = 0,
        fill.histogram = color,
        ylim = c(0, trunc(max_score, digit = 4)),
        yTicksAt = c(trunc(max_score, digit = 4)),
        rotation.title = 0,
        lwd = 0,
        alpha = 0.8
      )
    }
  )
  track <- OverlayTrack(gTrack_list)
  return(track)
}


plot_tracks <- function(peak_gr, TxDb, gene2symbol, plot_list, locus, window) {
  locus <- gene_coord[gene_coord$name %in% locus, ]

  full_locus <- locus
  start(full_locus) <- start(full_locus) - window[1]
  end(full_locus) <- end(full_locus) + window[2]

  locus_gr <- full_locus
  peak_gr_locus <- subsetByOverlaps(peak_gr, locus_gr)

  genome_track <- GenomeAxisTrack(
    col = "#333333",
    fontsize = 13,
    lwd = 1,
    distFromAxis = 0.5,
    labelPos = "alternating",
    sizes = 0.3,
    scale = 0.5,
    labelPos = "below"
  )

  gene_track <- GeneRegionTrack(
    TxDb,
    collapseTranscripts = "meta",
    name = "Genes",
    chromosome = as.character(seqnames(locus_gr)),
    start = start(locus_gr),
    end = end(locus_gr),
    fontface.group = "italic",
    col = 0,
    col.line = NULL,
    fill = "#585858",
    # lwd=0.3,
    fontcolor.group = "#333333",
    fontsize.group = 18,
    sizes = 0.4,
    rotation.title = 0,
    thinBoxFeature = "UTR"
    # just.group="above"
  )
  ranges(gene_track)$symbol <- gene2symbol[ranges(gene_track)$gene, "gene_name"]

  peak_track <- AnnotationTrack(
    peak_gr_locus,
    chromosome = as.character(seqnames(peak_gr_locus)),
    start = start(peak_gr_locus),
    end = end(peak_gr_locus),
    strand = as.character(strand(peak_gr_locus)),
    name = "OCRs",
    col.line = NULL,
    col = 0,
    fill = "#5f4780",
    sizes = 0.4,
    rotation.title = 0
  )


  ht <- HighlightTrack(
    trackList = c(plot_list, peak_track, gene_track),
    start = start(peak_gr_locus) - 2,
    end = end(peak_gr_locus) + 2,
    col = 0,
    fill = "#666666",
    chromosome = as.character(seqnames(peak_gr_locus)),
    alpha = 0.1,
    lwd = 0.5
  )

  plot <- plotTracks(
    c(ht, genome_track),
    chromosome = as.character(seqnames(locus_gr)),
    from = start(locus_gr),
    to = end(locus_gr),
    transcriptAnnotation = "symbol",
    background.title = "transparent",
    col.border.title = "transparent",
    col.title = "#333333",
    col.axis = "#333333",
    alpha.title = 1,
    alpha.axis = 1,
    sizes = c(rep(0.4, length(plot_list)), 0.3, 0.4, 0.3),
    cex.title = 0.9,
    cex.id = 0.5,
    cex.axis = 0.7,
    title.width = 1.4
  )


  return(plot)
}

###########################################
###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

bw_folder <- snakemake@params[["bw_folder"]]
genome <- snakemake@input[["genome"]]
peaks <- snakemake@input[["peaks"]]
gene_bed <- snakemake@input[["gene_bed"]]
gene_list <- snakemake@input[["gene_list"]]
save_folder <- snakemake@params[["save_folder"]]

# bw_folder <- "results/processed_data/mm10/ATAC_bigwig"
# # genome <- "workflow/data/mm10/iGenome_mm10_ucsc_genes.gtf.gz"
# genome <- "workflow/data/mm10/gencode.vM25.annotation.gtf.gz"
# peaks <- "results/processed_data/mm10/ATAC_norm_counts.csv"
# gene_bed <- "workflow/data/mm10/gene_standard.bed"
# gene_list <- "workflow/data/gTrack_gene_examples.tsv"

bw_files <- list.files(path = bw_folder, pattern = "REP..bw")
gene_list <- read.table(gene_list, header = TRUE)
window <- c(gene_list$Start, gene_list$End)
genes <- gene_list$Gene

peak_gr <- GRanges(rownames(read.csv(peaks, header = TRUE, row.names = 1)))
gene_coord <- GRanges(rtracklayer::import(gene_bed))

conditions <- paste(
  sapply(strsplit(bw_files, "_"), `[`, 1),
  sapply(strsplit(bw_files, "_"), `[`, 2),
  sep = "_"
)
conditions <- unique(conditions)
conditions <- c(conditions[c(TRUE, FALSE)], conditions[c(FALSE, TRUE)])
names(conditions_color) <- conditions

genome_gtf <- rtracklayer::import(genome)
gene2symbol <- mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
gene2symbol$gene_name <- paste0(gene2symbol$gene_name, " ")
rownames(gene2symbol) <- gene2symbol$gene_id

TxDb <- GenomicFeatures::makeTxDbFromGFF(genome)

# pdf(file=snakemake@output[['pdf']], width=8, height=8)
# pdf(file="test.pdf", width=7.5, height=5)

# plot <- lapply(genes, function(gene){
plot <- foreach(gene = genes) %dopar% {
  locus <- gene_list[gene_list$Gene == gene, ]
  window <- c(locus$Start, locus$End) * 1000
  max_score <- locus$MaxScore

  atac_coverage <- import_bw_files(
    bw_folder,
    bw_files,
    gene,
    window
  )

  if (max_score == 0) {
    max_score <- do.call("c", do.call("c", do.call("c", atac_coverage)))
    max_score <- max(unlist(lapply(max_score, function(x) x$score)))
  }

  XX_track_list <- lapply(
    atac_coverage[["XX"]],
    function(stage) {
      generate_genomic_tracks(
        stage,
        gene,
        window,
        max_score,
        conditions_color
      )
    }
  )

  XY_track_list <- lapply(
    atac_coverage[["XY"]],
    function(stage) {
      generate_genomic_tracks(
        stage,
        gene,
        window,
        max_score,
        conditions_color
      )
    }
  )

  png(file = paste0(save_folder, "/", gene, "_peak_tracks.png"), width = 6.5, height = 3.7, units = "in", res = 150)

  print(paste("Plot", gene))
  plot <- plot_tracks(
    peak_gr,
    TxDb,
    gene2symbol,
    c(XX_track_list, XY_track_list),
    gene,
    window
  )

  dev.off()
  is_plot_ok <- "y"
  write.csv(is_plot_ok, file = snakemake@output[["log"]])
  return(plot)
}
# })
