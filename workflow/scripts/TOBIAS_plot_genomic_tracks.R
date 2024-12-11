source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

# This script is not so optimal for Snakemake but it works...

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
  library("ggplot2")
})

doParallel::registerDoParallel(cores = 12)
options(ucscChromosomeNames = FALSE)


###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

genome <- snakemake@input[["genome"]]
peaks <- snakemake@input[["peaks"]]
peak_list <- snakemake@input[["peak_list"]]
TPM <-  snakemake@input[["TPM"]]
links <- snakemake@input[["linkage"]]
gene_bed <- snakemake@input[["gene_bed"]]

ChIP_folder <- snakemake@params[["ChIP_folder"]]

bw_folder <- snakemake@params[["bw_folder"]]
save_folder <- snakemake@params[["save_folder"]]

XX_footprint_bed <- snakemake@input[["XX_bound"]]
XY_footprint_bed <- snakemake@input[["XY_bound"]]

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Import the coverage data for a given genomic locus from bigwig files.
#' @param folder Path to the bigwig files.
#' @param files Bigwig file names.
#' @param full_locus GRanges coordnate of the extended locus of interest.
#' @return Return a list of GRanges objects.
import_bw_files <- function(folder, files, full_locus) {
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
    # Add space to correct the track name position
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
  return(gr_list)
}

#' Create the coverage tracks with replicates overlay.
#' @param gr_list List containing the extended locus coverage per condition.
#' @param locus GRanges object with the coordinate of the locus of interest.
#' @param max_score Maximal score to set the Y axis limit.
#' @param colors vector of hexadecimal colors.
#' @return Return a Gviz object.
generate_genomic_tracks <- function(gr_list, locus, max_score, colors) {
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

ChIP_tracks <- function(ChIP_gr, locus, name){
  print("1")
  ChIP_gr_locus <- subsetByOverlaps(ChIP_gr, locus)
  ChIP_track <- AnnotationTrack(
    ChIP_gr_locus,
    chromosome = as.character(seqnames(ChIP_gr_locus)),
    start = start(ChIP_gr_locus),
    end = end(ChIP_gr_locus),
    strand = as.character(strand(ChIP_gr_locus)),
    name = name,
    col.line = NULL,
    col = 0,
    # fill = colors,
    sizes = 0.6,
    rotation.title = 0
  )
}

Footprint_tracks <- function(footprint_gr, locus, name){
  print("2")
  footprint_gr_locus <- subsetByOverlaps(footprint_gr, locus)
  footprint_track <- AnnotationTrack(
    footprint_gr_locus,
    chromosome = as.character(seqnames(footprint_gr_locus)),
    start = start(footprint_gr_locus),
    end = end(footprint_gr_locus),
    strand = as.character(strand(footprint_gr_locus)),
    name = name,
    id = footprint_gr_locus$name,
    featureAnnotation = "id",
    showFeatureId = TRUE, 
    just.feature = "above",
    col.line = NULL,
    col = 0,
    fontcolor.feature = "#333333",
    fontsize = 11,
    fill = "#000080",
    sizes = 0.6,
    rotation.title = 0
  )
}


#' Plot the genomic tracks with the peaks and the genome annotation.
#' @param peak_gr GRanges object with the coordinate of the ATAC peaks.
#' @param TxDb Genome annotation.
#' @param gene2symbol Data frame containing the conversion of the Ensembl gene IDs to gene symbols.
#' @param plot_list Gviz objects generated using generate_genomic_tracks().
#' @param locus GRanges coordnate of the extended locus of interest.
#' @return Return a Grid object.
plot_tracks <- function(peak_gr, ChIP_gr, footprint_gr, TxDb, gene2symbol, link, plot_list, locus) {
  locus_gr <- locus
  print("3")
  peak_gr_locus <- subsetByOverlaps(peak_gr, locus_gr)
  gene_links <- links[links$Gene %in% locus$Gene, ]
  print(paste(nrow(gene_links), "links found."))

  ChIP_tracks <- lapply(seq_along(ChIP_gr), function(ChIP) ChIP_tracks(ChIP_gr[[ChIP]], locus_gr, names(ChIP_gr)[ChIP]))
  print("ChIP tracks OK")

  footprint_tracks <- lapply(seq_along(footprint_gr), function(foot) Footprint_tracks(footprint_gr[[foot]], locus_gr, names(footprint_gr)[foot]))
  print("Footprint tracks OK")
  # ChIP_tracks <- lapply(ChIP_gr, function(ChIP) ChIP_tracks(ChIP, GRanges("chr1:1245678-1255678")))

  # If links exist, plot link track
  if (nrow(gene_links) >= 1) {
    cor_type <- gene_links$correlations
    cor_type[cor_type > 0] <- "Positive"
    cor_type[cor_type < 0] <- "Negative"
    colors <- cor_type
    colors[colors == "Positive"] <- "#EF6351"
    colors[colors == "Negative"] <- "#2191FB"

    link_peak_gr <- GRanges(gene_links$Peak)
    link_peak_gr$type <- cor_type
    link_peak_gr$colors <- colors
    link_gene_TSS <- GRanges(seqname = seqnames(link_peak_gr), ranges = IRanges::IRanges(start = gene_links$TSS, end = gene_links$TSS + 1))
    link_gene_TSS$type <- rep("Promoter", length(link_gene_TSS))
    interactions <- GenomicInteractions(link_peak_gr, link_gene_TSS, counts = abs(gene_links$correlations))
    annotateInteractions(
      interactions,
      list(
        positive = link_peak_gr[link_peak_gr$type == "Positive", ],
        negative = link_peak_gr[link_peak_gr$type == "Negative", ]
      )
    )
    interaction_track <- InteractionTrack(interactions, name = "Linkage")

    displayPars(interaction_track) <- list(
      plot.anchors = FALSE,
      plot.outside = TRUE,
      col.interaction.types = c("positive-distal" = "#EF6351", "negative-distal" = "#2191FB"),
      rotation.title = 0,
      legend = TRUE
    )
    print("4")
    visible_linked_peaks <- subsetByOverlaps(link_peak_gr, locus_gr)

  } else {
    visible_linked_peaks <- NULL
  }

  genome_track <- GenomeAxisTrack(
    col = "#333333",
    fontsize = 13,
    lwd = 1,
    distFromAxis = 0.5,
    labelPos = "alternating",
    sizes = 0.5,
    scale = 0.8,
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
    fontcolor.group = "#333333",
    fontsize.group = 18,
    # sizes = 0.9,
    rotation.title = 0,
    thinBoxFeature = "UTR",
    just.group = "above"
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
    sizes = 0.2,
    rotation.title = 0
  )

  title <- paste0(seqnames(locus_gr), ":", start(locus_gr), "-", end(locus_gr))

  if (length(visible_linked_peaks) >= 1) {

    ht <- HighlightTrack(
      trackList = c(gene_track, plot_list, peak_track, footprint_tracks),
      # start = c(start(visible_linked_peaks), start(promoter)),
      start = start(visible_linked_peaks),
      # end = c(end(visible_linked_peaks), end(promoter)),
      end = end(visible_linked_peaks),
      col = c(visible_linked_peaks$colors),
      fill = 0,
      chromosome = as.character(unique(seqnames(visible_linked_peaks))),
      alpha = 0.4,
      lwd = 0.5
    )

    plot <- plotTracks(
      c(interaction_track, ht, genome_track),
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
      sizes = c(0.4, 0.4, rep(0.4, length(plot_list)), 0.2, 0.25, 0.25, 0.3),
      cex.title = 0.9,
      cex.id = 0.5,
      title.width = 1.4,
      main = title,
      cex.main = 1,
      col.main = "#333333"
    )

  } else {
    # ht <- HighlightTrack(
    #   trackList = c(peak_track, gene_track, plot_list),
    #   start = start(promoter),
    #   end = end(promoter),
    #   col = "#000000",
    #   fill = "#000000",
    #   chromosome = as.character(seqnames(promoter)),
    #   alpha = 0.2,
    #   lwd = 0.5
    # )

    plot <- plotTracks(
      c(gene_track, plot_list, peak_track, footprint_tracks, genome_track),
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
      sizes = c(0.4, rep(0.4, length(plot_list)), 0.2, 0.25, 0.25, 0.3),
      cex.title = 0.9,
      cex.id = 0.5,
      title.width = 1.4,
      main = title,
      cex.main = 1,
      col.main = "#333333"
    )
  }
  return(plot)
}

#' Plot gene expression as horizontal TPM to plot side by side with the genomic tracks.
#' @param TPM GRanges object with the coordinate of the ATAC peaks.
#' @param gene Genome annotation.
#' @return Return a Ggplot object.
plot_gene_expression <- function(TPM, gene, sex, conditions_color){
  TPM <- read.csv(TPM, header=TRUE, row.names=1)
  TPM <- TPM[,grep("E13", colnames(TPM))]

  if (sex=="all") {
    exp <- as.numeric(TPM[gene, ])
    sex <- sapply(strsplit(colnames(TPM), "_"), `[`, 2)
  } else {
    TPM <- TPM[,grep(sex, colnames(TPM))]
    exp <- as.numeric(TPM[gene, ])
    sex <- sapply(strsplit(colnames(TPM), "_"), `[`, 2)
  }

  names(conditions_color) <- paste0(
    sapply(strsplit(names(conditions_color), "_"), `[`, 2),
    "_",
    sapply(strsplit(names(conditions_color), "_"), `[`, 1)
  )

    gene_exp <- data.frame(
      sex = sex,
      stages = sapply(strsplit(colnames(TPM), "_"), `[`, 1),
      cond = paste0(sex, "_",sapply(strsplit(colnames(TPM), "_"), `[`, 1)),
      exp = exp
    )

    df <- dplyr::group_by(gene_exp, cond, sex)
    options(dplyr.summarise.inform = FALSE)
    df.summary2 <- dplyr::summarise(
      df,
      sd = sd(exp),
      len = mean(exp)
    )

    plot <- ggplot(df.summary2, aes(x = len, y = cond, color = cond, fill = cond)) +
      geom_bar(stat="identity", width=0.7) +
      geom_errorbar(
        aes(
          xmin = len - sd,
          xmax = len + sd
        ),
        width = .3
      ) +
      scale_color_manual(
        values = conditions_color
      ) +
      scale_fill_manual(
        values = conditions_color
      ) +
      ggtitle(paste0("\n\n\n\n\n\n\n\n\n", gene, "\n")) +
      xlab("TPM\n\n\n\n\n\n\n") +
      scale_x_continuous(labels = scales::comma, n.breaks = 3) +
      scale_y_discrete(limits=rev, expand = c(0,0)) +
      # coord_cartesian(ylim = c(NA, 0)) +
      theme_light() +
      theme(
        plot.title = element_text(size = 12, face = "bold.italic", hjust = 0.5),
        axis.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        # aspect.ratio = 3.5,
        legend.position = "none"
      )
    return(plot)
}

###########################################
#                                         #
#           Plot genomic tracks           #
#                                         #
###########################################

all_bw_files <- list.files(path = bw_folder, pattern = "REP..bw")
bw_files <- list.files(path = bw_folder, pattern = "E13.*REP..bw")
# bw_files <- list.files(path = bw_folder, pattern = "REP..bw")

links <- read.table(links, header = TRUE)

ChIP_files <- list.files(path = ChIP_folder, pattern = "*.bed")
ChIP_gr <- lapply(ChIP_files, function(bed) rtracklayer::import(paste0(ChIP_folder, bed)))
names(ChIP_gr) <- ChIP_files


footprint_gr <- lapply(c(XX_footprint_bed, XY_footprint_bed), rtracklayer::import)
names(footprint_gr) <- c("XX_footprints", "XY_footprints")



all_conditions <- paste(
  sapply(strsplit(all_bw_files, "_"), `[`, 1),
  sapply(strsplit(all_bw_files, "_"), `[`, 2),
  sep = "_"
)

conditions <- unique(all_conditions)
conditions <- c(conditions[c(TRUE, FALSE)], conditions[c(FALSE, TRUE)])
names(conditions_color) <- conditions

conditions_color <- conditions_color[grep("E13", names(conditions_color))]

peak_list <- read.table(peak_list, header = TRUE)
peak2plot <- peak_list$Locus


peak_gr <- GRanges(rownames(read.csv(peaks, header = TRUE, row.names = 1)))


genome_gtf <- rtracklayer::import(genome)
gene2symbol <- mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
gene2symbol$gene_name <- paste0(gene2symbol$gene_name, " ")
rownames(gene2symbol) <- gene2symbol$gene_id

TxDb <- GenomicFeatures::makeTxDbFromGFF(genome)

print("GTF file load: OK")


plot <- foreach(peak = peak2plot) %dopar% {

  sex <- peak_list[peak_list$Locus == peak, "Sex"]

  if (sex != "all") {
    bw_files <- grep(sex, bw_files, value = TRUE)
    conditions_color <- conditions_color[grep(sex, names(conditions_color))]
  }

  print("Defining conditions: OK")

  locus <- peak_list[peak_list$Locus == peak, ]
  peak_list_gr <- GRanges(locus$Locus)
  peak_list_gr$Gene <- locus$Gene
  window <- peak_list_gr

  print("Defining loci: OK")

  atac_coverage <- import_bw_files(
    bw_folder,
    bw_files,
    window
  )

  print("Import bw files: OK")

  max_score <- peak_list[peak_list$Locus == peak, "Max"]

  if (max_score == 0){
    max_score <- do.call("c", do.call("c", do.call("c", atac_coverage)))
    max_score <- max(unlist(lapply(max_score, function(x) x$score)))
  } 

  if (sex == "all") {
    XX_track_list <- lapply(
      atac_coverage[["XX"]],
      function(stage) {
        generate_genomic_tracks(
          stage,
          peak_list_gr,
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
          peak_list_gr,
          max_score,
          conditions_color
        )
      }
    )

    bw_track <- c(XX_track_list, XY_track_list)
    pdf_height <- 4

  } else {
    bw_track <- lapply(
      atac_coverage[[sex]],
      function(stage) {
        generate_genomic_tracks(
          stage,
          peak_list_gr,
          max_score,
          conditions_color
        )
      }
    )
    pdf_height <- 3
  }

  print("Make ATAC tracks: OK")


  gene <- peak_list[peak_list$Locus == peak, "Gene"]

  plot <- grid::grid.grabExpr(
      plot_tracks(
        peak_gr,
        ChIP_gr,
        footprint_gr,
        TxDb,
        gene2symbol,
        links,
        bw_track,
        window
      )
    )

  gene_expr_plot <- plot_gene_expression(TPM, gene, sex, conditions_color)

  pdf(
    file = paste0(
      save_folder, 
      "/", 
      "Footprints_", 
      sex, 
      "_", 
      gene,
      "_", 
      locus$Locus,
      "_peak_tracks.pdf"
      ), 
    width = 7, 
    height = pdf_height
  )


  gridExtra::grid.arrange(
    plot
    # ncol = 2, 
    # widths = c(4,0.5)
  )


  dev.off()

  # Ugly code:
  # If the code ran ok, create a file to let Snakemake know
  is_plot_ok <- "y"
  write.csv(is_plot_ok, file = snakemake@output[["log"]])
  return(plot)
}
