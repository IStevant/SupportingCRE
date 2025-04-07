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
#               Load data                 #
#                                         #
###########################################

bw_folder <- snakemake@params[["bw_folder"]]
genome <- snakemake@input[["genome"]]
peaks <- snakemake@input[["peaks"]]
links <- snakemake@input[["linkage"]]
gene_bed <- snakemake@input[["gene_bed"]]
peak_list <- snakemake@input[["peak_list"]]

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Import the coverage data for a given genomic locus from bigwig files.
#' @param folder Path to the bigwig files.
#' @param files Bigwig file names.
#' @param locus GRanges coordinate of the locus of interest.
#' @return Return a list of GRanges objects.
import_bw_files <- function(folder, files, locus) {
  # locus <- resize(gene_TSS[gene_TSS$name %in% locus, ], width = 2, fix = "start")
  full_locus <- locus

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
    gr$cond <- rep(paste("         ", sex, stage), length(gr))
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

#' Create the coverage tracks with replicates overlay.
#' @param gr_list List containing the extended locus coverage per condition.
#' @param locus GRanges object with the coordinate of the locus of interest.
#' @param max_score Maximal score to set the Y axis limit.
#' @param colors vector of hexadecimal colors.
#' @return Return a Gviz object.
generate_genomic_tracks <- function(gr_list, locus, max_score, colors) {
  # locus <- resize(gene_TSS[gene_TSS$name %in% locus, ], width = 2, fix = "start")

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
        window = res,
        chromosome = as.character(seqnames(locus)),
        name = unique(gr$cond),
        col.baseline = color,
        col.histogram = color,
        fill.histogram = color,
        ylim = c(0, trunc(max_score, digit = 4)),
        yTicksAt = c(0, trunc(max_score, digit = 4)),
        rotation.title = 0,
        lwd = 0,
        alpha = 0.8
      )
    }
  )
  track <- OverlayTrack(gTrack_list)
  return(track)
}

#' Plot the genomic tracks with the peaks and the genome annotation.
#' @param peak_gr GRanges object with the coordinate of the ATAC peaks.
#' @param TxDb Genome annotation.
#' @param gene2symbol Data frame containing the conversion of the Ensembl gene IDs to gene symbols.
#' @param link Dataframe with the link coordinates.
#' @param plot_list Gviz objects generated using generate_genomic_tracks().
#' @param locus GRanges coordinate of the locus of interest.
#' @return Return a Gviz plot.
plot_tracks <- function(peak_gr, TxDb, gene2symbol, link, plot_list, locus) {
  # locus_TSS <- resize(gene_TSS[gene_TSS$name %in% locus, ], width = 2, fix = "start")

  full_locus <- locus_TSS

  locus_gr <- full_locus
  peak_gr_locus <- subsetByOverlaps(peak_gr, locus_gr)
  gene_links <- links[links$Gene %in% locus, ]
  print(paste(nrow(gene_links), "interactions found."))

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
    interaction_track <- InteractionTrack(interactions, name = "Links")

    displayPars(interaction_track) <- list(
      plot.anchors = FALSE,
      plot.outside = TRUE,
      col.interaction.types = c("positive-distal" = "#EF6351", "negative-distal" = "#2191FB"),
      rotation.title = 0,
      legend = TRUE
    )
  }

  genome_track <- GenomeAxisTrack(
    col = "#333333",
    fontsize = 13,
    lwd = 1,
    distFromAxis = 0.5,
    labelPos = "alternating",
    sizes = 0.3
  )

  gene_track <- GeneRegionTrack(
    TxDb,
    collapseTranscripts = "meta",
    name = "Genes",
    chromosome = as.character(seqnames(locus_gr)),
    start = start(locus_gr),
    end = end(locus_gr),
    fontface.group = "italic",
    col = NULL,
    col.line = NULL,
    fill = "#585858",
    fontcolor.group = "#333333",
    fontsize.group = 18,
    sizes = 0.3,
    rotation.title = 0,
    thinBoxFeature = "UTR"
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
    col = NULL,
    fill = "#5f4780",
    sizes = 0.2,
    rotation.title = 0
  )


  visible_linked_peaks <- subsetByOverlaps(link_peak_gr, locus_gr)
  promoter <- locus_TSS + 500

  if (length(visible_linked_peaks) >= 1) {
    ht <- HighlightTrack(
      trackList = c(peak_track, gene_track, plot_list),
      # start = c(start(visible_linked_peaks), start(promoter)),
      start = start(visible_linked_peaks),
      # end = c(end(visible_linked_peaks), end(promoter)),
      end = end(visible_linked_peaks),
      col = c(visible_linked_peaks$colors, "#000000"),
      fill = c(visible_linked_peaks$colors, "#000000"),
      chromosome = as.character(unique(seqnames(visible_linked_peaks))),
      alpha = 0.2,
      lwd = 0.5
    )

    title <- paste0(seqnames(locus_gr), ":", start(locus_gr), "-", end(locus_gr))
    
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
      sizes = c(0.3, 0.2, 0.4, rep(0.4, length(plot_list)), 0.3),
      cex.title = 0.9,
      cex.id = 0.5,
      title.width = 1.4,
      main = title,
      cex.main = 1,
      col.main = "#333333"
    )

  } else if (nrow(gene_links) < 1) {
    ht <- HighlightTrack(
      trackList = c(peak_track, gene_track, plot_list),
      start = start(promoter),
      end = end(promoter),
      col = "#000000",
      fill = "#000000",
      chromosome = as.character(seqnames(promoter)),
      alpha = 0.2,
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
      sizes = c(0.3, 0.2, 0.4, rep(0.4, length(plot_list)), 0.3),
      cex.title = 0.9,
      cex.id = 0.5,
      title.width = 1.4
    )
  } else {
    ht <- HighlightTrack(
      trackList = c(peak_track, gene_track, plot_list),
      start = start(promoter),
      end = end(promoter),
      col = "#000000",
      fill = "#000000",
      chromosome = as.character(seqnames(promoter)),
      alpha = 0.2,
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
      sizes = c(0.3, 0.2, 0.4, rep(0.4, length(plot_list)), 0.3),
      cex.title = 0.9,
      cex.id = 0.5,
      title.width = 1.4
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
  names(conditions_color) <- paste0(
    sapply(strsplit(names(conditions_color), "_"), `[`, 2),
    "_",
    sapply(strsplit(names(conditions_color), "_"), `[`, 1)
  )

  print(conditions_color)

  if (sex=="all") {
    exp <- as.numeric(TPM[gene, ])
    sex <- sapply(strsplit(colnames(TPM), "_"), `[`, 2)
  } else {
    TPM <- TPM[,grep(sex, colnames(TPM))]
    exp <- as.numeric(TPM[gene, ])
    sex <- sapply(strsplit(colnames(TPM), "_"), `[`, 2)
  }

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
      ggtitle(paste0("\n", gene, "\n")) +
      xlab("TPM\n") +
      scale_x_continuous(labels = scales::comma, n.breaks = 3) +
      scale_y_discrete(limits=rev) +
      coord_cartesian(ylim = c(0, NA)) +
      theme_light() +
      theme(
        plot.title = element_text(size = 12, face = "bold.italic", hjust = 0.5),
        axis.text = element_text(size = 10),
        # axis.title.x = element_blank(),
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

bw_files <- list.files(path = bw_folder, pattern = "REP..bw")
links <- read.table(links, header = TRUE)
peak_list <- read.table(peak_list, header = TRUE)
genes <- peak_list$Gene
sex <-  peak_list$Sex

peak_gr <- GRanges(rownames(read.csv(peaks, header = TRUE, row.names = 1)))
gene_TSS <- GRanges(rtracklayer::import(gene_bed))

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
rownames(gene2symbol) <- gene2symbol$gene_id

TxDb <- GenomicFeatures::makeTxDbFromGFF(genome)


# pdf(file = snakemake@output[["pdf"]], width = 8, height = 8)

peak2plot <- peak_list$Locus


plot <- foreach(peak = peak2plot) %dopar% {
# plot <- lapply(genes, function(gene) {

  sex <- peak_list[peak_list$Locus == peak, "Sex"]

  if (sex != "all") {
    bw_files <- grep(sex, bw_files, value = TRUE)
    conditions_color <- conditions_color[grep(sex, names(conditions_color))]
  }

  locus <- peak_list[peak_list$Locus == peak, ]
  max_score <- locus$MaxScore

  atac_coverage <- import_bw_files(
    bw_folder,
    bw_files,
    gene
  )

  if (max_score == 0) {
    max_score <- do.call("c", do.call("c", do.call("c", atac_coverage)))
    max_score <- max(unlist(lapply(max_score, function(x) x$score)))
  }

  if (sex == "all") {
    XX_track_list <- lapply(
      atac_coverage[["XX"]],
      function(stage) {
        generate_genomic_tracks(
          stage,
          gene,
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
          max_score,
          conditions_color
        )
      }
    )

    bw_track <- c(XX_track_list, XY_track_list)
    pdf_height <- 6.5

  } else {

    plot <- grid::grid.grabExpr(
      plot_tracks(
        atac_coverage[[sex]],
        TxDb,
        gene2symbol,
        links,
        c(XX_track_list, XY_track_list),
        gene,
        window
      )
    )

  pdf_height <- 4

  }

  gene_expr_plot <- plot_gene_expression(TPM, gene, sex, conditions_color)

  pdf(
    file = paste0(
      save_folder, 
      "/", 
      "ATAC_", 
      sex, 
      "_", 
      gene, 
      "_peak_tracks.pdf"
      ), 
    width = 8, 
    height = pdf_height
  )

  gridExtra::grid.arrange(
    plot, gene_expr_plot, 
    ncol = 2, 
    widths = c(4,0.5)
  )

  dev.off()

  # Ugly code:
  # If the code ran ok, create a file to let Snakemake know
  is_plot_ok <- "y"
  write.csv(is_plot_ok, file = snakemake@output[["log"]])

  return(plot)
}

