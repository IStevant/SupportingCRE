source(".Rprofile")

source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("cowplot")
  library("grid")
  library("ggplot2")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# filtered_StageDARs
XX_filtered_StageDARs <- read.csv(file = snakemake@input[["XX_peak_list"]], header = TRUE, row.names = 1)
colnames(XX_filtered_StageDARs) <- "cluster"
XY_filtered_StageDARs <- read.csv(file = snakemake@input[["XY_peak_list"]], header = TRUE, row.names = 1)
colnames(XY_filtered_StageDARs) <- "cluster"

# Promoter region, i.e. distance to TSS
promoter <- snakemake@params[["promoter"]]

# Load mouse genome
genome_file <- snakemake@input[["genome"]]
Genes <- rtracklayer::import(genome_file)

txdb <- GenomicFeatures::makeTxDbFromGRanges(
  Genes,
  drop.stop.codons = FALSE
)

###########################################
#                                         #
#          ChIPseeker options             #
#                                         #
###########################################

# Ignore unnecessary annotation
options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Generate the read count matrix
#' @param anno annotation list containing outputs from ChIPseeker::annotatePeak.
#' @return Return a dataframe.
plot_anno_sex <- function(anno) {
  anno_XX <- anno[["XX"]]
  anno_XY <- anno[["XY"]]
  clusters <- names(anno_XX)
  anno_anno_XX <- lapply(1:length(anno_XX), function(x) {
    stat <- anno_XX[[x]]@annoStat
    stat$cluster <- rep(clusters[x], nrow(stat))
    return(stat)
  })
  data_XX <- data.table::rbindlist(anno_anno_XX)
  data_XX$Sex <- rep("Pre-gran.", nrow(data_XX))

  anno_anno_XY <- lapply(1:length(anno_XY), function(x) {
    stat <- anno_XY[[x]]@annoStat
    stat$cluster <- rep(clusters[x], nrow(stat))
    return(stat)
  })
  data_XY <- data.table::rbindlist(anno_anno_XY)
  data_XY$Sex <- rep("Sertoli", nrow(data_XY))

  data <- rbind(data_XX, data_XY)

  labels <- round(data$Frequency, digit = 2)
  labels[labels < 4] <- 0
  labels <- paste0(labels, "%")
  labels[labels == "0%"] <- " "

  colors <- c(
    "#f8777c",
    "#0e1b47",
    "#4461a8",
    "#21a3ea",
    "#3bc9d6",
    "#b886da"
  )

  plot <- ggplot(data, aes(fill = Feature, x = cluster, y = Frequency)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = labels, color = Feature), size = 5, position = position_stack(vjust = 0.5)) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = alpha(colors, 0.8)) +
    scale_color_manual(values = c("black", "white", "white", "black", "black", "black"), guide = "none") +
    facet_wrap(~Sex, nrow = 1) +
    theme_light() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 14, margin = margin(r = 10, unit = "pt")),
      strip.text.x = element_text(size = 14, face = "bold"),
      legend.box.spacing = unit(0, "mm"),
    )
  return(plot)
}

###########################################
#                                         #
#           Get peak annotation           #
#                                         #
###########################################

# XX_peaks <- unlist(peak_list[grep("XX", names(peak_list))])
# XY_peaks <- unlist(peak_list[grep("XY", names(peak_list))])

XX_clusters <- unique(XX_filtered_StageDARs$cluster)
XX_peaks <- lapply(
  XX_clusters,
  function(cl){
    GenomicRanges::GRanges(rownames(XX_filtered_StageDARs[XX_filtered_StageDARs$cluster==cl, , drop=FALSE]))
  }
)
names(XX_peaks) <- XX_clusters

# print(head(XX_peaks))

XY_clusters <- unique(XY_filtered_StageDARs$cluster)
XY_peaks <- lapply(
  XY_clusters,
  function(cl){
    GenomicRanges::GRanges(rownames(XY_filtered_StageDARs[XY_filtered_StageDARs$cluster==cl, , drop=FALSE]))
  }
)
names(XY_peaks) <- XY_clusters


XX_anno <- lapply(
  XX_peaks,
  function(cl) {
    ChIPseeker::annotatePeak(
      cl,
      tssRegion = c(-promoter, 0),
      TxDb = txdb,
      overlap = "all"
    )
  }
)

names(XX_anno) <- XX_clusters

XY_anno <- lapply(
  XY_peaks,
  function(cl) {
    ChIPseeker::annotatePeak(
      cl,
      tssRegion = c(-promoter,0),
      TxDb = txdb
    )
  }
)

names(XY_anno) <- XY_clusters

peak_anno_list <- list(
  XX = XX_anno,
  XY = XY_anno
)

# save(peak_anno_list, file = snakemake@output[["anno_list"]])

###########################################
#                                         #
#           Plot peak annotation          #
#                                         #
###########################################

plot_list <- plot_anno_sex(peak_anno_list)

figures <- plot_grid(
  plotlist = list(plot_list)
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  figures,
  base_width = 30,
  base_height = 10,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  figures,
  base_width = 30,
  base_height = 10,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
