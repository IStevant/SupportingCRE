source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("ComplexHeatmap")
  library("cowplot")
  library("gtable")
  library("grid")
  library("ggplot2")
  library("clusterProfiler")
  library("JLutils")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

norm_counts <- read.csv(file = snakemake@input[["norm_counts"]], row.names = 1)
load(snakemake@input[["sig_DARs"]])
samplesheet <- read.csv(file = snakemake@input[["samplesheet"]], row.names = 1)

sex <- snakemake@params[["sex"]]
clusters <- snakemake@params[["clusters"]]

names(conditions_color) <- sort(unique(samplesheet$conditions))
sex_colors <- conditions_color[grepl(sex, names(conditions_color))]

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw the heatmap of the z-scores from the regions differentially accessible along developmental stages
#' @param data Full expression matrix (TPM or normalized counts).
#' @param de_feature Vector containing the names of the regions found as differentially accessible.
#' @param colors Vector containing the hexadecimal colours corresponding to each conditions.
#' @param clusters Number of clusters to separate on the heatmap (k for k-mean clustering).
#' @param res_file Name of the output file containing the list of regions per cluster.
#' @return Pheatmap object.
plot_simple_heatmap <- function(data, de_feature, colors, clusters, res_file) {
  matrix_DAR <- data[rownames(data) %in% de_feature, ]

  # Calculate z-scores
  matrix <- t(scale(t(matrix_DAR)))

  set.seed(654)

  # Cluster matrix
  row_dend <- hclust(dist(matrix), method = "ward.D")

  # Calculate the optimal number of clusters
  clusters <- best.cutree(row_dend, min = 4, max = 15)

  clustering <- cutree(row_dend, k = clusters)

  clustering <- as.factor(cutree(row_dend, k = clusters))
  levels(clustering) <- letters[1:clusters]

  if (length(levels(clustering)) == 4){
    if (sex == "XX") {
      clustering <- factor(clustering, levels = c("c", "b", "a", "d"))
    } else if (sex == "XY") {
      clustering <- factor(clustering, levels = c("a", "c", "b", "d"))
    }  
  }

  levels(clustering) <- letters[1:clusters]

  # Save clustering
  write.csv(clustering, file = res_file, quote = FALSE)

  # Limit zscore to |2|
  matrix[matrix > 2] <- 2
  matrix[matrix < (-2)] <- (-2)

  # Prepare top annotation
  conditions <- sapply(strsplit(colnames(matrix), "_"), `[`, 1)
  annotation_col <- data.frame(
    Stages = conditions
  )
  rownames(annotation_col) <- colnames(matrix)

  # Stages palette
  col_stage <- colors
  names(col_stage) <- unique(conditions)

  # peak cluster color palette
  cluster_colors <- as.vector(MetBrewer::met.brewer("Hokusai1", n = clusters))

  # Color palette for the heatmap
  # blue-yellow-red
  cold <- colorRampPalette(c("#4677b7", "#709eca", "#9ac4dd", "#cce1e3", "#fffee8"))
  warm <- colorRampPalette(c("#fffee8", "#fdd4ab", "#fbab70", "#e96e33", "#d83329"))
  BYR <- c(cold(12), warm(12))
  # green-yellow-brown
  cold <- colorRampPalette(c("#138586", "#44aa9b", "#74cfb1", "#bbe7ce", "#fffee8"))
  warm <- colorRampPalette(c("#fffee8", "#f5da9f", "#ebb655", "#d18f43", "#b56832"))
  GYB <- c(cold(12), warm(12))
  # green-yellow-purple
  cold <- colorRampPalette(c("#138586", "#44aa9b", "#74cfb1", "#bbe7ce", "#fffee8"))
  warm <- colorRampPalette(c("#fffee8", "#f3c2d3", "#d981cf", "#b452c1", "#9030b4"))
  GYP <- c(cold(12), warm(12))

  mypalette <- GYP

  # Top annotation (stages)
  stage_anno <- HeatmapAnnotation(
    Stages = anno_block(
      gp = gpar(fill = col_stage, col = 0),
      labels = names(col_stage),
      labels_gp = gpar(col = "white", fontsize = 17, fontface = "bold"),
      height = unit(6.5, "mm")
    )
  )

  # Row annotation (peak clusters)
  cluster_anno <- rowAnnotation(
    Clusters = anno_empty(
      border = FALSE,
      width = unit(10, "mm")
    )
  )

  # Make the heatmap
  ht_list <- Heatmap(
    matrix,
    name = "z-score",
    top_annotation = stage_anno,
    left_annotation = cluster_anno,
    row_title_rot = 0,
    row_split = clustering,
    column_split = conditions,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    cluster_row_slices = FALSE,
    col = mypalette,
    heatmap_legend_param = list(direction = "vertical"),
    row_title = NULL,
    column_title = NULL,
    column_gap = unit(0.4, "mm")
  )

  height <- 9
  width <- 5

  gTree <- grid.grabExpr(
    {
      # Draw the heatmap
      draw(
        ht_list,
        row_title = paste(scales::comma(nrow(matrix)), "differentially accessible regions"),
        row_title_gp = gpar(fontsize = 19),
        # annotation_legend_list=cluster_legend,
        merge_legend = TRUE,
        use_raster = TRUE,
        raster_quality = 5
      )

      # Add peak cluster annotation as extra-rectangles to avoid bad rasterization (faded colors)
      for (i in 1:length(levels(clustering))) {
        decorate_annotation(
          "Clusters",
          slice = i,
          {
            grid.rect(x = 0.9, width = unit(0.7, "mm"), gp = gpar(fill = "black", col = NA), just = "right")
            grid.circle(x = 0.3, r = unit(3.8, "mm"), gp = gpar(fill = "black"))
            grid.text(x = 0.3, levels(clustering)[i], just = "center", gp = gpar(fontsize = 17, col = "white"))
          }
        )
      }
    },
    height = height,
    width = width
  )

  return(gTree)
}

###########################################
#                                         #
#              Draw heatmap               #
#                                         #
###########################################

matrix <- norm_counts[, grep(sex, colnames(norm_counts))]
sex_colors <- conditions_color[grepl(sex, names(conditions_color))]


dynamic_peaks <- plot_simple_heatmap(
  matrix,
  filtered_StageDARs,
  sex_colors,
  clusters,
  snakemake@output[["clusters"]]
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  dynamic_peaks,
  base_width = 14,
  base_height = 20,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  dynamic_peaks,
  base_width = 14,
  base_height = 20,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
