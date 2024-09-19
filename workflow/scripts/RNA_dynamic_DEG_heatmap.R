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

TF_genes <- snakemake@input[["TF_genes"]]
TF_pheno <- snakemake@input[["TF_pheno"]]

norm_counts <- read.csv(file = snakemake@input[["norm_counts"]], row.names = 1)
load(snakemake@input[["sig_DEGs"]])
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

#' Draw the heatmap of the z-scores from the genes differentially expressed along developmental stages
#' @param data Full expression matrix (TPM or normalized counts).
#' @param de_feature Vector containing the names of the genes found as differentially expressed.
#' @param colors Vector containing the hexadecimal colours corresponding to each conditions.
#' @param clusters Number of clusters to separate on the heatmap (k for k-mean clustering).
#' @param res_file Name of the output file containing the list of genes per cluster.
#' @return Pheatmap object.
draw_heatmap <- function(data, de_feature, colors, clusters, res_file) {
  matrix_DEG <- data[rownames(data) %in% de_feature, ]
  # Calculate z-scores
  matrix <- t(scale(t(matrix_DEG)))

  # Cluster matrix
  row_dend <- hclust(dist(matrix), method = "ward.D2")

  # Calculate the optimal number of clusters
  clusters <- best.cutree(row_dend, min = 7, max = 15)

  clustering <- cutree(row_dend, k = clusters)

  # Order clusters
  mean_per_cluster <- lapply(
    unique(clustering),
    function(cluster) {
      genes <- names(clustering[clustering == cluster])
      matrix_cluster <- matrix[rownames(matrix) %in% genes, ]
      mean <- as.data.frame(t(colMeans(matrix_cluster)))
      return(mean)
    }
  )
  mean_per_cluster <- data.table::rbindlist(mean_per_cluster)
  o1 <- seriation::seriate(dist(mean_per_cluster), method = "GW")
  levels <- seriation::get_order(o1)
  clustering <- factor(clustering, levels = levels)
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
  # Gene cluster color palette
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

  mypalette <- BYR
  # Top annotation (stages)
  stage_anno <- HeatmapAnnotation(
    Stages = anno_block(
      gp = gpar(fill = col_stage, col = 0),
      labels = names(col_stage),
      labels_gp = gpar(col = "white", fontsize = 17, fontface = "bold"),
      height = unit(6.5, "mm")
    )
  )
  # Row annotation (gene clusters)
  cluster_anno <- rowAnnotation(
    Clusters = anno_empty(
      border = FALSE,
      width = unit(10, "mm")
    )
  )
  # Prepare the gene clusters legend manually
  cluster_legend <- Legend(
    at = levels(clustering),
    title = "Clusters",
    legend_gp = gpar(fill = cluster_colors)
  )

  TF_list <- read.csv(TF_genes, header = FALSE)
  gonad_pheno_genes <- read.table(TF_pheno, header = FALSE, sep = "\t")
  TF_list <- as.vector(TF_list[, 1])
  total_TFs <- unlist(lapply(TF_list, function(TF) which(rownames(matrix) %in% TF)))

  gonad_pheno_genes <- as.vector(gonad_pheno_genes[, 1])
  TFs <- TF_list[TF_list %in% gonad_pheno_genes]
  matrix_TF_indexes <- unlist(lapply(TFs, function(TF) which(rownames(matrix) %in% TF)))

  if (length(matrix_TF_indexes) > 25) {
    split_indexes <- split(matrix_TF_indexes, cut(seq_along(matrix_TF_indexes), 25, labels = FALSE))
    set.seed(123)
    matrix_TF_indexes <- unlist(lapply(split_indexes, function(x) sample(x, 1)))
  }
  TF_names <- rownames(matrix)[matrix_TF_indexes]
  # Show transcription factors
  TFs <- rowAnnotation(
    TFs = anno_mark(
      at = matrix_TF_indexes, # TF row indexes
      labels = TF_names, # Gene names
      labels_gp = gpar(fontsize = 17, fontface = "italic"),
      padding = unit(2, "mm")
    )
  )
  # Make the heatmap
  ht_list <- Heatmap(
    matrix,
    name = "z-score",
    top_annotation = stage_anno,
    left_annotation = cluster_anno,
    right_annotation = TFs,
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
  width <- 7

  gTree <- grid.grabExpr(
    {
      # Draw the heatmap
      draw(
        ht_list,
        row_title = paste(scales::comma(nrow(matrix)), "differentially expressed genes"),
        row_title_gp = gpar(fontsize = 19),
        merge_legend = TRUE,
        use_raster = TRUE,
        raster_quality = 5
      )

      # Add gene cluster annotation as a black line, white letters over a black circle
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

  p_fix <- panel_fix(plot_grid(gTree), width = width, height = height, units = "inch")
  p <- plot_grid(p_fix)
  return(p_fix)
}

#' Get enriched Biological Process GO terms and simplify them to avoid term redundancy.
#' @param de_genes Vector containing the names of the genes found as differentially expressed.
#' @param res_file Name of the output file containing the list of genes per cluster.
#' @return Pheatmap object.
GO_term_per_cluster <- function(de_genes, res_file) {
  print("Calculate GO term over-representation...")
  formula_res <- compareCluster(
    gene ~ cluster,
    data = de_genes,
    fun = "enrichGO",
    keyType = "SYMBOL",
    OrgDb = "org.Mm.eg.db",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  print("Calculate GO term semantic similarities...")
  lineage1_ego <- simplify(
    formula_res,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )

  write.table(lineage1_ego, file = res_file, quote = FALSE, sep = "\t", row.names = FALSE)

  # # KEGG pathway enrichment
  # entrez_genes <- bitr(de_genes$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # de_gene_clusters <- de_genes[de_genes$gene %in% entrez_genes$SYMBOL,c("gene", "cluster")]
  # de_gene_clusters <- data.frame(
  # 	ENTREZID=entrez_genes$ENTREZID[which(entrez_genes$SYMBOL==de_gene_clusters$gene)],
  # 	cluster=de_gene_clusters$cluster
  # )

  # formula_res <- compareCluster(
  # 	ENTREZID~cluster,
  # 	data=de_gene_clusters,
  # 	fun="enrichKEGG",
  # 	pAdjustMethod = "BH",
  # 	pvalueCutoff  = 0.01,
  # 	qvalueCutoff  = 0.05
  # )

  return(lineage1_ego)
}

#' Plot the top enriched GO terms.
#' @param go_res GO term enrichment results dataframe from GO_term_per_cluster().
#' @param nb_terms Number of top GO terms to plot per condition.
#' @return Pheatmap object.
go_plot <- function(go_res, nb_terms = 5) {
  # Make the first GO term letter as capital letter
  go_res@compareClusterResult[, 4] <- gsub("^([a-z])", "\\U\\1", go_res[, 4], perl = TRUE)
  options(enrichplot.colours = c("#77BFA3", "#98C9A3", "#BFD8BD", "#DDE7C7", "#EDEEC9"))
  plot <- clusterProfiler::dotplot(go_res, showCategory = nb_terms)
  x_labels <- unique(as.data.frame(go_res)$Cluster)
  # Reload ggplot2 to apply new theme
  library(ggplot2)
  plot <- plot +
    geom_point(shape = 21) +
    ggtitle(paste0("Enriched biological process GO terms (Top ", nb_terms, ")")) +
    scale_x_discrete(labels = x_labels) +
    labs(color = "Adj. p-value", size = "Gene ratio") +
    guides(color = guide_colorbar(reverse = TRUE), size = guide_legend(reverse = TRUE)) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.title.x = element_blank()
    )
  return(plot)
}


###############################################
# Fix complexeheatmap gene annotation

panel_fix <- function(p = NULL, grob = NULL,
                      width = NULL, height = NULL, margin = 1, units = "cm",
                      filename = NULL) {
  if (is.null(p) & is.null(grob)) {
    stop("'p' or 'grob' must be provided with at least one.")
  }
  if (is.null(width) & is.null(height)) {
    stop("'width' or 'height' must be provided with at least one.")
  }
  if (is.null(grob)) {
    grob <- ggplotGrob(p)
  }

  panels <- grep("panel", grob[["layout"]][["name"]])
  panel_index_h <- sort(unique(grob[["layout"]][["t"]][panels]))
  panel_index_w <- sort(unique(grob[["layout"]][["l"]][panels]))
  nw <- length(panel_index_w)
  nh <- length(panel_index_h)
  raw_w <- as.numeric(grob[["widths"]][panel_index_w])
  raw_h <- as.numeric(grob[["heights"]][panel_index_h])
  raw_aspect <- raw_h / raw_w
  if (is.null(width)) {
    width <- height / raw_aspect
  }
  if (is.null(height)) {
    height <- width * raw_aspect
  }
  if (!length(width) %in% c(1, length(raw_aspect)) | !length(height) %in% c(1, length(raw_aspect))) {
    stop("The length of 'width' and 'height' must be 1 or the length of panels.")
  }
  if (length(width) == 1) {
    width <- rep(width, nw)
  }
  if (length(height) == 1) {
    height <- rep(height, nh)
  }
  grob[["widths"]][panel_index_w] <- unit(width, units = units)
  grob[["heights"]][panel_index_h] <- unit(height, units = units)
  grob <- gtable_add_padding(grob, unit(margin, units = units))
  plot_width <- convertWidth(sum(grob[["widths"]]), unitTo = units, valueOnly = TRUE)
  plot_height <- convertHeight(sum(grob[["heights"]]), unitTo = units, valueOnly = TRUE)

  if (!is.null(filename)) {
    ggsave(filename = filename, plot = grob, units = units, width = plot_width, height = plot_height)
  }
  attr(grob, "size") <- list(units = units, width = plot_width, height = plot_height)
  return(grob)
}


###########################################
#                                         #
#              Drax heatmap               #
#                                         #
###########################################
matrix <- norm_counts[, grep(sex, colnames(norm_counts))]
sex_colors <- conditions_color[grepl(sex, names(conditions_color))]


dynamic_genes <- draw_heatmap(
  matrix,
  filtered_StageDEGs,
  sex_colors,
  clusters,
  snakemake@output[["cluster_file"]]
)

de_genes <- read.csv(snakemake@output[["cluster_file"]])
colnames(de_genes) <- c("gene", "cluster")

go_analysis <- GO_term_per_cluster(
  de_genes,
  snakemake@output[["GO"]]
)

go_plot <- go_plot(go_analysis, nb_terms = 4)

go_plot <- go_plot +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

figure <- plot_grid(
  dynamic_genes, go_plot,
  labels = "AUTO",
  ncol = 2
)

##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################

save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 39,
  base_height = 23.2,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 39,
  base_height = 23.2,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
