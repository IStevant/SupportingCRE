source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("ggrepel")
  library("ggplot2")
  library("doParallel")
  library("foreach")
  library("clusterProfiler")
  library("org.Mm.eg.db")
  library("cowplot")
})



###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

load(snakemake@input[["all_DEGs"]])
samplesheet <- read.csv(file = snakemake@input[["samplesheet"]], row.names = 1)
stages <- unique(samplesheet$stages)
names(conditions_color) <- sort(unique(samplesheet$conditions))

adj.pval <- snakemake@params[["adjpval"]]
log2FC <- snakemake@params[["log2FC"]]
path <- snakemake@params[["path"]]

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 12))
doParallel::registerDoParallel(cores = n_cores)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Prepare the volcano plot.
#' @param dds DESeq2 analysis result object.
#' @param stage Embryonic stage.
#' @param p.adj adjusted p-value threshold used in the DE analysis.
#' @param log2FC Log2(Fold change) threshold used in the DE analysis.
#' @param colors Vector of hexadecimal stage colors for the current embryonic stage.
#' @return Ggplot object.
plot_volcano_sex <- function(dds, stage, p.adj, log2FC, colors) {
  colors <- conditions_color[grep(stage, names(conditions_color))]
  res <- DESeq2::results(dds, contrast = c("conditions", paste("XX", stage), paste("XY", stage)))
  de_res <- as.data.frame(res)

  res <- mutate(
    de_res,
    Diff.Exp. = dplyr::case_when(
      log2FoldChange >= log2FC & padj <= p.adj ~ "Up in XX",
      log2FoldChange <= (-log2FC) & padj <= p.adj ~ "Up in XY",
      TRUE ~ "non sig."
    )
  )

  res$log10_pval_noInf <- -log10(res$padj)
  res$log10_pval_noInf[res$log10_pval_noInf == Inf] <- 400

  sig.DE <- subset(res, padj < p.adj)
  sig.DE <- subset(sig.DE, abs(log2FoldChange) > log2FC)

  nb_sig_DE <- scales::comma(nrow(sig.DE))

  pval_order <- res[with(res, order(log10_pval_noInf, decreasing = TRUE)), ]

  top_p_up <- rownames(head(pval_order[pval_order$Diff.Exp. == "Up in XX", ], 10))
  top_p_down <- rownames(head(pval_order[pval_order$Diff.Exp. == "Up in XY", ], 10))
  top_genes <- c(top_p_up, top_p_down)
  top_genes <- c("Sox9", "Amh", "Dhh", "Fgf9", "Ptgds", "Sry", "Nr0b1", "Dmrt1", "Wt1", "Nr5a1", "Foxl2", "Fst", "Wnt4", "Lgr5", "Ctnnb1", "Bmp2", "Runx1", "Cdkn1a")
  sig_genes <- sig.DE[rownames(sig.DE) %in% top_genes, ]

  up_in_cond1 <- nrow(pval_order[pval_order$Diff.Exp. == "Up in XX", ])
  up_in_cond2 <- nrow(pval_order[pval_order$Diff.Exp. == "Up in XY", ])

  de_res <- na.omit(res)
  max_FC <- max(abs(res$log2FoldChange))

  volcano <- ggplot(res, aes(x = -log2FoldChange, y = log10_pval_noInf)) +
    geom_point(
      aes(
        colour = Diff.Exp.,
        fill = Diff.Exp.
      ),
      alpha = 0.5,
      shape = 21,
      size = 2
    ) +
    geom_point(
      data = sig_genes,
      shape = 21,
      size = 3,
      colour = "black"
    ) +
    geom_hline(
      yintercept = -log10(p.adj),
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c(-log2FC, log2FC),
      linetype = "dashed"
    ) +
    geom_label_repel(
      data = sig_genes,
      mapping = aes(
        label = rownames(sig_genes)
      ),
      color = "black",
      size = 4,
      min.segment.length = unit(0, "lines"),
      fill = alpha("white", 0.5),
      fontface = "italic",
      force = 3,
      max.overlaps = 25,
    ) +
    expand_limits(
      x = c(-25, 25),
      y = c(0, 450)
    ) +
    ggtitle(paste0(stage, " XX vs ", stage, " XY \n", nb_sig_DE, " significant DE genes")) +
    theme_bw() +
    xlab("log2(fold change)") +
    ylab("-log10(adjusted p-values)") +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, vjust = 0.5, hjust = 0.5),
      legend.position = "bottom",
      aspect.ratio = 0.8
    ) +
    scale_color_manual(
      values = c("#5D576B", as.vector(colors)),
      labels = c("non-sig.", paste0("Up in XX\n(", scales::comma(up_in_cond1), " genes)"), paste0("Up in XY\n(", scales::comma(up_in_cond2), " genes)"))
    ) +
    scale_fill_manual(
      values = c("#5D576B", as.vector(colors)),
      labels = c("non-sig.", paste0("Up in XX\n(", scales::comma(up_in_cond1), " genes)"), paste0("Up in XY\n(", scales::comma(up_in_cond2), " genes)"))
    )
  return(volcano)
}


#' Prepare the volcano plot.
#' @param dds DESeq2 analysis result object.
#' @param stage Embryonic stage.
#' @param p.adj adjusted p-value threshold used in the DE analysis.
#' @param log2FC Log2(Fold change) threshold used in the DE analysis.
#' @param colors Vector of hexadecimal stage colors for the current embryonic stage.
#' @param path Path to the GO term result table file.
#' @return Ggplot object.
plot_DEG <- function(dds, stage, p.adj, log2FC, colors, path) {
  volcano <- plot_volcano_sex(dds, stage, p.adj, log2FC, colors)

  res <- DESeq2::results(dds, contrast = c("conditions", paste("XX", stage), paste("XY", stage)))

  de_res <- as.data.frame(res)
  res <- mutate(
    de_res,
    Diff.Exp. = dplyr::case_when(
      log2FoldChange >= log2FC & padj <= p.adj ~ "Up in XX",
      log2FoldChange <= (-log2FC) & padj <= p.adj ~ "Up in XY",
      TRUE ~ "non sig."
    )
  )
  sig.DE <- subset(res, padj < p.adj)
  sig.DE <- subset(sig.DE, abs(log2FoldChange) > log2FC)

  de_genes <- data.frame(
    gene = rownames(sig.DE),
    cluster = sig.DE$Diff.Exp.
  )

  res_file <- paste0(path, "/RNA_", stage, "_GO_sex_DEG.csv")
  GO_terms <- GO_term_per_cluster(de_genes, res_file)
  go_term_plot <- go_plot(GO_terms, nb_terms = 5)
  # KEGG_terms <- KEGG_term_per_cluster(de_genes)
  # if (!is.null(KEGG_terms)){
  # 	KEGG_term_plot <- kegg_plot(KEGG_terms, nb_terms=5)
  # } else {
  # 	KEGG_term_plot <- NULL
  # }

  figure <- plot_grid(plotlist = list(volcano, go_term_plot), ncol = 2)

  return(figure)
}

#' Get GO term enrichment.
#' @param de_genes Vector of gene names.
#' @param res_file Path and name of the result file.
#' @return dataframe.
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
    cutoff = 0.6,
    by = "p.adjust",
    select_fun = min
  )

  write.table(lineage1_ego, file = res_file, quote = FALSE, sep = "\t", row.names = FALSE)

  return(lineage1_ego)
}

#' Get KEGG pathway enrichment.
#' @param de_genes Vector of gene names.
#' @return dataframe.
KEGG_term_per_cluster <- function(de_genes) {
  print("Calculate KEGG term over-representation...")
  entrez_genes <- bitr(de_genes$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  de_gene_clusters <- de_genes[de_genes$gene %in% entrez_genes$SYMBOL, c("gene", "cluster")]
  de_gene_clusters <- data.frame(
    ENTREZID = entrez_genes$ENTREZID[which(entrez_genes$SYMBOL == de_gene_clusters$gene)],
    cluster = de_gene_clusters$cluster
  )

  formula_res <- compareCluster(
    ENTREZID ~ cluster,
    data = de_gene_clusters,
    fun = "enrichKEGG",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
  )
  return(formula_res)
}

#' Plot GO term enrichment.
#' @param go_res Dataframe containing the GO term enrichment result.
#' @param nb_terms Number of top GO term per cluster to show.
#' @return Ggplot object
go_plot <- function(go_res, nb_terms = 5) {
  # Make the first GO term letter as capital letter
  go_res@compareClusterResult[, 4] <- gsub("^([a-z])", "\\U\\1", go_res[, 4], perl = TRUE)
  options(enrichplot.colours = c("#77BFA3", "#98C9A3", "#BFD8BD", "#DDE7C7", "#EDEEC9"))
  plot <- clusterProfiler::dotplot(go_res, showCategory = nb_terms)
  x_labels <- levels(as.data.frame(go_res)$Cluster)
  # Important: Reload ggplot2 to apply new theme
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

#' Plot GO term enrichment.
#' @param kegg_res Dataframe containing the KEGG enrichment result.
#' @return Ggplot object
kegg_plot <- function(kegg_res, nb_terms = 5) {
  plot <- dotplot(kegg_res)
  return(plot)
}
#################################################################################################################################

###########################################
#                                         #
#         Volcano + GO sex DEGs           #
#                                         #
###########################################
sex_DEG_plots <- foreach(stg = stages) %dopar% {
  plot_DEG(SexDEGs, stg, adj.pval, log2FC, conditions_color, path)
}

figure <- plot_grid(
  plotlist = sex_DEG_plots,
  labels = "AUTO",
  ncol = 1
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 28,
  base_height = 48,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 28,
  base_height = 48,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
