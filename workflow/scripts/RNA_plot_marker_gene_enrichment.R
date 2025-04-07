source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("cowplot")
  library("ggplot2")
  library("ComplexHeatmap")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

TPM <- read.csv(file = snakemake@input[["tpm"]], row.names = 1)
whole_gonad <- read.csv(file = snakemake@input[["whole_gonad"]])

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw a dotplot with the enrich$ent and expression level of marketr genes in sorted cells compated to whole gonads
#' @param TPM TPM expression matrix.
#' @param sex Genetic sex of the cells, either "XX" or "XY".
#' @return Ggplot object.
mCherry_dotplot <- function(TPM, sex) {
  if (sex == "XX") {
    genes <- c("Nr5a1", "Wt1", "Gata4", "Foxl2", "Fst", "Lgr5", "Runx1", "Irx3", "Wnt6", "Amhr2", "Bmp2", "Lef1", "Sp5", "Nr2f2", "Tcf21", "Pdgfra", "Wnt5a", "Arx", "Maf", "Stra8", "Mael", "Ddx4", "Dazl", "Figla")
    title <- "Enrichment of gonadal cell marker genes in Enh8-mCherry+ cells\ncompared to whole gonads"
  } else {
    genes <- c("Nr5a1", "Wt1", "Gata4", "Sry", "Gadd45g", "Nr0b1", "Sox9", "Fgf9", "Ptgds", "Amh", "Dmrt1", "Amhr2", "Nr2f2", "Tcf21", "Pdgfra", "Wnt5a", "Arx", "Maf", "Stra8", "Mael", "Ddx4", "Dazl", "Figla")
    title <- "Enrichment of gonadal cell marker genes in Sox9-IRES-GFP+ cells\ncompared to whole gonads"
  }
  expr_sex <- TPM[genes, grep(sex, colnames(TPM))]
  stages <- unique(sapply(strsplit(grep("whole", colnames(expr_sex), value = TRUE), "_"), `[`, 2))
  mean_expr <- lapply(stages, function(stg) {
    stage <- rep(stg, nrow(expr_sex))
    mCherry <- expr_sex[, grep("supporting", colnames(expr_sex))]
    median_expr_mCherry <- rowMeans(as.matrix(mCherry[, grep(stg, colnames(mCherry))]))
    whole <- expr_sex[, grep("whole", colnames(expr_sex))]
    median_expr_Zhao <- rowMeans(as.matrix(whole[, grep(stg, colnames(whole))]))

    ratio_expr <- log2(median_expr_mCherry / median_expr_Zhao)

    ratio_expr[ratio_expr > 1] <- 1
    ratio_expr[ratio_expr < (-1)] <- (-1)

    data <- data.frame(
      stage = stage,
      genes = names(median_expr_mCherry),
      exp = median_expr_mCherry,
      Enrichment = ratio_expr
    )
    return(data)
  })

  data <- do.call(rbind, mean_expr)
  data$genes <- factor(data$genes, levels = unique(data$genes))

  plot <- ggplot(data, aes(stage, genes, size = exp)) +
    geom_point(shape = 21, aes(fill = Enrichment), color = "#666666") +
    scale_fill_gradient2(low = "#04bbc6", mid = "#fffee8", high = "#eb2d62") +
    scale_size(range = c(0, 10)) +
    scale_y_discrete(limits = rev) +
    labs(size = "Expression (TPM)", fill = "Log2 enrichment") +
    ggtitle(title) +
    theme_light() +
    theme(
      strip.text.x = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 14),
      axis.text.y = element_text(face = "italic"),
      axis.title = element_blank(),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 11, margin = margin(r = 10, unit = "pt")),
      legend.position = "right"
    ) 
    # guides(
    #   size = guide_legend(title.position = "top", title.hjust = 0.5),
    #   fill = guide_colourbar(title.position = "top", title.hjust = 0.5)
    # )

  return(plot)
}

##########################################
#                                        #
#              Prepare data              #
#                                        #
##########################################

# Rename RNA-seq data to label if whole gonad or supporting cells
rownames(whole_gonad) <- make.names(whole_gonad$X, unique = TRUE)
whole_gonad <- whole_gonad[, -1]
colnames(whole_gonad) <- paste0("whole.gonad_", colnames(whole_gonad))

colnames(TPM) <- paste0("supporting_", colnames(TPM))

# Merge supporting and whole gonad RNA-seq data
TPM <- merge(TPM, whole_gonad, by = "row.names", all.x = TRUE)
rownames(TPM) <- TPM[, 1]
TPM <- TPM[, -1]

##########################################
#                                        #
#          Plot gene expression          #
#                                        #
##########################################

plot_mCherry_XX <- mCherry_dotplot(TPM, "XX")
plot_mCherry_XY <- mCherry_dotplot(TPM, "XY")

figure <- plot_grid(
  plotlist = list(plot_mCherry_XX, plot_mCherry_XY),
  labels = "AUTO",
  ncol = 2,
  align = "hv"
)

##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################

# As PDF
save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 25,
  base_height = 22,
  units = c("cm"),
  dpi = 300
)

# As PNG
save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 25,
  base_height = 22,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
