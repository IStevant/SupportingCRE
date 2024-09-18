source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("cowplot")
  library("ggplot2")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

TPM <- read.csv(file = snakemake@input[["tpm"]], row.names = 1)
exampleGenes <- as.vector(read.csv(file = snakemake@input[["genes"]], header=FALSE)[,1])

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw the expression plots for each example gene of a given cell type
#' @param genes Vector containing the example gene names.
#' @param TPM TPM expression matrix.
#' @param title Cell type name.
#' @return Ggplot object.
gene_exp <- function(genes, TPM, title) {
  plotlist <- list()
  for (gene in genes) {
    exp <- as.numeric(TPM[gene, ])
    gene_exp <- data.frame(
      sex = sapply(strsplit(colnames(TPM), "_"), `[`, 2),
      stages = sapply(strsplit(colnames(TPM), "_"), `[`, 1),
      exp = exp
    )
    df <- dplyr::group_by(gene_exp, stages, sex)
    options(dplyr.summarise.inform = FALSE)
    df.summary2 <- dplyr::summarise(
      df,
      sd = sd(exp),
      len = mean(exp)
    )
    plot <- ggplot(df.summary2, aes(x = stages, y = len, group = sex, color = sex)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      geom_errorbar(
        aes(
          ymin = len - sd,
          ymax = len + sd
        ),
        width = .1
      ) +
      scale_color_manual(
        values = c("#FFB100", "#339989"),
        labels = c("Pre-gran.", "Sertoli")
      ) +
      scale_y_continuous(labels = scales::comma) +
      coord_cartesian(ylim = c(0, NA)) +
      labs(title = gene, x = "Embryonic stages", y = "Expression (TPM)") +
      theme_light() +
      theme(
        plot.title = element_text(size = 13, face = "bold.italic", hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        aspect.ratio = 0.5
      )
    legend <- get_legend(
      plot + theme(legend.direction = "horizontal", legend.box.margin = margin(0, 0, 0, 0))
    )
    plot <- plot + theme(legend.position = "bottom")
    plotlist[[gene]] <- plot
  }

  return(plotlist)
}

##########################################
#                                        #
#          Plot gene expression          #
#                                        #
##########################################

# Plot expression of the example genes
plot_genes <- gene_exp(exampleGenes, TPM)

figure <- plot_grid(
  plotlist = plot_genes,
  labels = "AUTO",
  ncol = 5,
  nrow=5,
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
  base_width = 38,
  base_height = 30,
  units = c("cm"),
  dpi = 300
)

# As PNG
save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 38,
  base_height = 30,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
