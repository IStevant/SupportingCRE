source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("tidyr")
  library("cowplot")
  library("grid")
  library("ggplot2")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

bindetect_res <- read.table(snakemake@input[["bindetect"]], header = TRUE, sep="\t")

# Load the matrice used to plot sex differential TF-motif enrichment
# matrix
load(snakemake@input[["heatmap_matrice"]])

###########################################
#                                         #
#              Prepare plot               #
#                                         #
###########################################

bindetect_FC <- bindetect_res[,grep("change", colnames(bindetect_res))]
rownames(bindetect_FC) <- bindetect_res$motif_id

stages <- unique(sapply(strsplit(colnames(bindetect_FC), "_"), `[`, 1))

get_FC <- lapply(stages, function(stg){
  FC_tab <- bindetect_FC[,grep(paste0(stg,".*",stg), colnames(bindetect_FC)), drop=FALSE]
  return(FC_tab)
})

bindetect_FC <- do.call(cbind, get_FC)
bindetect_FC$motif <- rownames(bindetect_FC)

bindetect_FC <- bindetect_FC %>%  gather(stage, FC, -motif) 
bindetect_FC$stage <- sapply(strsplit(bindetect_FC$stage, "_"), `[`, 1)
bindetect_FC$FC <- bindetect_FC$FC*-1
bindetect_FC$color <- paste0(
  ifelse(bindetect_FC$FC<0, "XX", "XY"),
  "_",
  bindetect_FC$stage
  )

bindetect_FC$motif <- factor(bindetect_FC$motif, levels=rev(rownames(matrix)))

plot <- ggplot(bindetect_FC, aes(x = FC, y = motif)) +
  # geom_rect(aes(xmin = -0.15, xmax = 0.15, ymin = -Inf, ymax = Inf), fill = '#eeeeee') +
  geom_point(shape=21, size=6, aes(fill=color)) +
  geom_vline(
    xintercept = 0, 
    linetype="dashed", 
    color = "#666666", 
    linewidth=1
  ) +
  scale_fill_manual(values = conditions_color) +
  theme_bw() +
  ylab("TF-binding motifs") +
  xlab("Differential binding score") +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    # aspect.ratio = 1.5,
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################

# As PDF
save_plot(
  snakemake@output[["pdf"]],
  plot,
  base_width = 38,
  base_height = 22,
  units = c("cm"),
  dpi = 300
)

# As PNG
save_plot(
  snakemake@output[["png"]],
  plot,
  base_width = 38,
  base_height = 22,
  units = c("cm"),
  dpi = 300
)
