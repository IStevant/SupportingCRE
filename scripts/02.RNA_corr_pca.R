source("scripts/00.functions.R")
source("scripts/00.color_palettes.R")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

TPM <- read.csv(file=snakemake@input[['tpm']], row.names=1)
# print(head(TPM))

conditions <- paste(
	sapply(strsplit(colnames(TPM), "_"), `[`, 2), 
	sapply(strsplit(colnames(TPM), "_"), `[`, 1), 
	sep=" "
)

names(conditions_color) <- sort(unique(conditions))
# print(conditions)

corr_method <- snakemake@params[['corr_method']]

###########################################
#                                         #
#        Plot correlation and PCA         #
#                                         #
###########################################

corr_plot <- grid::grid.grabExpr(correlation(TPM, conditions, corr_method, conditions_color))

pca_plot <- plot.pca(TPM, conditions, conditions_color, c("PC1", "PC2"))

figure <- cowplot::plot_grid(
	plotlist=list(corr_plot, pca_plot[[1]]),
	labels = "AUTO",
	ncol=2
)

cowplot::save_plot(
	snakemake@output[['pdf']], 
	figure, 
	base_width=35,
	base_height=13.5,
	units = c("cm"), 
	dpi=300
)

cowplot::save_plot(
	snakemake@output[['png']], 
	figure, 
	base_width=35,
	base_height=13.5,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)