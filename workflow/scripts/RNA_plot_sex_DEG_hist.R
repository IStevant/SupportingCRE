source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

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

load(snakemake@input[['sig_DEGs']])
samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)
stages <- unique(samplesheet$stages)


names(conditions_color) <- sort(unique(samplesheet$conditions))
XX_colors <- conditions_color[grepl("XX" , names(conditions_color))]
XY_colors <- conditions_color[grepl("XY" , names(conditions_color))]

#################################################################################################################################

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw the histogram of the sexually dimorphic genes.
#' @param de_genes Dataframe with the number of differentially expressed genes in each sex and for each stage.
#' @param XX_colors Vector of hexadecimal stage colors for the XX samples.
#' @param XY_colors Vector of hexadecimal stage colors for the XY samples.
#' @return Ggplot object.
plot_dymorphic_genes <- function(de_genes, XX_colors, XY_colors){
	axis_margin <- 5.5
	xx_sex_de_genes <- de_genes[de_genes$sex=="XX",]
	p1 <- ggplot(xx_sex_de_genes, aes(x=DE, y=stage)) +
		geom_col(fill=alpha(XX_colors, 0.8), color=XX_colors) +
		geom_text(aes(label = scales::comma(abs(DE))), x=-500, size=4.5, color=c("#333333", rep("white", 3))) +
		scale_x_reverse(labels = scales::comma, limits=c(max(de_genes$DE)+150,0)) +
		scale_y_discrete(position = "right", limits=rev) +
		ggtitle("XX overexpressed genes") +
		xlab("Gene count") +
		theme_light() +
		theme(
			axis.text.y = element_blank(),
			axis.title.y = element_blank(),
			plot.margin = margin(axis_margin, 0, axis_margin, axis_margin),
			plot.title = element_text(size=12, hjust = 0.5, color="#333333"),
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			aspect.ratio=0.5
		)

	xy_sex_de_genes <- de_genes[de_genes$sex=="XY",]
	p2 <- ggplot(xy_sex_de_genes, aes(x=DE, y=stage)) +
		geom_col(fill=alpha(XY_colors, 0.8), color=XY_colors) +
		geom_text(aes(label = scales::comma(abs(DE))), x=500, size=4.5, color=c("#333333", rep("white", 3))) +
		scale_x_continuous(labels = scales::comma, limits=c(0,max(de_genes$DE)+150)) +
		scale_y_discrete(limits=rev) +
		ggtitle("XY overexpressed genes") +
		xlab("Gene count") +
		theme_light() +
		theme(
			axis.title.y = element_blank(),
			plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
			axis.text.y.left = element_text(margin = margin(0, 9, 0, 0)),
			plot.title = element_text(size=12, hjust = 0.5, color="#333333"),
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			aspect.ratio=0.5
		)

	gtitle <- ggdraw() + 
		draw_label(
			"Sexually dimorphic genes per stages",
			fontface = 'bold',
			size = 14,
			hjust = 0.5
		)

	plot_graphs <- plot_grid(p1, p2, align="hv", axis="bt", nrow=1)

	plot <- plot_grid(gtitle, plot_graphs, nrow=2, rel_heights = c(0.1, 1))

	return(plot)

}

#################################################################################################################################

###########################################
#                                         #
#           Histogram sex DEGs            #
#                                         #
###########################################

# Count the number of DEG in each condition
nb_DEGs <- lapply(filtered_SexDEGs, function(stg_Sex_DEG) table(stg_Sex_DEG$Diff.Exp.))

# Prepare the dataframe
sex_de_genes <- data.frame(
	stage=c("E11.5", "E11.5", "E12.5", "E12.5", "E13.5", "E13.5", "E15.5", "E15.5"),
	sex=c("XX", "XY", "XX", "XY", "XX", "XY", "XX", "XY"),
	DE=unlist(nb_DEGs, use.names = FALSE)
)

sex_dymorphic_plot <- plot_dymorphic_genes(sex_de_genes, XX_colors, XY_colors)

save_plot(
	snakemake@output[['pdf']], 
	sex_dymorphic_plot, 
	base_width=20,
	base_height=8,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)

save_plot(
	snakemake@output[['png']], 
	sex_dymorphic_plot, 
	base_width=20,
	base_height=8,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)