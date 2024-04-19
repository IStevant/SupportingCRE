source("scripts/00.color_palettes.R")

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
#               Functions                 #
#                                         #
###########################################

#' Draw the histogram of the sexually dimorphic accessible regions.
#' @param dar Dataframe with the number of differentially accessible peaks in each sex and for each stage.
#' @param XX_colors Vector of hexadecimal stage colors for the XX samples.
#' @param XY_colors Vector of hexadecimal stage colors for the XY samples.
#' @return Ggplot object.
plot_dymorphic_genes <- function(dar, XX_colors, XY_colors){
	axis_margin <- 5.5
	xx_sex_dar <- dar[dar$sex=="XX",]
	p1 <- ggplot(xx_sex_dar, aes(x=DA, y=stage)) +
		geom_col(fill=XX_colors) +
		geom_text(aes(label = scales::comma(abs(DA))), x=-5000, size=4.5, color="white") +
		scale_x_reverse(labels = scales::comma, limits=c(max(dar$DA)+150,0)) +
		scale_y_discrete(position = "right", limits=rev) +
		ggtitle("Regions more accessible in XX") +
		xlab("Region count") +
		theme_light() +
		theme(
			axis.text.y = element_blank(),
			axis.title.y = element_blank(),
			plot.margin = margin(axis_margin, 0, axis_margin, axis_margin),
			plot.title = element_text(size=12, hjust = 0.5, face="bold"),
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			aspect.ratio=0.5
		)

	xy_sex_dar <- dar[dar$sex=="XY",]
	p2 <- ggplot(xy_sex_dar, aes(x=DA, y=stage)) +
		geom_col(fill=XY_colors) +
		geom_text(aes(label = scales::comma(abs(DA))), x=5000, size=4.5, color="white") +
		scale_x_continuous(labels = scales::comma, limits=c(0,max(dar$DA)+150)) +
		scale_y_discrete(limits=rev) +
		ggtitle("Regions more accessible in XY") +
		xlab("Region count") +
		theme_light() +
		theme(
			axis.title.y = element_blank(),
			plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
			axis.text.y.left = element_text(margin = margin(0, 9, 0, 0)),
			plot.title = element_text(size=12, hjust = 0.5, face="bold"),
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			aspect.ratio=0.5
		)

	gtitle <- ggdraw() + 
		draw_label(
			"Sexually dimorphic open chromatin regions",
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
#               Load data                 #
#                                         #
###########################################

load(snakemake@input[['sig_DARs']])
samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)
stages <- unique(samplesheet$stages)


names(conditions_color) <- sort(unique(samplesheet$conditions))
XX_colors <- conditions_color[grepl("XX" , names(conditions_color))]
XY_colors <- conditions_color[grepl("XY" , names(conditions_color))]

###########################################
#                                         #
#           Histogram sex DARs            #
#                                         #
###########################################

# Count the number of DAR in each condition
nb_DARs <- lapply(filtered_SexDARs, function(stg_Sex_DAR) table(stg_Sex_DAR$Diff.Acc.))

# Prepare the dataframe
sex_dar <- data.frame(
	stage=c("E11.5", "E11.5", "E12.5", "E12.5", "E13.5", "E13.5", "E15.5", "E15.5"),
	sex=c("XX", "XY", "XX", "XY", "XX", "XY", "XX", "XY"),
	DA=unlist(nb_DARs, use.names = FALSE)
)

sex_dymorphic_plot <- plot_dymorphic_genes(sex_dar, XX_colors, XY_colors)

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