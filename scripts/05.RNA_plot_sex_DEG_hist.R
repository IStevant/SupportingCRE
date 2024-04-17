source("scripts/00.functions.R")
source("scripts/00.color_palettes.R")


load(snakemake@input[['sig_DEGs']])
samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)
stages <- unique(samplesheet$stages)


names(conditions_color) <- sort(unique(samplesheet$conditions))
XX_colors <- conditions_color[grepl("XX" , names(conditions_color))]
XY_colors <- conditions_color[grepl("XY" , names(conditions_color))]

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


# sex_dymorphic_plot_2 <- plot_dymorphic_genes_2(sex_de_genes, XX_colors, XY_colors)

# save_plot(
# 	"../graphs/240306_sex_dymorphic_genes_per_stage_2.pdf", 
# 	sex_dymorphic_plot_2, 
# 	base_width=20,
# 	base_height=8,
# 	units = c("cm"), 
# 	dpi=300, 
# 	bg = "white"
# )

# save_plot(
# 	"../graphs/240306_sex_dymorphic_genes_per_stage_2.png", 
# 	sex_dymorphic_plot_2, 
# 	base_width=20,
# 	base_height=8,
# 	units = c("cm"), 
# 	dpi=300, 
# 	bg = "white"
# )
