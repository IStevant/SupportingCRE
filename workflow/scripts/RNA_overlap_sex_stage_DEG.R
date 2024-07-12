source(".Rprofile")
###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
	library("cowplot")
})


###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

draw_venn_sex_dyn <- function(sexID, sex, dyn, title) {
	data <- c(
		sex=length(setdiff(sex,dyn)),
		dyn=length(setdiff(dyn,sex)),
		"sex&dyn"=length(intersect(sex,dyn))
	)

	if (sexID=="XX"){
		colours <- c( "#FCBF49", "#D62828")
	} else {
		colours <- c( "#94d574", "#1e8bd1")
	}

	venn <- eulerr::euler(data)
	title <- title
	venn_plot <- plot(
		venn,
		labels=c(
			paste0("Sexually dymorphic\n(", prettyNum(length(sex), big.mark = ","), ")"),
			paste0("Dynamic\n(", prettyNum(length(dyn), big.mark = ","), ")")
		),
		quantities = list(fontsize=10),
		edges = list(
			col = as.vector(colours), 
			lex = 2
		),
		fills = list(
			fill=as.vector(colours),
			alpha=0.35
		)
	)
}


#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# filtered_SexDEGs
load(snakemake@input[['sex_DEGs']])

# filtered_StageDEGs
load(snakemake@input[['XY_stage_DEGs']])
XY_filtered_StageDEGs <- filtered_StageDEGs

# filtered_StageDEGs
load(snakemake@input[['XX_stage_DEGs']])
XX_filtered_StageDEGs <- filtered_StageDEGs

###########################################
#                                         #
#               Draw Venn                 #
#                                         #
###########################################
XX_spe_genes <- unique(unlist(lapply(filtered_SexDEGs, function(x) rownames(x[x$Diff.Exp. == "Up in XX",]))))

XX_dyn_genes <- XX_filtered_StageDEGs
XX_venn <- draw_venn_sex_dyn(sexID="XX", XX_spe_genes, XX_dyn_genes, "XX genes")

XY_spe_genes <- unique(unlist(lapply(filtered_SexDEGs, function(x) rownames(x[x$Diff.Exp. == "Up in XY",]))))
XY_dyn_genes <- XY_filtered_StageDEGs
XY_venn <- draw_venn_sex_dyn(sexID="XY", XY_spe_genes, XY_dyn_genes, "XY genes")


figure <- plot_grid(
	XX_venn, XY_venn,
	labels = c("XX", "XY"),
	ncol=1
)

save_plot(
	snakemake@output[['pdf']], 
	figure, 
	base_width=12,
	base_height=15,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)

save_plot(
	snakemake@output[['png']], 
	figure, 
	base_width=12,
	base_height=15,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)

