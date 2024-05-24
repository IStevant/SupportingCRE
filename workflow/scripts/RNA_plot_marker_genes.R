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

TPM <- read.csv(file=snakemake@input[['tpm']], row.names=1)
markerGenes <- read.csv(file=snakemake@input[['marker_genes']])
whole_gonad <- read.csv(file=snakemake@input[['whole_gonad']])

#################################################################################################################################

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw the expression plots marker genes for each cell type
#' @param TPM TPM expression matrix.
#' @param markers Dataframe containing the marker genes and their corresponding cell types.
#' @return Ggplot object.
plot_expression_scatter <- function(TPM, markers){
	cellTypes <- unique(markers$Cell_type)
	plot_list <- lapply(cellTypes, function(cellType){
		genes <- markers[markers$Cell_type %in% cellType, "Gene"]
		gene_exp(genes, TPM, cellType)
	})
}

#' Draw the expression plots for each marker gene of a given cell type
#' @param genes Vector containing the marker gene names.
#' @param TPM TPM expression matrix.
#' @param title Cell type name.
#' @return Ggplot object.
gene_exp <- function(genes, TPM, title){
	plotlist <- list()
	for (gene in genes){
		exp <- as.numeric(TPM[gene, ])
		gene_exp <- data.frame(
			cells=sapply(strsplit(colnames(TPM), "_"), `[`, 1),
			sex=sapply(strsplit(colnames(TPM), "_"), `[`, 3),
			stages=sapply(strsplit(colnames(TPM), "_"), `[`, 2),
			exp=exp
		)
		df <- dplyr::group_by(gene_exp , stages, sex, cells)
		options(dplyr.summarise.inform = FALSE)
		df.summary2 <- dplyr::summarise(
			df,
			sd = sd(exp),
			len = mean(exp)
		)
		plot <- ggplot(df.summary2, aes(x=stages, y=len, group=interaction(sex, cells), color=sex, linetype=cells))+
		geom_line(linewidth=1)+
		geom_point(size=2.5)+
		geom_errorbar(
			aes(
				ymin=len-sd, 
				ymax=len+sd
			), 
			width=.1
		) +
		scale_color_manual(
			values=c("#FFB100", "#339989"),
			labels=c('Gran.', 'Sert.')
		) +
		scale_y_continuous(labels = scales::comma) +
		coord_cartesian(ylim = c(0, NA)) +
		labs(title=gene, x="Embryonic stages", y = "Expression (TPM)")+
		theme_light() +
		theme(
			plot.title = element_text(size=13, face="bold.italic", hjust = 0.5),
			axis.text = element_text(size=12),
			axis.title.x=element_blank(),
			legend.text = element_text(size =12),
			legend.title = element_blank()
		)
		legend <- get_legend(
			plot + theme(legend.direction="horizontal", legend.box.margin = margin(0, 0, 0, 0))
		)
		plot <- plot + theme(legend.position="none")
		plotlist[[gene]] <- plot
	}

	gtitle <- ggdraw() + 
		draw_label(
			title,
			fontface = 'bold',
			size = 14,
			hjust = 0,
			vjust = 0.5,
			x = 0
		) +
		theme(plot.background = element_rect(fill="#ECF0F5", color = NA),
			plot.margin = margin(0, 0, 0, 25)
		)
	y_axis <- ggdraw() + 
		draw_label(
			"Expression (TPM) vs embryonic stages",
			size = 12,
			hjust = 0.5
		)

	plot_graphs <- plot_grid(plotlist=plotlist, ncol = 5, align = "hv")

	plots <- plot_grid(
		gtitle, plot_graphs,
		ncol = 1,
		rel_heights = c(0.2, 1)
	)
	return(plots)
}

#################################################################################################################################

##########################################
#                                        #
#              Prepare data              #
#                                        #
##########################################

# Rename RNA-seq data to label if whole gonad or supporting cells
rownames(whole_gonad) <- make.names(whole_gonad$X, unique = TRUE)
whole_gonad <- whole_gonad[,-1]
colnames(whole_gonad) <- paste0("whole.gonad_", colnames(whole_gonad))

colnames(TPM) <- paste0("supporting_", colnames(TPM))

# Merge supporting and whole gonad RNA-seq data
TPM <- merge(TPM,whole_gonad,by="row.names",all.x=TRUE)
rownames(TPM) <- TPM[,1]
TPM <- TPM[,-1]

##########################################
#                                        #
#          Plot gene expression          #
#                                        #
##########################################

# Plot expression of the marker genes
plot_list <- plot_expression_scatter(TPM, markerGenes)

figure <- plot_grid(
	plotlist=plot_list,
	labels = "AUTO",
	ncol=1
)

##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################

# As PDF
save_plot(
	snakemake@output[['pdf']], 
	figure, 
	base_width=38,
	base_height=30,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)

# As PNG
save_plot(
	snakemake@output[['png']], 
	figure, 
	base_width=38,
	base_height=30,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)
