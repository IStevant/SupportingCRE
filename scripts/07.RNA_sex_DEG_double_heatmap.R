# source("scripts/00.functions.R")
source("scripts/00.color_palettes.R")

# install.packages('doParallel', repos="https://mirror.ibcp.fr/pub/CRAN/")
# install.packages('foreach', repos="https://mirror.ibcp.fr/pub/CRAN/")

# install.packages('MetBrewer', repos="https://mirror.ibcp.fr/pub/CRAN/")

# BiocManager::install('org.Mm.eg.db')
# BiocManager::install("ComplexHeatmap", force = TRUE)
# renv::snapshot()

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
	library("ComplexHeatmap")
	library("cowplot")
	library("gtable")
	library("grid")
	library("ggplot2")
})

# renv::snapshot()

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

plot_double_heatmap <- function(data, de_feature, colors, clusters, res_file){
	matrix_DEG <- data[rownames(data) %in% de_feature, ]
	# Split data by sex and invert XX data to prepare the double heatmap
	XX <- matrix_DEG[,grep("XX", colnames(matrix_DEG))]
	XX <- XX[,order(ncol(XX):1)]
	XY <- matrix_DEG[,grep("XY", colnames(matrix_DEG))]
	matrix_DEG <- cbind(XX, XY)

	# Calculate z-scores
	matrix <- t(scale(t(matrix_DEG)))

	# Cluster matrix
	row_dend <- hclust(dist(matrix), method= "ward.D2")
	clustering <- cutree(row_dend, k=clusters)

	# Order clusters
	mean_per_cluster <- lapply(
		unique(clustering),
		function(cluster){
			genes <- names(clustering[clustering==cluster])
			matrix_cluster <- matrix[rownames(matrix) %in% genes,]
			mean <- as.data.frame(t(colMeans(matrix_cluster)))
			return(mean)
		}
	)

	mean_per_cluster <- data.table::rbindlist(mean_per_cluster)
	o1 <- seriation::seriate(dist(mean_per_cluster), method="GW")
	levels <- seriation::get_order(o1)

	clustering <- factor(clustering, levels=levels)
	levels(clustering) <- letters[1:clusters]

	# Save clustering
	write.csv(clustering, file=res_file, quote=FALSE)

	# Limit zscore to |2|
	matrix[matrix>2] <- 2
	matrix[matrix<(-2)] <- (-2)

	# Prepare top annotation
	conditions <-paste(sapply(strsplit(colnames(matrix), "_"), `[`, 2), sapply(strsplit(colnames(matrix), "_"), `[`, 1))
	annotation_col <- data.frame(
		Stages=conditions
	)
	rownames(annotation_col) <- colnames(matrix)

	# Stages palette
	col_stage <- colors
	# names(col_stage) <- c(unique(conditions), rev(unique(conditions)))

	# Gene cluster color palette
	cluster_colors <- as.vector(MetBrewer::met.brewer("Hokusai1", n=clusters))

	# Color palette for the heatmap
	cold <- colorRampPalette(c('#e9f6e6','#e0f3f7','#8ab8d7','#4575b4'))
	warm <- colorRampPalette(c('#f4fbd2','#feeda3','#fa8a57','#d73027'))
	mypalette <- c(rev(cold(12)), warm(12))

	# set.seed(654)

	Sex <- sapply(strsplit(colnames(matrix), "_"), `[`, 2)

	# Top annotation (stages)
	stage_anno <- HeatmapAnnotation(
		Stages=conditions,
		col=list(
			Stages=col_stage
		),
		annotation_name_side = "left",
		simple_anno_size = unit(0.5, "cm")
	)
 
	# Row annotation (gene clusters)
	cluster_anno <- rowAnnotation(
		Clusters = anno_empty(
			border = FALSE,
			width = unit(10, "mm")
		)
	)

	# Prepare the gene clusters legend manually
	cluster_legend = Legend(
		at = levels(clustering), 
		title = "Clusters", 
		legend_gp = gpar(fill = cluster_colors)
	)

	TF_list <- read.csv("data/mouse_transcription_factors.txt", header=FALSE)
	gonad_pheno_genes <- read.csv("data/genes_gonad_asociated_pheno_MGI.txt", header=FALSE)

	TF_list <- as.vector(TF_list[,1])
	total_TFs <- unlist(lapply(TF_list, function(TF) which(rownames(matrix) %in% TF)))
	# print(paste(length(total_TFs), "TFs found in total."))

	gonad_pheno_genes <- as.vector(gonad_pheno_genes[,1])

	TFs <- TF_list[TF_list %in% gonad_pheno_genes]

	matrix_TF_indexes <- unlist(lapply(TFs, function(TF) which(rownames(matrix) %in% TF)))
	# print(paste(length(matrix_TF_indexes), "TFs found associated with gonadal phenoypes."))

	if(length(matrix_TF_indexes)>30){
		matrix_TF_indexes <- sample(matrix_TF_indexes,30)
	}

	TF_names <- rownames(matrix[matrix_TF_indexes,])

	# Show transcription factors
	TFs = rowAnnotation(
		TFs = anno_mark(
			# at = matrix_TF_indexes, # TF row indexes
			at = which(rownames(matrix) %in% TF_names),
			labels = rownames(matrix)[rownames(matrix)%in%TF_names],    # Gene names
			labels_gp = gpar(fontsize=10, fontface = "italic"),
			padding=unit(2, "mm")
		)
	)

	# Make the heatmap
	ht_list <- Heatmap(
		matrix, 
		name = "z-score",
		top_annotation = stage_anno,
		left_annotation = cluster_anno,
		right_annotation = TFs,
		column_split = Sex,
		row_title_rot = 0,
		row_split = clustering,
		cluster_columns = FALSE,
		show_column_names = FALSE,
		show_row_names = FALSE,
		show_row_dend = FALSE,
		cluster_row_slices = FALSE, 
		col = mypalette,
		heatmap_legend_param = list(direction = "vertical"),
		row_title=NULL
	)

	height <- 8
	width <- 6

	gTree <- grid.grabExpr(
		{

		# Draw the heatmap
		draw(
			ht_list,
			row_title = paste(scales::comma(nrow(matrix)), "differentially expressed genes"),
			annotation_legend_list=cluster_legend,
			merge_legend = TRUE,
			use_raster = TRUE, 
			raster_quality = 5
		)

		# Add gene cluster annotation as extra-rectangles to avoid bad rasterization (faded colors)
		for(i in 1:length(levels(clustering))) {
			decorate_annotation(
				"Clusters", 
				slice = i, {
					grid.rect(x=1, width = unit(3, "mm"), gp = gpar(fill = cluster_colors[i], col = NA), just = "right")
					grid.text(x=0.15, levels(clustering)[i], just = "left")
				}
			)
		}

	},
		height = height,
		width = width
	)


	p_fix <- panel_fix(plot_grid(gTree), width = width, height = height, units="inch")
	p <- plot_grid(p_fix)

	return(p_fix)
}


###############################################
# Fix complexeheatmap gene annotation

panel_fix <- function(p = NULL, grob = NULL,
						width = NULL, height = NULL, margin = 1, units = "cm",
						filename = NULL) {
	if (is.null(p) & is.null(grob)) {
			stop("'p' or 'grob' must be provided with at least one.")
	}
	if (is.null(width) & is.null(height)) {
			stop("'width' or 'height' must be provided with at least one.")
	}
	if (is.null(grob)) {
			grob <- ggplotGrob(p)
	}

	panels <- grep("panel", grob[["layout"]][["name"]])
	panel_index_h <- sort(unique(grob[["layout"]][["t"]][panels]))
	panel_index_w <- sort(unique(grob[["layout"]][["l"]][panels]))
	nw <- length(panel_index_w)
	nh <- length(panel_index_h)
	raw_w <- as.numeric(grob[["widths"]][panel_index_w])
	raw_h <- as.numeric(grob[["heights"]][panel_index_h])
	raw_aspect <- raw_h / raw_w
	if (is.null(width)) {
			width <- height / raw_aspect
	}
	if (is.null(height)) {
			height <- width * raw_aspect
	}
	if (!length(width) %in% c(1, length(raw_aspect)) | !length(height) %in% c(1, length(raw_aspect))) {
			stop("The length of 'width' and 'height' must be 1 or the length of panels.")
	}
	if (length(width) == 1) {
			width <- rep(width, nw)
	}
	if (length(height) == 1) {
			height <- rep(height, nh)
	}
	grob[["widths"]][panel_index_w] <- unit(width, units = units)
	grob[["heights"]][panel_index_h] <- unit(height, units = units)
	grob <- gtable_add_padding(grob, unit(margin, units = units))
	plot_width <- convertWidth(sum(grob[["widths"]]), unitTo = units, valueOnly = TRUE)
	plot_height <- convertHeight(sum(grob[["heights"]]), unitTo = units, valueOnly = TRUE)

	if (!is.null(filename)) {
			ggsave(filename = filename, plot = grob, units = units, width = plot_width, height = plot_height)
	}
	attr(grob, "size") <- list(units = units, width = plot_width, height = plot_height)
	return(grob)
}

#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

load(snakemake@input[['sig_DEGs']])
norm_counts <- read.csv(file=snakemake@input[['norm_counts']], row.names=1)
samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)
stages <- unique(samplesheet$stages)
names(conditions_color) <- sort(unique(samplesheet$conditions))

clusters <- snakemake@params[['clusters']]

###########################################
#                                         #
#             Double heatmap              #
#                                         #
###########################################

de_feature <- unique(unlist(lapply(filtered_SexDEGs, rownames)))
sex_DEG_heatmap <- plot_double_heatmap(
	norm_counts, 
	de_feature, 
	conditions_color, 
	clusters, 
	snakemake@output[['clusters']]
)

# plot_grid(sex_DEG_heatmap)

save_plot(
	snakemake@output[['pdf']], 
	sex_DEG_heatmap, 
	base_width=15,
	base_height=20.5,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)

save_plot(
	snakemake@output[['png']], 
	sex_DEG_heatmap, 
	base_width=15,
	base_height=20.5,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)
