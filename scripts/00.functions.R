library("ggplot2")
library("ggrepel")
library("viridis")
library("pheatmap")
library("dplyr")
library("grid")
library("cowplot")
library("DESeq2")
library("clusterProfiler")

###########################################
#                                         #
#              Get matrices               #
#                                         #
###########################################

#' Generate the read count matrix
#' @param csv_file Path to the read count matrix.
#' @param genes Keep all the genes or only the protein coding genes. Values can be "all" or "protein".
#' @return Return a dataframe.
get_read_counts <- function(csv_file, genes){
	# load file
	raw_counts <- read.csv(file=csv_file, row.names=1)
	# Transform values as integers for DESeq2 that does not support floats
	raw_counts <- round(raw_counts, digits = 0)
	# Select protein coding genes
	if (genes=="protein"){
		protein_coding_genes <- read.csv(file="data/mart_prot_coding_genes.txt")
		raw_counts <- raw_counts[rownames(raw_counts) %in% protein_coding_genes$external_gene_name,]
		rm(protein_coding_genes)
	}
	return(raw_counts)
}

#' Generate the tpm matrix
#' @param csv_file Path to the read count matrix.
#' @param genes Keep all the genes or only the protein coding genes. Values can be "all" or "protein".
#' @return Return a dataframe.
get_TPM_counts <- function(csv_file, genes){
	# load file
	tpm <- read.csv(file=csv_file, row.names=1)
	# Select protein coding genes
	if (genes=="protein"){
		protein_coding_genes <- read.csv(file="../data/mart_prot_coding_genes.txt")
		tpm <- tpm[rownames(tpm) %in% protein_coding_genes$external_gene_name,]
		rm(protein_coding_genes)
	}
	return(tpm)
}

###########################################
#                                         #
#          Normalize read counts          #
#                                         #
###########################################

#' Normalize the read count using the size factor normalization from DESeq2
#' @param raw_counts Read count matrix.
#' @param samplesheet Samplesheet for DESeq2.
#' @return Return a dataframe.
get_normalized_counts <- function(raw_counts, samplesheet){
	dds <- DESeq2::DESeqDataSetFromMatrix(
		countData = raw_counts,
		colData = samplesheet,
		design = ~conditions
	)
	dds <- estimateSizeFactors(dds)
	norm_counts <- counts(dds, normalized=TRUE)
	# norm_counts <- assay(vst(dds, blind=FALSE))

	return(norm_counts)
}

###########################################
#                                         #
#            Filter low counts            #
#                                         #
###########################################

filter_low_counts <- function(row, col_names, minReads) {
	if (max(row) < minReads) {
		return(setNames(rep(0, length(row)), col_names))
	} else {
		return(setNames(row, col_names))
	}
}

run_filter_low_counts <- function(data, minReads) {
	col_names <- colnames(data)
	data <- t(apply(data, 1, filter_low_counts, col_names = col_names, minReads = minReads))
	return(as.data.frame(data))
}

###################################################################################
#																				  #
#								Correlation Functions							  #
#																				  #
###################################################################################

correlation <- function(matrix, conditions, method, colours){
	sample_names <- paste(
		sapply(strsplit(colnames(matrix), "_"), `[`, 2), 
		sapply(strsplit(colnames(matrix), "_"), `[`, 1), 
		sapply(strsplit(colnames(matrix), "_"), `[`, 4), 
		sep="_"
	)
	colnames(matrix) <- sample_names
	matrix <- matrix[ ,order(names(matrix))]
	cor_data <- cor(matrix, method=method)
	cor_data[cor_data==1.000]<-NA

	col_annotation <- data.frame(
		Samples=conditions[order(conditions)]
	)

	rownames(col_annotation) <- names(matrix)

	col <- list(Samples=colours)

	pheatmap::pheatmap(
		mat = cor_data,
		cluster_cols = FALSE,
		cluster_rows = FALSE,
		color = viridis(15),
		annotation_colors = col,
		border_color = NA,
		cellwidth = 12,
		cellheight = 12,
		annotation_col = col_annotation,
		annotation_row = col_annotation
	)
}

run_pca <- function(rpkm_table){
	print("Calculating the PCA...")
	log.rpkm <- t(rpkm_table)
	log.rpkm.no0 <- log.rpkm[, colSums(log.rpkm)!=0]
	pca <- prcomp(
		log.rpkm.no0,
		center=TRUE,
		scale.=TRUE
	)
	return(pca)
}

roundUp <- function(x, nice=seq(1,10, 0.5)) {
	if(length(x) != 1) stop("'x' must be of length 1")
	10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

plot.pca <- function(
	count_table,
	conditions,
	colours,
	PCs
	)
	{
		pca <- run_pca(count_table)
		print("Generating the plots...")
		percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2))*100
		# print(percent_var_explained)
		cond <- factor(conditions)
		col <- factor(conditions)
		levels(col) <- colours
		col <- as.vector(col)
		conditions_2 <- conditions
		conditions_2[duplicated(conditions_2)] <- ""
		scores <- as.data.frame(pca$x)
		PCs.combinations <- combn(PCs,2)
		# print(PCs.combinations)
		# max.score <- roundUp(max(abs(scores[,PCs.combinations])))
		# max.score[max.score<150] <- 150
		# print(head(scores))
		plots <- apply(
			PCs.combinations,
			2,
			function(combination)
			{
				data <- scores[,c(combination[1], combination[2])]
				data$cond <- conditions
				colnames(data) <- c("PC_x", "PC_y", "cond")
				# print(head(data))
				plot <- ggplot(data, aes(x=PC_x, y=PC_y, fill=cond)) +
					geom_point(shape = 21, size = 6, stroke=0.5) +
					scale_fill_manual(values = conditions_color) +
					# geom_label_repel(
					# 	aes(label = conditions_2),
					# 	size = 4,
					# 	segment.size  = 0.8,
					# 	segment.color = "grey50",
					# 	color=col,
					# 	fill = alpha("white",0.5),
					# 	box.padding = 0.7,
					# 	max.overlaps = 15
					# ) +
					# ggtitle("PCA on TPM") +
					theme_bw() +
					xlab(paste(combination[1], " ", "(",round(percent_var_explained[as.numeric(gsub("PC", "", combination[1]))], digit=2),"%)", sep=""))+
					ylab(paste(combination[2], " ", "(",round(percent_var_explained[as.numeric(gsub("PC", "", combination[2]))], digit=2),"%)", sep=""))+
					theme(
						plot.title = element_text(size=10, face="bold", hjust = 0.5),
						axis.text=element_text(size=10),
						axis.title=element_text(size=10),
						aspect.ratio=1,
						legend.title=element_blank(),
						legend.text=element_text(size=10)
					)
				return(plot)
			}
		)
		print("Done.")
		return(plots)
}


plot_volcano <- function(de_res, ad.pval, FC, cond1, cond2, nb_sig_DE, colors){
	pval_order <- de_res[with(de_res, order(padj,decreasing = FALSE)),]
	# print(head(pval_order))
	top5_p_up <- rownames(head(pval_order[pval_order$Diff.Exp.==paste("Up in", cond1),],20))
	top5_p_down <- rownames(head(pval_order[pval_order$Diff.Exp.==paste("Up in", cond2),],20))
	sig_genes <- de_res[rownames(de_res) %in% c(top5_p_up, top5_p_down),]

	up_in_cond1 <- nrow(pval_order[pval_order$Diff.Exp.==paste("Up in", cond1),])
	up_in_cond2 <- nrow(pval_order[pval_order$Diff.Exp.==paste("Up in", cond2),])

	# print(sig_genes)
	de_res <- na.omit(de_res)
	max_FC <- max(abs(de_res$log2FoldChange))
	# print(max_FC)
	# print(head(de_res))
	volcano <- ggplot(de_res, aes(
		x = log2FoldChange*(-1), 
		y = -log10(padj)
		)
	) +
		geom_point(
			# shape=21,
			# colour="black"
			aes(
				colour = Diff.Exp.,
				fill = Diff.Exp.
			),
			alpha = 0.5,
			shape = 21,
			size = 2
		) +
		geom_point(
			data = sig_genes,
			shape = 21,
			size = 3,  
			colour = "black"
		) +
		geom_hline(
			yintercept = -log10(ad.pval),
			linetype = "dashed"
		) + 
		geom_vline(
			xintercept = c(FC*(-1), FC),
			linetype = "dashed"
		) +
		geom_label_repel(
			data = sig_genes,
			aes(
				label = rownames(sig_genes)
			),
			fill = alpha(c("white"),0.5),
			size = 4,
			fontface="italic",
			force = 2,
			nudge_y = 1,
			max.overlaps = 18
		) +
		ggtitle(paste0(cond1, " vs ", cond2, "\n", nb_sig_DE, " significant DE genes")) +
		# xlim(-max_FC, max_FC) +
		theme_bw() +
		xlab("log2(fold change)") +
		ylab("-log10(adjusted p-values)") +
		theme(
			# axis.text.x = element_text(angle = 45, hjust=1),
			plot.title = element_text(size=12, hjust = 0.5, face="bold"),
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			legend.title=element_blank(),
			legend.text=element_text(size=12, vjust = 0.5, hjust = 0.5),
			legend.position="bottom",
			plot.margin = unit(c(0.42,0,0,0), "cm"),
			aspect.ratio=0.8
		) +
		scale_color_manual(	
			values=c("#5D576B", as.vector(colors)),
			labels = c("non-sig.", paste0("Up in ", cond1, "\n(", up_in_cond1," genes)"), paste0("Up in ", cond2, "\n(", up_in_cond2," genes)"))
		) +
		scale_fill_manual(
			values=c("#5D576B", as.vector(colors)),
			labels = c("non-sig.", paste0("Up in ", cond1, "\n(", up_in_cond1," genes)"), paste0("Up in ", cond2, "\n(", up_in_cond2," genes)"))
		)
	return(volcano)
}


plot_dymorphic_genes <- function(de_genes, XX_colors, XY_colors){
	axis_margin <- 5.5
	xx_sex_de_genes <- de_genes[de_genes$sex=="XX",]
	p1 <- ggplot(xx_sex_de_genes, aes(x=DE, y=stage)) +
		geom_col(fill=XX_colors) +
		geom_text(aes(label = scales::comma(abs(DE))), x=-500, size=4.5, color="white") +
		scale_x_reverse(labels = scales::comma, limits=c(4100,0)) +
		scale_y_discrete(position = "right", limits=rev) +
		ggtitle("XX overexpressed genes") +
		xlab("Gene count") +
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

	xy_sex_de_genes <- de_genes[de_genes$sex=="XY",]
	p2 <- ggplot(xy_sex_de_genes, aes(x=DE, y=stage)) +
		geom_col(fill=XY_colors) +
		geom_text(aes(label = scales::comma(abs(DE))), x=500, size=4.5, color="white") +
		scale_x_continuous(labels = scales::comma, limits=c(0,4200)) +
		scale_y_discrete(limits=rev) +
		ggtitle("XY overexpressed genes") +
		xlab("Gene count") +
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
			"Sexually dimorphic genes per stages",
			fontface = 'bold',
			size = 14,
			hjust = 0.5
		)

	plot_graphs <- plot_grid(p1, p2, align="hv", axis="bt", nrow=1)

	plot <- plot_grid(gtitle, plot_graphs, nrow=2, rel_heights = c(0.1, 1))

	return(plot)

}


plot_dymorphic_genes_2 <- function(de_genes, XX_colors, XY_colors){
	axis_margin <- 5.5
	xx_sex_de_genes <- de_genes[de_genes$sex=="XX",]
	p1 <- ggplot(xx_sex_de_genes, aes(x=stage, y=DE)) +
		geom_col(fill=XX_colors) +
		geom_text(aes(label = scales::comma(abs(DE))), y=500, size=4.5, color="white") +
		scale_y_continuous(position = "right", labels = scales::comma, limits=c(0,4200)) +
		scale_x_discrete(limits=rev) +
		ggtitle("XX overexpressed genes") +
		xlab("Gene count") +
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

	xy_sex_de_genes <- de_genes[de_genes$sex=="XY",]
	p2 <- ggplot(xy_sex_de_genes, aes(x=stage, y=DE)) +
		geom_col(fill=XY_colors) +
		geom_text(aes(label = scales::comma(abs(DE))), y=500, size=4.5, color="white") +
		scale_y_continuous(labels = scales::comma, limits=c(0,4200)) +
		# scale_x_discrete(limits=rev) +
		ggtitle("XY overexpressed genes") +
		xlab("Gene count") +
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
			"Sexually dimorphic genes per stages",
			fontface = 'bold',
			size = 14,
			hjust = 0.5
		)

	plot_graphs <- plot_grid(p1, p2, align="hv", axis="bt", nrow=1)

	plot <- plot_grid(gtitle, plot_graphs, nrow=2, rel_heights = c(0.1, 1))

	return(plot)

}


plot_simple_heatmap <- function(data, de_feature, colors, clusters, res_file){
	matrix_DEG <- data[rownames(data) %in% de_feature, ]

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
	conditions <-sapply(strsplit(colnames(matrix), "_"), `[`, 1)
	annotation_col <- data.frame(
		Stages=conditions
	)
	rownames(annotation_col) <- colnames(matrix)

	# Stages palette
	col_stage <- colors
	names(col_stage) <- unique(conditions)

	# Gene cluster color palette
	cluster_colors <- as.vector(MetBrewer::met.brewer("Hokusai1", n=clusters))

	# Color palette for the heatmap
	cold <- colorRampPalette(c('#e9f6e6','#e0f3f7','#8ab8d7','#4575b4'))
	warm <- colorRampPalette(c('#f4fbd2','#feeda3','#fa8a57','#d73027'))
	mypalette <- c(rev(cold(12)), warm(12))

	# set.seed(654)

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

	TF_list <- read.csv("../data/mouse_transcription_factors.txt", header=FALSE)
	gonad_pheno_genes <- read.csv("../data/genes_gonad_asociated_pheno_MGI.txt", header=FALSE)

	TF_list <- as.vector(TF_list[,1])
	total_TFs <- unlist(lapply(TF_list, function(TF) which(rownames(matrix) %in% TF)))
	print(paste(length(total_TFs), "TFs found in total."))

	gonad_pheno_genes <- as.vector(gonad_pheno_genes[,1])

	TFs <- TF_list[TF_list %in% gonad_pheno_genes]

	matrix_TF_indexes <- unlist(lapply(TFs, function(TF) which(rownames(matrix) %in% TF)))
	print(paste(length(matrix_TF_indexes), "TFs found associated with gonadal phenoypes."))

	if(length(matrix_TF_indexes)>30){
		matrix_TF_indexes <- sample(matrix_TF_indexes,30)
	}

	TF_names <- rownames(matrix)[matrix_TF_indexes]

	# Show transcription factors
	TFs = rowAnnotation(
		TFs = anno_mark(
			at = matrix_TF_indexes, # TF row indexes
			labels = TF_names,    # Gene names
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
	height <- 9
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


plot_expression_heatmap <- function(data, markers, colors, clusters){
	matrix_DEG <- data[rownames(data) %in% markers$Gene, ]

	matrix <- matrix_DEG[match(markers$Gene, rownames(matrix_DEG)),]

	# matrix <- log10(matrix+1)

	# matrix <- matrix_DEG
	matrix[matrix>100] <- 100


	# Prepare heatmap
	conditions <-paste(sapply(strsplit(colnames(matrix), "_"), `[`, 2), sapply(strsplit(colnames(matrix), "_"), `[`, 1))
	matrix <- matrix[,order(conditions)]
	conditions <-paste(sapply(strsplit(colnames(matrix), "_"), `[`, 2), sapply(strsplit(colnames(matrix), "_"), `[`, 1))

	annotation_col <- data.frame(
		Stages=conditions,
		Sex=paste(sapply(strsplit(colnames(matrix), "_"), `[`, 2))
	)
	rownames(annotation_col) <- colnames(matrix)

	col_stage <- colors
	names(col_stage) <- unique(conditions)

	Sex <- paste(sapply(strsplit(colnames(matrix), "_"), `[`, 2))

	# set.seed(654)

	ha1 <- HeatmapAnnotation(
		Stages=conditions,
		col=list(
			Stages=col_stage
		),
		annotation_legend_param = list(
			Stages = list(direction = "horizontal")
		),
		annotation_name_side = "left",
		simple_anno_size = unit(1, "cm")
	)

	Cell.Type <- markers$Cell_type
	names(Cell.Type) <- markers$Gene

	ha2 <- rowAnnotation(
		Cell.Type = Cell.Type,
		annotation_legend_param = list(
			Cell.Type = list(direction = "horizontal")
		),
		simple_anno_size = unit(0.5, "cm")
	)

	# Color palette for the heatmap
	cold <- colorRampPalette(c('#e9f6e6','#e0f3f7','#8ab8d7','#4575b4'))
	warm <- colorRampPalette(c('#f4fbd2','#feeda3','#fa8a57','#d73027'))
	mypalette <- c(rev(cold(12)), warm(12))

	ht_list <- Heatmap(
				matrix, 
				name = "TPM",
				top_annotation = ha1,
				right_annotation = ha2,
				row_names_gp = gpar(fontface = "italic"),
				row_dend_reorder = FALSE,
				column_split = Sex,
				cluster_columns = FALSE,
				cluster_rows = FALSE,
				show_column_names = FALSE,
				col = viridis::viridis(50),
				heatmap_legend_param = list(direction = "horizontal")
			)
	draw(
		ht_list, 
		merge_legend = TRUE, 
		heatmap_legend_side = "right", 
		annotation_legend_side = "right"
	)
}


plot_expression_scatter <- function(data, markers){
	cellTypes <- unique(markers$Cell_type)
	plot_list <- lapply(cellTypes, function(cellType){
		genes <- markers[markers$Cell_type %in% cellType, "Gene"]
		gene_exp(genes, data, cellType)
	})
}

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
		df <- group_by(gene_exp , stages, sex, cells)
		options(dplyr.summarise.inform = FALSE)
		df.summary2 <- summarise(
			df,
			sd = sd(exp),
			len = mean(exp)
		)

		plot <- ggplot(df.summary2, aes(x=stages, y=len, group=interaction(sex, cells), color=sex, linetype=cells))+
		geom_line(size=1)+
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
		# coord_fixed() +
		# theme_bw()+
		# ggtitle(gene)+
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

GO_term_per_cluster <- function(de_genes, res_file){

	print("Calculate GO term over-representation...")
	formula_res <- compareCluster(
		gene~cluster, 
		data=de_genes, 
		fun="enrichGO", 
		keyType = "SYMBOL",
		OrgDb="org.Mm.eg.db",
		ont		   = "BP",
		pAdjustMethod = "BH",
		pvalueCutoff  = 0.01,
		qvalueCutoff  = 0.05,
		readable = TRUE
	)

	print("Calculate GO term semantic similarities...")
	lineage1_ego <- simplify(
		formula_res, 
		cutoff=0.6,
		by="p.adjust", 
		select_fun=min
	)

	write.csv(lineage1_ego, file=res_file, quote=FALSE)
	
	return(lineage1_ego)
}

go_plot <- function(go_res, nb_terms=5){
	# Make the first GO term letter as capital letter
	go_res@compareClusterResult[,4] <- gsub("^([a-z])", "\\U\\1", go_res[,4], perl=TRUE)
	options(enrichplot.colours = c("#77BFA3", "#98C9A3", "#BFD8BD", "#DDE7C7", "#EDEEC9"))
	plot <- clusterProfiler::dotplot(go_res, showCategory=nb_terms)
	x_labels <-  levels(as.data.frame(go_res)$Cluster)
	# Reload ggplot2 to apply new theme
	library(ggplot2)
	plot <- plot +
		geom_point(shape=21) +
		ggtitle(paste0("Enriched biological process GO terms (Top ", nb_terms, ")")) +
		scale_x_discrete(labels=x_labels) +
		labs(color = "Adj. p-value", size="Gene ratio") +
		guides(color = guide_colorbar(reverse=TRUE), size = guide_legend(reverse=TRUE)) +
		theme(
			plot.title = element_text(size=12, hjust = 0.5, face="bold"),
			axis.title.x = element_blank()
		)
	return(plot)
}

get_sex_DEG_per_stage <- function(dds, stage, p.adj, log2FC){
	res <- DESeq2::results(dds, contrast=c("conditions", paste("XX", stage), paste("XY", stage)))

	de_res <- as.data.frame(res)
	res <- mutate(
		de_res, 
		Diff.Exp. = dplyr::case_when(
			log2FoldChange >= log2FC & padj <= p.adj ~ "Up in XX",
			log2FoldChange <= (-log2FC) & padj <= p.adj ~ "Up in XY",
			TRUE ~ "non sig."
		)
	)
	sig.DE <- subset(res, padj < p.adj)
	sig.DE <- subset(sig.DE, abs(log2FoldChange) > log2FC)

	return(sig.DE)
}


plot_volcano_sex <- function(dds, stage, p.adj, log2FC, colors){

	colors <- conditions_color[grep(stage, names(conditions_color))]

	res <- results(dds, contrast=c("conditions", paste("XX", stage), paste("XY", stage)))
	de_res <- as.data.frame(res)

	res <- mutate(
		de_res, 
		Diff.Exp. = dplyr::case_when(
			log2FoldChange >= log2FC & padj <= p.adj ~ "Up in XX",
			log2FoldChange <= (-log2FC) & padj <= p.adj ~ "Up in XY",
			TRUE ~ "non sig."
		)
	)

	res$test <- -log10(res$padj)
	res$test[res$test==Inf] <- 400

	# print(head(res))

	sig.DE <- subset(res, padj < p.adj)
	sig.DE <- subset(sig.DE, abs(log2FoldChange) > log2FC)

	nb_sig_DE <- scales::comma(nrow(sig.DE))

	# print(nb_sig_DE)

	pval_order <- res[with(res, order(test, decreasing = TRUE)),]

	# print(head(pval_order))

	top_p_up <- rownames(head(pval_order[pval_order$Diff.Exp.=="Up in XX",],10))
	top_p_down <- rownames(head(pval_order[pval_order$Diff.Exp.=="Up in XY",],10))
	top_genes <- c(top_p_up, top_p_down)
	top_genes <- c("Sox9", "Amh", "Dhh", "Fgf9", "Ptgds", "Sry", "Nr0b1", "Dmrt1", "Wt1", "Nr5a1", "Foxl2", "Fst", "Wnt4", "Lgr5", "Ctnnb1", "Bmp2", "Runx1", "Cdkn1a")
	sig_genes <- sig.DE[rownames(sig.DE) %in% top_genes,]

	# print(sig_genes)

	up_in_cond1 <- nrow(pval_order[pval_order$Diff.Exp.=="Up in XX",])
	up_in_cond2 <- nrow(pval_order[pval_order$Diff.Exp.=="Up in XY",])

	# print(up_in_cond1)

	de_res <- na.omit(res)
	max_FC <- max(abs(res$log2FoldChange))
	# print(max_FC)
	# print(head(de_res))

	volcano <- ggplot(
		res, 
		aes(
			x = -log2FoldChange, 
			y = test
		)
	) +
		geom_point(
			# shape=21,
			# colour="black"
			aes(
				colour = Diff.Exp.,
				fill = Diff.Exp.
			),
			alpha = 0.5,
			shape = 21,
			size = 2
		) +
		geom_point(
			data = sig_genes,
			shape = 21,
			size = 3,  
			colour = "black"
		) +
		geom_hline(
			yintercept = -log10(p.adj),
			linetype = "dashed"
		) + 
		geom_vline(
			xintercept = c(-log2FC, log2FC),
			linetype = "dashed"
		) +
		geom_label_repel(
			data = sig_genes,
			mapping = aes(
				label = rownames(sig_genes)
			),
			color = 'black',
			size = 4,
			min.segment.length = unit(0, 'lines'),
			fill = alpha("white", 0.5),
			fontface="italic",
			force = 3,
			# nudge_y = 1,
			max.overlaps = 25,
			# seed=123
		) +
		expand_limits(
			x = c(-25, 25),
			y = c(0, 450)
		) +
		ggtitle(paste0(stage, " XX vs ", stage," XY \n", nb_sig_DE, " significant DE genes")) +
		# xlim(-max_FC, max_FC) +
		theme_bw() +
		xlab("log2(fold change)") +
		ylab("-log10(adjusted p-values)") +
		theme(
			# axis.text.x = element_text(angle = 45, hjust=1),
			plot.title = element_text(size=12, hjust = 0.5, face="bold"),
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			legend.title=element_blank(),
			legend.text=element_text(size=12, vjust = 0.5, hjust = 0.5),
			legend.position="bottom",
			# plot.margin = unit(c(0.42,0,0,0), "cm"),
			aspect.ratio=0.8
		) +
		scale_color_manual(	
			values=c("#5D576B", as.vector(colors)),
			labels = c("non-sig.", paste0("Up in XX\n(", scales::comma(up_in_cond1)," genes)"), paste0("Up in XY\n(", scales::comma(up_in_cond2)," genes)"))
		) +
		scale_fill_manual(
			values=c("#5D576B", as.vector(colors)),
			labels = c("non-sig.", paste0("Up in XX\n(", scales::comma(up_in_cond1)," genes)"), paste0("Up in XY\n(", scales::comma(up_in_cond2)," genes)"))
		)
	return(volcano)
}

plot_DEG <- function(dds, stage, p.adj, log2FC, colors){
	volcano <- plot_volcano_sex(dds, stage, p.adj, log2FC, colors)

	res <- DESeq2::results(dds, contrast=c("conditions", paste("XX", stage), paste("XY", stage)))

	de_res <- as.data.frame(res)
	res <- mutate(
		de_res, 
		Diff.Exp. = dplyr::case_when(
			log2FoldChange >= log2FC & padj <= p.adj ~ "Up in XX",
			log2FoldChange <= (-log2FC) & padj <= p.adj ~ "Up in XY",
			TRUE ~ "non sig."
		)
	)
	sig.DE <- subset(res, padj < p.adj)
	sig.DE <- subset(sig.DE, abs(log2FoldChange) > log2FC)

	de_genes <- data.frame(
		gene=rownames(sig.DE),
		cluster=sig.DE$Diff.Exp.
	)

	res_file <- paste0("../results/240306_", stage, "_GO_sex_DEG.csv")
	GO_terms <- GO_term_per_cluster(de_genes, res_file)
	go_term_plot <- go_plot(GO_terms, nb_terms=5)

	figure <- plot_grid(plotlist=list(volcano, go_term_plot), ncol=2)

	return(figure)
}

upset_plots_sex <- function(DEGs, sex){

	current_sex <- paste("Up in", sex)
	sex_DEG <- lapply(DEGs, function(DEG) rownames(DEG[DEG$Diff.Exp.==current_sex,,drop=FALSE]))
	venn <- eulerr::euler(sex_DEG)
	if (sex=="XX"){
		color <- "#FFB100"
		text_col <- c("black", "black")
	} else {
		color="#339989"
		text_col <- c("black", "white")
	}

	theme_upset <- function(){
		theme(
			axis.text=element_text(size=12),
			axis.title=element_text(size=12),
			legend.position="none"
		)
	}

	venn <- eulerr::euler(sex_DEG)
	data <- fromExpression(venn$original.values)
	colnames(data) <- names(DEGs)
	data <- data==1
	data <- as.data.frame(data)

	plot <- ComplexUpset::upset(
		data=data,
		intersect=names(DEGs),
		name="Sex-specific DE genes",
		sort_sets='ascending',
		stripes=alpha('white', 0),
		set_sizes=FALSE,
		# wrap=TRUE
		base_annotations=list(
			'Intersection size'=ComplexUpset::intersection_size(
				text_colors=text_col,
				mapping=aes(fill='bars_color')
			) +
			scale_fill_manual(values=c('bars_color'=color)) +
			scale_y_continuous(labels = scales::comma) +
			theme_upset()
		),
		matrix=ComplexUpset::intersection_matrix(
			segment=geom_segment(color=color),
			outline_color=list(active=color, inactive="grey70"),
		) + scale_color_manual(
			values=c('TRUE'=color, 'FALSE'='grey'),
			na.value='transparent',
			guide='none'
		),
		themes=ComplexUpset::upset_default_themes(
			axis.text=element_text(size=12),
			axis.title=element_text(size=12)
		)
	) 

	return(plot)
}


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

	TF_list <- read.csv("../data/mouse_transcription_factors.txt", header=FALSE)
	gonad_pheno_genes <- read.csv("../data/genes_gonad_asociated_pheno_MGI.txt", header=FALSE)

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



draw_venn_sex_dyn <- function(sex, dyn, title) {
	data <- c(
		sex=length(setdiff(sex,dyn)),
		dyn=length(setdiff(dyn,sex)),
		"sex&dyn"=length(intersect(sex,dyn))
	)

	colours <- c( "#FCBF49", "#D62828")
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


###############################################
# Fix complexeheatmap gene annotation

# library(ggplot2)
library(gtable)
library(grid)
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


