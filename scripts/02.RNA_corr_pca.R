# source("scripts/00.functions.R")
source("scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################
library("cowplot")
library("grid")
library("viridis")
library("ggplot2")

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw the correlation matrix between the samples
#' @param matrix Expression matrix.
#' @param conditions Vector containing the names of the conditions of the samples. The length should be the same as the number of samples.
#' @param method Method for the correlation (Example: "Spearman" or "Pearson"). Default is "Spearman".
#' @param colours Vector containing the hexadecimal colours corresponding to each conditions.
#' @return Pheatmap object.
correlation <- function(matrix, conditions, method="Spearman", colours){
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

#' Proceed to the PCA using prcomp.
#' @param matrix Expression matrix.
#' @return prcomp object.
run_pca <- function(matrix){
	print("Calculating the PCA...")
	t.matrix <- t(matrix)
	t.matrix.no0 <- t.matrix[, colSums(t.matrix)!=0]
	pca <- prcomp(
		t.matrix.no0,
		center=TRUE,
		scale.=TRUE
	)
	return(pca)
}

#' Plot the PCA.
#' @param matrix Expression matrix.
#' @param conditions Vector containing the names of the conditions of the samples. The length should be the same as the number of samples.
#' @param colours Vector containing the hexadecimal colours corresponding to each conditions.
#' @param PCs Vector PCs to plot. Default is "PC1" vs "PC2". If more than two PCs provided, all the possible combinations are drawn.
#' @return Pheatmap object.
plot.pca <- function(matrix, conditions, colours, PCs=c("PC1", "PC2")){
		pca <- run_pca(matrix)
		print("Generating the plots...")
		percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2))*100
		cond <- factor(conditions)
		col <- factor(conditions)
		levels(col) <- colours
		col <- as.vector(col)
		conditions_2 <- conditions
		conditions_2[duplicated(conditions_2)] <- ""
		scores <- as.data.frame(pca$x)
		PCs.combinations <- combn(PCs,2)
		plots <- apply(
			PCs.combinations,
			2,
			function(combination)
			{
				data <- scores[,c(combination[1], combination[2])]
				data$cond <- conditions
				colnames(data) <- c("PC_x", "PC_y", "cond")
				plot <- ggplot(data, aes(x=PC_x, y=PC_y, fill=cond)) +
					geom_point(shape = 21, size = 6, stroke=0.5) +
					scale_fill_manual(values = conditions_color) +
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

corr_plot <- grid.grabExpr(correlation(TPM, conditions, corr_method, conditions_color))

pca_plot <- plot.pca(TPM, conditions, conditions_color, c("PC1", "PC2"))

figure <- plot_grid(
	plotlist=list(corr_plot, pca_plot[[1]]),
	labels = "AUTO",
	ncol=2
)

save_plot(
	snakemake@output[['pdf']], 
	figure, 
	base_width=35,
	base_height=13.5,
	units = c("cm"), 
	dpi=300
)

save_plot(
	snakemake@output[['png']], 
	figure, 
	base_width=35,
	base_height=13.5,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)