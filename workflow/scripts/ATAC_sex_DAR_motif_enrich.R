source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
	library("monaLisa")
	library("BSgenome.Mmusculus.UCSC.mm10")
	library('doParallel')
	library('foreach')
	library('dplyr')
	library("cowplot")
	library("grid")
	library("ggplot2")
	library("gtable")
	library("SummarizedExperiment")
	library("ComplexHeatmap")
	library("circlize")
})

doParallel::registerDoParallel(cores=8)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' When the maximum value (read count) of a peak between samples is under a certain threshold, we considere it is not relevant and the values are set to 0.
#' @param data Read count matrix.
#' @param minExp Minimum value. Default is 5.
#' @return Return a dataframe.
run_filter_low_counts <- function(data, minExp=5) {
	col_names <- colnames(data)
	data <- t(apply(data, 1, filter_low_counts, col_names = col_names, minExp = minExp))
	return(as.data.frame(data))
}
filter_low_counts <- function(row, col_names, minExp) {
	if (max(row) < minExp) {
		return(setNames(rep(0, length(row)), col_names))
	} else {
		return(setNames(row, col_names))
	}
}

# Overide monaLisa heatmap plots
.assertScalar <- function(x,
						  type = NULL,
						  rngIncl = NULL,
						  rngExcl = NULL,
						  validValues = NULL) {
	args <- lapply(sys.call()[-1], as.character)
	xname <- if ("x" %in% names(args)) args$x else "argument"

	if (length(x) != 1L) {
		stop("'", xname, "' must be a scalar value (length one)")
	}
	
	if (is.null(type) && (!is.null(rngIncl) || !is.null(rngExcl))) {
		type <- "numeric"
	}

	if (!is.null(type) && !is(x, type)) {
		stop("'", xname, "' must be of type '", type, "'")
	}
	
	if (!is.null(rngIncl) && is.numeric(rngIncl) && length(rngIncl) == 2L &&
		(x < rngIncl[1] || x > rngIncl[2])) {
		stop("'", xname, "' must be within [", rngIncl[1], ",", 
			 rngIncl[2], "] (inclusive)")
	}
	
	if (!is.null(rngExcl) && is.numeric(rngExcl) && length(rngExcl) == 2L &&
		(x <= rngExcl[1] || x >= rngExcl[2])) {
		stop("'", xname, "' must be within (", rngExcl[1], ",", 
			 rngExcl[2], ") (exclusive)")
	}
	
	if (!is.null(validValues) && !(x %in% validValues)) {
		stop("'", xname, "' must be one of: ", paste(validValues, 
													 collapse = ", "))
	}
	
	return(invisible(TRUE))
}



plotMotifHeatmaps_exp <- function(x,
							  which.plots = c("negLog10P", "pearsonResid", 
											  "negLog10Padj", "log2enr", "expr"),
							  width = 4,
							  col.enr = c("#053061", "#2166AC", "#4393C3",
										  "#92C5DE", "#D1E5F0", "#F7F7F7",
										  "#FDDBC7", "#F4A582", "#D6604D",
										  "#B2182B", "#67001F"),
							  col.sig = RColorBrewer::brewer.pal(8, "BuGn"),
							  col.gc = c("#F7FCF5", "#E5F5E0", "#C7E9C0",
										 "#A1D99B", "#74C476", "#41AB5D",
										 "#238B45", "#006D2C", "#00441B"),
							  bincols = NULL,
							  maxEnr = NULL,
							  maxSig = NULL,
							  highlight = NULL,
							  cluster = FALSE,
							  show_dendrogram = FALSE,
							  show_motif_GC = FALSE,
							  show_seqlogo = FALSE,
							  width.seqlogo = 1.5,
							  use_raster = FALSE,
							  na_col = "white", 
							  doPlot = TRUE,
							  ...) {
	stopifnot(exprs = {
		is(x, "SummarizedExperiment")
		all(which.plots %in% assayNames(x))
		"bins" %in% names(metadata(x))
		(!show_motif_GC || "motif.percentGC" %in% colnames(rowData(x)))
	})
	b <- metadata(x)$bins
	.assertScalar(x = width, type = "numeric", rngExcl = c(0, Inf))
	.assertScalar(x = show_dendrogram, type = "logical")
	.assertScalar(x = show_motif_GC, type = "logical")
	.assertScalar(x = show_seqlogo, type = "logical")
	.assertScalar(x = width.seqlogo, type = "numeric", rngExcl = c(0, Inf))
	.assertScalar(x = use_raster, type = "logical")
	.assertScalar(x = na_col, type = "character")
	.assertScalar(x = doPlot, type = "logical")
	stopifnot(exprs = {
		ncol(x) == nlevels(b)
		all(which.plots %in% c("negLog10P", "negLog10Padj", 
							   "pearsonResid", "log2enr", "expr"))
		is.null(highlight) || (is.logical(highlight) && 
								 length(highlight) == nrow(x))
	})
	# bincols <- attr(getColsByBin(b), "cols")
	if (identical(cluster, TRUE)) {
		clAssayName <- "pearsonResid"
		clAssay <- assay(x, clAssayName)
		allNA <- rowSums(is.na(clAssay)) == ncol(clAssay)
		if (any(allNA)) {
			warning("removing motifs without finite values in '",
					clAssayName, "': ",
					paste(rownames(clAssay)[allNA], collapse = ", "))
			x <- x[!allNA, ]
			clAssay <- clAssay[!allNA, ]
		}
		clres <- hclust(dist(clAssay))
	} else if (identical(cluster, FALSE)) {
		clres <- FALSE
	} else if (is(cluster, "hclust")) {
		clres <- cluster
	} else {
		stop("'cluster' must be either TRUE, FALSE or an hclust-object.")
	}
	hmBin <- HeatmapAnnotation(df = data.frame(bin = colnames(x)), name = "bin",
							   col = list(bin = bincols),
							   show_annotation_name = FALSE,
							   which = "column", width = unit(width,"inch"),
							   annotation_height = unit(width / 16, "inch"),
							   show_legend = FALSE)
	tmp <- matrix(if (!is.null(highlight)) {
		as.character(highlight) 
	} else {
		rep(NA, nrow(x))
	},
	ncol = 1, dimnames = list(unname(rowData(x)$motif.name), NULL))
	hmSeqlogo <- NULL
	if (show_seqlogo) {
		pfms <- rowData(x)$motif.pfm
		maxwidth <- max(vapply(TFBSTools::Matrix(pfms), ncol, 0L))
		grobL <- lapply(pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")
		hmSeqlogo <- HeatmapAnnotation(
			logo = annoSeqlogo(grobL = grobL, which = "row",
							   space = unit(0.5, "mm"),
							   width = unit(width.seqlogo, "inch")),
			# logo=rownames(x),
			show_legend = FALSE, show_annotation_name = FALSE, which = "row")
	}
	hmMotifs <- Heatmap(
		matrix = tmp, name = "names",
		width = unit(if (!is.null(highlight)) .2 else 0, "inch"),
		na_col = NA, col = c("TRUE" = "green3", "FALSE" = "white"),
		cluster_rows = clres, show_row_dend = show_dendrogram,
		cluster_columns = FALSE, show_row_names = TRUE,
		row_names_side = "left", show_column_names = FALSE,
		show_heatmap_legend = FALSE, left_annotation = hmSeqlogo,
		...
	)

	assayNameMap1 <- c(negLog10P = "P value",
					   negLog10Padj = "adj. P value",
					   pearsonResid = "Pearson residual",
					   log2enr = "log2 enrichment",
					   expr = "log10 TPM")
	assayNameMap2 <- c(negLog10P = "P value (-log10)",
					   negLog10Padj = "adj. P value (-log10)",
					   pearsonResid = "Pearson residual (o-e)/sqrt(e)",
					   log2enr = "Motif \nenrichment",
					   expr = "TF\nexpression")
	L <- list(labels = hmMotifs)
	if (show_motif_GC) {
		tmp <- as.matrix(rowData(x)[, "motif.percentGC", drop = FALSE])
		hmPercentGC <- Heatmap(
			matrix = tmp, name = "Percent G+C",
			width = unit(0.2, "inch"), na_col = NA,
			col = colorRamp2(breaks = c(0, seq(20, 80, length.out = 254), 100),
							 colors = colorRampPalette(col.gc)(256)),
			cluster_rows = FALSE, cluster_columns = FALSE,
			show_row_names = FALSE, show_column_names = FALSE,
			show_heatmap_legend = TRUE,
			heatmap_legend_param = list(color_bar = "continuous"),
			use_raster = use_raster,
			...
		)
		L <- c(L, list("percentGC" = hmPercentGC))
	}
	ret <- c(L, lapply(which.plots, function(w) {
		dat <- assay(x, w)
		if ((w == "pearsonResid") | (w == "log2enr") ) {
			rng <- c(-1, 1) * if (is.null(maxEnr)) {
				quantile(abs(dat), .995, na.rm = TRUE) 
			} else {
				maxEnr
			}
			cols <- col.enr
		} else {
			rng <- c(0,4)
			cols <- col.sig
		}
		Heatmap(
			matrix = dat,
			name = assayNameMap1[w],
			width = unit(width,"inch"),
			column_title = assayNameMap2[w],
			col = colorRamp2(breaks = seq(rng[1], rng[2], length.out = 256),
							 colors = colorRampPalette(cols)(256)),
			cluster_rows = FALSE, cluster_columns = FALSE,
			show_row_names = FALSE, show_column_names = FALSE,
			# column_names_side = "bottom", 
			# column_names_max_height = unit(1.5,"inch"),
			top_annotation = hmBin, show_heatmap_legend = TRUE,
			heatmap_legend_param = list(color_bar = "continuous"),
			use_raster = use_raster,
			na_col = na_col, 
			...
		)
	}))
	names(ret)[seq(length(ret) - length(which.plots) + 1L, length(ret))] <- 
		which.plots
	if (doPlot) {
		show(Reduce(ComplexHeatmap::add_heatmap, ret))
	}
	invisible(ret)
}

get_entiched_TFs <- function(stg){
	# Select the genes expressed at a specific stage for both sexes
	genes <- TPM[,grep(stg, colnames(TPM))]
	# genes <- TPM[,grep("XX", colnames(genes))]
	# Discard lowly expressed genes
	genes <- run_filter_low_counts(genes, minTPM)
	genes <- rownames(genes[rowSums(genes)>0,])
	# Select only the TFs
	stg_TFs <- genes[which(genes %in% TFs)]

	# Get human TF matrices
	pwms_human <- TFBSTools::getMatrixSet(
		JASPAR,
		opts = list(
			name=toupper(stg_TFs),
			matrixtype = "PWM",
			# tax_group = "vertebrates"
			species = "Homo sapiens"
		)
	)
	# Get mouse TF matrices
	pwms_mouse <- TFBSTools::getMatrixSet(
		JASPAR,
		opts = list(
			name=stg_TFs,
			matrixtype = "PWM",
			# tax_group = "vertebrates"
			species = "Mus musculus"
		)
	)

	# merge human and mouse lists of TF matrices
	pwms <- c(pwms_human, pwms_mouse)


	# Get sex specific peaks for the given stage
	sex_peaks <- filtered_SexDARs[[stg]]

	# generate GRanges objects
	female <- GenomicRanges::GRanges(rownames(sex_peaks[sex_peaks$Diff.Acc.=="More in XX",]))
	male <- GenomicRanges::GRanges(rownames(sex_peaks[sex_peaks$Diff.Acc.=="More in XY",]))
	all <- c(female, male)

	# Get the peak sequences
	sequences <-  Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, all)

	# Define which sequences are male or female specific
	bins <- rep(c("XX", "XY"), c(length(female), length(male)))
	bins <- factor(bins)
	table(bins)

	# Calculate motif enrichments
	# If background is "genome", run the analysis against radom genomic regions
	# Else, compare the two sexes
	if (background=="genome"){
		se2 <- calcBinnedMotifEnrR(
			seqs = sequences, 
			bins = bins,
			pwmL = pwms,
			background = "genome",
			genome = BSgenome.Mmusculus.UCSC.mm10,
			genome.regions = NULL, # sample from full genome
			genome.oversample = 2,
			BPPARAM = BiocParallel::MulticoreParam(4)
		)
	} else {
		se2 <- calcBinnedMotifEnrR(
			seqs = sequences, 
			bins = bins,
			pwmL = pwms,
			BPPARAM = BiocParallel::MulticoreParam(4)
		)
	}

	# Get gene expression of the TFs
	TF_names <- tolower(as.vector(rowData(se2)$motif.name))
	names(TF_names) <- rowData(se2)$motif.id
	genes <- TPM[,grep(stg, colnames(TPM))]
	rownames(genes) <- tolower(rownames(genes))

	results <- list()

	for (i in seq_along(TF_names)) {
	  result <- genes[rownames(genes) == TF_names[i],]
	  rownames(result) <- names(TF_names[i])
	  results[[names(TF_names[i])]] <- result
	}
	names(results) <- names(TF_names)

	expr <- do.call(rbind.data.frame, results)

	TF_expression <- data.frame(
		XX=rowMeans(expr[,grep("XX",colnames(expr))]),
		XY=rowMeans(expr[,grep("XY",colnames(expr))])
	)

	assays(se2)$expr <- log10(TF_expression)


	# Select the motifs enriched with a -Log10Padj > 10
	sel2 <- apply(SummarizedExperiment::assay(se2, "negLog10Padj"), 1, 
				function(x) max(abs(x), 0, na.rm = TRUE)) > 10.0


	seSel <- se2[sel2, ]

	top_XX <- names(sort(assays(seSel)$negLog10Padj[,"XX"], decreasing=TRUE)[1:nbTFs])
	top_XX <- names(assays(seSel[top_XX,])$log2enr[,"XX"][assays(seSel[top_XX,])$log2enr[,"XX"] > 0])

	top_XY <- names(sort(assays(seSel)$negLog10Padj[,"XY"], decreasing=TRUE)[1:nbTFs])
	top_XY <- names(assays(seSel[top_XY,])$log2enr[,"XY"][assays(seSel[top_XY,])$log2enr[,"XY"] > 0])

	tops <- unique(c(top_XX, top_XY))

	seSel <- seSel[tops,]

	return(seSel)
}

plot_TF_heatmap <- function(seSel, stg){
	# Cluster the motifs by similarity
	# SimMatSel <- motifSimilarity(SummarizedExperiment::rowData(seSel)$motif.pfm)
	# hcl <- hclust(as.dist(1 - SimMatSel), method = "ward.D2")
	hcl <- hclust(dist(assays(seSel)$log2enr), method = "ward.D2")

	bincols <- c(
		XX=XX_colors[grep(stg, names(XX_colors))],
		XY=XY_colors[grep(stg, names(XY_colors))]
	)
	names(bincols) <- c("XX", "XY")

	# If the logo is TRUE, print the heatmaps in a pfd, one heatmap per page
	# If the logo is FALSE, print the heatmap side by side using cowplot
	if (logos=="TRUE") {
		# pdf(file=snakemake@output[['pdf']], width = 20,  height = nbTFs,  units = "cm")
		heatmap <- plotMotifHeatmaps_exp(
			x = seSel,
			bincols = bincols,
			which.plots = c("log2enr", "expr"), 
			show_seqlogo = TRUE,
			width = 1.3, 
			cluster = hcl,
			maxEnr = 1.5, 
			maxSig = 100,
			width.seqlogo = 1.8
		)
	} else {
		heatmap <- grid.grabExpr(
			plotMotifHeatmaps_exp(
				x = seSel,
				bincols = bincols,
				which.plots = c("log2enr", "expr"), 
				width = 1.3, 
				cluster = hcl,
				maxEnr = 1.5, 
				maxSig = 100
			)
		)
		figure <- plot_grid(
			heatmap,
			labels = stg,
			ncol=1
		)
		return(figure)	
	}
}
#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# filtered_SexDARs
load(snakemake@input[['sig_DARs']])

# Load RNA-seq TPM matrix to filter the TFs that are expressed in the gonads
TPM <- read.csv(file=snakemake@input[['TPM']], header=TRUE, row.names=1)

# Load the mouse TFs list
TFs <- as.vector(read.csv(snakemake@input[['TF_genes']], header=FALSE)[,1])

# Minimum TPM value to considere a gene expressed
minTPM <- snakemake@params[['minTPM']]

# Run analysis using the genome bakground or calculating the enrichment compared to the conditions
background <- snakemake@params[['background']]

# Print the logos of the TFs on the heatmap
logos <- snakemake@params[['logos']]

if (logos=="TRUE") {
	width <- 80
} else {
	width <- 50
}

# Nb of top TFs per sex to print 
nbTFs <- snakemake@params[['nbTFs']]

stage <- sapply(strsplit(colnames(TPM), "_"), `[`, 1)
sex <- sapply(strsplit(colnames(TPM), "_"), `[`, 2)
conditions <- unique(paste(sex, stage, sep=" "))
names(conditions_color) <- conditions[order(conditions)]
XX_colors <- conditions_color[grepl("XX" , names(conditions_color))]
XY_colors <- conditions_color[grepl("XY" , names(conditions_color))]

# Load Jaspar 2024 database
JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db

###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

stages <- unique(sapply(strsplit(colnames(TPM), "_"), `[`, 1))

# For each stage, get TFBS motid enrichments
enrichments <- foreach(stg=stages) %dopar% {
	get_entiched_TFs(stg)
}

# stg <- "E15.5"
# enrichments <- list(get_entiched_TFs(stg))

if (logos=="TRUE") {
	pdf(snakemake@output[['pdf']], width=8, height=nbTFs/2)
	# pdf("test.pdf", width=8, height=nbTFs/2)
		for (i in seq_along(enrichments)){ plot_TF_heatmap(enrichments[[i]], stages[i]) }
	dev.off()

} else {

	heatmap_list <- lapply(seq_along(enrichments), function(i) plot_TF_heatmap(enrichments[[i]], stages[i]))

	figure <- plot_grid(
		plotlist=heatmap_list,
		labels = "AUTO",
		ncol=4
	)

	save_plot(
		snakemake@output[['pdf']],
		figure,
		base_width=width,
		base_height=nbTFs,
		units = c("cm"), 
		dpi=300
	)

	save_plot(
		snakemake@output[['png']],
		figure,
		base_width=width,
		base_height=nbTFs,
		units = c("cm"), 
		dpi=300, 
		bg = "white"
	)
}
