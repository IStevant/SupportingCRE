source(".Rprofile")
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
})

doParallel::registerDoParallel(cores=4)

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

# Overide monaLisa default colors fir bins, i.e. for conditions
my.getColsByBin <- function(b,
                         col1 = c("#FFB100", "#FFB100", "#FFB100", "#FFB100", "#FFB100"),
                         col2 = c("#339989", "#339989", "#339989", "#339989", "#339989"),
                         col0 = "#F5F5F5") {
	print("###########################")
    if (!is.factor(b)) {
        b <- factor(b, levels = unique(b))
        b <- setZeroBin(b, NA)
    }

    if (!is.null(getZeroBin(b)) && !is.na(getZeroBin(b))) {
        bin0 <- getZeroBin(b)
        cols <- c(colorRampPalette(col1)(bin0 - 1L),
                  "#AAAAAA33",
                  colorRampPalette(col2)(nlevels(b) - bin0))
    } else {
        nh <- round(nlevels(b) / 2)
        cols <- c(colorRampPalette(col1)(nh),
                  colorRampPalette(col2)(nlevels(b) - nh))
    }

    res <- cols[b]
    names(cols) <- levels(b)
    attr(res, "cols") <- cols
    return(res)
}

unlockBinding("getColsByBin", as.environment("package:monaLisa"))
assign("getColsByBin", my.getColsByBin, "package:monaLisa")
lockBinding("getColsByBin", as.environment("package:monaLisa"))

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
TFs <- as.vector(read.csv("data/mouse_transcription_factors.txt", header=FALSE)[,1])

# Minimum TPM value to considere a gene expressed
minTPM <- snakemake@params[['minTPM']]

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
	# Select the genes expressed at a specific stage for both sexes
	genes <- TPM[,grep(stg, colnames(TPM))]
	# Discard lowly expressed genes
	genes <- run_filter_low_counts(genes, minTPM)
	genes <- rownames(genes[rowSums(genes)>0,])
	# Select only the TFs
	stg_TFs <- genes[which(genes %in% TFs)]

	write.csv(stg_TFs, file=paste0("results/RNA_", stg, "_expressed_TFs.csv"))

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
	se2 <- calcBinnedMotifEnrR(
		seqs = sequences, 
		bins = bins,
		pwmL = pwms
	)
	# Select the motifs enriched with a -Log10Padj > 10
	sel2 <- apply(SummarizedExperiment::assay(se2, "negLog10Padj"), 1, 
				function(x) max(abs(x), 0, na.rm = TRUE)) > 10.0

	# Calculate similarities between the motifs to cluster them
	seSel <- se2[sel2, ]
	SimMatSel <- motifSimilarity(SummarizedExperiment::rowData(seSel)$motif.pfm)
	# Cluster the motifs by similarity
	hcl <- hclust(as.dist(1 - SimMatSel), method = "average")

	heatmap <- grid.grabExpr(
		plotMotifHeatmaps(
			x = seSel, 
			which.plots = c("log2enr"), 
			width = 1, 
			cluster = hcl, 
			maxEnr = 2, 
			maxSig = 100,
			show_dendrogram = TRUE, 
			show_seqlogo = FALSE
			# width.seqlogo = 1.2
		)
	)

	figure <- plot_grid(
		heatmap,
		labels = stg,
		ncol=1
	)

	return(figure)
}

figure <- plot_grid(
	plotlist=enrichments,
	labels = "AUTO",
	ncol=4
)

save_plot(
	snakemake@output[['pdf']], 
	figure,
	base_width=48,
	base_height=68,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)
