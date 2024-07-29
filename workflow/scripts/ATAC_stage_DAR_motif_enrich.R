source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
	library("monaLisa")
	if (snakemake@params[['genome']]=="mm10"){
		library("BSgenome.Mmusculus.UCSC.mm10")
	} else {
		library("BSgenome.Mmusculus.UCSC.mm39")
	}
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
	library("motifStack")
})

doParallel::registerDoParallel(cores=12)

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

get_enriched_TFs <- function(save_folder){
	# sx <- "E11.5"
	# save_folder <- "."
	# Select the genes expressed at a specific stage for both sexes
	genes <- TPM
	# genes <- TPM[,grep("XX", colnames(genes))]
	# Discard lowly expressed genes
	genes <- run_filter_low_counts(genes, minTPM)
	genes <- rownames(genes[rowSums(genes)>0,])
	# Select only the TFs
	# stg_TFs <- genes[which(genes %in% TFs)]

	pwms <- TFBSTools::getMatrixSet(
		JASPAR,
		opts = list(
			matrixtype = "PWM",
			tax_group = "vertebrates"
			# species = "Mus musculus"
		)
	)

	peaks <- filtered_StageDARs
	clusters <-peaks$x

	# generate GRanges objects
	peak_gr <- GRanges(rownames(peaks)) 

	# Get the peak sequences
	if (genome_version=="mm10"){
		sequences <-  Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, peak_gr)
	} else {
		sequences <-  Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm39, peak_gr)
	}

	# # generate GRanges objects
	# female <- GenomicRanges::GRanges(rownames(peaks[peaks$Diff.Acc.=="More in XX",]))
	# male <- GenomicRanges::GRanges(rownames(peaks[peaks$Diff.Acc.=="More in XY",]))
	# all <- c(female, male)

	# Define which sequences are male or female specific
	bins <- clusters
	bins <- factor(bins)
	table(bins)

	# Calculate motif enrichments
	# If background is "genome", run the analysis against radom genomic regions
	# Else, compare the two sexes
	if (background=="genome"){
		if (genome_version=="mm10"){
			se2 <- calcBinnedMotifEnrR(
				seqs = sequences, 
				bins = bins,
				pwmL = pwms,
				background = "genome",
				genome = BSgenome.Mmusculus.UCSC.mm10,
				genome.regions = NULL, # sample from full genome
				genome.oversample = 2,
				BPPARAM = BiocParallel::MulticoreParam(12)
			)
		} else {
			se2 <- calcBinnedMotifEnrR(
				seqs = sequences, 
				bins = bins,
				pwmL = pwms,
				background = "genome",
				genome = BSgenome.Mmusculus.UCSC.mm39,
				genome.regions = NULL, # sample from full genome
				genome.oversample = 2,
				BPPARAM = BiocParallel::MulticoreParam(12)
			)
		}
	} else {
		se2 <- calcBinnedMotifEnrR(
			seqs = sequences, 
			bins = bins,
			pwmL = pwms,
			BPPARAM = BiocParallel::MulticoreParam(12)
		)
	}

	# Get gene expression of the TFs
	TF_names <- tolower(as.vector(rowData(se2)$motif.name))
	names(TF_names) <- rowData(se2)$motif.id
	genes <- TPM[,grep(sex,colnames(TPM))]
	genes <- run_filter_low_counts(genes, minExp=5)
	rownames(genes) <- tolower(rownames(genes))

	TF_exp <- matrix(0, nrow=length(TF_names), ncol=ncol(genes))
	rownames(TF_exp) <- TF_names
	colnames(TF_exp) <- colnames(genes)

	for (TF in TF_names) {
		if (TF %in% rownames(genes)){
			TF_exp[TF,] <- unlist(genes[TF,])
		}
	}

	TF_expression <- data.frame(
		a=rowMeans(TF_exp[,grep(sex,colnames(TF_exp))]),
		b=rowMeans(TF_exp[,grep(sex,colnames(TF_exp))]),
		c=rowMeans(TF_exp[,grep(sex,colnames(TF_exp))]),
		d=rowMeans(TF_exp[,grep(sex,colnames(TF_exp))])
	)

	rownames(TF_expression) <- rownames(SummarizedExperiment::assay(se2, "negLog10Padj"))

	assays(se2)$expr <- TF_expression

	# Select the motifs enriched with a -Log10Padj > 10
	sel2 <- apply(SummarizedExperiment::assay(se2, "negLog10Padj"), 1, 
				function(x) max(abs(x), 0, na.rm = TRUE)) > 10.0
	seSel <- se2[sel2, ]

	# Select TFs that are expressed
	sel2 <- apply(SummarizedExperiment::assay(seSel, "expr"), 1, 
			function(x) max(x, 0, na.rm = TRUE)) > 5
	seSel <- seSel[sel2, ]


	# sel2 <- apply(SummarizedExperiment::assay(seSel, "expr"), 1, 
	# 		function(x) max(x, 0, na.rm = TRUE)) > 5
	# seSel <- seSel[sel2, ]
	# seSel@elementMetadata$motif.name <- toupper(seSel@elementMetadata$motif.name)
	# assays(seSel)$expr <- log10(assays(seSel)$expr)

	TF_summary <- data.frame(
		SummarizedExperiment::assay(seSel, "log2enr"),
		TF.name=seSel@elementMetadata$motif.name
	)

	write.csv(TF_summary, file=paste0(save_folder, "/ATAC_stage_DAR_TF_", sex, "_", background,"_bg.csv"))


	return(seSel)
}


merge_TF_motifs <- function(seSel){
	# seSel <- enrichments
	# Cluster motifs by enrichment
	TF_enrichment <- SummarizedExperiment::assay(seSel, "log2enr")
	# rownames(TF_enrichment) <- seSel@elementMetadata$motif.name

	nbCluster=4

	hcl <- hclust(dist(TF_enrichment))
	clustering <- cutree(hcl, k=nbCluster)

	motif_sig <- lapply(1:nbCluster, function(cl){
		TFs <- names(clustering[clustering==cl])
		# print(TFs)
		if(length(TFs)>1){
			# TFs <- names(clustering[clustering==cl])
			# print(TFs)
			matrices <- seSel@elementMetadata$motif.pfm[TFs]
			# print(length(matrices))
			pfms <- universalmotif::convert_motifs(matrices, class="motifStack-pcm")
			hc <- motifStack::clusterMotifs(pfms)
			phylog <- ade4::hclust2phylog(hc)
			# extract the motif signatures
			motifSig <- motifSignature(pfms, phylog, cutoffPval = 0.005, min.freq=1)
			## get the signatures from object of motifSignature
			sig <- signatures(motifSig)	
		} else {
			sig <- universalmotif::convert_motifs(seSel@elementMetadata$motif.pfm[TFs], class="motifStack-pcm")
		}
		return(sig)
	})

	motif_sig_enr <- lapply(motif_sig, function(cl){
		lapply(cl, function(motif){
			names <- motif@name
			lapply(names, function(TF){
				tf_vector <- toupper(unlist(strsplit(TF, ";")))
				enrichment <- TF_enrichment
				rownames(enrichment) <- toupper(seSel@elementMetadata$motif.name)
				tf_enrichment <- enrichment[tf_vector, , drop=FALSE]
				if(nrow(tf_enrichment)>1){

					extract_prefix <- function(gene) {
						prefix <- sub("([a-zA-Z]+).*", "\\1", gene)
						if(prefix=="NR"){
							return(gene)
						} else if(grepl("^HOX", prefix)){
							prefix <- stringr::str_sub(prefix, end=-2)
							return(prefix)
						} else {
							return(prefix)
						}
					}

					grouped_genes <- split(tf_vector, sapply(tf_vector, extract_prefix))

					new_tf_vector <- sapply(grouped_genes, function(gene_group) {
						common_prefix <- extract_prefix(gene_group[1])
						if(common_prefix!="HOX"){
							suffixes <- gsub(paste0("^", common_prefix), "", gene_group)
							numeric_suffixes <- sort(as.numeric(suffixes[grepl("^\\d+$", suffixes)]), na.last = TRUE)
							non_numeric_suffixes <- sort(suffixes[!grepl("^\\d+$", suffixes)])
							all_suffixes <- c(non_numeric_suffixes, numeric_suffixes)
							concatenated_suffixes <- paste(all_suffixes, collapse = "/")
							paste0(common_prefix, concatenated_suffixes)
						} else {
							paste0(common_prefix,"s")
						}

					})

					new_tf_vector <- new_tf_vector[order(new_tf_vector)]
					paste(new_tf_vector, collapse=";")
					tf_names <- paste(new_tf_vector, collapse=";")

					tf_enrichment <- as.data.frame(colMedians(tf_enrichment))
					colnames(tf_enrichment) <- tf_names
					tf_enrichment <- as.data.frame(t(tf_enrichment))
					tf_enrichment$mat_names <- TF	
					# print(tf_enrichment)
				} else {
					tf_enrichment <- as.data.frame(tf_enrichment)
					tf_enrichment$mat_names <- TF
				}

				return(tf_enrichment)
			})
		})
	})

	enrichment <- do.call("rbind", do.call("rbind", do.call("rbind",motif_sig_enr)))
	enrichment <- enrichment[!duplicated(enrichment), ]
	
	merged_pfms <- do.call("c", do.call("c", do.call("c",motif_sig)))
	names(merged_pfms) <- unlist(lapply(merged_pfms, function(mat) mat@name))

# }

# stg <- "E15.5"

# plot_heatmap <- function(enrichment, merged_pfms, stg){

	# if(background=="genome"){
	# 	nbCluster=3
	# } else {
	# 	nbCluster=2
	# }

	motifs <- merged_pfms[enrichment$mat_names]
	motifs_pfms <- universalmotif::convert_motifs(motifs, class="TFBSTools-PFMatrix")

	maxwidth <- max(unlist(lapply(motifs_pfms, function(x) ncol(x@profileMatrix))))

	grobL <- lapply(motifs_pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")

	names(grobL) <- rownames(enrichment)

	matrix <- enrichment[,-5]
	matrix[matrix>1] <- 1
	matrix[matrix<(-1)] <- (-1)

	hmSeqlogo <- HeatmapAnnotation(
		logo = annoSeqlogo(
			grobL = grobL, which = "row",
			space = unit(0.5, "mm"),
			width = unit(1.5, "inch")
		),
		show_legend = FALSE, 
		show_annotation_name = FALSE, 
		which = "row"
	)

	bincols <- rep("black", ncol(matrix))
	names(bincols) <- colnames(matrix)
	conditions <- colnames(matrix)

	stage_anno <- HeatmapAnnotation(
		Stages = anno_block(
			gp = gpar(fill = bincols, col = 0), 
			labels = names(bincols),
			labels_gp = gpar(col = "white", fontsize = 14),
			height = unit(7, "mm")
		)
	)

	cold <- colorRampPalette(c('#04bbc6','#52d0cf','#8addd8','#c1f0e0',"#fffee8"))
	warm <- colorRampPalette(c("#fffee8",'#ffd9cb','#ffb1ad','#f5808d','#eb2d62'))
	TYP <- c(cold(12), warm(12))

	mypalette <- colorRamp2(breaks = seq(-1, 1, length.out = 24),
							 colors = TYP)

	# mypalette <- RColorBrewer::brewer.pal(11,"BuPu")

	if(nrow(matrix)<10) {
		cell_height <-10
	} else {
		cell_height <- 6
	}

	ht_list <-	Heatmap(
			matrix,
			clustering_method_rows = 'ward.D2',
			name = "Log2 enrichment",
			# row_km = nbCluster,
			right_annotation = hmSeqlogo,
			top_annotation = stage_anno,
			column_split = conditions,
			show_column_names = FALSE,
			show_row_dend = FALSE,
			# cluster_row_slices = FALSE,
			cluster_columns = FALSE,
			column_title = NULL,
			row_title = NULL,
			col = mypalette,
			width = ncol(matrix)*unit(15, "mm"),
			height = nrow(matrix)*unit(cell_height, "mm"),
			row_names_max_width = unit(17, "cm"),
			heatmap_legend_param = list(direction = "horizontal")
		)

	ht_list_2 <- grid.grabExpr(draw(ht_list, heatmap_legend_side = "bottom"), wrap.grobs = TRUE)

	# pdf("test.pdf", height=5, width=15)
	# 	plot(ht_list)
	# dev.off()

	return(ht_list_2)
}



#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# filtered_StageDARs
filtered_StageDARs <- read.csv(file=snakemake@input[['sig_DARs']], header=TRUE, row.names=1)
# Load RNA-seq TPM matrix to filter the TFs that are expressed in the gonads
TPM <- read.csv(file=snakemake@input[['TPM']], header=TRUE, row.names=1)
# Load the mouse TFs list
TFs <- as.vector(read.csv(snakemake@input[['TF_genes']], header=FALSE)[,1])
# Minimum TPM value to considere a gene expressed
minTPM <- snakemake@params[['minTPM']]
# Run analysis using the genome bakground or calculating the enrichment compared to the conditions
background <- snakemake@params[['background']]
save_folder <- snakemake@params[['save_folder']]
# Print the logos of the TFs on the heatmap
logos <- snakemake@params[['logos']]
genome_version <- snakemake@params[['genome']]
# Nb of top TFs per sex to print 
nbTFs <- snakemake@params[['nbTFs']]
sex <- snakemake@params[['sex']]



# filtered_StageDARs <- read.csv(file="results/tables/mm10/ATAC_XX_DAR_stage_heatmap_clusters.csv", header=TRUE, row.names=1)
# TPM <- read.csv(file="results/processed_data/mm10/RNA_TPM.csv", header=TRUE, row.names=1)
# TFs <- as.vector(read.csv("workflow/data/mouse_transcription_factors.txt", header=FALSE)[,1])
# minTPM <- 5
# background <- "conditions"
# # background <- "conditions"
# save_folder <- "."
# logos <- TRUE
# genome_version <- "mm10"
# nbTFs <- 5
# sex="XX"


# if (logos=="TRUE") {
# 	width <- 70
# } else {
# 	width <- 50
# }

stage <- sapply(strsplit(colnames(TPM), "_"), `[`, 1)
sx <- sapply(strsplit(colnames(TPM), "_"), `[`, 2)
conditions <- unique(paste(sex, stage, sep=" "))
names(conditions_color) <- conditions[order(conditions)]
XX_colors <- conditions_color[grepl("XX" , names(conditions_color))]
XY_colors <- conditions_color[grepl("XY" , names(conditions_color))]

TPM <- TPM[,grep(sex,colnames(TPM))]
conditions_color <- conditions_color[grep(sex, names(conditions_color))]

# Load Jaspar 2024 database
JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db

###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

stages <- unique(stage)
sexes <- unique(sex)
# stages <- "E11.5"

# For each stage, get TFBS motif enrichments
enrichments <- get_enriched_TFs(save_folder)


# stg <- "E15.5"
# enrichments <- list(get_enriched_TFs(stg))

figure <- merge_TF_motifs(enrichments)

# figure <- plot_grid(
# 	plotlist=heatmap_list,
# 	labels = "AUTO",
# 	ncol=4
# )

save_plot(
	snakemake@output[['pdf']],
	# "test.pdf",
	figure,
	base_width=30,
	base_height=30,
	units = c("cm"), 
	dpi=300
)

save_plot(
	snakemake@output[['png']],
	figure,
	base_width=30,
	base_height=30,
	units = c("cm"), 
	dpi=300
)
