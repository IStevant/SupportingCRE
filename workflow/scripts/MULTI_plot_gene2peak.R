source(".Rprofile")
source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################
suppressPackageStartupMessages({
	library("Gviz")
	library("GenomicInteractions")
	library("InteractionSet")
	library("dplyr")
	library('doParallel')
	library('foreach')
	library("cowplot")
	library("grid")
	library("ggplotify")
})

doParallel::registerDoParallel(cores=12)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

import_bw_files <- function(folder, files, locus, window){

	locus <- resize(gene_TSS[gene_TSS$name %in% locus,], width=2, fix='start')
	full_locus <- locus
	start(full_locus) <- start(full_locus) - window[1]
	end(full_locus) <- end(full_locus) + window[2]

	imports <- foreach(file=files, .final = function(file) setNames(file, files)) %dopar% {
		sex <- sapply(strsplit(file, "_"), `[`, 2)
		stage <- sapply(strsplit(file, "_"), `[`, 1)
		gr <- rtracklayer::import(
			paste(folder, file, sep="/"),
			which=full_locus
		)
		if(length(gr)<1){
			gr <- locus
			gr$score <- 0.1
		}
		gr$cond <- rep(paste0(sex, " ", stage), length(gr))
		return(gr)
	}
	sex <- unique(sapply(strsplit(files, "_"), `[`, 2))
	stages <- unique(sapply(strsplit(files, "_"), `[`, 1))
	gr_list <- lapply(sex, function(sx){
			sex_gr <- imports[grep(sx, names(imports))]
			lapply(stages, function(stg){
				stg_gr <- sex_gr[grep(stg, names(sex_gr))]
			})
		})
	names(gr_list) <- sex
	names(gr_list[[1]]) <- stages
	names(gr_list[[2]]) <- stages
	return(gr_list)
}

generate_genomic_tracks <- function(gr_list, locus, window, max_score, colors){

	locus <- resize(gene_TSS[gene_TSS$name %in% locus,], width=2, fix='start')

	cond <- stringr::str_extract(names(gr_list), "^.{8}")
	color <- colors[cond]
	res=2000
	gTrack_list <- lapply(
		gr_list,
		function(gr){
			DataTrack(
				range = gr, 
				type = "hist",
				baseline=0,
				window=res,
				chromosome = as.character(seqnames(locus)),
				name = unique(gr$cond),
				col.baseline=color,
				col.histogram=color,
				fill.histogram=color,
				ylim=c(0, trunc(max_score, digit=4)),
				yTicksAt=c(0,trunc(max_score, digit=4)),
				lwd=0,
				alpha=0.8
			)
		}
	)
	track <- OverlayTrack(gTrack_list)
	return(track)
}


plot_tracks <- function(peak_gr, TxDb, gene2symbol, link, plot_list, locus, window){

	locus_TSS <- resize(gene_TSS[gene_TSS$name %in% locus,], width=2, fix='start')

	full_locus <- locus_TSS
	start(full_locus) <- start(full_locus) - window[1]
	end(full_locus) <- end(full_locus) + window[2]

	locus_gr <- full_locus
	peak_gr_locus <- subsetByOverlaps(peak_gr, locus_gr)
	gene_links <- links[links$Gene %in% locus,]
	print(paste(nrow(gene_links), "interactions found."))

	# If links exist, plot link track
	if(nrow(gene_links)>=1){
		cor_type <- gene_links$correlations
		cor_type[cor_type>0] <- "Positive"
		cor_type[cor_type<0] <- "Negative"
		colors <- cor_type
		colors[colors=="Positive"] <- "#EF6351"
		colors[colors=="Negative"] <- "#2191FB"

		link_peak_gr <- GRanges(gene_links$Peak)
		link_peak_gr$type <- cor_type
		link_peak_gr$colors <- colors
		link_gene_TSS <- GRanges(seqname=seqnames(link_peak_gr), ranges = IRanges::IRanges(start=gene_links$TSS, end=gene_links$TSS+1))
		link_gene_TSS$type <- rep("Promoter", length(link_gene_TSS))
		interactions <- GenomicInteractions(link_peak_gr, link_gene_TSS, counts=abs(gene_links$correlations))
		annotateInteractions(
			interactions, 
			list(
				positive=link_peak_gr[link_peak_gr$type=="Positive",], 
				negative=link_peak_gr[link_peak_gr$type=="Negative",]
			)
		)
		# interactions$correlation <- gene_links$correlations
		interaction_track <- InteractionTrack(interactions, name = "Links")

		displayPars(interaction_track) <- list(
			plot.anchors=FALSE,
			# col.interactions="red",
			# interaction.measure = "counts",
			# interaction.dimension="width",
			# plot.trans=FALSE,
			plot.outside = TRUE,
			# plot.outside = FALSE,
			col.interaction.types= c('positive-distal'="#EF6351", 'negative-distal'="#2191FB"),
			# col.outside="grey",
			rotation.title=0,
			legend=TRUE
		)
	}

	genome_track <- GenomeAxisTrack(
		col = "#333333", 
		fontsize=13,
		lwd=1,
		distFromAxis=0.5,
		labelPos = "alternating",
		sizes=0.3
	)

	gene_track <- GeneRegionTrack(
		TxDb,
		collapseTranscripts = "meta",
		name = "Genes",
		chromosome = as.character(seqnames(locus_gr)),
		start = start(locus_gr),
		end = end(locus_gr),
		fontface.group="italic",
		col = NULL,
		col.line = NULL,
		fill = "#585858",
		# lwd=0.3,
		fontcolor.group= "#333333",
		fontsize.group=18,
		sizes=0.3,
		rotation.title=0,
		thinBoxFeature="UTR"
		# just.group="above"
	)
	ranges(gene_track)$symbol <- gene2symbol[ranges(gene_track)$gene, "gene_name"]

	peak_track <- AnnotationTrack(
		peak_gr_locus,
		chromosome = as.character(seqnames(peak_gr_locus)),
		start = start(peak_gr_locus),
		end = end(peak_gr_locus),
		strand = as.character(strand(peak_gr_locus)),
		name = "OCRs",
		col.line = NULL,
		col = NULL,
		fill = "#5f4780",
		sizes=0.2,
		rotation.title=0
	)


	visible_linked_peaks <- subsetByOverlaps(link_peak_gr, locus_gr)
	promoter <- locus_TSS+500

	if(length(visible_linked_peaks)>=1){
		ht <- HighlightTrack(
			trackList = c(peak_track, gene_track, plot_list),
			start = c(start(visible_linked_peaks), start(promoter)),
			end = c(end(visible_linked_peaks), end(promoter)),
			col = c(visible_linked_peaks$colors, "#000000"),
			fill = c(visible_linked_peaks$colors, "#000000"),
			chromosome = as.character(unique(seqnames(visible_linked_peaks))),
			alpha=0.2,
			lwd=0.5
		)

		plot <- plotTracks(
			c(interaction_track, ht, genome_track),
			chromosome = as.character(seqnames(locus_gr)),
			from = start(locus_gr),
			to = end(locus_gr),
			transcriptAnnotation="symbol",
			background.title = "transparent",
			col.border.title="transparent",
			col.title = "#333333",
			col.axis = "#333333",
			alpha.title=1,
			alpha.axis=1,
			sizes=c(0.3, 0.2, 0.4, rep(0.4, length(plot_list)), 0.3),
			cex.title=0.9,
			title.width=1
		)

	} else if(nrow(gene_links)<1) {
		ht <- HighlightTrack(
			trackList = c(peak_track, gene_track, plot_list),
			start = start(promoter),
			end = end(promoter),
			col = "#000000",
			fill = "#000000",
			chromosome = as.character(seqnames(promoter)),
			alpha=0.2,
			lwd=0.5
		)

		plot <- plotTracks(
			c(ht, genome_track),
			chromosome = as.character(seqnames(locus_gr)),
			from = start(locus_gr),
			to = end(locus_gr),
			transcriptAnnotation="symbol",
			background.title = "transparent",
			col.border.title="transparent",
			col.title = "#333333",
			col.axis = "#333333",
			alpha.title=1,
			alpha.axis=1,
			sizes=c(0.3, 0.2, 0.4, rep(0.4, length(plot_list)), 0.3),
			cex.title=0.9,
			title.width=1
		)
	} else {
		ht <- HighlightTrack(
			trackList = c(peak_track, gene_track, plot_list),
			start = start(promoter),
			end = end(promoter),
			col = "#000000",
			fill = "#000000",
			chromosome = as.character(seqnames(promoter)),
			alpha=0.2,
			lwd=0.5
		)

		plot <- plotTracks(
			c(interaction_track, ht, genome_track),
			chromosome = as.character(seqnames(locus_gr)),
			from = start(locus_gr),
			to = end(locus_gr),
			transcriptAnnotation="symbol",
			background.title = "transparent",
			col.border.title="transparent",
			col.title = "#333333",
			col.axis = "#333333",
			alpha.title=1,
			alpha.axis=1,
			sizes=c(0.3, 0.2, 0.4, rep(0.4, length(plot_list)), 0.3),
			cex.title=0.9,
			title.width=1
		)
	}

	return(plot)
}

###########################################
###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

bw_folder <- snakemake@params[['bw_folder']]
genome <- snakemake@input[['genome']]
peaks <- snakemake@input[['peaks']]
links <- snakemake@input[['linkage']]
gene_bed <- snakemake@input[['gene_bed']]
gene_list <- snakemake@input[['gene_list']]

# bw_folder <- "results/processed_data/mm10/bigWig"
# # genome <- "workflow/data/mm10/iGenome_mm10_ucsc_genes.gtf.gz"
# genome <- "workflow/data/mm10/gencode.vM25.annotation.gtf.gz"
# peaks <- "results/processed_data/mm10/ATAC_norm_counts.csv"
# links <- "results/tables/mm10/all_sig_gene2peak_linkage.csv"
# gene_bed <- "workflow/data/mm10/gene_standard.bed"
# gene_list <- "workflow/data/peak2gene_query.tsv"

bw_files <- list.files(path=bw_folder, pattern = "REP..bw")
links <- read.table(links, header=TRUE, row.names=1)
gene_list <- read.table(gene_list, header=TRUE)
window <- c(gene_list$Start, gene_list$End)
genes <- gene_list$Gene

peak_gr <- GRanges(rownames(read.csv(peaks, header=TRUE, row.names=1)))
gene_TSS <- GRanges(rtracklayer::import(gene_bed))

conditions <- paste(
	sapply(strsplit(bw_files, "_"), `[`, 1),
	sapply(strsplit(bw_files, "_"), `[`, 2),
	sep="_"
)
conditions <- unique(conditions)
conditions <- c(conditions[c(TRUE, FALSE)], conditions[c(FALSE, TRUE)])
names(conditions_color) <- conditions

genome_gtf <- rtracklayer::import(genome)
gene2symbol <- mcols(genome_gtf)[,c("gene_id","gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

TxDb <- GenomicFeatures::makeTxDbFromGFF(genome)


pdf(file=snakemake@output[['pdf']], width=6, height=8)
# pdf(file="test.pdf", width=6, height=8)

plot <- lapply(genes, function(gene){

	locus <- gene_list[gene_list$Gene == gene,]
	window <- c(locus$Start, locus$End)*1000
	max_score <- locus$MaxScore

	atac_coverage <- import_bw_files(
		bw_folder, 
		bw_files,
		gene, 
		window
	)

	if(max_score==0){
		max_score <- do.call("c", do.call("c", do.call("c", atac_coverage)))
		max_score <- max(unlist(lapply(max_score, function(x) x$score)))
	}

	XX_track_list <- lapply(
		atac_coverage[["XX"]], 
		function(stage) 
			generate_genomic_tracks(
				stage, 
				gene, 
				window, 
				max_score, 
				conditions_color
			)
		)

	XY_track_list <- lapply(
		atac_coverage[["XY"]], 
		function(stage) 
			generate_genomic_tracks(
				stage, 
				gene, 
				window, 
				max_score, 
				conditions_color
			)
		)


	plot <- plot_tracks(
		peak_gr, 
		TxDb,
		gene2symbol,
		links, 
		c(XX_track_list, XY_track_list), 
		gene, 
		window
	)

	return(plot)

})

dev.off()
