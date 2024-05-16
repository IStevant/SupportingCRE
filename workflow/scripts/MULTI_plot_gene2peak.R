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
	# library("BSgenome.Mmusculus.UCSC.mm10")
	# library("TxDb.Mmusculus.UCSC.mm10.knownGene")
	# library("Mus.musculus")
	library("dplyr")
	library('doParallel')
	library('foreach')
})

doParallel::registerDoParallel(cores=12)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

import_bw_files <- function(folder, files, locus, window){

	locus <- resize(gene_TSS[gene_TSS$name %in% locus,], width=2, fix='start')

	imports <- foreach(file=files, .final = function(file) setNames(file, files)) %dopar% {
		sex <- sapply(strsplit(file, "_"), `[`, 2)
		stage <- sapply(strsplit(file, "_"), `[`, 1)
		gr <- rtracklayer::import(
			paste(folder, file, sep="/"),
			which=locus+window
		)
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
				genome = "mm10",
				type = "hist",
				baseline=0,
				window=res,
				chromosome = as.character(seqnames(locus)),
				name = unique(gr$cond),
				col.baseline=color,
				col.histogram=color,
				fill.histogram=color,
				ylim=c(0, trunc(max_score, digit=4)),
				yTicksAt=c(0,trunc(max_score, digit=4))
				# alpha=0.8
			)
		}
	)
	track <- OverlayTrack(gTrack_list)
	return(track)
}

plot_tracks <- function(peak_gr, genome, link, plot_list, locus, window){

	locus_TSS <- resize(gene_TSS[gene_TSS$name %in% locus,], width=2, fix='start')

	locus_gr <- locus_TSS+window
	peak_gr_locus <- subsetByOverlaps(peak_gr, locus_gr)


	gene_links <- links[links$Gene %in% locus,]
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



	genome_track <- GenomeAxisTrack(
		col = "#333333", 
		fontsize=8,
		lwd=1,
		distFromAxis=0.5,
		labelPos = "alternating",
		sizes=0.3
	)

	gene_track <- GeneRegionTrack(
		genome,
		collapseTranscripts = "meta",
		name = "Genes",
		chromosome = as.character(seqnames(locus_gr)),
		start = start(locus_gr),
		end = end(locus_gr),
		fontface.group="bold.italic",
		col = NULL,
		fill = "#585858", 
		fontsize.group=18,
		sizes=0.3,
		rotation.title=0
		# just.group="above"
	)

	peak_track <- AnnotationTrack(
		peak_gr_locus,
		chromosome = as.character(seqnames(peak_gr_locus)),
		start = start(peak_gr_locus),
		end = end(peak_gr_locus),
		strand = as.character(strand(peak_gr_locus)),
		name = "Peaks",
		col = "#5f4780",
		fill = "#5f4780",
		sizes=0.2,
		rotation.title=0
	)


	visible_linked_peaks <- subsetByOverlaps(link_peak_gr, locus_gr)

	ht <- HighlightTrack(
		trackList = c(peak_track, gene_track, plot_list),
		start = start(visible_linked_peaks),
		end = end(visible_linked_peaks),
		col=visible_linked_peaks$colors,
		fill=visible_linked_peaks$colors,
		chromosome=as.character(seqnames(visible_linked_peaks)),
		alpha=0.1
	)

	# plotTracks(list(itrack, gtrack, ht), from = lim[1], to = lim[2])


	plot <- plotTracks(
		c(link_track, ht, genome_track),
		chromosome = as.character(seqnames(locus_gr)),
		from = start(locus_gr),
		to = end(locus_gr),
		fontsize=10,
		transcriptAnnotation="symbol",
		background.title = "transparent",
		col.border.title="transparent",
		col.title = "#333333",
		col.axis = "#333333",
		# rotation.title=0,
		alpha.title=1,
		alpha.axis=1,
		sizes=c(0.5, 0.2, 0.3, rep(0.5, length(plot_list)), 0.4),
		cex.title=0.8,
		title.width=0.7,
		stackHeight=2
		# rotation.title=0
	)

	return(plot)
}

###########################################
###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################
# samplesheet <- snakemake@input[['samplesheet']]
# bw_folder <- snakemake@params[['bw_folder']]
# window <- snakemake@params[['window']]
# genome <- snakemake@input[['genome']]
# peaks <- snakemake@input[['peaks']]
# links <- snakemake@input[['links']]
# gene_bed <- snakemake@input[['gene_bed']]


bw_folder <- "results/processed_data/mm10/bigWig"



bw_files <- list.files(path=bw_folder, pattern = "REP..bw")
genome <- "workflow/data/mm10/iGenome_mm10_ucsc_genes.gtf.gz"

peaks <- "results/processed_data/mm10/ATAC_norm_counts.csv"
links <- "results/tables/mm10/all_sig_gene2peak_linkage.csv"
gene_bed <- "workflow/data/mm10/gene_standard.bed"

links <- read.table(links, header=TRUE, row.names=1)

peak_gr <- GRanges(rownames(read.csv(peaks, header=TRUE, row.names=1)))
gene_TSS <- GRanges(rtracklayer::import(gene_bed))



window <- 85000
locus_gene <- "Amh"



conditions <- paste(
	sapply(strsplit(bw_files, "_"), `[`, 1),
	sapply(strsplit(bw_files, "_"), `[`, 2),
	sep="_"
)

conditions <- unique(conditions)

conditions <- c(conditions[c(TRUE, FALSE)], conditions[c(FALSE, TRUE)])

names(conditions_color) <- conditions


atac__coverage <- import_bw_files(
	bw_folder, 
	bw_files,
	locus_gene, 
	window
)

max_score <- do.call("c", do.call("c", do.call("c", atac__coverage)))
max_score <- max(unlist(lapply(max_score, function(x) x$score)))


link_track <- generate_link_track(locus_gene, links)
XX_track_list <- lapply(atac__coverage[["XX"]], function(stage) generate_genomic_tracks(stage, locus_gene, window, max_score, conditions_color))
XY_track_list <- lapply(atac__coverage[["XY"]], function(stage) generate_genomic_tracks(stage, locus_gene, window, max_score, conditions_color))

pdf(file="test.pdf", width=8, height=8)
	plot <- plot_tracks(peak_gr, genome, links, c(XX_track_list, XY_track_list), locus_gene, window)
dev.off()

