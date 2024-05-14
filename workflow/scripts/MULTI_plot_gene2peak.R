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
	# library("BSgenome.Mmusculus.UCSC.mm10")
	library("TxDb.Mmusculus.UCSC.mm10.knownGene")
	library("Mus.musculus")
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
	cond <- stringr::str_extract(names(gr_list), "^.{8}")
	color <- colors[cond]
	res=1000
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


plot_tracks <- function(plot_list, locus, window){

	locus_gr <- locus+window

	genome_track <- GenomeAxisTrack(
		col = "#333333", 
		fontsize=8,
		lwd=1,
		distFromAxis=0.5,
		labelPos = "alternating",
		sizes=0.3
	)

	gene_track <- GeneRegionTrack(
		TxDb.Mmusculus.UCSC.mm10.knownGene,
		collapseTranscripts = "meta",
		name = "Genes",
		chromosome = as.character(seqnames(locus_gr)),
		start = start(locus_gr),
		end = end(locus_gr),
		fontface.group="bold.italic",
		col = "#585858",
		fill = "#585858", 
		fontsize.group=18,
		sizes=0.3
		# just.group="above"
	)

	z <- ranges(gene_track)
	z$symbol <- mapIds(Mus.musculus, z$symbol, "SYMBOL", "TXNAME")
	z$symbol <- paste0(z$symbol, " ")
	ranges(gene_track) <- z

	genom_plot_list <- list(
		gene_track,
		genome_track
	)


	plot <- plotTracks(
		c(gene_track, plot_list, genome_track),
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
		sizes=c(0.3, rep(0.5, length(plot_list)), 0.3),
		cex.title=0.8,
		title.width=0.7
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

bw_folder <- "results/processed_data/mm10/bigWig"
window <- 85000
locus_test <- GRanges("chr2:105124529-105124530")
bw_files <- list.files(path=bw_folder, pattern = "REP..bw")

conditions <- paste(
	sapply(strsplit(bw_files, "_"), `[`, 1),
	sapply(strsplit(bw_files, "_"), `[`, 2),
	sep="_"
)

conditions <- unique(conditions)

conditions <- c(conditions[c(TRUE, FALSE)], conditions[c(FALSE, TRUE)])

names(conditions_color) <- conditions


test <- import_bw_files(
	bw_folder, 
	bw_files,
	locus_test, 
	window
)

max_score <- do.call("c", do.call("c", do.call("c", test)))
max_score <- max(unlist(lapply(max_score, function(x) x$score)))

XX_track_list <- lapply(test[["XX"]], function(stage) generate_genomic_tracks(stage, locus_test, window, max_score, conditions_color))
XY_track_list <- lapply(test[["XY"]], function(stage) generate_genomic_tracks(stage, locus_test, window, max_score, conditions_color))


pdf(file="test.pdf", width=8, height=6)
plot_tracks(c(XX_track_list, XY_track_list), locus_test, window)
dev.off()

