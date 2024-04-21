source(".Rprofile")

source("scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
	library("cowplot")
	library("grid")
# 	library("viridis")
	library("ggplot2")
})

# install.packages("rtracklayer", repos="https://mirror.ibcp.fr/pub/CRAN/")

# renv::snapshot()
###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

plot_anno_sex <- function(anno){
	anno_XX <- anno[["XX"]]
	anno_XY <- anno[["XY"]]
	stages <- names(anno_XX)
	anno_anno_XX <- lapply(1:length(anno_XX), function(x) {
		stat <- anno_XX[[x]]@annoStat
		stat$stage <- rep(stages[x], nrow(stat))
		return(stat)
		}
	)
	data_XX <- data.table::rbindlist(anno_anno_XX)
	data_XX$Sex <- rep("XX", nrow(data_XX))

	anno_anno_XY <- lapply(1:length(anno_XY), function(x) {
		stat <- anno_XY[[x]]@annoStat
		stat$stage <- rep(stages[x], nrow(stat))
		return(stat)
		}
	)
	data_XY <- data.table::rbindlist(anno_anno_XY)
	data_XY$Sex <- rep("XY", nrow(data_XY))

	data <- rbind(data_XX, data_XY)

	labels <- round(data$Frequency, digit=2)
	labels[labels<3.5] <- 0
	labels <- paste0(labels, "%")
	labels[labels=="0%"] <- " "

	plot <- ggplot(data, aes(fill=Feature, x=Sex, y=Frequency)) + 
		geom_bar(stat="identity") +
		geom_text(aes(label=labels, color=Feature), size = 5, position = position_stack(vjust = 0.5)) +
		# geom_text(aes(label = scales::comma(after_stat(y)), group = stage), size = 4, stat = 'summary', fun = sum, vjust = -0.3) +
		scale_y_continuous(labels = scales::comma) + 
		# scale_fill_paletteer_d("MetBrewer::Hiroshige") +
		scale_fill_manual(values=alpha(rev(MetBrewer::met.brewer("Hiroshige", n=length(unique(data$Feature)))), 0.8)) +
		scale_color_manual(values=c(rep("white", 3), rep("black", 8)), guide="none") +
		# ggtitle("Genomic features overlapped by expressed TEs") +
		facet_wrap(~stage, nrow=1) +
		theme_light() +
			theme(
			plot.title = element_text(size=14, hjust = 0.5, face="bold"),
			axis.text=element_text(size=14),
			axis.title=element_text(size=14),
			legend.title=element_blank(),
			legend.text=element_text(size=14, margin = margin(r = 10, unit = "pt")),
			# legend.position="bottom",
			# aspect.ratio=1.75,
			strip.text.x = element_text(size = 14, face="bold"),
			legend.box.spacing = unit(0, "mm"),
		)
	return(plot)
}



#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# peak_list
load(file=snakemake@input[['peak_list']])

# Promoter region, i.e. distance to TSS
promoter <- snakemake@params[['promoter']]

# Load mouse genome
mm10Genes <- rtracklayer::import("data/iGenome_mm10_ucsc_genes.gtf.gz")
# the column symbol was present in other version of the script so I add it again but ultimately change the $symbol to $gene_name
mm10Genes$symbol <- mm10Genes$gene_name

txdb <- GenomicFeatures::makeTxDbFromGRanges(
	mm10Genes,
	drop.stop.codons=FALSE
)

###########################################
#                                         #
#           Get peak annotation           #
#                                         #
###########################################

XX_peaks <- peak_list$XX
XY_peaks <- peak_list$XY

XX_anno <- lapply(
	XX_peaks, 
	function(stage) 
		ChIPseeker::annotatePeak(stage, tssRegion=c(-promoter, promoter), TxDb=txdb)
)

names(XX_anno) <- names(XX_peaks)


XY_anno <- lapply(
	XY_peaks, 
	function(stage) 
		ChIPseeker::annotatePeak(stage, tssRegion=c(-promoter, promoter), TxDb=txdb)
)

names(XY_anno) <- names(XY_peaks)

peak_anno_list <- list(
	XX=XX_anno,
	XY=XY_anno
)

save(peak_anno_list, file=snakemake@output[['anno_list']])

###########################################
#                                         #
#           Plot peak annotation          #
#                                         #
###########################################

plot_list <- plot_anno_sex(peak_anno_list)

figures <- plot_grid(
	plotlist=list(plot_list)
)


save_plot(
	snakemake@output[['pdf']], 
	figures, 
	base_width=30,
	base_height=20,
	units = c("cm"), 
	dpi=300
)

save_plot(
	snakemake@output[['png']], 
	figures, 
	base_width=30,
	base_height=20,
	units = c("cm"), 
	dpi=300
)