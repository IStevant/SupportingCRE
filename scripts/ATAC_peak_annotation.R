
source("scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

# suppressPackageStartupMessages({
# 	library("cowplot")
# 	library("grid")
# 	library("viridis")
# 	library("ggplot2")
# })

install.packages("rtracklayer", repos="https://mirror.ibcp.fr/pub/CRAN/")

renv::snapshot()
###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################




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
