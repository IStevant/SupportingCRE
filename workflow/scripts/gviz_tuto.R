# Load libraries
library("Gviz")
library("dplyr")
library('doParallel')
library('foreach')


# Load libraries
library("GenomicFeatures")
library("Gviz")


my_genome <- "path/to/mouse/gencode.vM34.annotation.gtf.gz"

my_genome <- "/media/zazooo/9C33-6BBD1/gencode.vM34.annotation.gtf.gz"

# Load GTF file as a TxDb object
TxDb <- GenomicFeatures::makeTxDbFromGFF(my_genome)

geneTrack <- GeneRegionTrack(
	rtacklayer::import(my_genome),
	# TxDb, 
	chromosome = "chr2",
	start=133372559,
	end=133426311,
	name = "Genes"
)

plotTracks(geneTrack)



head(gene(grtrack))
head(transcript(grtrack))
head(exon(grtrack))
head(symbol(grtrack))
plotTracks(grtrack)







genes_coord <- GenomicFeatures::genes(TxDb)
genome_gtf <- rtracklayer::import("/media/zazooo/9C33-6BBD1/gencode.vM34.annotation.gtf.gz")
gene2symbol <- mcols(genome_gtf)[,c("gene_id","gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

# Select the genomic window you want to see
my_gene <- "Bmp2"


columns(TxDb)


select(TxDb, keys = keys, columns="TXNAME", keytype="GENEID")


locus <- resize(gene_TSS[gene_TSS$name %in% locus,], width=2, fix='start')


gene_track <- GeneRegionTrack(
	TxDb,
	collapseTranscripts = "meta",
	name = "Genes",
	chromosome = as.character(seqnames(locus_gr)),
	start = start(locus_gr),
	end = end(locus_gr),
	# Mouse genes are in italic
	fontface.group="italic",
	# Exon color
	fill = "#585858",
	# Apply the same color for exon borders
	col = NULL,
	# Apply the same color to the intron lines
	col.line = NULL,
	# Gene name color
	fontcolor.group= "#333333",
	# Gene name fontsize
	fontsize.group=18
)

plotTrack(gene_track)






ranges(gene_track)$symbol <- gene2symbol[ranges(gene_track)$gene, "gene_name"]

genome_track <- GenomeAxisTrack(
	col = "#333333", 
	fontsize=13,
	lwd=1,
	distFromAxis=0.5,
	labelPos = "alternating",
)












TxDb <- GenomicFeatures::makeTxDbFromGFF("/media/zazooo/9C33-6BBD1/gencode.vM34.annotation.gtf.gz")
