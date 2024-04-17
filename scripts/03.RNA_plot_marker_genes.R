source("scripts/00.functions.R")
source("scripts/00.color_palettes.R")

TPM <- read.csv(file=snakemake@input[['tpm']], row.names=1)

##########################################
#                                        #
#              Marker Genes              #
#                                        #
##########################################

# Load marker gene list
markerGenes <- read.csv("data/gonad_marker_genes.csv")

# Load whole gonad RNA-seq data
whole_gonad <- read.csv(
	"data/Zhao_tpm_matrix_for_analysis.csv", 
)

# Rename RNA-seq data to label if whole gonad or supporting cells
rownames(whole_gonad) <- make.names(whole_gonad$X, unique = TRUE)
whole_gonad <- whole_gonad[,-1]
colnames(whole_gonad) <- paste0("whole.gonad_", colnames(whole_gonad))

colnames(TPM) <- paste0("supporting_", colnames(TPM))

# Merge supporting and whole gonad RNA-seq data
TPM <- merge(TPM,whole_gonad,by="row.names",all.x=TRUE)
rownames(TPM) <- TPM[,1]
TPM <- TPM[,-1]

# Plot expression of the marker genes
plot_list <- plot_expression_scatter(TPM, markerGenes)

figure <- plot_grid(
	plotlist=plot_list,
	labels = "AUTO",
	ncol=1
)

save_plot(
	snakemake@output[['pdf']], 
	figure, 
	base_width=38,
	base_height=30,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)

save_plot(
	snakemake@output[['png']], 
	figure, 
	base_width=38,
	base_height=30,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)
