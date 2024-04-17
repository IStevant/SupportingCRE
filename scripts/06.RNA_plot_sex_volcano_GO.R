# source("scripts/00.functions.R")
source("scripts/00.color_palettes.R")

# install.packages('doParallel', repos="https://mirror.ibcp.fr/pub/CRAN/")
# install.packages('foreach', repos="https://mirror.ibcp.fr/pub/CRAN/")
BiocManager::install('clusterProfiler')

renv::snapshot()

library('doParallel')
library('foreach')


doParallel::registerDoParallel(cores=4)

load(snakemake@input[['all_DEGs']])
samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)
stages <- unique(samplesheet$stages)
names(conditions_color) <- sort(unique(samplesheet$conditions))


adj.pval <- snakemake@params[['adjpval']]
log2FC <- snakemake@params[['log2FC']]

###########################################
#                                         #
#         Volcano + GO sex DEGs           #
#                                         #
###########################################

sex_DEG_plots <- foreach(stg=stages) %dopar% {
	plot_DEG(SexDEGs, stg, adj.pval, log2FC, conditions_color)
}

figure <- plot_grid(
	plotlist=sex_DEG_plots,
	labels = "AUTO",
	ncol=1
)

save_plot(
	snakemake@output[['pdf']], 
	figure, 
	base_width=28,
	base_height=48,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)

save_plot(
	snakemake@output[['png']], 
	figure, 
	base_width=28,
	base_height=48,
	units = c("cm"), 
	dpi=300, 
	bg = "white"
)
