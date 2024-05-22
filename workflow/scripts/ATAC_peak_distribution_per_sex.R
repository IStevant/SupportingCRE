source(".Rprofile")

source("workflow/scripts/00.color_palettes.R")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
	library("cowplot")
	library("grid")
	library("ggplot2")
})

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

plot_distribution <- function(peaks, conditions_color){

	# peaks <- peak_list

	peaks_XX <- peaks[["XX"]]
	peaks_XY <- peaks[["XY"]]

	stages_XX <- names(peaks_XX)
	stages_XY <- names(peaks_XY)

	nb_peaks_XX <-	unlist(lapply(stages_XX, function(stg) length(peaks_XX[[stg]])))
	nb_peaks_XY <- unlist(lapply(stages_XY, function(stg) length(peaks_XY[[stg]])))

	stage <- c(stages_XX, stages_XY)
	sex <- c(rep("XX", length(stages_XX)), rep("XY", length(stages_XY)))

	data <- data.frame(
		Stage = stage,
		Sex = sex,
		Counts = c(nb_peaks_XX, nb_peaks_XY), 
		condition = paste(sex, stage)
	)

	plot <- ggplot(data, aes(x=Sex, y=Counts)) +
		geom_bar(stat="identity", aes(fill=condition)) +
		geom_text(aes(label=scales::comma(Counts)), size = 5, position = position_stack(vjust = 0.5), color="white") +
		scale_y_continuous(labels = scales::comma) + 
		# scale_fill_paletteer_d("MetBrewer::Hiroshige") +
		scale_fill_manual(values=conditions_color, 0.8) +
		# scale_color_manual(values=c(rep("white", 3), rep("black", 8)), guide="none") +
		# ggtitle("Genomic features overlapped by expressed TEs") +
		facet_wrap(~Stage, nrow=1) +
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
# load(file="workflow/data/mm10/ATAC_all_consensus_peaks_2rep_list.Robj")
# samplesheet <- read.csv(file="results/processed_data/mm10/ATAC_samplesheet.csv", row.names=1)
# filtered_ATAC <- read.csv(file="results/processed_data/mm10/ATAC_norm_counts.csv", row.names=1)

samplesheet <- read.csv(file=snakemake@input[['samplesheet']], row.names=1)

names(conditions_color) <- sort(unique(samplesheet$conditions))

###########################################
#                                         #
#           Plot peak peakstation          #
#                                         #
###########################################

plot_list <- plot_distribution(peak_list, conditions_color)

figures <- plot_grid(
	plotlist=list(plot_list)
)

save_plot(
	snakemake@output[['pdf']],
	# "test.pdf",
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
