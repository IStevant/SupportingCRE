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
#               Load data                 #
#                                         #
###########################################

load(file = snakemake@input[["peak_list"]])
samplesheet <- read.csv(file = snakemake@input[["samplesheet"]], row.names = 1)
names(conditions_color) <- sort(unique(samplesheet$conditions))

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Plot the number of peaks per condition
#' @param peaks list of peaks from each sex and stages.
#' @param conditions_color vector of hexadecimal colors.
#' @return Return a ggplot object.
plot_distribution <- function(peaks, conditions_color) {
  peaks_XX <- unlist(peaks[grep("XX", names(peaks))])
  peaks_XY <- unlist(peaks[grep("XY", names(peaks))])

  stages_XX <- sapply(strsplit(names(peaks_XX), "_"), `[`, 1)
  stages_XY <- sapply(strsplit(names(peaks_XY), "_"), `[`, 1)

  nb_peaks_XX <- unlist(lapply(peaks_XX, length))
  nb_peaks_XY <- unlist(lapply(peaks_XY, length))

  stage <- c(stages_XX, stages_XY)
  sex <- c(rep("Pre-gran.", length(stages_XX)), rep("Sertoli", length(stages_XY)))

  data <- data.frame(
    Stage = stage,
    Sex = sex,
    Counts = c(nb_peaks_XX, nb_peaks_XY),
    condition = paste(sex, stage)
  )

  plot <- ggplot(data, aes(x = Sex, y = Counts)) +
    geom_bar(stat = "identity", aes(fill = condition)) +
    geom_text(aes(label = scales::comma(Counts)), size = 5, position = position_stack(vjust = 0.5), color = "white") +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = conditions_color, 0.8) +
    facet_wrap(~Stage, nrow = 1) +
    theme_light() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.position = "none"
    )
  return(plot)
}

###########################################
#                                         #
#         Plot peak distribution          #
#                                         #
###########################################

# There is an unpredicted behaviour here. If I run the function once, I get 1 peak only for one sample, which is not correct
# If I run the function a second time the peak number is correct...
plot_list <- plot_distribution(peak_list, conditions_color)
plot_list <- plot_distribution(peak_list, conditions_color)

figures <- plot_grid(
  plotlist = list(plot_list)
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  figures,
  base_width = 20,
  base_height = 10,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  figures,
  base_width = 20,
  base_height = 10,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
