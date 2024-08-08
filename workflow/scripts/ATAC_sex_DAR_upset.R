source(".Rprofile")
###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("cowplot")
  library("ggplot2")
})


###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

upset_plots_sex <- function(DARs, sex) {
  current_sex <- paste("More in", sex)
  filtered_SexDARs <- lapply(DARs, function(DAR) rownames(DAR[DAR$Diff.Acc. == current_sex, , drop = FALSE]))
  venn <- eulerr::euler(filtered_SexDARs)
  if (sex == "XX") {
    color <- "#FFB100"
    text_col <- c("black", "black")
  } else {
    color <- "#339989"
    text_col <- c("black", "white")
  }

  theme_upset <- function() {
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )
  }

  venn <- eulerr::euler(filtered_SexDARs)
  data <- UpSetR::fromExpression(venn$original.values)
  colnames(data) <- names(DARs)
  data <- data == 1
  data <- as.data.frame(data)

  plot <- ComplexUpset::upset(
    data = data,
    intersect = names(DARs),
    name = "Sex differentially accessible regions",
    sort_sets = "ascending",
    stripes = alpha("white", 0),
    set_sizes = FALSE,
    # wrap=TRUE
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size(
        text_colors = text_col,
        mapping = aes(fill = "bars_color")
      ) +
        scale_fill_manual(values = c("bars_color" = color)) +
        scale_y_continuous(labels = scales::comma) +
        theme_upset()
    ),
    matrix = ComplexUpset::intersection_matrix(
      segment = geom_segment(color = color),
      outline_color = list(active = color, inactive = "grey70"),
    ) + scale_color_manual(
      values = c("TRUE" = color, "FALSE" = "grey"),
      na.value = "transparent",
      guide = "none"
    ),
    themes = ComplexUpset::upset_default_themes(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
  )

  return(plot)
}
#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# filtered_SexDARs
load(snakemake@input[["sig_DARs"]])

###########################################
#                                         #
#               Plot Upset                #
#                                         #
###########################################

XX_upset <- upset_plots_sex(filtered_SexDARs, "XX")
XY_upset <- upset_plots_sex(filtered_SexDARs, "XY")


figure <- plot_grid(
  plotlist = list(XX_upset, XY_upset),
  labels = "AUTO",
  ncol = 1
)
figure

save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 28,
  base_height = 20,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)

save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 28,
  base_height = 20,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
