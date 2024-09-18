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
#               Load data                 #
#                                         #
###########################################

# filtered_SexDEGs
load(snakemake@input[["sig_DEGs"]])
output_folder <- snakemake@params[["output_folder"]]

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Get genes contained in each intersections.
#' @param sets_list Sex DESeq2 analysis result table
get_intersections <- function(sets_list) {
  intersections <- list()
  
  all_sets <- names(sets_list)
  comb <- expand.grid(lapply(all_sets, function(x) c(TRUE, FALSE)))
  comb <- comb[-nrow(comb), ]
  
  for (i in 1:nrow(comb)) {
    included_sets <- all_sets[comb[i, ] == TRUE]
    excluded_sets <- all_sets[comb[i, ] == FALSE]
    
    intersection_genes <- Reduce(intersect, sets_list[included_sets])
    
    if (length(excluded_sets) > 0) {
      excluded_genes <- unlist(sets_list[excluded_sets])
      intersection_genes <- setdiff(intersection_genes, excluded_genes)
    }
    
    intersections[[paste(included_sets, collapse = "+")]] <- intersection_genes
  }
  
  return(intersections)
}


#' Draw upset plot.
#' @param DEGs Sex DESeq2 analysis result table
#' @param sex Genetic sex of the cells, either "XX" or "XY".
#' @return Grig object.
upset_plots_sex <- function(DEGs, sex, output_folder) {
  current_sex <- paste("Up in", sex)
  filtered_SexDEGs <- lapply(DEGs, function(DEG) rownames(DEG[DEG$Diff.Exp. == current_sex, , drop = FALSE]))
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


  venn <- eulerr::euler(filtered_SexDEGs)


  intersections <- get_intersections(filtered_SexDEGs)

  intersection_df <- do.call(
    rbind, 
    lapply(
      names(intersections), 
      function(set) {
        data.frame(Intersection = set, Gene = intersections[[set]])
      }
    )
  )

  write.table(
    intersection_df, 
    file=paste0(output_folder, "RNA_upset_sex_DEGs_", sex, ".tsv"), 
    row.name=FALSE, 
    quote=FALSE,
    sep="\t"
  )


  data <- UpSetR::fromExpression(venn$original.values)
  colnames(data) <- names(DEGs)
  rownames(data) 
  data <- data == 1
  data <- as.data.frame(data)


  plot <- ComplexUpset::upset(
    data = data,
    intersect = names(DEGs),
    name = "Sex-specific DE genes",
    sort_sets = "ascending",
    stripes = alpha("white", 0),
    set_sizes = FALSE,
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

###########################################
#                                         #
#               Plot Upset                #
#                                         #
###########################################

XX_upset <- upset_plots_sex(filtered_SexDEGs, "XX", output_folder)
XY_upset <- upset_plots_sex(filtered_SexDEGs, "XY", output_folder)


figure <- plot_grid(
  plotlist = list(XX_upset, XY_upset),
  labels = "AUTO",
  ncol = 1
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 28,
  base_height = 20,
  units = c("cm"),
  dpi = 300
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
