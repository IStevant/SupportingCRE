source(".Rprofile")
###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("cowplot")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# filtered_StageDARs
load(snakemake@input[["XY_stage_DARs"]])
XY_filtered_StageDARs <- filtered_StageDARs

# filtered_StageDARs
load(snakemake@input[["XX_stage_DARs"]])
XX_filtered_StageDARs <- filtered_StageDARs

output_folder <- snakemake@params[["output_folder"]]

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Get elements contained in each intersections.
#' @param sets_list Sex DESeq2 analysis result table
get_intersections <- function(sets_list) {
  intersections <- list()
  
  all_sets <- names(sets_list)
  comb <- expand.grid(lapply(all_sets, function(x) c(TRUE, FALSE)))
  comb <- comb[-nrow(comb), ]
  
  for (i in 1:nrow(comb)) {
    included_sets <- all_sets[comb[i, ] == TRUE]
    excluded_sets <- all_sets[comb[i, ] == FALSE]
    
    intersection_elements <- Reduce(intersect, sets_list[included_sets])
    
    if (length(excluded_sets) > 0) {
      excluded_elements <- unlist(sets_list[excluded_sets])
      intersection_elements <- setdiff(intersection_elements, excluded_elements)
    }
    
    intersections[[paste(included_sets, collapse = "+")]] <- intersection_elements
  }
  
  return(intersections)
}

draw_venn_sex_dyn <- function(set1, set2, title, output_folder) {
  data <- list(XX_dynamic=set1, XY_dynamic=set2)

  intersections <- get_intersections(data)

  intersection_df <- do.call(
    rbind, 
    lapply(
      names(intersections), 
      function(set) {
        data.frame(Intersection = set, OCRs = intersections[[set]])
      }
    )
  )

  write.table(
    intersection_df, 
    file=paste0(output_folder, "ATAC_common_dynamic_DARs.tsv"), 
    row.name=FALSE, 
    quote=FALSE,
    sep="\t"
  )

  colours <- c("#D62828", "#1e8bd1")
  venn <- eulerr::euler(data)
  title <- title
  venn_plot <- plot(
    venn,
    labels = c(
      paste0("Dynamic in XX\n(", prettyNum(length(set1), big.mark = ","), ")"),
      paste0("Dynamic in XY\n(", prettyNum(length(set2), big.mark = ","), ")")
    ),
    quantities = list(fontsize = 10),
    edges = list(
      col = as.vector(colours),
      lex = 2
    ),
    fills = list(
      fill = as.vector(colours),
      alpha = 0.35
    )
  )
}

###########################################
#                                         #
#               Draw Venn                 #
#                                         #
###########################################

XX_dyn_OCRs <- XX_filtered_StageDARs
XY_dyn_OCRs <- XY_filtered_StageDARs

venn <- draw_venn_sex_dyn(XX_dyn_OCRs, XY_dyn_OCRs, "Dynamic OCRs", output_folder)


###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  venn,
  base_width = 12,
  base_height = 8,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  venn,
  base_width = 12,
  base_height = 8,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
