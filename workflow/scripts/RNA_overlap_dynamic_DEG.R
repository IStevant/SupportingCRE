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

# filtered_StageDEGs
load(snakemake@input[["XY_stage_DEGs"]])
XY_filtered_StageDEGs <- filtered_StageDEGs

# filtered_StageDEGs
load(snakemake@input[["XX_stage_DEGs"]])
XX_filtered_StageDEGs <- filtered_StageDEGs

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

draw_venn_sex_dyn <- function(set1, set2, title, output_folder) {
  data <- list(XX_dynamic=set1, XY_dynamic=set2)

  intersections <- get_intersections(data)

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
    file=paste0(output_folder, "RNA_common_dynamic_DEGs.tsv"), 
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
      paste0("Dynamic in XX\n(", prettyNum(length(set2), big.mark = ","), ")"),
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

XX_dyn_genes <- XX_filtered_StageDEGs
XY_dyn_genes <- XY_filtered_StageDEGs

venn <- draw_venn_sex_dyn(XX_dyn_genes, XY_dyn_genes, "Dynamic genes", output_folder)


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
