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

# filtered_SexDEGs
load(snakemake@input[["sex_DEGs"]])

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


draw_venn_sex_dyn <- function(sexID, sex, dyn, title, output_folder) {
  data <- list(Sex.biased=sex, Dynamic=dyn)

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
    file=paste0(output_folder, "RNA_common_sex_dynamic_DEGs.tsv"), 
    row.name=FALSE, 
    quote=FALSE,
    sep="\t"
  )


  if (sexID == "XX") {
    colours <- c("#FCBF49", "#D62828")
  } else {
    colours <- c("#94d574", "#1e8bd1")
  }

  venn <- eulerr::euler(data)
  title <- title
  venn_plot <- plot(
    venn,
    labels = c(
      paste0("Sexually dymorphic\n(", prettyNum(length(sex), big.mark = ","), ")"),
      paste0("Dynamic\n(", prettyNum(length(dyn), big.mark = ","), ")")
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
XX_spe_genes <- unique(unlist(lapply(filtered_SexDEGs, function(x) rownames(x[x$Diff.Exp. == "Up in XX", ]))))

XX_dyn_genes <- XX_filtered_StageDEGs
XX_venn <- draw_venn_sex_dyn(sexID = "XX", XX_spe_genes, XX_dyn_genes, "XX genes", output_folder)

XY_spe_genes <- unique(unlist(lapply(filtered_SexDEGs, function(x) rownames(x[x$Diff.Exp. == "Up in XY", ]))))
XY_dyn_genes <- XY_filtered_StageDEGs
XY_venn <- draw_venn_sex_dyn(sexID = "XY", XY_spe_genes, XY_dyn_genes, "XY genes", output_folder)


figure <- plot_grid(
  XX_venn, XY_venn,
  labels = c("XX", "XY"),
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
  base_width = 12,
  base_height = 15,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 12,
  base_height = 15,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
