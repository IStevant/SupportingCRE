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

# filtered_SexDARs
load(snakemake@input[["sex_DARs"]])

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

#' Get regions contained in each intersections.
#' @param sets_list Sex DESeq2 analysis result table
get_intersections <- function(sets_list) {
  intersections <- list()
  
  all_sets <- names(sets_list)
  comb <- expand.grid(lapply(all_sets, function(x) c(TRUE, FALSE)))
  comb <- comb[-nrow(comb), ]
  
  for (i in 1:nrow(comb)) {
    included_sets <- all_sets[comb[i, ] == TRUE]
    excluded_sets <- all_sets[comb[i, ] == FALSE]
    
    intersection_regions <- Reduce(intersect, sets_list[included_sets])
    
    if (length(excluded_sets) > 0) {
      excluded_regions <- unlist(sets_list[excluded_sets])
      intersection_regions <- setdiff(intersection_regions, excluded_regions)
    }
    
    intersections[[paste(included_sets, collapse = "+")]] <- intersection_regions
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
    file=paste0(output_folder, "ATAC_common_sex_dynamic_DARs.tsv"), 
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
XX_spe_regions <- unique(unlist(lapply(filtered_SexDARs, function(x) rownames(x[x$Diff.Acc. == "More in XX", ]))))

XX_dyn_regions <- XX_filtered_StageDARs
XX_venn <- draw_venn_sex_dyn(sexID = "XX", XX_spe_regions, XX_dyn_regions, "XX regions", output_folder)

XY_spe_regions <- unique(unlist(lapply(filtered_SexDARs, function(x) rownames(x[x$Diff.Acc. == "More in XY", ]))))
XY_dyn_regions <- XY_filtered_StageDARs
XY_venn <- draw_venn_sex_dyn(sexID = "XY", XY_spe_regions, XY_dyn_regions, "XY regions", output_folder)


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
