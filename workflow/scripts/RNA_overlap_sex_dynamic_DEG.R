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

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

draw_venn_sex_dyn <- function(sexID, sex, dyn, title) {
  data <- c(
    sex = length(setdiff(sex, dyn)),
    dyn = length(setdiff(dyn, sex)),
    "sex&dyn" = length(intersect(sex, dyn))
  )
  write.csv(setdiff(sex, dyn), file = "XX_dyn_only_DE_genes.csv")
  write.csv(setdiff(dyn, sex), file = "XY_dyn_only_DE_genes.csv")
  write.csv(intersect(dyn, sex), file = "both_dyb_DE_genes.csv")


  colours <- c("#D62828", "#1e8bd1")

  venn <- eulerr::euler(data)
  title <- title
  venn_plot <- plot(
    venn,
    labels = c(
      paste0("Dynamic in XX\n(", prettyNum(length(dyn), big.mark = ","), ")"),
      paste0("Dynamic in XY\n(", prettyNum(length(dyn), big.mark = ","), ")")
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

venn <- draw_venn_sex_dyn(sexID = "XX", XX_dyn_genes, XY_dyn_genes, "Dynamic genes")


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
