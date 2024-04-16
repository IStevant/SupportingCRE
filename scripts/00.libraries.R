# Install and load packages
packages <- c(
    "bioc::DESeq2",
    "bioc::clusterProfiler",
    "viridis",
    "ggrepel",
    "RColorBrewer",
    "bioc::ComplexHeatmap",
    "cowplot",
    "UpSetR",
    "eulerr",
    "grid",
    "pheatmap",
    "seriation",
    "MetBrewer",
    "bioc::org.Mm.eg.db",
    "foreach",
    "doParallel",
    "ggplot2",
    "ComplexUpset"
)

renv::install(packages)

# Update Renv
renv::snapshot()