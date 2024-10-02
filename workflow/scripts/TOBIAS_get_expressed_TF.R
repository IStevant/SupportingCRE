source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
	library("TFBSTools")
	library("motifStack")
	library("JASPAR2024")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

TPM <- read.csv(file = snakemake@input[["TPM"]], header = TRUE, row.names = 1)
TFs <- as.vector(read.csv(snakemake@input[["TF_genes"]], header = FALSE)[, 1])
minTPM <- snakemake@params[["minTPM"]]

# Load Jaspar 2024 database
JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' When the maximum value (read count) of a peak between samples is under a certain threshold, we considere it is not relevant and the values are set to 0.
#' @param data Read count matrix.
#' @param minExp Minimum value. Default is 5.
#' @return Return a dataframe.
run_filter_low_counts <- function(data, minExp = 5) {
  col_names <- colnames(data)
  data <- t(apply(data, 1, filter_low_counts, col_names = col_names, minExp = minExp))
  return(as.data.frame(data))
}
filter_low_counts <- function(row, col_names, minExp) {
  if (max(row) < minExp) {
    return(setNames(rep(0, length(row)), col_names))
  } else {
    return(setNames(row, col_names))
  }
}

###########################################
#                                         #
#             Get TF matrices             #
#                                         #
###########################################

genes <- rownames(run_filter_low_counts(TPM, minTPM))
# Select only the TFs
exp_TFs <- genes[which(genes %in% TFs)]

# Get all vertebrate TF matrices
pfms <- TFBSTools::getMatrixSet(
  JASPAR,
  opts = list(
    matrixtype = "PFM",
    tax_group = "vertebrates"
  )
)


tf_matrix_names <- sapply(pfms, function(x) name(x))

selected_pfms <- pfms[tolower(tf_matrix_names) %in% tolower(exp_TFs)]

pcms <- universalmotif::convert_motifs(selected_pfms, class = "motifStack-pcm")


motif_names <- sapply(pcms, function(x) x$name)
duplicated_names <- make.unique(motif_names, sep = "ZZZ")
duplicated_names <- gsub("-", "XXX", duplicated_names)
duplicated_names <- toupper(duplicated_names)
duplicated_names <- make.unique(duplicated_names, sep = "ZZZ")


pcms_uniq <- lapply(seq_along(pcms), function(x) {
	pcms[[x]]$name <- duplicated_names[x]
	return(pcms[[x]])
})

names(pcms_uniq) <- duplicated_names

pcms_uniq <- pcms_uniq[order(names(pcms_uniq))]

hc <- motifStack::clusterMotifs(pcms_uniq)
phylog <- ade4::hclust2phylog(hc)
motifSig <- motifStack::motifSignature(pcms_uniq, phylog, cutoffPval = 0.001, min.freq = 1)
## get the signatures from object of motifSignature
sig <- motifStack::signatures(motifSig)

sig <- lapply(sig, function(x) {
	x$name <- gsub("ZZZ", ".", x$name)
	x$name <- gsub("XXX", "-", x$name)
	return(x)
})

sig <- universalmotif::convert_motifs(sig, class = "motifStack-pcm")


output_file <- snakemake@output[["matrices"]]

fileConn <- file(output_file, "w")

for (i in seq_along(sig)) {
  pfm <- sig[[i]]
  
  tf_name <- pfm$name

  writeLines(paste0(">", tf_name), fileConn)
  
  pfm_matrix <- pfm$mat
  
  bases <- rownames(pfm_matrix)
  
  for (j in 1:nrow(pfm_matrix)) {
    pfm_values <- paste(pfm_matrix[j, ], collapse = "   ")  # Espaces entre les valeurs
    writeLines(paste(bases[j], "[", pfm_values, "]"), fileConn)
  }
}

close(fileConn)

