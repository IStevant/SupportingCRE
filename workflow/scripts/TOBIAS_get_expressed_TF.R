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

simplify_genes <- function(genes) {
  gene_list <- sort(unlist(strsplit(genes, ";")))
    simplified_genes <- c()
  
  i <- 1
  while (i <= length(gene_list)) {
    gene <- gene_list[i]
    prefix <- sub("[0-9]+.*$", "", gene)
    
    matching_genes <- grep(paste0("^", prefix), gene_list[i:length(gene_list)], value = TRUE)
    
    if (length(matching_genes) > 1) {
      suffixes <- sub("^.*?([0-9]+.*$)", "\\1", matching_genes)
      simplified_gene <- paste0(prefix, paste(suffixes, collapse = "/"))
      simplified_genes <- c(simplified_genes, simplified_gene)
      i <- i + length(matching_genes) - 1
    } else {
      simplified_genes <- c(simplified_genes, gene)
    }
    
    i <- i + 1
  }
  return(paste(simplified_genes, collapse = ";"))
}

extract_prefix <- function(gene) {
  prefix <- sub("([a-zA-Z]+).*", "\\1", gene)
  if (prefix == "NR") {
    return(gene)
  } else if (grepl("^HOX", prefix)) {
      prefix <- stringr::str_sub(prefix, end = -2)
      return(prefix)
  } else {
        return(prefix)
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

merged_TF_motifs <- lapply(sig, function(x) {
  x$name <- gsub("ZZZ.", "", x$name)
  x$name <- gsub("XXX", "-", x$name)
  x$name <- extract_prefix(x$name)
  if (common_prefix != "HOX") {
    suffixes <- gsub(paste0("^", common_prefix), "", gene_group)
    numeric_suffixes <- sort(as.numeric(suffixes[grepl("^\\d+$", suffixes)]), na.last = TRUE)
    non_numeric_suffixes <- sort(suffixes[!grepl("^\\d+$", suffixes)])
    all_suffixes <- c(non_numeric_suffixes, numeric_suffixes)
    concatenated_suffixes <- paste(all_suffixes, collapse = "/")
    paste0(common_prefix, concatenated_suffixes)
  } else {
    paste0(common_prefix, "s")
  }
  return(x)
})

merged_TF_motifs <- universalmotif::convert_motifs(merged_TF_motifs, class = "motifStack-pcm")

output_file <- snakemake@output[["matrices_txt"]]

fileConn <- file(output_file, "w")

for (i in seq_along(merged_TF_motifs)) {
  pfm <- merged_TF_motifs[[i]]
    tf_name <- pfm$name
  writeLines(paste0(">", tf_name), fileConn)
  pfm_matrix <- pfm$mat
  bases <- rownames(pfm_matrix)
  for (j in 1:nrow(pfm_matrix)) {
    pfm_values <- paste(pfm_matrix[j, ], collapse = "   ")
    writeLines(paste(bases[j], "[", pfm_values, "]"), fileConn)
  }
}

close(fileConn)

