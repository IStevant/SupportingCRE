source(".Rprofile")


###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Generate the read count matrix
#' @param csv_file Path to the read count matrix.s
#' @return Return a dataframe.
get_peak_matrix <- function(csv_file, peaks) {
  # load file
  raw_counts <- read.csv(file = csv_file, header = TRUE, sep = "\t")
  # Define peak names
  peak_names <- paste0(raw_counts$Chr, ":", raw_counts$Start, "-", raw_counts$End)
  # Remove the extra columns of the file
  raw_counts <- raw_counts[, -c(1:6)]
  # Apply peak names
  rownames(raw_counts) <- peak_names
  # Transform values as integers for DESeq2 that does not support floats
  raw_counts <- round(raw_counts, digits = 0)
  # Rename samples
  sample_names <- colnames(raw_counts)
  sample_names <- sapply(strsplit(sample_names, ".mLb"), `[`, 1)
  stage <- sapply(strsplit(sample_names, "_"), `[`, 1)
  sex <- sapply(strsplit(sample_names, "_"), `[`, 2)
  transgene <- sapply(strsplit(sample_names, "_"), `[`, 3)
  transgene <- gsub("enh8", "enh8.mcherry", transgene)
  transgene <- gsub("sox9", "SOX9.IRES.GFP", transgene)
  rep <- sapply(strsplit(sample_names, "_"), `[`, 4)
  rep <- gsub("REP", "rep", rep)
  sample_names <- paste(stage, sex, transgene, rep, sep = "_")
  colnames(raw_counts) <- sample_names
  raw_counts <- raw_counts[, order(names(raw_counts))]
  raw_counts <- raw_counts[grep("^chr[0-9|X|Y]{1,2}[:]", rownames(raw_counts)), ]
  return(raw_counts)
}

#' Normalize the read count using the size factor normalization from DESeq2
#' @param raw_counts Read count matrix.
#' @param samplesheet Samplesheet for DESeq2.
#' @return Return a dataframe.
get_normalized_counts <- function(raw_counts, samplesheet) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  # norm_counts <- DESeq2::counts(dds, normalized=TRUE)
  norm_counts <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))

  return(norm_counts)
}

get_size_factors <- function(raw_counts, samplesheet) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  size_factors <- DESeq2::sizeFactors(dds)
  # norm_counts <- DESeq2::counts(dds, normalized=TRUE)
  # norm_counts <- SummarizedExperiment::assay(DESeq2::vst(dds, blind=FALSE))
  size_factors <- data.frame(
    sample = names(size_factors),
    SizeFactor = size_factors
  )
  return(size_factors)
}



#' When the maximum value (read count) of a peak between samples is under a certain threshold, we considere it is not relevant and the values are set to 0.
#' @param data Read count matrix.
#' @param minExp Minimum value. Default is 5.
#' @return Return a dataframe.
run_filter_low_counts <- function(data, minExp = 5) {
  col_names <- colnames(data)
  data <- t(apply(data, 1, filter_low_counts, col_names = col_names, minExp = minExp))
  data <- as.data.frame(data)
  data <- data[rowSums(data[]) > 0, ]
  return(data)
}
filter_low_counts <- function(row, col_names, minExp) {
  if (max(row) < minExp) {
    return(setNames(rep(0, length(row)), col_names))
  } else {
    return(setNames(row, col_names))
  }
}

#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

count_file <- snakemake@input[["counts"]]

minReads <- snakemake@params[["minReads"]]

###########################################
#                                         #
#              Get matrices               #
#                                         #
###########################################

raw_counts <- get_peak_matrix(
  count_file
)

###########################################
#                                         #
#             Filter low peaks            #
#                                         #
###########################################

# Remove peaks if max value < x reads
raw_counts <- run_filter_low_counts(raw_counts, minReads)

###########################################
#                                         #
#             Get samplesheet             #
#                                         #
###########################################

stage <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 1)
sex <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 2)
conditions <- paste(sex, stage, sep = " ")

replicate <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 4)
transgene <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 3)

samplesheet <- data.frame(
  sample = colnames(raw_counts),
  stages = stage,
  sex = sex,
  conditions = conditions,
  replicate = replicate,
  transgene = transgene
)

###########################################
#                                         #
#           Normalize read counts         #
#                                         #
###########################################

norm_counts <- get_normalized_counts(
  raw_counts,
  samplesheet
)

size_factors <- get_size_factors(
  raw_counts,
  samplesheet
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

write.csv(raw_counts, snakemake@output[["counts"]])
write.csv(norm_counts, snakemake@output[["norm_counts"]])
write.csv(samplesheet, snakemake@output[["samplesheet"]])
write.csv(size_factors, snakemake@output[["size_factors"]])
