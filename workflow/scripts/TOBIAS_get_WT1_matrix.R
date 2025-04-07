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

# Load Jaspar 2024 database
JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db

###########################################
#                                         #
#             Get TF matrices             #
#                                         #
###########################################

# Get all vertebrate TF matrices
pfms <- TFBSTools::getMatrixByName(JASPAR, "Wt1")

pcms <- universalmotif::convert_motifs(pfms, class = "motifStack-pcm")

output_file <- snakemake@output[["matrix"]]

fileConn <- file(output_file, "w")

pfm <- pcms
tf_name <- pfm$name
writeLines(paste0(">", tf_name), fileConn)
pfm_matrix <- pfm$mat
bases <- rownames(pfm_matrix)
for (j in 1:nrow(pfm_matrix)) {
  pfm_values <- paste(pfm_matrix[j, ], collapse = "   ")
  writeLines(paste(bases[j], "[", pfm_values, "]"), fileConn)
}


close(fileConn)
