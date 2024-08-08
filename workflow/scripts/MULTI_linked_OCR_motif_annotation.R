source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("monaLisa")
  if (snakemake@params[["genome"]] == "mm10") {
    library("BSgenome.Mmusculus.UCSC.mm10")
  } else {
    library("BSgenome.Mmusculus.UCSC.mm39")
  }
  library("dplyr")
  library("TFBSTools")
  library("motifStack")
})

#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

linked_OCRs <- "results/tables/mm10/all_sig_gene2peak_linkage.csv"

TPM <- read.csv(file = snakemake@input[["TPM"]], header = TRUE, row.names = 1)
# TPM <- read.csv(file="results/processed_data/mm10/RNA_TPM.csv", header=TRUE, row.names=1)


genome_version <- snakemake@params[["genome"]]
# genome_version <- "mm10"


JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db

#################################################################################################################################

peaks <- read.table(linked_OCRs, header = TRUE)

gr <- unique(GenomicRanges::GRanges(peaks$Peak))

# Get the peak sequences
if (genome_version == "mm10") {
  seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)
} else {
  seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm39, gr)
}

names(seqs) <- paste0(seqnames(gr), ":", start(gr), "-", end(gr))

# Get all vertebrate TF matrices
pfms <- TFBSTools::getMatrixSet(
  JASPAR,
  opts = list(
    matrixtype = "PFM",
    tax_group = "vertebrates"
  )
)

filtered_pfm_list <- pfms[sapply(pfms, function(pfm) toupper(pfm@name) %in% toupper(rownames(TPM)))]

unique_names <- make.unique(sapply(filtered_pfm_list, function(pfm) toupper(pfm@name)))
for (i in seq_along(filtered_pfm_list)) {
  filtered_pfm_list[[i]]@name <- unique_names[i]
  filtered_pfm_list[[i]]@name <- gsub("[.]", "_", filtered_pfm_list[[i]]@name)
  filtered_pfm_list[[i]]@name <- gsub("-", "_", filtered_pfm_list[[i]]@name)
}

pcms <- universalmotif::convert_motifs(filtered_pfm_list, class = "motifStack-pcm")
pcms_names <- unlist(lapply(pcms, function(pcm) pcm@name))
names(pcms) <- pcms_names
pcms <- universalmotif::convert_motifs(pcms, class = "motifStack-pcm")


hc <- motifStack::clusterMotifs(pcms)
phylog <- ade4::hclust2phylog(hc)

leaves <- names(phylog$leaves)
pcms <- as.list(pcms[leaves])


motifSig <- motifStack::motifSignature(pcms, phylog, cutoffPval = 0.005, min.freq = 1)

## get the signatures from object of motifSignature
sig <- signatures(motifSig)

motifs_pwms <- universalmotif::convert_motifs(sig, class = "TFBSTools-PWMatrix")
cl_vertebrate_pwm <- do.call(PWMatrixList, motifs_pwms)

save(cl_vertebrate_pwm, file = "cl_vertebrate_pwm.Rds")

res <- findMotifHits(
  query = cl_vertebrate_pwm,
  subject = seqs,
  min.score = 8.0,
  method = "matchPWM",
  BPPARAM = BiocParallel::MulticoreParam(12)
)

m <- table(
  factor(seqnames(res), levels = names(seqs)),
  factor(res$pwmname)
)

head(m)
