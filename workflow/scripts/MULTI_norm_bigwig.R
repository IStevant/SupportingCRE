source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################
suppressPackageStartupMessages({
	library("dplyr")
	library('doParallel')
	library('foreach')
})

doParallel::registerDoParallel(cores=12)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

bigwig_folder <- snakemake@param[['bigwig_folder']]
size_factors <- snakemake@input[['size_factors']]
new_bigwig_folder <- snakemake@param[['new_bigwig_folder']]

# bigwig_folder <- "/home/istevant/work/data/ATACseq/240521_XX-enh8_XY-sox9_mm10_sub/mapping_nf-core_atacseq_v2.1.2_GRCm39/bwa/merged_library/bigwig"
# size_factors <- "/home/istevant/work/project/SupportingCRE/results/processed_data/mm39/ATAC_size_factors.csv"
# new_bigwig_folder <- "/home/istevant/work/project/SupportingCRE/results/processed_data/mm39/bigwig"

pattern <- "*.mLb.clN.bigWig"

bigwig_files <- list.files(path=bigwig_folder, pattern = pattern)
samples <- gsub(".mLb.clN.bigWig","",bigwig_files)

size_factors <- read.csv(size_factors, row.names=1)
rownames(size_factors) <- gsub("Sox9.IRES.GFP","Sox9-IRES-GFP", rownames(size_factors))
rownames(size_factors) <- gsub("Enh8.","Enh8-", rownames(size_factors))

norm_bigwig <- foreach(sample=samples) %dopar% {
	bw <- paste0(sample, ".mLb.clN.bigWig")
	# print(paste(bigwig_folder, bw, sep="/"))
	bw_data <- rtracklayer::import.bw(paste(bigwig_folder, bw, sep="/"))
	# print(bw_data)
	scale_folder <- paste(bigwig_folder, "scale", sep="/")
	scale_file <- list.files(path=scale_folder, pattern = paste0(sample, ".*"))
	scale <- read.csv(paste(scale_folder, scale_file, sep="/"), header=FALSE)[1,1]
	denorm_bw <- bw_data
	denorm_bw$score <- bw_data$score/scale
	sizeFactor <- size_factors[grep(tolower(sample), tolower(rownames(size_factors))), "SizeFactor"]
	norn_bw <- denorm_bw
	norn_bw$score <- denorm_bw$score*(1/sizeFactor)
	new_bw <- paste0(sample, ".bw")
	rtracklayer::export.bw(norn_bw, paste(new_bigwig_folder, new_bw, sep="/"))
}
