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
	# library("tidyverse")
	# library("MASS")
	# library("cisDynet")
})

doParallel::registerDoParallel(cores=4)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

# Modified version of the function from cisDynet

my_getPeak2Gene <- function(atac_matrix, rna_matrix, peak_annotation,
                          max_distance=20000, N_permutation=10000, save_path=NA){
  checkGeAnno()
  atac_paired_norm <- atac_matrix
  rna_paired <- rna_matrix
  logfile("Remove the gene with all expression value is 0.")
  rna_paired <- rna_paired[rowSums(rna_paired) > 0,]
  tss <- CATAnno$tss %>% as.data.frame()
  rownames(tss) <- tss$name
  tss_df <- tss[,c(1,2,3,4,6)]
  colnames(tss_df) <- c("chr","start","end","gene","strand")
  gr.tss <- GenomicRanges::makeGRangesFromDataFrame(tss_df,keep.extra.columns=TRUE)
  peak_bed <- peak_annotation
  peak_bed$Peak <- sprintf("%s:%s-%s",peak_bed$Chromosome, peak_bed$Start, peak_bed$End)
  rownames(peak_bed) <- NULL
  peak_gr <- GenomicRanges::makeGRangesFromDataFrame(peak_bed, keep.extra.columns=TRUE)

  ## Make the pseudo data for permutation
  logfile("Make the pseudo data for permutation...")
  test_random <- atac_paired_norm[sample(nrow(atac_paired_norm), N_permutation), ]


  ## Calculate the correlation coefficient and p-value for a given gene.
  calculate_pvalue <- function(row_data, df, gene, test_random=test_random){
    corr_random <- apply(test_random, 1, function(row){
      corr <- cor(as.numeric(row_data), row)
    })
    correlations <- apply(df, 1, function(row){
      corr <- cor(as.numeric(row_data), row)
      if (is.na(corr)){
        corr <- 0
        p <- 1
      } else {
        p <- BSDA::z.test(corr_random, mu = corr, sigma.x = 15)$p.value
      }
      return(list(corr = corr, p = p))
    })
    result <- data.frame(Peak = rownames(df), Gene = gene)
    result$correlations <- sapply(correlations, function(x) x$corr)
    result$p.value <- sapply(correlations, function(x) x$p)
    return(result)
  }

  ## get the gene2peak links function
  getGenePeaks <- function(gene, max_distance){
    gr.bm_gs <- unique(gr.tss[gr.tss$gene %in% gene])
    GenomicRanges::start(gr.bm_gs) <- GenomicRanges::start(gr.bm_gs) - max_distance
    GenomicRanges::end(gr.bm_gs) <- GenomicRanges::end(gr.bm_gs) + max_distance
    ol <- GenomicRanges::findOverlaps(peak_gr, gr.bm_gs)
    peak_in_window <- peak_gr[unique(IRanges::from(ol))]
    if(length(peak_in_window)!=0){
      df1 <- data.frame(chr=GenomicRanges::seqnames(peak_in_window), start=GenomicRanges::start(peak_in_window), end=GenomicRanges::end(peak_in_window))
      df1$peak <- sprintf("%s:%s-%s", df1$chr, df1$start, df1$end)
      ATAC <- atac_paired_norm[df1$peak,]
      RNA <- rna_paired[gene,]
      res <- calculate_pvalue(RNA, ATAC, gene, test_random=test_random)
      return(res)
    }else{
      return(NA)
    }
  }
  ##  calculate the all gene peak2gene links
  gene_list <- rownames(rna_paired)
  # pb <- progress::progress_bar$new(total = length(gene_list))
  result_list <- lapply(c(1:length(gene_list)), function(x) {
    result <- getGenePeaks(gene_list[x],max_distance=max_distance)
    # pb$tick()
    return(result)
  })
  na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
  res_list <- na.omit.list(result_list)
  all_res <- dplyr::bind_rows(res_list)
  final <- merge(all_res,peak_bed,by="Peak",all.x=T)
  final_res <- merge(final,tss_df, by.x="Gene", by.y="gene", all.x=T)
  final_res$Summit2TSS <- final_res$summit - final_res$start
  final_res$orientation <- ifelse((final_res$strand=="+" & final_res$Summit2TSS>=0)|(final_res$strand=="-" & final_res$Summit2TSS<0), "Downstream", "Upstream")
  final_res <- final_res[,c("Peak","Gene","correlations","p.value","Type","summit","start","Summit2TSS","strand","orientation")]
  colnames(final_res)[6:7] <- c("PeakSummit","TSS")
  final_res <- tibble::add_column(final_res, FDR = p.adjust(final_res$p.value, method = "fdr"), .after = 4)
  if(!is.na(save_path)){
    saveRDS(final_res, sprintf("%s/Peak2Gene_All_Links.rds",save_path))
  }
  return(final_res)
}

addAnnotation <- function(gene_bed, gtf, genome_size) {
  # gene <- valr::read_bed(gene_bed, n_fields = 6)
  gene <- valr::read_bed(gene_bed)  ## valr 0.7.0 update
  tss <- gene
  tss$start <- ifelse(tss$strand=="+", tss$start, tss$end)
  tss$end <- ifelse(tss$strand=="+", tss$start + 1, tss$end+1)
  genome_size <- valr::read_genome(genome_size)
  CATAnno <- list(gene, tss, genome_size, gtf)
  names(CATAnno) <- c("gene", "tss", "genome","gtf")
  assign("CATAnno", CATAnno, envir = .GlobalEnv)
}

annoMergedPeaks <- function(quant_data, tss_flank, cutoff, save_path=NA, save_name=NA){
  checkGeAnno()
  peak <- read.table(quant_data, header = T, row.names = 1)
  merged<- data.frame(peak=rownames(peak))
  merged_peaks <- merged %>% tidyr::separate(peak, c("chrom", "tem1"), ":") %>% tidyr::separate(tem1, c("start", "end"), "-")
  merged_peaks$start <- as.integer(merged_peaks$start)
  merged_peaks$end <- as.integer(merged_peaks$end)

  gene <- CATAnno$gene
  tss <- CATAnno$tss
  gene$start <- ifelse(gene$strand == "+", gene$start + tss_flank, gene$start)
  gene$end <- ifelse(gene$strand == "+", gene$end, gene$end - tss_flank)
  genome <- CATAnno$genome
  promoter <- valr::bed_flank(tss, genome, left = cutoff, right = tss_flank, strand = TRUE)

  merged_peaks$Start <- merged_peaks$start
  merged_peaks$End <- merged_peaks$end
  merged_peaks$start <- as.integer((merged_peaks$start + merged_peaks$end) / 2)
  merged_peaks$end <- merged_peaks$start + 1
  pro <- valr::bed_intersect(merged_peaks, promoter)[,c(1:5)]
  pro <-  pro[!duplicated(pro[,c(1,2,3)]),]
  pro$Type <- "Proximal"
  colnames(pro)[2:5] <- c("start","end","Start","End")
  intra <- valr::bed_intersect(merged_peaks, promoter,invert = TRUE)
  intra <-  intra[!duplicated(intra[,c(1,2,3)]),]
  intra_res <- valr::bed_intersect(intra, gene)[,c(1:5)]
  intra_res <- intra_res[!duplicated(intra_res[,c(1,2,3)]),]
  intra_res$Type <- "Intragenic"
  colnames(intra_res)[2:5] <- c("start","end","Start","End")
  re <- valr::bed_intersect(intra, intra_res,invert = TRUE)[,c(1:5)]
  re <-  re[!duplicated(re[,c(1,2,3)]),]
  re$Type <- "Distal"
  intra_pe <- rbind(re, intra_res, pro)
  intra_pe <- as.data.frame(intra_pe[!duplicated(intra_pe[,c(1,2,3)]), ])
  intra_pe <- intra_pe[,c(1,4,5,6,2)]
  colnames(intra_pe) <- c("Chromosome","Start","End","Type","summit")
  rownames(intra_pe) <- sprintf("%s:%s-%s", intra_pe$Chromosome, intra_pe$Start, intra_pe$End)
  if(!is.na(save_path)){
  write.table(sprintf("%s/%s_Merged_Peaks_Annotations.tsv", save_path, save_name), sep='\t', quote=F, col.names=T, row.names=T)
  }
  return(intra_pe)
}

###########################################
#                                         #
#         Prepare cisDynet files          #
#                                         #
###########################################

# install.packages("devtools")
# devtools::install_github("tzhu-bio/cisDynet")
# renv::install("tidyverse")
# renv::snapshot()
# genome_file <- "workflow/data/iGenome_mm10_ucsc_genes.gtf.gz"
# mm10Genes <- rtracklayer::import(genome_file)
# rtracklayer::export(mm10Genes, format="bed")
# genes <- subset(mm10Genes, type == "gene")
# txdb <- GenomicFeatures::makeTxDbFromGRanges(
# 	mm10Genes,
# 	drop.stop.codons=FALSE
# )
# genes <- genes(txdb)
# mm10_genes <- data.frame(
# 	seqnames(genes),
# 	start(genes),
# 	end(genes),
# 	genes$gene_id,
# 	rep(1,length(genes)),
# 	strand(genes)
# )
# test <- mm10_genes[grep("^chr[0-9|X|Y]{1,2}$", mm10_genes$seqnames.genes.),]
# write.table(test, file="mm10_gene_standard.bed", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

addAnnotation(
	gene_bed = "workflow/data/mm10_gene_standard.bed",
	gtf = "workflow/data/iGenome_mm10_ucsc_genes.gtf.gz",
	genome_size = "workflow/data/mm10.chrom.size"
)

anno <- annoMergedPeaks(
	quant_data = "results/processed_data/ATAC_norm_counts.csv", 
	cutoff=3000,
	tss_flank=1000
)

#################################################################################################################################

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

rna_samplesheet <- read.csv(snakemake@input[['RNA_samplesheet']], row.names=1)
atac_samplesheet <- read.csv(snakemake@input[['ATAC_samplesheet']], row.names=1)
# rna_samplesheet <- read.csv("results/processed_data/RNA_samplesheet.csv", row.names=1)
# atac_samplesheet <- read.csv("results/processed_data/ATAC_samplesheet.csv", row.names=1)

rna <- read.csv(snakemake@input[['RNA_norm_counts']], row.names=1)
# rna <- read.csv("results/processed_data/RNA_norm_counts.csv", row.names=1)

atac <- read.csv(snakemake@input[['ATAC_norm_counts']], row.names=1)
# atac <- read.csv("results/processed_data/ATAC_norm_counts.csv", row.names=1)

sex <- snakemake@params[['sex']]
# sex <- "XX"

distance <- snakemake@params[['distance']]

min_cor <- snakemake@params[['min_cor']]

max_FDR <- snakemake@params[['FDR']]


###########################################
#                                         #
#        Get gene-peak correlation        #
#                                         #
###########################################

# Mean expression of the replicates
rna <- rna[,grep(sex, colnames(rna))]
conditions <- unique(rna_samplesheet$conditions)
conditions <- grep(sex, conditions, value=TRUE)
samples <- lapply(conditions, function(cond) rna_samplesheet[rna_samplesheet$conditions==cond,1])
means <- lapply(samples, function(cond) rowMeans(rna[,cond]))
names(means) <- conditions
rna_means <- as.data.frame(do.call(cbind, means))
# split_RNA data into 4 tables to parallelize the work
split_rna <- split(rna_means, factor(sort(rank(row.names(rna_means)) %% 4)))

# split_rna <- split(rna, factor(sort(rank(row.names(rna)) %% 8)))

# Mean expression of the replicates
atac <- atac[,grep(sex, colnames(atac))]
conditions <- unique(atac_samplesheet$conditions)
conditions <- grep(sex, conditions, value=TRUE)
samples <- lapply(conditions, function(cond) atac_samplesheet[atac_samplesheet$conditions==cond,1])
means <- lapply(samples, function(cond) rowMeans(atac[,cond]))
names(means) <- conditions
atac_means <- as.data.frame(do.call(cbind, means))


p2g_split <- foreach(split=split_rna) %dopar% {
	# print(head(split))
	my_getPeak2Gene(
		atac_matrix = atac_means,
		rna_matrix = split,
		peak_annotation = anno,
		max_distance = distance,
		N_permutation = 10000
	)
}

p2g_res <- do.call(rbind.data.frame, p2g_split)

sig_p2g_res <- p2g_res[p2g_res$FDR<max_FDR,]
sig_p2g_res <- sig_p2g_res[abs(sig_p2g_res$correlations)>min_cor,]

write.table(sig_p2g_res, snakemake@output[['linkage']])
