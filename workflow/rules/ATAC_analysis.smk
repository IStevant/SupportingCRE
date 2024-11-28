'''
Author: Isabelle Stévant
Affiliation: Mammalian sex determination lab, University of Bar Ilan
Date: 23/09/2024
Licence: MIT


Snakemake rules concerning the ATAC-seq analysis

The pipeline support mm10 and mm39 mouse genome versions
The parameters of the analsysis are defined in the analysis_parameters.yaml configuration file.

Pipeline created to analyse the paired time series RNA and ATAC-seq of the supporting cells from mouse fetal gonads.
Stévant et al. 2024
doi:

'''

configfile: "analysis_parameters.yaml"

sexes = config["sexes"]
TFBS_background = ["genome", "conditions"]

# Get file path according to the genome version
input_data = f'{config["path_to_data"]}{config["genome_version"]}'
processed_data = f'{config["path_to_process"]}{config["genome_version"]}/ATAC'
output_png = f'{config["path_to_graphs"]}{config["genome_version"]}/ATAC'
output_pdf = f'{config["path_to_graphs"]}{config["genome_version"]}/ATAC/PDF'
output_tables = f'{config["path_to_tables"]}{config["genome_version"]}/ATAC'
RNA_tables = f'{config["path_to_tables"]}{config["genome_version"]}/RNA'
RNA_processed_data = f'{config["path_to_process"]}{config["genome_version"]}/RNA'



if config["genome_version"] == "mm10":
    genome = f'{input_data}/gencode.vM25.annotation.gtf.gz'
    ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm10"]
    ATAC_norm_bigwig_folder = config["ATAC_norm_bigwig_folder_mm10"]
else :
    genome = f'{input_data}/gencode.vM34.annotation.gtf.gz'
    ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm39"]
    ATAC_norm_bigwig_folder = config["ATAC_norm_bigwig_folder_mm39"]

# List of output files
rule_ATAC_input_list = [
    f"{output_png}/ATAC_corr_pca_all_samples.png",
    f"{output_png}/ATAC_consensus_peak_distribution.png",
    f"{output_png}/ATAC_all_consensus_peak_annotation.png",
    f"{output_png}/ATAC_sex_DAR_histograms.png",
    f"{output_png}/ATAC_sig_sex_DARs_annotation.png",
    f"{output_png}/ATAC_sex_DAR_upset.png",
    # expand(f"{output_png}/ATAC_sex_DAR_TF_motifs_{{background}}_bg.png", background = TFBS_background),
    f"{output_png}/ATAC_sex_DAR_TF_motifs_genome_bg.png",
    f"{output_png}/ATAC_sex_DAR_TF_motifs_merged.png",
    expand(f"{output_png}/ATAC_{{sex}}_DAR_stage_heatmap.png", sex=sexes),
    f"{output_png}/ATAC_common_dynamic_DARs.png",
    f"{output_png}/ATAC_common_sex_dynamic_DARs.png",
    f"{output_png}/ATAC_sig_stage_DARs_annotation.png",
    expand(f"{output_png}/ATAC_{{sex}}_dyn_DAR_TF_motifs_{{background}}_bg.png", sex=sexes, background = TFBS_background),
    f"{processed_data}/plot_example_1.log"
]

## Uncomment if you want to run only this pipeline
# rule all:
#   input:
#       rule_ATAC_input_list


###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

# Generate the expression matrices that will be used for the rest of the analysis
rule ATAC_Get_matrices:
    input:
        counts=f'{input_data}/{config["ATAC_counts"]}', # Raw read count per peak comming from the nf-core/atacseq pipeline
    params:
        minReads=config["ATAC_minReads"],  # Minimal number of raw reads from which we consider a peak accessible
    output:
        counts=f"{output_tables}/ATAC_raw_counts.csv",          # Filtered raw read counts
        norm_counts=f"{output_tables}/ATAC_norm_counts.csv",    # Filtered normalized read counts (normalization by the library size)
        samplesheet=f"{output_tables}/ATAC_samplesheet.csv",    # Description of the samples for downstream analysis
        size_factors=f"{output_tables}/ATAC_size_factors.csv",  # Size factors used to normalize the data (was used to normalize the BigWig files)
        bed=f"{output_tables}/ATAC_all_OCRs.bed"                # Bed files with the merged OCRs coordinates
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/ATAC_clean_matrices.R"

# BigWig files from the nf-core/atacseq pipeline are normalized by read per million.
# This normalization is not appropriate when samples have different FRiP scores (i.e. different level of background noise).
# We decided to normalize the BigWig using the same normalization applied for the peaks in the differential accessibility analysis.
# The rational is that most of the peaks are unchanged between conditions, so the signal can be normalized using the peak signals only.
# We provide the normalized bigwig files on GEO.
# Force the execution if you need to re-normalize the bigwig files (not recommanded, apply only on output nf-core/atacseq files normalized by read per million).
rule ATAC_Normalize_bigwig:
    input:
        size_factors=f"{output_tables}/ATAC_size_factors.csv"   # Size factors used to normalize the data
    params:
        bigwig_folder=f"{ATAC_bigwig_folder}",             # Path to the original nf-core/atacseq BigWig files
        new_bigwig_folder=f"{ATAC_norm_bigwig_folder}"  # Path to the normalized files
    output:
        output_file=f"{processed_data}/ATAC_bigwig/scale/ATAC_size_factors.csv"  # Save the size factors with the new BigWig (to signify snakemake that the rule executed well)
    threads: 12
    resources:
        mem_mb=16000
    script:
        "../scripts/MULTI_norm_bigwig.R"

# Get differentially accessible regions between XX and XY at each embryonic stage using DESeq2 Wald test
rule ATAC_Get_sex_DARs:
    input:
        counts=f"{output_tables}/ATAC_raw_counts.csv",        # Filtered raw read counts
        samplesheet=f"{output_tables}/ATAC_samplesheet.csv",  # Description of the samples
        gtf=f"{genome}"                                       # GTF file
    params:
        adjpval=config["ATAC_adjpval"],             # Adjusted p-value threshold to condider a region differentially accessible
        log2FC=config["ATAC_log2FC"],               # Log2 fold change threshold to condider a region differentially accessible
        promoter=config["ATAC_promoter_distance"],  # distance to the TSS defining the promoter region
        save_folder=f"{output_tables}"              # Location where the result tables are saved
    output:
        sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",      # R object containing the filtered DESeq2 results
        sig_DARs_GR=f"{processed_data}/ATAC_sig_SexDARs_GR.Robj",# R object containing DAR regions as a GRanges list
        sig_DARs_bed=f"{output_tables}/ATAC_sig_SexDARs.bed"
    threads: 12
    resources:
        mem_mb=16000
    script:
        "../scripts/ATAC_sex_DAR.R"


# Get Jaspar2024 TF motifs, select the TFs expressed in the gonad using our RNA-seq, and merge ythe motofs by similarity
rule ATAC_Merge_TF_motifs:
    input:
        TPM=f"{RNA_tables}/TPM.csv",
        TF_genes=config["TF_genes"]
    output:
        matrices_txt=f"{processed_data}/merged_exp_TF_matrices_pfms.txt",
        matrices_obj=f"{processed_data}/merged_exp_TF_matrices_pfms.RData"
    params:
        minTPM=config["RNA_minTPM"]
    threads: 1
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_merge_TF_motifs.R"

# Get TFBS enrichment in sex-biased OCRs against random genomic background and against the opposite sex as background.
# Returns the plots and the enrichment scores as tables.
rule ATAC_TFBS_motifs_sex_DAR:
    input:
        TF_genes=config["TF_genes"],                         # List of known mouse transcripotion factors
        sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",  # R object containing the filtered DESeq2 results
        TPM=f"{RNA_tables}/TPM.csv",                          # Gene expression matrix
        DEG=f"{RNA_processed_data}/sig_SexDEGs.Robj" 
    output:
        # pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_{{background}}_bg.pdf",  # Figure as PDF
        # png=f"{output_png}/ATAC_sex_DAR_TF_motifs_{{background}}_bg.png"   # Figure as PNG
        pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_genome_bg.pdf",  # Figure as PDF
        png=f"{output_png}/ATAC_sex_DAR_TF_motifs_genome_bg.png"   # Figure as PNG
    params:
        minTPM=config["RNA_minTPM"],                        # Minimal number of TPM from which we consider a gene expressed
        # background=lambda wildcards: wildcards.background,  # Background to use to do the enrichment analysis (either "genome" or "conditions")
        background="genome",
        genome=config["genome_version"],                    # Version of the genome ("mm10" or "mm39")
        save_folder=f"{output_tables}"                      # Location where the result tables are saved
    threads: 48
    resources:
        mem_mb=84000
    script:
        "../scripts/ATAC_sex_DAR_motif_enrich_bis.R"

# Get TFBS enrichment in sex-biased OCRs from all stages together against the opposite sex as background.
rule ATAC_TFBS_motifs_sex_DAR_merged_stages:
    input:
        TF_genes=config["TF_genes"],                         # List of known mouse transcripotion factors
        sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",  # R object containing the filtered DESeq2 results
        TPM=f"{RNA_tables}/TPM.csv"                          # Gene expression matrix
    params:
        minTPM=config["RNA_minTPM"],                         # Minimal number of TPM from which we consider a gene expressed
        background="conditions",                             # Background to use to do the enrichment analysis (either "genome" or "conditions")
        genome=config["genome_version"],                     # Version of the genome ("mm10" or "mm39")
        save_folder=f"{output_tables}"                       # Location where the result tables are saved
    output:
        matrices_txt=f"{processed_data}/sex_merged_motifs_pfms.txt",
        pdf=f"{output_pdf}/ATAC_sex_DAR_TF_motifs_merged.pdf",
        png=f"{output_png}/ATAC_sex_DAR_TF_motifs_merged.png",
        heatmap_matrice=f"{processed_data}/merged_enr_matrix.RData"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_sex_DAR_motif_enrich_merge_stages.R"


rule ATAC_Get_XX_dynamic_DARs:
    input:
        counts=f"{output_tables}/ATAC_raw_counts.csv",
        samplesheet=f"{output_tables}/ATAC_samplesheet.csv",
        peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj",
        gtf=f"{genome}"
    params:
        adjpval=config["ATAC_adjpval"],
        log2FC=config["ATAC_log2FC"],
        promoter=config["ATAC_promoter_distance"],
        sex="XX"
    output:
        tsv=f"{output_tables}/ATAC_XX_DAR_stage.tsv",
        sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XX.Robj"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_stage_DAR.R"

rule ATAC_Get_XY_dynamic_DARs:
    input:
        counts=f"{output_tables}/ATAC_raw_counts.csv",
        samplesheet=f"{output_tables}/ATAC_samplesheet.csv",
        peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj",
        gtf=f"{genome}"
    params:
        adjpval=config["ATAC_adjpval"],
        log2FC=config["ATAC_log2FC"],
        promoter=config["ATAC_promoter_distance"],
        sex="XY"
    output:
        tsv=f"{output_tables}/ATAC_XY_DAR_stage.tsv",
        sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XY.Robj"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_stage_DAR.R"


# Get TFBS enrichment in sex-biased OCRs against random genomic background and against the opposite sex as background.
# Returns the plots and the enrichment scores as tables.
rule ATAC_TFBS_motifs_dynamic_DAR:
    input:
        TF_genes=config["TF_genes"],                         # List of known mouse transcripotion factors
        sig_DARs=f"{output_tables}/ATAC_{{sex}}_DAR_stage_heatmap_clusters.csv",  # R object containing the filtered DESeq2 results
        TPM=f"{RNA_tables}/TPM.csv"                          # Gene expression matrix
    output:
        pdf=f"{output_pdf}/ATAC_{{sex}}_dyn_DAR_TF_motifs_{{background}}_bg.pdf",  # Figure as PDF
        png=f"{output_png}/ATAC_{{sex}}_dyn_DAR_TF_motifs_{{background}}_bg.png"   # Figure as PNG
    params:
        minTPM=config["RNA_minTPM"],                        # Minimal number of TPM from which we consider a gene expressed
        background=lambda wildcards: wildcards.background,  # Background to use to do the enrichment analysis (either "genome" or "conditions")
        logos="TRUE",
        genome=config["genome_version"],                    # Version of the genome ("mm10" or "mm39")
        save_folder=f"{output_tables}",                     # Location where the result tables are saved
        sex=lambda wildcards: wildcards.sex
    threads: 48
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_stage_DAR_motif_enrich.R"


###########################################
#                                         #
#                 Plots                   #
#                                         #
###########################################

rule ATAC_Plot_nb_consensus_peak:
    input:
        peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj",
        samplesheet=f"{output_tables}/ATAC_samplesheet.csv",
        norm_data=f"{output_tables}/ATAC_norm_counts.csv"
    output:
        pdf=f"{output_pdf}/ATAC_consensus_peak_distribution.pdf",
        png=f"{output_png}/ATAC_consensus_peak_distribution.png"
    threads: 12
    resources:
        mem_mb=16000
    script:
        "../scripts/ATAC_peak_distribution_per_sex.R"

rule ATAC_Plot_consensus_peak_annotation:
    input:
        genome=f"{genome}",
        peak_list=f"{input_data}/ATAC_all_consensus_peaks_2rep_list.Robj"
    params:
        promoter=config["ATAC_promoter_distance"]
    output:
        anno_list=f"{processed_data}/ATAC_all_consensus_peak_annotation.Robj",
        pdf=f"{output_pdf}/ATAC_all_consensus_peak_annotation.pdf",
        png=f"{output_png}/ATAC_all_consensus_peak_annotation.png"
    threads: 12
    resources:
        mem_mb=16000
    script:
        "../scripts/ATAC_peak_annotation_per_sex.R"

rule ATAC_corr_PCA:
    input:
        norm_data=f"{output_tables}/ATAC_norm_counts.csv"
    params:
        corr_method=config["ATAC_corr_met"]
    output:
        pdf=f"{output_pdf}/ATAC_corr_pca_all_samples.pdf",
        png=f"{output_png}/ATAC_corr_pca_all_samples.png"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/Corr_pca.R"

rule ATAC_Plot_sex_DAR_histogram:
    input:
        sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",
        samplesheet=f"{output_tables}/ATAC_samplesheet.csv"
    output:
        pdf=f"{output_pdf}/ATAC_sex_DAR_histograms.pdf",
        png=f"{output_png}/ATAC_sex_DAR_histograms.png"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_plot_sex_DAR_hist.R"

rule ATAC_Plot_sex_DAR_peak_annotation:
    input:
        genome=f"{genome}",
        peak_list=f"{processed_data}/ATAC_sig_SexDARs_GR.Robj"
    params:
        promoter=config["ATAC_promoter_distance"]
    output:
        anno_list=f"{processed_data}/ATAC_sig_SexDARs_annotation.Robj",
        pdf=f"{output_pdf}/ATAC_sig_sex_DARs_annotation.pdf",
        png=f"{output_png}/ATAC_sig_sex_DARs_annotation.png"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_peak_annotation_per_sex.R"

rule ATAC_Plot_sex_DAR_upset:
    input:
        sig_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj"
    params:
        output_folder=f"{output_tables}/"      # Location where the result tables are saved
    output:
        pdf=f"{output_pdf}/ATAC_sex_DAR_upset.pdf",
        png=f"{output_png}/ATAC_sex_DAR_upset.png",
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_sex_DAR_upset.R"

rule ATAC_Plot_heatmap_dyn_DARs:
    input:
        sig_DARs=f"{processed_data}/ATAC_sig_stage_DARs_{{sex}}.Robj",
        norm_counts=f"{output_tables}/ATAC_norm_counts.csv",
        samplesheet=f"{output_tables}/ATAC_samplesheet.csv"
    params:
        sex=lambda wildcards: wildcards.sex,
        clusters=config["ATAC_stage_DAR_clusters"]
    output:
        clusters=f"{output_tables}/ATAC_{{sex}}_DAR_stage_heatmap_clusters.csv",
        pdf=f"{output_pdf}/ATAC_{{sex}}_DAR_stage_heatmap.pdf",
        png=f"{output_png}/ATAC_{{sex}}_DAR_stage_heatmap.png"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_stage_DAR_heatmap.R"

# Draw venn diagrams of the comparison of dynamic DARs of both sexes
rule RNA_Plot_common_dynamic_DARs:
    input:
        XX_stage_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XX.Robj",    # Robj with the filtered XX dynamic genes
        XY_stage_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XY.Robj",    # Robj with the filtered XY dynamic genes
        samplesheet=f"{output_tables}/ATAC_samplesheet.csv"               # Description of the samples
    params:
        output_folder=f"{output_tables}/"      # Location where the result tables are saved
    output:
        pdf=f"{output_pdf}/ATAC_common_dynamic_DARs.pdf",
        png=f"{output_png}/ATAC_common_dynamic_DARs.png"
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/ATAC_overlap_dynamic_DAR.R"


# Draw venn diagrams of the comparison of dynamic and sex-specific DARs for each sex
rule ATAC_Plot_common_sex_dynamic_DARs:
    input:
        sex_DARs=f"{processed_data}/ATAC_sig_SexDARs.Robj",               # Robj with the filtered sex DARs per stage
        XX_stage_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XX.Robj",    # Robj with the filtered XX dynamic regions
        XY_stage_DARs=f"{processed_data}/ATAC_sig_stage_DARs_XY.Robj",     # Robj with the filtered XY dynamic regions
        samplesheet=f"{output_tables}/ATAC_samplesheet.csv"               # Description of the samples
    params:
        output_folder=f"{output_tables}/"      # Location where the result tables are saved
    output:
        pdf=f"{output_pdf}/ATAC_common_sex_dynamic_DARs.pdf",
        png=f"{output_png}/ATAC_common_sex_dynamic_DARs.png"
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/ATAC_overlap_sex_dynamic_DAR.R"

rule ATAC_Plot_dynamic_DAR_peak_annotation:
    input:
        genome=f"{genome}",
        XX_peak_list=f"{output_tables}/ATAC_XX_DAR_stage_heatmap_clusters.csv",
        XY_peak_list=f"{output_tables}/ATAC_XY_DAR_stage_heatmap_clusters.csv",
    params:
        promoter=config["ATAC_promoter_distance"]
    output:
        pdf=f"{output_pdf}/ATAC_sig_stage_DARs_annotation.pdf",
        png=f"{output_png}/ATAC_sig_stage_DARs_annotation.png"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/ATAC_dynamic_peak_annotation.R"

rule ATAC_Plot_peak_examples:
    input:
        genome=f"{genome}",
        peak_list=f"{input_data}/gTrack_DAR_peak_examples.tsv",
        peaks=f"{output_tables}/ATAC_norm_counts.csv",
        TPM=f"{RNA_tables}/TPM.csv"
    params:
        bw_folder=f"{ATAC_norm_bigwig_folder}",
        save_folder=f"{output_pdf}"
    output: 
        log= f"{processed_data}/plot_example_1.log"
    threads: 12
    resources:
        cpus_per_task=12,
        mem_mb=64000
    script:
        "../scripts/ATAC_plot_DAR_peak_examples.R"