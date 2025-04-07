'''
Author: Isabelle Stévant
Affiliation: Mammalian sex determination lab, University of Bar Ilan
Date: 24/09/2024
Licence: MIT


Snakemake rules concerning the Multi-omics analysis

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
processed_data = f'{config["path_to_process"]}{config["genome_version"]}/MULTI'
ATAC_processed_data = f'{config["path_to_process"]}{config["genome_version"]}/ATAC'
output_png = f'{config["path_to_graphs"]}{config["genome_version"]}/MULTI'
output_pdf = f'{config["path_to_graphs"]}{config["genome_version"]}/MULTI/PDF'
output_tables = f'{config["path_to_tables"]}{config["genome_version"]}/MULTI'
RNA_tables = f'{config["path_to_tables"]}{config["genome_version"]}/RNA'
ATAC_tables = f'{config["path_to_tables"]}{config["genome_version"]}/ATAC'

if config["genome_version"] == "mm10":
    genome = f'{input_data}/gencode.vM25.annotation.gtf.gz'
    ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm10"]
    ATAC_norm_bigwig_folder = config["ATAC_norm_bigwig_folder_mm10"]
else :
    genome = f'{input_data}/gencode.vM34.annotation.gtf.gz'
    ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm39"]
    ATAC_norm_bigwig_folder = config["ATAC_norm_bigwig_folder_mm39"]


# List of output figures
rule_MULTI_input_list = [
    f"{output_tables}/all_sig_gene2peak_linkage.csv",
    f"{output_png}/MULTI_linkage_stats.png",
    f"{processed_data}/plot_example_2.log",
    f"{processed_data}/plot_example_3.log"
#   f"{output_pdf}/MULTI_gene2peak_plots.pdf",
#   f"{output_png}/MULTI_TFBS_motifs_peak_XX_genes.png",
#   f"{output_png}/MULTI_TFBS_motifs_peak_XY_genes.png",
#   # f"{processed_data}/plot_example_1.log",
#   # f"{processed_data}/plot_example_2.log",
#   f"{processed_data}/plot_example_3.log"
]

## Uncomment if you want to run only this pipeline
# rule all:
#   input:
#       rule_MULTI_input_list


###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

# Calculate the correlation between gene expressin and ATAC peaks located within 500kb upstream and downstream of the gene TSS.
# Peaks on gene promoters are excluded to avoid kinking promoters of a neighbour gene as a putative enhancer.
rule MULTI_Get_all_gene_peak_correlation:
    input:
        RNA_samplesheet=f"{RNA_tables}/samplesheet.csv",
        RNA_norm_counts=f"{RNA_tables}/TPM.csv",
        ATAC_samplesheet=f"{ATAC_tables}/ATAC_samplesheet.csv",
        ATAC_norm_counts=f"{ATAC_tables}/ATAC_norm_counts.csv",
        chrom_size=f"{input_data}/chrom.size",
        genes=f"{input_data}/gene_standard.bed",
        gtf=f"{genome}"
    params:
        sex="all",
        distance=config["MULTI_peak_gene_distance"],
        min_cor=config["MULTI_peak_gene_min_cor"],
        FDR=config["MULTI_peak_gene_FDR"]
    output:
        linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
        bedpe=f"{output_tables}/all_sig_gene2peak_linkage.bedpe"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/MULTI_gene2peak.R"

###########################################
#                                         #
#                 Plots                   #
#                                         #
###########################################


rule MULTI_plot_linkage_stat:
    input:
        linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
    output:
        pdf=f"{output_pdf}/MULTI_linkage_stats.pdf",
        png=f"{output_png}/MULTI_linkage_stats.png"
    threads: 2
    resources:
        cpus_per_task=2,
        mem_mb=12000
    script:
        "../scripts/MULTI_plot_stat_linkage.R"


rule MULTI_Plot_peak_examples:
    input:
        genome=f"{genome}",
        gene_bed=f"{input_data}/gene_standard.bed",
        peaks=f"{ATAC_tables}/ATAC_norm_counts.csv",
        linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
        peak_list=f"{input_data}/gTrack_DAR_link_examples.tsv",
        TPM=f"{RNA_tables}/TPM.csv",
        bw=f"{ATAC_processed_data}/dl_bigwig.log"
    params:
        bw_folder=f"{ATAC_norm_bigwig_folder}",
        save_folder=f"{output_pdf}"
    output: 
        log= f"{processed_data}/plot_example_2.log"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/MULTI_plot_gene2peak_bis.R"


rule MULTI_Plot_PCHiC_examples:
    input:
        genome=f"{genome}",
        gene_bed=f"{input_data}/gene_standard.bed",
        peaks=f"{ATAC_tables}/ATAC_norm_counts.csv",
        linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
        PCHiC=f"{input_data}/PCHiC_5kb_score3_merged.bedpe",
        peak_list=f"{input_data}/gTrack_DAR_PCHiC_examples.tsv",
        TPM=f"{RNA_tables}/TPM.csv",
        bw=f"{ATAC_processed_data}/dl_bigwig.log"
    params:
        bw_folder=f"{ATAC_norm_bigwig_folder}",
        save_folder=f"{output_pdf}"
    output: 
        log= f"{processed_data}/plot_example_3.log"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/MULTI_plot_genomic_tracks_PCHiC.R"

# rule MULTI_TFBS_motifs_peak_XX_genes:
#     input:
#         linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
#         sex_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj",
#         TPM=f"{output_tables}/RNA_TPM.csv"
#     params:
#         background="conditions",
#         genome=config["genome_version"],
#         save_folder=f"{output_tables}",
#         sex="XX"
#     output:
#         pdf=f"{output_pdf}/MULTI_TFBS_motifs_peak_XX_genes.pdf",
#         png=f"{output_png}/MULTI_TFBS_motifs_peak_XX_genes.png"
#     resources:
#         cpus_per_task=12,
#         mem_mb=64000
#     script:
#         "../scripts/MULTI_linked_OCR_motif_enrichment.R"

# rule MULTI_TFBS_motifs_peak_XY_genes:
#     input:
#         linkage=f"{output_tables}/all_sig_gene2peak_linkage.csv",
#         sex_DEGs=f"{processed_data}/RNA_sig_SexDEGs.Robj",
#         TPM=f"{output_tables}/RNA_TPM.csv"
#     params:
#         background="conditions",
#         genome=config["genome_version"],
#         save_folder=f"{output_tables}",
#         sex="XY"
#     output:
#         pdf=f"{output_pdf}/MULTI_TFBS_motifs_peak_XY_genes.pdf",
#         png=f"{output_png}/MULTI_TFBS_motifs_peak_XY_genes.png"
#     resources:
#         cpus_per_task=12,
#         mem_mb=64000
#     script:
#         "../scripts/MULTI_linked_OCR_motif_enrichment.R"

