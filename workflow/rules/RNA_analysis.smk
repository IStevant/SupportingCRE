'''
Author: Isabelle Stévant
Affiliation: Mammalian sex determination lab, University of Bar Ilan
Date: 19/09/2024
Licence: MIT


Snakemake rules concerning the RNA-seq analysis

The parameters of the analsysis are defined in the analysis_parameters.yaml configuration file.

Pipeline created to analyse the paired time series RNA and ATAC-seq of the supporting cells from mouse fetal gonads.
Stévant et al. 2024
doi:

'''
configfile: "analysis_parameters.yaml"

sexes = config["sexes"]

# Get file path according to the genome version
input_data = f'{config["path_to_data"]}{config["genome_version"]}'
processed_data = f'{config["path_to_process"]}{config["genome_version"]}/RNA'
output_png = f'{config["path_to_graphs"]}{config["genome_version"]}/RNA'
output_pdf = f'{config["path_to_graphs"]}{config["genome_version"]}/RNA/PDF'
output_tables = f'{config["path_to_tables"]}{config["genome_version"]}/RNA'

# List of output figures
rule_RNA_input_list = [
    f"{output_png}/RNA_corr_pca_all_samples.png",
    f"{output_png}/RNA_corr_pca.png",
    f"{output_png}/RNA_marker_genes.png",
    f"{output_png}/RNA_marker_genes_enrichment.png",
    f"{output_png}/RNA_sex_DEG_histograms.png",
    f"{output_png}/RNA_sex_DEG_volcano.png",
    f"{output_png}/RNA_sex_DEG_upset.png",
    expand(f"{output_png}/RNA_{{sex}}_DEG_stage_heatmap.png", sex=sexes),
    f"{output_png}/RNA_common_dynamic_DEGs.png",
    f"{output_png}/RNA_common_sex_dynamic_DEGs.png"
]

# If there is no outliers, do not run the analysis that discard them
if len(config["RNA_outliers"])<1:
    rule_RNA_input_list.remove(f"{output_png}/RNA_corr_pca_all_samples.png")


## Uncomment if you want to run only this pipeline
# rule all:
#   input:
#       rule_RNA_input_list


###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

# Generate the expression matrices that will be used for the rest of the analysis
rule RNA_Get_matrices:
    input:
        counts=f'{input_data}/{config["RNA_counts"]}', # Raw read count per gene matrix comming from the nf-core/rnaseq pipeline
        tpm=f'{input_data}/{config["RNA_TPM"]}',       # TPM per gene matrix comming from the nf-core/rnaseq pipeline
        protein_genes=config["protein_genes"]          # List of the mouse protein coding genes
    params:
        minReads=config["RNA_minReads"],     # Minimal number of raw reads from which we consider a gene expressed
        minTPM=config["RNA_minTPM"],         # Minimal number of TPM from which we consider a gen eexpressed
        RNA_outliers=config["RNA_outliers"]  # If we know in advance we have outlier samples, generate the matrices with and without the outliers
    output:
        tpm_all=f"{output_tables}/TPM_all_samples.csv",           # TPM matrix with all the samples
        tpm=f"{output_tables}/TPM.csv",                           # TPM matrix without the outliers
        counts=f"{output_tables}/raw_counts.csv",                 # Filtered raw read counts
        norm_counts=f"{output_tables}/norm_counts.csv",           # Filtered normalized read counts (normalization by the library size)
        norm_counts_all=f"{output_tables}/norm_counts_all_samples.csv",   # Filtered normalized read counts without the outliers (normalization by the library size)
        samplesheet=f"{output_tables}/samplesheet.csv"            # Description of the samples for downstream analysis
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/RNA_clean_matrices.R"

# Get genes differentially expressed between XX and XY at each embryonic stage using DESeq2 Wald test
rule RNA_Get_sex_DEGs:
    input:
        counts=f"{output_tables}/raw_counts.csv",        # Filtered raw read counts
        samplesheet=f"{output_tables}/samplesheet.csv",  # Description of the samples
        TF_genes=config["TF_genes"],                     # List of known mouse transcription factors
        TF_pheno=config["TF_pheno"]                      # List ig genes with gonadal phenotypes from MGI OBO database
    params:
        adjpval=config["RNA_adjpval"],    # Adjusted p-value threshold to condider a gene differentially expressed
        log2FC=config["RNA_log2FC"],      # Log2 fold change threshold to condider a gene differentially expressed
        save_folder=f"{output_tables}"    # Location where the result tables are saved
    output:
        all_DEGs=f"{processed_data}/all_SexDEGs.Robj",  # R object containing the full DESeq2 results
        sig_DEGs=f"{processed_data}/sig_SexDEGs.Robj"   # R object containing the only the significant DESeq2 results
    threads: 12
    resources:
        mem_mb=12000
    script:
        "../scripts/RNA_sex_DEG.R"

# Get genes differentially expressed between embryonic stages using DESeq2 LTR test
rule RNA_Get_dynamic_DEGs:
    input:
        counts=f"{output_tables}/raw_counts.csv",        # Filtered raw read counts
        samplesheet=f"{output_tables}/samplesheet.csv",  # Description of the samples
        TF_genes=config["TF_genes"],                     # List of known mouse transcription factors
        TF_pheno=config["TF_pheno"]                      # List ig genes with gonadal phenotypes from MGI OBO database
    output:
        tsv=f"{output_tables}/{{sex}}_DEG_stage.tsv",             # Filtered results exported as a tab separated value table
        sig_DEGs=f"{processed_data}/sig_stage_DEGs_{{sex}}.Robj"  # Filtered results exported as an R object
    params:
        adjpval=config["RNA_adjpval"],        # Adjusted p-value threshold to condider a gene differentially expressed
        log2FC=config["RNA_log2FC"],          # Log2 fold change threshold to condider a gene differentially expressed
        sex=lambda wildcards: wildcards.sex   # Run the analysis for each sex
    threads: 12
    resources:
        mem_mb=12000
    script:
        "../scripts/RNA_dynamic_DEG.R"


###########################################
#                                         #
#                 Plots                   #
#                                         #
###########################################

# Draw the correlation matrix between samples and the PCA for all samples, including the outliers
rule RNA_corr_PCA_with_outliers:
    input:
        norm_data=f"{output_tables}/norm_counts_all_samples.csv"   # Filtered normalized read counts (normalization by the library size)
    params:
        corr_method=config["RNA_corr_met"]     # Correlation method (can be either "Pearson" or "Spearman")
    output:
        pdf=f"{output_pdf}/RNA_corr_pca_all_samples.pdf",  # Figure as PDF
        png=f"{output_png}/RNA_corr_pca_all_samples.png"   # Figure as PNG
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/Corr_pca.R"

# Draw the correlation matrix between samples and the PCA without the outliers
rule RNA_corr_PCA:
    input:
        norm_data=f"{output_tables}/norm_counts.csv"   # Filtered normalized read counts (normalization by the library size)
    params:
        corr_method=config["RNA_corr_met"]     # Correlation method (can be either "Pearson" or "Spearman")
    output:
        pdf=f"{output_pdf}/RNA_corr_pca.pdf",  # Figure as PDF
        png=f"{output_png}/RNA_corr_pca.png"   # Figure as PNG
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/Corr_pca.R"

# Draw marker gene expression over time compared to whole gonad RNA-seq data
rule RNA_Plot_marker_genes:
    input:
        marker_genes=config["marker_genes"],                          # List of gonadal marker genes and their associated cell type
        whole_gonad=f'{input_data}/{config["whole_gonad_RNAseq"]}',   # Expression matrix of whole gonad RNA-seq (TPMs)
        tpm=f"{output_tables}/TPM.csv"                                # TPM matrix without the outliers
    output:
        pdf=f"{output_pdf}/RNA_marker_genes.pdf",
        png=f"{output_png}/RNA_marker_genes.png"
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/RNA_plot_marker_genes.R"

# Dot plot showing the expression enrichment of marker genes compared to whole gonads
rule RNA_Plot_marker_gene_enrichment:
    input:
        whole_gonad=f'{input_data}/{config["whole_gonad_RNAseq"]}',   # Expression matrix of whole gonad RNA-seq (TPMs)
        tpm=f"{output_tables}/TPM.csv"                                # TPM matrix without the outliers
    output:
        pdf=f"{output_pdf}/RNA_marker_genes_enrichment.pdf",
        png=f"{output_png}/RNA_marker_genes_enrichment.png"
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/RNA_plot_marker_gene_enrichment.R"

# Draw the horizontal histogram showing the number of XX and XY overexpressed genes at each stage
rule RNA_Plot_sex_DEG_histogram:
    input:
        sig_DEGs=f"{processed_data}/sig_SexDEGs.Robj",   # Robj containing the filtered DEGs
        samplesheet=f"{output_tables}/samplesheet.csv"   # Description of the samples
    output:
        pdf=f"{output_pdf}/RNA_sex_DEG_histograms.pdf",
        png=f"{output_png}/RNA_sex_DEG_histograms.png"
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/RNA_plot_sex_DEG_hist.R"

# Draw the sex differentially expressed genes as Volcano plots with the top enriched GO terms side by side
rule RNA_Plot_sex_DEG_volcano_GO:
    input:
        all_DEGs=f"{processed_data}/all_SexDEGs.Robj",   # Robj containing the unfiltered DEGs
        samplesheet=f"{output_tables}/samplesheet.csv"   # Description of the samples
    params:
        adjpval=config["RNA_adjpval"],    # Adjusted p-value threshold to condider a gene differentially expressed
        log2FC=config["RNA_log2FC"],      # Log2 fold change threshold to condider a gene differentially expressed
        path=f"{output_tables}"           # Location where the result tables are saved
    output:
        pdf=f"{output_pdf}/RNA_sex_DEG_volcano.pdf",
        png=f"{output_png}/RNA_sex_DEG_volcano.png"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/RNA_plot_sex_volcano_GO.R"

# Get how many sex DEGs are found in other embryonic stages and draw an upset plot.
# Return plots and tables with the genes found in each intersection.
rule RNA_Plot_sex_DEG_upset:
    input:
        sig_DEGs=f"{processed_data}/sig_SexDEGs.Robj"   # Robj containing the filtered DEGs
    params:
        output_folder=f"{output_tables}/"      # Location where the result tables are saved
    output:
        pdf=f"{output_pdf}/RNA_sex_DEG_upset.pdf",
        png=f"{output_png}/RNA_sex_DEG_upset.png"
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/RNA_sex_DEG_upset.R"

# Draw the heatmap of the dynamically expressed genes and top GO term enrichment for each cluster side by side.
# The heatmap is annotated with 25 TFs with gonadal phenotypes.
# Return the genes in each clusters and the GO term enrichment result table.
rule RNA_Plot_heatmap_GO:
    input:
        TF_genes=config["TF_genes"],                          # List of known mouse transcription factors
        TF_pheno=config["TF_pheno"],                          # List ig genes with gonadal phenotypes from MGI OBO database
        sig_DEGs=f"{processed_data}/sig_stage_DEGs_{{sex}}.Robj",  # Robj with the filtered dynamic genes
        norm_counts=f"{output_tables}/norm_counts.csv",       # Filtered normalized read counts (normalization by the library size)
        samplesheet=f"{output_tables}/samplesheet.csv"        # Description of the samples
    params:
        sex=lambda wildcards: wildcards.sex,               # Current sex
        clusters=config["RNA_stage_DEG_clusters"]            # Minimal number of clusters (decided by visual inspection of the heatmap)
    output:
        GO=f"{output_tables}/{{sex}}_GO_DEG_stage.csv",                          # Simplified GO term enrichment result table
        cluster_file=f"{output_tables}/{{sex}}_DEG_stage_heatmap_clusters.csv",  # Genes per clusters
        pdf=f"{output_pdf}/RNA_{{sex}}_DEG_stage_heatmap.pdf",
        png=f"{output_png}/RNA_{{sex}}_DEG_stage_heatmap.png"
    threads: 12
    resources:
        mem_mb=16000
    script:
        "../scripts/RNA_dynamic_DEG_heatmap.R"

# Draw venn diagrams of the comparison of dynamic DEGs of both sexes
rule RNA_Plot_common_dynamic_DEGs:
    input:
        XX_stage_DEGs=f"{processed_data}/sig_stage_DEGs_XX.Robj",    # Robj with the filtered XX dynamic genes
        XY_stage_DEGs=f"{processed_data}/sig_stage_DEGs_XY.Robj",    # Robj with the filtered XY dynamic genes
        samplesheet=f"{output_tables}/samplesheet.csv"               # Description of the samples
    params:
        output_folder=f"{output_tables}/"      # Location where the result tables are saved
    output:
        pdf=f"{output_pdf}/RNA_common_dynamic_DEGs.pdf",
        png=f"{output_png}/RNA_common_dynamic_DEGs.png"
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/RNA_overlap_dynamic_DEG.R"

# Draw venn diagrams of the comparison of dynamic and sex-specific DEGs for each sex
rule RNA_Plot_common_sex_dynamic_DEGs:
    input:
        sex_DEGs=f"{processed_data}/sig_SexDEGs.Robj",               # Robj with the filtered sex DEGs per stage
        XX_stage_DEGs=f"{processed_data}/sig_stage_DEGs_XX.Robj",    # Robj with the filtered XX dynamic genes
        XY_stage_DEGs=f"{processed_data}/sig_stage_DEGs_XY.Robj",    # Robj with the filtered XY dynamic genes
        samplesheet=f"{output_tables}/samplesheet.csv"               # Description of the samples
    params:
        output_folder=f"{output_tables}/"      # Location where the result tables are saved
    output:
        pdf=f"{output_pdf}/RNA_common_sex_dynamic_DEGs.pdf",
        png=f"{output_png}/RNA_common_sex_dynamic_DEGs.png"
    threads: 12
    resources:
        mem_mb=4000
    script:
        "../scripts/RNA_overlap_sex_dynamic_DEG.R"

