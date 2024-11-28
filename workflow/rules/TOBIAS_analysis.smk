'''
Author: Isabelle Stévant
Affiliation: Mammalian sex determination lab, University of Bar Ilan
Date: 01/10/2024
Licence: MIT


Snakemake rules concerning the ATAC footprint analysis using TOBIAS

The pipeline support mm10 and mm39 mouse genome versions
The parameters of the analsysis are defined in the analysis_parameters.yaml configuration file.

Pipeline created to analyse the paired time series RNA and ATAC-seq of the supporting cells from mouse fetal gonads.
Stévant et al. 2024
doi:

'''

configfile: "analysis_parameters.yaml"

sexes = config["sexes"]
stages = config["stages"]
transgenes = ["XX_Enh8-mCherry", "XY_Sox9-IRES-GFP"]
TFBS_background = ["genome", "conditions"]
TFBS_to_plot = ["_ARID3BEMX2HOXsISXLHX9MSX1", "_DMRT1SOX4568910111315"]

# Get file path according to the genome version
input_data = f'{config["path_to_data"]}{config["genome_version"]}'
processed_data = f'{config["path_to_process"]}{config["genome_version"]}/TOBIAS'
output_png = f'{config["path_to_graphs"]}{config["genome_version"]}/TOBIAS'
output_pdf = f'{config["path_to_graphs"]}{config["genome_version"]}/TOBIAS/PDF'
output_tables = f'{config["path_to_tables"]}{config["genome_version"]}/TOBIAS'
RNA_tables = f'{config["path_to_tables"]}{config["genome_version"]}/RNA'
ATAC_tables = f'{config["path_to_tables"]}{config["genome_version"]}/ATAC'
ATAC_processed_data = f'{config["path_to_process"]}{config["genome_version"]}/ATAC'


if config["genome_version"] == "mm10":
    genome = f'{input_data}/gencode.vM25.annotation.gtf.gz'
    # ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm10"]
    # ATAC_norm_bigwig_folder = config["ATAC_norm_bigwig_folder_mm10"]
    fasta_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz"
    fasta_gz = f"{processed_data}/GRCm38.p6.genome.fa.gz"
    fasta = f"{processed_data}/GRCm38.p6.genome.fa"
else :
    genome = f'{input_data}/gencode.vM34.annotation.gtf.gz'
    # ATAC_bigwig_folder = config["ATAC_bigwig_folder_mm39"]
    # ATAC_norm_bigwig_folder = config["ATAC_norm_bigwig_folder_mm39"]
    fasta_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.genome.fa.gz"
    fasta_gz = f"{processed_data}/GRCm39.genome.fa.gz"
    fasta = f"{processed_data}/GRCm39.genome.fa"

# List of output figures
rule_TOBIAS_input_list = [
    f"{output_png}/TOBIAS_sex_DAR_bindiff.png",
    # "test.pdf"
    expand(f"{output_pdf}/{{TFs}}_footprint.pdf", TFs=TFBS_to_plot)
]

## Uncomment if you want to run only this pipeline
# rule all:
#   input:
#       rule_TOBIAS_input_list


###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

rule TOBIAS_download_genome:
    output:
        fa = f"{fasta}"
    shell:
        "curl {fasta_url} --output {fasta_gz} && gzip -d  {fasta_gz}"

rule TOBIAS_download_merged_bam:
    output: f"{processed_data}/{{stage}}_{{transgene}}_merged_sorted.bam"
    shell:
        "cp /home/istevant/work/data/ATACseq/merged_bam/*_merged_sorted.bam {processed_data}"

rule TOBIAS_ATACorrect:
    input:
        bam=f"{processed_data}/{{stage}}_{{transgene}}_merged_sorted.bam",
        genome=f"{fasta}",
        peaks=f'{input_data}/ATAC_all_OCR.bed'
    output:
        output_bw=f"{processed_data}/footprint/{{stage}}_{{transgene}}_merged_sorted_corrected.bw"
    threads: 12
    resources:
        mem_mb=64000
    shell:
        "TOBIAS ATACorrect \
        --bam {input.bam} \
        --genome {input.genome} \
        --peaks {input.peaks} \
        --cores 12 \
        --outdir {processed_data}/footprint"


rule TOBIAS_ScoreBigwig:
    input:
        signal = f"{processed_data}/footprint/{{stage}}_{{transgene}}_merged_sorted_corrected.bw",
        bam=f"{processed_data}/{{stage}}_{{transgene}}_merged_sorted.bam",
        genome=f"{fasta}",
        # peaks=f'{input_data}/ATAC_all_OCR.bed'
        peaks=f"{ATAC_tables}/ATAC_sig_SexDARs.bed"
    output:
        output_bw=f"{processed_data}/footprint/{{stage}}_{{transgene}}_footprints.bw"
    threads: 12
    resources:
        mem_mb=64000
    shell:
        "TOBIAS ScoreBigwig \
        --signal {input.signal} \
        --regions {input.peaks} \
        --cores 12 \
        --output {output.output_bw}"

rule TOBIAS_get_expressed_TF_matrices:
    input:
        TPM=f"{RNA_tables}/TPM.csv",
        TF_genes=config["TF_genes"]
    output:
        matrices=f"{output_tables}/Merged_exp_TF_matrices_pfms.txt"
    params:
        minTPM=config["RNA_minTPM"]
    threads: 1
    resources:
        mem_mb=64000
    script:
        "../scripts/TOBIAS_get_expressed_TF.R"

rule TOBIAS_BINDetect:
    input:
        motifs = f"{ATAC_processed_data}/sex_merged_motifs_pfms.txt",
        signal = expand(f"{processed_data}/footprint/{{stage}}_{{transgene}}_footprints.bw", stage=stages, transgene=transgenes),
        genome=f"{fasta}",
        # peaks=f'{input_data}/ATAC_all_OCR.bed'
        peaks=f"{ATAC_tables}/ATAC_sig_SexDARs.bed"
    output:
        output_file=f"{processed_data}/BINDdetect/bindetect_results.txt",
    params:
        output_dir=f"{processed_data}/BINDdetect",
        conditions=expand(f"{{stage}}_{{transgene}}", stage=stages, transgene=transgenes)
    threads: 12
    resources:
        mem_mb=64000
    shell:
        "TOBIAS BINDetect \
        --motifs {input.motifs}\
        --signal {input.signal}\
        --genome {input.genome}\
        --peaks {input.peaks}\
        --cores 12 \
        --cond-names {params.conditions}\
        --outdir {params.output_dir}"


###########################################
#                                         #
#                 Plots                   #
#                                         #
###########################################

rule TOBIAS_plot_sex_DAR_BINDetect:
    input:
        bindetect=f"{processed_data}/BINDdetect/bindetect_results.txt",
        heatmap_matrice=f"{ATAC_processed_data}/merged_enr_matrix.RData"
    output:
        pdf=f"{output_pdf}/TOBIAS_sex_DAR_bindiff.pdf",
        png=f"{output_png}/TOBIAS_sex_DAR_bindiff.png"
    threads: 2
    resources:
        cpus_per_task=2,
        mem_mb=12000
    script:
        "../scripts/TOBIAS_sex_DAR_bindiff.R"

rule TOBIAS_PlotAggregate_XX:
    input:
        TFBS = f"{processed_data}/BINDdetect/_ARID3BEMX2HOXsISXLHX9MSX1/beds/",
        signal = expand(f"{processed_data}/footprint/", transgene=transgenes),
    output:
        output_file=f"{output_pdf}/_ARID3BEMX2HOXsISXLHX9MSX1_footprint.pdf"
        # output_file="test.pdf"
    params:
        TFs = f"{processed_data}/BINDdetect/_ARID3BEMX2HOXsISXLHX9MSX1/beds/_ARID3BEMX2HOXsISXLHX9MSX1_all.bed",
        signal = expand(f"{processed_data}/footprint/E13.5*_merged_sorted_corrected.bw", transgene=transgenes),
    threads: 2
    resources:
        mem_mb=12000
    shell:
        "TOBIAS PlotAggregate \
        --TFBS {params.TFs} \
        --signals {params.signal} \
        --output {output.output_file} \
        --flank 75 \
        --log-transform \
        --signal-on-x \
        --plot_boundaries"

rule TOBIAS_PlotAggregate_XY:
    input:
        TFBS = f"{processed_data}/BINDdetect/_DMRT1SOX4568910111315/beds/",
        signal = expand(f"{processed_data}/footprint/", transgene=transgenes),
    output:
        output_file=f"{output_pdf}/_DMRT1SOX4568910111315_footprint.pdf"
        # output_file="test.pdf"
    params:
        TFs = f"{processed_data}/BINDdetect/_DMRT1SOX4568910111315/beds/_DMRT1SOX4568910111315_all.bed",
        signal = expand(f"{processed_data}/footprint/E13.5*_merged_sorted_corrected.bw", transgene=transgenes),
    threads: 2
    resources:
        mem_mb=12000
    shell:
        "TOBIAS PlotAggregate \
        --TFBS {params.TFs} \
        --signals {params.signal} \
        --output {output.output_file} \
        --flank 75 \
        --log-transform \
        --signal-on-x \
        --plot_boundaries"