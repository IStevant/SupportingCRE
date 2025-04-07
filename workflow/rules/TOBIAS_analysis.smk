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
MULTI_tables = f'{config["path_to_tables"]}{config["genome_version"]}/MULTI'


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
    f"{processed_data}/BINDdetect/bindetect_results.txt",
    f"{output_png}/TOBIAS_sex_DAR_bindiff.png",
    expand(f"{output_pdf}/{{TFs}}_footprint.pdf", TFs=TFBS_to_plot),
    f"{processed_data}/plot_example_4.log",
    f"{processed_data}/BINDdetect/XX_all_bound_TFs.bed",
    f"{processed_data}/BINDdetect_WT1/XY_bound_WT1.bed"
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
        "cp /home/istevant/work/data/ATACseq/SupportingCRE_merged_bam/*_merged_sorted.bam {processed_data}"

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

# rule TOBIAS_get_expressed_TF_matrices:
#     input:
#         TPM=f"{RNA_tables}/TPM.csv",
#         TF_genes=config["TF_genes"]
#     output:
#         matrices=f"{output_tables}/Merged_exp_TF_matrices_pfms.txt"
#     params:
#         minTPM=config["RNA_minTPM"]
#     threads: 1
#     resources:
#         mem_mb=64000
#     script:
#         "../scripts/TOBIAS_get_expressed_TF.R"

rule TOBIAS_get_WT1_matrix:
    output:
        matrix=f"{output_tables}/WT1_pfm.txt"
    threads: 1
    resources:
        mem_mb=64000
    script:
        "../scripts/TOBIAS_get_WT1_matrix.R"

rule TOBIAS_BINDetect:
    input:
        motifs = f"{ATAC_processed_data}/sex_merged_motifs_pfms.txt",
        signal = expand(f"{processed_data}/footprint/{{stage}}_{{transgene}}_footprints.bw", stage=stages, transgene=transgenes),
        genome = f"{fasta}",
        peaks = f"{ATAC_tables}/ATAC_sig_SexDARs.bed"
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



rule TOBIAS_BINDetect_WT1:
    input:
        motifs = f"{output_tables}/WT1_pfm.txt",
        signal = expand(f"{processed_data}/footprint/{{stage}}_{{transgene}}_footprints.bw", stage=stages, transgene=transgenes),
        genome = f"{fasta}",
        peaks = f"{ATAC_tables}/ATAC_sig_SexDARs.bed"
    output:
        output_file=f"{processed_data}/BINDdetect_WT1/bindetect_results.txt",
    params:
        output_dir=f"{processed_data}/BINDdetect_WT1",
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


rule TOBIAS_merge_bed:
    input:
        bindetect=f"{processed_data}/BINDdetect/bindetect_results.txt"
    params:
        bed_folders = f"{processed_data}/BINDdetect/_*",
        wt1_bed_folders = f"{processed_data}/BINDdetect_WT1/_*"
    output:
        XX_bound = f"{processed_data}/BINDdetect/XX_bound_TFs.bed",
        XY_bound = f"{processed_data}/BINDdetect/XY_bound_TFs.bed"
    threads: 1
    resources:
        mem_mb=12000
    shell:
        "cat {params.bed_folders}/beds/*XX*_bound.bed {params.wt1_bed_folders}/beds/*XX*_bound.bed | grep -E 'EMX2|GATA|FOX|RUNX1|NR5A1|DMRT1|Wt1' | awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$4,$11,$6}}' | sort | uniq | sed 's/_Wt1/WT1/g' | sed 's/_.*EMX2.*\t0\t/EMX2-LHX9\t0\t/g' | sed 's/_.*GATA.*\t0\t/GATAs\t0\t/g' | sed 's/_FOXs.*\t0\t/FOXs\t0\t/g' | sed 's/_.*RUNX.*\t0\t/RUNX1\t0\t/g' | sed 's/_NR5A1.*\t0\t/NR5A1\t0\t/g' | sed 's/_.*SOX.*\t0\t/DMRT1-SOXs\t0\t/g' > {output.XX_bound} && \
        cat {params.bed_folders}/beds/*XY*_bound.bed {params.wt1_bed_folders}/beds/*XY*_bound.bed | grep -E 'EMX2|GATA|FOX|RUNX1|NR5A1|DMRT1|Wt1' | awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$4,$11,$6}}' | sort | uniq  | sed 's/_Wt1/WT1/g' | sed 's/_.*EMX2.*\t0\t/EMX2-LHX9\t0\t/g' | sed 's/_.*GATA.*\t0\t/GATAs\t0\t/g' | sed 's/_FOXs.*\t0\t/FOXs\t0\t/g' | sed 's/_.*RUNX.*\t0\t/RUNX1\t0\t/g' | sed 's/_NR5A1.*\t0\t/NR5A1\t0\t/g' | sed 's/_.*SOX.*\t0\t/DMRT1-SOXs\t0\t/g' > {output.XY_bound}"


rule TOBIAS_merge_all_TF_bed:
    input:
        bindetect=f"{processed_data}/BINDdetect/bindetect_results.txt"
    params:
        bed_folders = f"{processed_data}/BINDdetect/_*"
    output:
        XX_bound = f"{processed_data}/BINDdetect/XX_all_bound_TFs.bed",
        XY_bound = f"{processed_data}/BINDdetect/XY_all_bound_TFs.bed"
    threads: 2
    resources:
        mem_mb=12000
    shell:
        "cat {params.bed_folders}/beds/*XX*_bound.bed | awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$4,$6}}' | sort | uniq > {output.XX_bound} && \
        cat {params.bed_folders}/beds/*XY*_bound.bed | awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$4,$6}}' | sort | uniq > {output.XY_bound}"


rule TOBIAS_merge_bed_WT1:
    input:
        bindetect=f"{processed_data}/BINDdetect_WT1/bindetect_results.txt"
    params:
        bed_folders = f"{processed_data}/BINDdetect_WT1/_*"
    output:
        XX_bound = f"{processed_data}/BINDdetect_WT1/XX_bound_WT1.bed",
        XY_bound = f"{processed_data}/BINDdetect_WT1/XY_bound_WT1.bed"
    threads: 2
    resources:
        mem_mb=12000
    shell:
        "cat {params.bed_folders}/beds/*XX*_bound.bed | awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$4,$6}}' | sort | uniq > {output.XX_bound} && \
        cat {params.bed_folders}/beds/*XY*_bound.bed | awk -v FS='\t' -v OFS='\t' '{{print $1,$2,$3,$4,$6}}' | sort | uniq > {output.XY_bound}"


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
        bindDetect=f"{output_png}/TOBIAS_sex_DAR_bindiff.png"
    output:
        output_file=f"{output_pdf}/_ARID3BEMX2HOXsISXLHX9MSX1_footprint.pdf"
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
        # TFBS = f"{processed_data}/BINDdetect/_DMRT1SOX4568910111315/beds/",
        # signal = expand(f"{processed_data}/footprint/", transgene=transgenes),
        bindDetect=f"{output_png}/TOBIAS_sex_DAR_bindiff.png"
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

rule TOBIAS_Plot_footprints_examples:
    input:
        genome=f"{genome}",
        gene_bed=f"{input_data}/gene_standard.bed",
        peaks=f"{ATAC_tables}/ATAC_norm_counts.csv",
        linkage=f"{MULTI_tables}/all_sig_gene2peak_linkage.csv",
        peak_list=f"{input_data}/gTrack_footprint_examples.tsv",
        TPM=f"{RNA_tables}/TPM.csv",
        XX_bound = f"{processed_data}/BINDdetect/XX_bound_TFs.bed",
        XY_bound = f"{processed_data}/BINDdetect/XY_bound_TFs.bed"
    params:
        bw_folder=f"{ATAC_norm_bigwig_folder}",
        ChIP_folder=f"{input_data}/ChIP/",
        save_folder=f"{output_pdf}"
    output: 
        log= f"{processed_data}/plot_example_4.log"
    threads: 12
    resources:
        mem_mb=64000
    script:
        "../scripts/TOBIAS_plot_genomic_tracks.R"
