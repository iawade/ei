# [TITLE]]

from typing import List

# with open(config["datasets"]) as f:
#     datasets = [line.strip() for line in f]
#
# config["datasets"] = datasets
#
# with open(config["dataset_name_files"]) as f1:
#     dataset_name_files = [line.strip() for line in f1]

# Target Rule for Completion of Pipeline
rule all:
    input:
        "pipeline_complete.txt"

# Get Regions
rule biomart_extract:
    output:
        "{gene}_exons_raw.bed".format(gene=config["gene"])
    params:
        gene=config["gene"]
    shell:
        "python ../../scripts/biomart_gene_query.py --gene {params.gene} --output {output}"

rule bedtools_process:
    input:
        "{gene}_exons_raw.bed"
    output:
        "{gene}_processed_exons.bed"
    params:
        gene=config["gene"],
        chromosome_sizes=config["chromosome_sizes"]
    shell:
        "bash ../../scripts/bedtools_process_exons.sh {input} {params.chromosome_sizes} {output}"

rule identify_pvcf_chunk:
    input:
        "{gene}_processed_exons.bed"
    output:
        "{gene}_pvcf_chunk_name.txt"
    params:
        gene=config["gene"],
        pvcf_blocks=config["pvcf_blocks"]
    shell:
        """
        python ../../scripts/UKB_exome_file_cross-ref.py {input} {output} {params.pvcf_blocks}
        """

rule extract_variants_and_download:
    input:
        "{gene}_processed_exons.bed",
        "{gene}_pvcf_chunk_name.txt"
    output:
        "{gene}_ukb_variants.vcf.gz",
        "{gene}_ukb_variants.vcf.gz.tbi"
    params:
        gene=config["gene"],
        ukb_bed_path=config["ukb_bed_path"]
    shell:
        """
        sbatch --output="test_slurm/slurm-%x-%A_%a.out" --time=5-00:00:00 ../../scripts/extract_regions_from_ukb.sh {input[0]} {params.ukb_bed_path} {input[1]} {output[0]}
        """

# Catch all due to snakemake quirks re wildcards in target rules and input/output
rule final_rule:
    input:
        "{gene}_ukb_variants.vcf.gz".format(gene=config["gene"]),
        "{gene}_ukb_variants.vcf.gz.tbi".format(gene=config["gene"])
    output:
        "pipeline_complete.txt"
    shell:
        "touch {output}"
