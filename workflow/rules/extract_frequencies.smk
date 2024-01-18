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
        "{gene}_ukb_variants.vcf.gz.tbi".format(gene=config["gene"]),
        "{gene}_ukb_variants_participants_filtered.vcf.gz.tbi".format(gene=config["gene"]),
        "{gene}_ukb_variants_participants_filtered_normalised.vcf.gz.tbi".format(gene=config["gene"]),
        "{gene}_ukb_variants_participants_filtered_normalised_missingness_filter.vcf.gz.tbi".format(gene=config["gene"]),
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes.vcf.gz.tbi".format(gene=config["gene"]),
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes_annotated.vcf.gz.csi".format(gene=config["gene"]),
        "{gene}_variant_lists_ptv_only.tsv".format(gene=config["gene"]),
        "{gene}_variant_lists_ptv_clinvar_one.tsv".format(gene=config["gene"]),
        "{gene}_variant_lists_ptv_clinvar_two.tsv".format(gene=config["gene"]),
        "{gene}_variant_lists_rare.tsv".format(gene=config["gene"]),
        expand("{gene}_ukb_variants_participants_filtered_normalised_{variant_tranche}.vcf.gz.tbi", gene=config["gene"], variant_tranche=glob_wildcards("{gene}_variant_lists_{variant_tranche}.tsv").variant_tranche),
        expand("{gene}_ukb_{variant_tranche}_tally.tsv", gene=config["gene"], variant_tranche=glob_wildcards("{gene}_variant_lists_{variant_tranche}.tsv").variant_tranche),
    output:
        "pipeline_complete.txt"
    shell:
        "touch {output}"

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

# TODO : I'm not sure if this script could handle if pvcf_chunk > 1; Check and Fix
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

rule filter_participants:
    input:
        "{gene}_ukb_variants.vcf.gz",
    output:
        "{gene}_ukb_variants_participants_filtered.vcf.gz",
        "{gene}_ukb_variants_participants_filtered.vcf.gz.tbi",
    params:
        gene=config["gene"],
        participants_list=config["participants_list"]
    shell:
        """
        sbatch --output="test_slurm/slurm-%x-%A_%a.out" --time=5-00:00:00 --cpus-per-task=2 \
        ../../scripts/filter_participants.sh {input} {params.participants_list} {output[0]}
        """

rule normalisation:
    input:
        "{gene}_ukb_variants_participants_filtered.vcf.gz",
    output:
        "{gene}_ukb_variants_participants_filtered_normalised.vcf.gz",
        "{gene}_ukb_variants_participants_filtered_normalised.vcf.gz.tbi",
    params:
        gene=config["gene"],
        singularity_sif=config["singularity_sif"],
        reference_genome_path=config["reference_genome_path"]
    shell:
        """
       sbatch --output="test_slurm/slurm-%x-%A_%a.out" --time=5-00:00:00 --cpus-per-task=2 \
       ../../scripts/variant_normalisation.sh {input} {params.singularity_sif} {params.reference_genome_path} {output[0]}
       """

rule missingness_filter:
    input:
        "{gene}_ukb_variants_participants_filtered_normalised.vcf.gz",
    output:
        "{gene}_ukb_variants_participants_filtered_normalised_missingness_filter.vcf.gz",
        "{gene}_ukb_variants_participants_filtered_normalised_missingness_filter.vcf.gz.tbi"
    params:
        missingness_cutoff=config["missingness_cutoff"],
        sbatch_job_name = "--job-name=missingness_filter",
        sbatch_params = config["sbatch_high_cpu_low_mem"],
    shell:
        """
        sbatch {params.sbatch_job_name} {params.sbatch_params} ../../scripts/missingness_filter.sh {input} {params.missingness_cutoff} {output}
        """

rule drop_genotypes:
    input:
        "{gene}_ukb_variants_participants_filtered_normalised_missingness_filter.vcf.gz",
    output:
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes.vcf.gz",
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes.vcf.gz.tbi",
    params:
        gene=config["gene"],
    shell:
        """
       sbatch --output="test_slurm/slurm-%x-%A_%a.out" --time=5-00:00:00 \
       ../../scripts/drop_genotypes.sh {input} {output[0]}
       """

# Annotation
rule annotation:
    input:
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes.vcf.gz",
    output:
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes_annotated.vcf.gz",
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes_annotated.vcf.gz.csi",
    params:
        sbatch_job_name="--job-name=annotation_vep",
        sbatch_params=config["sbatch_low_cpu"],
        assembly=config["assembly"],
        fasta_reference=config["reference_genome_path"],
        vep_cache_version=config["vep_cache_version"],
        revel=config["revel"],
        clinvar=config["clinvar"],
        loftee=config["loftee"],
        cadd=config["cadd"]
    shell:
        """
        sbatch {params.sbatch_job_name} {params.sbatch_params} \
         ../../scripts/annotation_vep.sh {input} {output[0]} "/data/scratch/DGE/DUDGE/PREDIGEN/iwade/WES_pipeline/workflow/envs/ensembl_vep_loftee_v110.sif" \
         {params.assembly} {params.fasta_reference} "/data/scratch/DGE/DUDGE/PREDIGEN/iwade/WES_pipeline/reference_data" {params.vep_cache_version} "/data/scratch/DGE/DUDGE/PREDIGEN/iwade/WES_pipeline/reference_data" {params.revel} \
         {params.clinvar} {params.loftee} {params.cadd}
        """

rule remove_vcf_header:
    input:
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes_annotated.vcf.gz"
    output:
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes_annotated.tsv"
    params:
        gene=config["gene"],
        sbatch_job_name="--job-name=remove_header",
        sbatch_params=config["sbatch_low_cpu"],
    shell:
        """
        sbatch {params.sbatch_job_name} {params.sbatch_params} ../../scripts/remove_header.sh {input} {output}
        """

rule identify_variants:
    input:
        "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes_annotated.tsv"
    output:
        "{gene}_variant_lists_ptv_only.tsv",
        "{gene}_variant_lists_ptv_clinvar_one.tsv",
        "{gene}_variant_lists_ptv_clinvar_two.tsv",
        "{gene}_variant_lists_rare.tsv",
    params:
        gene=config["gene"],
        clinvar_phenotypes=config["clinvar_phenotypes"],
        maf_cutoff=config["maf_cutoff"]
    shell:
        """
        Rscript ../../scripts/generate_variant_lists.R {input} {params.clinvar_phenotypes} {params.maf_cutoff} \
        {output[0]} {output[1]} {output[2]} {output[3]} 
        """

rule filter_variants:
    input:
        "{gene}_ukb_variants_participants_filtered_normalised_missingness_filter.vcf.gz",
        "{gene}_variant_lists_{variant_tranche}.tsv"
    output:
        "{gene}_ukb_variants_participants_filtered_normalised_{variant_tranche}.vcf.gz",
        "{gene}_ukb_variants_participants_filtered_normalised_{variant_tranche}.vcf.gz.tbi"
    params:
        gene = config["gene"],
        sbatch_job_name="--job-name=filter_variants",
        sbatch_params=config["sbatch_low_cpu"],
    shell:
        """
        sbatch {params.sbatch_job_name} {params.sbatch_params} ../../scripts/filter_variants.sh \
        {input[0]} {input[1]} {output}
        """

rule tally_variants:
    input:
        "{gene}_ukb_variants_participants_filtered_normalised_{variant_tranche}.vcf.gz"
    output:
        "{gene}_ukb_{variant_tranche}_tally.tsv"
    params:
        gene = config["gene"]
    resources:
        mem_mb=8000
    shell:
        """
        python ../../scripts/tally_variants.py {input} {output}
        """
