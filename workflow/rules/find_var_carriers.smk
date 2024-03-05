# [TITLE]]

# Target Rule for Completion of Pipeline
rule all:
    input:
        "breast_and_lynch_genes_exons_processed.bed",
        "query_vcf_regions_filtered.vcf.gz.tbi",
        "query_vcf_regions_filtered_normalised.vcf.gz.tbi",
        "query_vcf_regions_filtered_normalised_missingness_filter.vcf.gz.tbi",
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes.vcf.gz.tbi",
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes_annotated.vcf.gz.csi",
        "query_vcf_regions_filtered_normalised_missingness_filter_ptv_clinvar_2.vcf.gz.tbi",
        "query_vcf_ptv_clinvar_2_tally.tsv"
    output:
        "pipeline_complete.txt"
    shell:
        "touch {output}"

# Get Regions
rule biomart_extract:
    output:
        "breast_and_lynch_genes_exons_raw.bed"
    params:
        gene=config["gene"]
    shell:
        "python ../../scripts/biomart_gene_query.py --gene {params.gene} --output {output}"

rule bedtools_process:
    input:
        "breast_and_lynch_genes_exons_raw.bed"
    output:
        "breast_and_lynch_genes_exons_processed.bed"
    params:
        chromosome_sizes=config["chromosome_sizes"]
    shell:
        "bash ../../scripts/bedtools_process_exons.sh {input} {params.chromosome_sizes} {output}"

rule filter_pvcf_to_regions:
    input:
        "breast_and_lynch_genes_exons_processed.bed"
    output:
        "query_vcf_regions_filtered.vcf.gz",
        "query_vcf_regions_filtered.vcf.gz.tbi",
    params:
        query_vcf=config["query_vcf"]
    shell:
        """
        sbatch --output="test_slurm/slurm-%x-%A_%a.out" --time=5-00:00:00 --cpus-per-task=2 \
        ../../scripts/filter_regions.sh {params.query_vcf} {input} {output[0]}
        """

rule normalisation:
    input:
        "query_vcf_regions_filtered.vcf.gz",
    output:
        "query_vcf_regions_filtered_normalised.vcf.gz",
        "query_vcf_regions_filtered_normalised.vcf.gz.tbi",
    params:
        singularity_sif=config["singularity_sif"],
        reference_genome_path=config["reference_genome_path"]
    shell:
        """
       sbatch --output="test_slurm/slurm-%x-%A_%a.out" --time=5-00:00:00 --cpus-per-task=2 \
       ../../scripts/variant_normalisation.sh {input} {params.singularity_sif} {params.reference_genome_path} {output[0]}
       """

rule missingness_filter:
    input:
        "query_vcf_regions_filtered_normalised.vcf.gz",
    output:
        "query_vcf_regions_filtered_normalised_missingness_filter.vcf.gz",
        "query_vcf_regions_filtered_normalised_missingness_filter.vcf.gz.tbi"
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
        "query_vcf_regions_filtered_normalised_missingness_filter.vcf.gz",
    output:
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes.vcf.gz",
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes.vcf.gz.tbi",
    shell:
        """
       sbatch --output="test_slurm/slurm-%x-%A_%a.out" --time=5-00:00:00 \
       ../../scripts/drop_genotypes.sh {input} {output[0]}
       """

# Annotation
rule annotation:
    input:
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes.vcf.gz",
    output:
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes_annotated.vcf.gz",
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes_annotated.vcf.gz.csi",
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
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes_annotated.vcf.gz",
    output:
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes_annotated.tsv"
    params:
        sbatch_job_name="--job-name=remove_header",
        sbatch_params=config["sbatch_low_cpu"],
    shell:
        """
        sbatch {params.sbatch_job_name} {params.sbatch_params} ../../scripts/remove_header.sh {input} {output}
        """

rule identify_variants:
    input:
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes_annotated.tsv"
    output:
        "query_vcf_variant_lists_ptv_clinvar_two.tsv",
    params:
        gene=config["gene"]
    shell:
        """
        Rscript ../../scripts/generate_variant_lists.R {input} {params.gene} {output[0]}
        """

rule filter_variants:
    input:
        "query_vcf_regions_filtered_normalised_missingness_filter.vcf.gz",
        "query_vcf_variant_lists_ptv_clinvar_two.tsv",
    output:
        "query_vcf_regions_filtered_normalised_missingness_filter_ptv_clinvar_2.vcf.gz",
        "query_vcf_regions_filtered_normalised_missingness_filter_ptv_clinvar_2.vcf.gz.tbi"
    params:
        gene = config["gene"],
        sbatch_job_name="--job-name=filter_variants",
        sbatch_params=config["sbatch_low_cpu"],
    shell:
        """
        sbatch {params.sbatch_job_name} {params.sbatch_params} ../../scripts/filter_variants.sh \
        {input[0]} {input[1]} {output[0]}
        """

rule tally_variants:
    input:
        "query_vcf_regions_filtered_normalised_missingness_filter_ptv_clinvar_2.vcf.gz",
    output:
        "query_vcf_ptv_clinvar_2_tally.tsv"
    resources:
        mem_mb=8000
    shell:
        """
        sbatch -J tallies \
            -o test_slurm/tallies.out \
            --wrap='python ../../scripts/tally_variants.py {input} {output}_int && mv {output}_int {output}'
        """
