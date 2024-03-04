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
        "query_vcf_regions_filtered_normalised_missingness_filter_drop_genotypes_annotated.tsv",
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

# rule identify_variants:
#     input:
#         "{gene}_ukb_variants_participants_filtered_normalised_drop_genotypes_annotated.tsv"
#     output:
#         "{gene}_variant_lists_ptv_only.tsv",
#         "{gene}_variant_lists_ptv_clinvar_one.tsv",
#         "{gene}_variant_lists_ptv_clinvar_two.tsv",
#         "{gene}_variant_lists_rare.tsv",
#         "{gene}_variant_lists_ptv_clinvar_one_rare.tsv",
#         "{gene}_variant_lists_ptv_clinvar_two_rare.tsv",
#     params:
#         gene=config["gene"],
#         clinvar_phenotypes=config["clinvar_phenotypes"],
#         maf_cutoff=config["maf_cutoff"]
#     shell:
#         """
#         Rscript ../../scripts/generate_variant_lists.R {input} {params.gene} {params.clinvar_phenotypes} {params.maf_cutoff} \
#         {output[0]} {output[1]} {output[2]} {output[3]} {output[4]} {output[5]}
#         """
#
# rule filter_variants:
#     input:
#         "{gene}_ukb_variants_participants_filtered_normalised_missingness_filter.vcf.gz",
#         "{gene}_variant_lists_{variant_tranche}.tsv"
#     output:
#         "{gene}_ukb_variants_participants_filtered_normalised_{variant_tranche}.vcf.gz",
#         "{gene}_ukb_variants_participants_filtered_normalised_{variant_tranche}.vcf.gz.tbi"
#     params:
#         gene = config["gene"],
#         sbatch_job_name="--job-name=filter_variants",
#         sbatch_params=config["sbatch_low_cpu"],
#     shell:
#         """
#         sbatch {params.sbatch_job_name} {params.sbatch_params} ../../scripts/filter_variants.sh \
#         {input[0]} {input[1]} {output}
#         """
#
# rule tally_variants:
#     input:
#         "{gene}_ukb_variants_participants_filtered_normalised_{variant_tranche}.vcf.gz"
#     output:
#         "{gene}_ukb_{variant_tranche}_tally.tsv"
#     params:
#         gene = config["gene"]
#     resources:
#         mem_mb=8000
#     shell:
#         """
#         sbatch -J {wildcards.gene}_{wildcards.variant_tranche}_tallies \
#             -o test_slurm/{wildcards.gene}_{wildcards.variant_tranche}_tallies.out \
#             --wrap='python ../../scripts/tally_variants.py {input} {output}_int && mv {output}_int {output}'
#         """
