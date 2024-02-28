# [TITLE]]

# Target Rule for Completion of Pipeline
rule all:
    input:
        "{gene}_exons_raw.bed".format(gene=config["gene"])
    output:
        "pipeline_complete.txt"
    shell:
        "touch {output}"

# Get Regions
rule rule:
    input:
        input.txt	
    output:
        "{gene}_exons_raw.bed".format(gene=config["gene"])
    params:
        gene=config["gene"]
    shell:
        "python ../../scripts/biomart_gene_query.py --gene {params.gene} --output {output}"
