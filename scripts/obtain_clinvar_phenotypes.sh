#!/bin/bash

module load BCFtools

GENE="$1"

bcftools query -f '%INFO/CLNDN;%INFO/GENEINFO\n' /data/scratch/DGE/DUDGE/PREDIGEN/iwade/ref/clinvar/clinvar_20230626.vcf.gz | \
       	grep "$GENE" | cut -d ";" -f 1 | tr '|' '\n' | sort | uniq > "$GENE"_clinvar_phenotypes_for_review.txt

