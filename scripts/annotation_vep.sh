#!/bin/bash

# Load Variables
INPUT="$1"
OUTPUT="$2"
VEP_SIF="$3"
ASSEMBLY="$4"
FASTA_REFERENCE="$5"
DIR_CACHE="$6"
CACHE_VERSION="$7"
DIR_PLUGINS="$8"
REVEL="$9"
CLINVAR=${10}
LOFTEE=${11}
CADD=${12}

module load BCFtools

singularity exec  -B/data:/data "$VEP_SIF"  vep \
    --assembly "$ASSEMBLY" \
    --offline \
    -i "$INPUT" \
    -o "${OUTPUT//.vcf.gz/.vcf}" \
    --symbol \
    --gene \
    --format vcf \
    --no_stats \
    --fasta "$FASTA_REFERENCE" \
    --cache \
    --dir_cache "$DIR_CACHE" \
    --cache_version "$CACHE_VERSION" \
    --hgvs \
    --af_gnomad \
    --mane \
    --vcf \
    --pick_allele \
    --fork 1 \
    --buffer_size 500 \
    --dir_plugins "$DIR_PLUGINS" \
    --plugin REVEL,"$REVEL" \
    --custom "$CLINVAR",ClinVar,vcf,exact,0,CLNDN,CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNSIGINCL,CLNVC,CLNVCSO,CLNVI \
    --plugin LoF,"$LOFTEE" \
    --plugin CADD,"$CADD" \
    --force_overwrite \
    --vcf_info_field vep

## LOFTEE runs into errors when used with fork - minor - but for completion run with only one thread
## But some kind of error where container won't put samtools into path unless forking is used
### --fork 1 - solid workaround

# Fix names of vars
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' "${OUTPUT//.vcf.gz/.vcf}" -Oz -o "$OUTPUT"

rm "${OUTPUT//.vcf.gz/.vcf}"

# Create index file
tabix -f -C -p vcf "$OUTPUT"

rm "${OUTPUT//.vcf.gz/.vcf}"_warnings.txt
