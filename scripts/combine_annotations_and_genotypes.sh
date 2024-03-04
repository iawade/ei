#!/bin/bash

module load BCFtools

# Variables
no_genotypes_vcf="$1"
with_genotypes_anno_vcf="$2"
output="$3"
uncompressed_output=${output//.vcf.gz/.vcf}

# Extract header lines from the header file
# zgrep "^##" "$no_genotypes_vcf" | grep -v "drop-genotypes" > "$uncompressed_output" #bcftool commands in header for posterity
zgrep "^#CHROM" "$with_genotypes_anno_vcf" > "$uncompressed_output"

# Combine columns 1-8 from columns_1_8_file with columns 9-end from columns_9_end_file
paste <(zcat "$no_genotypes_vcf" | awk '{OFS="\t"; if($0 !~ /^#/) print $1,$2,$3,$4,$5,$6,$7,$8}') \
      <(zcat "$with_genotypes_anno_vcf" | awk '{OFS="\t"; if($0 !~ /^#/) {for(i=9; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}}') \
      >> "$uncompressed_output"

# Compress and Index Output
bgzip "$uncompressed_output"
tabix -p vcf "$output"
