#!/bin/bash

input="$1"
bcftools_sif="$2"
reference_genome="$3"
output="$4"

if zgrep -m 1 -q "##INFO=<ID=AC" "$input"; then
  singularity exec  -B/data:/data "$bcftools_sif" bcftools norm -m-any --multi-overlaps 0 -f "$reference_genome" "$input" -Ou --threads 2 \
  | singularity exec  -B/data:/data "$bcftools_sif" bcftools view --threads 2 -i "AC>0" -Oz -o "$output"
else
  singularity exec  -B/data:/data "$bcftools_sif" bcftools norm -m-any --multi-overlaps 0 -f "$reference_genome" "$input" -Ou --threads 2 \
  | singularity exec  -B/data:/data "$bcftools_sif" bcftools view --threads 2 -i "AF>0" -Oz -o "$output"
fi

module load BCFtools
tabix -f -p vcf "$output"
