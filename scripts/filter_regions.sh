#!/bin/bash

INPUT="$1"
BED_FILE="$2"
OUTPUT="$3"

module load BCFtools

# Remove Sites with Missingness < Specified Amount
bcftools view -R "$BED_FILE" "$INPUT" --threads 8 -Oz -o "$OUTPUT"
tabix -f -p vcf "$OUTPUT"
