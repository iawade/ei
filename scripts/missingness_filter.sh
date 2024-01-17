#!/bin/bash

INPUT="$1"
MISSINGNESS_CUTOFF="$2"
OUTPUT="$3"

module load BCFtools

# Remove Sites with Missingness < Specified Amount
bcftools filter "$INPUT" --threads 8 -Oz -o "$OUTPUT" --exclude "F_MISSING > $MISSINGNESS_CUTOFF"
tabix -f -p vcf "$OUTPUT"
