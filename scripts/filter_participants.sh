#!/bin/bash

module load BCFtools

# Set Variables
INPUT="$1"
DATASET_NAMES="$2"
OUTPUT="$3"

bcftools view \
  -S "$DATASET_NAMES"  --force-samples \
        "$INPUT" \
        -Oz -o "$OUTPUT" \
        --threads 2

tabix -f -p vcf "$OUTPUT"
