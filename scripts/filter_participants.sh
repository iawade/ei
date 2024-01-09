#!/bin/bash

. "/opt/software/applications/anaconda/3/etc/profile.d/conda.sh" && conda activate /data/scratch/DGE/DUDGE/PREDIGEN/iwade/SMARCA4_related/workflow/envs/SMARCA4_related

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
