#!/bin/bash

INPUT_BED="$1"
PVCF_CHUNK="$2"
OUTPUT="$3"

UKB_EXOME_DIR="/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release"

echo "$INPUT_BED" "$UKB_EXOME_DIR/$PVCF_CHUNK" "${OUTPUT}"
echo ""$UKB_EXOME_DIR/$PVCF_CHUNK""

bcftools view -R "$INPUT_BED" "$PVCF_CHUNK" --threads 8 -Oz -o "$OUTPUT"

