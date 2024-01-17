#!/bin/bash

# Variables
INPUT="$1" 
GENOME="$2"
OUTPUT="$3"

# Sort, merge and slop 
# (simplest representation of exonic regions of gene from bed produced by biomart)
sort -k1,1 -k2,2n $INPUT > "${INPUT//.bed/_sorted.bed}"
bedtools merge -i "${INPUT//.bed/_sorted.bed}" > "${INPUT//.bed/_sorted_merged.bed}"
bedtools slop -i "${INPUT//.bed/_sorted_merged.bed}" -g "$GENOME" -b 20  > "$OUTPUT"

# Clean up
rm "${INPUT//.bed/_sorted.bed}" "${INPUT//.bed/_sorted_merged.bed}"
