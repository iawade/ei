#!/bin/bash

# submit as
# sbatch <script> <input_vcf> <target_intervals_tsv>

# Load packages
module load BCFtools

# Set Variables
INPUT_VCF="$1"
TARGET_INTERVALS_TSV="$2"
OUTPUT="$3"

# Extract comment lines from VCF file
zgrep "^#" "$INPUT_VCF" > "comment_lines_${TARGET_INTERVALS_TSV}.txt"

# Match non-comment lines using the first two fields of TSV file
awk -F'\t' 'NR==FNR{a[$1,$2,$4,$5]; next} ($1,$2,$4,$5) in a' "$TARGET_INTERVALS_TSV" \
  <(zgrep -v "^#" "$INPUT_VCF") > "${TARGET_INTERVALS_TSV}_matched_lines.txt"


# Combine comment lines and matched lines
cat "comment_lines_${TARGET_INTERVALS_TSV}.txt" "${TARGET_INTERVALS_TSV}_matched_lines.txt" > "${TARGET_INTERVALS_TSV}_combined.vcf"

bgzip -f "${TARGET_INTERVALS_TSV}_combined.vcf"
tabix -f -p vcf "${TARGET_INTERVALS_TSV}_combined.vcf.gz"

mv "${TARGET_INTERVALS_TSV}_combined.vcf.gz" "$OUTPUT"
mv "${TARGET_INTERVALS_TSV}"_combined.vcf.gz.tbi "$OUTPUT".tbi

# shellcheck disable=SC2035
rm *_matched_lines* comment_lines_*
