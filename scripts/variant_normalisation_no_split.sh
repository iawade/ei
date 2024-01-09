#!/bin/bash

. "/opt/software/applications/anaconda/3/etc/profile.d/conda.sh" && conda activate /data/scratch/DGE/DUDGE/PREDIGEN/iwade/SMARCA4_related/workflow/envs/SMARCA4_related

input="$1"
reference_genome="$2"
output="$3"

if zgrep -m 1 -q "##INFO=<ID=AC" "$input"; then
  bcftools norm -f "$reference_genome" "$input" -Ou --threads 8 \
  | bcftools view --threads 8 -i "AC>0" -Oz -o "$output"
else
  bcftools norm -f "$reference_genome" "$input" -Ou --threads 8 \
  | bcftools view --threads 8 -i "AF>0" -Oz -o "$output"
fi
