#!/bin/bash

module load BCFtools

input="$1"
output="$2"

bcftools view "$input" --drop-genotypes -Oz -o "$output"

tabix -f -p vcf "$output"