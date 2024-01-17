#!/bin/bash

INPUT="$1"
OUTPUT="$2"

# Remove header lines from compressed VCF
zgrep -v "^##" "$INPUT" | sed 's/^#//' > "$OUTPUT"
