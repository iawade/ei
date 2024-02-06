#!/bin/bash

<<'COMMENT'
This script relies on my access and directory structure
Also dependent on files/scripts already on my UKB project
So not portable to another user
COMMENT

module load BCFtools

# Variables
INPUT_BED="$1"
UKB_BED_PATH="$2"
PVCF_CHUNK=$(cat "$3")
OUTPUT_FILENAME="$4"

export project=$(dx pwd)
UKB_EXOME_DIR="/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release"

# Upload Bed File if Needed
if dx ls "${project}/${UKB_BED_PATH}/${INPUT_BED}" &>/dev/null; then
    echo "File ${INPUT_BED} already exists. Skipping upload."
else
    # Upload Bed File
    dx upload "$INPUT_BED" --path "${project}/${UKB_BED_PATH}/"
fi

# Loop over PVCF_CHUNK if it's greater than 1
if [[ $(echo "$PVCF_CHUNK" | wc -l) -gt 1 ]]; then
    # Create a directory to store intermediate outputs

    for i in $(echo $PVCF_CHUNK); do
        # Extract Region
        dx run app-swiss-army-knife \
            -iin="${project}/scripts/query_UKB.sh" \
            -iin="${project}/$UKB_BED_PATH/$INPUT_BED" \
            -iin="$UKB_EXOME_DIR/${i}" \
            -iin="$UKB_EXOME_DIR/${i}.tbi" \
            -icmd="sh query_UKB.sh $INPUT_BED ${i} ${OUTPUT_FILENAME}_${i}" \
            --name="UKB_variant_extractor_${i}" \
            --instance-type="mem2_ssd1_v2_x8" \
            --destination="${project}/$UKB_BED_PATH/" \
            -y

        # Wait for the output file to become available
        while ! dx ls "${project}/$UKB_BED_PATH/${OUTPUT_FILENAME}_${i}" &>/dev/null; do
            echo "Waiting for the output file ${OUTPUT_FILENAME}_${i} to become available..."
            sleep 60  # Adjust the sleep duration as needed
        done

        dx download "${project}/$UKB_BED_PATH/${OUTPUT_FILENAME}_${i}"
        tabix -p vcf "${OUTPUT_FILENAME}_${i}"
    done

    # Concatenate the output files into a single file
    # TODO double-check this ; can't ls with quotes
    bcftools concat $(ls ${OUTPUT_FILENAME}_*vcf.gz) -Oz -o "${OUTPUT_FILENAME}"
    tabix -p vcf "${OUTPUT_FILENAME}"

    # Clean up intermediate files


else
    # Extract Region
    dx run app-swiss-army-knife \
        -iin="${project}/scripts/query_UKB.sh" \
        -iin="${project}/$UKB_BED_PATH/$INPUT_BED" \
        -iin="$UKB_EXOME_DIR/${PVCF_CHUNK}" \
        -iin="$UKB_EXOME_DIR/${PVCF_CHUNK}.tbi" \
        -icmd="sh query_UKB.sh $INPUT_BED $PVCF_CHUNK $OUTPUT_FILENAME" \
        --name="UKB_variant_extractor" \
        --instance-type="mem2_ssd1_v2_x8" \
        --destination="${project}/$UKB_BED_PATH/" \
        -y

    # Wait for the output file to become available
    while ! dx ls "${project}/$UKB_BED_PATH/$OUTPUT_FILENAME" &>/dev/null; do
        echo "Waiting for the output file to become available..."
        sleep 60  # Adjust the sleep duration as needed
    done

    # Download Output
    dx download "${project}/$UKB_BED_PATH/$OUTPUT_FILENAME"
    tabix -p vcf "$OUTPUT_FILENAME"
fi
