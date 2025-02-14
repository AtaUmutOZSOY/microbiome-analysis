#!/bin/bash

# Hata yönetimi
set -e
set -o pipefail

# Metadata dosyasının yolu
METADATA="projects/selin_ugrakli/metadata_files/metadata_qiime.tsv"
READS_DIR="projects/selin_ugrakli/reads"

# Sample counter for S[X] in filename
sample_counter=1

# Metadata'dan SampleID ve internal_barcode eşleştirmelerini al
while IFS=$'\t' read -r sampleid barcode internal_barcode rest; do
    # Skip header and empty lines
    if [[ $sampleid == "#SampleID" ]] || [[ -z $sampleid ]]; then
        continue
    fi

    # Look for files matching the internal barcode pattern
    for file in projects/selin_ugrakli/reads/${internal_barcode}*_R1_001.fastq.gz; do
        if [ -f "$file" ]; then
            # Get the R1 file
            r1_file="$file"
            # Construct the R2 file name
            r2_file="${file/_R1_/_R2_}"
            
            # Construct new filenames in Casava format
            new_r1="${file%/*}/${sampleid}_S${sample_counter}_L001_R1_001.fastq.gz"
            new_r2="${file%/*}/${sampleid}_S${sample_counter}_L001_R2_001.fastq.gz"
            
            # Rename the files
            mv "$r1_file" "$new_r1"
            mv "$r2_file" "$new_r2"
            
            # Increment sample counter
            ((sample_counter++))
        else
            echo "Warning: No matching file found for internal barcode $internal_barcode R1"
        fi
    done
done < "$METADATA"

echo "Renaming complete!" 