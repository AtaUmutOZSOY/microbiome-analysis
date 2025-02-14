#!/bin/bash

# Error handling
set -e
set -o pipefail

# Help function
usage() {
    echo "Usage: $0 
    -i INPUT_DIR [QIIME2 output directory]
    -o OUTPUT_DIR [directory for desktop-friendly outputs]"
    exit 1
}

# Parse command line arguments
while getopts "i:o:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        h) usage;;
        ?) usage;;
    esac
done

# Check if required arguments are present
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Required arguments missing"
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/taxonomy"
mkdir -p "$OUTPUT_DIR/diversity"
mkdir -p "$OUTPUT_DIR/abundance"
mkdir -p "$OUTPUT_DIR/feature_table"

# Function to activate QIIME2 environment
activate_qiime2() {
    eval "$(conda shell.bash hook)"
    conda activate qiime2-2023.5
}

# Activate QIIME2
activate_qiime2

echo "Converting QIIME2 artifacts to desktop-friendly formats..."

# Export feature table
if [ -f "$INPUT_DIR/table.qza" ]; then
    echo "Exporting feature table..."
    qiime tools export \
        --input-path "$INPUT_DIR/table.qza" \
        --output-path "$OUTPUT_DIR/feature_table"
    
    # Convert biom to TSV
    biom convert \
        -i "$OUTPUT_DIR/feature_table/feature-table.biom" \
        -o "$OUTPUT_DIR/feature_table/feature-table.tsv" \
        --to-tsv
fi

# Export representative sequences
if [ -f "$INPUT_DIR/rep-seqs.qza" ]; then
    echo "Exporting representative sequences..."
    qiime tools export \
        --input-path "$INPUT_DIR/rep-seqs.qza" \
        --output-path "$OUTPUT_DIR/feature_table"
fi

# Export taxonomy
if [ -f "$INPUT_DIR/taxonomy.qza" ]; then
    echo "Exporting taxonomy..."
    qiime tools export \
        --input-path "$INPUT_DIR/taxonomy.qza" \
        --output-path "$OUTPUT_DIR/taxonomy"
fi

# Export phylogenetic tree
if [ -f "$INPUT_DIR/rooted-tree.qza" ]; then
    echo "Exporting phylogenetic tree..."
    qiime tools export \
        --input-path "$INPUT_DIR/rooted-tree.qza" \
        --output-path "$OUTPUT_DIR"
fi

# Export diversity metrics
for metric in observed_features shannon evenness faith_pd; do
    if [ -f "$INPUT_DIR/diversity/${metric}_vector.qza" ]; then
        echo "Exporting ${metric} diversity..."
        qiime tools export \
            --input-path "$INPUT_DIR/diversity/${metric}_vector.qza" \
            --output-path "$OUTPUT_DIR/diversity/${metric}"
    fi
done

# Export beta diversity distance matrices
for metric in unweighted_unifrac weighted_unifrac bray_curtis jaccard; do
    if [ -f "$INPUT_DIR/diversity/core-metrics-results/${metric}_distance_matrix.qza" ]; then
        echo "Exporting ${metric} distance matrix..."
        qiime tools export \
            --input-path "$INPUT_DIR/diversity/core-metrics-results/${metric}_distance_matrix.qza" \
            --output-path "$OUTPUT_DIR/diversity/${metric}"
    fi
done

# Export taxonomic level tables
for level in {2..7}; do
    if [ -f "$INPUT_DIR/taxonomy/table-l${level}.qza" ]; then
        echo "Exporting taxonomic level ${level} table..."
        qiime tools export \
            --input-path "$INPUT_DIR/taxonomy/table-l${level}.qza" \
            --output-path "$OUTPUT_DIR/taxonomy/level${level}"
        
        # Convert biom to TSV
        biom convert \
            -i "$OUTPUT_DIR/taxonomy/level${level}/feature-table.biom" \
            -o "$OUTPUT_DIR/taxonomy/level${level}/feature-table.tsv" \
            --to-tsv
    fi
done

# Export relative frequency tables
for level in {2..7}; do
    if [ -f "$INPUT_DIR/taxonomy/relative-table-l${level}.qza" ]; then
        echo "Exporting relative frequency table for level ${level}..."
        qiime tools export \
            --input-path "$INPUT_DIR/taxonomy/relative-table-l${level}.qza" \
            --output-path "$OUTPUT_DIR/taxonomy/level${level}_relative"
        
        # Convert biom to TSV
        biom convert \
            -i "$OUTPUT_DIR/taxonomy/level${level}_relative/feature-table.biom" \
            -o "$OUTPUT_DIR/taxonomy/level${level}_relative/feature-table.tsv" \
            --to-tsv
    fi
done

echo "Conversion complete! Desktop-friendly files are available in $OUTPUT_DIR" 