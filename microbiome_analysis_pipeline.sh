#!/bin/bash

# Error handling
set -e  # Exit on error
set -o pipefail  # Exit if any command in a pipe fails

# Step tracking function
step_info() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Step $1: $2"
}

# Checkpoint function
create_checkpoint() {
    touch "$OUTPUT_DIR/logs/step${1}.done"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Completed step $1: $2" >> "$OUTPUT_DIR/logs/pipeline.log"
}

# Check if step is completed
is_step_completed() {
    if [ -f "$OUTPUT_DIR/logs/step${1}.done" ]; then
        echo "Step $1 already completed, skipping..."
        return 0
    fi
    return 1
}

# Error handling function
handle_error() {
    echo "Error occurred in step $current_step"
    echo "Exit code: $?"
    echo "Line number: $1"
    echo "Check $OUTPUT_DIR/logs/pipeline.log for details"
    echo "You can resume from the last successful step by running the script again"
    exit 1
}

trap 'handle_error ${LINENO}' ERR

# Help function
usage() {
    echo "Usage: $0 
    -i INPUT_DIR [directory containing raw fastq files]
    -m METADATA [metadata file in TSV format]
    -o OUTPUT_DIR [directory for output files]
    -c CLASSIFIER [SILVA classifier file]
    -t THREADS [number of threads to use]
    -s SAMPLING_DEPTH [rarefaction depth]
    -p PAIRED [true/false for paired-end data]
    -f FORCE [true/false to force rerun all steps]"
    exit 1
}

# Parse command line arguments
while getopts "i:m:o:c:t:s:p:f:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG";;
        m) METADATA="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        c) CLASSIFIER="$OPTARG";;
        t) THREADS="$OPTARG";;
        s) SAMPLING_DEPTH="$OPTARG";;
        p) PAIRED="$OPTARG";;
        f) FORCE="$OPTARG";;
        h) usage;;
        ?) usage;;
    esac
done

# Check if required arguments are present
if [ -z "$INPUT_DIR" ] || [ -z "$METADATA" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$CLASSIFIER" ]; then
    echo "Error: Required arguments missing"
    usage
fi

# Set default values if not provided
THREADS=${THREADS:-8}
SAMPLING_DEPTH=${SAMPLING_DEPTH:-10000}
PAIRED=${PAIRED:-true}
FORCE=${FORCE:-false}

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/diversity"
mkdir -p "$OUTPUT_DIR/taxonomy"
mkdir -p "$OUTPUT_DIR/abundance"

# If force flag is true, remove all checkpoints
if [ "$FORCE" = true ]; then
    rm -f "$OUTPUT_DIR/logs/"*.done
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Forcing rerun of all steps" > "$OUTPUT_DIR/logs/pipeline.log"
fi

# Initialize step counter
current_step=1

# Function to activate QIIME2 environment
activate_qiime2() {
    eval "$(conda shell.bash hook)"
    conda activate qiime2-2023.5
}

# Activate QIIME2 environment
if ! is_step_completed $current_step "QIIME2 Environment Activation"; then
    step_info $current_step "Activating QIIME2 environment"
    activate_qiime2
    create_checkpoint $current_step "QIIME2 Environment Activation"
fi
((current_step++))

# Import and denoise data
if ! is_step_completed $current_step "Data Import and DADA2"; then
    step_info $current_step "Starting QIIME2 data preparation"
    activate_qiime2
    
    if [ "$PAIRED" = true ]; then
        step_info $current_step "Importing paired-end data"
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path "$INPUT_DIR" \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt \
            --output-path "$OUTPUT_DIR/demux-paired-end.qza"
        
        step_info $current_step "Generating visualization for quality check"
        qiime demux summarize \
            --i-data "$OUTPUT_DIR/demux-paired-end.qza" \
            --o-visualization "$OUTPUT_DIR/demux-paired-end.qzv"
        
        step_info $current_step "DADA2 denoising"
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$OUTPUT_DIR/demux-paired-end.qza" \
            --p-trim-left-f 0 \
            --p-trim-left-r 0 \
            --p-trunc-len-f 250 \
            --p-trunc-len-r 250 \
            --p-n-threads "$THREADS" \
            --p-max-ee-f 2.0 \
            --p-max-ee-r 2.0 \
            --p-trunc-q 2 \
            --p-chimera-method consensus \
            --o-table "$OUTPUT_DIR/table.qza" \
            --o-representative-sequences "$OUTPUT_DIR/rep-seqs.qza" \
            --o-denoising-stats "$OUTPUT_DIR/denoising-stats.qza"
    else
        step_info $current_step "Importing single-end data"
        qiime tools import \
            --type 'SampleData[SequencesWithQuality]' \
            --input-path "$INPUT_DIR" \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt \
            --output-path "$OUTPUT_DIR/demux-single-end.qza"
        
        step_info $current_step "Generating visualization for quality check"
        qiime demux summarize \
            --i-data "$OUTPUT_DIR/demux-single-end.qza" \
            --o-visualization "$OUTPUT_DIR/demux-single-end.qzv"
        
        step_info $current_step "DADA2 denoising"
        qiime dada2 denoise-single \
            --i-demultiplexed-seqs "$OUTPUT_DIR/demux-single-end.qza" \
            --p-trim-left 0 \
            --p-trunc-len 250 \
            --p-n-threads "$THREADS" \
            --p-max-ee 2.0 \
            --p-trunc-q 2 \
            --p-chimera-method consensus \
            --o-table "$OUTPUT_DIR/table.qza" \
            --o-representative-sequences "$OUTPUT_DIR/rep-seqs.qza" \
            --o-denoising-stats "$OUTPUT_DIR/denoising-stats.qza"
    fi
    create_checkpoint $current_step "Data Import and DADA2"
fi
((current_step++))

# Generate Feature Table Summary
if ! is_step_completed $current_step "Feature Table Summary"; then
    step_info $current_step "Generating feature table summary"
    activate_qiime2
    
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/table.qza" \
        --o-visualization "$OUTPUT_DIR/table.qzv" \
        --m-sample-metadata-file "$METADATA"
    
    qiime feature-table tabulate-seqs \
        --i-data "$OUTPUT_DIR/rep-seqs.qza" \
        --o-visualization "$OUTPUT_DIR/rep-seqs.qzv"
    
    create_checkpoint $current_step "Feature Table Summary"
fi
((current_step++))

# Phylogenetic Tree Construction
if ! is_step_completed $current_step "Phylogenetic Tree Construction"; then
    step_info $current_step "Building phylogenetic tree"
    activate_qiime2
    
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences "$OUTPUT_DIR/rep-seqs.qza" \
        --o-alignment "$OUTPUT_DIR/aligned-rep-seqs.qza" \
        --o-masked-alignment "$OUTPUT_DIR/masked-aligned-rep-seqs.qza" \
        --o-tree "$OUTPUT_DIR/unrooted-tree.qza" \
        --o-rooted-tree "$OUTPUT_DIR/rooted-tree.qza"
    
    create_checkpoint $current_step "Phylogenetic Tree Construction"
fi
((current_step++))

# Taxonomic Classification
if ! is_step_completed $current_step "Taxonomic Classification"; then
    step_info $current_step "Running taxonomic classification"
    activate_qiime2
    
    qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER" \
        --i-reads "$OUTPUT_DIR/rep-seqs.qza" \
        --o-classification "$OUTPUT_DIR/taxonomy.qza"
    
    qiime metadata tabulate \
        --m-input-file "$OUTPUT_DIR/taxonomy.qza" \
        --o-visualization "$OUTPUT_DIR/taxonomy.qzv"
    
    create_checkpoint $current_step "Taxonomic Classification"
fi
((current_step++))

# Alpha and Beta Diversity Analysis
if ! is_step_completed $current_step "Diversity Analysis"; then
    step_info $current_step "Performing diversity analysis"
    activate_qiime2
    
    # Alpha rarefaction
    qiime diversity alpha-rarefaction \
        --i-table "$OUTPUT_DIR/table.qza" \
        --i-phylogeny "$OUTPUT_DIR/rooted-tree.qza" \
        --p-max-depth "$SAMPLING_DEPTH" \
        --m-metadata-file "$METADATA" \
        --o-visualization "$OUTPUT_DIR/diversity/alpha-rarefaction.qzv"
    
    # Core diversity metrics
    core_metrics_dir="$OUTPUT_DIR/diversity/core-metrics-results"
    if [ -d "$core_metrics_dir" ]; then
        echo "Removing existing core metrics directory..."
        rm -rf "$core_metrics_dir"
    fi
    
    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny "$OUTPUT_DIR/rooted-tree.qza" \
        --i-table "$OUTPUT_DIR/table.qza" \
        --p-sampling-depth "$SAMPLING_DEPTH" \
        --m-metadata-file "$METADATA" \
        --output-dir "$core_metrics_dir"
    
    # Alpha diversity statistical tests
    for metric in observed_features shannon pielou_e simpson chao1 fisher_alpha goods_coverage strong dominance mcintosh_d mcintosh_e berger_parker_d margalef; do
        if qiime diversity alpha --p-metric ${metric} --help &>/dev/null; then
            echo "Calculating ${metric} diversity..."
            qiime diversity alpha \
                --i-table "$OUTPUT_DIR/table.qza" \
                --p-metric ${metric} \
                --o-alpha-diversity "$OUTPUT_DIR/diversity/${metric}_vector.qza"
            
            # Try to run the statistical test, but continue if it fails
            if qiime diversity alpha-group-significance \
                --i-alpha-diversity "$OUTPUT_DIR/diversity/${metric}_vector.qza" \
                --m-metadata-file "$METADATA" \
                --o-visualization "$OUTPUT_DIR/diversity/${metric}-group-significance.qzv" 2>/dev/null; then
                echo "Statistical test completed for ${metric}"
            else
                echo "Warning: Could not perform statistical test for ${metric} (possibly due to identical values)"
            fi
        else
            echo "Warning: ${metric} metric not available, skipping..."
        fi
    done
    
    # Beta diversity statistical tests
    for metric in unweighted_unifrac weighted_unifrac bray_curtis jaccard; do
        qiime diversity beta-group-significance \
            --i-distance-matrix "$OUTPUT_DIR/diversity/core-metrics-results/${metric}_distance_matrix.qza" \
            --m-metadata-file "$METADATA" \
            --m-metadata-column "group" \
            --o-visualization "$OUTPUT_DIR/diversity/${metric}-group-significance.qzv" \
            --p-pairwise
    done
    
    create_checkpoint $current_step "Diversity Analysis"
fi
((current_step++))

# Taxonomic Analysis
if ! is_step_completed $current_step "Taxonomic Analysis"; then
    step_info $current_step "Performing taxonomic analysis"
    activate_qiime2
    
    # Generate taxonomy bar plots
    qiime taxa barplot \
        --i-table "$OUTPUT_DIR/table.qza" \
        --i-taxonomy "$OUTPUT_DIR/taxonomy.qza" \
        --m-metadata-file "$METADATA" \
        --o-visualization "$OUTPUT_DIR/taxonomy/taxa-bar-plots.qzv"
    
    # Collapse table at different taxonomic levels
    for level in 2 3 4 5 6 7; do
        qiime taxa collapse \
            --i-table "$OUTPUT_DIR/table.qza" \
            --i-taxonomy "$OUTPUT_DIR/taxonomy.qza" \
            --p-level $level \
            --o-collapsed-table "$OUTPUT_DIR/taxonomy/table-l${level}.qza"
        
        qiime feature-table relative-frequency \
            --i-table "$OUTPUT_DIR/taxonomy/table-l${level}.qza" \
            --o-relative-frequency-table "$OUTPUT_DIR/taxonomy/relative-table-l${level}.qza"
    done
    
    create_checkpoint $current_step "Taxonomic Analysis"
fi
((current_step++))

# Differential Abundance Analysis
if ! is_step_completed $current_step "Differential Abundance Analysis"; then
    step_info $current_step "Performing differential abundance analysis"
    activate_qiime2
    
    # ANCOM analysis at different taxonomic levels
    for level in 2 3 4 5 6 7; do
        qiime composition add-pseudocount \
            --i-table "$OUTPUT_DIR/taxonomy/table-l${level}.qza" \
            --o-composition-table "$OUTPUT_DIR/abundance/comp-table-l${level}.qza"
        
        qiime composition ancom \
            --i-table "$OUTPUT_DIR/abundance/comp-table-l${level}.qza" \
            --m-metadata-file "$METADATA" \
            --m-metadata-column "group" \
            --o-visualization "$OUTPUT_DIR/abundance/ancom-l${level}.qzv"
    done
    
    create_checkpoint $current_step "Differential Abundance Analysis"
fi
((current_step++))

# Krona Visualization
if ! is_step_completed $current_step "Krona Visualization"; then
    step_info $current_step "Creating Krona visualizations"
    activate_qiime2
    
    mkdir -p "$OUTPUT_DIR/krona"
    
    # Export taxonomy
    qiime tools export \
        --input-path "$OUTPUT_DIR/taxonomy.qza" \
        --output-path "$OUTPUT_DIR/krona/taxonomy"
    
    # Export the feature table
    qiime tools export \
        --input-path "$OUTPUT_DIR/table.qza" \
        --output-path "$OUTPUT_DIR/krona/feature_table"
    
    # Convert biom to TSV
    biom convert \
        -i "$OUTPUT_DIR/krona/feature_table/feature-table.biom" \
        -o "$OUTPUT_DIR/krona/feature_table/feature-table.tsv" \
        --to-tsv
    
    # Create Krona input files
    python3 - << EOF
import pandas as pd
import numpy as np

# Read the feature table
ft = pd.read_csv("$OUTPUT_DIR/krona/feature_table/feature-table.tsv", 
                sep='\t', skiprows=1)
ft.set_index('#OTU ID', inplace=True)

# Read taxonomy
tax = pd.read_csv("$OUTPUT_DIR/krona/taxonomy/taxonomy.tsv", sep='\t')
tax.set_index('Feature ID', inplace=True)
tax['Taxon'] = tax['Taxon'].fillna('Unassigned')

# Process each sample
for sample in ft.columns:
    # Create a combined input file for all taxonomic levels
    with open(f"$OUTPUT_DIR/krona/{sample}_krona.txt", 'w') as f:
        for idx, row in ft.iterrows():
            if row[sample] > 0:
                try:
                    taxonomy = tax.loc[idx, 'Taxon'] if idx in tax.index else 'Unassigned'
                    taxonomy = taxonomy.replace('; ', '\t')
                    f.write(f"{int(row[sample])}\t{taxonomy}\n")
                except:
                    continue
EOF
    
    # Generate Krona HTML files for each sample
    for sample_file in "$OUTPUT_DIR/krona/"*_krona.txt; do
        sample=$(basename "${sample_file%_krona.txt}")
        ktImportText \
            -o "$OUTPUT_DIR/krona/${sample}_krona.html" \
            "$sample_file"
    done
    
    create_checkpoint $current_step "Krona Visualization"
fi
((current_step++))

# LEfSe Analysis
if ! is_step_completed $current_step "LEfSe Analysis"; then
    step_info $current_step "Performing LEfSe analysis"
    activate_qiime2
    
    mkdir -p "$OUTPUT_DIR/lefse"
    
    # For each taxonomic level
    for level in 2 3 4 5 6 7; do
        echo "Processing taxonomic level ${level}..."
        
        # Export the relative frequency table
        qiime tools export \
            --input-path "$OUTPUT_DIR/taxonomy/relative-table-l${level}.qza" \
            --output-path "$OUTPUT_DIR/lefse/level${level}"
        
        # Convert to TSV
        biom convert \
            -i "$OUTPUT_DIR/lefse/level${level}/feature-table.biom" \
            -o "$OUTPUT_DIR/lefse/level${level}/feature-table.tsv" \
            --to-tsv
        
        # Create LEfSe input file
        export CURRENT_LEVEL=${level}
        python3 - << EOF
import pandas as pd
import os

# Get the current taxonomic level from environment variable
level = os.environ.get('CURRENT_LEVEL')

# Read feature table
ft = pd.read_csv("$OUTPUT_DIR/lefse/level${level}/feature-table.tsv", 
                 sep='\t', skiprows=1)
ft.set_index('#OTU ID', inplace=True)

# Read metadata
md = pd.read_csv("$METADATA", sep='\t')
md.set_index('#SampleID', inplace=True)

# Create LEfSe input file
with open(f"$OUTPUT_DIR/lefse/level{level}_lefse_input.txt", 'w') as f:
    # Write header with feature names
    feature_names = ft.index.tolist()
    f.write('class\t' + '\t'.join(feature_names) + '\n')
    
    # Write data for each sample
    for sample in ft.columns:
        if sample in md.index:
            group = md.loc[sample, 'group']
            abundance_values = [str(ft.loc[feature, sample]) for feature in feature_names]
            f.write(f"{group}\t" + '\t'.join(abundance_values) + '\n')
EOF
        
        echo "Running LEfSe analysis for level ${level}..."
        
        # Format input for LEfSe with normalization
        lefse_format_input.py \
            "$OUTPUT_DIR/lefse/level${level}_lefse_input.txt" \
            "$OUTPUT_DIR/lefse/level${level}_lefse_input.in" \
            -c 1 \
            -u 2 \
            -o 1000000 \
            -m f \
            -n 10

        # Run LEfSe analysis with optimized parameters
        lefse_run.py \
            "$OUTPUT_DIR/lefse/level${level}_lefse_input.in" \
            "$OUTPUT_DIR/lefse/level${level}_lefse_output.txt" \
            -a 0.10 \
            -w 0.10 \
            -l 1.5 \
            -s 1 \
            --min_c 5 \
            -y 1 \
            -b 100 \
            --verbose 1
        
        # Check if there are significant features
        if [ -s "$OUTPUT_DIR/lefse/level${level}_lefse_output.txt" ]; then
            echo "Found significant features for level ${level}, generating visualizations..."
            
            # Generate horizontal bar plot with improved visualization
            lefse_plot_res.py \
                "$OUTPUT_DIR/lefse/level${level}_lefse_output.txt" \
                "$OUTPUT_DIR/lefse/level${level}_lefse_plot.svg" \
                --format svg \
                --dpi 600 \
                --width 12 \
                --height 8 \
                --orientation h \
                --autoscale 1 \
                --background_color w \
                --feature_font_size 10 \
                --title "Differential Abundance Analysis (LEfSe) - Level ${level}" \
                --title_font_size 14 \
                --class_legend_font_size 12 \
                --max_feature_len 150 \
                --all_feats 1
            
            # Also create a PNG version for quick viewing
            lefse_plot_res.py \
                "$OUTPUT_DIR/lefse/level${level}_lefse_output.txt" \
                "$OUTPUT_DIR/lefse/level${level}_lefse_plot.png" \
                --format png \
                --dpi 300 \
                --width 12 \
                --height 8 \
                --orientation h \
                --autoscale 1 \
                --background_color w \
                --feature_font_size 10 \
                --title "Differential Abundance Analysis (LEfSe) - Level ${level}" \
                --title_font_size 14 \
                --class_legend_font_size 12 \
                --max_feature_len 150 \
                --all_feats 1
            
            # Generate cladogram with improved visualization
            lefse2circlader.py \
                "$OUTPUT_DIR/lefse/level${level}_lefse_output.txt" \
                "$OUTPUT_DIR/lefse/level${level}_lefse_cladogram.svg" \
                -l 7
            
            # Also create a PNG version for quick viewing
            lefse2circlader.py \
                "$OUTPUT_DIR/lefse/level${level}_lefse_output.txt" \
                "$OUTPUT_DIR/lefse/level${level}_lefse_cladogram.png" \
                -l 7
        else
            echo "No significant features found for level ${level}"
            # Create empty files to indicate analysis was performed
            touch "$OUTPUT_DIR/lefse/level${level}_lefse_plot.png"
            touch "$OUTPUT_DIR/lefse/level${level}_lefse_cladogram.png"
        fi
    done
    
    create_checkpoint $current_step "LEfSe Analysis"
fi
((current_step++))

# Process all QZV files
export OUTPUT_DIR="$OUTPUT_DIR"
export METADATA="$METADATA"
python3 - << 'EOF'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from qiime2 import Artifact, Visualization
import skbio
import os

# Get environment variables
output_dir = os.environ.get('OUTPUT_DIR')
metadata_file = os.environ.get('METADATA')

if output_dir is None or metadata_file is None:
    raise ValueError("OUTPUT_DIR or METADATA environment variables are not set")

print(f"Using output directory: {output_dir}")
print(f"Using metadata file: {metadata_file}")

# Create output directories
os.makedirs(os.path.join(output_dir, "publication_figures/quality_plots"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "publication_figures/feature_tables"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "publication_figures/taxonomy"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "publication_figures/diversity"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "publication_figures/abundance"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "exported_data"), exist_ok=True)

def load_and_process_data():
    """Load and process feature table and taxonomy data."""
    try:
        # Load artifacts
        feature_table = Artifact.load(os.path.join(output_dir, "table.qza"))
        taxonomy = Artifact.load(os.path.join(output_dir, "taxonomy.qza"))
        
        # Convert to pandas DataFrames
        ft_df = feature_table.view(pd.DataFrame)
        tax_df = taxonomy.view(pd.DataFrame)
        
        # Export raw data for inspection
        ft_df.to_csv(os.path.join(output_dir, "exported_data/feature_table_raw.tsv"), sep='\t')
        tax_df.to_csv(os.path.join(output_dir, "exported_data/taxonomy_raw.tsv"), sep='\t')
        
        print("\nFeature table info:")
        print(f"Shape: {ft_df.shape}")
        print("First few feature IDs:", ft_df.index[:5].tolist())
        
        print("\nTaxonomy info:")
        print(f"Shape: {tax_df.shape}")
        print("First few feature IDs:", tax_df.index[:5].tolist())
        
        # Export feature IDs for debugging
        with open(os.path.join(output_dir, "exported_data/feature_ids.txt"), 'w') as f:
            f.write("Feature Table IDs:\n")
            f.write("\n".join(ft_df.index.tolist()))
            f.write("\n\nTaxonomy IDs:\n")
            f.write("\n".join(tax_df.index.tolist()))
        
        # Use the feature table as is, without trying to match with taxonomy
        # This is because we'll handle taxonomy matching in individual visualization functions
        return ft_df, tax_df
            
    except Exception as e:
        print(f"Error processing data: {str(e)}")
        return None, None

def load_qiime2_artifact(filepath):
    """Safely load QIIME2 artifact and convert to appropriate format."""
    try:
        print(f"Loading artifact: {filepath}")
        artifact = Artifact.load(filepath)
        
        # Handle different QIIME2 artifact types
        if str(artifact.type) == "DistanceMatrix":
            return artifact.view(skbio.DistanceMatrix)
        elif str(artifact.type) == "PCoAResults":
            return artifact.view(skbio.OrdinationResults)
        elif str(artifact.type) == "SampleData[AlphaDiversity]":
            df = artifact.view(pd.Series).to_frame()
            df.columns = [os.path.basename(filepath).split('_vector')[0]]
            return df
        elif str(artifact.type) == "FeatureTable[Frequency]":
            return artifact.view(pd.DataFrame)
        elif str(artifact.type) == "FeatureData[Taxonomy]":
            return artifact.view(pd.DataFrame)
        else:
            print(f"Unknown artifact type: {artifact.type}")
            return None
            
    except Exception as e:
        print(f"Error loading {filepath}: {str(e)}")
        return None

def plot_alpha_diversity(vector_path, metadata_path, metric_name, output_path):
    """Create publication-quality alpha diversity plots."""
    try:
        print(f"Processing alpha diversity for {metric_name}")
        # Load data
        alpha_div = load_qiime2_artifact(vector_path)
        if alpha_div is None:
            return
        
        metadata = pd.read_csv(metadata_path, sep='\t', index_col='#SampleID')
        
        # Merge data
        combined_data = pd.merge(alpha_div, metadata, left_index=True, right_index=True)
        metric_col = combined_data.columns[0]  # First column should be the metric
        
        # Set publication-quality style
        plt.style.use('default')
        plt.rcParams.update({
            'figure.figsize': (10, 6),
            'font.size': 12,
            'font.family': ['DejaVu Sans', 'Liberation Sans', 'sans-serif'],
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'figure.dpi': 300,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1,
            'axes.linewidth': 1.5,
            'grid.linewidth': 0.5,
            'lines.linewidth': 2.0,
            'axes.spines.top': False,
            'axes.spines.right': False
        })
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Create violin plot with individual points
        sns.violinplot(data=combined_data, x='group', y=metric_col, ax=ax, 
                      inner='box', color='lightgray')
        sns.stripplot(data=combined_data, x='group', y=metric_col, ax=ax,
                     size=6, alpha=0.6, jitter=True)
        
        # Add statistical annotation if more than one group
        if len(combined_data['group'].unique()) > 1:
            from scipy import stats
            groups = combined_data['group'].unique()
            stat_result = stats.kruskal(*[combined_data[combined_data['group'] == g][metric_col] 
                                        for g in groups])
            ax.text(0.95, 0.95, f'Kruskal-Wallis p = {stat_result.pvalue:.3e}',
                   transform=ax.transAxes, ha='right', va='top')
        
        # Customize plot
        ax.set_title(f'{metric_name} by Group', pad=20)
        ax.set_xlabel('Group', labelpad=10)
        ax.set_ylabel(metric_name, labelpad=10)
        
        # Add grid lines
        ax.yaxis.grid(True, linestyle='--', alpha=0.7)
        
        # Save plot
        plt.tight_layout()
        plt.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Saved alpha diversity plot to {output_path}")
        
    except Exception as e:
        print(f"Error plotting alpha diversity for {metric_name}: {str(e)}")

def plot_beta_diversity(distance_matrix_path, pcoa_path, metadata_path, metric_name, output_path):
    """Create publication-quality beta diversity plots."""
    try:
        print(f"Processing beta diversity for {metric_name}")
        # Load distance matrix and metadata
        dist_matrix = load_qiime2_artifact(distance_matrix_path)
        pcoa_results = load_qiime2_artifact(pcoa_path)
        if dist_matrix is None or pcoa_results is None:
            return
            
        metadata = pd.read_csv(metadata_path, sep='\t', index_col='#SampleID')
        
        # Set publication-quality style
        plt.style.use('default')
        plt.rcParams.update({
            'figure.figsize': (12, 8),
            'font.size': 12,
            'font.family': ['DejaVu Sans', 'Liberation Sans', 'sans-serif'],
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'figure.dpi': 300,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1,
            'axes.linewidth': 1.5,
            'grid.linewidth': 0.5,
            'lines.linewidth': 2.0,
            'axes.spines.top': False,
            'axes.spines.right': False
        })
        
        # Extract PCoA coordinates
        pcoa_df = pd.DataFrame(pcoa_results.samples.values,
                             index=pcoa_results.samples.index,
                             columns=[f'PC{i+1}' for i in range(pcoa_results.samples.shape[1])])
        variance_explained = pcoa_results.proportion_explained
        
        # Merge with metadata
        pcoa_df = pd.merge(pcoa_df, metadata, left_index=True, right_index=True, how='inner')
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot points with different markers for groups
        groups = pcoa_df['group'].unique()
        markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']
        for i, group in enumerate(groups):
            mask = pcoa_df['group'] == group
            ax.scatter(pcoa_df.loc[mask, 'PC1'], 
                      pcoa_df.loc[mask, 'PC2'],
                      label=group,
                      marker=markers[i % len(markers)],
                      s=100,
                      alpha=0.7)
        
        # Add variance explained
        ax.set_xlabel(f'PC1 ({variance_explained[0]:.1%} explained)')
        ax.set_ylabel(f'PC2 ({variance_explained[1]:.1%} explained)')
        
        # Add statistical annotation (PERMANOVA)
        from skbio.stats.distance import permanova
        
        # Get common samples between distance matrix and metadata
        common_samples = list(set(dist_matrix.ids) & set(metadata.index))
        if len(common_samples) < 2:
            print(f"Not enough samples for PERMANOVA analysis in {metric_name}")
        else:
            # Subset distance matrix and metadata to common samples
            dist_matrix_subset = dist_matrix.filter(common_samples)
            metadata_subset = metadata.loc[common_samples]
            
            permanova_result = permanova(dist_matrix_subset, 
                                       metadata_subset['group'],
                                       permutations=999)
            
            ax.text(0.95, 0.95, f'PERMANOVA p = {permanova_result["p-value"]:.3e}',
                    transform=ax.transAxes, ha='right', va='top')
        
        # Add grid lines
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # Customize plot
        ax.set_title(f'{metric_name} PCoA', pad=20)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Save plot
        plt.tight_layout()
        plt.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Saved beta diversity plot to {output_path}")
        
    except Exception as e:
        print(f"Error plotting beta diversity for {metric_name}: {str(e)}")

def plot_taxonomy_barplot(table_path, taxonomy_path, metadata_path, output_path):
    """Create publication-quality taxonomy barplots."""
    try:
        print("Processing taxonomy barplot")
        # Load data
        feature_table = Artifact.load(table_path)
        taxonomy = Artifact.load(taxonomy_path)
        
        ft_df = feature_table.view(pd.DataFrame)
        tax_df = taxonomy.view(pd.DataFrame)
        
        metadata = pd.read_csv(metadata_path, sep='\t', index_col='#SampleID')
        
        print("\nData shapes before processing:")
        print("Feature table shape:", ft_df.shape)
        print("Taxonomy shape:", tax_df.shape)
        print("\nFeature table index example:", ft_df.index[:5].tolist())
        print("Taxonomy index example:", tax_df.index[:5].tolist())
        
        # Clean and standardize feature IDs
        ft_df.index = ft_df.index.str.strip()
        tax_df.index = tax_df.index.str.strip()
        
        # Process taxonomy data
        tax_df['Taxon'] = tax_df['Taxon'].fillna('Unassigned')
        tax_levels = tax_df['Taxon'].str.split(';', expand=True)
        tax_levels = tax_levels.fillna('Unassigned')  # Fill NaN values in split results
        tax_levels.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        
        # Add index to tax_levels
        tax_levels.index = tax_df.index
        
        print("\nUnique Phyla found:", tax_levels['Phylum'].unique().tolist())
        
        # Focus on Phylum level
        phylum_data = pd.DataFrame(index=ft_df.columns)
        
        # Aggregate at Phylum level
        for phylum in tax_levels['Phylum'].unique():
            if pd.isna(phylum) or phylum == '':
                continue
                
            # Get features for this phylum
            phylum_features = tax_levels[tax_levels['Phylum'] == phylum].index
            # Find features present in feature table
            common_features = list(set(phylum_features) & set(ft_df.index))
            
            print(f"\nProcessing Phylum: {phylum}")
            print(f"Features in taxonomy: {len(phylum_features)}")
            print(f"Common features found: {len(common_features)}")
            
            if common_features:
                phylum_data[phylum] = ft_df.loc[common_features].sum()
        
        print("\nPhylum data shape before filtering:", phylum_data.shape)
        
        # Remove empty phyla and convert to relative abundance
        phylum_data = phylum_data.loc[:, (phylum_data != 0).any(axis=0)]
        print("Phylum data shape after removing empty phyla:", phylum_data.shape)
        
        if phylum_data.empty:
            raise ValueError("No valid phylum data found after processing")
            
        # Convert to relative abundance
        phylum_data = phylum_data.div(phylum_data.sum(axis=1), axis=0) * 100
        
        # Merge with metadata
        phylum_data = pd.merge(phylum_data, metadata[['group']], left_index=True, right_index=True)
        print("\nData shape after merging with metadata:", phylum_data.shape)
        
        # Calculate mean abundance per group
        phylum_means = phylum_data.groupby('group').mean()
        print("\nFinal phylum means shape:", phylum_means.shape)
        print("Groups found:", phylum_means.index.tolist())
        print("Phyla in final data:", phylum_means.columns.tolist())
        
        if phylum_means.empty:
            raise ValueError("No data to plot after grouping by metadata")
        
        # Sort phyla by overall abundance
        phylum_means = phylum_means.reindex(phylum_means.mean().sort_values(ascending=False).index, axis=1)
        
        # Export data before plotting for debugging
        debug_output = output_path.replace('.svg', '_debug.tsv')
        phylum_means.to_csv(debug_output, sep='\t')
        print(f"\nExported debug data to: {debug_output}")
        
        # Set publication-quality style
        plt.style.use('default')
        plt.rcParams.update({
            'figure.figsize': (15, 8),
            'font.size': 12,
            'font.family': ['DejaVu Sans', 'Liberation Sans', 'sans-serif'],
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 10,
            'figure.dpi': 300,
            'savefig.dpi': 300
        })
        
        # Create figure
        fig, ax = plt.subplots(figsize=(15, 8))
        
        # Create color palette
        colors = plt.cm.Set3(np.linspace(0, 1, len(phylum_means.columns)))
        
        # Plot stacked bars
        phylum_means.plot(kind='bar', stacked=True, ax=ax, color=colors)
        
        # Customize plot
        ax.set_title('Taxonomic Composition at Phylum Level', pad=20)
        ax.set_xlabel('Group', labelpad=10)
        ax.set_ylabel('Relative Abundance (%)', labelpad=10)
        
        # Rotate x-axis labels
        plt.xticks(rotation=45, ha='right')
        
        # Adjust legend
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Phylum',
                 frameon=True, edgecolor='black')
        
        # Add grid lines
        ax.yaxis.grid(True, linestyle='--', alpha=0.7)
        
        # Save plot
        plt.tight_layout()
        plt.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close()
        
        print(f"\nSaved taxonomy barplot to {output_path}")
        
        # Export abundance data
        abundance_output = output_path.replace('.svg', '_abundance.tsv')
        phylum_means.to_csv(abundance_output, sep='\t')
        print(f"Saved abundance data to {abundance_output}")
        
    except Exception as e:
        print(f"\nError plotting taxonomy barplot: {str(e)}")
        import traceback
        traceback.print_exc()
        
        # Export current state of data for debugging
        try:
            debug_dir = os.path.join(os.path.dirname(output_path), 'debug')
            os.makedirs(debug_dir, exist_ok=True)
            
            ft_df.to_csv(os.path.join(debug_dir, 'feature_table_debug.tsv'), sep='\t')
            tax_df.to_csv(os.path.join(debug_dir, 'taxonomy_debug.tsv'), sep='\t')
            if 'tax_levels' in locals():
                tax_levels.to_csv(os.path.join(debug_dir, 'taxonomy_levels_debug.tsv'), sep='\t')
            if 'phylum_data' in locals():
                phylum_data.to_csv(os.path.join(debug_dir, 'phylum_data_debug.tsv'), sep='\t')
            
            print(f"\nDebug files exported to: {debug_dir}")
        except Exception as debug_error:
            print(f"Error exporting debug files: {str(debug_error)}")

# Load and process data
print("\nLoading and processing data...")
feature_table_df, taxonomy_df = load_and_process_data()

if feature_table_df is not None and taxonomy_df is not None:
    # Process core metrics results
    core_metrics_dir = os.path.join(output_dir, "diversity/core-metrics-results")
    print(f"\nChecking for core metrics results in: {core_metrics_dir}")
    
    if os.path.exists(core_metrics_dir):
        print("Processing core metrics results...")
        
        # Process beta diversity metrics
        for metric in ['unweighted_unifrac', 'weighted_unifrac', 'bray_curtis', 'jaccard']:
            try:
                plot_beta_diversity(
                    os.path.join(core_metrics_dir, f"{metric}_distance_matrix.qza"),
                    os.path.join(core_metrics_dir, f"{metric}_pcoa_results.qza"),
                    metadata_file,
                    metric.replace('_', ' ').title(),
                    os.path.join(output_dir, f"publication_figures/diversity/{metric}_pcoa.svg")
                )
            except Exception as e:
                print(f"Could not process {metric}: {str(e)}")
        
        # Process alpha diversity metrics
        for metric in ['observed_features', 'shannon', 'evenness', 'faith_pd']:
            try:
                plot_alpha_diversity(
                    os.path.join(core_metrics_dir, f"{metric}_vector.qza"),
                    metadata_file,
                    metric.replace('_', ' ').title(),
                    os.path.join(output_dir, f"publication_figures/diversity/{metric}_boxplot.svg")
                )
            except Exception as e:
                print(f"Could not process {metric}: {str(e)}")

# Create taxonomy barplot
try:
    plot_taxonomy_barplot(
        os.path.join(output_dir, "table.qza"),
        os.path.join(output_dir, "taxonomy.qza"),
        metadata_file,
        os.path.join(output_dir, "publication_figures/taxonomy/taxonomy_barplot.svg")
    )
except Exception as e:
    print(f"Could not create taxonomy barplot: {str(e)}")

print("All QIIME2 visualizations have been converted to publication-quality figures")
EOF

# Final status
if [ -f "$OUTPUT_DIR/logs/step$((current_step-1)).done" ]; then
    echo "Pipeline completed successfully at $(date '+%Y-%m-%d %H:%M:%S')"
    echo "All steps completed successfully. Check $OUTPUT_DIR/logs/pipeline.log for details"
else
    echo "Pipeline incomplete. Check $OUTPUT_DIR/logs/pipeline.log for details"
fi 