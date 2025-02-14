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

# QIIME2 Data Preparation Pipeline
if ! is_step_completed $current_step "Data Import and DADA2"; then
    step_info $current_step "Starting QIIME2 data preparation"
    
    # Ensure QIIME2 is activated
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

# 2. Generate Feature Table Summary
if ! is_step_completed $current_step "Feature Table Summary"; then
    step_info $current_step "Generating feature table summary"
    
    # Ensure QIIME2 is activated
    activate_qiime2
    
    qiime feature-table summarize \
        --i-table "$OUTPUT_DIR/table.qza" \
        --o-visualization "$OUTPUT_DIR/table.qzv" \
        --m-sample-metadata-file "$METADATA"
    create_checkpoint $current_step "Feature Table Summary"
fi
((current_step++))

# 3. Phylogenetic Tree Construction
if ! is_step_completed $current_step "Phylogenetic Tree Construction"; then
    step_info $current_step "Building phylogenetic tree"
    
    # Ensure QIIME2 is activated
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

# 4. Taxonomic Classification
if ! is_step_completed $current_step "Taxonomic Classification"; then
    step_info $current_step "Running taxonomic classification"
    
    # Ensure QIIME2 is activated
    activate_qiime2
    
    qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER" \
        --i-reads "$OUTPUT_DIR/rep-seqs.qza" \
        --o-classification "$OUTPUT_DIR/taxonomy.qza"
    create_checkpoint $current_step "Taxonomic Classification"
fi
((current_step++))

# 5. Metadata and Taxonomy Validation
if ! is_step_completed $current_step "Metadata and Taxonomy Validation"; then
    step_info $current_step "Validating metadata and taxonomy"
    
    # Ensure QIIME2 is activated
    activate_qiime2
    
    # Validate metadata against feature table
    step_info $current_step "Checking metadata and feature table compatibility"
    qiime metadata tabulate \
        --m-input-file "$OUTPUT_DIR/table.qza" \
        --o-visualization "$OUTPUT_DIR/feature-table-metadata.qzv"
    
    # Examine taxonomy assignments
    step_info $current_step "Examining taxonomy assignments"
    qiime metadata tabulate \
        --m-input-file "$OUTPUT_DIR/taxonomy.qza" \
        --o-visualization "$OUTPUT_DIR/taxonomy.qzv"
    
    # Merge feature table with taxonomy
    step_info $current_step "Merging feature table with taxonomy"
    qiime feature-table merge-taxa \
        --i-data "$OUTPUT_DIR/taxonomy.qza" \
        --o-merged-data "$OUTPUT_DIR/merged-taxonomy.qza"
    
    create_checkpoint $current_step "Metadata and Taxonomy Validation"
fi
((current_step++))

# 6. Export files for Phyloseq
if ! is_step_completed $current_step "Exporting files for Phyloseq analysis"; then
    step_info $current_step "Exporting files for Phyloseq analysis"
    
    # Ensure QIIME2 is activated
    activate_qiime2
    
    mkdir -p "$OUTPUT_DIR/exported"

    # Export feature table
    step_info $current_step "Exporting feature table"
    qiime tools export \
        --input-path "$OUTPUT_DIR/table.qza" \
        --output-path "$OUTPUT_DIR/exported"

    # Export taxonomy
    step_info $current_step "Exporting taxonomy"
    qiime tools export \
        --input-path "$OUTPUT_DIR/taxonomy.qza" \
        --output-path "$OUTPUT_DIR/exported"

    # Export tree
    step_info $current_step "Exporting tree"
    qiime tools export \
        --input-path "$OUTPUT_DIR/rooted-tree.qza" \
        --output-path "$OUTPUT_DIR/exported"

    # Convert biom to TSV and fix sample names
    step_info $current_step "Converting biom to TSV and fixing sample names"
    biom convert \
        -i "$OUTPUT_DIR/exported/feature-table.biom" \
        -o "$OUTPUT_DIR/exported/feature-table_temp.txt" \
        --to-tsv

    # Fix sample names in feature table
    step_info $current_step "Fixing sample names in feature table"
    sed 's/^# Constructed from biom file/# QIIME2 Feature Table/' "$OUTPUT_DIR/exported/feature-table_temp.txt" | \
        sed 's/^# OTU ID/# Feature ID/' | \
        sed 's/^X//g' | \
        sed 's/\.0$//g' | \
        sed 's/\.0\.[0-9]*$//g' > "$OUTPUT_DIR/exported/feature-table.txt"

    # Clean up temporary files
    step_info $current_step "Cleaning up temporary files"
    rm "$OUTPUT_DIR/exported/feature-table_temp.txt"

    create_checkpoint $current_step "Exporting files for Phyloseq analysis"
fi
((current_step++))

# 7. Create and run R script for Phyloseq analysis
if ! is_step_completed $current_step "Creating R script for Phyloseq analysis"; then
    step_info $current_step "Creating R script for Phyloseq analysis"
    cat > "$OUTPUT_DIR/phyloseq_analysis.R" << 'EOF'
#!/usr/bin/env Rscript

# Function to check and install packages only if they're not already installed
install_if_missing <- function(package_name) {
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
        message(paste("Installing package:", package_name))
        if (package_name %in% c("phyloseq", "DESeq2", "microbiome")) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager", repos = "http://cran.us.r-project.org", quiet = TRUE)
            }
            BiocManager::install(package_name, quiet = TRUE)
        } else {
            install.packages(package_name, repos = "http://cran.us.r-project.org", quiet = TRUE)
        }
    }
}

# List of required packages
required_packages <- c("phyloseq", "ggplot2", "vegan", "DESeq2", "microbiome", 
                      "data.table", "dplyr", "tidyr", "RColorBrewer", "reshape2")

# Install missing packages
invisible(sapply(required_packages, install_if_missing))

# Load all required packages
invisible(sapply(required_packages, library, character.only = TRUE))

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript phyloseq_analysis.R <exported_dir> <metadata_file>")
}

exported_dir <- args[1]
metadata_file <- args[2]

# Create output directories
dir.create(file.path(exported_dir, "..", "alpha_diversity"), showWarnings = FALSE)
dir.create(file.path(exported_dir, "..", "beta_diversity"), showWarnings = FALSE)
dir.create(file.path(exported_dir, "..", "abundance_analysis"), showWarnings = FALSE)
dir.create(file.path(exported_dir, "..", "differential_abundance"), showWarnings = FALSE)

# Import data
otu_table <- read.table(file.path(exported_dir, "feature-table.txt"), header=T, row.names=1, skip=1)
taxonomy <- read.table(file.path(exported_dir, "taxonomy.tsv"), header=T, sep="\t", row.names=1)
metadata <- read.table(metadata_file, header=T, sep="\t", row.names=1, check.names=FALSE)
tree <- read_tree(file.path(exported_dir, "tree.nwk"))

# Print data summary
cat("\nData Summary:\n")
cat("Number of samples in feature table:", ncol(otu_table), "\n")
cat("Number of samples in metadata:", nrow(metadata), "\n")
cat("Number of features:", nrow(otu_table), "\n")
cat("Number of taxonomic assignments:", nrow(taxonomy), "\n")

# Create phyloseq object
OTU <- otu_table(as.matrix(otu_table), taxa_are_rows=TRUE)
TAX <- tax_table(as.matrix(taxonomy))
META <- sample_data(metadata)
physeq <- phyloseq(OTU, TAX, META, tree)

# Print phyloseq object summary
cat("\nPhyloseq object summary:\n")
print(physeq)

# 1. Alpha Diversity Analysis
cat("\nCalculating alpha diversity metrics...\n")
alpha_div <- estimate_richness(physeq, measures=c("Observed", "Chao1", "Shannon", "Simpson", "Fisher"))
# Calculate Pielou's evenness
alpha_div$Pielou <- alpha_div$Shannon / log(alpha_div$Observed)
write.csv(alpha_div, file.path(exported_dir, "..", "alpha_diversity/alpha_diversity_metrics.csv"))

# Alpha diversity plots
pdf(file.path(exported_dir, "..", "alpha_diversity/alpha_diversity_boxplots.pdf"), width=12, height=8)
par(mfrow=c(2,3))
for(metric in c("Observed", "Chao1", "Shannon", "Simpson", "Fisher", "Pielou")) {
    boxplot(alpha_div[[metric]] ~ sample_data(physeq)$Group,
            main=metric,
            xlab="Group",
            ylab=metric)
}
dev.off()

# 2. Beta Diversity Analysis
cat("\nCalculating beta diversity metrics...\n")

# Calculate distance matrices
dist_methods <- c("bray", "jaccard", "unifrac", "wunifrac")
dist_matrices <- list()
for(method in dist_methods) {
    dist_matrices[[method]] <- phyloseq::distance(physeq, method=method)
}

# Ordination methods
ord_methods <- c("NMDS", "PCoA")
for(dist_method in names(dist_matrices)) {
    for(ord_method in ord_methods) {
        tryCatch({
            ord <- ordinate(physeq, method=ord_method, distance=dist_matrices[[dist_method]])
            
            # Create ordination plot
            p <- plot_ordination(physeq, ord, color="Group") +
                ggtitle(paste(ord_method, "-", dist_method)) +
                theme_bw()
            
            ggsave(file.path(exported_dir, "..", paste0("beta_diversity/", tolower(ord_method), "_", 
                         dist_method, ".pdf")), p, width=8, height=6)
            
            # Save stress plot for NMDS
            if(ord_method == "NMDS") {
                pdf(file.path(exported_dir, "..", paste0("beta_diversity/stress_", dist_method, ".pdf")))
                stressplot(ord)
                dev.off()
            }
        }, error=function(e) {
            cat("Error in ordination:", dist_method, ord_method, "\n")
            print(e)
        })
    }
}

# Statistical tests
cat("\nPerforming statistical tests...\n")
sink(file.path(exported_dir, "..", "beta_diversity/statistical_tests.txt"))

for(dist_method in names(dist_matrices)) {
    cat("\n=== Results for", dist_method, "distance ===\n")
    
    # ANOSIM
    ano <- anosim(dist_matrices[[dist_method]], metadata$Group)
    cat("\nANOSIM Results:\n")
    print(ano)
    
    # ADONIS (PERMANOVA)
    ado <- adonis2(dist_matrices[[dist_method]] ~ Group, data=metadata)
    cat("\nADONIS Results:\n")
    print(ado)
    
    # MRPP
    mrpp_result <- mrpp(dist_matrices[[dist_method]], metadata$Group)
    cat("\nMRPP Results:\n")
    print(mrpp_result)
}
sink()

# 3. Relative Abundance Analysis
cat("\nPerforming relative abundance analysis...\n")

# Transform to relative abundance
rel_abund <- transform_sample_counts(physeq, function(x) x/sum(x)*100)

# Function to get top N taxa at a given taxonomic level
get_top_taxa <- function(physeq, level, n=10) {
    tax_table <- tax_table(physeq)
    otu_table <- otu_table(physeq)
    
    # Sum abundances for each taxon
    tax_sum <- tapply(colSums(otu_table), tax_table[,level], sum)
    top_taxa <- names(sort(tax_sum, decreasing=TRUE)[1:n])
    
    return(top_taxa)
}

# Create abundance tables and plots for each taxonomic level
tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
for(level in tax_levels) {
    # Get top 10 taxa
    top_taxa <- get_top_taxa(rel_abund, level)
    
    # Create abundance table
    abund_table <- tax_glom(rel_abund, taxrank=level)
    tax_mat <- tax_table(abund_table)[,level]
    otu_mat <- otu_table(abund_table)
    abund_table_df <- data.frame(t(otu_mat))
    colnames(abund_table_df) <- tax_mat
    
    # Save abundance table
    write.csv(abund_table_df, 
              file.path(exported_dir, "..", paste0("abundance_analysis/abundance_table_", tolower(level), ".csv")))
    
    # Create stacked bar plot
    mdf <- melt(abund_table_df)
    colnames(mdf) <- c("Sample", "Taxon", "Abundance")
    
    p <- ggplot(mdf, aes(x=Sample, y=Abundance, fill=Taxon)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        ggtitle(paste("Relative Abundance -", level)) +
        scale_fill_brewer(palette="Set3")
    
    ggsave(file.path(exported_dir, "..", paste0("abundance_analysis/abundance_plot_", tolower(level), ".pdf")),
           p, width=12, height=6)
}

# 4. Differential Abundance Analysis with DESeq2
cat("\nPerforming differential abundance analysis with DESeq2...\n")

# Convert to DESeq2 object
dds <- phyloseq_to_deseq2(physeq, ~ Group)

# Estimate size factors and dispersions
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# Perform Wald test
dds <- DESeq(dds, test="Wald", fitType="parametric")
res <- results(dds)

# Add taxonomy information to results
tax_table <- tax_table(physeq)
res_df <- as.data.frame(res)
res_df$Feature_ID <- rownames(res_df)
res_df <- merge(res_df, as.data.frame(tax_table), by="row.names")
rownames(res_df) <- res_df$Row.names
res_df$Row.names <- NULL

# Save results
write.csv(res_df, file.path(exported_dir, "..", "differential_abundance/deseq2_results.csv"))

# Create volcano plot
pdf(file.path(exported_dir, "..", "differential_abundance/volcano_plot.pdf"))
plot(res$log2FoldChange, -log10(res$padj),
     main="Volcano Plot",
     xlab="log2 Fold Change",
     ylab="-log10 adjusted p-value",
     pch=20)
abline(h=-log10(0.05), col="red", lty=2)
abline(v=c(-1,1), col="red", lty=2)
dev.off()

# Save session info
writeLines(capture.output(sessionInfo()), file.path(exported_dir, "..", "session_info.txt"))

cat("\nAnalysis complete! Check the output directories for results.\n")
EOF

    chmod +x "$OUTPUT_DIR/phyloseq_analysis.R"

    # Run Phyloseq analysis
    step_info $current_step "Running Phyloseq analysis"
    Rscript "$OUTPUT_DIR/phyloseq_analysis.R" "$OUTPUT_DIR/exported" "$METADATA"

    step_info $current_step "Analysis complete! Results are in $OUTPUT_DIR"

    create_checkpoint $current_step "Phyloseq Analysis"
fi
((current_step++))

# Final status
if [ -f "$OUTPUT_DIR/logs/step$((current_step-1)).done" ]; then
    echo "Pipeline completed successfully at $(date '+%Y-%m-%d %H:%M:%S')" > "$OUTPUT_DIR/pipeline_status.txt"
    echo "All steps completed successfully. Check $OUTPUT_DIR/logs/pipeline.log for details"
else
    echo "Pipeline incomplete. Check $OUTPUT_DIR/logs/pipeline.log for details"
fi 