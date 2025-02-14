#!/usr/bin/env python3

import os
import sys
import argparse
import logging
from datetime import datetime
import qiime2
from qiime2 import Artifact, Visualization
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, anosim, permdisp
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import confusion_matrix, classification_report, roc_curve, auc
from scipy import stats
import subprocess
import networkx as nx
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import shutil

class MicrobiomeAnalysisPipeline:
    def __init__(self, input_dir, metadata, output_dir, classifier, threads=8, 
                 sampling_depth=10000, paired=True, force=False):
        """Initialize the pipeline with parameters."""
        self.input_dir = input_dir
        self.metadata = metadata
        self.output_dir = output_dir
        self.classifier = classifier
        self.threads = threads
        self.sampling_depth = sampling_depth
        self.paired = paired
        self.force = force
        
        # Setup logging
        self.setup_logging()
        
        # Create output directories
        self.create_directories()
        
        # Initialize step tracking
        self.current_step = 1
        
        # Store artifacts for later use
        self.artifacts = {}
        
    def setup_logging(self):
        """Setup logging configuration."""
        log_dir = os.path.join(self.output_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)
        
        log_file = os.path.join(log_dir, "pipeline.log")
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def create_directories(self):
        """Create necessary output directories."""
        dirs = [
            "logs",
            "diversity",
            "taxonomy",
            "abundance",
            "exported_data"
        ]
        for dir_name in dirs:
            os.makedirs(os.path.join(self.output_dir, dir_name), exist_ok=True)
            
    def create_checkpoint(self, step_name):
        """Create a checkpoint file for completed steps."""
        checkpoint_file = os.path.join(self.output_dir, "logs", f"step{self.current_step}.done")
        with open(checkpoint_file, 'w') as f:
            f.write(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        self.logger.info(f"Completed step {self.current_step}: {step_name}")
        
    def is_step_completed(self, step_name):
        """Check if a step has been completed."""
        checkpoint_file = os.path.join(self.output_dir, "logs", f"step{self.current_step}.done")
        if os.path.exists(checkpoint_file) and not self.force:
            self.logger.info(f"Step {self.current_step} ({step_name}) already completed, skipping...")
            return True
        return False
        
    def process_sequence_statistics(self):
        """Process sequence statistics and generate quality visualizations."""
        step_name = "Sequence Statistics"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            from qiime2.plugins import demux
            
            # Import sequences based on whether data is paired or single-end
            if self.paired:
                self.logger.info("Importing paired-end data")
                imported_seqs = Artifact.import_data(
                    'SampleData[PairedEndSequencesWithQuality]',
                    self.input_dir,
                    'CasavaOneEightSingleLanePerSampleDirFmt'
                )
                demux_file = os.path.join(self.output_dir, "demux-paired-end.qza")
                viz_file = os.path.join(self.output_dir, "demux-paired-end.qzv")
            else:
                self.logger.info("Importing single-end data")
                imported_seqs = Artifact.import_data(
                    'SampleData[SequencesWithQuality]',
                    self.input_dir,
                    'CasavaOneEightSingleLanePerSampleDirFmt'
                )
                demux_file = os.path.join(self.output_dir, "demux-single-end.qza")
                viz_file = os.path.join(self.output_dir, "demux-single-end.qzv")
            
            # Save imported sequences
            imported_seqs.save(demux_file)
            
            # Store artifact for later use
            self.artifacts['demux'] = imported_seqs
            
            # Generate visualization using QIIME2 demux plugin
            demux_vis = demux.visualizers.summarize(data=imported_seqs)
            demux_vis.visualization.save(viz_file)
            
            # Export statistics to a more accessible format
            stats_dir = os.path.join(self.output_dir, "exported_data", "sequence_stats")
            os.makedirs(stats_dir, exist_ok=True)
            
            # Extract and save basic statistics from the demux visualization
            per_sample_fastq_counts = os.path.join(stats_dir, "per_sample_fastq_counts.tsv")
            
            # Use qiime tools export to get the statistics
            export_dir = os.path.join(stats_dir, "temp_export")
            os.makedirs(export_dir, exist_ok=True)
            
            demux_vis.visualization.export_data(export_dir)
            
            # Read the per-sample counts
            counts_file = os.path.join(export_dir, "per-sample-fastq-counts.tsv")
            stats_df = pd.read_csv(counts_file, sep='\t', index_col=0)
            
            # Save the statistics
            stats_df.to_csv(per_sample_fastq_counts, sep='\t')
            
            # Calculate and save summary statistics
            summary_stats = {
                'total_sequences': stats_df['forward sequence count'].sum(),
                'mean_sequences_per_sample': stats_df['forward sequence count'].mean(),
                'min_sequences': stats_df['forward sequence count'].min(),
                'max_sequences': stats_df['forward sequence count'].max(),
                'std_sequences': stats_df['forward sequence count'].std()
            }
            
            with open(os.path.join(stats_dir, "summary_stats.txt"), 'w') as f:
                for stat, value in summary_stats.items():
                    f.write(f"{stat}: {value}\n")
            
            # Clean up temporary export directory
            shutil.rmtree(export_dir)
            
            self.logger.info("Sequence statistics processing completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def run_dada2_denoising(self):
        """Run DADA2 denoising on the sequence data."""
        step_name = "DADA2 Denoising"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            from qiime2.plugins import dada2
            
            # Get the demux artifact
            demux = self.artifacts.get('demux')
            if demux is None:
                raise ValueError("Demultiplexed sequences not found. Run process_sequence_statistics first.")
            
            # Set DADA2 parameters
            dada2_params = {
                'n_threads': self.threads,
                'max_ee_f': 2.0,
                'max_ee_r': 2.0 if self.paired else None,
                'trunc_q': 2,
                'chimera_method': 'consensus',
                'min_fold_parent_over_abundance': 1.0
            }
            
            # Add paired-end specific parameters
            if self.paired:
                self.logger.info("Running paired-end DADA2 denoising")
                dada2_params.update({
                    'trim_left_f': 0,
                    'trim_left_r': 0,
                    'trunc_len_f': 250,
                    'trunc_len_r': 250
                })
                
                # Run DADA2 paired-end
                dada2_results = dada2.methods.denoise_paired(
                    demultiplexed_seqs=demux,
                    **dada2_params
                )
            else:
                self.logger.info("Running single-end DADA2 denoising")
                dada2_params.update({
                    'trim_left': 0,
                    'trunc_len': 250
                })
                
                # Run DADA2 single-end
                dada2_results = dada2.methods.denoise_single(
                    demultiplexed_seqs=demux,
                    **dada2_params
                )
            
            # Save DADA2 results
            dada2_results.table.save(os.path.join(self.output_dir, "table.qza"))
            dada2_results.representative_sequences.save(os.path.join(self.output_dir, "rep-seqs.qza"))
            dada2_results.denoising_stats.save(os.path.join(self.output_dir, "denoising-stats.qza"))
            
            # Store artifacts for later use
            self.artifacts['table'] = dada2_results.table
            self.artifacts['rep_seqs'] = dada2_results.representative_sequences
            self.artifacts['denoising_stats'] = dada2_results.denoising_stats
            
            # Generate visualizations using QIIME2 plugins
            from qiime2.plugins import feature_table, metadata
            from qiime2.plugins import dada2
            
            # Create feature table visualization
            metadata_obj = qiime2.Metadata.load(self.metadata)
            table_viz = feature_table.visualizers.summarize(
                table=dada2_results.table,
                sample_metadata=metadata_obj
            )
            table_viz.visualization.save(os.path.join(self.output_dir, "table.qzv"))
            
            # Create representative sequences visualization using feature-table plugin
            rep_seqs_viz = feature_table.visualizers.tabulate_seqs(
                data=dada2_results.representative_sequences
            )
            rep_seqs_viz.visualization.save(os.path.join(self.output_dir, "rep-seqs.qzv"))
            
            # Create denoising stats visualization using metadata plugin
            stats_viz = metadata.visualizers.tabulate(
                metadata=qiime2.Metadata.load(os.path.join(self.output_dir, "denoising-stats.qza"))
            )
            stats_viz.visualization.save(os.path.join(self.output_dir, "denoising-stats.qzv"))
            
            self.logger.info("DADA2 denoising completed successfully")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def generate_feature_table_summary(self):
        """Generate feature table summary and visualizations."""
        step_name = "Feature Table Summary"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            from qiime2.plugins import feature_table
            
            # Get the feature table artifact
            table = self.artifacts.get('table')
            if table is None:
                raise ValueError("Feature table not found. Run DADA2 denoising first.")
            
            # Get the representative sequences artifact
            rep_seqs = self.artifacts.get('rep_seqs')
            if rep_seqs is None:
                raise ValueError("Representative sequences not found. Run DADA2 denoising first.")
            
            # Generate feature table summary
            self.logger.info("Generating feature table summary")
            metadata = qiime2.Metadata.load(self.metadata)
            
            table_summary = feature_table.visualizers.summarize(
                table=table,
                sample_metadata=metadata
            )
            table_summary.visualization.save(
                os.path.join(self.output_dir, "table.qzv")
            )
            
            # Generate representative sequences summary
            self.logger.info("Generating representative sequences summary")
            seq_summary = feature_table.visualizers.tabulate_seqs(
                data=rep_seqs
            )
            seq_summary.visualization.save(
                os.path.join(self.output_dir, "rep-seqs.qzv")
            )
            
            # Export feature table statistics
            stats_dir = os.path.join(self.output_dir, "exported_data", "feature_table_stats")
            os.makedirs(stats_dir, exist_ok=True)
            
            # Export feature table to biom format
            feature_table_biom = os.path.join(stats_dir, "feature-table.biom")
            table.export_data(feature_table_biom)
            
            # Convert to TSV for easier viewing
            feature_table_tsv = os.path.join(stats_dir, "feature-table.tsv")
            table_df = table.view(pd.DataFrame)
            table_df.to_csv(feature_table_tsv, sep='\t')
            
            # Calculate and save summary statistics
            summary_stats = {
                'total_features': table_df.shape[0],
                'total_samples': table_df.shape[1],
                'total_frequency': table_df.sum().sum(),
                'mean_frequency_per_sample': table_df.sum(axis=0).mean(),
                'mean_features_per_sample': (table_df > 0).sum(axis=0).mean()
            }
            
            with open(os.path.join(stats_dir, "summary_stats.txt"), 'w') as f:
                for stat, value in summary_stats.items():
                    f.write(f"{stat}: {value}\n")
            
            self.logger.info("Feature table summary completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1
    
    def build_phylogenetic_tree(self):
        """Build phylogenetic tree using MAFFT and FastTree."""
        step_name = "Phylogenetic Tree Construction"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            from qiime2.plugins import phylogeny
            
            # Get representative sequences
            rep_seqs = self.artifacts.get('rep_seqs')
            if rep_seqs is None:
                raise ValueError("Representative sequences not found. Run DADA2 denoising first.")
            
            self.logger.info("Performing multiple sequence alignment with MAFFT")
            aligned_seqs = phylogeny.methods.align_to_tree_mafft_fasttree(
                sequences=rep_seqs,
                n_threads=self.threads
            )
            
            # Save all outputs
            aligned_seqs.alignment.save(
                os.path.join(self.output_dir, "aligned-rep-seqs.qza")
            )
            aligned_seqs.masked_alignment.save(
                os.path.join(self.output_dir, "masked-aligned-rep-seqs.qza")
            )
            aligned_seqs.tree.save(
                os.path.join(self.output_dir, "unrooted-tree.qza")
            )
            aligned_seqs.rooted_tree.save(
                os.path.join(self.output_dir, "rooted-tree.qza")
            )
            
            # Store artifacts for later use
            self.artifacts['aligned_seqs'] = aligned_seqs.alignment
            self.artifacts['masked_aligned_seqs'] = aligned_seqs.masked_alignment
            self.artifacts['unrooted_tree'] = aligned_seqs.tree
            self.artifacts['rooted_tree'] = aligned_seqs.rooted_tree
            
            # Export tree in Newick format for visualization in other tools
            tree_dir = os.path.join(self.output_dir, "exported_data", "phylogenetic_tree")
            os.makedirs(tree_dir, exist_ok=True)
            
            # Export unrooted tree
            unrooted_tree_file = os.path.join(tree_dir, "unrooted_tree.nwk")
            with open(unrooted_tree_file, 'w') as f:
                f.write(str(aligned_seqs.tree.view(scikit_bio.TreeNode)))
            
            # Export rooted tree
            rooted_tree_file = os.path.join(tree_dir, "rooted_tree.nwk")
            with open(rooted_tree_file, 'w') as f:
                f.write(str(aligned_seqs.rooted_tree.view(scikit_bio.TreeNode)))
            
            self.logger.info("Phylogenetic tree construction completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def run_taxonomic_classification(self):
        """Run taxonomic classification and generate visualizations."""
        step_name = "Taxonomic Classification"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            from qiime2.plugins import feature_classifier
            from qiime2.plugins import taxa
            
            # Get representative sequences
            rep_seqs = self.artifacts.get('rep_seqs')
            if rep_seqs is None:
                raise ValueError("Representative sequences not found. Run DADA2 denoising first.")
            
            # Get feature table
            table = self.artifacts.get('table')
            if table is None:
                raise ValueError("Feature table not found. Run DADA2 denoising first.")
            
            # Load the classifier
            self.logger.info("Loading pre-trained classifier")
            classifier = Artifact.load(self.classifier)
            
            # Run taxonomic classification
            self.logger.info("Running taxonomic classification")
            taxonomy = feature_classifier.methods.classify_sklearn(
                reads=rep_seqs,
                classifier=classifier,
                n_jobs=self.threads
            )
            
            # Save taxonomy artifact
            taxonomy_file = os.path.join(self.output_dir, "taxonomy.qza")
            taxonomy.classification.save(taxonomy_file)
            
            # Store artifact for later use
            self.artifacts['taxonomy'] = taxonomy.classification
            
            # Generate taxonomy visualization
            taxonomy_viz = taxa.visualizers.tabulate(
                data=taxonomy.classification
            )
            taxonomy_viz.visualization.save(
                os.path.join(self.output_dir, "taxonomy.qzv")
            )
            
            # Generate taxonomy barplots
            barplots = taxa.visualizers.barplot(
                table=table,
                taxonomy=taxonomy.classification,
                metadata=qiime2.Metadata.load(self.metadata)
            )
            barplots.visualization.save(
                os.path.join(self.output_dir, "taxa-bar-plots.qzv")
            )
            
            # Export taxonomy data
            tax_dir = os.path.join(self.output_dir, "exported_data", "taxonomy")
            os.makedirs(tax_dir, exist_ok=True)
            
            # Convert taxonomy data to DataFrame
            tax_df = taxonomy.classification.view(pd.DataFrame)
            tax_df.to_csv(os.path.join(tax_dir, "taxonomy.tsv"), sep='\t')
            
            # Generate custom visualizations
            self.logger.info("Generating custom taxonomy visualizations")
            
            # Create publication figures directory
            pub_dir = os.path.join(self.output_dir, "publication_figures", "taxonomy")
            os.makedirs(pub_dir, exist_ok=True)
            
            # Process taxonomy data for visualization
            tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            
            # Split taxonomy strings into separate columns
            tax_split = tax_df['Taxon'].str.split(';', expand=True)
            tax_split.columns = tax_levels
            
            # Clean up taxonomy names
            for col in tax_levels:
                tax_split[col] = tax_split[col].str.replace(r'^[a-z]__', '', regex=True)
                tax_split[col] = tax_split[col].fillna('Unassigned')
            
            # Get feature table as DataFrame
            table_df = table.view(pd.DataFrame)
            
            # Calculate relative abundance
            rel_abundance = table_df.div(table_df.sum(axis=0), axis=1) * 100
            
            # Load metadata
            metadata = pd.read_csv(self.metadata, sep='\t', index_col='#SampleID')
            
            # Generate stacked barplots for each taxonomic level
            for level in tax_levels:
                self.logger.info(f"Generating {level} level visualization")
                
                # Aggregate abundances at the current taxonomic level
                level_abundance = pd.DataFrame(index=rel_abundance.columns)
                
                for taxon in tax_split[level].unique():
                    if taxon == 'Unassigned':
                        continue
                    mask = tax_split[level] == taxon
                    if mask.any():
                        features = tax_split[mask].index
                        level_abundance[taxon] = rel_abundance.loc[features].sum()
                
                # Sort taxa by mean abundance
                taxa_means = level_abundance.mean()
                top_taxa = taxa_means.sort_values(ascending=False).head(10).index
                
                # Create figure
                plt.figure(figsize=(15, 8))
                
                # Create color palette
                colors = plt.cm.Set3(np.linspace(0, 1, len(top_taxa)))
                
                # Plot stacked bars
                ax = level_abundance[top_taxa].plot(
                    kind='bar',
                    stacked=True,
                    color=colors,
                    width=0.8
                )
                
                # Customize plot
                plt.title(f'Top 10 Most Abundant {level}s', pad=20, fontsize=14)
                plt.xlabel('Samples', labelpad=10, fontsize=12)
                plt.ylabel('Relative Abundance (%)', labelpad=10, fontsize=12)
                plt.xticks(rotation=45, ha='right')
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=level)
                plt.grid(axis='y', linestyle='--', alpha=0.7)
                
                # Adjust layout and save
                plt.tight_layout()
                plt.savefig(
                    os.path.join(pub_dir, f"{level.lower()}_barplot.svg"),
                    bbox_inches='tight',
                    dpi=300
                )
                plt.close()
                
                # Save abundance data
                level_abundance.to_csv(
                    os.path.join(tax_dir, f"{level.lower()}_abundance.tsv"),
                    sep='\t'
                )
            
            # Generate heatmap of top 20 taxa at genus level
            genus_abundance = pd.DataFrame(index=rel_abundance.columns)
            for taxon in tax_split['Genus'].unique():
                if taxon == 'Unassigned':
                    continue
                mask = tax_split['Genus'] == taxon
                if mask.any():
                    features = tax_split[mask].index
                    genus_abundance[taxon] = rel_abundance.loc[features].sum()
            
            # Get top 20 genera
            top_genera = genus_abundance.mean().sort_values(ascending=False).head(20).index
            
            # Create heatmap
            plt.figure(figsize=(12, 8))
            sns.heatmap(
                genus_abundance[top_genera].T,
                cmap='YlOrRd',
                robust=True,
                cbar_kws={'label': 'Relative Abundance (%)'}
            )
            plt.title('Top 20 Most Abundant Genera', pad=20)
            plt.xlabel('Samples')
            plt.ylabel('Genus')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(
                os.path.join(pub_dir, "genus_heatmap.svg"),
                bbox_inches='tight',
                dpi=300
            )
            plt.close()
            
            # Generate Krona visualization if available
            try:
                import krona
                self.logger.info("Generating Krona visualization")
                
                # Prepare data for Krona
                krona_data = []
                for idx, row in tax_df.iterrows():
                    if idx in table_df.index:
                        abundance = table_df.loc[idx].sum()
                        taxonomy = row['Taxon'].split(';')
                        taxonomy = [t.split('__')[-1] for t in taxonomy]
                        krona_data.append([abundance] + taxonomy)
                
                # Create Krona HTML
                krona_file = os.path.join(pub_dir, "krona_chart.html")
                krona.create(krona_data, krona_file)
                
            except ImportError:
                self.logger.warning("Krona package not available, skipping Krona visualization")
            
            self.logger.info("Taxonomic classification and visualization completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def analyze_alpha_diversity(self):
        """Analyze alpha diversity using multiple metrics and generate visualizations."""
        step_name = "Alpha Diversity Analysis"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            from qiime2.plugins import diversity
            
            # Get feature table
            table = self.artifacts.get('table')
            if table is None:
                raise ValueError("Feature table not found. Run DADA2 denoising first.")
            
            # Get rooted tree for phylogenetic metrics
            rooted_tree = self.artifacts.get('rooted_tree')
            if rooted_tree is None:
                raise ValueError("Rooted tree not found. Run phylogenetic tree construction first.")
            
            # Create directory for alpha diversity results
            alpha_dir = os.path.join(self.output_dir, "diversity", "alpha")
            os.makedirs(alpha_dir, exist_ok=True)
            
            # Create directory for publication-quality figures
            pub_dir = os.path.join(self.output_dir, "publication_figures", "diversity")
            os.makedirs(pub_dir, exist_ok=True)
            
            # Load metadata
            metadata = qiime2.Metadata.load(self.metadata)
            
            # Calculate alpha diversity metrics
            self.logger.info("Calculating alpha diversity metrics")
            
            # Define metrics to calculate
            metrics = [
                ('observed_features', 'Observed Features'),
                ('chao1', 'Chao1'),
                ('fisher_alpha', 'Fisher Alpha'),
                ('pielou_e', 'Pielou\'s Evenness'),
                ('shannon', 'Shannon Diversity'),
                ('simpson', 'Simpson Diversity')
            ]
            
            # Store results for each metric
            alpha_results = {}
            
            # Calculate each metric
            for metric_id, metric_name in metrics:
                self.logger.info(f"Calculating {metric_name}")
                
                # Calculate alpha diversity
                alpha_div = diversity.methods.alpha_diversity(
                    table=table,
                    metric=metric_id
                )
                
                # Save artifact
                alpha_div.alpha_diversity.save(
                    os.path.join(alpha_dir, f"{metric_id}_vector.qza")
                )
                
                # Store for later use
                alpha_results[metric_id] = alpha_div.alpha_diversity
                
                # Generate statistical visualization
                alpha_vis = diversity.visualizers.alpha_group_significance(
                    alpha=alpha_div.alpha_diversity,
                    metadata=metadata
                )
                alpha_vis.visualization.save(
                    os.path.join(alpha_dir, f"{metric_id}_group_significance.qzv")
                )
            
            # Generate publication-quality visualizations
            self.logger.info("Generating publication-quality visualizations")
            
            # Set plot style
            plt.style.use('default')
            plt.rcParams.update({
                'figure.figsize': (10, 6),
                'font.size': 12,
                'font.family': 'sans-serif',
                'axes.labelsize': 14,
                'axes.titlesize': 16,
                'xtick.labelsize': 12,
                'ytick.labelsize': 12,
                'legend.fontsize': 12
            })
            
            # Load metadata as DataFrame
            metadata_df = pd.read_csv(self.metadata, sep='\t', index_col='#SampleID')
            
            # Create combined visualization for all metrics
            fig, axes = plt.subplots(2, 3, figsize=(20, 12))
            axes = axes.ravel()
            
            for idx, (metric_id, metric_name) in enumerate(metrics):
                # Get alpha diversity data
                alpha_data = alpha_results[metric_id].view(pd.Series)
                
                # Combine with metadata
                plot_data = pd.DataFrame({
                    'value': alpha_data,
                    'group': metadata_df.loc[alpha_data.index, 'group']
                })
                
                # Create violin plot with individual points
                ax = axes[idx]
                sns.violinplot(data=plot_data, x='group', y='value', ax=ax, 
                             inner='box', color='lightgray')
                sns.stripplot(data=plot_data, x='group', y='value', ax=ax,
                            size=6, alpha=0.6, jitter=True)
                
                # Add statistical annotation
                groups = plot_data['group'].unique()
                if len(groups) > 1:
                    stat_result = stats.kruskal(*[plot_data[plot_data['group'] == g]['value'] 
                                                for g in groups])
                    ax.text(0.95, 0.95, f'p = {stat_result.pvalue:.2e}',
                           transform=ax.transAxes, ha='right', va='top')
                
                # Customize plot
                ax.set_title(metric_name)
                ax.set_xlabel('Group')
                ax.set_ylabel(metric_name)
                ax.tick_params(axis='x', rotation=45)
                ax.grid(True, linestyle='--', alpha=0.7)
                
                # Save individual plot
                plt.figure(figsize=(10, 6))
                sns.violinplot(data=plot_data, x='group', y='value', inner='box', color='lightgray')
                sns.stripplot(data=plot_data, x='group', y='value', size=6, alpha=0.6, jitter=True)
                
                if len(groups) > 1:
                    plt.text(0.95, 0.95, f'p = {stat_result.pvalue:.2e}',
                           transform=plt.gca().transAxes, ha='right', va='top')
                
                plt.title(metric_name, pad=20)
                plt.xlabel('Group', labelpad=10)
                plt.ylabel(metric_name, labelpad=10)
                plt.xticks(rotation=45)
                plt.grid(True, linestyle='--', alpha=0.7)
                plt.tight_layout()
                plt.savefig(
                    os.path.join(pub_dir, f"{metric_id}_boxplot.svg"),
                    bbox_inches='tight',
                    dpi=300
                )
                plt.close()
            
            # Adjust layout of combined plot
            fig.suptitle('Alpha Diversity Metrics', y=1.02, fontsize=18)
            plt.tight_layout()
            plt.savefig(
                os.path.join(pub_dir, "alpha_diversity_combined.svg"),
                bbox_inches='tight',
                dpi=300
            )
            plt.close()
            
            # Export alpha diversity data
            alpha_dir = os.path.join(self.output_dir, "exported_data", "alpha_diversity")
            os.makedirs(alpha_dir, exist_ok=True)
            
            # Combine all metrics into one DataFrame
            combined_data = pd.DataFrame()
            for metric_id, metric_name in metrics:
                combined_data[metric_name] = alpha_results[metric_id].view(pd.Series)
            
            # Add metadata
            combined_data = pd.merge(
                combined_data,
                metadata_df[['group']],
                left_index=True,
                right_index=True
            )
            
            # Save combined data
            combined_data.to_csv(
                os.path.join(alpha_dir, "alpha_diversity_metrics.tsv"),
                sep='\t'
            )
            
            # Calculate and save summary statistics
            summary_stats = combined_data.groupby('group').agg(['mean', 'std', 'min', 'max'])
            summary_stats.to_csv(
                os.path.join(alpha_dir, "alpha_diversity_summary.tsv"),
                sep='\t'
            )
            
            self.logger.info("Alpha diversity analysis completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def analyze_beta_diversity(self):
        """Analyze beta diversity using multiple metrics and generate visualizations."""
        step_name = "Beta Diversity Analysis"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            from qiime2.plugins import diversity
            
            # Get necessary artifacts
            table = self.artifacts.get('table')
            if table is None:
                raise ValueError("Feature table not found. Run DADA2 denoising first.")
            
            rooted_tree = self.artifacts.get('rooted_tree')
            if rooted_tree is None:
                raise ValueError("Rooted tree not found. Run phylogenetic tree construction first.")
            
            # Create directories
            beta_dir = os.path.join(self.output_dir, "diversity", "beta")
            pub_dir = os.path.join(self.output_dir, "publication_figures", "diversity")
            os.makedirs(beta_dir, exist_ok=True)
            os.makedirs(pub_dir, exist_ok=True)
            
            # Load metadata
            metadata = qiime2.Metadata.load(self.metadata)
            metadata_df = metadata.to_dataframe()
            
            # Define distance metrics
            metrics = [
                ('braycurtis', 'Bray-Curtis'),
                ('jaccard', 'Jaccard'),
                ('unweighted_unifrac', 'Unweighted UniFrac'),
                ('weighted_unifrac', 'Weighted UniFrac')
            ]
            
            # Store distance matrices and ordination results
            distance_matrices = {}
            pcoa_results = {}
            
            # Calculate distance matrices and ordinations
            for metric_id, metric_name in metrics:
                self.logger.info(f"Processing {metric_name} distances")
                
                # Calculate distance matrix
                if metric_id in ['unweighted_unifrac', 'weighted_unifrac']:
                    dm = diversity.methods.beta_phylogenetic(
                        table=table,
                        phylogeny=rooted_tree,
                        metric=metric_id
                    )
                    distance_matrix = dm.distance_matrix
                else:
                    dm = diversity.methods.beta(
                        table=table,
                        metric=metric_id
                    )
                    distance_matrix = dm.distance_matrix
                
                # Save distance matrix
                distance_matrix.save(
                    os.path.join(beta_dir, f"{metric_id}_distance_matrix.qza")
                )
                
                # Store for later use
                distance_matrices[metric_id] = distance_matrix
                
                # Perform statistical tests
                self.logger.info(f"Performing statistical tests for {metric_name}")
                
                # Convert to skbio DistanceMatrix
                dm_array = distance_matrix.view(skbio.DistanceMatrix)
                
                # PERMANOVA
                permanova_results = permanova(
                    dm_array,
                    metadata_df['group'],
                    permutations=999
                )
                
                # ANOSIM
                anosim_results = anosim(
                    dm_array,
                    metadata_df['group'],
                    permutations=999
                )
                
                # Save statistical results
                stats_file = os.path.join(beta_dir, f"{metric_id}_statistical_tests.txt")
                with open(stats_file, 'w') as f:
                    f.write(f"PERMANOVA Results for {metric_name}:\n")
                    f.write(f"R-squared: {permanova_results['test-statistic']:.4f}\n")
                    f.write(f"p-value: {permanova_results['p-value']:.4f}\n\n")
                    
                    f.write(f"ANOSIM Results for {metric_name}:\n")
                    f.write(f"R-value: {anosim_results['test-statistic']:.4f}\n")
                    f.write(f"p-value: {anosim_results['p-value']:.4f}\n")
                
                # Perform PCoA
                pcoa_result = pcoa(dm_array)
                pcoa_results[metric_id] = pcoa_result
                
                # Generate PCoA plot
                self.logger.info(f"Generating PCoA plot for {metric_name}")
                
                # Create PCoA plot
                fig, ax = plt.subplots(figsize=(10, 8))
                
                # Get coordinates
                coords = pcoa_result.samples.copy()
                
                # Add group information
                coords['group'] = metadata_df.loc[coords.index, 'group']
                
                # Create scatter plot
                groups = coords['group'].unique()
                colors = plt.cm.Set3(np.linspace(0, 1, len(groups)))
                
                for i, group in enumerate(groups):
                    mask = coords['group'] == group
                    ax.scatter(
                        coords.loc[mask, 0],
                        coords.loc[mask, 1],
                        c=[colors[i]],
                        label=group,
                        s=100,
                        alpha=0.7
                    )
                
                # Add explained variance
                explained_var = pcoa_result.proportion_explained
                ax.set_xlabel(f'PC1 ({explained_var[0]:.1%} explained)')
                ax.set_ylabel(f'PC2 ({explained_var[1]:.1%} explained)')
                
                # Add statistical results
                ax.text(0.05, 0.95,
                       f'PERMANOVA p = {permanova_results["p-value"]:.2e}\n'
                       f'ANOSIM p = {anosim_results["p-value"]:.2e}',
                       transform=ax.transAxes,
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                # Customize plot
                ax.set_title(f'{metric_name} PCoA')
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                ax.grid(True, linestyle='--', alpha=0.7)
                
                # Save plot
                plt.tight_layout()
                plt.savefig(
                    os.path.join(pub_dir, f"{metric_id}_pcoa.svg"),
                    bbox_inches='tight',
                    dpi=300
                )
                plt.close()
                
                # Generate NMDS plot if available
                if nmds_results[metric_id] is not None:
                    self.logger.info(f"Generating NMDS plot for {metric_name}")
                    
                    fig, ax = plt.subplots(figsize=(10, 8))
                    
                    # Get coordinates
                    nmds_coords = pd.DataFrame(
                        nmds_results[metric_id].samples.copy(),
                        columns=['NMDS1', 'NMDS2'],
                        index=metadata_df.index
                    )
                    
                    # Add group information
                    nmds_coords['group'] = metadata_df['group']
                    
                    # Create scatter plot
                    for i, group in enumerate(groups):
                        mask = nmds_coords['group'] == group
                        ax.scatter(
                            nmds_coords.loc[mask, 'NMDS1'],
                            nmds_coords.loc[mask, 'NMDS2'],
                            c=[colors[i]],
                            label=group,
                            s=100,
                            alpha=0.7
                        )
                    
                    # Add stress value
                    stress = nmds_results[metric_id].stress
                    ax.text(0.05, 0.95,
                           f'Stress = {stress:.3f}',
                           transform=ax.transAxes,
                           verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                    
                    # Customize plot
                    ax.set_title(f'{metric_name} NMDS')
                    ax.set_xlabel('NMDS1')
                    ax.set_ylabel('NMDS2')
                    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                    ax.grid(True, linestyle='--', alpha=0.7)
                    
                    # Save plot
                    plt.tight_layout()
                    plt.savefig(
                        os.path.join(pub_dir, f"{metric_id}_nmds.svg"),
                        bbox_inches='tight',
                        dpi=300
                    )
                    plt.close()
                    
                    # Generate stress plot
                    fig, ax = plt.subplots(figsize=(8, 8))
                    observed_distances = dm_array.condensed_form()
                    nmds_distances = pdist(nmds_results[metric_id].samples)
                    
                    ax.scatter(observed_distances, nmds_distances, alpha=0.5)
                    ax.plot([min(observed_distances), max(observed_distances)],
                           [min(observed_distances), max(observed_distances)],
                           'r--', label='1:1 line')
                    
                    ax.set_xlabel('Observed Distances')
                    ax.set_ylabel('Ordinated Distances')
                    ax.set_title(f'{metric_name} NMDS Stress Plot')
                    ax.legend()
                    
                    plt.tight_layout()
                    plt.savefig(
                        os.path.join(pub_dir, f"{metric_id}_stress_plot.svg"),
                        bbox_inches='tight',
                        dpi=300
                    )
                    plt.close()
            
            # Perform PCA on normalized abundance data
            self.logger.info("Performing PCA analysis")
            
            # Get feature table as DataFrame
            table_df = table.view(pd.DataFrame)
            
            # Normalize and scale data
            rel_abundance = table_df.div(table_df.sum(axis=0), axis=1) * 100
            scaled_data = StandardScaler().fit_transform(rel_abundance.T)
            
            # Perform PCA
            pca = PCA()
            pca_result = pca.fit_transform(scaled_data)
            
            # Create PCA plot
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Create DataFrame with PCA results
            pca_df = pd.DataFrame(
                pca_result[:, :2],
                columns=['PC1', 'PC2'],
                index=metadata_df.index
            )
            
            # Add group information
            pca_df['group'] = metadata_df['group']
            
            # Plot PCA
            for i, group in enumerate(groups):
                mask = pca_df['group'] == group
                ax.scatter(
                    pca_df.loc[mask, 'PC1'],
                    pca_df.loc[mask, 'PC2'],
                    c=[colors[i]],
                    label=group,
                    s=100,
                    alpha=0.7
                )
            
            # Add explained variance
            ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} explained)')
            ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} explained)')
            
            # Customize plot
            ax.set_title('PCA of Normalized Abundance Data')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(True, linestyle='--', alpha=0.7)
            
            # Save plot
            plt.tight_layout()
            plt.savefig(
                os.path.join(pub_dir, "abundance_pca.svg"),
                bbox_inches='tight',
                dpi=300
            )
            plt.close()
            
            # Export ordination results
            ord_dir = os.path.join(self.output_dir, "exported_data", "ordination")
            os.makedirs(ord_dir, exist_ok=True)
            
            # Save PCA results
            pca_df.to_csv(os.path.join(ord_dir, "pca_coordinates.tsv"), sep='\t')
            
            # Save coordinates for each metric
            for metric_id, metric_name in metrics:
                # Save PCoA coordinates
                pcoa_coords = pd.DataFrame(
                    pcoa_results[metric_id].samples,
                    columns=[f'PC{i+1}' for i in range(pcoa_results[metric_id].samples.shape[1])]
                )
                pcoa_coords['group'] = metadata_df.loc[pcoa_coords.index, 'group']
                pcoa_coords.to_csv(
                    os.path.join(ord_dir, f"{metric_id}_pcoa_coordinates.tsv"),
                    sep='\t'
                )
                
                # Save NMDS coordinates if available
                if nmds_results[metric_id] is not None:
                    nmds_coords = pd.DataFrame(
                        nmds_results[metric_id].samples,
                        columns=['NMDS1', 'NMDS2']
                    )
                    nmds_coords['group'] = metadata_df.loc[nmds_coords.index, 'group']
                    nmds_coords.to_csv(
                        os.path.join(ord_dir, f"{metric_id}_nmds_coordinates.tsv"),
                        sep='\t'
                    )
            
            self.logger.info("Beta diversity analysis completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def run_lefse_analysis(self):
        """Run LEfSe analysis to identify biomarkers between groups."""
        step_name = "LEfSe Analysis"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            # Get necessary artifacts
            table = self.artifacts.get('table')
            taxonomy = self.artifacts.get('taxonomy')
            
            if table is None or taxonomy is None:
                raise ValueError("Feature table and taxonomy artifacts required for LEfSe analysis")
            
            # Create LEfSe directory
            lefse_dir = os.path.join(self.output_dir, "lefse")
            os.makedirs(lefse_dir, exist_ok=True)
            
            # Load feature table and taxonomy
            table_df = table.view(pd.DataFrame)
            tax_df = taxonomy.view(pd.DataFrame)
            
            # Load metadata
            metadata_df = pd.read_csv(self.metadata, sep='\t', index_col='#SampleID')
            
            # Calculate relative abundance
            rel_abundance = table_df.div(table_df.sum(axis=0), axis=1) * 100
            
            # Process taxonomy data
            tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            tax_split = tax_df['Taxon'].str.split(';', expand=True)
            tax_split.columns = tax_levels
            
            # Clean taxonomy names
            for col in tax_levels:
                tax_split[col] = tax_split[col].str.replace(r'^[a-z]__', '', regex=True)
                tax_split[col] = tax_split[col].fillna('Unassigned')
            
            # Create LEfSe input files for each taxonomic level
            for level in tax_levels:
                self.logger.info(f"Processing {level} level for LEfSe analysis")
                
                # Aggregate abundances at the current taxonomic level
                level_abundance = pd.DataFrame(index=rel_abundance.columns)
                
                for taxon in tax_split[level].unique():
                    if taxon == 'Unassigned':
                        continue
                    mask = tax_split[level] == taxon
                    if mask.any():
                        features = tax_split[mask].index
                        level_abundance[taxon] = rel_abundance.loc[features].sum()
                
                # Prepare LEfSe input file
                lefse_input = level_abundance.T
                
                # Add class information
                lefse_input['class'] = metadata_df.loc[lefse_input.index, 'group']
                
                # Save LEfSe input
                input_file = os.path.join(lefse_dir, f"level{tax_levels.index(level)+1}_lefse_input.txt")
                lefse_input.to_csv(input_file, sep='\t')
                
                # Run LEfSe analysis using command line tools
                self.logger.info(f"Running LEfSe for {level} level")
                
                # Format input file
                format_cmd = f"lefse-format_input.py {input_file} {input_file}.in -c 1"
                subprocess.run(format_cmd, shell=True, check=True)
                
                # Run LEfSe
                run_cmd = f"run_lefse.py {input_file}.in {input_file}.res"
                subprocess.run(run_cmd, shell=True, check=True)
                
                # Generate plots
                # Plot LDA scores
                plot_cmd = f"lefse-plot_res.py {input_file}.res {lefse_dir}/level{tax_levels.index(level)+1}_lefse_lda.png"
                subprocess.run(plot_cmd, shell=True, check=True)
                
                # Generate cladogram
                clado_cmd = f"lefse-plot_cladogram.py {input_file}.res {lefse_dir}/level{tax_levels.index(level)+1}_lefse_cladogram.png --format png"
                subprocess.run(clado_cmd, shell=True, check=True)
                
                # Also save as SVG for publication quality
                clado_svg_cmd = f"lefse-plot_cladogram.py {input_file}.res {lefse_dir}/level{tax_levels.index(level)+1}_lefse_cladogram.svg --format svg"
                subprocess.run(clado_svg_cmd, shell=True, check=True)
                
                # Read and parse LEfSe results
                results_file = f"{input_file}.res"
                if os.path.exists(results_file):
                    results_df = pd.read_csv(results_file, sep='\t', header=None)
                    results_df.columns = ['Feature', 'LogMaxMean', 'Class', 'LDA', 'pValue']
                    
                    # Save formatted results
                    results_df.to_csv(
                        os.path.join(lefse_dir, f"level{tax_levels.index(level)+1}_lefse_results.tsv"),
                        sep='\t',
                        index=False
                    )
                    
                    # Create custom visualization of LEfSe results
                    if not results_df.empty:
                        plt.figure(figsize=(12, 8))
                        
                        # Sort by LDA score
                        results_df = results_df.sort_values('LDA', ascending=True)
                        
                        # Create color map for different classes
                        unique_classes = results_df['Class'].unique()
                        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_classes)))
                        class_colors = dict(zip(unique_classes, colors))
                        
                        # Create horizontal bar plot
                        bars = plt.barh(
                            range(len(results_df)),
                            results_df['LDA'],
                            color=[class_colors[c] for c in results_df['Class']]
                        )
                        
                        # Customize plot
                        plt.yticks(
                            range(len(results_df)),
                            results_df['Feature'],
                            fontsize=10
                        )
                        plt.xlabel('LDA Score (log10)')
                        plt.title(f'LEfSe Results - {level} Level')
                        
                        # Add legend
                        legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color)
                                         for color in colors]
                        plt.legend(legend_elements, unique_classes,
                                 loc='center left', bbox_to_anchor=(1, 0.5))
                        
                        # Save plot
                        plt.tight_layout()
                        plt.savefig(
                            os.path.join(lefse_dir, f"level{tax_levels.index(level)+1}_lefse_custom.svg"),
                            bbox_inches='tight',
                            dpi=300
                        )
                        plt.close()
            
            self.logger.info("LEfSe analysis completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def run_random_forest_analysis(self):
        """Run Random Forest analysis to identify important features and create classification models."""
        step_name = "Random Forest Analysis"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            # Get necessary artifacts
            table = self.artifacts.get('table')
            taxonomy = self.artifacts.get('taxonomy')
            
            if table is None or taxonomy is None:
                raise ValueError("Feature table and taxonomy artifacts required for Random Forest analysis")
            
            # Create Random Forest directory
            rf_dir = os.path.join(self.output_dir, "random_forest")
            os.makedirs(rf_dir, exist_ok=True)
            
            # Load feature table and taxonomy
            table_df = table.view(pd.DataFrame)
            tax_df = taxonomy.view(pd.DataFrame)
            
            # Load metadata
            metadata_df = pd.read_csv(self.metadata, sep='\t', index_col='#SampleID')
            
            # Calculate relative abundance
            rel_abundance = table_df.div(table_df.sum(axis=0), axis=1) * 100
            
            # Process taxonomy data
            tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            tax_split = tax_df['Taxon'].str.split(';', expand=True)
            tax_split.columns = tax_levels
            
            # Clean taxonomy names
            for col in tax_levels:
                tax_split[col] = tax_split[col].str.replace(r'^[a-z]__', '', regex=True)
                tax_split[col] = tax_split[col].fillna('Unassigned')
            
            # Perform Random Forest analysis for each taxonomic level
            for level in tax_levels:
                self.logger.info(f"Running Random Forest analysis for {level} level")
                
                # Aggregate abundances at the current taxonomic level
                level_abundance = pd.DataFrame(index=rel_abundance.columns)
                
                for taxon in tax_split[level].unique():
                    if taxon == 'Unassigned':
                        continue
                    mask = tax_split[level] == taxon
                    if mask.any():
                        features = tax_split[mask].index
                        level_abundance[taxon] = rel_abundance.loc[features].sum()
                
                # Prepare data for Random Forest
                X = level_abundance.values
                y = metadata_df.loc[level_abundance.index, 'group']
                
                # Encode labels
                le = LabelEncoder()
                y_encoded = le.fit_transform(y)
                
                # Initialize Random Forest classifier
                rf = RandomForestClassifier(
                    n_estimators=1000,
                    max_features='sqrt',
                    n_jobs=self.threads,
                    random_state=42
                )
                
                # Perform cross-validation
                cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
                cv_scores = cross_val_score(rf, X, y_encoded, cv=cv, scoring='accuracy')
                
                # Fit the model on the entire dataset
                rf.fit(X, y_encoded)
                
                # Get feature importance
                importance = pd.DataFrame({
                    'Feature': level_abundance.columns,
                    'Importance': rf.feature_importances_
                })
                importance = importance.sort_values('Importance', ascending=False)
                
                # Save feature importance
                importance.to_csv(
                    os.path.join(rf_dir, f"{level.lower()}_feature_importance.tsv"),
                    sep='\t',
                    index=False
                )
                
                # Create feature importance plot
                plt.figure(figsize=(12, 8))
                sns.barplot(
                    data=importance.head(20),
                    x='Importance',
                    y='Feature',
                    palette='viridis'
                )
                plt.title(f'Top 20 Important Features - {level} Level')
                plt.xlabel('Feature Importance')
                plt.tight_layout()
                plt.savefig(
                    os.path.join(rf_dir, f"{level.lower()}_feature_importance.svg"),
                    bbox_inches='tight',
                    dpi=300
                )
                plt.close()
                
                # Generate confusion matrix
                y_pred = rf.predict(X)
                cm = confusion_matrix(y_encoded, y_pred)
                
                # Plot confusion matrix
                plt.figure(figsize=(8, 6))
                sns.heatmap(
                    cm,
                    annot=True,
                    fmt='d',
                    cmap='Blues',
                    xticklabels=le.classes_,
                    yticklabels=le.classes_
                )
                plt.title(f'Confusion Matrix - {level} Level')
                plt.xlabel('Predicted')
                plt.ylabel('True')
                plt.tight_layout()
                plt.savefig(
                    os.path.join(rf_dir, f"{level.lower()}_confusion_matrix.svg"),
                    bbox_inches='tight',
                    dpi=300
                )
                plt.close()
                
                # Generate ROC curves for each class
                if len(le.classes_) == 2:
                    # Binary classification
                    y_prob = rf.predict_proba(X)[:, 1]
                    fpr, tpr, _ = roc_curve(y_encoded, y_prob)
                    roc_auc = auc(fpr, tpr)
                    
                    plt.figure(figsize=(8, 6))
                    plt.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc:.2f})')
                    plt.plot([0, 1], [0, 1], 'k--')
                    plt.xlim([0.0, 1.0])
                    plt.ylim([0.0, 1.05])
                    plt.xlabel('False Positive Rate')
                    plt.ylabel('True Positive Rate')
                    plt.title(f'ROC Curve - {level} Level')
                    plt.legend(loc="lower right")
                    plt.tight_layout()
                    plt.savefig(
                        os.path.join(rf_dir, f"{level.lower()}_roc_curve.svg"),
                        bbox_inches='tight',
                        dpi=300
                    )
                    plt.close()
                
                # Save classification report
                report = classification_report(y_encoded, y_pred, target_names=le.classes_)
                with open(os.path.join(rf_dir, f"{level.lower()}_classification_report.txt"), 'w') as f:
                    f.write(f"Random Forest Classification Report - {level} Level\n\n")
                    f.write(f"Cross-validation scores: {cv_scores}\n")
                    f.write(f"Mean CV accuracy: {cv_scores.mean():.3f} (+/- {cv_scores.std() * 2:.3f})\n\n")
                    f.write(report)
            
            self.logger.info("Random Forest analysis completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def run_network_analysis(self):
        """Run network analysis to identify microbial interactions."""
        step_name = "Network Analysis"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            # Get necessary artifacts
            table = self.artifacts.get('table')
            taxonomy = self.artifacts.get('taxonomy')
            
            if table is None or taxonomy is None:
                raise ValueError("Feature table and taxonomy artifacts required for network analysis")
            
            # Create network directory
            network_dir = os.path.join(self.output_dir, "network")
            os.makedirs(network_dir, exist_ok=True)
            
            # Load feature table and taxonomy
            table_df = table.view(pd.DataFrame)
            tax_df = taxonomy.view(pd.DataFrame)
            
            # Calculate relative abundance
            rel_abundance = table_df.div(table_df.sum(axis=0), axis=1) * 100
            
            # Process taxonomy data
            tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            tax_split = tax_df['Taxon'].str.split(';', expand=True)
            tax_split.columns = tax_levels
            
            # Clean taxonomy names
            for col in tax_levels:
                tax_split[col] = tax_split[col].str.replace(r'^[a-z]__', '', regex=True)
                tax_split[col] = tax_split[col].fillna('Unassigned')
            
            # Perform network analysis at genus level
            self.logger.info("Performing network analysis at genus level")
            
            # Aggregate abundances at genus level
            genus_abundance = pd.DataFrame(index=rel_abundance.columns)
            for taxon in tax_split['Genus'].unique():
                if taxon == 'Unassigned':
                    continue
                mask = tax_split['Genus'] == taxon
                if mask.any():
                    features = tax_split[mask].index
                    genus_abundance[taxon] = rel_abundance.loc[features].sum()
            
            # Calculate Spearman correlations
            self.logger.info("Calculating Spearman correlations")
            corr_matrix = pd.DataFrame(index=genus_abundance.columns, columns=genus_abundance.columns)
            pval_matrix = pd.DataFrame(index=genus_abundance.columns, columns=genus_abundance.columns)
            
            for i in genus_abundance.columns:
                for j in genus_abundance.columns:
                    if i != j:
                        corr, pval = spearmanr(genus_abundance[i], genus_abundance[j])
                        corr_matrix.loc[i, j] = corr
                        pval_matrix.loc[i, j] = pval
            
            # Multiple testing correction
            pvals_flat = pval_matrix.values[~np.eye(pval_matrix.shape[0], dtype=bool)]
            _, pvals_adj, _, _ = multipletests(pvals_flat, method='fdr_bh')
            
            # Reshape adjusted p-values back to matrix
            pval_matrix_adj = pd.DataFrame(
                np.eye(pval_matrix.shape[0]),
                index=pval_matrix.index,
                columns=pval_matrix.columns
            )
            pval_matrix_adj.values[~np.eye(pval_matrix.shape[0], dtype=bool)] = pvals_adj
            
            # Create network
            self.logger.info("Creating correlation network")
            G = nx.Graph()
            
            # Add nodes (genera)
            for genus in genus_abundance.columns:
                G.add_node(genus)
            
            # Add edges (significant correlations)
            for i in genus_abundance.columns:
                for j in genus_abundance.columns:
                    if i < j:  # Avoid duplicate edges
                        if pval_matrix_adj.loc[i, j] < 0.05:  # Significant correlation
                            corr = corr_matrix.loc[i, j]
                            if abs(corr) >= 0.3:  # Correlation strength threshold
                                G.add_edge(i, j, weight=corr)
            
            # Calculate network metrics
            self.logger.info("Calculating network metrics")
            metrics = {
                'Number of nodes': G.number_of_nodes(),
                'Number of edges': G.number_of_edges(),
                'Average degree': np.mean([d for n, d in G.degree()]),
                'Network density': nx.density(G),
                'Average clustering coefficient': nx.average_clustering(G),
                'Average path length': nx.average_shortest_path_length(G) if nx.is_connected(G) else 'NA'
            }
            
            # Save network metrics
            with open(os.path.join(network_dir, "network_metrics.txt"), 'w') as f:
                for metric, value in metrics.items():
                    f.write(f"{metric}: {value}\n")
            
            # Save correlation and p-value matrices
            corr_matrix.to_csv(os.path.join(network_dir, "correlation_matrix.tsv"), sep='\t')
            pval_matrix_adj.to_csv(os.path.join(network_dir, "adjusted_pvalue_matrix.tsv"), sep='\t')
            
            # Generate network visualization
            self.logger.info("Generating network visualization")
            
            # Calculate node sizes based on mean abundance
            node_sizes = {genus: np.mean(genus_abundance[genus]) for genus in G.nodes()}
            
            # Calculate layout
            pos = nx.spring_layout(G, k=1/np.sqrt(G.number_of_nodes()), iterations=50)
            
            # Create figure
            plt.figure(figsize=(15, 15))
            
            # Draw nodes
            nx.draw_networkx_nodes(
                G, pos,
                node_size=[node_sizes[n]*100 for n in G.nodes()],
                node_color='lightblue',
                alpha=0.7
            )
            
            # Draw edges with different colors for positive and negative correlations
            edges_pos = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > 0]
            edges_neg = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] < 0]
            
            nx.draw_networkx_edges(
                G, pos,
                edgelist=edges_pos,
                edge_color='green',
                alpha=0.5
            )
            nx.draw_networkx_edges(
                G, pos,
                edgelist=edges_neg,
                edge_color='red',
                alpha=0.5
            )
            
            # Add labels
            nx.draw_networkx_labels(
                G, pos,
                font_size=8,
                font_weight='bold'
            )
            
            plt.title("Microbial Co-occurrence Network\n(Green: positive correlation, Red: negative correlation)")
            plt.axis('off')
            
            # Save network visualization
            plt.savefig(
                os.path.join(network_dir, "correlation_network.svg"),
                bbox_inches='tight',
                dpi=300
            )
            plt.close()
            
            # Generate heatmap visualization
            plt.figure(figsize=(12, 10))
            mask = np.triu(np.ones_like(corr_matrix))
            sns.heatmap(
                corr_matrix,
                mask=mask,
                cmap='RdBu_r',
                center=0,
                vmin=-1,
                vmax=1,
                square=True,
                annot=False,
                fmt='.2f'
            )
            plt.title('Correlation Heatmap')
            plt.tight_layout()
            plt.savefig(
                os.path.join(network_dir, "correlation_heatmap.svg"),
                bbox_inches='tight',
                dpi=300
            )
            plt.close()
            
            # Export network in GraphML format for external visualization
            nx.write_graphml(G, os.path.join(network_dir, "network.graphml"))
            
            self.logger.info("Network analysis completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def analyze_relative_abundance(self):
        """Analyze and generate detailed relative abundance tables and visualizations."""
        step_name = "Relative Abundance Analysis"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            # Get necessary artifacts
            table = self.artifacts.get('table')
            taxonomy = self.artifacts.get('taxonomy')
            
            if table is None or taxonomy is None:
                raise ValueError("Feature table and taxonomy artifacts required for relative abundance analysis")
            
            # Create relative abundance directory
            abundance_dir = os.path.join(self.output_dir, "abundance")
            os.makedirs(abundance_dir, exist_ok=True)
            
            # Create publication figures directory
            pub_dir = os.path.join(self.output_dir, "publication_figures", "abundance")
            os.makedirs(pub_dir, exist_ok=True)
            
            # Load feature table and taxonomy
            table_df = table.view(pd.DataFrame)
            tax_df = taxonomy.view(pd.DataFrame)
            
            # Load metadata
            metadata_df = pd.read_csv(self.metadata, sep='\t', index_col='#SampleID')
            
            # Calculate relative abundance
            rel_abundance = table_df.div(table_df.sum(axis=0), axis=1) * 100
            
            # Process taxonomy data
            tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            tax_split = tax_df['Taxon'].str.split(';', expand=True)
            tax_split.columns = tax_levels
            
            # Clean taxonomy names
            for col in tax_levels:
                tax_split[col] = tax_split[col].str.replace(r'^[a-z]__', '', regex=True)
                tax_split[col] = tax_split[col].fillna('Unassigned')
            
            # Create relative abundance tables for each taxonomic level
            for level in tax_levels:
                self.logger.info(f"Processing {level} level relative abundance")
                
                # Aggregate abundances at the current taxonomic level
                level_abundance = pd.DataFrame(index=rel_abundance.columns)
                
                for taxon in tax_split[level].unique():
                    if taxon == 'Unassigned':
                        continue
                    mask = tax_split[level] == taxon
                    if mask.any():
                        features = tax_split[mask].index
                        level_abundance[taxon] = rel_abundance.loc[features].sum()
                
                # Save full abundance table
                level_abundance.to_csv(
                    os.path.join(abundance_dir, f"{level.lower()}_abundance_table.tsv"),
                    sep='\t'
                )
                
                # Calculate summary statistics
                summary_stats = pd.DataFrame({
                    'Mean': level_abundance.mean(),
                    'Std': level_abundance.std(),
                    'Min': level_abundance.min(),
                    'Max': level_abundance.max(),
                    'Prevalence': (level_abundance > 0).mean() * 100
                })
                
                # Save summary statistics
                summary_stats.to_csv(
                    os.path.join(abundance_dir, f"{level.lower()}_abundance_summary.tsv"),
                    sep='\t'
                )
                
                # Calculate group-wise statistics if groups exist
                if 'group' in metadata_df.columns:
                    group_stats = pd.DataFrame()
                    for group in metadata_df['group'].unique():
                        group_samples = metadata_df[metadata_df['group'] == group].index
                        group_abundance = level_abundance.loc[group_samples]
                        
                        group_stats[f'{group}_Mean'] = group_abundance.mean()
                        group_stats[f'{group}_Std'] = group_abundance.std()
                        group_stats[f'{group}_Prevalence'] = (group_abundance > 0).mean() * 100
                    
                    # Save group-wise statistics
                    group_stats.to_csv(
                        os.path.join(abundance_dir, f"{level.lower()}_abundance_by_group.tsv"),
                        sep='\t'
                    )
                
                # Generate visualizations
                # Top 20 taxa abundance boxplot
                plt.figure(figsize=(15, 10))
                top_taxa = level_abundance.mean().sort_values(ascending=False).head(20).index
                
                # Create boxplot data
                plot_data = pd.melt(
                    level_abundance[top_taxa],
                    var_name=level,
                    value_name='Relative Abundance (%)'
                )
                
                # Create boxplot
                sns.boxplot(
                    data=plot_data,
                    x=level,
                    y='Relative Abundance (%)',
                    color='lightblue'
                )
                sns.stripplot(
                    data=plot_data,
                    x=level,
                    y='Relative Abundance (%)',
                    color='black',
                    size=4,
                    alpha=0.3
                )
                
                plt.xticks(rotation=45, ha='right')
                plt.title(f'Top 20 Most Abundant {level}s')
                plt.tight_layout()
                plt.savefig(
                    os.path.join(pub_dir, f"{level.lower()}_top20_boxplot.svg"),
                    bbox_inches='tight',
                    dpi=300
                )
                plt.close()
                
                # Generate group-wise boxplots if groups exist
                if 'group' in metadata_df.columns:
                    plt.figure(figsize=(15, 10))
                    plot_data['Group'] = metadata_df.loc[level_abundance.index, 'group'].values
                    
                    sns.boxplot(
                        data=plot_data,
                        x=level,
                        y='Relative Abundance (%)',
                        hue='Group',
                        palette='Set3'
                    )
                    
                    plt.xticks(rotation=45, ha='right')
                    plt.title(f'Top 20 Most Abundant {level}s by Group')
                    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                    plt.tight_layout()
                    plt.savefig(
                        os.path.join(pub_dir, f"{level.lower()}_top20_by_group_boxplot.svg"),
                        bbox_inches='tight',
                        dpi=300
                    )
                    plt.close()
                
                # Generate abundance heatmap
                plt.figure(figsize=(15, 10))
                sns.heatmap(
                    level_abundance[top_taxa].T,
                    cmap='YlOrRd',
                    robust=True,
                    cbar_kws={'label': 'Relative Abundance (%)'}
                )
                plt.title(f'Top 20 Most Abundant {level}s Heatmap')
                plt.tight_layout()
                plt.savefig(
                    os.path.join(pub_dir, f"{level.lower()}_top20_heatmap.svg"),
                    bbox_inches='tight',
                    dpi=300
                )
                plt.close()
            
            # Generate overall abundance summary
            self.logger.info("Generating overall abundance summary")
            
            # Calculate total abundance for each taxonomic level
            total_abundance = pd.DataFrame()
            for level in tax_levels:
                total_abundance[level] = tax_split[level].value_counts()
            
            # Save total abundance summary
            total_abundance.to_csv(
                os.path.join(abundance_dir, "total_abundance_summary.tsv"),
                sep='\t'
            )
            
            # Create summary visualization
            plt.figure(figsize=(12, 6))
            total_abundance.sum().plot(kind='bar')
            plt.title('Total Number of Taxa at Each Taxonomic Level')
            plt.xlabel('Taxonomic Level')
            plt.ylabel('Number of Taxa')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(
                os.path.join(pub_dir, "taxonomic_level_summary.svg"),
                bbox_inches='tight',
                dpi=300
            )
            plt.close()
            
            self.logger.info("Relative abundance analysis completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

    def run_deseq2_analysis(self):
        """Run DESeq2 differential abundance analysis."""
        step_name = "DESeq2 Analysis"
        
        if self.is_step_completed(step_name):
            self.current_step += 1
            return
            
        self.logger.info(f"Starting {step_name}")
        
        try:
            # Get necessary artifacts
            table = self.artifacts.get('table')
            taxonomy = self.artifacts.get('taxonomy')
            
            if table is None or taxonomy is None:
                raise ValueError("Feature table and taxonomy artifacts required for DESeq2 analysis")
            
            # Create DESeq2 directory
            deseq2_dir = os.path.join(self.output_dir, "deseq2")
            os.makedirs(deseq2_dir, exist_ok=True)
            
            # Create publication figures directory
            pub_dir = os.path.join(self.output_dir, "publication_figures", "deseq2")
            os.makedirs(pub_dir, exist_ok=True)
            
            # Load feature table and taxonomy
            table_df = table.view(pd.DataFrame)
            tax_df = taxonomy.view(pd.DataFrame)
            
            # Load metadata
            metadata_df = pd.read_csv(self.metadata, sep='\t', index_col='#SampleID')
            
            # Process taxonomy data
            tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            tax_split = tax_df['Taxon'].str.split(';', expand=True)
            tax_split.columns = tax_levels
            
            # Clean taxonomy names
            for col in tax_levels:
                tax_split[col] = tax_split[col].str.replace(r'^[a-z]__', '', regex=True)
                tax_split[col] = tax_split[col].fillna('Unassigned')
            
            # Initialize R packages
            self.logger.info("Initializing R packages")
            pandas2ri.activate()
            deseq = importr('DESeq2')
            base = importr('base')
            stats = importr('stats')
            
            # Perform DESeq2 analysis for each taxonomic level
            for level in tax_levels:
                self.logger.info(f"Running DESeq2 analysis for {level} level")
                
                # Aggregate counts at the current taxonomic level
                level_counts = pd.DataFrame(index=table_df.columns)
                
                for taxon in tax_split[level].unique():
                    if taxon == 'Unassigned':
                        continue
                    mask = tax_split[level] == taxon
                    if mask.any():
                        features = tax_split[mask].index
                        level_counts[taxon] = table_df.loc[features].sum()
                
                # Convert to integer counts
                level_counts = level_counts.astype(int)
                
                # Create sample information DataFrame
                sample_info = metadata_df[['group']].copy()
                
                # Convert DataFrames to R objects
                count_matrix = pandas2ri.py2rpy(level_counts.T)
                coldata = pandas2ri.py2rpy(sample_info)
                
                # Create DESeqDataSet
                dds = deseq.DESeqDataSetFromMatrix(
                    countData=count_matrix,
                    colData=coldata,
                    design=robjects.Formula('~ group')
                )
                
                # Run DESeq2
                dds = deseq.DESeq(dds)
                
                # Get results
                groups = sorted(metadata_df['group'].unique())
                contrast = ('group', groups[1], groups[0])  # Compare second group to first group
                res = deseq.results(dds, contrast=contrast)
                
                # Convert results to pandas DataFrame
                res_df = pd.DataFrame(
                    base.as_data_frame(res),
                    columns=['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
                )
                
                # Add taxonomy information
                res_df.index = level_counts.columns
                
                # Save results
                res_df.to_csv(
                    os.path.join(deseq2_dir, f"{level.lower()}_deseq2_results.tsv"),
                    sep='\t'
                )
                
                # Generate visualizations
                # MA plot
                plt.figure(figsize=(10, 8))
                plt.scatter(
                    np.log10(res_df['baseMean']),
                    res_df['log2FoldChange'],
                    c='grey',
                    alpha=0.5
                )
                
                # Highlight significant features
                significant = res_df['padj'] < 0.05
                plt.scatter(
                    np.log10(res_df.loc[significant, 'baseMean']),
                    res_df.loc[significant, 'log2FoldChange'],
                    c='red',
                    alpha=0.7
                )
                
                plt.axhline(y=0, color='blue', linestyle='--')
                plt.xlabel('log10(Base Mean)')
                plt.ylabel('log2(Fold Change)')
                plt.title(f'MA Plot - {level} Level')
                
                plt.tight_layout()
                plt.savefig(
                    os.path.join(pub_dir, f"{level.lower()}_ma_plot.svg"),
                    bbox_inches='tight',
                    dpi=300
                )
                plt.close()
                
                # Volcano plot
                plt.figure(figsize=(10, 8))
                plt.scatter(
                    res_df['log2FoldChange'],
                    -np.log10(res_df['pvalue']),
                    c='grey',
                    alpha=0.5
                )
                
                # Highlight significant features
                plt.scatter(
                    res_df.loc[significant, 'log2FoldChange'],
                    -np.log10(res_df.loc[significant, 'pvalue']),
                    c='red',
                    alpha=0.7
                )
                
                plt.axvline(x=0, color='blue', linestyle='--')
                plt.axhline(y=-np.log10(0.05), color='blue', linestyle='--')
                plt.xlabel('log2(Fold Change)')
                plt.ylabel('-log10(p-value)')
                plt.title(f'Volcano Plot - {level} Level')
                
                plt.tight_layout()
                plt.savefig(
                    os.path.join(pub_dir, f"{level.lower()}_volcano_plot.svg"),
                    bbox_inches='tight',
                    dpi=300
                )
                plt.close()
                
                # Heatmap of significant features
                if significant.any():
                    significant_taxa = res_df.index[significant]
                    sig_counts = level_counts[significant_taxa]
                    
                    # Normalize counts
                    normalized_counts = sig_counts.div(sig_counts.sum(axis=1), axis=0) * 100
                    
                    # Create heatmap
                    plt.figure(figsize=(12, len(significant_taxa) * 0.3 + 2))
                    sns.heatmap(
                        normalized_counts.T,
                        cmap='YlOrRd',
                        robust=True,
                        cbar_kws={'label': 'Relative Abundance (%)'}
                    )
                    plt.title(f'Significant {level}s Heatmap')
                    plt.tight_layout()
                    plt.savefig(
                        os.path.join(pub_dir, f"{level.lower()}_significant_heatmap.svg"),
                        bbox_inches='tight',
                        dpi=300
                    )
                    plt.close()
                
                # Generate summary statistics
                summary_stats = {
                    'Total features': len(res_df),
                    'Significant features (FDR < 0.05)': significant.sum(),
                    'Upregulated': ((res_df['log2FoldChange'] > 0) & significant).sum(),
                    'Downregulated': ((res_df['log2FoldChange'] < 0) & significant).sum()
                }
                
                with open(os.path.join(deseq2_dir, f"{level.lower()}_summary.txt"), 'w') as f:
                    f.write(f"DESeq2 Analysis Summary - {level} Level\n\n")
                    for stat, value in summary_stats.items():
                        f.write(f"{stat}: {value}\n")
            
            self.logger.info("DESeq2 analysis completed")
            self.create_checkpoint(step_name)
            
        except Exception as e:
            self.logger.error(f"Error in {step_name}: {str(e)}")
            raise
            
        self.current_step += 1

def main():
    parser = argparse.ArgumentParser(description='Microbiome Analysis Pipeline')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory with fastq files')
    parser.add_argument('-m', '--metadata', required=True, help='Metadata file (TSV format)')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-c', '--classifier', required=True, help='SILVA classifier file')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads')
    parser.add_argument('-s', '--sampling_depth', type=int, default=10000, help='Rarefaction depth')
    parser.add_argument('-p', '--paired', type=str, default='true', help='Paired-end data (true/false)')
    parser.add_argument('-f', '--force', type=str, default='false', help='Force rerun all steps (true/false)')
    
    args = parser.parse_args()
    
    # Convert string boolean arguments to actual booleans
    paired = args.paired.lower() == 'true'
    force = args.force.lower() == 'true'
    
    # Initialize and run pipeline
    pipeline = MicrobiomeAnalysisPipeline(
        input_dir=args.input_dir,
        metadata=args.metadata,
        output_dir=args.output_dir,
        classifier=args.classifier,
        threads=args.threads,
        sampling_depth=args.sampling_depth,
        paired=paired,
        force=force
    )
    
    # Run sequence statistics
    pipeline.process_sequence_statistics()
    
    # Run DADA2 denoising
    pipeline.run_dada2_denoising()
    
    # Generate feature table summary
    pipeline.generate_feature_table_summary()
    
    # Build phylogenetic tree
    pipeline.build_phylogenetic_tree()
    
    # Run taxonomic classification
    pipeline.run_taxonomic_classification()
    
    # Run relative abundance analysis
    pipeline.analyze_relative_abundance()
    
    # Run alpha diversity analysis
    pipeline.analyze_alpha_diversity()
    
    # Run beta diversity analysis
    pipeline.analyze_beta_diversity()
    
    # Run LEfSe analysis
    pipeline.run_lefse_analysis()
    
    # Run Random Forest analysis
    pipeline.run_random_forest_analysis()
    
    # Run Network analysis
    pipeline.run_network_analysis()
    
    # Run DESeq2 analysis
    pipeline.run_deseq2_analysis()

if __name__ == "__main__":
    main() 