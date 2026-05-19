#!/bin/bash

# =============================================================================
# QIIME2 Metabarcoding Pipeline
# =============================================================================
# Adapted from:
# QIIME2 "Moving Pictures" tutorial:
# https://docs.qiime2.org/2019.7/tutorials/moving-pictures/
#
# Australian BioCommons workshop:
# https://zenodo.org/records/6350808
#
# This workflow processes paired-end amplicon sequencing data using QIIME2.
# Sequences obtained from this study are publicly available in the National Centre for Biotechnology # Sequence Read Archive BioProjectID: PRJNA1432314
# =============================================================================

# -----------------------------------------------------------------------------
# Activate QIIME2 environment
# -----------------------------------------------------------------------------

conda env list
conda activate qiime2-2021.8

# -----------------------------------------------------------------------------
# Project directory structure
# -----------------------------------------------------------------------------

mkdir -p ~/qiime2_project/raw_data
mkdir -p ~/qiime2_project/seqs
mkdir -p ~/qiime2_project/trim
mkdir -p ~/qiime2_project/dada2out
mkdir -p ~/qiime2_project/taxonomy
mkdir -p ~/qiime2_project/downstream
mkdir -p ~/qiime2_project/export

# -----------------------------------------------------------------------------
# STEP 1: Rename FASTQ files to Casava 1.8 format
# -----------------------------------------------------------------------------
# Example format:
# SampleID_15_L001_R1_001.fastq.gz
# SampleID_15_L001_R2_001.fastq.gz
#
# Example rename command:
# mv old_sample_name.fastq new_sample_name.fastq
#
# After renaming, compress FASTQ files:
# gzip *.fastq

cd ~/qiime2_project/raw_data

# Verify files
ls

# Compress FASTQ files
gzip *.fastq

#Import the data files (two per sample i.e. forward and reverse reads ‘R1’ and ‘R2’ respectively) #and export them as a single qiime2 artefact file.
#If samples are already demultiplexed (i.e. sequences from each sample have been written to #separate files), a metadata file is not initially required.

# -----------------------------------------------------------------------------
# STEP 2: Import demultiplexed paired-end sequences
# -----------------------------------------------------------------------------

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ~/qiime2_project/raw_data \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path ~/qiime2_project/seqs/16S_demuxed.qza

# -----------------------------------------------------------------------------
# STEP 3: Remove primers using Cutadapt
# -----------------------------------------------------------------------------
# 16S primers (515F/806R)

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ~/qiime2_project/seqs/16S_demuxed.qza \
  --p-front-f GTGCCAGCMGCCGCGGTAA \
  --p-front-r GGACTACHVGGGTWTCTAAT \
  --p-error-rate 0.20 \
  --output-dir ~/qiime2_project/trim \
  --verbose

# ITS2 primers (fITS7/ITS4)
# Uncomment if processing ITS2 data

# qiime cutadapt trim-paired \
#   --i-demultiplexed-sequences ~/qiime2_project/seqs/ITS_demuxed.qza \
#   --p-front-f GTGARTCATCGAATCTTTG \
#   --p-front-r TCCTCCGCTTATTGATATGC \
#   --p-error-rate 0.20 \
#   --output-dir ~/qiime2_project/trim_ITS2 \
#   --verbose

# -----------------------------------------------------------------------------
# STEP 4: Summarize sequence quality
# -----------------------------------------------------------------------------

qiime demux summarize \
  --i-data ~/qiime2_project/trim/trimmed_sequences.qza \
  --o-visualization ~/qiime2_project/trim/trimmed_sequences.qzv

# View .qzv files at:
# https://view.qiime2.org

# -----------------------------------------------------------------------------
# STEP 5: Denoise sequences using DADA2
# -----------------------------------------------------------------------------

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ~/qiime2_project/trim/trimmed_sequences.qza \
  --p-trunc-len-f 220 \ #These are for 16S, 200 used for ITS2
  --p-trunc-len-r 150 \
  --p-n-threads 3 \
  --output-dir ~/qiime2_project/dada2out \
  --verbose

# Outputs:
# - table.qza
# - representative_sequences.qza
# - denoising_stats.qza

# -----------------------------------------------------------------------------
# STEP 6: Summarize denoising outputs
# -----------------------------------------------------------------------------

qiime feature-table summarize \
  --i-table ~/qiime2_project/dada2out/table.qza \
  --m-sample-metadata-file ~/qiime2_project/metadata.tsv \
  --o-visualization ~/qiime2_project/dada2out/16s_table.qzv \
  --verbose

qiime feature-table tabulate-seqs \
  --i-data ~/qiime2_project/dada2out/representative_sequences.qza \
  --o-visualization ~/qiime2_project/dada2out/16s_rep_seqs.qzv \
  --verbose

qiime metadata tabulate \
  --m-input-file ~/qiime2_project/dada2out/denoising_stats.qza \
  --o-visualization ~/qiime2_project/dada2out/16s_denoising_stats.qzv \
  --verbose

# -----------------------------------------------------------------------------
# STEP 7: Export feature table
# -----------------------------------------------------------------------------

qiime tools export \
  --input-path ~/qiime2_project/dada2out/table.qza \
  --output-path ~/qiime2_project/export/table_export

# -----------------------------------------------------------------------------
# STEP 8: Assign taxonomy
# -----------------------------------------------------------------------------
# Use ONE classifier appropriate for your marker gene:
# - SILVA for 16S
# - UNITE for ITS #both are provided

qiime feature-classifier classify-sklearn \
  --i-classifier ~/qiime2_project/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads ~/qiime2_project/dada2out/representative_sequences.qza \
  --p-n-jobs 1 \
  --output-dir ~/qiime2_project/taxonomy \
  --verbose

# -----------------------------------------------------------------------------
# STEP 9: Visualize taxonomy assignments
# -----------------------------------------------------------------------------

qiime metadata tabulate \
  --m-input-file ~/qiime2_project/taxonomy/classification.qza \
  --o-visualization ~/qiime2_project/taxonomy/taxonomy.qzv \
  --verbose

# -----------------------------------------------------------------------------
# STEP 10: Filter mitochondria and chloroplasts
# -----------------------------------------------------------------------------

qiime taxa filter-table \
  --i-table ~/qiime2_project/dada2out/table.qza \
  --i-taxonomy ~/qiime2_project/taxonomy/classification.qza \
  --p-exclude Mitochondria,Chloroplast \
  --o-filtered-table ~/qiime2_project/16s_table_filtered.qza \
  --verbose

qiime feature-table summarize \
  --i-table ~/qiime2_project/16s_table_filtered.qza \
  --m-sample-metadata-file ~/qiime2_project/metadata.tsv \
  --o-visualization ~/qiime2_project/16s_table_filtered.qzv \
  --verbose

# -----------------------------------------------------------------------------
# STEP 11: Sequence alignment
# -----------------------------------------------------------------------------

qiime alignment mafft \
  --i-sequences ~/qiime2_project/dada2out/representative_sequences.qza \
  --p-n-threads 15 \
  --o-alignment ~/qiime2_project/aligned_16s_rep_seqs.qza \
  --verbose

# -----------------------------------------------------------------------------
# STEP 12: Mask highly variable regions
# -----------------------------------------------------------------------------

qiime alignment mask \
  --i-alignment ~/qiime2_project/aligned_16s_rep_seqs.qza \
  --o-masked-alignment ~/qiime2_project/masked_aligned_16s_rep_seqs.qza \
  --verbose

# -----------------------------------------------------------------------------
# STEP 13: Build phylogenetic tree
# -----------------------------------------------------------------------------

qiime phylogeny fasttree \
  --i-alignment ~/qiime2_project/masked_aligned_16s_rep_seqs.qza \
  --p-n-threads 1 \
  --o-tree ~/qiime2_project/16s_unrooted_tree.qza \
  --verbose

# -----------------------------------------------------------------------------
# STEP 14: Root phylogenetic tree
# -----------------------------------------------------------------------------

qiime phylogeny midpoint-root \
  --i-tree ~/qiime2_project/16s_unrooted_tree.qza \
  --o-rooted-tree ~/qiime2_project/16s_rooted_tree.qza \
  --verbose

# -----------------------------------------------------------------------------
# STEP 15: Alpha rarefaction. Use metadata.tsv file for correct target -16S or ITS2
# -----------------------------------------------------------------------------

qiime diversity alpha-rarefaction \
  --i-table ~/qiime2_project/dada2out/table.qza \
  --i-phylogeny ~/qiime2_project/16s_rooted_tree.qza \
  --p-max-depth 10000 \
  --m-metadata-file ~/qiime2_project/metadata.tsv \
  --o-visualization ~/qiime2_project/downstream/16s_alpha_rarefaction.qzv \
  --verbose

# -----------------------------------------------------------------------------
# STEP 16: Taxonomy barplots. Use metadata.tsv file for correct target -16S or ITS2
# -----------------------------------------------------------------------------

qiime taxa barplot \
  --i-table ~/qiime2_project/16s_table_filtered.qza \
  --i-taxonomy ~/qiime2_project/taxonomy/classification.qza \
  --m-metadata-file ~/qiime2_project/metadata.tsv \
  --o-visualization ~/qiime2_project/downstream/barchart.qzv \
  --verbose

# -----------------------------------------------------------------------------
# STEP 17: Export phylogenetic tree. N.B. this output was not used in further analyses
# -----------------------------------------------------------------------------

qiime tools export \
  --input-path ~/qiime2_project/16s_unrooted_tree.qza \
  --output-path ~/qiime2_project/export/tree

# -----------------------------------------------------------------------------
# STEP 18: Export feature table
# -----------------------------------------------------------------------------

qiime tools export \
  --input-path ~/qiime2_project/16s_table_filtered.qza \
  --output-path ~/qiime2_project/export/feature_table

# Convert BIOM to TSV

biom convert \
  -i ~/qiime2_project/export/feature_table/feature-table.biom \
  -o ~/qiime2_project/export/feature-table.tsv \
  --to-tsv

# -----------------------------------------------------------------------------
# STEP 19: Export taxonomy table
# -----------------------------------------------------------------------------

qiime tools export \
  --input-path ~/qiime2_project/taxonomy/classification.qza \
  --output-path ~/qiime2_project/export/taxonomy

# =============================================================================
# End of QIIME2 workflow
# =============================================================================
