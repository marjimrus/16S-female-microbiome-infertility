#!/bin/bash
# ==============================================================================
# 16S rRNA Analysis Pipeline for Ion Torrent Data using QIIME2
# ==============================================================================
#
# Description: Complete workflow for 16S rRNA gene sequencing analysis from
#              Ion Torrent platform using QIIME2 and DADA2.
#
# Author: Maria J. Rus
# Date: 2023
#
# Requirements:
#   - QIIME2 >= 2022.8
#   - SAMtools
#   - Greengenes 13_8 99% classifier
#
# Input:
#   - BAM files from Ion Torrent sequencer (zipped)
#   - Metadata file (TSV format)
#   - Manifest file (sample-id, filepath, direction)
#
# Output:
#   - Feature table (ASVs)
#   - Taxonomy assignments
#   - Phylogenetic tree
#   - Alpha/Beta diversity metrics
#
# Usage: bash qiime2-16s-iontorrent-pipeline.sh
#
# Reference: Marcos et al. (2024) Microorganisms 12:989
#            https://doi.org/10.3390/microorganisms12050989
# ==============================================================================

# ------------------------------------------------------------------------------
# CONFIGURATION - Modify these parameters according to your analysis
# ------------------------------------------------------------------------------

# Input/Output directories
WORK_DIR="."
OUTPUT_DIR="./results"
MANIFEST_FILE="manifest.tsv"
METADATA_FILE="metadata.tsv"

# DADA2 parameters
TRIM_LEFT=15          # Bases to trim from 5' end
TRUNC_LEN=0           # Truncation length (0 = no truncation, recommended for Ion Torrent)
TRUNC_Q=20            # Quality score threshold

# Diversity analysis parameters
SAMPLING_DEPTH=3400   # Rarefaction depth (adjust based on your data)

# Reference database (Greengenes)
CLASSIFIER="path/to/gg-13-8-99-515-806-nb-classifier.qza"

# Number of threads
N_THREADS=4

# ------------------------------------------------------------------------------
# STEP 0: SETUP AND DATA PREPARATION
# ------------------------------------------------------------------------------

echo "======================================"
echo "Setting up directories..."
echo "======================================"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/denoising"
mkdir -p "${OUTPUT_DIR}/taxonomy"
mkdir -p "${OUTPUT_DIR}/diversity"
mkdir -p "${OUTPUT_DIR}/exports"

# ------------------------------------------------------------------------------
# STEP 1: CONVERT BAM TO FASTQ (if needed)
# ------------------------------------------------------------------------------

echo "======================================"
echo "Converting BAM files to FASTQ..."
echo "======================================"

# Unzip BAM files if compressed
for f in *.zip; do
    if [ -f "$f" ]; then
        unzip -p "$f" > "${f%.zip}.bam"
    fi
done

# Convert BAM to FASTQ using SAMtools
for bam_file in *.bam; do
    if [ -f "$bam_file" ]; then
        samtools bam2fq "${bam_file}" > "${bam_file%.bam}.fastq"
        echo "Converted: ${bam_file}"
    fi
done

# ------------------------------------------------------------------------------
# STEP 2: IMPORT DATA INTO QIIME2
# ------------------------------------------------------------------------------

echo "======================================"
echo "Importing sequences into QIIME2..."
echo "======================================"

# Import demultiplexed single-end sequences
# Note: Manifest file format (TSV with columns: sample-id, absolute-filepath, direction)
qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "${MANIFEST_FILE}" \
    --output-path "${OUTPUT_DIR}/demux-seqs.qza" \
    --input-format SingleEndFastqManifestPhred33V2

# Generate demultiplexing summary
qiime demux summarize \
    --i-data "${OUTPUT_DIR}/demux-seqs.qza" \
    --o-visualization "${OUTPUT_DIR}/demux-seqs.qzv"

echo "View quality plots: qiime tools view ${OUTPUT_DIR}/demux-seqs.qzv"

# ------------------------------------------------------------------------------
# STEP 3: DENOISING WITH DADA2 (Pyrosequencing mode)
# ------------------------------------------------------------------------------

echo "======================================"
echo "Denoising sequences with DADA2..."
echo "======================================"

# DADA2 denoise-pyro is specifically designed for Ion Torrent data
# It handles the homopolymer errors characteristic of pyrosequencing
qiime dada2 denoise-pyro \
    --i-demultiplexed-seqs "${OUTPUT_DIR}/demux-seqs.qza" \
    --p-trunc-len ${TRUNC_LEN} \
    --p-trunc-q ${TRUNC_Q} \
    --p-trim-left ${TRIM_LEFT} \
    --p-n-threads ${N_THREADS} \
    --o-table "${OUTPUT_DIR}/denoising/feature-table.qza" \
    --o-representative-sequences "${OUTPUT_DIR}/denoising/rep-seqs.qza" \
    --o-denoising-stats "${OUTPUT_DIR}/denoising/stats.qza" \
    --verbose

# Generate feature table summary
qiime feature-table summarize \
    --i-table "${OUTPUT_DIR}/denoising/feature-table.qza" \
    --m-sample-metadata-file "${METADATA_FILE}" \
    --o-visualization "${OUTPUT_DIR}/denoising/feature-table-summary.qzv"

# Visualize representative sequences
qiime feature-table tabulate-seqs \
    --i-data "${OUTPUT_DIR}/denoising/rep-seqs.qza" \
    --o-visualization "${OUTPUT_DIR}/denoising/rep-seqs.qzv"

# Denoising statistics
qiime metadata tabulate \
    --m-input-file "${OUTPUT_DIR}/denoising/stats.qza" \
    --o-visualization "${OUTPUT_DIR}/denoising/stats.qzv"

# ------------------------------------------------------------------------------
# STEP 4: TAXONOMIC CLASSIFICATION
# ------------------------------------------------------------------------------

echo "======================================"
echo "Classifying taxonomy with Greengenes..."
echo "======================================"

# Classify ASVs using pre-trained Naive Bayes classifier
qiime feature-classifier classify-sklearn \
    --i-reads "${OUTPUT_DIR}/denoising/rep-seqs.qza" \
    --i-classifier "${CLASSIFIER}" \
    --p-n-jobs ${N_THREADS} \
    --o-classification "${OUTPUT_DIR}/taxonomy/taxonomy.qza"

# Visualize taxonomy
qiime metadata tabulate \
    --m-input-file "${OUTPUT_DIR}/taxonomy/taxonomy.qza" \
    --o-visualization "${OUTPUT_DIR}/taxonomy/taxonomy.qzv"

# Generate taxonomy barplot
qiime taxa barplot \
    --i-table "${OUTPUT_DIR}/denoising/feature-table.qza" \
    --i-taxonomy "${OUTPUT_DIR}/taxonomy/taxonomy.qza" \
    --m-metadata-file "${METADATA_FILE}" \
    --o-visualization "${OUTPUT_DIR}/taxonomy/taxa-barplot.qzv"

# ------------------------------------------------------------------------------
# STEP 5: PHYLOGENETIC TREE CONSTRUCTION
# ------------------------------------------------------------------------------

echo "======================================"
echo "Building phylogenetic tree..."
echo "======================================"

# Multiple sequence alignment with MAFFT
qiime alignment mafft \
    --i-sequences "${OUTPUT_DIR}/denoising/rep-seqs.qza" \
    --p-n-threads ${N_THREADS} \
    --o-alignment "${OUTPUT_DIR}/denoising/aligned-rep-seqs.qza"

# Mask highly variable positions
qiime alignment mask \
    --i-alignment "${OUTPUT_DIR}/denoising/aligned-rep-seqs.qza" \
    --o-masked-alignment "${OUTPUT_DIR}/denoising/masked-aligned-rep-seqs.qza"

# Build unrooted tree with FastTree
qiime phylogeny fasttree \
    --i-alignment "${OUTPUT_DIR}/denoising/masked-aligned-rep-seqs.qza" \
    --p-n-threads ${N_THREADS} \
    --o-tree "${OUTPUT_DIR}/denoising/unrooted-tree.qza"

# Root tree at midpoint
qiime phylogeny midpoint-root \
    --i-tree "${OUTPUT_DIR}/denoising/unrooted-tree.qza" \
    --o-rooted-tree "${OUTPUT_DIR}/denoising/rooted-tree.qza"

# ------------------------------------------------------------------------------
# STEP 6: ALPHA AND BETA DIVERSITY ANALYSIS
# ------------------------------------------------------------------------------

echo "======================================"
echo "Calculating diversity metrics..."
echo "======================================"

# Alpha rarefaction curves (to determine appropriate sampling depth)
qiime diversity alpha-rarefaction \
    --i-table "${OUTPUT_DIR}/denoising/feature-table.qza" \
    --i-phylogeny "${OUTPUT_DIR}/denoising/rooted-tree.qza" \
    --m-metadata-file "${METADATA_FILE}" \
    --p-min-depth 10 \
    --p-max-depth 5000 \
    --o-visualization "${OUTPUT_DIR}/diversity/alpha-rarefaction.qzv"

# Core diversity metrics (alpha + beta)
qiime diversity core-metrics-phylogenetic \
    --i-table "${OUTPUT_DIR}/denoising/feature-table.qza" \
    --i-phylogeny "${OUTPUT_DIR}/denoising/rooted-tree.qza" \
    --m-metadata-file "${METADATA_FILE}" \
    --p-sampling-depth ${SAMPLING_DEPTH} \
    --p-n-jobs-or-threads ${N_THREADS} \
    --output-dir "${OUTPUT_DIR}/diversity/core-metrics"

# Faith's Phylogenetic Diversity (alpha)
qiime diversity alpha-phylogenetic \
    --i-phylogeny "${OUTPUT_DIR}/denoising/rooted-tree.qza" \
    --i-table "${OUTPUT_DIR}/denoising/feature-table.qza" \
    --p-metric faith_pd \
    --o-alpha-diversity "${OUTPUT_DIR}/diversity/faith-pd.qza"

# Weighted UniFrac (beta)
qiime diversity beta-phylogenetic \
    --i-phylogeny "${OUTPUT_DIR}/denoising/rooted-tree.qza" \
    --i-table "${OUTPUT_DIR}/denoising/feature-table.qza" \
    --p-metric weighted_unifrac \
    --o-distance-matrix "${OUTPUT_DIR}/diversity/weighted-unifrac.qza"

# Unweighted UniFrac (beta)
qiime diversity beta-phylogenetic \
    --i-phylogeny "${OUTPUT_DIR}/denoising/rooted-tree.qza" \
    --i-table "${OUTPUT_DIR}/denoising/feature-table.qza" \
    --p-metric unweighted_unifrac \
    --o-distance-matrix "${OUTPUT_DIR}/diversity/unweighted-unifrac.qza"

# ------------------------------------------------------------------------------
# STEP 7: EXPORT DATA FOR R ANALYSIS
# ------------------------------------------------------------------------------

echo "======================================"
echo "Exporting data for downstream analysis..."
echo "======================================"

# Export feature table
qiime tools export \
    --input-path "${OUTPUT_DIR}/denoising/feature-table.qza" \
    --output-path "${OUTPUT_DIR}/exports/feature-table"

# Convert BIOM to TSV
biom convert \
    -i "${OUTPUT_DIR}/exports/feature-table/feature-table.biom" \
    -o "${OUTPUT_DIR}/exports/feature-table/feature-table.tsv" \
    --to-tsv

# Export taxonomy
qiime tools export \
    --input-path "${OUTPUT_DIR}/taxonomy/taxonomy.qza" \
    --output-path "${OUTPUT_DIR}/exports/taxonomy"

# Export distance matrices
qiime tools export \
    --input-path "${OUTPUT_DIR}/diversity/core-metrics/bray_curtis_distance_matrix.qza" \
    --output-path "${OUTPUT_DIR}/exports/bray-curtis"

qiime tools export \
    --input-path "${OUTPUT_DIR}/diversity/weighted-unifrac.qza" \
    --output-path "${OUTPUT_DIR}/exports/weighted-unifrac"

qiime tools export \
    --input-path "${OUTPUT_DIR}/diversity/core-metrics/jaccard_distance_matrix.qza" \
    --output-path "${OUTPUT_DIR}/exports/jaccard"

# Export Faith's PD
qiime tools export \
    --input-path "${OUTPUT_DIR}/diversity/faith-pd.qza" \
    --output-path "${OUTPUT_DIR}/exports/faith-pd"

echo "======================================"
echo "Pipeline completed successfully!"
echo "======================================"
echo ""
echo "Output files are located in: ${OUTPUT_DIR}"
echo ""
echo "Next steps:"
echo "  1. Review quality plots in ${OUTPUT_DIR}/demux-seqs.qzv"
echo "  2. Check rarefaction curves in ${OUTPUT_DIR}/diversity/alpha-rarefaction.qzv"
echo "  3. Explore taxonomy in ${OUTPUT_DIR}/taxonomy/taxa-barplot.qzv"
echo "  4. Use exported files in ${OUTPUT_DIR}/exports/ for R analysis"
echo ""
