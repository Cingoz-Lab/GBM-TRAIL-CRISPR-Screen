#!/bin/bash
# =============================================================================
# CRISPR Pooled Screen Analysis — nf-core/crisprseq v2.3.0
# =============================================================================
# This script runs the nf-core/crisprseq pipeline in screening mode.
# It performs:
#   1. Quality control (FastQC)
#   2. sgRNA counting (MAGeCK count)
#   3. Differential abundance analysis (MAGeCK test)
#   4. MultiQC report aggregation
#
# Requirements:
#   - Nextflow >= 22.10
#   - Conda environment named 'nf-env' with nf-core installed
#     (or Docker — change -profile below)
#   - FASTQ files placed in the Fastq/ directory
#   - samplesheet.csv, contrasts.csv, library TSV in the same directory
#
# Usage:
#   bash CRISPR_nfcore.sh
# =============================================================================

set -e

# Always run from the directory where this script lives
cd "$(dirname "$0")"

# Initialize conda so 'conda activate' works inside scripts
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate the environment that has nextflow + nf-core installed
conda activate nf-env

# Run the pipeline
# -profile: use 'docker' if Docker is running, or 'conda' otherwise
# -resume: safely restart from the last successful step if interrupted
nextflow run nf-core/crisprseq -r 2.3.0 \
  -c nf_fix.config \
  --analysis screening \
  --input samplesheet.csv \
  --library library_sabatini_metko.tsv \
  --contrasts contrasts.csv \
  --outdir crisprseq_out \
  -profile conda \
  -resume
