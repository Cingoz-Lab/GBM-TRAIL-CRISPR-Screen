#!/bin/bash
# =============================================================================
# 02_mageck_test.sh
# =================
# Runs MAGeCK test for all three CRISPR screen comparisons using either:
#   (a) The de novo count table from 01_denovo_count.py, or
#   (b) The nf-core/crisprseq-generated count table (if preferred)
#
# Comparisons:
#   - Spost_vs_S0    : Susceptible post-treatment vs Day-0 baseline
#   - Rpost_vs_R0    : Resistant post-treatment vs Day-0 baseline
#   - Spost_vs_Rpost : Susceptible vs Resistant (both post-treatment)
#
# MAGeCK test algorithm (RRA — Robust Rank Aggregation):
#   1. Normalises read counts (median normalisation by default)
#   2. Calculates sgRNA-level fold change and p-values (negative binomial)
#   3. Aggregates sgRNA scores to gene level using RRA
#   4. Applies BH FDR correction for multiple testing
#
# Requirements:
#   conda activate mageck-env   (conda install -c bioconda mageck)
#
# Usage:
#   bash scripts/02_mageck_test.sh
#   bash scripts/02_mageck_test.sh --count-table results/denovo_count_table.txt
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")/.."          # always run from project root

# ── Defaults ─────────────────────────────────────────────────────────────────
COUNT_TABLE="${1:-results/denovo_count_table.txt}"
OUT_DIR="results/mageck_test"
mkdir -p "$OUT_DIR"

echo "============================================================"
echo "MAGeCK RRA Test"
echo "Count table : $COUNT_TABLE"
echo "Output dir  : $OUT_DIR"
echo "============================================================"

# Helper function
run_mageck() {
    local treatment="$1"
    local control="$2"
    local name="$3"
    echo ""
    echo "▶  $name  (treatment=$treatment  control=$control)"
    mageck test \
        -k "$COUNT_TABLE" \
        -t "$treatment" \
        -c "$control" \
        -n "$OUT_DIR/$name" \
        --gene-lfc-method median \
        --remove-zero both \
        --remove-zero-threshold 0 \
        2>&1 | grep -E "^(INFO|WARNING|ERROR)" | sed 's/INFO  @ [^:]*: /  /; s/WARNING @ [^:]*: /  [WARN] /; s/ERROR @ [^:]*: /  [ERROR] /'
    echo "  → $OUT_DIR/$name.gene_summary.txt"
}

# ── Run all comparisons ───────────────────────────────────────────────────────
# Column labels in the count table (comma-separated for multiple replicates)

run_mageck "Spost_1,Spost_2,Spost_3" \
           "S0_1,S0_2,S0_3" \
           "Spost_vs_S0"

run_mageck "Rpost_1,Rpost_2,Rpost_3,Rpost_4" \
           "R0_1,R0_2,R0_3" \
           "Rpost_vs_R0"

run_mageck "Spost_1,Spost_2,Spost_3" \
           "Rpost_1,Rpost_2,Rpost_3,Rpost_4" \
           "Spost_vs_Rpost"

echo ""
echo "============================================================"
echo "MAGeCK test complete. Files in: $OUT_DIR"
ls "$OUT_DIR"/*.gene_summary.txt
echo "============================================================"
