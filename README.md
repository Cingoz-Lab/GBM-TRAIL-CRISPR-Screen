# GBM TRAIL Metabolic CRISPR/Cas9 Screen — Analysis Pipeline

A Python pipeline for analysing genome-scale CRISPR/Cas9 pooled screens with MAGeCK, including sgRNA counting, statistical testing, and publication-quality figure generation.

---

## Overview

This repository contains analysis scripts for a pooled CRISPR metabolic knockout screen. The pipeline handles paired-end FASTQ input, runs MAGeCK count and test, and produces a standardised set of publication figures.

**Key feature:** Automatically concatenates R1 and R2 reads before MAGeCK counting to recover sgRNAs sequenced in both orientations — a common issue with paired-end CRISPR screen libraries.

---

## Requirements

### Conda environment

```bash
conda create -n mageck-env -c bioconda -c conda-forge \
    mageck matplotlib seaborn scipy adjustText -y
conda activate mageck-env
```

---

## Repository Structure

```
├── analysis/
│   └── scripts/
│       ├── last.py                    # Main pipeline (recommended entry point)
│       ├── 01_denovo_count.py         # De novo sgRNA extraction from raw FASTQ
│       ├── 02_mageck_test.sh          # MAGeCK test wrapper (legacy)
│       └── 05_publication_figures.py  # Standalone figure generator
├── .gitignore
├── LICENSE
└── README.md
```

---

## Usage

### Main pipeline (`last.py`)

Runs the full analysis from raw paired-end FASTQs to figures in one command:

```bash
conda activate mageck-env
python analysis/scripts/last.py
```

**Steps performed automatically:**
1. Concatenate R1 + R2 FASTQ per sample → `results/last_combined_fastq/`
2. `mageck count` on combined reads → `results/last_count/`
3. Fix duplicate column names in count table
4. `mageck test` (3 comparisons) → `results/last_test/`
5. Generate 9 publication figures → `figures_last/`

**Skip recounting** if count results already exist:

```bash
python analysis/scripts/last.py --skip-count
```

### Expected inputs

Edit the `SAMPLES` list and `LIB_FILE` path at the top of `last.py` to match your data:

```python
# Path to sgRNA library TSV (sgRNA, sequence, gene columns)
LIB_FILE = Path("path/to/your_library.tsv")

# FASTQ directory
FASTQ_DIR = Path("Fastq/")

# Sample list: (label, condition, R1_filename, R2_filename)
SAMPLES = [
    ("R0_1", "R0", "sample_R01_R1.fq.gz", "sample_R01_R2.fq.gz"),
    ...
]
```

### Comparisons

Edit the `COMP_GROUPS` dict to define your treatment vs control groups:

```python
COMP_GROUPS = {
    "Treated_vs_Control": ("treated_1,treated_2", "control_1,control_2"),
}
```

---

## Output Files

| Path | Description |
|------|-------------|
| `figures_last/Fig1_volcano.*` | Volcano plots (LFC vs −log₁₀p) |
| `figures_last/Fig2_rank_lfc.*` | Gene rank by fold change |
| `figures_last/Fig3_heatmap.*` | LFC heatmap across comparisons |
| `figures_last/Fig4_barplot.*` | Top genes per comparison |
| `figures_last/Fig5_qc_metrics.*` | Mapping rate, Gini index, depth |
| `figures_last/Fig6_sgrna_profiles.*` | Per-sgRNA count profiles for key genes |
| `figures_last/Fig7_replicate_corr.*` | Replicate Spearman correlations |
| `figures_last/Fig8_cross_comparison.*` | Bubble plot across comparisons |
| `figures_last/Fig9_summary_table.*` | Table of all significant genes |
| `figures_last/significant_genes_last.csv` | Significant genes (p < 0.05) |

All figures are saved as both **PDF** and **PNG** at 300 dpi.

---

## Library QC Thresholds

The pipeline flags samples against standard CRISPR screen benchmarks:

| Metric | Ideal | Poor |
|--------|-------|------|
| Mapping rate | ≥ 60% | < 10% |
| Gini index | < 0.2 | > 0.9 |
| sgRNA dropout | < 20% | > 80% |
| Reads per sample | ≥ 9M | < 3M |

---

## De Novo sgRNA Counting (`01_denovo_count.py`)

If the reference library does not match the sequenced library (low mapping rate), extract sgRNA sequences directly from reads:

```bash
python analysis/scripts/01_denovo_count.py \
    --fastq-dir Fastq/ \
    --outdir results/ \
    --scaffold GTTTTAGAGCTAGAAATAGC \
    --min-reads 5
```

Outputs a count matrix and a de-novo library TSV for downstream BLAST annotation.

---

## Standalone Figure Generation (`05_publication_figures.py`)

Regenerate figures only from existing MAGeCK results:

```bash
python analysis/scripts/05_publication_figures.py
```

Reads from `mageck_test/` by default. Edit `MT_DIR` at the top of the script to point to a different results directory.

---

## Citation

If you use this pipeline, please cite:

**MAGeCK:**
> Li W, et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. *Genome Biology*, 15, 554.

---

## License

MIT — see [LICENSE](LICENSE)
