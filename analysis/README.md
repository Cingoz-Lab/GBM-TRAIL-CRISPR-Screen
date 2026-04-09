# GBM TRAIL Metabolic CRISPR/Cas9 Screen — Analysis
**Sabatini Human Metabolic Gene Knockout Library (Addgene #110066)**

A genome-scale CRISPR-Cas9 pooled screen investigating metabolic gene essentiality in TRAIL-resistant versus TRAIL-sensitive GBM cell lines, with pre- and post-treatment comparisons.

---

## Quick Start

```bash
# 1. Activate the MAGeCK conda environment
conda activate mageck-env

# 2. Run the full pipeline (R1+R2 concat → MAGeCK count → test → figures)
python scripts/last.py

# Or skip re-counting if count results already exist:
python scripts/last.py --skip-count

# Full report is in:
# reports/analysis_report.md
# reports/correctness_audit.md
```

---

## Project Structure

```
analysis/
├── scripts/
│   ├── last.py                    # Main pipeline: R1+R2 concat → count → test → figures
│   ├── 01_denovo_count.py         # De novo sgRNA extraction from FASTQ
│   ├── 02_mageck_test.sh          # MAGeCK RRA statistical test (legacy)
│   ├── 05_publication_figures.py  # Standalone figure generation (figures_v2/)
├── mageck_test/
│   ├── count_table_fixed.txt      # Count matrix from nf-core run (original)
│   ├── Spost_vs_S0.*              # MAGeCK output: sensitive comparison
│   ├── Rpost_vs_R0.*              # MAGeCK output: resistant comparison
│   └── Spost_vs_Rpost.*           # MAGeCK output: sensitive vs resistant
├── results/
│   ├── last_count/                # MAGeCK count output (R1+R2 combined)
│   └── last_test/                 # MAGeCK test output (R1+R2 combined)
├── figures_last/                  # Final publication figures (PDF + PNG)
│   ├── Fig1_volcano.*
│   ├── Fig2_rank_lfc.*
│   ├── Fig3_heatmap.*
│   ├── Fig4_barplot.*
│   ├── Fig5_qc_metrics.*
│   ├── Fig6_sgrna_profiles.*
│   ├── Fig7_replicate_corr.*
│   ├── Fig8_cross_comparison.*
│   ├── Fig9_summary_table.*
│   └── significant_genes_last.csv
└── reports/
    ├── analysis_report.md         # Full analysis report with methods and results
    └── correctness_audit.md       # Comparison with company analysis report
```

---

## Background

### Library
- **Addgene #110066** — Human CRISPR Metabolic Gene Knockout Library (Sabatini Lab)
- 29,790 sgRNAs targeting 2,981 human metabolic genes
- ~10 sgRNAs per gene
- Vector: lentiCRISPR v1 (Addgene #49535)
- Reference: Birsoy et al., *Cell* 2015

### Experimental Design

| Condition | Label | Description |
|-----------|-------|-------------|
| R0 | Resistant Day-0 | Baseline, pre-treatment (3 replicates) |
| Rpost | Resistant Post | Post-treatment, resistant line (4 replicates) |
| S0 | Sensitive Day-0 | Baseline, pre-treatment (3 replicates) |
| Spost | Sensitive Post | Post-treatment, TRAIL-sensitive line (3 replicates) |

### Comparisons
| ID | Treatment → Control | Question |
|----|---------------------|---------|
| Spost_vs_S0 | Spost → S0 | Genes essential for sensitive cells under treatment |
| Rpost_vs_R0 | Rpost → R0 | Genes essential for resistant cells under treatment |
| Spost_vs_Rpost | Spost → Rpost | Differential essentiality: sensitive vs resistant |

---

## Pipeline Details

### Step 1 — sgRNA Counting (`nf-core/crisprseq`)

```bash
# Run nf-core/crisprseq (from nf-core_crisprseq/ directory)
bash CRISPR_nfcore.sh
```

The nf-core pipeline performs:
1. **FastQC** — raw read quality assessment
2. **TrimGalore** — adapter and quality trimming
3. **MAGeCK count** — align reads to sgRNA library, generate count table

Output: `nf-core_crisprseq/crisprseq_out/mageck/count/count_table.count.txt`

**Known issue:** The nf-core count table has duplicate column names (all R0 columns labeled "R0"). The analysis pipeline fixes this automatically.

### Step 2 — MAGeCK Test (`scripts/02_mageck_test.sh`)

Uses MAGeCK RRA (Robust Rank Aggregation) algorithm:
- Normalizes read counts per million (or no normalization for sparse data)
- Calculates per-sgRNA fold change and negative-binomial p-values
- Aggregates sgRNA-level scores to gene level using RRA
- Applies Benjamini-Hochberg FDR correction

```bash
# Dependencies
conda install -c bioconda mageck

# Run
bash scripts/02_mageck_test.sh
```

### Step 3 — De Novo Counting (Alternative, `scripts/01_denovo_count.py`)

If the reference library file does not match the sequenced library:

```bash
python scripts/01_denovo_count.py \
    --fastq-dir ../Fastq/ \
    --outdir results/ \
    --scaffold GTTTTAGAGCTAGAAATAGC \
    --min-reads 5
```

This extracts all 20-mers upstream of the CRISPR scaffold without requiring a reference, then outputs:
- `results/denovo_count_table.txt` — count matrix for MAGeCK
- `results/denovo_library.tsv` — discovered sgRNA sequences (gene names require BLAST annotation)

### Step 4 — Figure Generation (`scripts/04_generate_figures.py`)

```bash
# Dependencies
conda install -c conda-forge matplotlib seaborn scipy adjustText

# Run
python scripts/04_generate_figures.py
```

Generates 9 figures (PDF + PNG):

| Figure | Content |
|--------|---------|
| Fig1 | Volcano plots: LFC vs −log₁₀(p) for all comparisons |
| Fig2 | Gene rank plots: genes sorted by LFC |
| Fig3 | Heatmap: LFC of significant genes across comparisons |
| Fig4 | Barplots: top 20 genes per comparison |
| Fig5 | Library QC: read depth, dropout rate, Gini index |
| Fig6 | Bubble plot: significance across comparisons |
| Fig7 | sgRNA count profiles for key biological genes |
| Fig8 | Replicate correlation scatter plots |
| Fig9 | Summary table of all significant genes |

---

## Key Results

### Significant Genes (p < 0.05)

| Comparison | Gene | LFC | p-value | Function |
|------------|------|-----|---------|---------|
| Spost_vs_S0 | **CYP2E1** | −13.5 | 0.038 | Drug metabolism, ROS; essential in sensitive cells |
| Spost_vs_S0 | NAT2 | +5.6 | 0.038 | Drug acetylation; enriched in sensitive post |
| Rpost_vs_R0 | **ASL** | −11.5 | 0.002 | Urea cycle; essential in resistant cells |
| Rpost_vs_R0 | **CYB5A** | −4.2 | 0.043 | Electron transfer; essential in resistant cells |
| Spost_vs_Rpost | **PLA2G4E** | −12.1 | 0.003 | More essential in sensitive than resistant (post-treatment) |

---

## Data Quality Warning

This screen has severe library bottlenecking:
- Only **1,006 / 30,197 sgRNAs** (3.3%) detected in sequencing data
- **Gini index ~0.997** (ideal: <0.2) — reads concentrated in very few sgRNAs
- **Low mapping rate** (1–10%; ideal: >60%)

**Results are exploratory.** Key hits (CYP2E1, ASL, PLA2G4E) are biologically plausible and consistent with published metabolic screen data, but must be validated.

### Why Bottlenecking Occurs

Library bottlenecking happens when:
1. Low MOI during lentiviral transduction (fewer cells receive the library)
2. Insufficient cell numbers at transduction (<500× coverage per sgRNA)
3. Library amplification issues during viral production

### Recommendations for Re-screening

1. **Target MOI 0.3** during transduction (ensures single-copy integration)
2. **Minimum cell number:** 500 × 30,197 = **15.1 million cells** at transduction
3. **Sequencing depth:** >300 reads/sgRNA → **>9M reads** per sample
4. **Library QC:** Sequence the plasmid library to verify sgRNA representation before use
5. **Use lentiCRISPR v2** (Addgene #52961) — updated vector recommended by Addgene

---

## Environment Setup

```bash
# Create and activate environment
conda create -n mageck-env -c bioconda -c conda-forge \
    mageck matplotlib seaborn scipy adjustText -y
conda activate mageck-env
```

Or install dependencies individually:
```bash
conda install -c bioconda mageck
conda install -c conda-forge matplotlib seaborn scipy adjustText
```

---

## References

1. Birsoy K, et al. (2015). An Essential Role of the Mitochondrial Electron Transport Chain in Cell Proliferation Is to Enable Aspartate Synthesis. *Cell*, 162(3):540-551. [PMC4501826]
2. Li W, et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. *Genome Biology*, 15, 554.
3. Ewels PA, et al. (2020). The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology*.

---

## Contact

For questions about this analysis, refer to `reports/analysis_report.md` for the full methods and results description.
