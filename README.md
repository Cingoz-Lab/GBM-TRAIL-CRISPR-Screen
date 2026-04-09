# GBM TRAIL Metabolic CRISPR/Cas9 Screen

Genome-scale CRISPR-Cas9 pooled screen to identify metabolic gene dependencies in TRAIL-resistant versus TRAIL-sensitive glioblastoma (GBM) cell lines.

---

## Overview

Glioblastoma (GBM) is the most common and aggressive primary brain tumor. TRAIL (Tumor Necrosis Factor-Related Apoptosis-Inducing Ligand) selectively induces apoptosis in tumor cells, but resistance frequently occurs. To understand the metabolic basis of TRAIL resistance, this screen applied the Sabatini Human Metabolic Gene Knockout Library to perform a genome-scale *in vivo* CRISPR/Cas9 screen comparing TRAIL-resistant (R) and TRAIL-sensitive (S) GBM cell lines before and after treatment.

**Library:** [Sabatini Human CRISPR Metabolic Knockout Library (Addgene #110066)](https://www.addgene.org/pooled-library/sabatini-human-crispr-metabolic-knockout/)  
29,790 sgRNAs targeting 2,981 metabolic genes (~10 sgRNAs/gene)

**Pipeline:** nf-core/crisprseq v2.3.0 → MAGeCK v0.5.9.5 (Robust Rank Aggregation)

---

## Key Findings

| Comparison | Gene | LFC | p-value | Interpretation |
|------------|------|-----|---------|----------------|
| Sensitive: Post vs Day-0 | **CYP2E1** | −13.5 | 0.038 | Essential in TRAIL-sensitive cells; CYP450/ROS metabolism |
| Resistant: Post vs Day-0 | **ASL** | −11.5 | 0.002 | Essential in resistant cells; arginine/urea cycle |
| Resistant: Post vs Day-0 | **CYB5A** | −4.2 | 0.043 | Electron transfer; depleted in resistant cells |
| Sensitive vs Resistant (Post) | **PLA2G4E** | −12.1 | 0.003 | More essential in sensitive cells post-treatment |

> **Data quality note:** Severe library bottlenecking — only ~3.3% of sgRNAs detected (Gini ~0.997). Results are exploratory and require independent validation. No genes pass FDR < 0.25 owing to extreme sgRNA dropout; p < 0.05 is used as the primary threshold with appropriate caveats.

---

## Repository Structure

```
├── analysis/
│   ├── scripts/
│   │   ├── last.py                    # Main pipeline: R1+R2 concat → count → test → figures
│   │   ├── 01_denovo_count.py         # De novo sgRNA counting from FASTQ
│   │   ├── 02_mageck_test.sh          # MAGeCK test (legacy, nf-core count input)
│   │   └── 05_publication_figures.py  # Standalone figure generator (figures_v2/)
│   ├── figures_last/                  # Final publication figures (PDF + PNG)
│   ├── results/
│   │   ├── last_count/                # MAGeCK count output (R1+R2 combined run)
│   │   └── last_test/                 # MAGeCK test output (R1+R2 combined run)
│   ├── mageck_test/                   # MAGeCK results from original nf-core run
│   └── reports/
│       ├── analysis_report.md         # Full methods and results
│       └── correctness_audit.md       # Cross-analysis quality audit
├── nf-core_crisprseq/
│   ├── samplesheet.csv                # Sample sheet for nf-core pipeline
│   └── CRISPR_nfcore.sh               # nf-core/crisprseq run command
├── library_sabatini_metko.tsv         # sgRNA reference library (Addgene #110066)
├── LICENSE
└── README.md
```

---

## Quick Start

### Requirements

```bash
conda create -n mageck-env -c bioconda -c conda-forge \
    mageck matplotlib seaborn scipy adjustText -y
conda activate mageck-env
```

### Run the full pipeline

```bash
conda activate mageck-env

# Step 1–4 in one command:
#   Concatenate R1+R2 FASTQs per sample
#   → MAGeCK count (combined reads)
#   → MAGeCK test (3 comparisons)
#   → 9 publication figures
python analysis/scripts/last.py

# Skip re-counting if count results already exist:
python analysis/scripts/last.py --skip-count
```

### Outputs

- `analysis/figures_last/` — 9 figures in PDF and PNG
- `analysis/figures_last/significant_genes_last.csv` — all significant genes (p < 0.05)
- `analysis/results/last_count/` — MAGeCK count table and QC summary
- `analysis/results/last_test/` — per-comparison gene summary files

---

## Pipeline Details

### Why R1+R2 concatenation?

In this paired-end dataset, ~35% of R1 reads contain the sgRNA in reverse-complement orientation (detected by RC scaffold scan). MAGeCK count only processes the forward strand. Concatenating R1 and R2 before counting gives MAGeCK access to both orientations, improving the mapping rate.

```bash
# Handled automatically by last.py — equivalent to:
cat sample_R1.fq.gz sample_R2.fq.gz > sample_combined.fq.gz
mageck count -l library.tsv --fastq sample_combined.fq.gz ...
```

### Library bottlenecking

Despite the orientation fix, the Gini index remains ~0.997 (ideal: <0.2). The dominant cause is library coverage — only ~1,006/30,197 sgRNAs were detected across all samples. This is likely due to suboptimal transduction efficiency (MOI and cell number) during viral library delivery.

**Recommendations for re-screening:**
- ≥15 million cells at transduction (500× coverage for 30K sgRNAs)
- Target MOI 0.3 for single-copy integration
- ≥300 reads/sgRNA sequencing depth (>9M reads per sample)
- QC the plasmid library before transduction

### Comparisons

| Comparison ID | Description |
|---------------|-------------|
| `Spost_vs_S0` | Essential genes in sensitive cells under treatment |
| `Rpost_vs_R0` | Essential genes in resistant cells under treatment |
| `Spost_vs_Rpost` | Differential essentiality: sensitive vs resistant (post-treatment) |

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| MAGeCK | 0.5.9.5 | sgRNA counting and RRA statistical test |
| nf-core/crisprseq | 2.3.0 | Nextflow pipeline wrapper |
| Python | ≥3.10 | Analysis scripts |
| matplotlib | — | Figure generation |
| seaborn | — | Figure styling |
| scipy | — | Spearman correlation |
| adjustText | — | Non-overlapping gene labels |

---

## Citation

If you use or adapt this analysis pipeline, please cite:

**MAGeCK:**
> Li W, et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. *Genome Biology*, 15, 554.

**Library:**
> Birsoy K, et al. (2015). An Essential Role of the Mitochondrial Electron Transport Chain in Cell Proliferation Is to Enable Aspartate Synthesis. *Cell*, 162(3):540–551.

**nf-core/crisprseq:**
> Ewels PA, et al. (2020). The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology*.

---

## License

MIT — see [LICENSE](LICENSE)
