# Downstream R Analysis — GBM TRAIL Metabolic CRISPR Screen

R-based downstream analysis of MAGeCK output. Three scripts cover QC, visualisation, and pathway enrichment.

## Requirements

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork",
                   "pheatmap", "RColorBrewer", "scales",
                   "ggrepel", "cowplot", "stringr"))

if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db",
                       "enrichplot", "DESeq2"))
```

## Run

```r
# In RStudio — open reanalysis/ as project, then:
source("R/run_all.R")

# Or from terminal:
Rscript R/run_all.R
```

## Scripts

| Script | Input | Output |
|--------|-------|--------|
| `01_qc.R` | `last_count.countsummary.txt`, count matrix | Mapping rate, Gini, depth, dropout, replicate correlation heatmap |
| `02_mageck_viz.R` | MAGeCK gene summary files | Volcano, rank, heatmap, barplot, sgRNA profiles |
| `03_enrichment.R` | MAGeCK gene summary files | GO BP + KEGG bubble plots (ORA), GSEA dotplots |

## Output

All figures saved to `figures/` as PDF + PNG (300 dpi).

## Data Quality Note

This dataset has severe library bottlenecking (Gini ~0.997, ~3% sgRNA coverage).
Enrichment results should be treated as hypothesis-generating only.
