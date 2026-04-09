## ============================================================
## run_all.R — Run the complete downstream R analysis
## GBM TRAIL Metabolic CRISPR/Cas9 Screen
## ============================================================
## Usage (from reanalysis/ directory):
##   Rscript R/run_all.R
## Or in RStudio: open and Source this file.
## ============================================================


message("==========================================================")
message("GBM TRAIL Metabolic CRISPR Screen — Downstream R Analysis")
message("==========================================================")

message("\n[1/3] QC ...")
source("R/01_qc.R")

message("\n[2/3] MAGeCK visualisation ...")
source("R/02_mageck_viz.R")

message("\n[3/3] Pathway enrichment ...")
source("R/03_enrichment.R")

message("\n==========================================================")
message("Done. All figures written to: ", file.path(getwd(), "figures"))
message("==========================================================")
