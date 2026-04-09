## ============================================================
## 03_enrichment.R — Pathway & GO Enrichment Analysis
## GBM TRAIL Metabolic CRISPR/Cas9 Screen
## ============================================================
## Two approaches:
##   (A) ORA  — Over-Representation Analysis on significant genes
##   (B) GSEA — Gene Set Enrichment on all genes ranked by LFC
## Runs GO (BP) and KEGG for each comparison.
## Bubble plot style adapted from RSeqPlots/GO_BP_bubble.R.
## ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(ggrepel)
})

# ── Paths ──────────────────────────────────────────────────────
PROJ_ROOT   <- "/home/emrebora/Desktop/AC_Lab/CRISPR-Seq"
RESULTS_DIR <- file.path(PROJ_ROOT, "analysis/results/last_test")
FIG_DIR     <- file.path(PROJ_ROOT, "reanalysis/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

PVAL_THRESH  <- 0.05
ENRICH_PVAL  <- 0.10     # relaxed for small gene lists
TOP_TERMS    <- 15

COMPARISONS <- list(
  Spost_vs_S0    = "Sensitive: Post vs Day-0",
  Rpost_vs_R0    = "Resistant: Post vs Day-0",
  Spost_vs_Rpost = "Sensitive vs Resistant (Post)"
)

save_plot <- function(p, name, w = 10, h = 7) {
  tryCatch({
    ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h)
    ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
    message("  ✓  ", name)
  }, error = function(e) message("  [WARN] Could not save ", name, ": ", e$message))
}

# ── Helper: symbol → Entrez ID ────────────────────────────────
sym_to_entrez <- function(symbols) {
  ids <- mapIds(org.Hs.eg.db, keys = symbols,
                keytype = "SYMBOL", column = "ENTREZID",
                multiVals = "first")
  ids[!is.na(ids)]
}

# ── Helper: custom bubble plot ────────────────────────────────
# Follows the style from RSeqPlots/GO_BP_bubble.R / GO_KEGG_bubble.R
bubble_plot <- function(enrich_result, title, top_n = TOP_TERMS) {
  df <- as.data.frame(enrich_result)
  if (nrow(df) == 0) return(NULL)

  df <- df |>
    arrange(p.adjust) |>
    slice_head(n = top_n) |>
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]; as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      neg_log10_padj = -log10(p.adjust),
      Description    = stringr::str_wrap(Description, width = 45)
    ) |>
    arrange(neg_log10_padj)

  ggplot(df, aes(x = GeneRatio_num,
                 y = reorder(Description, neg_log10_padj),
                 size = Count, colour = neg_log10_padj)) +
    geom_point(alpha = 0.9) +
    scale_colour_gradient(low = "#AED6F1", high = "#1A5276",
                          name = expression(-log[10](adj.p))) +
    scale_size_continuous(name = "Gene count", range = c(3, 10)) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y  = element_text(size = 8.5),
      plot.title   = element_text(face = "bold", size = 10),
      legend.position = "right"
    )
}

# ── Load all comparisons ───────────────────────────────────────
load_comparison <- function(name) {
  path <- file.path(RESULTS_DIR, paste0(name, ".gene_summary.txt"))
  if (!file.exists(path)) return(NULL)
  read.delim(path, stringsAsFactors = FALSE) |>
    mutate(
      LFC   = neg.lfc,
      pval  = pmin(neg.p.value, pos.p.value),
      sig   = pval < PVAL_THRESH
    )
}

data_list <- lapply(names(COMPARISONS), load_comparison)
names(data_list) <- names(COMPARISONS)
data_list <- Filter(Negate(is.null), data_list)

# ── Run enrichment per comparison ────────────────────────────
message("── Running enrichment analyses …")

for (nm in names(data_list)) {
  df    <- data_list[[nm]]
  label <- COMPARISONS[[nm]]
  message("\n  Comparison: ", label)

  # ── A. ORA on significant genes ──────────────────────────────
  sig_symbols <- df |> filter(sig) |> pull(id)
  bg_symbols  <- df |> pull(id)          # all testable genes as background

  if (length(sig_symbols) < 3) {
    message("  [SKIP ORA] Too few significant genes (", length(sig_symbols), ")")
  } else {
    sig_ids <- sym_to_entrez(sig_symbols)
    bg_ids  <- sym_to_entrez(bg_symbols)
    message("  ORA: ", length(sig_ids), " genes / ", length(bg_ids), " background")

    # GO BP
    tryCatch({
      go_bp <- enrichGO(
        gene          = sig_ids,
        universe      = bg_ids,
        OrgDb         = org.Hs.eg.db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = ENRICH_PVAL,
        qvalueCutoff  = 0.2,
        readable      = TRUE,
        minGSSize     = 5,
        maxGSSize     = 500
      )
      if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
        p_go <- bubble_plot(go_bp, paste0("GO Biological Process\n", label))
        if (!is.null(p_go))
          save_plot(p_go, paste0("Fig_ORA_GO_BP_", nm), w = 11, h = 7)
      } else {
        message("  [no GO-BP terms at padj<", ENRICH_PVAL, "]")
      }
    }, error = function(e) message("  [WARN] GO-BP failed: ", e$message))

    # KEGG
    tryCatch({
      kegg <- enrichKEGG(
        gene          = sig_ids,
        universe      = bg_ids,
        organism      = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff  = ENRICH_PVAL,
        qvalueCutoff  = 0.2,
        minGSSize     = 5
      )
      if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
        kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        p_kegg <- bubble_plot(kegg, paste0("KEGG Pathways\n", label))
        if (!is.null(p_kegg))
          save_plot(p_kegg, paste0("Fig_ORA_KEGG_", nm), w = 11, h = 7)
      } else {
        message("  [no KEGG terms at padj<", ENRICH_PVAL, "]")
      }
    }, error = function(e) message("  [WARN] KEGG failed: ", e$message))
  }

  # ── B. GSEA on full ranked gene list ────────────────────────
  # Rank by signed LFC (negative = depleted = essential)
  ranked <- df |>
    filter(!is.na(LFC)) |>
    arrange(desc(LFC)) |>
    mutate(entrez = sym_to_entrez(id)[id])

  ranked <- ranked[!is.na(ranked$entrez), ]

  if (nrow(ranked) < 20) {
    message("  [SKIP GSEA] Too few ranked genes (", nrow(ranked), ")")
    next
  }

  gene_rank <- setNames(ranked$LFC, ranked$entrez)

  # GSEA — GO BP
  tryCatch({
    gsea_go <- gseGO(
      geneList      = gene_rank,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      minGSSize     = 5,
      maxGSSize     = 500,
      pvalueCutoff  = 0.15,
      pAdjustMethod = "BH",
      verbose       = FALSE
    )
    if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
      gsea_go <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

      # Dotplot top enriched + depleted
      p_gsea <- dotplot(gsea_go, showCategory = 12, split = ".sign") +
        facet_grid(. ~ .sign) +
        labs(title = paste0("GSEA — GO Biological Process\n", label)) +
        theme(plot.title = element_text(face = "bold", size = 10))
      save_plot(p_gsea, paste0("Fig_GSEA_GO_BP_", nm), w = 13, h = 7)
    } else {
      message("  [no GSEA-GO terms]")
    }
  }, error = function(e) message("  [WARN] GSEA-GO failed: ", e$message))

  # GSEA — KEGG
  tryCatch({
    gsea_kegg <- gseKEGG(
      geneList      = gene_rank,
      organism      = "hsa",
      minGSSize     = 5,
      maxGSSize     = 500,
      pvalueCutoff  = 0.15,
      pAdjustMethod = "BH",
      verbose       = FALSE
    )
    if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
      gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

      p_gsea_k <- dotplot(gsea_kegg, showCategory = 12, split = ".sign") +
        facet_grid(. ~ .sign) +
        labs(title = paste0("GSEA — KEGG Pathways\n", label)) +
        theme(plot.title = element_text(face = "bold", size = 10))
      save_plot(p_gsea_k, paste0("Fig_GSEA_KEGG_", nm), w = 13, h = 7)
    } else {
      message("  [no GSEA-KEGG terms]")
    }
  }, error = function(e) message("  [WARN] GSEA-KEGG failed: ", e$message))
}

message("\n── Enrichment complete. Figures in: ", FIG_DIR)
