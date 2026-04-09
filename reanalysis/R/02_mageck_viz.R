## ============================================================
## 02_mageck_viz.R — MAGeCK Results Visualisation
## GBM TRAIL Metabolic CRISPR/Cas9 Screen
## ============================================================
## Reads MAGeCK gene summary files (3 comparisons).
## Produces: volcano plots, rank plots, LFC heatmap,
##           top-genes barplots, sgRNA count profiles.
## ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(scales)
})

# ── Paths ──────────────────────────────────────────────────────
PROJ_ROOT   <- "/home/emrebora/Desktop/AC_Lab/CRISPR-Seq"
RESULTS_DIR <- file.path(PROJ_ROOT, "analysis/results/last_test")
COUNT_FILE  <- file.path(PROJ_ROOT, "analysis/results/last_count/last_count_fixed.txt")
FIG_DIR     <- file.path(PROJ_ROOT, "reanalysis/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

PVAL_THRESH <- 0.05
TOP_LABEL   <- 8     # max gene labels per panel

COMPARISONS <- list(
  Spost_vs_S0    = "Sensitive: Post vs Day-0",
  Rpost_vs_R0    = "Resistant: Post vs Day-0",
  Spost_vs_Rpost = "Sensitive vs Resistant (Post)"
)

COND_COLS <- c(R0 = "#2471A3", Rpost = "#CB4335",
               S0  = "#1E8449", Spost = "#CA6F1E")
C_DEPL <- "#1A5276"
C_ENRI <- "#922B21"
C_GREY <- "#BDC3C7"

save_plot <- function(p, name, w = 12, h = 5) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  message("  ✓  ", name)
}

# ── Load all comparisons ───────────────────────────────────────
load_comparison <- function(name) {
  path <- file.path(RESULTS_DIR, paste0(name, ".gene_summary.txt"))
  if (!file.exists(path)) {
    message("  [SKIP] not found: ", path); return(NULL)
  }
  df <- read.delim(path, stringsAsFactors = FALSE) |>
    mutate(
      LFC     = `neg.lfc`,
      pval    = pmin(`neg.p.value`, `pos.p.value`),
      FDR     = pmin(`neg.fdr`, `pos.fdr`),
      log10p  = -log10(pmax(pval, 1e-10)),
      sig     = pval < PVAL_THRESH,
      direction = ifelse(LFC < 0, "Depleted", "Enriched"),
      comparison = name,
      label   = COMPARISONS[[name]]
    )
  message("  ", name, ": ", nrow(df), " genes, ",
          sum(df$sig, na.rm = TRUE), " significant (p<", PVAL_THRESH, ")")
  df
}

data_list <- lapply(names(COMPARISONS), load_comparison)
names(data_list) <- names(COMPARISONS)
data_list <- Filter(Negate(is.null), data_list)

if (length(data_list) == 0) stop("No gene summary files found. Check RESULTS_DIR path.")

all_data <- bind_rows(data_list)

# ── Save significant genes CSV ─────────────────────────────────
sig_df <- all_data |>
  filter(sig) |>
  select(comparison, label, gene = id, pval, FDR, LFC, direction) |>
  arrange(comparison, pval)

write.csv(sig_df, file.path(FIG_DIR, "significant_genes.csv"), row.names = FALSE)
message("  Significant genes: ", nrow(sig_df))

# ── Figure 1: Volcano plots ────────────────────────────────────
message("── Fig 1: Volcano plots …")

volcano_list <- lapply(names(data_list), function(nm) {
  df    <- data_list[[nm]]
  label <- COMPARISONS[[nm]]

  top_sig <- df |> filter(sig) |> slice_min(pval, n = TOP_LABEL)

  ggplot(df, aes(x = LFC, y = log10p)) +
    geom_point(data = filter(df, !sig),
               colour = C_GREY, size = 1.5, alpha = 0.6) +
    geom_point(data = filter(df, sig & LFC < 0),
               colour = C_DEPL, size = 2.8, alpha = 0.9) +
    geom_point(data = filter(df, sig & LFC > 0),
               colour = C_ENRI, size = 2.8, alpha = 0.9) +
    geom_hline(yintercept = -log10(PVAL_THRESH),
               linetype = "dashed", colour = "#555555", linewidth = 0.6) +
    geom_vline(xintercept = 0,
               linetype = "dotted", colour = "#AAAAAA", linewidth = 0.5) +
    geom_text_repel(data = top_sig, aes(label = id),
                    size = 2.8, max.overlaps = 20,
                    colour = ifelse(top_sig$LFC < 0, C_DEPL, C_ENRI),
                    fontface = "bold",
                    box.padding = 0.4, point.padding = 0.3,
                    segment.colour = "#999999", segment.size = 0.3) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("▼ ", sum(df$sig & df$LFC < 0, na.rm=TRUE), " depleted"),
             hjust = -0.1, vjust = 1.5, size = 3, colour = C_DEPL, fontface = "bold") +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("▲ ", sum(df$sig & df$LFC > 0, na.rm=TRUE), " enriched"),
             hjust = -0.1, vjust = 3.2, size = 3, colour = C_ENRI) +
    labs(title = label,
         x = "Log₂ Fold Change", y = expression(-log[10](p-value))) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10))
})

volcano_combined <- wrap_plots(volcano_list, nrow = 1) +
  plot_annotation(
    title = "Volcano Plots — MAGeCK RRA (p < 0.05)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )
save_plot(volcano_combined, "Fig_Volcano", w = 15, h = 5)

# ── Figure 2: Gene rank plots ──────────────────────────────────
message("── Fig 2: Rank plots …")

rank_list <- lapply(names(data_list), function(nm) {
  df    <- data_list[[nm]] |> arrange(LFC) |> mutate(rank = seq_len(n()))
  label <- COMPARISONS[[nm]]
  top_sig <- df |> filter(sig) |> slice_min(pval, n = TOP_LABEL)

  ggplot(df, aes(x = rank, y = LFC)) +
    geom_point(data = filter(df, !sig),
               colour = C_GREY, size = 1.2, alpha = 0.6) +
    geom_point(data = filter(df, sig & LFC < 0),
               colour = C_DEPL, size = 2.5, alpha = 0.9) +
    geom_point(data = filter(df, sig & LFC > 0),
               colour = C_ENRI, size = 2.5, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted",
               colour = "#AAAAAA", linewidth = 0.5) +
    geom_text_repel(data = top_sig, aes(label = id),
                    size = 2.8, max.overlaps = 15,
                    colour = ifelse(top_sig$LFC < 0, C_DEPL, C_ENRI),
                    fontface = "bold",
                    box.padding = 0.4, segment.colour = "#999999", segment.size = 0.3) +
    labs(title = label,
         x = "Gene Rank (by LFC)", y = "Log₂ Fold Change") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10))
})

rank_combined <- wrap_plots(rank_list, nrow = 1) +
  plot_annotation(
    title = "Gene Rank Plots",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )
save_plot(rank_combined, "Fig_Rank", w = 15, h = 5)

# ── Figure 3: LFC heatmap ──────────────────────────────────────
message("── Fig 3: LFC heatmap …")

sig_genes <- unique(sig_df$gene)

if (length(sig_genes) > 0) {
  heat_mat <- sapply(names(data_list), function(nm) {
    df <- data_list[[nm]]
    lfc_vec <- setNames(df$LFC, df$id)
    lfc_vec[sig_genes]
  })
  rownames(heat_mat) <- sig_genes
  colnames(heat_mat) <- unlist(COMPARISONS[names(data_list)])
  heat_mat[is.na(heat_mat)] <- 0

  # Sort by most negative LFC
  heat_mat <- heat_mat[order(apply(heat_mat, 1, min, na.rm = TRUE)), , drop = FALSE]

  vmax <- min(ceiling(max(abs(heat_mat), na.rm = TRUE)), 16)
  bks  <- seq(-vmax, vmax, length.out = 101)

  pheatmap(
    heat_mat,
    color        = colorRampPalette(c("#154360", "#2E86C1", "#AED6F1",
                                      "white", "#F1948A", "#922B21", "#641E16"))(100),
    breaks       = bks,
    border_color = NA,
    display_numbers = TRUE,
    number_format   = "%.1f",
    fontsize_number = 8,
    fontsize_row    = 9,
    fontsize_col    = 9,
    angle_col       = 45,
    main         = paste0("Significant Genes — LFC Across Comparisons\n",
                          "(p<", PVAL_THRESH, "; 0 = not testable in that comparison)"),
    filename     = file.path(FIG_DIR, "Fig_Heatmap.pdf"),
    width = 8, height = max(4, length(sig_genes) * 0.35 + 2)
  )
  pheatmap(
    heat_mat,
    color        = colorRampPalette(c("#154360", "#2E86C1", "#AED6F1",
                                      "white", "#F1948A", "#922B21", "#641E16"))(100),
    breaks       = bks,
    border_color = NA,
    display_numbers = TRUE,
    number_format   = "%.1f",
    fontsize_number = 8,
    fontsize_row    = 9,
    fontsize_col    = 9,
    angle_col       = 45,
    main         = paste0("Significant Genes — LFC Across Comparisons\n",
                          "(p<", PVAL_THRESH, "; 0 = not testable in that comparison)"),
    filename     = file.path(FIG_DIR, "Fig_Heatmap.png"),
    width = 8, height = max(4, length(sig_genes) * 0.35 + 2)
  )
  message("  ✓  Fig_Heatmap")
} else {
  message("  [SKIP] No significant genes for heatmap")
}

# ── Figure 4: Top-genes barplots ───────────────────────────────
message("── Fig 4: Top genes barplot …")

bar_list <- lapply(names(data_list), function(nm) {
  df    <- data_list[[nm]]
  label <- COMPARISONS[[nm]]
  top15 <- df |> slice_min(pval, n = 15) |> arrange(log10p)

  ggplot(top15, aes(x = reorder(id, log10p), y = log10p,
                    fill = ifelse(LFC < 0, "Depleted", "Enriched"))) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = -log10(PVAL_THRESH),
               linetype = "dashed", colour = "#555555", linewidth = 0.5) +
    geom_text(aes(label = sprintf("p=%.4f  LFC=%+.1f", pval, LFC)),
              hjust = -0.05, size = 2.6, colour = "#444444") +
    scale_fill_manual(values = c(Depleted = C_DEPL, Enriched = C_ENRI),
                      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    coord_flip() +
    labs(title = label, x = NULL, y = expression(-log[10](p-value))) +
    theme_classic(base_size = 10) +
    theme(
      plot.title     = element_text(face = "bold", size = 9.5),
      legend.position = "bottom",
      panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = 0.4)
    )
})

bar_combined <- wrap_plots(bar_list, ncol = 1) +
  plot_annotation(
    title = "Top Genes per Comparison (ranked by p-value)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )
save_plot(bar_combined, "Fig_Barplot", w = 10, h = 4.5 * length(bar_list))

# ── Figure 5: sgRNA count profiles for key genes ───────────────
message("── Fig 5: sgRNA count profiles …")

KEY_GENES <- c("CYP2E1", "ASL", "PLA2G4E", "CYB5A", "MTF1", "TRPM7")

if (file.exists(COUNT_FILE)) {
  ct_raw <- read.delim(COUNT_FILE, stringsAsFactors = FALSE, row.names = 1)
  ct_raw <- ct_raw[!grepl("INTERGENIC", ct_raw$Gene, ignore.case = TRUE), ]
  mat_raw <- as.matrix(ct_raw[, -1])
  totals  <- colSums(mat_raw)
  rpm     <- sweep(mat_raw, 2, totals, "/") * 1e6

  sample_order <- colnames(rpm)
  cond_of      <- sub("_[0-9]+$", "", sample_order)

  genes_found <- intersect(KEY_GENES, ct_raw$Gene)
  if (length(genes_found) == 0) {
    message("  [SKIP] No key genes found in count table")
  } else {
    prof_plots <- lapply(genes_found, function(g) {
      rows  <- which(ct_raw$Gene == g)
      g_rpm <- as.data.frame(t(rpm[rows, , drop = FALSE]))
      g_rpm$sample    <- rownames(g_rpm)
      g_rpm$condition <- cond_of

      long <- g_rpm |>
        pivot_longer(cols = -c(sample, condition),
                     names_to = "sgRNA", values_to = "RPM") |>
        mutate(
          sample    = factor(sample, levels = sample_order),
          condition = factor(condition, levels = c("R0", "Rpost", "S0", "Spost"))
        )

      ggplot(long, aes(x = sample, y = RPM, colour = sgRNA, group = sgRNA)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2.5) +
        scale_colour_brewer(palette = "Set1") +
        labs(title = g, x = NULL, y = "Reads per million",
             colour = "sgRNA") +
        theme_classic(base_size = 10) +
        theme(
          axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
          plot.title  = element_text(face = "bold"),
          legend.position = "right"
        )
    })

    prof_combined <- wrap_plots(prof_plots, ncol = 1) +
      plot_annotation(
        title = "sgRNA Count Profiles — Key Hit Genes",
        theme = theme(plot.title = element_text(size = 13, face = "bold"))
      )
    save_plot(prof_combined, "Fig_sgRNA_Profiles",
              w = 11, h = 3.5 * length(genes_found))
  }
} else {
  message("  [SKIP] Count file not found: ", COUNT_FILE)
}

message("\n── Visualisation complete. Figures in: ", FIG_DIR)
