## ============================================================
## 01_qc.R — Library & Sequencing Quality Control
## GBM TRAIL Metabolic CRISPR/Cas9 Screen
## ============================================================
## Reads MAGeCK count summary + count matrix.
## Produces: mapping rate, Gini index, depth, dropout,
##           replicate Spearman correlation heatmap.
## ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(patchwork)
})

# ── Paths ──────────────────────────────────────────────────────
PROJ_ROOT <- "/home/emrebora/Desktop/AC_Lab/CRISPR-Seq"
COUNT_QC  <- file.path(PROJ_ROOT, "analysis/results/last_count/last_count.countsummary.txt")
COUNT_MT  <- file.path(PROJ_ROOT, "analysis/results/last_count/last_count_fixed.txt")
FIG_DIR   <- file.path(PROJ_ROOT, "reanalysis/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Colour palette (condition-consistent) ─────────────────────
COND_COLS <- c(R0 = "#2471A3", Rpost = "#CB4335",
               S0  = "#1E8449", Spost = "#CA6F1E")

save_plot <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  message("  ✓  ", name)
}

# ── 1. Read QC summary ─────────────────────────────────────────
qc <- read.delim(COUNT_QC, stringsAsFactors = FALSE) |>
  mutate(
    sample    = Label,
    condition = sub("_[0-9]+$", "", Label),
    map_pct   = Percentage * 100,
    depth_M   = Reads / 1e6,
    dropout   = Zerocounts / TotalsgRNAs * 100
  ) |>
  mutate(condition = factor(condition, levels = c("R0", "Rpost", "S0", "Spost")))

message("── QC summary loaded: ", nrow(qc), " samples")

# ── Panel A: Sequencing depth ──────────────────────────────────
pA <- ggplot(qc, aes(x = reorder(sample, depth_M), y = depth_M, fill = condition)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 5, linetype = "dashed", colour = "#555555", linewidth = 0.6) +
  annotate("text", x = 0.7, y = 5.3, label = "5M guideline",
           hjust = 0, size = 3, colour = "#555555") +
  scale_fill_manual(values = COND_COLS) +
  coord_flip() +
  labs(title = "A. Sequencing Depth", x = NULL, y = "Total reads (M)", fill = "Condition") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

# ── Panel B: Mapping rate ──────────────────────────────────────
pB <- ggplot(qc, aes(x = reorder(sample, map_pct), y = map_pct, fill = condition)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 60, linetype = "dashed", colour = "#E74C3C", linewidth = 0.6) +
  annotate("text", x = 0.7, y = 62, label = "Ideal ≥60%",
           hjust = 0, size = 3, colour = "#E74C3C") +
  scale_fill_manual(values = COND_COLS) +
  coord_flip() +
  labs(title = "B. sgRNA Mapping Rate", x = NULL, y = "Mapping rate (%)", fill = "Condition") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

# ── Panel C: Gini index ────────────────────────────────────────
pC <- ggplot(qc, aes(x = reorder(sample, GiniIndex), y = GiniIndex, fill = condition)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0.2, linetype = "dashed", colour = "#27AE60", linewidth = 0.6) +
  annotate("text", x = 0.7, y = 0.21, label = "Ideal <0.2",
           hjust = 0, size = 3, colour = "#27AE60") +
  scale_fill_manual(values = COND_COLS) +
  coord_flip() +
  labs(title = "C. Gini Index (Library Evenness)", x = NULL, y = "Gini index", fill = "Condition") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

# ── Panel D: sgRNA dropout ─────────────────────────────────────
pD <- ggplot(qc, aes(x = reorder(sample, dropout), y = dropout, fill = condition)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "#27AE60", linewidth = 0.6) +
  annotate("text", x = 0.7, y = 21, label = "Ideal <20%",
           hjust = 0, size = 3, colour = "#27AE60") +
  scale_fill_manual(values = COND_COLS) +
  coord_flip() +
  labs(title = "D. sgRNA Dropout Rate", x = NULL, y = "Zero-count sgRNAs (%)", fill = "Condition") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

# Shared legend
legend_p <- ggplot(qc, aes(x = sample, y = map_pct, fill = condition)) +
  geom_col() +
  scale_fill_manual(values = COND_COLS, name = "Condition") +
  theme_classic() +
  theme(legend.position = "bottom")
leg <- cowplot::get_legend(legend_p)

qc_combined <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title    = "Library & Sequencing Quality Control",
    subtitle = paste0(nrow(qc), " samples | R1+R2 combined FASTQs"),
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

save_plot(qc_combined, "Fig_QC_panels", w = 14, h = 10)

# ── 2. Replicate correlation heatmap ──────────────────────────
message("── Computing replicate correlations …")

ct <- read.delim(COUNT_MT, stringsAsFactors = FALSE, row.names = 1)
ct <- ct[ct$Gene != "INTERGENIC", ]
mat <- as.matrix(ct[, -1])                       # drop Gene column
mat <- mat[rowSums(mat) > 0, ]                    # remove all-zero rows
log_rpm <- log1p(sweep(mat, 2, colSums(mat), "/") * 1e6)

# Spearman correlation (pairwise complete)
cor_mat <- cor(log_rpm, method = "spearman", use = "pairwise.complete.obs")

# Annotation for heatmap columns
ann_df <- data.frame(
  Condition = sub("_[0-9]+$", "", colnames(cor_mat)),
  row.names  = colnames(cor_mat)
)
ann_colors <- list(Condition = COND_COLS[levels(factor(ann_df$Condition))])

pheatmap(
  cor_mat,
  annotation_col  = ann_df,
  annotation_colors = ann_colors,
  color           = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
  breaks          = seq(0, 1, length.out = 101),
  display_numbers = TRUE,
  number_format   = "%.2f",
  fontsize_number = 7,
  border_color    = NA,
  main            = "Replicate Spearman Correlation (log₁RPM)",
  filename        = file.path(FIG_DIR, "Fig_Replicate_Correlation.pdf"),
  width           = 9, height = 8
)
pheatmap(
  cor_mat,
  annotation_col  = ann_df,
  annotation_colors = ann_colors,
  color           = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
  breaks          = seq(0, 1, length.out = 101),
  display_numbers = TRUE,
  number_format   = "%.2f",
  fontsize_number = 7,
  border_color    = NA,
  main            = "Replicate Spearman Correlation (log₁RPM)",
  filename        = file.path(FIG_DIR, "Fig_Replicate_Correlation.png"),
  width           = 9, height = 8
)
message("  ✓  Fig_Replicate_Correlation")

message("\n── QC complete. Figures in: ", FIG_DIR)
