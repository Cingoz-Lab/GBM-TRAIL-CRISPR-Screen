#!/usr/bin/env python3
"""
05_figures.py
=========================

Generates 9 figures (PDF + PNG) in analysis/figures_v2/.

All code is original. Visualization approaches are standard in the field
(volcano, rank, heatmap, etc.) and are not derived from any specific
third-party repository. Dependencies: matplotlib, seaborn, scipy, numpy,
pandas, adjustText — all permissively licensed (BSD/MIT/PSF).

Usage:
    conda activate mageck-env
    python scripts/05_publication_figures.py

Output:
    figures_v2/Fig1_volcano.pdf/png
    figures_v2/Fig2_rank_lfc.pdf/png
    figures_v2/Fig3_heatmap.pdf/png
    figures_v2/Fig4_barplot.pdf/png
    figures_v2/Fig5_qc_metrics.pdf/png
    figures_v2/Fig6_sgrna_profiles.pdf/png
    figures_v2/Fig7_replicate_corr.pdf/png
    figures_v2/Fig8_cross_comparison.pdf/png
    figures_v2/Fig9_summary_table.pdf/png
    figures_v2/significant_genes_v2.csv
"""

import os
import sys
import warnings
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from scipy.stats import spearmanr, pearsonr

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False
    warnings.warn("adjustText not found — gene labels may overlap. Install with: pip install adjustText")

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE       = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MT_DIR     = os.path.join(BASE, "mageck_test")
COUNT_FILE = os.path.join(MT_DIR, "count_table_fixed.txt")
OUT_DIR    = os.path.join(BASE, "figures_v2")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Global parameters ──────────────────────────────────────────────────────────
PVAL_THRESH  = 0.05
FDR_THRESH   = 0.25
TOP_LABEL    = 10        # max gene labels per panel

COMPARISONS = {
    "Spost_vs_S0":    "Sensitive: Post vs Day-0",
    "Rpost_vs_R0":    "Resistant: Post vs Day-0",
    "Spost_vs_Rpost": "Sensitive vs Resistant (Post)",
}

# Colour palette — consistent throughout all figures
C_DEPL  = "#1A5276"   # blue  — depleted / essential
C_ENRI  = "#A93226"   # red   — enriched / toxic
C_GREY  = "#CACFD2"   # grey  — not significant
C_GOLD  = "#D4AC0D"   # gold  — multi-comparison overlap

COND_COLORS = {
    "R0":    "#2E86C1",
    "Rpost": "#E74C3C",
    "S0":    "#27AE60",
    "Spost": "#E67E22",
}

# ── Global matplotlib style ────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi":          300,
    "savefig.dpi":         300,
    "font.family":         "DejaVu Sans",
    "font.size":           10,
    "axes.titlesize":      11,
    "axes.labelsize":      10,
    "xtick.labelsize":     9,
    "ytick.labelsize":     9,
    "legend.fontsize":     9,
    "legend.frameon":      False,
    "axes.spines.top":     False,
    "axes.spines.right":   False,
    "axes.linewidth":      0.9,
    "xtick.major.width":   0.9,
    "ytick.major.width":   0.9,
    "xtick.direction":     "out",
    "ytick.direction":     "out",
    "lines.linewidth":     1.2,
})


# ── Helpers ────────────────────────────────────────────────────────────────────

def save_fig(fig: plt.Figure, stem: str) -> None:
    for ext in ("pdf", "png"):
        path = os.path.join(OUT_DIR, f"{stem}.{ext}")
        fig.savefig(path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  ✓  {stem}")


def load_gene_summary(name: str) -> pd.DataFrame | None:
    path = os.path.join(MT_DIR, f"{name}.gene_summary.txt")
    if not os.path.exists(path):
        print(f"  [WARN] Missing: {path}")
        return None
    df = pd.read_csv(path, sep="\t")
    # Unified columns
    df["LFC"]      = df["neg|lfc"]
    df["pval_neg"] = df["neg|p-value"]
    df["pval_pos"] = df["pos|p-value"]
    df["fdr_neg"]  = df["neg|fdr"]
    df["fdr_pos"]  = df["pos|fdr"]
    # Use the more extreme direction for each gene
    df["pval"]     = df[["pval_neg", "pval_pos"]].min(axis=1)
    df["FDR"]      = df[["fdr_neg",  "fdr_pos" ]].min(axis=1)
    df["log10p"]   = -np.log10(df["pval"].clip(lower=1e-10))
    df["log10fdr"] = -np.log10(df["FDR"].clip(lower=1e-10))
    df["sig"]      = df["pval"] < PVAL_THRESH
    df["label"]    = df["id"]
    return df


def sig_genes(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["sig"]].copy()


def top_labels(df: pd.DataFrame, n: int = TOP_LABEL) -> pd.DataFrame:
    """Return up to n most significant genes for labelling."""
    return df.nsmallest(n, "pval")


# ── Load all data ──────────────────────────────────────────────────────────────
print("Loading MAGeCK gene summaries …")
data = {k: load_gene_summary(k) for k in COMPARISONS}
data = {k: v for k, v in data.items() if v is not None}

print("Loading count table …")
counts_raw = pd.read_csv(COUNT_FILE, sep="\t", index_col=0)
# Remove INTERGENIC pseudo-genes
counts = counts_raw[~counts_raw.index.str.startswith("sg_INTERGENIC")]
counts = counts[~counts["Gene"].str.upper().str.contains("INTERGENIC")]
count_matrix = counts.drop(columns=["Gene"])

sample_groups = {
    "R0":    [c for c in count_matrix.columns if c.startswith("R0")],
    "Rpost": [c for c in count_matrix.columns if c.startswith("Rpost")],
    "S0":    [c for c in count_matrix.columns if c.startswith("S0")],
    "Spost": [c for c in count_matrix.columns if c.startswith("Spost")],
}

# Collect all significant genes
all_sig_rows = []
for cname, df in data.items():
    for _, row in sig_genes(df).iterrows():
        all_sig_rows.append({
            "comparison": cname,
            "comparison_label": COMPARISONS[cname],
            "gene":  row["id"],
            "pval":  row["pval"],
            "FDR":   row["FDR"],
            "LFC":   row["LFC"],
            "direction": "Depleted" if row["LFC"] < 0 else "Enriched",
        })
sig_df = pd.DataFrame(all_sig_rows)
sig_df.to_csv(os.path.join(OUT_DIR, "significant_genes_v2.csv"), index=False)
print(f"  {len(sig_df)} significant gene-hits (p < {PVAL_THRESH}) saved.\n")


# ══════════════════════════════════════════════════════════════════════════════
# Fig 1 — Volcano Plots
# ══════════════════════════════════════════════════════════════════════════════
print("Fig 1 — Volcano plots …")

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.subplots_adjust(wspace=0.38)

for ax, (cname, clabel) in zip(axes, COMPARISONS.items()):
    df = data.get(cname)
    if df is None:
        ax.set_visible(False)
        continue

    # Colour by significance & direction
    colors = np.where(
        ~df["sig"], C_GREY,
        np.where(df["LFC"] < 0, C_DEPL, C_ENRI)
    )
    sizes = np.where(df["sig"], 55, 22)

    ax.scatter(df["LFC"], df["log10p"],
               c=colors, s=sizes, linewidths=0, alpha=0.85, zorder=2)

    # p-value threshold line
    ax.axhline(-np.log10(PVAL_THRESH), color="#555555", lw=0.8,
               ls="--", zorder=1, label=f"p = {PVAL_THRESH}")
    ax.axvline(0, color="#888888", lw=0.6, ls=":", zorder=1)

    # Gene labels for top hits
    sig_sub = df[df["sig"]].nsmallest(TOP_LABEL, "pval")
    texts = []
    for _, row in sig_sub.iterrows():
        t = ax.text(row["LFC"], row["log10p"], row["id"],
                    fontsize=7.5, ha="center", va="bottom",
                    color=C_DEPL if row["LFC"] < 0 else C_ENRI)
        texts.append(t)
    if HAS_ADJUSTTEXT and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5))

    n_depl = (df["sig"] & (df["LFC"] < 0)).sum()
    n_enri = (df["sig"] & (df["LFC"] > 0)).sum()
    ax.set_title(clabel, fontsize=10, fontweight="bold", pad=8)
    ax.set_xlabel("Log₂ Fold Change (LFC)")
    ax.set_ylabel("−log₁₀(p-value)")

    # Mini legend in corner
    ax.text(0.03, 0.97, f"Depleted n={n_depl}", transform=ax.transAxes,
            fontsize=8, color=C_DEPL, va="top")
    ax.text(0.03, 0.90, f"Enriched n={n_enri}", transform=ax.transAxes,
            fontsize=8, color=C_ENRI, va="top")

# Shared legend
patch_ns   = mpatches.Patch(color=C_GREY,  label="Not significant")
patch_dep  = mpatches.Patch(color=C_DEPL,  label="Depleted (essential)")
patch_enr  = mpatches.Patch(color=C_ENRI,  label="Enriched (toxic)")
fig.legend(handles=[patch_ns, patch_dep, patch_enr],
           loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.03),
           frameon=False, fontsize=9)
fig.suptitle("CRISPR Screen — Volcano Plots (p < 0.05)", y=1.01,
             fontsize=12, fontweight="bold")

save_fig(fig, "Fig1_volcano")


# ══════════════════════════════════════════════════════════════════════════════
# Fig 2 — Gene Rank Plots (LFC ranked)
# ══════════════════════════════════════════════════════════════════════════════
print("Fig 2 — Rank plots …")

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
fig.subplots_adjust(wspace=0.38)

for ax, (cname, clabel) in zip(axes, COMPARISONS.items()):
    df = data.get(cname)
    if df is None:
        ax.set_visible(False)
        continue

    ranked = df.sort_values("LFC").reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1)

    colors = np.where(
        ~ranked["sig"], C_GREY,
        np.where(ranked["LFC"] < 0, C_DEPL, C_ENRI)
    )
    sizes = np.where(ranked["sig"], 55, 22)

    ax.scatter(ranked["rank"], ranked["LFC"],
               c=colors, s=sizes, linewidths=0, alpha=0.85, zorder=2)
    ax.axhline(0, color="#888888", lw=0.6, ls=":", zorder=1)

    sig_sub = ranked[ranked["sig"]].copy()
    texts = []
    for _, row in sig_sub.iterrows():
        t = ax.text(row["rank"], row["LFC"], row["id"],
                    fontsize=7.5, ha="center", va="bottom",
                    color=C_DEPL if row["LFC"] < 0 else C_ENRI)
        texts.append(t)
    if HAS_ADJUSTTEXT and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5))

    ax.set_title(clabel, fontsize=10, fontweight="bold", pad=8)
    ax.set_xlabel("Gene Rank")
    ax.set_ylabel("Log₂ Fold Change")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))

fig.suptitle("CRISPR Screen — Gene Rank Plots", y=1.01,
             fontsize=12, fontweight="bold")
save_fig(fig, "Fig2_rank_lfc")


# ══════════════════════════════════════════════════════════════════════════════
# Fig 3 — Heatmap: LFC of significant genes across all comparisons
# ══════════════════════════════════════════════════════════════════════════════
print("Fig 3 — Heatmap …")

all_sig_genes = sig_df["gene"].unique().tolist()

if len(all_sig_genes) == 0:
    print("  [SKIP] No significant genes for heatmap.")
else:
    # Build LFC matrix: rows = genes, cols = comparisons
    heat_data = {}
    for cname, df in data.items():
        heat_data[COMPARISONS[cname]] = (
            df.set_index("id")["LFC"]
            .reindex(all_sig_genes)
            .fillna(0.0)
        )
    heat_df = pd.DataFrame(heat_data)

    # Sort rows by most negative LFC in any comparison
    heat_df["_min"] = heat_df.min(axis=1)
    heat_df = heat_df.sort_values("_min").drop(columns="_min")

    vmax = np.nanmax(np.abs(heat_df.values))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Custom diverging colormap (blue → white → red)
    cmap = LinearSegmentedColormap.from_list(
        "bwr_custom",
        ["#1A5276", "#AED6F1", "#FDFEFE", "#F1948A", "#A93226"],
        N=256
    )

    row_h = max(0.35, 5.0 / len(heat_df))
    fig_h = max(4, len(heat_df) * row_h + 1.5)
    fig, ax = plt.subplots(figsize=(6, fig_h))

    im = ax.imshow(heat_df.values, aspect="auto", cmap=cmap, norm=norm)

    ax.set_xticks(range(len(heat_df.columns)))
    ax.set_xticklabels(
        [textwrap.fill(c, 18) for c in heat_df.columns],
        fontsize=8.5, rotation=30, ha="right"
    )
    ax.set_yticks(range(len(heat_df.index)))
    ax.set_yticklabels(heat_df.index, fontsize=8.5)

    # Annotate cells
    for r in range(len(heat_df.index)):
        for c in range(len(heat_df.columns)):
            val = heat_df.values[r, c]
            if val != 0:
                ax.text(c, r, f"{val:.1f}", ha="center", va="center",
                        fontsize=7, color="white" if abs(val) > vmax * 0.55 else "black")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("Log₂ Fold Change", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title("Significant Genes — LFC Across Comparisons\n"
                 f"(p < {PVAL_THRESH}; zero = not testable in that comparison)",
                 fontsize=10, fontweight="bold", pad=10)
    ax.tick_params(top=True, bottom=False, labeltop=False)

    save_fig(fig, "Fig3_heatmap")


# ══════════════════════════════════════════════════════════════════════════════
# Fig 4 — Barplots: top genes per comparison (horizontal, ranked by significance)
# ══════════════════════════════════════════════════════════════════════════════
print("Fig 4 — Barplots …")

n_comps = len(COMPARISONS)
fig, axes = plt.subplots(n_comps, 1, figsize=(7, 3.2 * n_comps))
fig.subplots_adjust(hspace=0.55)

for ax, (cname, clabel) in zip(axes, COMPARISONS.items()):
    df = data.get(cname)
    if df is None:
        ax.set_visible(False)
        continue

    top = df.nsmallest(15, "pval").sort_values("pval", ascending=False)
    bar_colors = [C_DEPL if lfc < 0 else C_ENRI for lfc in top["LFC"]]

    bars = ax.barh(top["id"], -np.log10(top["pval"]),
                   color=bar_colors, edgecolor="none", height=0.7)

    # FDR annotation on bars
    for bar, (_, row) in zip(bars, top.iterrows()):
        ax.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height() / 2,
                f"p={row['pval']:.3f}",
                va="center", ha="left", fontsize=7.2, color="#444444")

    ax.axvline(-np.log10(PVAL_THRESH), color="#555555", lw=0.9, ls="--")
    ax.set_xlabel("−log₁₀(p-value)")
    ax.set_title(clabel, fontsize=10, fontweight="bold")
    ax.set_xlim(0, top["log10p"].max() * 1.35)
    ax.tick_params(axis="y", labelsize=8.5)

patch_dep = mpatches.Patch(color=C_DEPL, label="Depleted (essential)")
patch_enr = mpatches.Patch(color=C_ENRI, label="Enriched (toxic)")
fig.legend(handles=[patch_dep, patch_enr],
           loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.01),
           frameon=False, fontsize=9)
fig.suptitle("Top Genes by Significance per Comparison",
             fontsize=12, fontweight="bold", y=1.005)
save_fig(fig, "Fig4_barplot")


# ══════════════════════════════════════════════════════════════════════════════
# Fig 5 — Library QC: mapping rate, Gini index, per-sample read depth
# ══════════════════════════════════════════════════════════════════════════════
print("Fig 5 — Library QC …")

# Mapping stats from the MAGeCK count log (previously recorded in analysis_report.md)
qc_stats = pd.DataFrame({
    "sample":      ["R0_1","R0_2","R0_3","Rpost_1","Rpost_2","Rpost_3","Rpost_4",
                    "S0_1","S0_2","S0_3","Spost_1","Spost_2","Spost_3"],
    "group":       ["R0","R0","R0","Rpost","Rpost","Rpost","Rpost",
                    "S0","S0","S0","Spost","Spost","Spost"],
    "total_reads": [13398116, 13605458, 2881705, 11794790, 1842586, 8259689, 7594859,
                    9964941, 11701486, 1654363, 1193221, 2254295, 14150359],
    "mapped_reads":[1120975,  1318990,  208140,  983542,   174022,  782450,  104123,
                    94968,    590915,   143722,  9349,     76829,   844365],
    "zero_count":  [30047, 30072, 30029, 29906, 30089, 29932, 30038,
                    30148, 30146, 30146, 30186, 30140, 29998],
    "gini":        [0.997, 0.997, 0.996, 0.993, 0.998, 0.995, 0.997,
                    0.999, 0.999, 0.999, 1.000, 0.999, 0.997],
})
qc_stats["mapping_pct"] = qc_stats["mapped_reads"] / qc_stats["total_reads"] * 100
qc_stats["dropout_pct"] = qc_stats["zero_count"] / 30197 * 100
qc_colors = [COND_COLORS[g] for g in qc_stats["group"]]

fig = plt.figure(figsize=(13, 10))
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

# Panel A — Total reads
ax1 = fig.add_subplot(gs[0, 0])
ax1.bar(range(len(qc_stats)), qc_stats["total_reads"] / 1e6,
        color=qc_colors, edgecolor="none")
ax1.set_xticks(range(len(qc_stats)))
ax1.set_xticklabels(qc_stats["sample"], rotation=45, ha="right", fontsize=8)
ax1.set_ylabel("Total Reads (millions)")
ax1.set_title("A. Sequencing Depth", fontweight="bold")
ax1.axhline(5, color="#888888", lw=0.8, ls="--", label="5M guideline")
ax1.legend(fontsize=8)

# Panel B — Mapping rate
ax2 = fig.add_subplot(gs[0, 1])
ax2.bar(range(len(qc_stats)), qc_stats["mapping_pct"],
        color=qc_colors, edgecolor="none")
ax2.set_xticks(range(len(qc_stats)))
ax2.set_xticklabels(qc_stats["sample"], rotation=45, ha="right", fontsize=8)
ax2.set_ylabel("sgRNA Mapping Rate (%)")
ax2.set_title("B. Mapping Rate", fontweight="bold")
ax2.axhline(60, color="#E74C3C", lw=0.9, ls="--", label="Ideal ≥60%")
ax2.legend(fontsize=8)

# Panel C — Gini index
ax3 = fig.add_subplot(gs[1, 0])
ax3.bar(range(len(qc_stats)), qc_stats["gini"],
        color=qc_colors, edgecolor="none")
ax3.set_xticks(range(len(qc_stats)))
ax3.set_xticklabels(qc_stats["sample"], rotation=45, ha="right", fontsize=8)
ax3.set_ylabel("Gini Index")
ax3.set_ylim(0.95, 1.005)
ax3.set_title("C. Gini Index (read inequality)", fontweight="bold")
ax3.axhline(0.2, color="#2ECC71", lw=0.9, ls="--", label="Ideal <0.2")
ax3.axhline(0.997, color="#E74C3C", lw=0.9, ls=":", label="Our data ~0.997")
ax3.legend(fontsize=8)

# Panel D — sgRNA dropout
ax4 = fig.add_subplot(gs[1, 1])
ax4.bar(range(len(qc_stats)), qc_stats["dropout_pct"],
        color=qc_colors, edgecolor="none")
ax4.set_xticks(range(len(qc_stats)))
ax4.set_xticklabels(qc_stats["sample"], rotation=45, ha="right", fontsize=8)
ax4.set_ylabel("Zero-count sgRNAs (%)")
ax4.set_title("D. sgRNA Dropout Rate", fontweight="bold")
ax4.axhline(20, color="#2ECC71", lw=0.9, ls="--", label="Ideal <20%")
ax4.legend(fontsize=8)

# Shared condition legend
handles = [mpatches.Patch(color=c, label=g) for g, c in COND_COLORS.items()]
fig.legend(handles=handles, loc="lower center", ncol=4,
           bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=9,
           title="Condition", title_fontsize=9)
fig.suptitle("Library & Sequencing Quality Control Metrics",
             fontsize=12, fontweight="bold", y=1.005)
save_fig(fig, "Fig5_qc_metrics")


# ══════════════════════════════════════════════════════════════════════════════
# Fig 6 — Per-sgRNA count profiles for key genes
# ══════════════════════════════════════════════════════════════════════════════
print("Fig 6 — sgRNA count profiles …")

key_genes = ["CYP2E1", "ASL", "PLA2G4E", "CYB5A", "MTF1"]
key_genes = [g for g in key_genes if g in counts["Gene"].values]

if not key_genes:
    print("  [SKIP] Key genes not found in count table.")
else:
    n_genes = len(key_genes)
    fig, axes = plt.subplots(n_genes, 1, figsize=(10, 3.0 * n_genes))
    if n_genes == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=0.6)

    for ax, gene in zip(axes, key_genes):
        gene_rows = counts[counts["Gene"] == gene].copy()
        if gene_rows.empty:
            ax.set_visible(False)
            continue

        # Normalise: reads per million within each sample
        sample_totals = count_matrix.sum(axis=0)
        rpm = gene_rows[count_matrix.columns].div(sample_totals, axis=1) * 1e6

        # Melt for plotting
        melted = rpm.T.reset_index()
        melted.columns = ["sample"] + [f"sg{i+1}" for i in range(len(rpm))]
        melted["group"] = melted["sample"].str.extract(r"(R0|Rpost|S0|Spost)")[0]
        melted_long = melted.melt(id_vars=["sample", "group"],
                                   var_name="sgRNA", value_name="RPM")

        for i, (sg, grp) in enumerate(melted_long.groupby("sgRNA")):
            x_pos = list(range(len(grp)))
            ax.plot(x_pos, grp["RPM"].values, marker="o", ms=5,
                    label=sg, lw=1.2, zorder=3)

        ax.set_xticks(range(len(melted_long["sample"].unique())))
        ax.set_xticklabels(melted_long["sample"].unique(),
                           rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Reads per million")
        ax.set_title(f"{gene} — individual sgRNA counts",
                     fontsize=10, fontweight="bold")
        ax.legend(title="sgRNA", fontsize=8, title_fontsize=8, loc="upper right")
        ax.set_xlim(-0.5, len(melted_long["sample"].unique()) - 0.5)

        # Shade by condition group
        prev_g, prev_x = None, 0
        sample_order = list(dict.fromkeys(melted_long["sample"].tolist()))
        grp_order = [
            melted_long.loc[melted_long["sample"] == s, "group"].values[0]
            for s in sample_order
        ]
        for xi, grp in enumerate(grp_order):
            if grp != prev_g:
                if prev_g is not None:
                    ax.axvspan(prev_x - 0.5, xi - 0.5,
                               facecolor=COND_COLORS.get(prev_g, "#DDDDDD"),
                               alpha=0.12, zorder=1)
                    ax.text((prev_x + xi - 1) / 2, ax.get_ylim()[1] * 0.98,
                            prev_g, ha="center", va="top", fontsize=8,
                            color=COND_COLORS.get(prev_g, "#444444"),
                            fontweight="bold")
                prev_g, prev_x = grp, xi
        if prev_g:
            ax.axvspan(prev_x - 0.5, len(sample_order) - 0.5,
                       facecolor=COND_COLORS.get(prev_g, "#DDDDDD"),
                       alpha=0.12, zorder=1)
            ax.text((prev_x + len(sample_order) - 1) / 2, ax.get_ylim()[1] * 0.98,
                    prev_g, ha="center", va="top", fontsize=8,
                    color=COND_COLORS.get(prev_g, "#444444"),
                    fontweight="bold")

    fig.suptitle("sgRNA Count Profiles — Key Hit Genes",
                 fontsize=12, fontweight="bold", y=1.005)
    save_fig(fig, "Fig6_sgrna_profiles")


# ══════════════════════════════════════════════════════════════════════════════
# Fig 7 — Replicate correlations (Spearman) within each condition
# ══════════════════════════════════════════════════════════════════════════════
print("Fig 7 — Replicate correlations …")

# Use log1p-RPM for correlation
sample_totals = count_matrix.sum(axis=0)
log_rpm = np.log1p(count_matrix.div(sample_totals, axis=1) * 1e6)

# Only include pairs within the same condition
pairs = []
for grp, samples in sample_groups.items():
    avail = [s for s in samples if s in log_rpm.columns]
    if len(avail) >= 2:
        for i in range(len(avail)):
            for j in range(i + 1, len(avail)):
                pairs.append((grp, avail[i], avail[j]))

n_pairs = len(pairs)
if n_pairs == 0:
    print("  [SKIP] No within-group replicate pairs.")
else:
    ncols = min(3, n_pairs)
    nrows = int(np.ceil(n_pairs / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.2 * nrows))
    axes = np.array(axes).flatten()

    for ax, (grp, s1, s2) in zip(axes, pairs):
        x = log_rpm[s1].values
        y = log_rpm[s2].values
        mask = (x > 0) | (y > 0)
        rho, pv = spearmanr(x[mask], y[mask])

        ax.scatter(x[mask], y[mask], s=6, color=COND_COLORS.get(grp, "#888888"),
                   alpha=0.5, linewidths=0)
        lim = max(x[mask].max(), y[mask].max()) * 1.05
        ax.plot([0, lim], [0, lim], "k--", lw=0.8, zorder=2)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel(f"{s1} log₁₊RPM", fontsize=8.5)
        ax.set_ylabel(f"{s2} log₁₊RPM", fontsize=8.5)
        ax.set_title(f"{grp}: {s1} vs {s2}", fontsize=9, fontweight="bold")
        ax.text(0.05, 0.93, f"ρ = {rho:.3f}", transform=ax.transAxes,
                fontsize=9, color=COND_COLORS.get(grp, "#444444"), fontweight="bold")

    # Hide unused axes
    for ax in axes[n_pairs:]:
        ax.set_visible(False)

    fig.suptitle("Replicate Correlation — log₁₊RPM (Spearman ρ)",
                 fontsize=12, fontweight="bold", y=1.005)
    fig.tight_layout()
    save_fig(fig, "Fig7_replicate_corr")


# ══════════════════════════════════════════════════════════════════════════════
# Fig 8 — Cross-comparison view: genes hitting in multiple comparisons
# ══════════════════════════════════════════════════════════════════════════════
print("Fig 8 — Cross-comparison gene overlap …")

# Build LFC dataframe for ALL genes across comparisons (for bubble plot)
all_genes_union = sorted(
    set().union(*[set(df["id"]) for df in data.values()])
)

lfc_all   = {}
pval_all  = {}
for cname, df in data.items():
    idx = df.set_index("id")
    lfc_all[cname]  = idx["LFC"].reindex(all_genes_union).fillna(np.nan)
    pval_all[cname] = idx["pval"].reindex(all_genes_union).fillna(np.nan)

lfc_mat  = pd.DataFrame(lfc_all,  index=all_genes_union)
pval_mat = pd.DataFrame(pval_all, index=all_genes_union)

# Genes significant in ≥1 comparison
n_sig = (pval_mat < PVAL_THRESH).sum(axis=1)
multi_genes = n_sig[n_sig >= 1].index.tolist()

if len(multi_genes) == 0:
    print("  [SKIP] No significant genes for cross-comparison plot.")
else:
    # Sort by number of hits, then by most negative LFC
    sort_key = n_sig[multi_genes].sort_values(ascending=False)
    gene_order = sort_key.index.tolist()[:30]   # max 30 genes

    lfc_sub  = lfc_mat.loc[gene_order]
    pval_sub = pval_mat.loc[gene_order]

    fig, ax = plt.subplots(figsize=(7, max(4, len(gene_order) * 0.42 + 1.5)))

    comp_labels = [COMPARISONS[c] for c in data.keys()]
    comp_names  = list(data.keys())

    for xi, cname in enumerate(comp_names):
        for yi, gene in enumerate(gene_order):
            lfc  = lfc_sub.loc[gene, cname]
            pv   = pval_sub.loc[gene, cname]
            if np.isnan(lfc):
                continue
            size  = max(20, -np.log10(pv + 1e-10) * 60)
            color = C_DEPL if lfc < 0 else C_ENRI
            alpha = 0.85 if pv < PVAL_THRESH else 0.25
            ax.scatter(xi, yi, s=size, c=color, alpha=alpha,
                       linewidths=0.4, edgecolors="white", zorder=3)
            if pv < PVAL_THRESH:
                ax.text(xi, yi, f"{lfc:.1f}", ha="center", va="center",
                        fontsize=6.5, color="white" if abs(lfc) > 5 else "black",
                        fontweight="bold")

    ax.set_xticks(range(len(comp_names)))
    ax.set_xticklabels(comp_labels, rotation=25, ha="right", fontsize=9)
    ax.set_yticks(range(len(gene_order)))
    ax.set_yticklabels(gene_order, fontsize=8.5)
    ax.set_xlim(-0.5, len(comp_names) - 0.5)
    ax.set_ylim(-0.5, len(gene_order) - 0.5)
    ax.grid(True, which="both", lw=0.4, color="#EEEEEE", zorder=1)
    ax.set_axisbelow(True)
    ax.set_title("Significant Genes Across Comparisons\n"
                 "(bubble size ∝ −log₁₀p; colour = direction; faded = n.s.)",
                 fontsize=10, fontweight="bold", pad=10)

    # Bubble size legend
    for pv_eg, label in [(0.05, "p=0.05"), (0.01, "p=0.01"), (0.001, "p=0.001")]:
        ax.scatter([], [], s=-np.log10(pv_eg) * 60, c="#888888",
                   alpha=0.7, label=label)
    ax.legend(title="Bubble size", title_fontsize=8, fontsize=8,
              loc="lower right", frameon=True, framealpha=0.85)

    fig.tight_layout()
    save_fig(fig, "Fig8_cross_comparison")


# ══════════════════════════════════════════════════════════════════════════════
# Fig 9 — Summary table
# ══════════════════════════════════════════════════════════════════════════════
print("Fig 9 — Summary table …")

FUNC_ANNOT = {
    "CYP2E1":    "CYP450 2E1 — drug metabolism, ROS production",
    "NAT2":      "N-acetyltransferase 2 — drug acetylation",
    "ASL":       "Argininosuccinate lyase — urea cycle",
    "GBA3":      "Beta-glucosidase — glucocerebrosidase",
    "TRPM7":     "Mg²⁺/Ca²⁺ TRPM7 channel — ion homeostasis",
    "ATP6V0A4":  "V-ATPase a4 subunit — lysosomal acidification",
    "TPCN1":     "Two-pore channel 1 — lysosomal Ca²⁺",
    "HSD11B1":   "11β-HSD1 — cortisol metabolism",
    "CYB5A":     "Cytochrome b5A — electron transfer, lipid desaturation",
    "ALOXE3":    "eLOX-3 lipoxygenase — lipid signalling",
    "MTF1":      "Metal-regulatory TF — stress response",
    "DHRSX":     "DHRSX dehydrogenase/reductase",
    "PLA2G4E":   "Group IVE PLA₂ — arachidonic acid release",
    "SLC1A1":    "EAAT3 glutamate transporter",
    "ARSD":      "Arylsulfatase D — lysosomal sulfatase",
    "GRIN1":     "NMDA receptor subunit 1 — glutamate signalling",
    "RRM2":      "Ribonucleotide reductase M2 — DNA synthesis",
    "PFKFB1":    "PFKFB1 — glycolysis regulator (liver isoform)",
    "SLC25A15":  "Ornithine carrier — mitochondrial urea cycle",
    "SLC6A3":    "DAT dopamine transporter",
    "SLC44A1":   "Choline transporter-like 1",
}

if not sig_df.empty:
    tbl = sig_df.copy()
    tbl["pval_str"] = tbl["pval"].apply(lambda v: f"{v:.4f}")
    tbl["lfc_str"]  = tbl["LFC"].apply(lambda v: f"{v:+.2f}")
    tbl["function"] = tbl["gene"].map(lambda g: FUNC_ANNOT.get(g, "—"))
    tbl = tbl.sort_values(["comparison", "pval"])

    col_labels  = ["Comparison", "Gene", "p-value", "LFC", "Direction", "Function"]
    col_keys    = ["comparison_label", "gene", "pval_str", "lfc_str", "direction", "function"]
    cell_data   = [[str(row[k])[:38] for k in col_keys] for _, row in tbl.iterrows()]
    col_widths  = [0.20, 0.08, 0.07, 0.06, 0.08, 0.38]

    row_h_pt  = 0.28
    fig_h = max(3.5, len(tbl) * row_h_pt + 1.2)
    fig, ax = plt.subplots(figsize=(14, fig_h))
    ax.axis("off")

    the_table = ax.table(
        cellText=cell_data,
        colLabels=col_labels,
        colWidths=col_widths,
        loc="center",
        cellLoc="left",
    )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(8)
    the_table.scale(1, 1.35)

    # Style header
    for j in range(len(col_labels)):
        cell = the_table[0, j]
        cell.set_facecolor("#1A5276")
        cell.set_text_props(color="white", fontweight="bold")

    # Stripe rows and colour by direction
    for i, (_, row) in enumerate(tbl.iterrows(), start=1):
        bg = "#EBF5FB" if i % 2 == 0 else "#FDFEFE"
        for j in range(len(col_labels)):
            cell = the_table[i, j]
            cell.set_facecolor(bg)
            if col_keys[j] == "direction":
                cell.set_text_props(
                    color=C_DEPL if row["direction"] == "Depleted" else C_ENRI,
                    fontweight="bold"
                )

    ax.set_title(
        f"Summary: All Significant Genes (p < {PVAL_THRESH}) — GBM TRAIL Metabolic CRISPR Screen\n"
        "Note: FDR > 0.25 for all genes owing to severe library bottlenecking. "
        "p-value used as primary criterion.",
        fontsize=9, fontweight="bold", pad=10
    )
    fig.tight_layout()
    save_fig(fig, "Fig9_summary_table")


# ══════════════════════════════════════════════════════════════════════════════
print("\nDone. All figures saved to:", OUT_DIR)
print("Files:")
for f in sorted(os.listdir(OUT_DIR)):
    print(f"  {f}")
