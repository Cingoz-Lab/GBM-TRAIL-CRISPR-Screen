#!/usr/bin/env python3
"""
04_generate_figures.py
======================
Generates publication-quality figures from MAGeCK results.

This script wraps 03_downstream_analysis.py with the correct paths
for this project and a p-value threshold suitable for sparse data.

Data quality note:
  The screen has severe library bottlenecking (~97% sgRNA dropout).
  Only 1,006 of 30,197 sgRNAs have reads. With only 18 testable genes
  in the Spost_vs_S0 comparison, FDR adjustment is too conservative.
  We therefore report raw p-value < 0.05 as the threshold in addition
  to the standard FDR < 0.25.

Usage:
  python scripts/04_generate_figures.py

Output (in figures/):
  Fig1_volcano_plots.pdf/png
  Fig2_rank_plots.pdf/png
  Fig3_heatmap_sig_genes.pdf/png      (if genes pass threshold)
  Fig4_top_genes_barplot.pdf/png
  Fig5_library_QC.pdf/png
  Fig6_bubble_plot.pdf/png            (if genes pass threshold)
  Fig7_multi_comparison_genes.pdf/png (if ≥2 comparisons have hits)
  Fig8_replicate_correlations.pdf/png
  Fig9_summary_table.pdf/png          (if genes pass threshold)
  significant_genes.csv
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy.stats import spearmanr
from adjustText import adjust_text

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE        = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAGECK_DIR  = os.path.join(BASE, "mageck_test")
COUNT_FILE  = os.path.join(BASE, "mageck_test", "count_table_fixed.txt")
FIG_DIR     = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
FDR_THRESH    = 0.25    # standard FDR
PVAL_THRESH   = 0.05    # use p-value when FDR is too conservative
TOP_N_LABEL   = 12

COMPARISONS = {
    "Spost_vs_S0":    "Susceptible: Post vs Day-0",
    "Rpost_vs_R0":    "Resistant: Post vs Day-0",
    "Spost_vs_Rpost": "Susceptible vs Resistant (Post)",
}

DEPL = "#1565C0"
ENRI = "#C62828"
GREY = "#BDBDBD"
COND_PAL = {"R0": "#1565C0", "Rpost": "#EF5350", "S0": "#66BB6A", "Spost": "#FFA726"}

plt.rcParams.update({
    "figure.dpi": 300, "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.2, "legend.frameon": False,
})


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  ✓ {name}")


def load(name):
    path = os.path.join(MAGECK_DIR, f"{name}.gene_summary.txt")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    df["FDR"]      = df[["neg|fdr", "pos|fdr"]].min(axis=1)
    df["pval"]     = df[["neg|p-value", "pos|p-value"]].min(axis=1)
    df["LFC"]      = df["neg|lfc"]
    df["log10FDR"] = -np.log10(df["FDR"].clip(lower=1e-10))
    df["log10p"]   = -np.log10(df["pval"].clip(lower=1e-10))
    # Significance: use pval if very few genes (< 50), else FDR
    n = len(df)
    df["sig_fdr"]  = df["FDR"] < FDR_THRESH
    df["sig_pval"] = df["pval"] < PVAL_THRESH
    df["sig"]      = df["sig_pval"] if n < 50 else df["sig_fdr"]
    df["threshold_used"] = "p<0.05" if n < 50 else f"FDR<{FDR_THRESH}"
    df["direction"] = np.where(df["neg|lfc"] < df["pos|lfc"], "depleted", "enriched")
    return df


# ════════════════════════════════════════════════════════════════════════════
# Figure 1 — Volcano plots
# ════════════════════════════════════════════════════════════════════════════
def fig1_volcano(dfs):
    print("Fig 1: Volcano plots")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("CRISPR Screen — Volcano Plots (Gene Level)", fontsize=14, fontweight="bold", y=1.01)

    for ax, (name, label) in zip(axes, COMPARISONS.items()):
        df = dfs.get(name)
        if df is None:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(label, fontweight="bold")
            continue

        # Use -log10(pval) for y-axis (more informative with sparse data)
        colors = df.apply(lambda r: (DEPL if r.direction == "depleted" else ENRI) if r.sig else GREY, axis=1)
        ax.scatter(df["LFC"], df["log10p"], c=colors, s=35, alpha=0.78, linewidths=0)

        thresh_line = -np.log10(PVAL_THRESH if df["threshold_used"].iloc[0].startswith("p") else FDR_THRESH)
        ax.axhline(thresh_line, color="black", lw=1, ls="--", alpha=0.5)
        ax.axvline(0, color="black", lw=0.8, alpha=0.35)

        sig = df[df["sig"]]
        top = sig.reindex(sig["LFC"].abs().sort_values(ascending=False).index).head(TOP_N_LABEL)
        texts = [ax.text(r["LFC"], r["log10p"], r["id"], fontsize=8, fontweight="bold")
                 for _, r in top.iterrows()]
        if texts:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.6))

        thresh_label = df["threshold_used"].iloc[0]
        n_dep = (sig["direction"] == "depleted").sum()
        n_enr = (sig["direction"] == "enriched").sum()
        patches = [
            mpatches.Patch(color=DEPL, label=f"Depleted n={n_dep}"),
            mpatches.Patch(color=ENRI, label=f"Enriched n={n_enr}"),
            mpatches.Patch(color=GREY, label="Not sig."),
        ]
        ax.legend(handles=patches, fontsize=8, loc="upper right")
        ax.set_title(f"{label}\n({thresh_label})", fontweight="bold", fontsize=10)
        ax.set_xlabel("Log₂ Fold Change", fontsize=10)
        ax.set_ylabel("−log₁₀(p-value)", fontsize=10)

    plt.tight_layout()
    save(fig, "Fig1_volcano_plots")


# ════════════════════════════════════════════════════════════════════════════
# Figure 2 — Rank plots
# ════════════════════════════════════════════════════════════════════════════
def fig2_rank(dfs):
    print("Fig 2: Rank plots")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("CRISPR Screen — Gene Rank Plots", fontsize=14, fontweight="bold", y=1.01)

    for ax, (name, label) in zip(axes, COMPARISONS.items()):
        df = dfs.get(name)
        if df is None:
            ax.set_title(label, fontweight="bold")
            continue

        df_s = df.sort_values("LFC").reset_index(drop=True)
        df_s["rank"] = np.arange(len(df_s))
        colors = df_s.apply(
            lambda r: (DEPL if r.direction == "depleted" else ENRI) if r.sig else GREY, axis=1)
        ax.scatter(df_s["rank"], df_s["LFC"], c=colors, s=20, alpha=0.78, linewidths=0)
        ax.axhline(0, color="black", lw=0.8, alpha=0.5)

        sig = df_s[df_s["sig"]]
        top_dep = sig[sig["direction"] == "depleted"].head(8)
        top_enr = sig[sig["direction"] == "enriched"].tail(8)
        texts = [ax.text(r["rank"], r["LFC"], r["id"], fontsize=8, fontweight="bold")
                 for _, r in pd.concat([top_dep, top_enr]).iterrows()]
        if texts:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.6))

        ax.set_title(label, fontweight="bold", fontsize=11)
        ax.set_xlabel("Gene Rank", fontsize=10)
        ax.set_ylabel("Log₂ Fold Change", fontsize=10)
        patches = [mpatches.Patch(color=DEPL, label="Depleted"), mpatches.Patch(color=ENRI, label="Enriched")]
        ax.legend(handles=patches, fontsize=8)

    plt.tight_layout()
    save(fig, "Fig2_rank_plots")


# ════════════════════════════════════════════════════════════════════════════
# Figure 3 — Heatmap of significant genes
# ════════════════════════════════════════════════════════════════════════════
def fig3_heatmap(dfs):
    print("Fig 3: Heatmap")
    sig_genes = set()
    idx = {}
    for name, df in dfs.items():
        if df is not None:
            idx[name] = df.set_index("id")
            sig_genes |= set(df[df["sig"]]["id"])

    if not sig_genes:
        print("  No significant genes — skipping")
        return

    cols = {n: COMPARISONS[n] for n in COMPARISONS if n in dfs}
    lfc_mat = pd.DataFrame(0.0, index=sorted(sig_genes), columns=list(cols.values()))
    pval_mat = pd.DataFrame(1.0, index=sorted(sig_genes), columns=list(cols.values()))

    for name, label in cols.items():
        if name not in idx:
            continue
        for gene in sig_genes:
            if gene in idx[name].index:
                lfc_mat.loc[gene, label] = idx[name].loc[gene, "LFC"]
                pval_mat.loc[gene, label] = idx[name].loc[gene, "pval"]

    lfc_mat["_abs"] = lfc_mat.abs().max(axis=1)
    lfc_mat = lfc_mat.sort_values("_abs", ascending=False).drop(columns="_abs")
    pval_mat = pval_mat.loc[lfc_mat.index]

    vmax = max(lfc_mat.abs().max().max(), 0.5)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    fig_h = max(5, len(sig_genes) * 0.5 + 2.5)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    im = ax.imshow(lfc_mat.values, cmap="RdBu_r", norm=norm, aspect="auto")

    for i, gene in enumerate(lfc_mat.index):
        for j, label in enumerate(lfc_mat.columns):
            p = pval_mat.loc[gene, label]
            marker = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            if marker:
                ax.text(j, i, marker, ha="center", va="center",
                        fontsize=10, color="white", fontweight="bold")
            # LFC value annotation
            lfc_val = lfc_mat.loc[gene, label]
            if abs(lfc_val) > 0.5:
                ax.text(j, i + 0.3, f"{lfc_val:.1f}", ha="center", va="center",
                        fontsize=6.5, color="white", alpha=0.8)

    ax.set_xticks(range(len(lfc_mat.columns)))
    ax.set_xticklabels(lfc_mat.columns, rotation=30, ha="right", fontsize=10)
    ax.set_yticks(range(len(lfc_mat.index)))
    ax.set_yticklabels(lfc_mat.index, fontsize=9)
    plt.colorbar(im, ax=ax, label="Log₂ Fold Change", shrink=0.6, pad=0.02)
    ax.set_title("LFC Heatmap — Significant Genes Across Comparisons\n(* p<0.05, ** p<0.01, *** p<0.001)",
                 fontsize=12, fontweight="bold", pad=12)
    plt.tight_layout()
    save(fig, "Fig3_heatmap_sig_genes")


# ════════════════════════════════════════════════════════════════════════════
# Figure 4 — Top genes barplot
# ════════════════════════════════════════════════════════════════════════════
def fig4_barplot(dfs):
    print("Fig 4: Top gene barplots")
    fig, axes = plt.subplots(3, 1, figsize=(10, 15))
    fig.suptitle("Top Genes by Comparison", fontsize=14, fontweight="bold")

    for ax, (name, label) in zip(axes, COMPARISONS.items()):
        df = dfs.get(name)
        if df is None:
            ax.set_title(label, fontweight="bold")
            continue

        # Show top 20 by p-value regardless of significance
        top = df.sort_values("pval").head(20).sort_values("log10p", ascending=True)
        sig_mask = top["sig"].values
        colors = [DEPL if (d == "depleted" and s) else (ENRI if (d == "enriched" and s) else GREY)
                  for d, s in zip(top["direction"], sig_mask)]
        bars = ax.barh(top["id"], top["log10p"], color=colors, edgecolor="none", height=0.7)

        for bar, (_, row) in zip(bars, top.iterrows()):
            ax.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height() / 2,
                    f"p={row['pval']:.3f}", va="center", fontsize=7.5)

        thresh = -np.log10(PVAL_THRESH)
        ax.axvline(thresh, color="black", ls="--", lw=1, label=f"p={PVAL_THRESH}")
        ax.set_xlabel("−log₁₀(p-value)", fontsize=10)
        ax.set_title(label, fontweight="bold", fontsize=11)
        patches = [
            mpatches.Patch(color=DEPL, label="Depleted (sig.)"),
            mpatches.Patch(color=ENRI, label="Enriched (sig.)"),
            mpatches.Patch(color=GREY, label="Not significant"),
        ]
        ax.legend(handles=patches, fontsize=8, loc="lower right")

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, "Fig4_top_genes_barplot")


# ════════════════════════════════════════════════════════════════════════════
# Figure 5 — Library QC
# ════════════════════════════════════════════════════════════════════════════
def fig5_library_qc():
    print("Fig 5: Library QC")
    if not os.path.exists(COUNT_FILE):
        print(f"  Count file not found: {COUNT_FILE}")
        return

    df = pd.read_csv(COUNT_FILE, sep="\t", index_col=0)
    sample_cols = [c for c in df.columns if c != "Gene"]

    def cond(s):
        for p in ["Rpost", "R0", "Spost", "S0"]:
            if s.startswith(p): return p
        return "?"

    conditions = {s: cond(s) for s in sample_cols}
    colors = [COND_PAL.get(conditions[s], "grey") for s in sample_cols]

    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

    # (a) Total read depth
    ax = fig.add_subplot(gs[0, 0])
    totals = df[sample_cols].sum()
    ax.bar(range(len(sample_cols)), totals / 1e6, color=colors, edgecolor="none")
    ax.set_xticks(range(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Total mapped reads (M)", fontsize=10)
    ax.set_title("(a) Total Read Depth", fontweight="bold")
    patches = [mpatches.Patch(color=v, label=k) for k, v in COND_PAL.items()]
    ax.legend(handles=patches, fontsize=7, ncol=2)

    # (b) Zero-count %
    ax = fig.add_subplot(gs[0, 1])
    zero_frac = (df[sample_cols] == 0).mean() * 100
    ax.bar(range(len(sample_cols)), zero_frac, color=colors, edgecolor="none")
    ax.axhline(90, color="red", ls="--", lw=1, label="90%", alpha=0.7)
    ax.set_xticks(range(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Zero-count sgRNAs (%)", fontsize=10)
    ax.set_title("(b) Library Dropout per Sample", fontweight="bold")
    ax.legend(fontsize=8)

    # (c) Library coverage (non-zero sgRNAs per sample)
    ax = fig.add_subplot(gs[0, 2])
    nonzero = (df[sample_cols] > 0).sum()
    ax.bar(range(len(sample_cols)), nonzero, color=colors, edgecolor="none")
    ax.set_xticks(range(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Non-zero sgRNAs", fontsize=10)
    ax.set_title("(c) Detected sgRNAs per Sample", fontweight="bold")
    ax.axhline(df[sample_cols].shape[0] * 0.5, color="grey", ls="--", lw=1, alpha=0.7, label="50% of library")
    ax.legend(fontsize=8)

    # (d) Log10 count distribution violin
    ax = fig.add_subplot(gs[1, :2])
    plot_data = []
    for col in sample_cols:
        vals = df[col][df[col] > 0]
        for v in np.log10(vals + 1):
            plot_data.append({"Sample": col, "log10(count+1)": v, "Condition": conditions[col]})
    if plot_data:
        pdf = pd.DataFrame(plot_data)
        sns.violinplot(data=pdf, x="Sample", y="log10(count+1)",
                       hue="Condition", palette=COND_PAL,
                       ax=ax, inner="box", linewidth=0.8, legend=False)
        ax.set_xticklabels(sample_cols, rotation=60, ha="right", fontsize=8)
    ax.set_title("(d) Count Distribution (non-zero sgRNAs)", fontweight="bold")
    ax.set_xlabel("")

    # (e) Gini index annotation
    ax = fig.add_subplot(gs[1, 2])
    # Calculate Gini for each sample
    def gini(x):
        x = np.sort(x[x > 0])
        if len(x) == 0: return 1.0
        n = len(x)
        idx = np.arange(1, n + 1)
        return (2 * np.sum(idx * x) / (n * np.sum(x))) - (n + 1) / n

    ginis = [gini(df[s].values) for s in sample_cols]
    ax.bar(range(len(sample_cols)), ginis, color=colors, edgecolor="none")
    ax.axhline(0.9, color="red", ls="--", lw=1, label="Gini=0.9 (poor)", alpha=0.7)
    ax.axhline(0.5, color="green", ls="--", lw=1, label="Gini=0.5 (good)", alpha=0.7)
    ax.set_xticks(range(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=60, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Gini Index", fontsize=10)
    ax.set_title("(e) Gini Index (library diversity)", fontweight="bold")
    ax.legend(fontsize=8)

    fig.suptitle("Library Quality Control Metrics\n"
                 "(Gini index near 1.0 = poor library diversity; ideal < 0.2 for pooled screens)",
                 fontsize=12, fontweight="bold")
    save(fig, "Fig5_library_QC")


# ════════════════════════════════════════════════════════════════════════════
# Figure 6 — Bubble plot
# ════════════════════════════════════════════════════════════════════════════
def fig6_bubble(dfs):
    print("Fig 6: Bubble plot")
    sig_genes = set()
    idx = {}
    for name, df in dfs.items():
        if df is not None:
            idx[name] = df.set_index("id")
            sig_genes |= set(df[df["sig"]]["id"])

    if not sig_genes:
        print("  No significant genes — skipping")
        return

    records = []
    for name, label in COMPARISONS.items():
        if name not in dfs or dfs[name] is None:
            continue
        for gene in sig_genes:
            row = idx[name].loc[gene] if gene in idx[name].index else None
            records.append({
                "Gene": gene, "Comparison": label,
                "LFC": row["LFC"] if row is not None else 0,
                "log10p": row["log10p"] if row is not None else 0,
                "sig": row["sig"] if row is not None else False,
            })

    pdf = pd.DataFrame(records)
    gene_order = sorted(sig_genes, key=lambda g: max(
        abs(idx[n].loc[g, "LFC"]) if (n in idx and g in idx[n].index) else 0
        for n in dfs if dfs[n] is not None
    ), reverse=True)

    comp_labels = [COMPARISONS[n] for n in COMPARISONS if n in dfs and dfs[n] is not None]
    gi = {g: i for i, g in enumerate(gene_order)}
    ci = {c: i for i, c in enumerate(comp_labels)}

    fig, ax = plt.subplots(figsize=(12, max(6, len(gene_order) * 0.6 + 2)))
    for _, row in pdf.iterrows():
        x = ci.get(row["Comparison"], -1)
        y = gi.get(row["Gene"], -1)
        if x < 0 or y < 0: continue
        ax.scatter(x, y, s=max(row["log10p"] * 80, 12),
                   c=DEPL if row["LFC"] < 0 else ENRI,
                   alpha=0.85 if row["sig"] else 0.2,
                   linewidths=0.4, edgecolors="white")

    ax.set_xticks(range(len(comp_labels)))
    ax.set_xticklabels(comp_labels, rotation=25, ha="right", fontsize=10)
    ax.set_yticks(range(len(gene_order)))
    ax.set_yticklabels(gene_order, fontsize=9)
    ax.set_xlim(-0.5, len(comp_labels) - 0.5)
    ax.set_ylim(-0.5, len(gene_order) - 0.5)
    ax.grid(axis="x", lw=0.5, alpha=0.4)
    patches = [mpatches.Patch(color=DEPL, label="Depleted"), mpatches.Patch(color=ENRI, label="Enriched")]
    ax.legend(handles=patches, fontsize=9)
    ax.set_title("Significant Genes Across Comparisons\n(bubble size ∝ −log₁₀(p-value))",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    save(fig, "Fig6_bubble_plot")


# ════════════════════════════════════════════════════════════════════════════
# Figure 7 — sgRNA-level counts for key genes
# ════════════════════════════════════════════════════════════════════════════
def fig7_sgrna_counts():
    print("Fig 7: sgRNA-level counts for key genes")
    if not os.path.exists(COUNT_FILE):
        return

    count_df = pd.read_csv(COUNT_FILE, sep="\t", index_col=0)
    sample_cols = [c for c in count_df.columns if c != "Gene"]
    conditions = {s: s.rsplit("_", 1)[0] for s in sample_cols}

    # Key genes of interest from the screen
    genes_of_interest = ["CYP2E1", "PTGS1", "CYB5A", "TRPM7", "MTF1",
                         "SLC25A4", "ACSM2B", "ASL", "PLA2G4E"]
    gene_col = count_df["Gene"]

    fig, axes = plt.subplots(3, 3, figsize=(16, 14))
    fig.suptitle("sgRNA-Level Count Profiles for Key Genes\n(points = individual sgRNAs)",
                 fontsize=13, fontweight="bold")

    for ax, gene in zip(axes.flat, genes_of_interest):
        gene_sgrnas = count_df[gene_col == gene]
        if gene_sgrnas.empty:
            ax.set_title(gene, fontweight="bold")
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue

        # CPM normalise within each sample
        sample_totals = count_df[sample_cols].sum()
        cpm = gene_sgrnas[sample_cols].div(sample_totals) * 1e6

        for i, (sgrna_id, row) in enumerate(cpm.iterrows()):
            jitter = np.random.uniform(-0.15, 0.15, len(sample_cols))
            ax.scatter(np.arange(len(sample_cols)) + jitter, row.values,
                       s=25, alpha=0.7, zorder=3)

        # Group means
        for j, col in enumerate(sample_cols):
            mean_val = cpm[col].mean()
            ax.plot([j - 0.3, j + 0.3], [mean_val, mean_val], "k-", lw=2, zorder=4)

        ax.set_xticks(range(len(sample_cols)))
        ax.set_xticklabels(sample_cols, rotation=60, ha="right", fontsize=7)
        ax.set_ylabel("CPM", fontsize=9)
        ax.set_title(gene, fontweight="bold", fontsize=11)

        # Add condition separators
        cond_boundaries = [0, 3, 7, 10, 13]
        for b in cond_boundaries[1:-1]:
            ax.axvline(b - 0.5, color="grey", lw=0.8, ls=":", alpha=0.6)

    for ax in axes.flat[len(genes_of_interest):]:
        ax.set_visible(False)

    plt.tight_layout()
    save(fig, "Fig7_sgrna_counts_key_genes")


# ════════════════════════════════════════════════════════════════════════════
# Figure 8 — Replicate correlations
# ════════════════════════════════════════════════════════════════════════════
def fig8_replicate_corr():
    print("Fig 8: Replicate correlations")
    if not os.path.exists(COUNT_FILE):
        return

    df = pd.read_csv(COUNT_FILE, sep="\t", index_col=0)
    groups = {
        "R0":    [c for c in df.columns if c.startswith("R0_")],
        "Rpost": [c for c in df.columns if c.startswith("Rpost_")],
        "S0":    [c for c in df.columns if c.startswith("S0_")],
        "Spost": [c for c in df.columns if c.startswith("Spost_")],
    }
    groups = {k: v for k, v in groups.items() if len(v) >= 2}

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle("Replicate Correlation (log₁₀ CPM)", fontsize=13, fontweight="bold")

    for ax, (cond, cols) in zip(axes.flat, groups.items()):
        sub = df[cols].astype(float)
        sub = (sub.div(sub.sum()) * 1e6)
        sub_log = np.log10(sub + 1)
        c1, c2 = cols[0], cols[1]
        r, p = spearmanr(sub_log[c1], sub_log[c2])
        ax.scatter(sub_log[c1], sub_log[c2], s=8, alpha=0.35,
                   color=COND_PAL.get(cond, "grey"), linewidths=0)
        lims = [sub_log[[c1, c2]].min().min(), sub_log[[c1, c2]].max().max()]
        ax.plot(lims, lims, "k--", lw=0.8, alpha=0.6)
        ax.set_xlabel(f"{c1} log₁₀(CPM+1)", fontsize=9)
        ax.set_ylabel(f"{c2} log₁₀(CPM+1)", fontsize=9)
        ax.set_title(f"{cond} replicates  (ρ = {r:.3f}, p = {p:.3f})", fontweight="bold")

        # Add note if correlation is poor
        if r < 0.5:
            ax.text(0.05, 0.05, "⚠ Low replicate correlation\n(library bottlenecking)",
                    transform=ax.transAxes, fontsize=8, color="red",
                    va="bottom", style="italic")

    for ax in axes.flat[len(groups):]:
        ax.set_visible(False)

    plt.tight_layout()
    save(fig, "Fig8_replicate_correlations")


# ════════════════════════════════════════════════════════════════════════════
# Figure 9 — Summary table
# ════════════════════════════════════════════════════════════════════════════
def fig9_summary_table(dfs):
    print("Fig 9: Summary table")
    rows = []
    for name, label in COMPARISONS.items():
        df = dfs.get(name)
        if df is None: continue
        for _, r in df[df["sig"]].iterrows():
            rows.append({
                "Comparison": label,
                "Gene": r["id"],
                "p-value": round(r["pval"], 4),
                "FDR": round(r["FDR"], 3),
                "LFC": round(r["LFC"], 2),
                "Direction": r["direction"].capitalize(),
            })

    if not rows:
        print("  No significant genes — skipping")
        return

    table_df = pd.DataFrame(rows).sort_values(["Comparison", "p-value"])
    table_df.to_csv(os.path.join(FIG_DIR, "significant_genes.csv"), index=False)
    print("  Saved: significant_genes.csv")

    col_labels = list(table_df.columns)
    cell_colors = []
    for _, row in table_df.iterrows():
        c = "#E3F2FD" if row["Direction"] == "Depleted" else "#FFEBEE"
        cval = "#BBDEFB" if row["Direction"] == "Depleted" else "#FFCDD2"
        cell_colors.append([c, c, c, c, cval, c])

    fig_h = max(4, len(table_df) * 0.38 + 1.8)
    fig, ax = plt.subplots(figsize=(14, fig_h))
    ax.axis("off")
    tbl = ax.table(cellText=table_df.values, colLabels=col_labels,
                   cellColours=cell_colors, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.auto_set_column_width(col=list(range(len(col_labels))))
    for j in range(len(col_labels)):
        tbl[(0, j)].set_facecolor("#37474F")
        tbl[(0, j)].set_text_props(color="white", fontweight="bold")

    ax.set_title("Significant Genes Summary\n(threshold: p<0.05 for sparse data; FDR<0.25 for full data)",
                 fontsize=12, fontweight="bold", pad=12)
    plt.tight_layout()
    save(fig, "Fig9_summary_table")


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 65)
    print("CRISPR Screen — Figure Generation")
    print(f"MAGeCK dir : {MAGECK_DIR}")
    print(f"Count file : {COUNT_FILE}")
    print(f"Output dir : {FIG_DIR}")
    print("=" * 65)

    # Load all gene summaries
    dfs = {}
    for name in COMPARISONS:
        df = load(name)
        if df is not None:
            dfs[name] = df
            n_genes = len(df)
            n_sig = df["sig"].sum()
            thresh = df["threshold_used"].iloc[0]
            print(f"  {name}: {n_genes} genes, {n_sig} significant ({thresh})")
        else:
            print(f"  [missing] {name}")

    print()
    np.random.seed(42)

    fig1_volcano(dfs)
    fig2_rank(dfs)
    fig3_heatmap(dfs)
    fig4_barplot(dfs)
    fig5_library_qc()
    fig6_bubble(dfs)
    fig7_sgrna_counts()
    fig8_replicate_corr()
    fig9_summary_table(dfs)

    print(f"\nAll figures saved to: {FIG_DIR}")
    print("\nData quality warning:")
    print("  Gini index ~0.997 and 97% library dropout indicate severe bottlenecking.")
    print("  Results are exploratory. Key hits (CYP2E1, ASL, PLA2G4E) are consistent")
    print("  with published findings but require validation.")


if __name__ == "__main__":
    main()
