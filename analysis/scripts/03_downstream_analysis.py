#!/usr/bin/env python3
"""
03_downstream_analysis.py
=========================
Publication-quality downstream analysis for CRISPR-Cas9 screen results.

Generates 9 figures from MAGeCK gene_summary output:
  Fig 1 — Volcano plots (all 3 comparisons)
  Fig 2 — Gene rank plots (LFC-ranked)
  Fig 3 — Heatmap: LFC of significant genes across comparisons
  Fig 4 — Horizontal barplots: top significant genes per comparison
  Fig 5 — Library QC: read depth, zero-count fraction, count distribution
  Fig 6 — Bubble plot: gene significance across comparisons
  Fig 7 — Grouped bar: genes significant in ≥2 comparisons
  Fig 8 — Replicate correlation scatter plots
  Fig 9 — Formatted summary table of all significant genes

Input files expected:
  results/mageck_test/Spost_vs_S0.gene_summary.txt
  results/mageck_test/Rpost_vs_R0.gene_summary.txt
  results/mageck_test/Spost_vs_Rpost.gene_summary.txt
  results/denovo_count_table.txt   (or any MAGeCK-compatible count table)

Usage:
  python scripts/03_downstream_analysis.py
  python scripts/03_downstream_analysis.py --mageck-dir results/mageck_test --fdr 0.25

Key parameters:
  --fdr         FDR threshold for significance (default: 0.25)
  --top-label   Number of genes to label in volcano/rank plots (default: 15)
  --outdir      Output directory for figures (default: figures/)
"""

import argparse
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns
from scipy.stats import spearmanr
from adjustText import adjust_text

warnings.filterwarnings("ignore")

# ── Colour palette ─────────────────────────────────────────────────────────────
DEPL_COLOR = "#1565C0"   # blue  — depleted (essential)
ENRI_COLOR = "#C62828"   # red   — enriched (toxic / growth-suppressor)
GREY_COLOR = "#BDBDBD"   # grey  — not significant

COND_PALETTE = {
    "R0":    "#1565C0",
    "Rpost": "#EF5350",
    "S0":    "#66BB6A",
    "Spost": "#FFA726",
}

COMPARISONS = {
    "Spost_vs_S0":    "Sensitive: Post vs Day-0",
    "Rpost_vs_R0":    "Resistant: Post vs Day-0",
    "Spost_vs_Rpost": "Sensitive vs Resistant (Post)",
}


# ── Style ─────────────────────────────────────────────────────────────────────
def set_style():
    plt.rcParams.update({
        "figure.dpi":        300,
        "font.family":       "DejaVu Sans",
        "font.size":         11,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.linewidth":    1.2,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "legend.frameon":    False,
    })


# ── I/O helpers ────────────────────────────────────────────────────────────────
def load_gene_summary(path: str, fdr_thresh: float) -> pd.DataFrame:
    """Load and annotate a MAGeCK gene_summary file."""
    df = pd.read_csv(path, sep="\t")
    df["FDR"]      = df[["neg|fdr", "pos|fdr"]].min(axis=1)
    df["neg_FDR"]  = df["neg|fdr"]
    df["pos_FDR"]  = df["pos|fdr"]
    df["LFC"]      = df["neg|lfc"]
    df["log10FDR"] = -np.log10(df["FDR"].clip(lower=1e-10))
    df["sig"]      = df["FDR"] < fdr_thresh
    df["direction"] = np.where(df["neg|lfc"] < df["pos|lfc"], "depleted", "enriched")
    return df


def save(fig, fig_dir: str, name: str):
    for ext in ("pdf", "png"):
        p = os.path.join(fig_dir, f"{name}.{ext}")
        fig.savefig(p, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  Saved: {name}.pdf / .png")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — Volcano plots
# ══════════════════════════════════════════════════════════════════════════════
def fig_volcano(dfs: dict, fig_dir: str, fdr_thresh: float, top_n: int):
    print("Figure 1: Volcano plots …")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("CRISPR Screen – Volcano Plots (Gene Level)",
                 fontsize=14, fontweight="bold", y=1.01)

    for ax, (name, label) in zip(axes, COMPARISONS.items()):
        if name not in dfs:
            ax.set_title(f"{label}\n(no data)", fontweight="bold")
            continue
        df = dfs[name]
        colors = df.apply(
            lambda r: (DEPL_COLOR if r.direction == "depleted" else ENRI_COLOR)
                       if r.sig else GREY_COLOR, axis=1)

        ax.scatter(df["LFC"], df["log10FDR"], c=colors, s=30,
                   alpha=0.75, linewidths=0)
        ax.axhline(-np.log10(fdr_thresh), color="black", lw=1, ls="--", alpha=0.6)
        ax.axvline(0, color="black", lw=0.8, alpha=0.4)

        sig = df[df["sig"]].copy()
        top = sig.reindex(sig["LFC"].abs().sort_values(ascending=False).index).head(top_n)
        texts = [ax.text(r["LFC"], r["log10FDR"], r["id"], fontsize=8, fontweight="bold")
                 for _, r in top.iterrows()]
        if texts:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="grey", lw=0.6))

        n_dep = (sig["direction"] == "depleted").sum()
        n_enr = (sig["direction"] == "enriched").sum()
        patches = [
            mpatches.Patch(color=DEPL_COLOR, label=f"Depleted (n={n_dep})"),
            mpatches.Patch(color=ENRI_COLOR, label=f"Enriched (n={n_enr})"),
            mpatches.Patch(color=GREY_COLOR, label="Not significant"),
        ]
        ax.legend(handles=patches, fontsize=8, loc="upper right")
        ax.set_title(label, fontweight="bold", fontsize=11)
        ax.set_xlabel("Log₂ Fold Change", fontsize=10)
        ax.set_ylabel("−log₁₀(FDR)", fontsize=10)

    plt.tight_layout()
    save(fig, fig_dir, "Fig1_volcano_plots")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — Rank plots
# ══════════════════════════════════════════════════════════════════════════════
def fig_rank(dfs: dict, fig_dir: str, top_n: int):
    print("Figure 2: Rank plots …")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("CRISPR Screen – Gene Rank Plots",
                 fontsize=14, fontweight="bold", y=1.01)

    for ax, (name, label) in zip(axes, COMPARISONS.items()):
        if name not in dfs:
            ax.set_title(f"{label}\n(no data)", fontweight="bold")
            continue
        df = dfs[name].sort_values("LFC").reset_index(drop=True)
        df["rank"] = np.arange(len(df))
        colors = df.apply(
            lambda r: (DEPL_COLOR if r.direction == "depleted" else ENRI_COLOR)
                       if r.sig else GREY_COLOR, axis=1)

        ax.scatter(df["rank"], df["LFC"], c=colors, s=18,
                   alpha=0.75, linewidths=0)
        ax.axhline(0, color="black", lw=0.8, alpha=0.5)

        sig = df[df["sig"]]
        top_dep = sig[sig["direction"] == "depleted"].head(8)
        top_enr = sig[sig["direction"] == "enriched"].tail(8)
        texts = [ax.text(r["rank"], r["LFC"], r["id"], fontsize=8, fontweight="bold")
                 for _, r in pd.concat([top_dep, top_enr]).iterrows()]
        if texts:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="grey", lw=0.6))

        ax.set_title(label, fontweight="bold", fontsize=11)
        ax.set_xlabel("Gene Rank", fontsize=10)
        ax.set_ylabel("Log₂ Fold Change", fontsize=10)
        patches = [
            mpatches.Patch(color=DEPL_COLOR, label="Depleted (essential)"),
            mpatches.Patch(color=ENRI_COLOR, label="Enriched (toxic)"),
        ]
        ax.legend(handles=patches, fontsize=8)

    plt.tight_layout()
    save(fig, fig_dir, "Fig2_rank_plots")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 — Heatmap
# ══════════════════════════════════════════════════════════════════════════════
def fig_heatmap(dfs: dict, fig_dir: str, fdr_thresh: float):
    print("Figure 3: Heatmap …")
    sig_genes: set = set()
    idx = {}
    for name, df in dfs.items():
        idx[name] = df.set_index("id")
        sig_genes |= set(df[df["sig"]]["id"])

    if not sig_genes:
        print("  No significant genes — skipping heatmap.")
        return

    col_labels = {n: COMPARISONS[n] for n in dfs}
    lfc_mat = pd.DataFrame(0.0, index=sorted(sig_genes),
                           columns=list(col_labels.values()))
    fdr_mat = pd.DataFrame(1.0, index=sorted(sig_genes),
                           columns=list(col_labels.values()))

    for name, label in col_labels.items():
        for gene in sig_genes:
            if gene in idx[name].index:
                lfc_mat.loc[gene, label] = idx[name].loc[gene, "LFC"]
                fdr_mat.loc[gene, label] = idx[name].loc[gene, "FDR"]

    lfc_mat["_abs"] = lfc_mat.abs().max(axis=1)
    lfc_mat = lfc_mat.sort_values("_abs", ascending=False).drop(columns="_abs")
    fdr_mat = fdr_mat.loc[lfc_mat.index]

    vmax = max(lfc_mat.abs().max().max(), 0.1)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    fig_h = max(5, len(sig_genes) * 0.45 + 2)

    fig, ax = plt.subplots(figsize=(9, fig_h))
    im = ax.imshow(lfc_mat.values, cmap="RdBu_r", norm=norm, aspect="auto")

    for i, gene in enumerate(lfc_mat.index):
        for j, label in enumerate(lfc_mat.columns):
            fdr_val = fdr_mat.loc[gene, label]
            marker = ("***" if fdr_val < 0.01
                      else "**" if fdr_val < 0.05
                      else "*" if fdr_val < fdr_thresh else "")
            if marker:
                ax.text(j, i, marker, ha="center", va="center",
                        fontsize=9, color="white", fontweight="bold")

    ax.set_xticks(range(len(lfc_mat.columns)))
    ax.set_xticklabels(lfc_mat.columns, rotation=30, ha="right", fontsize=10)
    ax.set_yticks(range(len(lfc_mat.index)))
    ax.set_yticklabels(lfc_mat.index, fontsize=9)
    plt.colorbar(im, ax=ax, label="Log₂ Fold Change", shrink=0.7)
    ax.set_title(
        f"Significant Genes: LFC Across Comparisons\n"
        f"(FDR < {fdr_thresh}; * p<{fdr_thresh}, ** p<0.05, *** p<0.01)",
        fontsize=12, fontweight="bold", pad=12)

    plt.tight_layout()
    save(fig, fig_dir, "Fig3_heatmap_sig_genes")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 4 — Top gene barplots
# ══════════════════════════════════════════════════════════════════════════════
def fig_barplot(dfs: dict, fig_dir: str, fdr_thresh: float):
    print("Figure 4: Top gene barplots …")
    fig, axes = plt.subplots(3, 1, figsize=(10, 14))
    fig.suptitle("Top Significant Genes by Comparison",
                 fontsize=14, fontweight="bold")

    for ax, (name, label) in zip(axes, COMPARISONS.items()):
        if name not in dfs:
            ax.set_title(f"{label}\n(no data)", fontweight="bold")
            continue
        df = dfs[name]
        sig = df[df["sig"]].copy()

        if sig.empty:
            ax.text(0.5, 0.5, f"No significant genes (FDR < {fdr_thresh})",
                    ha="center", va="center", transform=ax.transAxes, fontsize=11)
            ax.set_title(label, fontweight="bold")
            continue

        sig = sig.sort_values("log10FDR", ascending=True)
        colors = [DEPL_COLOR if d == "depleted" else ENRI_COLOR for d in sig["direction"]]
        bars = ax.barh(sig["id"], sig["log10FDR"], color=colors,
                       edgecolor="none", height=0.7)

        for bar, (_, row) in zip(bars, sig.iterrows()):
            ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                    f"FDR={row['FDR']:.3f}", va="center", fontsize=8)

        ax.set_xlabel("−log₁₀(FDR)", fontsize=10)
        ax.set_title(label, fontweight="bold", fontsize=11)
        ax.axvline(-np.log10(fdr_thresh), color="grey", ls="--", lw=1)
        patches = [
            mpatches.Patch(color=DEPL_COLOR, label="Depleted (Essential)"),
            mpatches.Patch(color=ENRI_COLOR, label="Enriched (Toxic)"),
        ]
        ax.legend(handles=patches, fontsize=8, loc="lower right")

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, fig_dir, "Fig4_top_genes_barplot")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 5 — Library QC
# ══════════════════════════════════════════════════════════════════════════════
def fig_library_qc(count_file: str, fig_dir: str):
    print("Figure 5: Library QC …")
    if not os.path.exists(count_file):
        print(f"  Count file not found: {count_file} — skipping QC figure.")
        return

    df = pd.read_csv(count_file, sep="\t", index_col=0)
    sample_cols = [c for c in df.columns if c != "Gene"]

    # Condition mapping from sample name prefix
    def cond(s):
        for prefix in ["Rpost", "R0", "Spost", "S0"]:
            if s.startswith(prefix):
                return prefix
        return "?"

    conditions = {s: cond(s) for s in sample_cols}
    colors = [COND_PALETTE.get(conditions[s], "grey") for s in sample_cols]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle("Library QC – Count Distributions",
                 fontsize=14, fontweight="bold")

    # (a) Total read depth
    ax = axes[0]
    totals = df[sample_cols].sum()
    ax.bar(range(len(sample_cols)), totals / 1e6, color=colors, edgecolor="none")
    ax.set_xticks(range(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Total reads (M)", fontsize=10)
    ax.set_title("(a) Total Read Depth per Sample", fontweight="bold")
    patches = [mpatches.Patch(color=v, label=k) for k, v in COND_PALETTE.items()]
    ax.legend(handles=patches, fontsize=8)

    # (b) Zero-count fraction
    ax = axes[1]
    zero_frac = (df[sample_cols] == 0).mean()
    ax.bar(range(len(sample_cols)), zero_frac * 100, color=colors, edgecolor="none")
    ax.set_xticks(range(len(sample_cols)))
    ax.set_xticklabels(sample_cols, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Zero-count sgRNAs (%)", fontsize=10)
    ax.set_title("(b) Zero-Count sgRNAs per Sample", fontweight="bold")
    ax.axhline(50, color="red", ls="--", lw=1, alpha=0.7, label="50% threshold")
    ax.legend(fontsize=8)

    # (c) Count distribution violin
    ax = axes[2]
    plot_data = []
    for col in sample_cols:
        vals = df[col][df[col] > 0]
        for v in np.log10(vals + 1):
            plot_data.append({"Sample": col, "Log10(count+1)": v,
                              "Condition": conditions[col]})
    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        sns.violinplot(data=plot_df, x="Sample", y="Log10(count+1)",
                       hue="Condition", palette=COND_PALETTE,
                       ax=ax, inner="quartile", linewidth=0.8, legend=False)
        ax.set_xticklabels(sample_cols, rotation=60, ha="right", fontsize=8)
    ax.set_title("(c) Count Distribution (non-zero sgRNAs)", fontweight="bold")
    ax.set_xlabel("")

    plt.tight_layout()
    save(fig, fig_dir, "Fig5_library_QC")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 6 — Bubble plot
# ══════════════════════════════════════════════════════════════════════════════
def fig_bubble(dfs: dict, fig_dir: str, fdr_thresh: float):
    print("Figure 6: Bubble plot …")
    all_sig = set()
    idx = {}
    for name, df in dfs.items():
        idx[name] = df.set_index("id")
        all_sig |= set(df[df["sig"]]["id"])

    if not all_sig:
        print("  No significant genes — skipping bubble plot.")
        return

    records = []
    for name, label in COMPARISONS.items():
        if name not in dfs:
            continue
        for gene in all_sig:
            row = idx[name].loc[gene] if gene in idx[name].index else None
            records.append({
                "Gene": gene,
                "Comparison": label,
                "LFC":  row["LFC"] if row is not None else 0,
                "FDR":  row["FDR"] if row is not None else 1.0,
                "log10FDR": row["log10FDR"] if row is not None else 0,
                "sig":  row["sig"] if row is not None else False,
            })

    plot_df = pd.DataFrame(records)
    gene_order = sorted(all_sig, key=lambda g: max(
        abs(idx[n].loc[g, "LFC"]) if (n in idx and g in idx[n].index) else 0
        for n in dfs
    ), reverse=True)

    comp_labels = [COMPARISONS[n] for n in COMPARISONS if n in dfs]
    gene_idx = {g: i for i, g in enumerate(gene_order)}
    comp_idx = {c: i for i, c in enumerate(comp_labels)}

    fig, ax = plt.subplots(figsize=(12, max(6, len(all_sig) * 0.55 + 2)))
    for _, row in plot_df.iterrows():
        x = comp_idx.get(row["Comparison"], -1)
        y = gene_idx.get(row["Gene"], -1)
        if x < 0 or y < 0:
            continue
        ax.scatter(x, y,
                   s=max(row["log10FDR"] * 60, 10),
                   c=DEPL_COLOR if row["LFC"] < 0 else ENRI_COLOR,
                   alpha=0.85 if row["sig"] else 0.25,
                   linewidths=0.5, edgecolors="white")

    ax.set_xticks(range(len(comp_labels)))
    ax.set_xticklabels(comp_labels, rotation=25, ha="right", fontsize=10)
    ax.set_yticks(range(len(gene_order)))
    ax.set_yticklabels(gene_order, fontsize=9)
    ax.set_xlim(-0.5, len(comp_labels) - 0.5)
    ax.set_ylim(-0.5, len(gene_order) - 0.5)
    ax.grid(axis="x", lw=0.5, alpha=0.4)
    patches = [
        mpatches.Patch(color=DEPL_COLOR, label="Depleted (LFC < 0)"),
        mpatches.Patch(color=ENRI_COLOR, label="Enriched (LFC > 0)"),
    ]
    ax.legend(handles=patches, fontsize=9, loc="upper right")
    ax.set_title("Significant Genes Across Comparisons\n(bubble size = −log₁₀(FDR))",
                 fontsize=12, fontweight="bold")

    plt.tight_layout()
    save(fig, fig_dir, "Fig6_bubble_plot")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 7 — Multi-comparison grouped bar
# ══════════════════════════════════════════════════════════════════════════════
def fig_multi_comparison_bar(dfs: dict, fig_dir: str):
    print("Figure 7: Multi-comparison gene bars …")
    idx = {n: df.set_index("id") for n, df in dfs.items()}

    all_sig = set()
    for df in dfs.values():
        all_sig |= set(df[df["sig"]]["id"])

    multi = [g for g in all_sig if sum(
        g in idx[n].index and idx[n].loc[g, "sig"] for n in dfs
    ) >= 2]
    if not multi:
        multi = sorted(all_sig)
    if not multi:
        print("  No genes to plot — skipping.")
        return

    bar_colors = ["#1565C0", "#C62828", "#2E7D32"]
    x = np.arange(len(multi))
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(8, len(multi) * 1.2), 6))
    for i, (name, color) in enumerate(zip(list(COMPARISONS.keys()), bar_colors)):
        if name not in dfs:
            continue
        lfcs = [idx[name].loc[g, "LFC"] if g in idx[name].index else 0 for g in multi]
        ax.bar(x + i * width, lfcs, width, label=COMPARISONS[name],
               color=color, alpha=0.8, edgecolor="none")

    ax.set_xticks(x + width)
    ax.set_xticklabels(multi, rotation=40, ha="right", fontsize=9)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("Log₂ Fold Change", fontsize=11)
    ax.set_title("Genes Significant in ≥2 Comparisons", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)

    plt.tight_layout()
    save(fig, fig_dir, "Fig7_multi_comparison_genes")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 8 — Replicate correlation
# ══════════════════════════════════════════════════════════════════════════════
def fig_replicate_correlation(count_file: str, fig_dir: str):
    print("Figure 8: Replicate correlations …")
    if not os.path.exists(count_file):
        print(f"  Count file not found — skipping.")
        return

    df = pd.read_csv(count_file, sep="\t", index_col=0)
    groups = {
        "R0":    [c for c in df.columns if c.startswith("R0_")],
        "Rpost": [c for c in df.columns if c.startswith("Rpost_")],
        "S0":    [c for c in df.columns if c.startswith("S0_")],
        "Spost": [c for c in df.columns if c.startswith("Spost_")],
    }
    groups = {k: v for k, v in groups.items() if len(v) >= 2}

    n = len(groups)
    if n == 0:
        print("  No groups with ≥2 replicates — skipping.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle("Replicate Correlation (log₁₀ CPM)", fontsize=13, fontweight="bold")

    for ax, (cond, cols) in zip(axes.flat, groups.items()):
        sub = df[cols].astype(float)
        sub = sub.div(sub.sum()) * 1e6
        sub = np.log10(sub + 1)
        c1, c2 = cols[0], cols[1]
        r, _ = spearmanr(sub[c1], sub[c2])
        ax.scatter(sub[c1], sub[c2], s=8, alpha=0.4,
                   color=COND_PALETTE.get(cond, "grey"), linewidths=0)
        lims = [sub[[c1, c2]].min().min(), sub[[c1, c2]].max().max()]
        ax.plot(lims, lims, "k--", lw=0.8, alpha=0.6)
        ax.set_xlabel(f"{c1} (log₁₀ CPM)", fontsize=9)
        ax.set_ylabel(f"{c2} (log₁₀ CPM)", fontsize=9)
        ax.set_title(f"{cond} replicates", fontweight="bold")
        ax.text(0.05, 0.93, f"Spearman ρ = {r:.3f}", transform=ax.transAxes, fontsize=9)

    for ax in axes.flat[len(groups):]:
        ax.set_visible(False)

    plt.tight_layout()
    save(fig, fig_dir, "Fig8_replicate_correlations")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 9 — Summary table
# ══════════════════════════════════════════════════════════════════════════════
def fig_summary_table(dfs: dict, fig_dir: str, fdr_thresh: float):
    print("Figure 9: Summary table …")
    rows = []
    for name, label in COMPARISONS.items():
        if name not in dfs:
            continue
        for _, r in dfs[name][dfs[name]["sig"]].iterrows():
            rows.append({
                "Comparison": label,
                "Gene": r["id"],
                "FDR": round(r["FDR"], 4),
                "LFC": round(r["LFC"], 2),
                "Direction": r["direction"].capitalize(),
            })

    if not rows:
        print("  No significant genes — skipping summary table.")
        return

    table_df = pd.DataFrame(rows).sort_values(["Comparison", "FDR"])
    table_df.to_csv(os.path.join(fig_dir, "significant_genes.csv"), index=False)
    print("  Saved: significant_genes.csv")

    col_labels = list(table_df.columns)
    cell_colors = []
    for _, row in table_df.iterrows():
        c = "#E3F2FD" if row["Direction"] == "Depleted" else "#FFEBEE"
        cval = "#BBDEFB" if row["Direction"] == "Depleted" else "#FFCDD2"
        cell_colors.append([c, c, c, cval, c])

    fig_h = max(4, len(table_df) * 0.38 + 1.5)
    fig, ax = plt.subplots(figsize=(13, fig_h))
    ax.axis("off")

    tbl = ax.table(cellText=table_df.values, colLabels=col_labels,
                   cellColours=cell_colors, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.auto_set_column_width(col=list(range(len(col_labels))))
    for j in range(len(col_labels)):
        tbl[(0, j)].set_facecolor("#37474F")
        tbl[(0, j)].set_text_props(color="white", fontweight="bold")

    ax.set_title(f"Significant Genes (FDR < {fdr_thresh})",
                 fontsize=13, fontweight="bold", pad=10)
    plt.tight_layout()
    save(fig, fig_dir, "Fig9_summary_table")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
def parse_args():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mageck-dir", default=os.path.join(base, "results", "mageck_test"),
                   help="Directory containing .gene_summary.txt files")
    p.add_argument("--count-table",
                   default=os.path.join(base, "results", "denovo_count_table.txt"),
                   help="Count table for QC plots")
    p.add_argument("--outdir", default=os.path.join(base, "figures"),
                   help="Output directory for figures")
    p.add_argument("--fdr", type=float, default=0.25,
                   help="FDR significance threshold (default: 0.25)")
    p.add_argument("--top-label", type=int, default=15,
                   help="Max genes to label in volcano/rank plots (default: 15)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    set_style()

    print("=" * 65)
    print("CRISPR Screen Downstream Analysis")
    print(f"MAGeCK dir  : {args.mageck_dir}")
    print(f"FDR thresh  : {args.fdr}")
    print(f"Output dir  : {args.outdir}")
    print("=" * 65)

    # Load gene summaries
    dfs = {}
    for name in COMPARISONS:
        path = os.path.join(args.mageck_dir, f"{name}.gene_summary.txt")
        if os.path.exists(path):
            dfs[name] = load_gene_summary(path, args.fdr)
            n_sig = dfs[name]["sig"].sum()
            print(f"  Loaded {name}: {len(dfs[name])} genes, {n_sig} significant")
        else:
            print(f"  [WARN] Not found: {path}")

    if not dfs:
        print("ERROR: No gene_summary files found. Run 02_mageck_test.sh first.")
        return

    fig_volcano(dfs, args.outdir, args.fdr, args.top_label)
    fig_rank(dfs, args.outdir, args.top_label)
    fig_heatmap(dfs, args.outdir, args.fdr)
    fig_barplot(dfs, args.outdir, args.fdr)
    fig_library_qc(args.count_table, args.outdir)
    fig_bubble(dfs, args.outdir, args.fdr)
    fig_multi_comparison_bar(dfs, args.outdir)
    fig_replicate_correlation(args.count_table, args.outdir)
    fig_summary_table(dfs, args.outdir, args.fdr)

    print(f"\nAll figures saved to: {args.outdir}")


if __name__ == "__main__":
    main()
