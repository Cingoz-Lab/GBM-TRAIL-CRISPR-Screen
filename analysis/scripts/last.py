#!/usr/bin/env python3
"""
last.py — Final CRISPR Screen Analysis
=======================================
GBM TRAIL Metabolic CRISPR/Cas9 Screen | Cingöz Lab

This script fixes the primary cause of the low mapping rate (1–10%) that
was identified by read-level inspection:

  Root cause: ~34.6% of R1 reads contain the sgRNA in REVERSE-COMPLEMENT
  orientation. MAGeCK count only detects the FORWARD scaffold, so these
  reads are silently discarded. R2 for those same fragments contains the
  guide in FORWARD orientation and was never used.

  Fix: Concatenate R1 + R2 for each sample before MAGeCK count.
       This nearly doubles the detectable reads per sample.

Pipeline steps
--------------
  1. Concatenate R1 + R2 FASTQ per sample → combined/ (gzipped)
  2. MAGeCK count on combined FASTQs → results/last_count/
  3. Fix duplicate column names in count table
  4. MAGeCK test (3 comparisons) → results/last_test/
  5. Generate publication-quality figures → figures_last/

Licensing note
--------------
  All Python code is original. Visual approach informed by MAGeCKFlute
  (Bioconductor, GPL>=3) and CRISPR-pipeline tools in Support/, but no
  code is copied. Dependencies: matplotlib, seaborn, scipy, numpy, pandas,
  adjustText — all MIT/BSD/PSF.

Usage
-----
  conda activate mageck-env
  python scripts/last.py

  # Skip recounting if already done (uses existing results/last_count/):
  python scripts/last.py --skip-count
"""

import argparse
import os
import subprocess
import sys
import textwrap
import warnings
from pathlib import Path

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
from scipy.stats import spearmanr

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False
    print("[WARN] adjustText not installed. Gene labels may overlap.")
    print("       pip install adjustText")

warnings.filterwarnings("ignore")

# ── Argument parsing ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--skip-count", action="store_true",
                    help="Skip MAGeCK count (reuse existing results/last_count/)")
args = parser.parse_args()

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE      = Path(__file__).resolve().parent.parent   # analysis/
FASTQ_DIR = BASE.parent / "Fastq"
LIB_FILE  = BASE.parent / "nf-core_crisprseq" / "library_sabatini_metko.tsv"

COMB_DIR  = BASE / "results" / "last_combined_fastq"
COUNT_DIR = BASE / "results" / "last_count"
TEST_DIR  = BASE / "results" / "last_test"
FIG_DIR   = BASE / "figures_last"

for d in [COMB_DIR, COUNT_DIR, TEST_DIR, FIG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Sample map ────────────────────────────────────────────────────────────────
# Each sample: (label, condition, R1_glob, R2_glob)
# Edit these filenames to match your sequencing data
# Format: (sample_label, condition, R1_filename, R2_filename)
SAMPLES = [
    ("R0_1",    "R0",    "R0_1_R1.fq.gz",    "R0_1_R2.fq.gz"),
    ("R0_2",    "R0",    "R0_2_R1.fq.gz",    "R0_2_R2.fq.gz"),
    ("R0_3",    "R0",    "R0_3_R1.fq.gz",    "R0_3_R2.fq.gz"),
    ("Rpost_1", "Rpost", "Rpost_1_R1.fq.gz", "Rpost_1_R2.fq.gz"),
    ("Rpost_2", "Rpost", "Rpost_2_R1.fq.gz", "Rpost_2_R2.fq.gz"),
    ("Rpost_3", "Rpost", "Rpost_3_R1.fq.gz", "Rpost_3_R2.fq.gz"),
    ("Rpost_4", "Rpost", "Rpost_4_R1.fq.gz", "Rpost_4_R2.fq.gz"),
    ("S0_1",    "S0",    "S0_1_R1.fq.gz",    "S0_1_R2.fq.gz"),
    ("S0_2",    "S0",    "S0_2_R1.fq.gz",    "S0_2_R2.fq.gz"),
    ("S0_3",    "S0",    "S0_3_R1.fq.gz",    "S0_3_R2.fq.gz"),
    ("Spost_1", "Spost", "Spost_1_R1.fq.gz", "Spost_1_R2.fq.gz"),
    ("Spost_2", "Spost", "Spost_2_R1.fq.gz", "Spost_2_R2.fq.gz"),
    ("Spost_3", "Spost", "Spost_3_R1.fq.gz", "Spost_3_R2.fq.gz"),
]

COMPARISONS = {
    "Spost_vs_S0":    "Sensitive: Post vs Day-0",
    "Rpost_vs_R0":    "Resistant: Post vs Day-0",
    "Spost_vs_Rpost": "Sensitive vs Resistant (Post)",
}

COMP_GROUPS = {
    "Spost_vs_S0":    ("Spost_1,Spost_2,Spost_3", "S0_1,S0_2,S0_3"),
    "Rpost_vs_R0":    ("Rpost_1,Rpost_2,Rpost_3,Rpost_4", "R0_1,R0_2,R0_3"),
    "Spost_vs_Rpost": ("Spost_1,Spost_2,Spost_3", "Rpost_1,Rpost_2,Rpost_3,Rpost_4"),
}

# ── Analysis parameters ────────────────────────────────────────────────────────
PVAL_THRESH = 0.05
FDR_THRESH  = 0.25
TOP_LABEL   = 7       # max gene labels on any single panel

# ── Colour palette ─────────────────────────────────────────────────────────────
C_DEPL = "#1A5276"
C_ENRI = "#922B21"
C_GREY = "#D5D8DC"
C_GOLD = "#B7950B"
COND_COLORS = {"R0": "#2471A3", "Rpost": "#CB4335",
               "S0": "#1E8449", "Spost": "#CA6F1E"}

# ── Global plot style ──────────────────────────────────────────────────────────
sns.set_theme(style="ticks", font_scale=1.0)
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
    "axes.linewidth":      0.8,
    "axes.spines.top":     False,
    "axes.spines.right":   False,
    "xtick.direction":     "out",
    "ytick.direction":     "out",
    "figure.facecolor":    "white",
    "axes.facecolor":      "white",
    "grid.color":          "#E8E8E8",
    "grid.linewidth":      0.6,
})


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Concatenate R1 + R2 per sample
# ═══════════════════════════════════════════════════════════════════════════════

def run(cmd: str, desc: str = "") -> int:
    """Run a shell command, streaming output. Returns exit code."""
    print(f"  $ {cmd[:120]}{'…' if len(cmd)>120 else ''}")
    result = subprocess.run(cmd, shell=True, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in result.stdout.splitlines():
        if any(k in line for k in ("INFO", "WARN", "ERROR", "Mapped", "Total")):
            print(f"    {line}")
    if result.returncode != 0:
        print(f"  [ERROR] exit {result.returncode}")
        print(result.stdout[-2000:])
    return result.returncode


if not args.skip_count:
    print("\n" + "="*65)
    print("STEP 1 — Concatenate R1 + R2 per sample")
    print("="*65)
    print("Rationale: 34.6% of R1 reads carry the sgRNA in RC orientation.")
    print("           R2 of those reads has it in forward orientation.")
    print("           Combining both files gives MAGeCK access to all reads.\n")

    combined_paths = {}
    for label, cond, r1_name, r2_name in SAMPLES:
        r1 = FASTQ_DIR / r1_name
        r2 = FASTQ_DIR / r2_name
        out = COMB_DIR / f"{label}_combined.fq.gz"
        combined_paths[label] = out
        if out.exists():
            print(f"  [skip] {label} combined FASTQ already exists.")
            continue
        if not r1.exists() or not r2.exists():
            print(f"  [WARN] FASTQ not found for {label}: {r1_name} or {r2_name}")
            continue
        # cat gzip files — gzip concatenation is valid (zcat decompresses both)
        cmd = f"cat '{r1}' '{r2}' > '{out}'"
        print(f"  Combining {label} …")
        ret = run(cmd)
        if ret == 0:
            size_mb = out.stat().st_size / 1e6
            print(f"    → {out.name}  ({size_mb:.0f} MB)")

    # ── STEP 2: MAGeCK count ──────────────────────────────────────────────────
    print("\n" + "="*65)
    print("STEP 2 — MAGeCK count (R1+R2 combined)")
    print("="*65)

    count_out = COUNT_DIR / "last_count"
    label_str  = ",".join(s[0] for s in SAMPLES)
    fastq_str  = " ".join(f"'{COMB_DIR / (s[0] + '_combined.fq.gz')}'"
                          for s in SAMPLES
                          if (COMB_DIR / (s[0] + "_combined.fq.gz")).exists())

    if not fastq_str:
        print("  [ERROR] No combined FASTQs found. Check FASTQ_DIR path.")
        sys.exit(1)

    mageck_count_cmd = (
        f"mageck count"
        f" -l '{LIB_FILE}'"
        f" -n '{count_out}'"
        f" --sample-label '{label_str}'"
        f" --fastq {fastq_str}"
        f" --sgrna-len 20"
        f" --norm-method median"
    )
    print("  Running MAGeCK count …")
    run(mageck_count_cmd)

    # Fix duplicate column names from MAGeCK output
    raw_count = count_out.with_suffix(".count.txt")
    fixed_count = COUNT_DIR / "last_count_fixed.txt"
    if raw_count.exists():
        df = pd.read_csv(raw_count, sep="\t")
        # Rename duplicate columns (MAGeCK sometimes duplicates sample names)
        cols = list(df.columns)
        seen: dict = {}
        new_cols = []
        for c in cols:
            if c in ("sgRNA", "Gene"):
                new_cols.append(c)
            else:
                n = seen.get(c, 0)
                new_cols.append(f"{c}_{n+1}" if n > 0 else c)
                seen[c] = n + 1
        df.columns = new_cols
        df.to_csv(fixed_count, sep="\t", index=False)
        print(f"  → Count table: {fixed_count} ({len(df):,} sgRNAs)")
    else:
        print(f"  [WARN] Count file not found: {raw_count}")
        fixed_count = None

    # ── STEP 3: MAGeCK test ───────────────────────────────────────────────────
    if fixed_count and fixed_count.exists():
        print("\n" + "="*65)
        print("STEP 3 — MAGeCK test (RRA)")
        print("="*65)

        for comp, (treatment, control) in COMP_GROUPS.items():
            out_prefix = TEST_DIR / comp
            cmd = (
                f"mageck test"
                f" -k '{fixed_count}'"
                f" -t {treatment}"
                f" -c {control}"
                f" -n '{out_prefix}'"
                f" --gene-lfc-method median"
                f" --remove-zero both"
                f" --remove-zero-threshold 0"
            )
            print(f"  {comp} …")
            run(cmd)
            gs = out_prefix.with_suffix(".gene_summary.txt")
            if gs.exists():
                n = len(pd.read_csv(gs, sep="\t"))
                print(f"    → {n} testable genes")
else:
    print("\n[--skip-count] Using existing results in results/last_test/")
    fixed_count = COUNT_DIR / "last_count_fixed.txt"


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Load results (new or fallback to original mageck_test/)
# ═══════════════════════════════════════════════════════════════════════════════

def load(name: str, test_dir: Path) -> pd.DataFrame | None:
    path = test_dir / f"{name}.gene_summary.txt"
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    df["LFC"]      = df["neg|lfc"]
    df["pval"]     = df[["neg|p-value", "pos|p-value"]].min(axis=1)
    df["FDR"]      = df[["neg|fdr", "pos|fdr"]].min(axis=1)
    df["log10p"]   = -np.log10(df["pval"].clip(lower=1e-10))
    df["log10fdr"] = -np.log10(df["FDR"].clip(lower=1e-10))
    df["sig"]      = df["pval"] < PVAL_THRESH
    return df


print("\n" + "="*65)
print("STEP 4 — Loading gene summaries")
print("="*65)

data: dict[str, pd.DataFrame] = {}
for cname in COMPARISONS:
    # Prefer new last_test results; fall back to original mageck_test
    df = load(cname, TEST_DIR)
    if df is None:
        fallback = BASE / "mageck_test"
        df = load(cname, fallback)
        if df is not None:
            print(f"  {cname}: using original mageck_test/ results")
    if df is not None:
        data[cname] = df
        n_sig = df["sig"].sum()
        print(f"  {cname}: {len(df)} genes, {n_sig} significant (p<{PVAL_THRESH})")

if not data:
    print("[ERROR] No gene summaries found. Run without --skip-count first.")
    sys.exit(1)

# Load count table (new or original)
ct_path = fixed_count if (fixed_count and fixed_count.exists()) \
          else BASE / "mageck_test" / "count_table_fixed.txt"
counts_raw = pd.read_csv(ct_path, sep="\t", index_col=0)
counts = counts_raw[~counts_raw["Gene"].str.upper().str.contains("INTERGENIC", na=False)]
count_matrix = counts.drop(columns=["Gene"])

sample_groups = {
    "R0":    [c for c in count_matrix.columns if c.startswith("R0")],
    "Rpost": [c for c in count_matrix.columns if c.startswith("Rpost")],
    "S0":    [c for c in count_matrix.columns if c.startswith("S0")],
    "Spost": [c for c in count_matrix.columns if c.startswith("Spost")],
}

# Collect significant genes
sig_rows = []
for cname, df in data.items():
    for _, row in df[df["sig"]].iterrows():
        sig_rows.append({
            "comparison":       cname,
            "comparison_label": COMPARISONS[cname],
            "gene":  row["id"],
            "pval":  row["pval"],
            "FDR":   row["FDR"],
            "LFC":   row["LFC"],
            "direction": "Depleted" if row["LFC"] < 0 else "Enriched",
        })
sig_df = pd.DataFrame(sig_rows)
sig_df.to_csv(FIG_DIR / "significant_genes_last.csv", index=False)
print(f"\n  Total significant hits: {len(sig_df)}")


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def save_fig(fig: plt.Figure, name: str) -> None:
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  ✓  {name}")


def label_genes(ax: plt.Axes, sub: pd.DataFrame,
                x_col: str, y_col: str,
                n: int = TOP_LABEL) -> None:
    """Add non-overlapping gene labels using adjustText."""
    top = sub.nsmallest(n, "pval")
    texts = []
    for _, row in top.iterrows():
        color = C_DEPL if row["LFC"] < 0 else C_ENRI
        t = ax.text(row[x_col], row[y_col], row["id"],
                    fontsize=7.5, color=color, fontweight="bold",
                    ha="center", va="bottom", zorder=5)
        texts.append(t)
    if HAS_ADJUSTTEXT and texts:
        adjust_text(
            texts, ax=ax,
            expand_text=(1.8, 1.8),
            expand_points=(2.0, 2.0),
            force_text=(0.5, 0.8),
            force_points=(0.3, 0.5),
            arrowprops=dict(arrowstyle="-", color="#888888",
                            lw=0.5, shrinkA=2, shrinkB=2),
        )


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Volcano plots
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*65)
print("STEP 5 — Generating figures")
print("="*65)
print("Fig 1 — Volcano plots …")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.subplots_adjust(wspace=0.42)

for ax, (cname, clabel) in zip(axes, COMPARISONS.items()):
    df = data.get(cname)
    if df is None:
        ax.set_visible(False)
        continue

    # Background points
    ns = df[~df["sig"]]
    ax.scatter(ns["LFC"], ns["log10p"], c=C_GREY, s=18,
               linewidths=0, alpha=0.6, zorder=1, rasterized=True)

    # Significant — depleted
    dep = df[df["sig"] & (df["LFC"] < 0)]
    ax.scatter(dep["LFC"], dep["log10p"], c=C_DEPL, s=50,
               linewidths=0.3, edgecolors="white", alpha=0.9, zorder=3)

    # Significant — enriched
    enr = df[df["sig"] & (df["LFC"] > 0)]
    ax.scatter(enr["LFC"], enr["log10p"], c=C_ENRI, s=50,
               linewidths=0.3, edgecolors="white", alpha=0.9, zorder=3)

    # Threshold lines
    ax.axhline(-np.log10(PVAL_THRESH), color="#555555",
               lw=0.9, ls="--", zorder=2)
    ax.axvline(0, color="#AAAAAA", lw=0.6, ls=":", zorder=1)

    # Gene labels — top significant only
    sig_sub = df[df["sig"]]
    if not sig_sub.empty:
        label_genes(ax, sig_sub, "LFC", "log10p", n=min(TOP_LABEL, len(sig_sub)))

    # Counters
    ax.text(0.03, 0.97, f"▼ {len(dep)} depleted",
            transform=ax.transAxes, fontsize=8, color=C_DEPL,
            va="top", fontweight="bold")
    ax.text(0.03, 0.90, f"▲ {len(enr)} enriched",
            transform=ax.transAxes, fontsize=8, color=C_ENRI, va="top")

    ax.set_title(clabel, fontsize=10, fontweight="bold", pad=10)
    ax.set_xlabel("Log₂ Fold Change", labelpad=6)
    ax.set_ylabel("−log₁₀(p-value)", labelpad=6)
    sns.despine(ax=ax)

legend_handles = [
    mpatches.Patch(color=C_GREY, label="Not significant"),
    mpatches.Patch(color=C_DEPL, label="Depleted (essential)"),
    mpatches.Patch(color=C_ENRI, label="Enriched (toxic/growth-promoting KO)"),
]
fig.legend(handles=legend_handles, loc="lower center", ncol=3,
           bbox_to_anchor=(0.5, -0.04), frameon=False, fontsize=9)
fig.suptitle("CRISPR Screen — Volcano Plots  (p < 0.05)",
             fontsize=12, fontweight="bold", y=1.02)
save_fig(fig, "Fig1_volcano")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Gene Rank Plots
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 2 — Rank plots …")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.subplots_adjust(wspace=0.42)

for ax, (cname, clabel) in zip(axes, COMPARISONS.items()):
    df = data.get(cname)
    if df is None:
        ax.set_visible(False)
        continue

    ranked = df.sort_values("LFC").reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1)

    # Grey background scatter
    ns_mask = ~ranked["sig"]
    ax.scatter(ranked.loc[ns_mask, "rank"], ranked.loc[ns_mask, "LFC"],
               c=C_GREY, s=16, linewidths=0, alpha=0.7, zorder=1, rasterized=True)

    dep_mask = ranked["sig"] & (ranked["LFC"] < 0)
    enr_mask = ranked["sig"] & (ranked["LFC"] > 0)
    ax.scatter(ranked.loc[dep_mask, "rank"], ranked.loc[dep_mask, "LFC"],
               c=C_DEPL, s=48, linewidths=0.3, edgecolors="white",
               alpha=0.9, zorder=3)
    ax.scatter(ranked.loc[enr_mask, "rank"], ranked.loc[enr_mask, "LFC"],
               c=C_ENRI, s=48, linewidths=0.3, edgecolors="white",
               alpha=0.9, zorder=3)

    ax.axhline(0, color="#AAAAAA", lw=0.7, ls=":", zorder=1)

    sig_sub = ranked[ranked["sig"]]
    if not sig_sub.empty:
        label_genes(ax, sig_sub, "rank", "LFC", n=min(TOP_LABEL, len(sig_sub)))

    ax.set_title(clabel, fontsize=10, fontweight="bold", pad=10)
    ax.set_xlabel("Gene Rank (by LFC)", labelpad=6)
    ax.set_ylabel("Log₂ Fold Change", labelpad=6)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    sns.despine(ax=ax)

fig.suptitle("CRISPR Screen — Gene Rank Plots",
             fontsize=12, fontweight="bold", y=1.02)
save_fig(fig, "Fig2_rank_lfc")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — LFC Heatmap
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 3 — Heatmap …")

all_sig_genes = sig_df["gene"].unique().tolist() if not sig_df.empty else []

if not all_sig_genes:
    print("  [SKIP] No significant genes.")
else:
    heat_data = {}
    for cname, df in data.items():
        heat_data[COMPARISONS[cname]] = (
            df.set_index("id")["LFC"]
            .reindex(all_sig_genes).fillna(0.0)
        )
    heat_df = pd.DataFrame(heat_data)
    heat_df = heat_df.sort_values(heat_df.min(axis=1).name
                                  if heat_df.min(axis=1).name else heat_df.columns[0])
    heat_df = heat_df.loc[heat_df.min(axis=1).sort_values().index]

    vmax = max(abs(heat_df.values[~np.isnan(heat_df.values)].max()),
               abs(heat_df.values[~np.isnan(heat_df.values)].min())) + 0.5
    vmax = min(vmax, 16)

    cmap = LinearSegmentedColormap.from_list(
        "bwr2",
        ["#154360", "#2E86C1", "#AED6F1", "#FDFEFE", "#F1948A", "#922B21", "#641E16"],
        N=256
    )
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    row_h  = max(0.32, 5.5 / max(len(heat_df), 1))
    fig_h  = max(4.5, len(heat_df) * row_h + 2.0)
    fig, ax = plt.subplots(figsize=(6.5, fig_h))

    im = ax.imshow(heat_df.values, aspect="auto", cmap=cmap, norm=norm)

    ax.set_xticks(range(len(heat_df.columns)))
    ax.set_xticklabels(
        [textwrap.fill(c, 20) for c in heat_df.columns],
        fontsize=9, rotation=30, ha="right"
    )
    ax.set_yticks(range(len(heat_df.index)))
    ax.set_yticklabels(heat_df.index, fontsize=8.5)

    for r in range(len(heat_df.index)):
        for c in range(len(heat_df.columns)):
            v = heat_df.values[r, c]
            if v != 0:
                txt_color = "white" if abs(v) > vmax * 0.5 else "#333333"
                ax.text(c, r, f"{v:.1f}", ha="center", va="center",
                        fontsize=7, color=txt_color)

    ax.tick_params(top=False, bottom=True, labelbottom=True)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.04, shrink=0.8)
    cbar.set_label("Log₂ Fold Change", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(
        "Significant Genes — LFC Across Comparisons\n"
        f"(p < {PVAL_THRESH}; 0 = not testable in that comparison)",
        fontsize=10, fontweight="bold", pad=12
    )
    fig.tight_layout()
    save_fig(fig, "Fig3_heatmap")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Barplots (top genes per comparison)
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 4 — Barplots …")

n_comps = len(data)
fig, axes = plt.subplots(n_comps, 1, figsize=(8, 3.6 * n_comps))
fig.subplots_adjust(hspace=0.65)
if n_comps == 1:
    axes = [axes]

for ax, (cname, clabel) in zip(axes, COMPARISONS.items()):
    df = data.get(cname)
    if df is None:
        ax.set_visible(False)
        continue

    top = df.nsmallest(15, "pval").sort_values("log10p")
    colors = [C_DEPL if lfc < 0 else C_ENRI for lfc in top["LFC"]]

    bars = ax.barh(range(len(top)), top["log10p"].values,
                   color=colors, edgecolor="none", height=0.65)

    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["id"].values, fontsize=8.5)
    ax.axvline(-np.log10(PVAL_THRESH), color="#555555", lw=0.9, ls="--")

    # p-value annotations
    for i, (bar, (_, row)) in enumerate(zip(bars, top.iterrows())):
        ax.text(bar.get_width() + 0.03, i,
                f"p={row['pval']:.4f}  LFC={row['LFC']:+.1f}",
                va="center", ha="left", fontsize=7.5, color="#444444")

    ax.set_xlim(0, top["log10p"].max() * 1.55)
    ax.set_xlabel("−log₁₀(p-value)", labelpad=6)
    ax.set_title(clabel, fontsize=10, fontweight="bold", pad=8)
    ax.grid(axis="x", lw=0.5, color="#EEEEEE")
    sns.despine(ax=ax, left=True)

patch_dep = mpatches.Patch(color=C_DEPL, label="Depleted (essential)")
patch_enr = mpatches.Patch(color=C_ENRI, label="Enriched (toxic)")
fig.legend(handles=[patch_dep, patch_enr], loc="lower center", ncol=2,
           bbox_to_anchor=(0.5, -0.01), frameon=False, fontsize=9)
fig.suptitle("Top Significant Genes per Comparison",
             fontsize=12, fontweight="bold", y=1.005)
save_fig(fig, "Fig4_barplot")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Library & Mapping QC
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 5 — QC metrics …")

qc = pd.DataFrame({
    "sample":       ["R0_1","R0_2","R0_3","Rpost_1","Rpost_2","Rpost_3","Rpost_4",
                     "S0_1","S0_2","S0_3","Spost_1","Spost_2","Spost_3"],
    "group":        ["R0","R0","R0","Rpost","Rpost","Rpost","Rpost",
                     "S0","S0","S0","Spost","Spost","Spost"],
    "total_M":      [13.40,13.61,2.88,11.79,1.84,8.26,7.59,
                     9.96,11.70,1.65,1.19,2.25,14.15],
    "mapped_M":     [1.12,1.32,0.21,0.98,0.17,0.78,0.10,
                     0.09,0.59,0.14,0.01,0.08,0.84],
    "gini":         [0.997,0.997,0.996,0.993,0.998,0.995,0.997,
                     0.999,0.999,0.999,1.000,0.999,0.997],
    "dropout_pct":  [99.50,99.58,99.44,98.97,99.64,99.78,99.47,
                     99.84,99.83,99.84,99.96,99.81,99.34],
})
qc["mapping_pct"] = qc["mapped_M"] / qc["total_M"] * 100
colors = [COND_COLORS[g] for g in qc["group"]]

fig = plt.figure(figsize=(14, 9))
gs2 = gridspec.GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.38)

panels = [
    (gs2[0, 0], qc["total_M"],      "Sequencing Depth",    "Total Reads (M)", 5,   "≥5M guideline", "#888888"),
    (gs2[0, 1], qc["mapping_pct"],  "sgRNA Mapping Rate",  "Mapping Rate (%)",60,  "Ideal ≥60%",    "#E74C3C"),
    (gs2[1, 0], qc["gini"],         "Gini Index",          "Gini Index",      0.2, "Ideal <0.2",    "#27AE60"),
    (gs2[1, 1], qc["dropout_pct"],  "sgRNA Dropout Rate",  "Zero-count (%)",  20,  "Ideal <20%",    "#27AE60"),
]

for i, (pos, vals, title, ylabel, ref_val, ref_label, ref_color) in enumerate(panels):
    ax = fig.add_subplot(pos)
    ax.bar(range(len(qc)), vals, color=colors, edgecolor="none", width=0.75)
    ax.axhline(ref_val, color=ref_color, lw=1.0, ls="--",
               label=ref_label, alpha=0.8)
    ax.set_xticks(range(len(qc)))
    ax.set_xticklabels(qc["sample"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(f"{'ABCD'[i]}. {title}", fontsize=10, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(axis="y", lw=0.5, color="#F0F0F0")
    sns.despine(ax=ax)

cond_handles = [mpatches.Patch(color=c, label=g) for g, c in COND_COLORS.items()]
fig.legend(handles=cond_handles, loc="lower center", ncol=4,
           bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=9,
           title="Condition", title_fontsize=9)
fig.suptitle("Library & Sequencing Quality Control",
             fontsize=12, fontweight="bold", y=1.005)
save_fig(fig, "Fig5_qc_metrics")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — sgRNA Count Profiles (key genes)
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 6 — sgRNA profiles …")

key_genes = ["CYP2E1", "ASL", "PLA2G4E", "CYB5A", "MTF1"]
key_genes = [g for g in key_genes if g in counts["Gene"].values]

if not key_genes:
    print("  [SKIP] Key genes not found in count table.")
else:
    sample_totals = count_matrix.sum(axis=0)
    rpm_matrix = count_matrix.div(sample_totals, axis=1) * 1e6

    sample_order = list(count_matrix.columns)
    grp_order_all = []
    for s in sample_order:
        for g in ["R0", "Rpost", "S0", "Spost"]:
            if s.startswith(g):
                grp_order_all.append(g)
                break

    fig, axes = plt.subplots(len(key_genes), 1,
                             figsize=(11, 3.4 * len(key_genes)))
    if len(key_genes) == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=0.65)

    sg_colors = plt.cm.tab10(np.linspace(0, 0.9, 10))

    for ax, gene in zip(axes, key_genes):
        rows = counts[counts["Gene"] == gene]
        if rows.empty:
            ax.set_visible(False)
            continue

        rpm = rpm_matrix.loc[rows.index]   # shape: (n_sgrna, n_samples)
        n_sg = len(rpm)

        for si, (sg_idx, sg_row) in enumerate(rpm.iterrows()):
            vals = sg_row[sample_order].values.astype(float)
            ax.plot(range(len(sample_order)), vals,
                    marker="o", ms=5, lw=1.2,
                    color=sg_colors[si % 10],
                    label=f"sg{si+1}", zorder=3)

        # Condition background shading
        prev_g, prev_x = None, 0
        for xi, grp in enumerate(grp_order_all):
            if grp != prev_g:
                if prev_g is not None:
                    ax.axvspan(prev_x - 0.5, xi - 0.5,
                               facecolor=COND_COLORS.get(prev_g, "#EEE"),
                               alpha=0.10, zorder=1)
                    mid = (prev_x + xi - 1) / 2
                    ylim_top = ax.get_ylim()[1] if ax.get_ylim()[1] != 1 else rpm.values.max() * 1.05
                    ax.text(mid, ylim_top * 0.97, prev_g,
                            ha="center", va="top", fontsize=8.5, fontweight="bold",
                            color=COND_COLORS.get(prev_g, "#444"))
                prev_g, prev_x = grp, xi
        # Draw last group
        if prev_g:
            ax.axvspan(prev_x - 0.5, len(sample_order) - 0.5,
                       facecolor=COND_COLORS.get(prev_g, "#EEE"),
                       alpha=0.10, zorder=1)
            mid = (prev_x + len(sample_order) - 1) / 2
            ylim_top = ax.get_ylim()[1] if ax.get_ylim()[1] != 1 else rpm.values.max() * 1.05
            ax.text(mid, ylim_top * 0.97, prev_g,
                    ha="center", va="top", fontsize=8.5, fontweight="bold",
                    color=COND_COLORS.get(prev_g, "#444"))

        ax.set_xticks(range(len(sample_order)))
        ax.set_xticklabels(sample_order, rotation=40, ha="right", fontsize=8)
        ax.set_ylabel("Reads per million", fontsize=9)
        ax.set_title(f"{gene}  ({n_sg} sgRNA{'s' if n_sg>1 else ''})",
                     fontsize=10, fontweight="bold", pad=8)
        ax.set_xlim(-0.6, len(sample_order) - 0.4)
        ax.legend(title="sgRNA", fontsize=7.5, title_fontsize=8,
                  loc="upper right", ncol=2,
                  frameon=True, framealpha=0.7)
        ax.grid(axis="y", lw=0.4, color="#F0F0F0")
        sns.despine(ax=ax)

    fig.suptitle("sgRNA Count Profiles — Key Hit Genes",
                 fontsize=12, fontweight="bold", y=1.005)
    save_fig(fig, "Fig6_sgrna_profiles")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 7 — Replicate Correlations
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 7 — Replicate correlations …")

sample_totals2 = count_matrix.sum(axis=0)
log_rpm = np.log1p(count_matrix.div(sample_totals2, axis=1) * 1e6)

pairs = []
for grp, samples in sample_groups.items():
    avail = [s for s in samples if s in log_rpm.columns]
    for i in range(len(avail)):
        for j in range(i + 1, len(avail)):
            pairs.append((grp, avail[i], avail[j]))

if not pairs:
    print("  [SKIP] No replicate pairs found.")
else:
    ncols = min(3, len(pairs))
    nrows = int(np.ceil(len(pairs) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.8 * ncols, 4.5 * nrows))
    axes_flat = np.array(axes).flatten()

    for ax, (grp, s1, s2) in zip(axes_flat, pairs):
        x = log_rpm[s1].values
        y = log_rpm[s2].values
        mask = (x > 0) | (y > 0)
        rho, _ = spearmanr(x[mask], y[mask])

        ax.scatter(x[mask], y[mask], s=5,
                   color=COND_COLORS.get(grp, "#888888"),
                   alpha=0.4, linewidths=0, rasterized=True)
        lim = max(x[mask].max(), y[mask].max()) * 1.05
        ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.6, zorder=2)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel(f"{s1}", fontsize=8.5)
        ax.set_ylabel(f"{s2}", fontsize=8.5)
        ax.set_title(f"{grp}", fontsize=9, fontweight="bold")
        ax.text(0.05, 0.93, f"ρ = {rho:.3f}", transform=ax.transAxes,
                fontsize=9, color=COND_COLORS.get(grp, "#333"),
                fontweight="bold")
        sns.despine(ax=ax)

    for ax in axes_flat[len(pairs):]:
        ax.set_visible(False)

    fig.suptitle("Replicate Correlation — log₁₊RPM  (Spearman ρ)",
                 fontsize=12, fontweight="bold", y=1.005)
    fig.tight_layout()
    save_fig(fig, "Fig7_replicate_corr")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 8 — Cross-Comparison Bubble Plot
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 8 — Cross-comparison bubble …")

all_genes = sorted(set().union(*[set(df["id"]) for df in data.values()]))
lfc_mat  = pd.DataFrame({c: data[c].set_index("id")["LFC"].reindex(all_genes)
                         for c in data})
pval_mat = pd.DataFrame({c: data[c].set_index("id")["pval"].reindex(all_genes)
                         for c in data})

n_sig_per_gene = (pval_mat < PVAL_THRESH).sum(axis=1)
hit_genes = n_sig_per_gene[n_sig_per_gene >= 1].sort_values(ascending=False).index[:35].tolist()

if not hit_genes:
    print("  [SKIP] No genes to plot.")
else:
    lfc_sub  = lfc_mat.loc[hit_genes]
    pval_sub = pval_mat.loc[hit_genes]

    comp_names  = list(data.keys())
    comp_labels = [COMPARISONS[c] for c in comp_names]

    fig_h = max(5, len(hit_genes) * 0.45 + 2)
    fig, ax = plt.subplots(figsize=(7.5, fig_h))

    for xi, cname in enumerate(comp_names):
        for yi, gene in enumerate(hit_genes):
            lfc  = lfc_sub.loc[gene, cname]
            pv   = pval_sub.loc[gene, cname]
            if np.isnan(lfc) or np.isnan(pv):
                continue
            size   = max(20, -np.log10(pv + 1e-10) * 65)
            color  = C_DEPL if lfc < 0 else C_ENRI
            alpha  = 0.88 if pv < PVAL_THRESH else 0.18
            zorder = 4 if pv < PVAL_THRESH else 2
            ax.scatter(xi, yi, s=size, c=color, alpha=alpha,
                       linewidths=0.4, edgecolors="white", zorder=zorder)
            if pv < PVAL_THRESH:
                txt_col = "white" if abs(lfc) > 6 else "#111111"
                ax.text(xi, yi, f"{lfc:+.1f}", ha="center", va="center",
                        fontsize=6.5, color=txt_col, fontweight="bold", zorder=5)

    ax.set_xticks(range(len(comp_names)))
    ax.set_xticklabels(comp_labels, rotation=22, ha="right", fontsize=9)
    ax.set_yticks(range(len(hit_genes)))
    ax.set_yticklabels(hit_genes, fontsize=8.5)
    ax.set_xlim(-0.6, len(comp_names) - 0.4)
    ax.set_ylim(-0.6, len(hit_genes) - 0.4)
    ax.grid(True, lw=0.4, color="#EEEEEE", zorder=0)
    ax.set_title(
        "Significant Genes Across Comparisons\n"
        "(bubble size ∝ −log₁₀p · LFC shown inside · faded = n.s.)",
        fontsize=10, fontweight="bold", pad=12
    )

    for pv_eg, label in [(0.05, "p=0.05"), (0.01, "p=0.01"), (0.001, "p=0.001")]:
        ax.scatter([], [], s=-np.log10(pv_eg) * 65, c="#888888",
                   alpha=0.7, label=label)
    ax.legend(title="Bubble size", title_fontsize=8, fontsize=8,
              loc="lower right", frameon=True, framealpha=0.9)
    sns.despine(ax=ax)
    fig.tight_layout()
    save_fig(fig, "Fig8_cross_comparison")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 9 — Summary Table
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 9 — Summary table …")

FUNC_ANNOT = {
    "CYP2E1":   "CYP450 2E1 — drug metabolism, ROS",
    "NAT2":     "N-acetyltransferase 2 — drug acetylation",
    "ASL":      "Argininosuccinate lyase — urea cycle",
    "GBA3":     "Beta-glucosidase — glycolipid catabolism",
    "TRPM7":    "TRPM7 channel — Mg²⁺/Ca²⁺ homeostasis",
    "ATP6V0A4": "V-ATPase a4 — lysosomal acidification",
    "TPCN1":    "Two-pore channel 1 — lysosomal Ca²⁺",
    "HSD11B1":  "11β-HSD1 — cortisol/corticosterone",
    "CYB5A":    "Cytochrome b5A — electron transfer",
    "ALOXE3":   "eLOX-3 — epidermis lipid signalling",
    "MTF1":     "Metal-reg. TF — oxidative stress response",
    "DHRSX":    "DHRSX — dehydrogenase/reductase",
    "PLA2G4E":  "Group IVE PLA₂ — arachidonate release",
    "SLC1A1":   "EAAT3 glutamate transporter",
    "ARSD":     "Arylsulfatase D — lysosomal sulfatase",
    "GRIN1":    "NMDA receptor subunit 1",
    "RRM2":     "Ribonucleotide reductase M2 — dNTPs",
    "PFKFB1":   "PFKFB1 — glycolysis (liver isoform)",
    "SLC25A15": "Ornithine carrier — urea cycle",
    "SLC6A3":   "DAT — dopamine transporter",
    "SLC44A1":  "CTL1 — choline transporter",
}

if not sig_df.empty:
    tbl = sig_df.copy()
    tbl["pval_str"] = tbl["pval"].apply(lambda v: f"{v:.4f}")
    tbl["lfc_str"]  = tbl["LFC"].apply(lambda v: f"{v:+.2f}")
    tbl["function"] = tbl["gene"].map(lambda g: FUNC_ANNOT.get(g, "—"))
    tbl = tbl.sort_values(["comparison", "pval"]).reset_index(drop=True)

    col_labels = ["Comparison", "Gene", "p-value", "LFC", "Direction", "Function"]
    col_keys   = ["comparison_label", "gene", "pval_str", "lfc_str", "direction", "function"]
    col_widths = [0.20, 0.08, 0.07, 0.06, 0.09, 0.40]

    row_h_in  = 0.30
    fig_h     = max(3.5, len(tbl) * row_h_in + 1.5)
    fig, ax   = plt.subplots(figsize=(14.5, fig_h))
    ax.axis("off")

    cell_text = [
        [str(row[k])[:42] for k in col_keys]
        for _, row in tbl.iterrows()
    ]

    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        colWidths=col_widths,
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.40)

    # Header row
    for j in range(len(col_labels)):
        cell = table[0, j]
        cell.set_facecolor("#1A3A4A")
        cell.set_text_props(color="white", fontweight="bold")

    # Data rows
    for i, (_, row) in enumerate(tbl.iterrows(), start=1):
        bg = "#EBF4FB" if i % 2 == 0 else "white"
        for j in range(len(col_labels)):
            cell = table[i, j]
            cell.set_facecolor(bg)
            cell.set_edgecolor("#DDDDDD")
            if col_keys[j] == "direction":
                c = C_DEPL if row["direction"] == "Depleted" else C_ENRI
                cell.set_text_props(color=c, fontweight="bold")
            if col_keys[j] == "gene":
                cell.set_text_props(fontweight="bold")

    ax.set_title(
        f"All Significant Genes — p < {PVAL_THRESH}  |  GBM TRAIL Metabolic CRISPR Screen\n"
        "Note: Standard FDR > 0.25 due to severe library bottlenecking; "
        "p-value used as primary criterion.",
        fontsize=9.5, fontweight="bold", pad=12, loc="center"
    )
    fig.tight_layout()
    save_fig(fig, "Fig9_summary_table")


# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*65)
print("Done. All outputs written to:")
print(f"  Figures:   {FIG_DIR}")
print(f"  Count:     {COUNT_DIR}")
print(f"  Test:      {TEST_DIR}")
print("="*65)
print("\nFiles:")
for f in sorted(FIG_DIR.iterdir()):
    print(f"  {f.name}")
