#!/usr/bin/env python3
"""
01_denovo_count.py
==================
De novo sgRNA counting from CRISPR screen FASTQ files.

This script does NOT require a pre-existing sgRNA library file.
It extracts all 20-nt sequences upstream of the scaffold motif
and builds a count matrix across all samples.

Why de novo counting?
---------------------
The physical CRISPR library used in the experiment may differ from
the designed reference library. A de novo approach extracts and
counts whatever sgRNA sequences are actually present in the data,
ensuring accurate quantification regardless of library version.

Workflow:
  1. Extract 20-nt sequences upstream of the CRISPR scaffold
  2. Pool across all samples to build a reference set
  3. Filter for sgRNAs present in ≥ MIN_SAMPLES samples
  4. Build and save a per-sample count matrix (TSV format)
  5. Output the discovered sgRNA reference library for annotation

Output:
  results/denovo_count_table.txt   — count matrix (MAGeCK-compatible)
  results/denovo_library.tsv       — discovered sgRNAs

Usage:
  python scripts/01_denovo_count.py
  python scripts/01_denovo_count.py --fastq-dir /path/to/fastq --outdir results/
"""

import argparse
import gzip
import os
import sys
import collections
from pathlib import Path

# ── Default configuration ────────────────────────────────────────────────────
SCAFFOLD     = "GTTTTAGAGCTAGAAATAGC"   # gRNA scaffold motif
SGRNA_LEN   = 20                         # expected guide length
MIN_READS    = 5                          # min total reads across all samples
MIN_SAMPLES  = 1                          # min samples with ≥1 read

FASTQ_DIR    = Path(__file__).parents[1] / "Fastq"
OUT_DIR      = Path(__file__).parents[1] / "results"

# Sample name → FASTQ R1 file mapping
# Edit these filenames to match your sequencing data
SAMPLES = {
    "R0_1":    "R0_1_R1.fq.gz",
    "R0_2":    "R0_2_R1.fq.gz",
    "R0_3":    "R0_3_R1.fq.gz",
    "Rpost_1": "Rpost_1_R1.fq.gz",
    "Rpost_2": "Rpost_2_R1.fq.gz",
    "Rpost_3": "Rpost_3_R1.fq.gz",
    "Rpost_4": "Rpost_4_R1.fq.gz",
    "S0_1":    "S0_1_R1.fq.gz",
    "S0_2":    "S0_2_R1.fq.gz",
    "S0_3":    "S0_3_R1.fq.gz",
    "Spost_1": "Spost_1_R1.fq.gz",
    "Spost_2": "Spost_2_R1.fq.gz",
    "Spost_3": "Spost_3_R1.fq.gz",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fastq-dir", type=Path, default=FASTQ_DIR,
                   help=f"Directory containing FASTQ files (default: {FASTQ_DIR})")
    p.add_argument("--outdir", type=Path, default=OUT_DIR,
                   help=f"Output directory (default: {OUT_DIR})")
    p.add_argument("--scaffold", default=SCAFFOLD,
                   help=f"sgRNA scaffold sequence (default: {SCAFFOLD})")
    p.add_argument("--sgrna-len", type=int, default=SGRNA_LEN,
                   help=f"sgRNA length (default: {SGRNA_LEN})")
    p.add_argument("--min-reads", type=int, default=MIN_READS,
                   help=f"Min total reads to keep an sgRNA (default: {MIN_READS})")
    return p.parse_args()


def count_sample(fastq_path: Path, scaffold: str, sgrna_len: int) -> collections.Counter:
    """Count sgRNA occurrences in one FASTQ file."""
    counts: collections.Counter = collections.Counter()
    total = found = 0

    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    with opener(fastq_path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().strip()
            fh.readline()   # +
            fh.readline()   # quality

            total += 1
            idx = seq.find(scaffold)
            if idx >= sgrna_len:
                guide = seq[idx - sgrna_len : idx]
                counts[guide] += 1
                found += 1

    pct = 100 * found / total if total else 0
    print(f"    {fastq_path.name}: {total:>10,} reads  →  "
          f"{found:>10,} matched ({pct:.1f}%)  "
          f"unique: {len(counts):,}")
    return counts


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 65)
    print("De novo sgRNA Counting")
    print(f"Scaffold : {args.scaffold}")
    print(f"sgRNA len: {args.sgrna_len}")
    print("=" * 65)

    # ── Step 1: Count each sample ────────────────────────────────────────────
    all_counts: dict[str, collections.Counter] = {}
    for label, fname in SAMPLES.items():
        fq = args.fastq_dir / fname
        if not fq.exists():
            print(f"  [WARN] {fq} not found — skipping {label}")
            continue
        print(f"\n  Counting {label} …")
        all_counts[label] = count_sample(fq, args.scaffold, args.sgrna_len)

    if not all_counts:
        sys.exit("ERROR: No samples found. Check --fastq-dir path.")

    # ── Step 2: Build union of sgRNAs ────────────────────────────────────────
    all_guides: collections.Counter = collections.Counter()
    for cnts in all_counts.values():
        all_guides += cnts

    print(f"\nTotal unique sgRNAs detected: {len(all_guides):,}")

    # ── Step 3: Filter ───────────────────────────────────────────────────────
    kept = {g: c for g, c in all_guides.items() if c >= args.min_reads}
    print(f"After min-reads filter (≥{args.min_reads}): {len(kept):,}")

    if not kept:
        sys.exit("ERROR: No sgRNAs passed the filter. Try lowering --min-reads.")

    guides_sorted = sorted(kept, key=lambda g: -all_guides[g])

    # ── Step 4: Write de novo library TSV ────────────────────────────────────
    lib_path = args.outdir / "denovo_library.tsv"
    with open(lib_path, "w") as f:
        f.write("id\tsequence\tgene\n")
        for i, g in enumerate(guides_sorted, 1):
            f.write(f"sgrna_{i:05d}\t{g}\tUNKNOWN\n")
    print(f"\nLibrary written: {lib_path}  ({len(guides_sorted):,} sgRNAs)")
    print("  NOTE: gene names are 'UNKNOWN' — annotate with BLAST/Bowtie2 against genome")

    # ── Step 5: Write count matrix ───────────────────────────────────────────
    samples = list(all_counts.keys())
    count_path = args.outdir / "denovo_count_table.txt"
    with open(count_path, "w") as f:
        f.write("\t".join(["sgRNA", "Gene"] + samples) + "\n")
        for i, g in enumerate(guides_sorted, 1):
            row = [f"sgrna_{i:05d}", "UNKNOWN"]
            row += [str(all_counts[s].get(g, 0)) for s in samples]
            f.write("\t".join(row) + "\n")

    print(f"Count table written: {count_path}")

    # ── Step 6: Print summary ────────────────────────────────────────────────
    print("\n── Top 20 most abundant sgRNAs ──────────────────────────────")
    print(f"{'Rank':>5}  {'Sequence':<22}  {'Total reads':>12}")
    for i, g in enumerate(guides_sorted[:20], 1):
        print(f"{i:>5}  {g:<22}  {all_guides[g]:>12,}")


if __name__ == "__main__":
    main()
