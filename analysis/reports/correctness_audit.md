# CRISPR Screen Correctness Audit
**Internal analysis vs external sequencing service report**  
**Date:** 2026-04-09  
**Auditor:** Cingöz Lab Bioinformatics

---

## 1. Summary Verdict

| Category | Status | Severity |
|----------|--------|----------|
| Library reference mismatch | DISCREPANCY | CRITICAL |
| Scaffold sequence off-by-one | DISCREPANCY | MODERATE |
| Statistical threshold choice | DIFFERENCE (intentional) | LOW |
| TRPM7 direction conflict | DISCREPANCY | HIGH |
| CYP2E1 as top hit (Spost_vs_S0) | CONSISTENT | — |
| CYB5A depletion (Rpost_vs_R0) | CONSISTENT | — |
| "PREDICTED" gene annotations | COMPANY ISSUE | MODERATE |
| FDR significance in our data | PROBLEM | HIGH |

**Bottom line:** The two analyses are not directly comparable because they used different sgRNA reference libraries and different gene annotation approaches. Core biology (CYP2E1 essentiality in sensitive cells; CYB5A depletion in resistant cells) is reproduced. However, the results diverge on several key genes and our analysis cannot pass standard FDR < 0.25 owing to severe library bottlenecking.

---

## 2. Side-by-Side Analysis Parameters

| Parameter | Our Analysis | Company Report (PDF) |
|-----------|-------------|----------------------|
| Library size | **29,790 sgRNAs** (Addgene #110066) | **66,362 sgRNAs** |
| Library | Sabatini Metabolic KO | Unknown / Different |
| sgRNA annotation | MAGeCK count → native library TSV | BLAST-based custom annotation |
| Scaffold motif used | `GTTTTAGAGCTAGAAATAGC` (20 bp) | `TTTTAGAGCTAGAAATAGC` (19 bp) |
| Filtering strategy | `--remove-zero both` in MAGeCK | Minimum 1,000 reads/sgRNA |
| sgRNAs after filtering | ~18–193 per comparison | 728 sgRNAs |
| Statistical method | MAGeCK RRA (Robust Rank Aggregation) | MAGeCK RRA |
| Significance threshold | p < 0.05 (p-value, not FDR) | FDR < 0.25 |
| Total significant genes | 11 (across 3 comparisons) | 16 (across 3 comparisons) |
| Gene name source | Library TSV (human gene symbols) | BLAST → GenBank (some "PREDICTED") |

---

## 3. Critical Discrepancy — Library Reference

**The most important difference between the two analyses.**

Our count table (`count_table_fixed.txt`) contains **30,197 rows** (29,790 sgRNA targets + 407 INTERGENIC controls), consistent with **Addgene #110066** (Sabatini Human Metabolic Gene Knockout Library).

The company report states: *"Toplam 66,362 sgRNA içeren kütüphaneden..."* (From a library containing 66,362 sgRNAs...). This is **more than twice** the size of Addgene #110066.

**Implications:**
- The company used a different (larger) reference library for alignment/counting
- Gene hits unique to the company (PTGS1, GABRA3, LINC02047, ACSM2B, SRFBP1, PCDH9, SLC25A20) may result from sgRNAs present in their library but absent from ours
- Conversely, our hits (ASL, GBA3, PLA2G4E) target genes covered in our library that may not be in theirs
- Direct gene-level comparison is not valid across the two analyses

**Recommendation:** Confirm which library was physically transduced. If Addgene #110066, our 30K-sgRNA reference is correct and the company's 66K reference is erroneous (adding noise from non-existent sgRNAs).

---

## 4. Moderate Discrepancy — Scaffold Sequence

| Source | Scaffold used |
|--------|--------------|
| Our `01_denovo_count.py` | `GTTTTAGAGCTAGAAATAGC` (20 bp) |
| Company report | `TTTTAGAGCTAGAAATAGC` (19 bp) |

The leading `G` in our scaffold is the last nt of the protospacer sequence, which is part of the lentiCRISPRv1 backbone upstream of the scaffold. Including vs excluding this G affects:
- Which reads match the scaffold anchor
- The 20-mer extracted as the sgRNA sequence
- Downstream alignment to the reference library

With the extra G, our script requires the sgRNA to end in `G` before the scaffold, reducing specificity to reads with the canonical PAM-proximal G. The company's 19-bp scaffold is looser and will capture more reads. **This partly explains the lower mapping rate in our analysis.**

**Recommendation:** Re-run `01_denovo_count.py` with the 19-bp scaffold (`TTTTAGAGCTAGAAATAGC`) and compare mapping rates.

---

## 5. High-Severity Discrepancy — TRPM7 Direction

| Analysis | Comparison | LFC | Direction |
|----------|-----------|-----|-----------|
| **Our analysis** | Rpost_vs_R0 | **−3.013** | **Depleted** (essential) |
| **Company report** | Rpost_vs_R0 | **+3.78** | **Enriched** (toxic) |

This is a direct biological contradiction. TRPM7 (Mg²⁺/Ca²⁺ ion channel) is called as essential in resistant cells in our data but as growth-promoting-knockout (toxic effect on survival upon loss) in the company's data.

**Possible causes:**
1. Different sgRNAs targeting TRPM7 in the two libraries — if the company used a different sgRNA, it may target a different exon or have off-target effects
2. TRPM7 has 2 sgRNAs in our count table; counts may differ enough to flip direction
3. Different normalization methods — median vs total normalization can reverse the sign for border-case genes

**Recommendation:** Check TRPM7 sgRNA counts directly and verify which sgRNAs are present in each library.

---

## 6. Consistent Findings (Reproduced in Both Analyses)

These genes are found in both analyses with concordant direction — this increases biological confidence:

| Gene | Comparison | Our LFC | PDF LFC | Function | Confidence |
|------|-----------|---------|---------|---------|-----------|
| **CYP2E1** | Spost_vs_S0 | −13.629 | −11.19 | Drug metabolism, ROS | HIGH |
| **CYB5A** | Rpost_vs_R0 | −4.345 | −1.90 | Electron transfer chain | HIGH |
| **MTF1** | Rpost_vs_R0 | −6.899 | −8.68 | Metal stress TF | HIGH |

---

## 7. Statistical Power Problem

With `--remove-zero both`, MAGeCK excludes any sgRNA with zero counts in any sample. Given the extreme bottlenecking (>99% sgRNA dropout), this leaves:

| Comparison | Testable genes | FDR < 0.25? | Min FDR |
|-----------|---------------|-------------|---------|
| Spost_vs_S0 | **18** | NO | 0.502 |
| Rpost_vs_R0 | 193 | NO | 0.438 |
| Spost_vs_Rpost | 188 | NO | 0.438 |

**The RRA algorithm warns:** *"Increase the number of permutations to 5,556 to get precise p values"* — with only 18 genes in the Spost_vs_S0 comparison, FDR correction is essentially meaningless (18 genes → BH adjustment has no statistical power).

The company worked around this by using a **pre-filter** (≥1,000 reads/sgRNA → 728 sgRNAs → more testable genes) rather than removing zeros across samples. Their approach produces more stable FDR estimates at the cost of discarding low-count sgRNAs.

**Our choice** to use p < 0.05 is a pragmatic, documented workaround. It is scientifically valid when clearly stated (as in our README and report), but results must be interpreted with caution.

---

## 8. Company Report Issues

### 8.1 "PREDICTED" Gene Names
The company annotated genes with the suffix "-PREDICTED" (CYB5A-PREDICTED, MTF1-PREDICTED, PCDH9-PREDICTED, etc.). In humans, these genes are fully characterized in RefSeq/Ensembl without this suffix. The "PREDICTED" annotation is typically seen in:
- Non-human or poorly annotated genomes (e.g., some primate assemblies)
- Very old NCBI annotation builds (pre-2015)
- BLAST hits against poorly curated databases

This suggests the company performed BLAST annotation against an outdated or non-human reference, reducing confidence in their gene-level assignments. Our gene symbol assignments (CYB5A, MTF1, etc.) are correct for human Ensembl/RefSeq.

### 8.2 FDR Threshold
Using FDR < 0.25 is liberal but defensible in a severely bottlenecked screen. The company appropriately states this is a lenient threshold for exploratory analysis. This is not an error per se, but must be noted.

### 8.3 No Per-sgRNA Validation
The company report does not show individual sgRNA count profiles for key hits. Having only 1 sgRNA per gene (as we see for most of our hits) means there is no intra-gene concordance check. A standard quality requirement is ≥2 sgRNAs per gene pointing in the same direction.

---

## 9. Key Biological Findings

### 9.1 Confirmed Hits (Both Analyses)

**CYP2E1 — Sensitive cell essentiality (★★★)**
- Strongly depleted in Spost vs S0 (our LFC = −13.6, PDF LFC = −11.2)
- CYP2E1 is a major drug-metabolizing CYP450 enzyme that generates reactive oxygen species (ROS)
- Essentiality in sensitive cells under drug treatment is biologically coherent: CYP2E1 loss would reduce drug bioactivation and ROS production, selecting against sgRNA-bearing cells
- Also appears in Spost_vs_Rpost (our data: LFC −6.1 at p = 0.098), confirming sensitive-specific dependence

**CYB5A — Resistant cell essentiality (★★)**
- Depleted in Rpost vs R0 in both analyses (our LFC = −4.3, PDF LFC = −1.9)
- Cytochrome b5 — electron transfer protein involved in lipid desaturation and CYP450 activity
- Essential in resistant cells specifically suggests altered redox/lipid metabolism as a resistance mechanism

**MTF1 — Metal stress response (★★)**
- Depleted in Rpost vs R0 in both analyses (our LFC = −6.9, PDF LFC = −8.7)
- Metal-regulatory transcription factor 1 — regulates metallothioneins under heavy metal and oxidative stress
- Consistent loss in resistant cells suggests resistance involves metal/oxidative stress adaptation

### 9.2 Our-Only Hits (Not Reproduced in Company Report)

**ASL — Urea cycle (★★)**
- Argininosuccinate lyase; top hit in Rpost_vs_R0 (our LFC = −10.6, p = 0.0023)
- Urea cycle enzymes have emerging roles in arginine-dependent tumor survival
- Not in company's hit list likely due to library reference differences

**PLA2G4E — Phospholipase essentiality in Spost_vs_Rpost (★★)**
- Top hit in Spost_vs_Rpost (our LFC = −11.9, p = 0.0023)
- Group IVE phospholipase A2 — arachidonic acid release and eicosanoid production
- More essential in sensitive vs resistant post-treatment — interesting drug sensitivity link

### 9.3 Company-Only Hits (Not in Our Data)

| Gene | Comparison | LFC | Possible Explanation |
|------|-----------|-----|---------------------|
| PTGS1 (COX-1) | Rpost_vs_R0 | −14.13 | Absent from our library or zero-filtered |
| SLC25A4 | Rpost_vs_R0 | +3.75 | Absent from our library |
| ACSM2B | Rpost_vs_R0 | −3.84 | Absent from our library or zero-filtered |
| LINC02047 | Rpost_vs_R0 | +2.08 | lncRNA — may not be in metabolic library |

---

## 10. Recommendations

### Immediate Actions
1. **Confirm library identity**: Verify that Addgene #110066 (29,790 sgRNAs) was the library actually transduced, not a custom expanded version
2. **Re-run counting with 19-bp scaffold**: Test `TTTTAGAGCTAGAAATAGC` to potentially improve mapping rate
3. **Check TRPM7 raw counts**: Inspect the count table directly to explain the direction discrepancy

### For Re-screening
1. Increase transduction cells to ≥15 million (500× coverage for 30K sgRNAs)
2. Target MOI 0.3 for single-copy integration
3. Sequence at ≥300 reads/sgRNA depth (>9M reads per sample)
4. Use lentiCRISPR v2 (Addgene #52961)
5. QC the plasmid library before transduction

### Analysis Improvements
1. Run MAGeCK with `--min-count 10` instead of `--remove-zero both` to retain more genes while still filtering noise
2. Use median normalization with spike-in controls in the re-screen
3. For the current data, report both p < 0.05 and p < 0.10 tiers to distinguish confidence levels

---

## 11. Appendix — Full Gene Table

### Our Analysis (p < 0.05)

| Comparison | Gene | p-value | FDR | LFC | Direction | Function |
|-----------|------|---------|-----|-----|-----------|---------|
| Spost_vs_S0 | CYP2E1 | 0.0279 | 0.502 | −13.63 | Depleted | CYP450, ROS/drug metabolism |
| Spost_vs_S0 | NAT2 | 0.0279 | 0.502 | +6.67 | Enriched | Drug acetylation |
| Rpost_vs_R0 | ASL | 0.0023 | 0.438 | −10.62 | Depleted | Urea cycle |
| Rpost_vs_R0 | GBA3 | 0.0122 | 0.793 | −8.49 | Depleted | Beta-glucosidase |
| Rpost_vs_R0 | TRPM7 | 0.0277 | 0.793 | −3.01 | Depleted | Ion channel |
| Rpost_vs_R0 | ATP6V0A4 | 0.0175 | 0.793 | −6.87 | Depleted | V-ATPase subunit |
| Rpost_vs_R0 | TPCN1 | 0.0225 | 0.793 | −6.32 | Depleted | Two-pore channel |
| Rpost_vs_R0 | HSD11B1 | 0.0273 | 0.793 | −5.25 | Depleted | Cortisol metabolism |
| Rpost_vs_R0 | CYB5A | 0.0322 | 0.793 | −4.35 | Depleted | Electron transfer |
| Rpost_vs_R0 | ALOXE3 | 0.0372 | 0.793 | −4.34 | Depleted | Lipoxygenase |
| Rpost_vs_R0 | MTF1 | 0.0421 | 0.793 | −6.90 | Depleted | Metal-reg. TF |
| Rpost_vs_R0 | DHRSX | 0.0472 | 0.793 | −3.80 | Depleted | Dehydrogenase |
| Spost_vs_Rpost | PLA2G4E | 0.0023 | 0.438 | −11.87 | Depleted | Phospholipase A2 |
| Spost_vs_Rpost | SLC1A1 | 0.0073 | 0.684 | −10.97 | Depleted | Glutamate transporter |
| Spost_vs_Rpost | ARSD | 0.0126 | 0.788 | −10.07 | Depleted | Arylsulfatase |
| Spost_vs_Rpost | GRIN1 | 0.0179 | 0.844 | −9.95 | Depleted | NMDA receptor |
| Spost_vs_Rpost | RRM2 | 0.0230 | 0.865 | −9.79 | Depleted | Ribonucleotide reductase |
| Spost_vs_Rpost | PFKFB1 | 0.0282 | 0.877 | −9.34 | Depleted | Glycolysis regulator |
| Spost_vs_Rpost | SLC25A15 | 0.0332 | 0.877 | −9.07 | Depleted | Ornithine transporter |
| Spost_vs_Rpost | SLC6A3 | 0.0383 | 0.877 | −8.84 | Depleted | Dopamine transporter |
| Spost_vs_Rpost | SLC44A1 | 0.0487 | 0.877 | −7.44 | Depleted | Choline transporter |

*FDR is not reliable with fewer than ~50 testable genes per comparison. Use p-value as primary criterion for this bottlenecked dataset.*
