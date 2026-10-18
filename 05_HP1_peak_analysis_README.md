Markdown · 05_HP1_peak_analysis_README.md

# 05 — HP1 ChIP-seq Peak Analysis

This section covers MACS3 peak calling on HP1α and HP1β ChIP-seq data across four conditions, consensus peak generation, promoter and H3K9me3 overlap analysis, and bigWig generation for genome browser visualization.

---

## Samples

Twenty libraries were sequenced across four conditions.
| Sample ID | GEO library name | Description |
| --- | --- | --- |
| TS1 | MCF10A_WT_ChIP-HP1a_rep1 | HP1α ChIP, WT, rep1 |
| TS2 | MCF10A_WT_ChIP-HP1a_rep2 | HP1α ChIP, WT, rep2 |
| TS7 | MCF10A_WT_ChIP-HP1b_rep1 | HP1β ChIP, WT, rep1 |
| TS8 | MCF10A_WT_ChIP-HP1b_rep2 | HP1β ChIP, WT, rep2 |
| TS3 | MCF10A_HP1BKD_ChIP-HP1a_rep1 | HP1α ChIP, shHP1β KD, rep1 |
| TS4 | MCF10A_HP1BKD_ChIP-HP1a_rep2 | HP1α ChIP, shHP1β KD, rep2 |
| TS9 | MCF10A_HP1AKD_ChIP-HP1b_rep1 | HP1β ChIP, shHP1α KD, rep1 |
| TS10 | MCF10A_HP1AKD_ChIP-HP1b_rep2 | HP1β ChIP, shHP1α KD, rep2 |
| TS5 | MCF10A_H2AZKD_ChIP-HP1a_rep1 | HP1α ChIP, shH2AZ KD, rep1 |
| TS6 | MCF10A_H2AZKD_ChIP-HP1a_rep2 | HP1α ChIP, shH2AZ KD, rep2 |
| TS11 | MCF10A_H2AZKD_ChIP-HP1b_rep1 | HP1β ChIP, shH2AZ KD, rep1 |
| TS12 | MCF10A_H2AZKD_ChIP-HP1b_rep2 | HP1β ChIP, shH2AZ KD, rep2 |
| TS13 | MCF10A_HP1AKD_ChIP-H2AZ_rep1 | H2AZ ChIP, shHP1α KD, rep1 |
| TS14 | MCF10A_HP1AKD_ChIP-H2AZ_rep2 | H2AZ ChIP, shHP1α KD, rep2 |
| TS15 | MCF10A_HP1AKD_Input_rep1 | Input, shHP1α KD, rep1 |
| TS16 | MCF10A_HP1AKD_Input_rep2 | Input, shHP1α KD, rep2 |
| TS17 | MCF10A_HP1BKD_Input_rep1 | Input, shHP1β KD, rep1 |
| TS18 | MCF10A_HP1BKD_Input_rep2 | Input, shHP1β KD, rep2 |
| TS19 | MCF10A_H2AZKD_Input_rep1 | Input, shH2AZ KD, rep1 |
| TS20 | MCF10A_H2AZKD_Input_rep2 | Input, shH2AZ KD, rep2 |


---

## Upstream processing

Trimming, alignment, BAM sorting, and duplicate removal followed the same shared pipeline described in Section 02 (ChIP-seq branch). All samples were aligned to hs1/CHM13v2.0.

---

## Step 1 — MACS3 peak calling

Peaks were called in broad mode for all HP1α, HP1β, and H2AZ ChIP samples using their condition-matched input controls. Each replicate was called independently.

```bash
macs3 callpeak \
    -t sample_rmdup.bam \
    -c input_rmdup.bam \
    -f BAMPE \
    -g 2948000000 \
    --broad \
    --outdir macs3_output/ \
    -n sample_name \
    --keep-dup all
```

**Input control assignments:**

| ChIP condition | Input used |
| --- | --- |
| WT HP1α (TS1, TS2) | TS19, TS20 |
| WT HP1β (TS7, TS8) | TS19, TS20 |
| shHP1β KD HP1α (TS3, TS4) | TS17, TS18 |
| shHP1α KD HP1β (TS9, TS10) | TS15, TS16 |
| shH2AZ KD HP1α (TS5, TS6) | TS19, TS20 |
| shH2AZ KD HP1β (TS11, TS12) | TS19, TS20 |
| shHP1α KD H2AZ (TS13, TS14) | TS15, TS16 |

Output: `*.broadPeak` files used for consensus peak generation.

---

## Step 2 — Consensus peak generation

Consensus peaks were generated per condition by intersecting the two replicates and then merging overlapping intervals. This was done for all six ChIP targets (HP1α WT, HP1β WT, HP1α shHP1β KD, HP1β shHP1α KD, HP1α shH2AZ KD, HP1β shH2AZ KD).

```bash
# Intersect two replicates — keep only peaks present in both
bedtools intersect \
    -a sample_rep1.broadPeak \
    -b sample_rep2.broadPeak \
    > sample_intersect.bed

# Merge overlapping intervals into consensus peaks
bedtools merge \
    -i sample_intersect.bed \
    > sample_consensus_peaks.bed
```

This produces one consensus peak BED file per condition and target.

---

## Step 3 — Promoter overlap analysis

Consensus peaks were intersected with a TSS ±1 kb window file to identify promoter-associated HP1 peaks. The `-u` flag reports each peak once regardless of how many TSS intervals it overlaps.

```bash
# Peaks overlapping TSS ±1 kb
bedtools intersect \
    -a sample_consensus_peaks.bed \
    -b /home/jbenoit/index/T2T_TSS_1kb_reduced_clean.bed \
    -u \
    > sample_promoter_peaks.bed

# Peaks not overlapping TSS ±1 kb
bedtools intersect \
    -a sample_consensus_peaks.bed \
    -b /home/jbenoit/index/T2T_TSS_1kb_reduced_clean.bed \
    -v \
    > sample_nonpromoter_peaks.bed
```

Peak counts from promoter and non-promoter files were used to generate the promoter HP1 peak count figure (Section 08).

---

## Step 4 — H3K9me3 overlap analysis

Promoter-associated HP1 consensus peaks were intersected with H3K9me3 CUT&RUN consensus peaks (generated in Section 02) to classify peaks as H3K9me3-positive or H3K9me3-negative.

```bash
# HP1 promoter peaks overlapping H3K9me3
bedtools intersect \
    -a sample_promoter_peaks.bed \
    -b H3K9me3_consensus_peaks.bed \
    -u \
    > sample_promoter_H3K9me3_positive.bed

# HP1 promoter peaks not overlapping H3K9me3
bedtools intersect \
    -a sample_promoter_peaks.bed \
    -b H3K9me3_consensus_peaks.bed \
    -v \
    > sample_promoter_H3K9me3_negative.bed
```

The proportion of HP1 promoter peaks overlapping H3K9me3 was compared between WT and KD_H2AZ conditions for both isoforms (Section 08 figures).

---

## Step 5 — bigWig generation

Log2 ratio bigWigs (ChIP/input) were generated using deepTools bamCompare for genome browser visualization and heatmap generation.

```bash
bamCompare \
    -b1 sample_rmdup.bam \
    -b2 input_rmdup.bam \
    --binSize 10 \
    --normalizeUsing RPKM \
    --effectiveGenomeSize 2948000000 \
    --operation log2 \
    --skipZeroOverZero \
    --extendReads \
    -o sample_log2ratio.bw \
    -p 20
```

---

## Output files

| File | Description |
| --- | --- |
| `*_consensus_peaks.bed` | Consensus peaks per condition and target |
| `*_promoter_peaks.bed` | Consensus peaks overlapping TSS ±1 kb |
| `*_nonpromoter_peaks.bed` | Consensus peaks outside TSS ±1 kb |
| `*_promoter_H3K9me3_positive.bed` | Promoter peaks overlapping H3K9me3 |
| `*_promoter_H3K9me3_negative.bed` | Promoter peaks not overlapping H3K9me3 |
| `*_log2ratio.bw` | log2(ChIP/input) bigWig tracks |

---

## GEO deposit

Raw FASTQ files for all 20 libraries were deposited to NCBI GEO. No condition-specific input was sequenced for the WT condition; TS19/TS20 serve as the input proxy and are listed under the shH2AZ KD condition.
