# 02 — ChIP-seq and CUT&RUN Processing Pipeline

This section covers the full processing pipeline for both ChIP-seq and CUT&RUN data. The upstream steps (trimming, alignment, BAM filtering) are identical for both assay types. The two pipelines diverge at the normalization and peak calling steps, which are described separately below.

All steps were run on pauper2 HPC. MACS3 peak calling was run in the dedicated `MACS3` conda environment.

---

## Shared upstream pipeline

### Samples

**ChIP-seq** — HP1α and HP1β ChIP across three cell line conditions:

| Sample ID | Cell line | Antibody | Replicate |
|---|---|---|---|
| TS1 | MCF10A WT | HP1α | rep1 |
| TS2 | MCF10A WT | HP1α | rep2 |
| TS3 | MCF10A shHP1β KD | HP1α | rep1 |
| TS4 | MCF10A shHP1β KD | HP1α | rep2 |
| TS5 | MCF10A shH2AZ KD | HP1α | rep1 |
| TS6 | MCF10A shH2AZ KD | HP1α | rep2 |
| TS7 | MCF10A WT | HP1β | rep1 |
| TS8 | MCF10A WT | HP1β | rep2 |
| TS9 | MCF10A shHP1α KD | HP1β | rep1 |
| TS10 | MCF10A shHP1α KD | HP1β | rep2 |
| TS11 | MCF10A shH2AZ KD | HP1β | rep1 |
| TS12 | MCF10A shH2AZ KD | HP1β | rep2 |
| TS15 | MCF10A shHP1α KD | Input | rep1 |
| TS16 | MCF10A shHP1α KD | Input | rep2 |
| TS17 | MCF10A shHP1β KD | Input | rep1 |
| TS18 | MCF10A shHP1β KD | Input | rep2 |
| TS19 | MCF10A shH2AZ KD | Input | rep1 |
| TS20 | MCF10A shH2AZ KD | Input | rep2 |

> **Note on input controls:** Condition-matched inputs were used where available. No WT-specific input was sequenced; the closest available condition-matched input was used as a proxy for WT ChIP samples (shHP1α KD input for WT HP1α ChIP; shHP1β KD input for WT HP1β ChIP).

**CUT&RUN** — histone marks and architectural proteins in MCF10A WT cells (2 biological replicates per target):

| Target | Peak type |
|---|---|
| H3K4me3 | narrow |
| H3K27me3 | broad |
| H3K9me3 | broad |
| HP1α | broad |
| HP1β | broad |
| CTCF | narrow |
| IgG | control |

---

### Reference files

- Human genome: T2T CHM13v2.0 (hs1), Bowtie2 index at `chm13v2.0/chm13v2.0`
- CUT&RUN spike-in genome: *S. cerevisiae* R64-1-1, Bowtie2 index at appropriate path
- Adapter sequences: `TruSeq3-PE.fa`

---

### Step 1 — Adapter trimming

```bash
for i in *_R1*; do
    java -jar /home/bioinformatics/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -threads 20 \
        -phred33 \
        $i ${i/R1/R2} \
        ${i/R1/R1_paired} \
        ${i/R1/R1_unpaired} \
        ${i/R1/R2_paired} \
        ${i/R1/R2_unpaired} \
        ILLUMINACLIP:/home/jbenoit/index/TruSeq3-PE.fa:2:30:10:1:TRUE \
        MINLEN:20 \
        LEADING:5 \
        SLIDINGWINDOW:4:15
done &> trimming_B.log &
```

---

### Step 2 — Alignment to hs1

```bash
for i in *_R1_paired*; do
    bowtie2 \
        --end-to-end \
        --no-mixed \
        --no-discordant \
        -p 20 \
        -x /home/jbenoit/index/chm13v2.0/chm13v2.0 \
        -1 $i \
        -2 ${i/R1/R2} \
        -S ${i/.fastq.gz/.sam}
done &> T2T_alignment.log &
```

---

### Step 3 — SAM to BAM conversion and quality filtering

```bash
for i in *sam; do
    samtools view -b -q 10 $i > ${i/.sam/.bam}
done &> bamToSam.log &
```

---

### Step 4 — Sorting, indexing, and mitochondrial read removal

```bash
# Sort
for i in *.bam; do
    samtools sort -@ 20 $i > ${i/.bam/_sorted.bam}
done &> bamSorting.log &

# Index
for i in *_sorted.bam; do
    samtools index -@ 40 $i
done

# Remove mitochondrial reads
for i in *_sorted.bam; do
    samtools idxstats -@ 40 $i | cut -f 1 | grep -v chrM | \
    xargs samtools view -b $i > ${i/.bam/_rmMit.bam}
done

# Index
for i in *_rmMit.bam; do
    samtools index -@ 40 $i
done
```

---

### Step 5 — Duplicate removal

```bash
for i in *rmMit.bam; do
    java -jar /home/jbenoit/index/picard/picard.jar MarkDuplicates \
        -I $i \
        -O ${i/.bam/_rmdup.bam} \
        -M ${i/.bam/_rmdup.txt} \
        --REMOVE_DUPLICATES true
done &

for i in *_rmdup.bam; do
    samtools index -@ 40 $i
done
```

Intermediate SAM, `*_sorted.bam`, and `*_rmMit.bam` files were removed after each step to conserve disk space.

---

## CUT&RUN-specific steps

### Step 6 (CUT&RUN) — Spike-in alignment and normalization

CUT&RUN samples were aligned to the *S. cerevisiae* R64-1-1 genome using the same Bowtie2 parameters as the human alignment. Unique yeast reads were counted per sample and used to calculate a normalization factor for downsampling.

```bash
# Align to S. cerevisiae R64-1-1
for i in *_R1_paired*; do
    bowtie2 \
        --end-to-end \
        --no-mixed \
        --no-discordant \
        -p 20 \
        -x /path/to/R64-1-1/index \
        -1 $i \
        -2 ${i/R1/R2} \
        -S ${i/.fastq.gz/_sacCer3.sam}
done &> sacCer3_alignment.log &

# Count unique yeast reads per sample (after BAM filtering as above)
for i in *sacCer3*rmdup.bam; do
    samtools view -c $i
done
```

**Normalization factor calculation:**
```
normalization_factor = min(yeast_reads) / sample_yeast_reads
```

The sample with the fewest unique yeast reads is the reference. All other samples are downsampled by their normalization factor using Picard DownsampleSam:

```bash
java -jar /home/jbenoit/index/picard/picard.jar DownsampleSam \
    -I sample_rmdup.bam \
    -O sample_sub.bam \
    -P <normalization_factor>
```

---

### Step 7 (CUT&RUN) — Peak calling

Peak calling was performed using MACS3 with matched IgG as the control per replicate. Targets were called in narrow or broad mode based on the expected mark distribution. Genome size was set to 3.0e9 for hs1.

```bash
#!/usr/bin/env bash
# Run with MACS3 conda environment active

set -uo pipefail

INDIR="subsampled"
OUTDIR="macs3_out_sub"
GENOME_SIZE="3.0e9"
QVAL="0.01"
BROAD_QVAL="0.05"
BROAD_CUTOFF="0.05"

declare -A MODE=(
    [CTCF]=narrow
    [H3K4me3]=narrow
    [H3K27me3]=broad
    [H3K9me3]=broad
    [HP1a]=broad
    [HP1b]=broad
)

mkdir -p "$OUTDIR"

call_narrow() {
    macs3 callpeak -t "$1" -c "$2" \
        -f BAMPE -g "$GENOME_SIZE" -q "$QVAL" \
        -n "$3" --outdir "$OUTDIR" \
        --keep-dup all --bdg --SPMR
}

call_broad() {
    macs3 callpeak -t "$1" -c "$2" \
        -f BAMPE -g "$GENOME_SIZE" -q "$BROAD_QVAL" \
        --broad --broad-cutoff "$BROAD_CUTOFF" \
        -n "$3" --outdir "$OUTDIR" \
        --keep-dup all --bdg --SPMR
}

shopt -s nullglob
for bam in "$INDIR"/10A_rep*_*_R1_paired_sorted_rmMit_sub.bam; do
    f="$(basename "$bam")"
    rep="${f#10A_}"; rep="${rep%%_*}"
    target="${f#10A_${rep}_}"; target="${target%%_*}"
    [[ "$target" == "IgG" ]] && continue

    ctrl="${INDIR}/10A_${rep}_IgG_R1_paired_sorted_rmMit_sub.bam"
    [[ ! -f "$ctrl" ]] && { echo "ERROR: missing control $ctrl" >&2; continue; }

    mode="${MODE[$target]:-}"
    [[ -z "$mode" ]] && { echo "WARNING: no mode for $target" >&2; continue; }

    base="10A_${rep}_${target}_macs3_sub"
    case "$mode" in
        narrow) call_narrow "$bam" "$ctrl" "$base" ;;
        broad)  call_broad  "$bam" "$ctrl" "${base}_broad" ;;
    esac
done
echo "=== done ==="
```

**Key parameters:**
- `-f BAMPE` — paired-end BAM input
- `--keep-dup all` — retain all reads (duplicates already removed by Picard)
- `--bdg --SPMR` — output bedGraph normalized by sequencing depth
- Narrow peaks: q-value 0.01; Broad peaks: q-value 0.05, broad-cutoff 0.05

---

### Step 8 (CUT&RUN) — Consensus peak sets

Reproducible peaks were identified as regions present in both replicates using bedtools intersect followed by merging of overlapping intervals. Results were validated by concordance with MACS3 bdgdiff common peak calls on replicate bedgraph pairs.

```bash
mk_consensus() {
    local t="$1" ext="$2"
    local a="peaks/hs1_10A_rep1_${t}_peaks.${ext}"
    local b="peaks/hs1_10A_rep2_${t}_peaks.${ext}"
    bedtools intersect -a "$a" -b "$b" -u \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        > "peaks/hs1_10A_${t}_consensus.bed"
    echo "${t} consensus: $(wc -l < "peaks/hs1_10A_${t}_consensus.bed") regions"
}

mk_consensus H3K4me3  narrowPeak
mk_consensus H3K27me3 broadPeak
mk_consensus H3K9me3  broadPeak
```

---

### Step 9 (CUT&RUN) — BigWig generation

CPM-normalized bigWig files at 25 bp resolution. For sparse targets (HP1α, HP1β, H3K9me3), additional log2(target/IgG) tracks were generated at 10 kb resolution for genome browser visualization.

```bash
INDIR="dedup"
OUTDIR="dedup_bw"
THREADS=8
BINSIZE=25

mkdir -p "$OUTDIR"

for bam in "$INDIR"/10A_rep*_dedup.bam; do
    base=$(basename "$bam" _dedup.bam)
    out="${OUTDIR}/${base}.bw"
    [[ -s "$out" ]] && { echo "skip: $out"; continue; }
    [[ -f "${bam}.bai" ]] || samtools index -@ "$THREADS" "$bam"

    bamCoverage -b "$bam" -o "$out" \
        --binSize "$BINSIZE" \
        --normalizeUsing CPM \
        --extendReads \
        --samFlagInclude 2 \
        --numberOfProcessors "$THREADS"
done
```

---

## ChIP-seq-specific steps

### Step 6 (ChIP) — Read depth normalization

ChIP-seq samples were normalized by read depth using Picard DownsampleSam, downsampling all samples to the lowest paired-end fragment count across the dataset (no spike-in normalization).

```bash
# Count paired-end fragments per sample
for i in *rmdup.bam; do
    samtools view -c -f 1 -F 12 $i
done

# Downsample to lowest fragment count
java -jar /home/jbenoit/index/picard/picard.jar DownsampleSam \
    -I sample_rmdup.bam \
    -O sample_downsampled.bam \
    -P <proportion>   # lowest_count / sample_count
```

---

### Step 7 (ChIP) — Peak calling

Broad peaks were called for all ChIP targets using MACS3 with condition-matched input as control. Note that no WT-specific input was available; the closest condition-matched input was used as a proxy (see sample table above).

```bash
# WT HP1α (TS1, TS2) vs shHP1α KD input (TS15, TS16)
macs3 callpeak \
    -t TS1_S1_R1_paired.fastq_q10_sorted_rmMit.bam \
       TS2_S2_R1_paired.fastq_q10_sorted_rmMit.bam \
    -c TS16_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS15_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -f BAMPE -g hs -n ChIP_alpha --broad

# shHP1β KD HP1α (TS3, TS4) vs shHP1β KD input (TS17, TS18)
macs3 callpeak \
    -t TS3_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS4_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -c TS17_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS18_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -f BAMPE -g hs -n KDbeta_ChIP_alpha --broad

# shH2AZ KD HP1α (TS5, TS6) vs shH2AZ KD input (TS19, TS20)
macs3 callpeak \
    -t TS5_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS6_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -c TS19_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS20_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -f BAMPE -g hs -n KDZ_KDB_ChIP_alpha --broad

# WT HP1β (TS7, TS8) vs shHP1β KD input (TS17, TS18)
macs3 callpeak \
    -t TS7_S7_R1_paired.fastq_q10_sorted_rmMit.bam \
       TS8_S8_R1_paired.fastq_q10_sorted_rmMit.bam \
    -c TS17_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS18_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -f BAMPE -g hs -n ChIP_beta --broad

# shHP1α KD HP1β (TS9, TS10) vs shHP1α KD input (TS15, TS16)
macs3 callpeak \
    -t TS9_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS10_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -c TS16_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS15_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -f BAMPE -g hs -n KDalpha_ChIP_beta --broad

# shH2AZ KD HP1β (TS11, TS12) vs shH2AZ KD input (TS19, TS20)
macs3 callpeak \
    -t TS11_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS12_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -c TS17_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
       TS18_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
    -f BAMPE -g hs -n KDZ_KDA_ChIP_beta --broad
```

> **Note:** `-g hs` uses MACS3's built-in human genome size estimate. This is equivalent to the explicit `3.0e9` value used in CUT&RUN peak calling.

---

### Step 8 (ChIP) — BigWig generation

Log2(ChIP/Input) bigWig tracks were generated at 10 bp resolution using deepTools bamCompare with RPKM normalization.

```bash
# WT HP1α
bamCompare -b1 TS1_S1_R1_paired.fastq_q10_sorted_rmMit.bam \
           -b2 TS17_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o HP1alpha_rep1_RPKM_log2ratio.bw -p 15

bamCompare -b1 TS2_S2_R1_paired.fastq_q10_sorted_rmMit.bam \
           -b2 TS17_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o HP1alpha_rep2_RPKM_log2ratio.bw -p 15

# shHP1β KD HP1α
bamCompare -b1 TS3_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -b2 TS17_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o KDB_HP1alpha_rep1_RPKM_log2ratio.bw -p 15

bamCompare -b1 TS4_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -b2 TS17_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o KDB_HP1alpha_rep2_RPKM_log2ratio.bw -p 15

# shH2AZ KD HP1α
bamCompare -b1 TS5_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -b2 TS19_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o KDZ_KDB_HP1alpha_rep1_RPKM_log2ratio.bw -p 15

bamCompare -b1 TS6_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -b2 TS19_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o KDZ_KDB_HP1alpha_rep2_RPKM_log2ratio.bw -p 15

# WT HP1β
bamCompare -b1 TS7_S7_R1_paired.fastq_q10_sorted_rmMit.bam \
           -b2 TS15_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o HP1beta_rep1_RPKM_log2ratio.bw -p 15

bamCompare -b1 TS8_S8_R1_paired.fastq_q10_sorted_rmMit.bam \
           -b2 TS15_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o HP1beta_rep2_RPKM_log2ratio.bw -p 15

# shHP1α KD HP1β
bamCompare -b1 TS9_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -b2 TS15_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o KDA_HP1beta_rep1_RPKM_log2ratio.bw -p 15

bamCompare -b1 TS10_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -b2 TS15_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o KDA_HP1beta_rep2_RPKM_log2ratio.bw -p 15

# shH2AZ KD HP1β
bamCompare -b1 TS11_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -b2 TS19_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o KDZ_KDA_HP1beta_rep1_RPKM_log2ratio.bw -p 15

bamCompare -b1 TS12_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -b2 TS19_R1_paired_merged.fastq_q10_sorted_rmMit.bam \
           -bs 10 --normalizeUsing RPKM --scaleFactorsMethod None \
           --operation log2 --centerReads \
           -o KDZ_KDA_HP1beta_rep2_RPKM_log2ratio.bw -p 15
```

**Key parameters:**
- `--normalizeUsing RPKM` — reads per kilobase per million normalization
- `--scaleFactorsMethod None` — normalization handled by RPKM, not by scaling
- `--operation log2` — log2(ChIP/Input) ratio
- `--centerReads` — center reads on fragment midpoint
- `-bs 10` — 10 bp bin size

---

## Output files used in manuscript

> **Note on genome build:** All files use T2T CHM13v2.0 (hs1) coordinates. Load the hs1 genome in IGV before visualizing tracks.

### CUT&RUN
| File | Type | Description |
|---|---|---|
| `10A_rep1/2_H3K4me3_uniq_peaks.narrowPeak` | narrowPeak | H3K4me3 peaks per replicate |
| `10A_rep1/2_H3K27me3_uniq_broad_peaks.broadPeak` | broadPeak | H3K27me3 peaks per replicate |
| `10A_rep1/2_H3K9me3_uniq_broad_peaks.broadPeak` | broadPeak | H3K9me3 peaks per replicate |
| `10A_rep1/2_HP1a_macs3_broad_peaks.broadPeak` | broadPeak | HP1α peaks per replicate |
| `10A_rep1/2_HP1b_macs3_broad_peaks.broadPeak` | broadPeak | HP1β peaks per replicate |
| `10A_rep1/2_CTCF_uniq_peaks.narrowPeak` | narrowPeak | CTCF peaks per replicate |
| `hs1_10A_*_CPM.bw` | bigWig | CPM-normalized, 25 bp bins |
| `hs1_10A_*_log2IgG_10kb.bw` | bigWig | log2(target/IgG), 10 kb bins |
| `hs1_10A_*_consensus.bed` | BED | Reproducible consensus peak sets |

> **Note on file naming:** The `_uniq_` suffix reflects a naming convention only and does not indicate a distinct filtering step.

### ChIP-seq
| File | Type | Description |
|---|---|---|
| `HP1alpha_rep1/2_RPKM_log2ratio.bw` | bigWig | WT HP1α log2(ChIP/Input) |
| `KDB_HP1alpha_rep1/2_RPKM_log2ratio.bw` | bigWig | shHP1β KD HP1α log2(ChIP/Input) |
| `KDZ_KDB_HP1alpha_rep1/2_RPKM_log2ratio.bw` | bigWig | shH2AZ KD HP1α log2(ChIP/Input) |
| `HP1beta_rep1/2_RPKM_log2ratio.bw` | bigWig | WT HP1β log2(ChIP/Input) |
| `KDA_HP1beta_rep1/2_RPKM_log2ratio.bw` | bigWig | shHP1α KD HP1β log2(ChIP/Input) |
| `KDZ_KDA_HP1beta_rep1/2_RPKM_log2ratio.bw` | bigWig | shH2AZ KD HP1β log2(ChIP/Input) |
