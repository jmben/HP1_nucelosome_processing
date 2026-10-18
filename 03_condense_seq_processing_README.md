Markdown · 03_condense_seq_processing_README.md

# 03 — Condense-seq and Genomic MNase Sensitivity Processing

This section covers processing of two related datasets that both serve as inputs to the chromatin polymer simulation (Section 06):

1. **Condense-seq** — native mononucleosome isolation followed by spermine-induced condensation, sequencing of input (pre-condensation) and supernatant fractions to derive per-nucleosome condensability scores
2. **Genomic MNase sensitivity** — genome-wide (non-sequence-captured) MNase sensitivity profiling in shScramble MCF10A cells, run in parallel with the TSS-restricted mTSS-seq pipeline described in Section 01

Both datasets were aligned to hs1/CHM13v2.0 on the corona2 HPC cluster.

---

## Samples

### Condense-seq

Three libraries were sequenced:

| Label | Sample ID | Description |
| --- | --- | --- |
| 20U_input | Dennis14_S846_L007 | Pre-condensation input, 20U MNase |
| 20U_supernatant_rep1 | Dennis7_S843_L007 | Post-condensation supernatant rep1, 20U MNase |
| 20U_supernatant_rep2 | Dennis8_S844_L007 | Post-condensation supernatant rep2, 20U MNase |

### Genomic MNase sensitivity

Two replicates of shScramble MCF10A cells processed without sequence capture, generating genome-wide nucleosome sensitivity profiles:

| Label | Description |
| --- | --- |
| MCF10A_genomic_rep1 | shScramble clone1, genome-wide MNase sensitivity |
| MCF10A_genomic_rep2 | shScramble clone2, genome-wide MNase sensitivity |

> **Note:** These are non-sequence-captured versions of the shScramble control samples described in Section 01. The same MNase light/heavy digest protocol was used, but without the Roche NimbleGen SeqCap EZ capture step, providing genome-wide coverage rather than TSS-restricted coverage.

---

## Shared upstream pipeline

The trimming, alignment, and BAM processing steps were identical for both condense-seq and genomic MNase sensitivity samples.

### Step 1 — Adapter trimming

Same parameters as mTSS-seq (Section 01):

```bash
for i in *_R1*; do
    java -jar /pauper/jbenoit/index/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -threads 30 \
        -phred33 \
        $i ${i/_R1/_R2} \
        ${i/_R1/_R1_paired} \
        ${i/_R1/_R1_unpaired} \
        ${i/_R1/_R2_paired} \
        ${i/_R1/_R2_unpaired} \
        ILLUMINACLIP:/pauper/jbenoit/index/TruSeq3-PE.fa:2:30:10:1:TRUE \
        MINLEN:25 \
        LEADING:4 \
        TRAILING:4 \
        SLIDINGWINDOW:4:20
done &> trimming.log &
```

---

### Step 2 — Alignment to hs1

```bash
for i in *_R1_paired*; do
    bowtie2 \
        -p 20 \
        -x /pauper/jbenoit/index/chm13v2.0/chm13v2.0 \
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

### Step 4 — Sorting and indexing

```bash
for i in *.bam; do
    samtools sort -@ 40 $i > ${i/.bam/_sorted.bam}
done &> bamSorting.log &

for i in *_sorted.bam; do
    samtools index -@ 40 $i
done &
```

---

### Step 5 — Duplicate removal

```bash
for i in *sorted.bam; do
    java -jar /home/jbenoit/index/picard/picard.jar MarkDuplicates \
        -I $i \
        -O ${i/.bam/_rmdup.bam} \
        -M ${i/.bam/_rmdup.txt} \
        --REMOVE_DUPLICATES true
done &> rmdup_picard.log &

for i in *_rmdup.bam; do
    samtools index -@ 40 $i
done &
```

---

## Condense-seq-specific steps

### Step 6 (condense-seq) — Coverage generation

Per-base coverage BEDGraphs (raw counts) and CPM-normalized bigWig files were generated from deduplicated BAMs using deepTools bamCoverage at 1 bp resolution. Both are required: raw BEDGraphs feed into the condensability score computation; CPM bigWigs are used for genome browser visualization.

```bash
# CPM-normalized bigWig (visualization)
bamCoverage \
    --bam sample_rmdup.bam \
    --outFileName sample.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --normalizeUsing CPM \
    --extendReads \
    --ignoreDuplicates \
    --numberOfProcessors 8

# Raw counts BEDGraph (condensability calculation)
bamCoverage \
    --bam sample_rmdup.bam \
    --outFileName sample.bedgraph \
    --outFileFormat bedgraph \
    --binSize 1 \
    --normalizeUsing None \
    --extendReads \
    --ignoreDuplicates \
    --numberOfProcessors 8
```

> **Note on `--ignoreDuplicates`:** Duplicates were already removed by Picard in Step 5; this flag was retained as a safety measure and has no effect on deduplicated BAMs.

This was run for all three samples via `generate_deposit_bigwigs.sh`:

```bash
BAM_DIR=/data2/jbenoit/condense-seq/bams \
CHROM_SIZES=/data2/jbenoit/condense-seq/hs1.chrom.sizes \
OUTDIR=/data2/jbenoit/condense-seq/deposit_bigwigs \
bash generate_deposit_bigwigs.sh
```

---

### Step 7 (condense-seq) — Nucleosome peak calling and condensability scoring

Condensability scores were computed using `condensability.py`, which:

1. Calls nucleosome peak positions from the input coverage profile using local maxima detection with a minimum inter-peak distance of 107 bp (147 bp footprint minus 40 bp maximum overlap)
2. Applies a greedy non-overlapping selection to enforce the 40 bp maximum overlap constraint
3. Computes per-nucleosome condensability as: `S = −ln((fraction_area + 1) / (input_area + 1))` where areas are summed over a 171 bp window centered on each dyad
4. Also computes sliding window (171 bp / 25 bp step), 1 kb, and 10 kb binned scores
5. Outputs BEDGraph and bigWig files at all four resolutions

```bash
python3 condensability.py \
    --input-bedgraph    condensability_output/coverage/20U_input.bedgraph \
    --fraction-bedgraph condensability_output/coverage/20U_supernatant_rep1.bedgraph \
    --chrom-sizes       chm13v2.0.chrom.sizes \
    --titration         20U \
    --fraction          supernatant_rep1 \
    --outdir            condensability_output \
    --threads           8
```

This was run for both 20U supernatant replicates vs the 20U input:

- 20U supernatant rep1 (Dennis7) vs 20U input (Dennis14)
- 20U supernatant rep2 (Dennis8, downsampled) vs 20U input (Dennis14)

**Downsampling (20U rep2):** Dennis8 had approximately 2× more fragments than the 20U input prior to condensability scoring. It was downsampled to match the input fragment count (136,334,240 fragments) using Picard DownsampleSam before BEDGraph generation to avoid a depth-driven sign artifact:

```bash
java -jar /home/jbenoit/index/picard/picard.jar DownsampleSam \
    -I Dennis8_S844_L007_R1_paired_001_sorted_rmdup.bam \
    -O Dennis8_S844_L007_R1_paired_001_sorted_rmdup_ds.bam \
    -P <proportion>   # 136334240 / Dennis8_fragment_count
```

Both replicates are handled automatically by `generate_deposit_bigwigs.sh`.

**Key parameters (from paper):**

- Nucleosome footprint: 147 bp
- Coverage window for area calculation: 171 bp
- Maximum inter-nucleosome overlap: 40 bp
- Sliding window step: 25 bp
- Bin sizes: 1 kb and 10 kb
- Pseudocount: 1 (added to both numerator and denominator)

---

### Step 8 (condense-seq) — Sign orientation

Condensability bigWig files were inverted (multiplied by −1) where needed to ensure consistent orientation across tracks for downstream correlation analysis and polymer simulation input. This was performed using `invert_bigwig.py`:

```bash
python3 invert_bigwig.py input.bw input_inverted.bw
```

> **Sign convention:** `condensability.py` computes `S = −ln(supernatant/input)`, which by default gives negative values in A compartment (euchromatin). All condensability bigWig files were inverted after computation using `invert_bigwig.py` to produce the final convention used throughout this study: **positive values = A compartment (euchromatin, nucleosomes enriched in supernatant = not condensed); negative values = B compartment (heterochromatin, nucleosomes depleted in supernatant = condensed)**. This convention was validated by correlation with MCF10A Hi-C compartment scores (GSE246947) and is consistent with the published GM12878 condense-seq reference track (GSE252941).

---

## Genomic MNase sensitivity-specific steps

### Step 6 (genomic sensitivity) — Read depth normalization and merging

Fragment counts were checked per sample, then light and heavy digest BAMs were downsampled to matching depth using Picard DownsampleSam (same approach as Section 01). Downsampled light and heavy BAMs were then merged using Picard MergeSamFiles to produce total occupancy files.

```bash
# Count paired-end fragments per sample
for i in *rmdup.bam; do
    samtools view -c -f 1 -F 12 $i
done

# Downsample light and heavy to matching fragment count
java -jar /home/jbenoit/index/picard/picard.jar DownsampleSam \
    -I sample_light_rmdup.bam \
    -O sample_light_downsampled.bam \
    -P <proportion>   # lowest_count / sample_count

# Merge downsampled light and heavy into total occupancy BAM
java -jar /home/jbenoit/index/picard/picard.jar MergeSamFiles \
    I=sample_light_downsampled.bam \
    I=sample_heavy_downsampled.bam \
    O=sample_merged_total.bam

# Index merged BAM
for i in *merged*.bam; do
    samtools index -@ 40 $i
done
```

> **Note:** Unlike the TSS-restricted pipeline (Section 01), no restriction to TSS windows was applied — all reads mapping to the hs1 genome were retained for genome-wide analysis.

---

### Step 7 (genomic sensitivity) — DANPOS3 sensitivity calling

Genome-wide nucleosome sensitivity was called using DANPOS3 following the same approach as the TSS-restricted pipeline (Section 01), but applied genome-wide without TSS restriction. Light digest BAM is the treatment and heavy digest BAM is the control.

```bash
# Example — shScramble clone1
python3 danpos.py dpos \
    path/to/shScr_clone1_light_downsampled.bam:path/to/shScr_clone1_heavy_downsampled.bam \
    -m 1 \
    --mifrsz 110 \
    --mafrsz 200 \
    -o MCF10A_genomic_rep1_sensitivity
```

DANPOS `.wig` output was converted to bigWig:

```bash
./wigToBigWig name_of_file.wig hs1.chrom.sizes \
    MCF10A_genomic_rep1_downsampled_sensitivity.bw
```

The log2 ratio of light to heavy digest signal was computed using deepTools bigwigCompare:

```bash
bigwigCompare \
    --binSize 1 \
    --outFileName MCF10A_genomic_rep1_sensitivity.bw \
    --outFileFormat bigwig \
    --bigwig1 light_downsampled.bw \
    --bigwig2 heavy_downsampled.bw \
    --skipZeroOverZero \
    --operation log2 \
    --skipNonCoveredRegions
```

Two replicates were processed (shScramble clone1 and clone2), producing:

- `MCF10A_genomic_rep1_downsampled_sensitivity.bw`
- `MCF10A_genomic_rep2_downsampled_sensitivity.bw`

---

## Output files

### Condense-seq (intermediate files)

| File | Description |
| --- | --- |
| `coverage/*_CPM.bw` | CPM-normalized coverage bigWig, 1 bp bins |
| `coverage/*.bedgraph` | Raw count coverage BEDGraph, 1 bp bins |
| `nucleosomes/*_nucleosomes.tsv` | Per-nucleosome peak table with coordinates and scores |
| `bigwig/*_per_nuc.bw` | Condensability score per nucleosome dyad |
| `bigwig/*_slide171_25.bw` | Sliding window condensability, 171 bp / 25 bp step |
| `bigwig/*_1kb.bw` | Condensability binned at 1 kb (pre-inversion) |
| `bigwig/*_10kb.bw` | Condensability binned at 10 kb (pre-inversion) |

### GEO deposit files

Three processed data files (CPM-normalized bigWigs) were deposited to NCBI GEO, one per library:

| GEO library name | Processed data file | Raw files |
| --- | --- | --- |
| MCF10A_condense-seq_input | `MCF10A_condenseseq_20U_input_CPM.bw` | Dennis14_S846_L007_R1/R2_001.fastq.gz |
| MCF10A_condense-seq_rep1 | `MCF10A_condenseseq_20U_supernatant_rep1_CPM.bw` | Dennis7_S843_L007_R1/R2_001.fastq.gz |
| MCF10A_condense-seq_rep2 | `MCF10A_condenseseq_20U_supernatant_rep2_CPM.bw` | Dennis8_S844_L007_R1/R2_001.fastq.gz |

Condensability score bigWigs (1 kb bins, inverted) are provided in the code repository rather than GEO, as they are derived analysis outputs rather than primary processed data.

### Genomic MNase sensitivity

| File | Description |
| --- | --- |
| `MCF10A_genomic_rep1_downsampled_sensitivity.bw` | Genome-wide sensitivity rep1 |
| `MCF10A_genomic_rep2_downsampled_sensitivity.bw` | Genome-wide sensitivity rep2 |

---

## Downstream use

Both condense-seq and genomic MNase sensitivity bigWig files were used as inputs to the chromatin polymer simulation pipeline (Section 06). The condensability φ vector for each genomic region was derived by binning these tracks at 25 kb resolution using `01_bin_bigwig.py` as described in the simulation README.

Genome-wide correlation analysis between condense-seq tracks, MNase sensitivity tracks, and MCF10A Hi-C compartment scores was performed using `correlation_analysis.py`, `pairwise_correlation.py`, and `plot_correlation_forest.py` (Section 06).