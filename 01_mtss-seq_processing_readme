# 01 — mTSS-seq Processing Pipeline

This section covers the full processing pipeline for micrococcal nuclease TSS-seq (mTSS-seq) data, from raw reads through nucleosome position calling with DANPOS3. All steps were run on pauper2 HPC.

---

## Overview of steps

1. Adapter trimming (Trimmomatic)
2. Alignment to hs1 (Bowtie2)
3. SAM → BAM conversion and quality filtering (Samtools)
4. Sorting, indexing, and mitochondrial read removal (Samtools)
5. Duplicate removal (Picard MarkDuplicates)
6. TSS restriction (Samtools + T2T TSS BED file)
7. Read depth normalization by downsampling (Samtools fragment counting + Picard DownsampleSam)
8. Nucleosome position calling (DANPOS3)

---

## Input files

- Paired-end FASTQ files (Illumina, PE50)
- Sequenced on Illumina NovaSeq 6000 at FSU College of Medicine
- Target depth: >20 million paired-end reads per sample (~20× coverage of ~42 Mb captured promoter space)
- Samples include: shScramble (clone1, clone2), shHP1α (clone5, clone6/cloneE), shHP1β (cloneD, clone15), shH2AZ

---

## Reference files

- Genome: T2T CHM13v2.0 (hs1), Bowtie2 index at `chm13v2.0/chm13v2.0`
- TSS BED file: `T2T_TSS_1kb_reduced_clean.bed` (±1 kb flanking all human gene TSSs)
- Adapter sequences: `TruSeq3-PE.fa`

---

## Step 1 — Adapter trimming

Paired-end adapter trimming with moderate stringency. Parameters were chosen to remove low-quality bases while retaining sufficient read length for accurate alignment.

```bash
for i in *_R1*; do
    java -jar /home/bioinformatics/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -threads 30 \
        -phred33 \
        $i ${i/_R1/_R2} \
        ${i/_R1/_R1_paired} \
        ${i/_R1/_R1_unpaired} \
        ${i/_R1/_R2_paired} \
        ${i/_R1/_R2_unpaired} \
        ILLUMINACLIP:/home/jbenoit/index/TruSeq3-PE.fa:2:30:10:1:TRUE \
        MINLEN:25 \
        LEADING:4 \
        TRAILING:4 \
        SLIDINGWINDOW:4:20
done &> trimming.log &
```

**Key parameters:**
- `LEADING:4 TRAILING:4` — remove bases with quality < 4 from read ends
- `SLIDINGWINDOW:4:20` — trim when average quality in a 4-base window falls below 20
- `MINLEN:25` — discard reads shorter than 25 bp after trimming
- Only paired output files (`*_paired`) are carried forward

---

## Step 2 — Alignment

Paired-end alignment to the T2T human genome (hs1/CHM13v2.0). Flags enforce concordant paired-end mapping only.

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

**Key flags:**
- `--end-to-end` — full read must align (no soft clipping)
- `--no-mixed` — both reads in a pair must align
- `--no-discordant` — discard pairs that align in discordant orientation

---

## Step 3 — SAM to BAM conversion and quality filtering

Convert SAM to BAM and discard reads with MAPQ < 10 (multi-mapping reads).

```bash
for i in *sam; do
    samtools view -b -q 10 $i > ${i/.sam/.bam}
done &> bamToSam.log
```

---

## Step 4 — Sorting, indexing, and mitochondrial read removal

```bash
# Sort BAM files
for i in *.bam; do
    samtools sort -@ 40 $i > ${i/.bam/_sorted.bam}
done &> bamSorting.log &

# Index sorted BAM files
for i in *_sorted.bam; do
    samtools index -@ 40 $i
done &

# Remove mitochondrial reads by excluding chrM from index
for i in *_sorted.bam; do
    samtools idxstats -@ 40 $i | cut -f 1 | grep -v chrM | \
    xargs samtools view -b $i > ${i/.bam/_rmMit.bam}
done &> rmMit.log &

# Index mitochondrial-filtered BAM files
for i in *_rmMit.bam; do
    samtools index -@ 40 $i
done &
```

---

## Step 5 — Duplicate removal

Duplicates were marked and removed using Picard MarkDuplicates. Three BAM file types were deduplicated independently:
- Light digest files (MNase 2U, 6 min)
- Heavy digest files (MNase 20U, 10 min)
- Merged files (light + heavy combined for total occupancy analysis)

```bash
for i in *rmMit.bam; do
    java -jar /home/jbenoit/index/picard/picard.jar MarkDuplicates \
        -I $i \
        -O ${i/.bam/_rmdup.bam} \
        -M ${i/.bam/_rmdup.txt} \
        --REMOVE_DUPLICATES true
done &> rmdup_picard.log &

# Index deduplicated files
for i in *_rmdup.bam; do
    samtools index -@ 40 $i
done &
```

---

## Step 6 — TSS restriction

Restrict BAM files to reads mapping within ±1 kb of all human gene TSSs (sequence-captured regions).

```bash
for i in *rmdup.bam; do
    samtools view -b -L T2T_TSS_1kb_reduced_clean.bed $i > ${i/.bam/_TSS.bam}
done &
```

---

## Step 7 — Read depth normalization

Samples were downsampled to the lowest fragment count across all samples to ensure uniform nucleosome calling depth.

**Count paired-end fragments per TSS-restricted BAM:**

```bash
for i in *TSS*bam; do
    samtools view -c -f 1 -F 12 $i
done
```

- `-f 1` — retain only paired reads
- `-F 12` — exclude unmapped reads and reads with unmapped mates

**Downsample to lowest fragment count using Picard DownsampleSam:**

The downsampling proportion (`-P`) was calculated per sample as:
`P = lowest_fragment_count / sample_fragment_count`

```bash
java -jar /home/jbenoit/index/picard/picard.jar DownsampleSam \
    -I sample_TSS.bam \
    -O sample_TSS_downsampled.bam \
    -P <proportion>   # calculated per sample as above
```

---

## Step 8 — Nucleosome position calling (DANPOS3)

DANPOS3 was used for two types of analysis:

### 8a — Nucleosome MNase sensitivity (light vs. heavy digest)

Light digest BAM (2U) is the treatment; heavy digest BAM (20U) is the control. DANPOS computes per-nucleosome summit score log2 fold change (light/heavy), which reflects MNase sensitivity. Fragment size cutoffs (`--mifrsz 110 --mafrsz 200`) were applied consistently across all samples to restrict analysis to nucleosome-sized fragments.

```bash
# Example — shScramble clone1
python3 danpos.py dpos \
    path/to/MCF10A_shScr_clone1_light_TSS_downsampled.bam:path/to/MCF10A_shScr_clone1_heavy_TSS_downsampled.bam \
    -m 1 \
    --mifrsz 110 \
    --mafrsz 200 \
    -o MCF10A_final/10A_shScr_clone1_L_vs_H_downsampled
```

Samples run with this approach:
- shScramble clone1, clone2
- shHP1α clone5, clone6
- shHP1β cloneD, clone15

### 8b — Nucleosome occupancy changes (merged KD vs. scramble)

Merged (total occupancy) KD BAM is the treatment; matched scramble BAM is the control. This identifies nucleosomes with significant occupancy changes between conditions.

```bash
# Example — shHP1α clone5 vs shScr clone2
python3 danpos.py dpos \
    path/to/MCF10A_shHP1a_clone5_merged_TSS_downsampled.bam:path/to/MCF10A_shScr_clone2_merged_TSS_downsampled.bam \
    -m 1 \
    --mifrsz 110 \
    --mafrsz 200 \
    -o MCF10A_final/MCF10A_shHP1a_clone5_vs_scr2_occupancy_downsampled
```

Comparisons run:
- shHP1α clone5 vs shScr clone2
- shHP1α clone6/cloneE vs shScr clone1/clone2
- shHP1β cloneD vs shScr clone1
- shHP1β clone15 vs shScr clone2

**Key DANPOS parameters (consistent across all runs):**
- `-m 1` — paired-end mode
- `--mifrsz 110` — minimum fragment size 110 bp (excludes sub-nucleosomal fragments)
- `--mafrsz 200` — maximum fragment size 200 bp (excludes di-nucleosomes and above)
- The colon (`:`) between BAM files directs DANPOS to compare the second file (control) to the first (treatment)

---

## Output files

DANPOS produces `.positions.integrative.xls` files containing per-nucleosome statistics including:
- Summit position and displacement (`treat2control_dis`)
- Summit score log2 fold change (`smt_log2FC`)
- Statistical significance (`smt_diff_log10pval`, `smt_diff_FDR`)
- Fuzziness metrics

These files are used as input for the nucleosome sensitivity classification scripts in `04_nucleosome_sensitivity/`.

---

## Downstream processing

DANPOS `.wig` output files were converted to BigWig format for genome browser visualization and deepTools analysis:

```bash
./wigToBigWig name_of_file.wig hs1.chrom.sizes name_of_output_file.bw
```

The log2 ratio of light to heavy digest signal was computed using deepTools bigwigCompare:

```bash
bigwigCompare \
    --binSize 1 \
    --outFileName output.bw \
    --outFileFormat bigwig \
    --bigwig1 light_digest.bw \
    --bigwig2 heavy_digest.bw \
    --skipZeroOverZero \
    --operation log2 \
    --skipNonCoveredRegions
```
