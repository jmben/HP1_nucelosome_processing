# HP1_nucelosome_processing
Scripts for processing MNase-seq and multi-omics related to HP1 and H2A.Z in MCF10A breast epithelial cells


# Code Repository: Chromatin Organization and HP1 Isoform Function in MCF10A Cells

**Associated manuscript:** Cooperative recruitment of HP1 by H2A.Z mediates gene silencing through MNase-sensitive higher-order chromatin architecture
**Authors:** Jane Benoit2,3, Yasmin Dijkwel1, Tatiana Soboleva1, Jonathan Dennis2# and David John Tremethick1#
---

## Overview

This repository contains all custom analysis scripts used to process, analyze, and visualize data in the associated manuscript. The study characterizes chromatin organization, nucleosome dynamics, and HP1 isoform (HP1α and HP1β) function in MCF10A human breast epithelial cells, with a focus on the relationship between H2A.Z, HP1, and H3K9me3 at gene promoters.

All sequencing data were aligned to the T2T human genome assembly (hs1/CHM13v2.0). Where publicly available datasets were reprocessed, identical pipelines were applied.

---

## Repository Structure

```
├── README.md                                           # This file
├── 01_mTSS_seq_processing_README.md                    # Detailed notes for mTSS-seq pipeline
├── 02_ChIP_CUT_and_RUN_processing_README.md            # Detailed notes for ChIP/CUT&RUN pipeline
├── 03_condense_seq_processing_README.md                # Detailed notes for condense-seq pipeline
├── 04_nucleosome_sensitivity_README.md                 # Detailed notes for sensitivity classification
├── 05_HP1_peak_analysis_README.md                      # Detailed notes for HP1/H3K9me3 intersection analysis
├── 06_chromatin_polymer_simulation_README.md           # Detailed notes for polymer simulation
├── 07_GO_enrichment_README.md                          # Detailed notes for GO analysis
└── scripts
    ├── individual scripts for each of the above processes 
```

---

## Dependencies

### Bioinformatics tools
All tools were run on an HPC cluster using the following versions:

| Tool | Version | Use |
|---|---|---|
| Trimmomatic | 0.39 | Adapter trimming |
| Bowtie2 | 2.5.2 | Alignment to hs1 |
| Samtools | 1.10 | BAM processing and filtering |
| Picard | 2.27.3 | Duplicate removal, downsampling, merging |
| DANPOS | 3.1.1 | Nucleosome position calling |
| deepTools | 3.5.1 | Coverage computation, log2 ratio, heatmaps |
| MACS3 | 3.0.4 | ChIP/CUT&RUN peak calling |
| bedtools | v2.31.1 | Genomic interval operations |
| UCSC kentUtils | — | wigToBigWig, liftOver |
| OpenMM | 8.2 | Chromatin polymer simulation (NVIDIA A100 GPU) |
| cooler | 0.10.4 | Hi-C contact matrix extraction |
| pyBigWig | 0.3.25 | BigWig file reading |
| SciPy | 1.15.2 | Nucleosome dyad peak calling |

### R 4.5.2 packages
| Package | Version | Use |
|---|---|---|
| ggplot2 | 4.0.3 | Figures |
| ggpattern | 1.3.1 | Hatched bar charts |
| dplyr / tidyr | 1.2.1/1.3.2 | Data manipulation |
| patchwork | 1.3.2 | Multi-panel figures |
| clusterProfiler | 4.18.4 | GO enrichment analysis |
| org.Hs.eg.db | 3.22.0 | Gene ID mapping |
| PlotGardener | 1.4.1 | Genomic loci visualization |
| VennDiagram | 1.8.2 | Venn diagrams |

### Python packages
| Package | Use |
|---|---|
| pyBigWig | BigWig signal extraction |
| SciPy | Peak finding |
| NumPy | Array operations |

---

## Reference Genome

All analyses were performed using the T2T human genome assembly **hs1 (CHM13v2.0)**. Publicly available datasets originally aligned to hg38 were lifted over to hs1 using UCSC liftOver prior to analysis.

Genome files used:
- Bowtie2 index: `chm13v2.0`
- Chromosome sizes: `hs1.chrom.sizes`
- Subcompartment annotations: hs1 A0–A3, B0–B3 (provided in repository)
- TSS reference: `T2T_TSS_1kb_reduced_clean.bed` (provided in repository)

---

## Publicly Available Datasets

| Dataset | Accession | Use |
|---|---|---|
| MCF10A Hi-C | GSE246947 | Compartment scores, TAD boundaries |
| GM12878 condense-seq | GSE252941 | Genome-wide correlation analysis |
| MCF10A histone ChIP-seq | PRJNA336352 | H3K27ac, H3K27me3, H3K9me3, H3K4me3 |
| MCF10A ATAC-seq | SRR12006157 | Genome-wide correlation analysis |

---

## Data Availability

Raw and processed sequencing data generated in this study are deposited in NCBI GEO

---

## Contact

For questions about the code or analysis, please contact:  
Jane Benoit (jane.benoit@nationwidechildrens.org)
Bishop Lab, Nationwide Children's Hospital / Formerly Dennis Lab, Florida State University
