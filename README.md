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
├── README.md                          # This file
├── 01_mTSS_seq_processing/
│   ├── README.md                      # Detailed notes for mTSS-seq pipeline
│   └── mTSS_seq_processing.sh         # Full bash pipeline: QC → alignment → DANPOS
├── 02_ChIP_CUT_and_RUN_processing/
│   ├── README.md                      # Detailed notes for ChIP/CUT&RUN pipeline
│   └── ChIP_CUTandRUN_processing.sh   # Full bash pipeline: QC → alignment → peak calling
├── 03_condense_seq_processing/
│   ├── README.md                      # Detailed notes for condense-seq pipeline
│   └── condense_seq_processing.sh     # Nucleosome isolation → condensation → library prep pipeline
├── 04_nucleosome_sensitivity/
│   ├── README.md                      # Detailed notes for sensitivity classification
│   ├── nucleosome_sensitivity_classification.R   # DANPOS output parsing and S/R/N classification
│   └── nucleosome_sensitivity_figures.R          # Sensitivity proportion and count figure scripts
├── 05_HP1_peak_analysis/
│   ├── README.md                      # Detailed notes for HP1/H3K9me3 intersection analysis
│   ├── HP1_H3K9me3_intersect.sh       # Genome-wide HP1 vs H3K9me3 overlap quantification
│   ├── HP1_H3K9me3_promoter_intersect.sh  # Promoter-restricted HP1 vs H3K9me3 overlap
│   ├── HP1_peak_signal_extract.sh     # Peak signal extraction by H3K9me3 overlap status
│   └── HP1_peak_figures.R             # HP1 overlap and signal distribution figures
├── 06_chromatin_polymer_simulation/
│   ├── README.md                      # Detailed notes for polymer simulation
│   ├── 01_build_condensability_vector.py   # Condense-seq/MNase signal binning and normalization
│   ├── 02_build_energy_matrix.py          # Intermediate energy matrix construction
│   ├── 03_run_simulation.py               # OpenMM simulation execution
│   └── 04_contact_matrix_analysis.py      # Contact map generation, TAD calling, compartment PCA
├── 07_GO_enrichment/
│   ├── README.md                      # Detailed notes for GO analysis
│   └── HP1_cluster_GO.R               # GO Biological Process enrichment across nucleosome clusters
└── 08_figure_scripts/
    ├── README.md                      # Figure generation notes
    ├── nucleosome_distribution_barplot.R         # Figure: nucleosome distribution changes
    ├── nucleosome_sensitivity_HP1_dualaxis.R     # Figure: sensitivity proportions shScr vs shHP1
    ├── nucleosome_sensitivity_H2AZ_combined.R    # Figure: sensitivity proportions shScr vs shH2AZ
    ├── HP1_promoter_H3K9me3_twopanel.R           # Figure: promoter HP1 peak counts and H3K9me3 overlap
    └── HP1_WT_KD_H3K9me3_venn.svg               # Figure: three-circle Venn HP1/H3K9me3 overlap
```

---

## Dependencies

### Bioinformatics tools
All tools were run on an HPC cluster (r1pl-hpcf, SLURM; pauper2; corona2) using the following versions:

| Tool | Version | Use |
|---|---|---|
| Trimmomatic | 0.39 | Adapter trimming |
| Bowtie2 | — | Alignment to hs1 |
| Samtools | 1.10 | BAM processing and filtering |
| Picard | 2.27.3 | Duplicate removal, downsampling, merging |
| DANPOS | 3.1.1 | Nucleosome position calling |
| deepTools | 3.5.1 | Coverage computation, log2 ratio, heatmaps |
| MACS3 | — | ChIP/CUT&RUN peak calling |
| bedtools | — | Genomic interval operations |
| UCSC kentUtils | — | wigToBigWig, liftOver |
| OpenMM | 8.2 | Chromatin polymer simulation (NVIDIA A100 GPU) |
| cooler | — | Hi-C contact matrix extraction |
| pyBigWig | — | BigWig file reading |
| SciPy | — | Nucleosome dyad peak calling |

### R packages
| Package | Version | Use |
|---|---|---|
| ggplot2 | — | Figures |
| ggpattern | — | Hatched bar charts |
| dplyr / tidyr | — | Data manipulation |
| patchwork | — | Multi-panel figures |
| clusterProfiler | — | GO enrichment analysis |
| org.Hs.eg.db | — | Gene ID mapping |
| PlotGardener | 1.4.1 | Genomic loci visualization |
| VennDiagram | — | Venn diagrams |

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

Raw and processed sequencing data generated in this study are deposited in NCBI GEO under accession **[GSE###### — to be assigned upon acceptance]**.

---

## Contact

For questions about the code or analysis, please contact:  
Jane Benoit (jane.benoit@nationwidechildrens.org)
Bishop Lab, Nationwide Children's Hospital / Formerly Dennis Lab, Florida State University
