Markdown · 06_chromatin_polymer_simulation_README_v2.md

# 06 — Chromatin Polymer Simulation and Condensability Analysis

Chromatin organization was modeled using a bead-spring polymer simulation in which each bead represents a 25 kb genomic bin. The condensability (φ) of each bead was derived from experimental epigenomic tracks — condense-seq supernatant replicates and genomic MNase sensitivity — and used to drive differential bead–bead attraction potentials. This section also covers the upstream computation of condensability scores from raw BEDGraph files, track orientation validation, and genome-wide correlation analyses linking condensability and nucleosome sensitivity to Hi-C compartmentalization. Simulations were run on NVIDIA A100 GPUs via SLURM job arrays on a GPU cluster.

---

## Scripts

### Upstream condensability processing

| Script | Purpose |
| --- | --- |
| `condensability.py` | Compute genome-wide condensability scores at multiple resolutions from raw BEDGraph files; call nucleosome peaks; output bigWig tracks |
| `invert_bigwig.py` | Negate all values in a bigWig file; used to produce `*_inverted.bw` files for tracks that must be sign-corrected before use |
| `check_condensability_sign.py` | Diagnostic: verify that the condensability track orientation is consistent with Hi-C compartment sign convention |

### Simulation pipeline

| Script | Purpose |
| --- | --- |
| `01_bin_bigwig.py` | Bin condense-seq and MNase sensitivity bigWig tracks at a chosen resolution (default 10 kb); output φ vectors as `.npz` |
| `02_build_energy_matrix.py` | Convert binned φ vectors to per-bead pairwise interaction energy matrices; output `.npz` energy files |
| `03_run_simulation.py` | Run bead-spring polymer simulation on GPU (OpenMM); output per-trajectory contact matrices |
| `04_analyse_results.py` | Compute average contact map, O/E normalisation, A/B compartments by PCA, and summary statistics from simulation trajectories |

### Simulation visualization and comparison

| Script | Purpose |
| --- | --- |
| `plot_triangle.py` | Plot triangle contact map from a single simulation `contact_matrix_avg.npy` |
| `alt_triangle_stacked.py` | Stack Hi-C and simulation triangle contact maps for visual comparison |
| `fancy_contacts.py` | Quantitative comparison of simulation vs Hi-C: P(s) curves, SCC, Pearson/Spearman on O/E and rank-normalised matrices, difference maps |
| `plot_split_triangle.py` | Split triangle figure (upper = condense-seq supernatant, lower = MNase sensitivity) with raw and O/E versions |
| `mcool_triangle.py` | Plot triangle map directly from `.mcool` Hi-C file |
| `TADS_cleaner.py` | TAD boundary calling via insulation scores (Crane et al. 2015) and overlap analysis between simulation and Hi-C |

### Track orientation validation

| Script | Purpose |
| --- | --- |
| `compare_tracks3.py` | Verify sign convention of each input track vs Hi-C compartment reference; uses `*_inverted.bw` tracks with flip=False for all; reference `MCF10A_subcompartments_20000_T2T.bw` (positive = B) |

### Genome-wide correlation analyses

| Script | Purpose |
| --- | --- |
| `correlation_analysis.py` | Per-chromosome Pearson r at 1 Mb bins; genome-wide forest plot (Pearson r + Spearman ρ + 95% CI); hexbin scatter plots per assay category |
| `pairwise_correlation.py` | Pairwise Pearson r between all tracks; hierarchically clustered heatmap; forest plot comparing each track to a common reference |
| `plot_correlation_forest.py` | Standalone forest plot from pre-computed correlation values; Pearson r + Spearman ρ vs Hi-C compartment score |
| `multi-track_correlation_clusters.py` | Genome-wide Spearman correlation matrix at 10 kb + TSS profiles stratified by deepTools cluster (7-cluster BED `HP1_HPTM_2kb_sen_matrix_og_sort_k7.bed`); per-cluster grid, per-track overlay, per-cluster correlation matrices |
| `replot_cluster_grid.py` | Re-renders the per-cluster TSS profile grid directly from bigWig files without recomputing intermediate matrices |
| `subcompartment_analysis.py` | Assigns 1 Mb bins to subcompartments A0–B3 by dominant overlap; violin plots, heatmap, and line profiles per track per subcompartment |

### Nucleosome sensitivity vs compartment analyses

| Script | Purpose |
| --- | --- |
| `genomic_sensitivity_vs_compartment.py` | Genome-wide nucleosome sensitivity counts per subcompartment A0–B3; chi-squared test per subcompartment vs genome-wide average; Spearman ρ across subcompartment rank; grouped bar, proportional stacked bar, and fraction profile plots |
| `nucleosome_sensitivity_vs_compartment.py` | TSS-level nucleosome sensitivity assigned to subcompartments; gene-level promoter class assignment (SR score = (S−R)/(S+R)); 4 classes (Sens/Sens, Sens/Res, Res/Res, Res/Sens); stacked bars, fraction lines, heatmap, and grouped bar plots |
| `promoter_class_proportions.py` | Promoter sensitivity class proportions (%) normalised within each subcompartment; chi-squared significance stars; Scr / HP1α KD / HP1β KD comparison |
| `nuc-summary_revised.R` | Nucleosome sensitivity counts per TSS file for HP1α/HP1β KD vs scramble; deduplicates to unique nucleosomes; stacked bar plots (raw counts, proportions, and combined mean ± SD across clones) |
| `H2AZ_nucleosome_sensitivity.R` | Same structure for H2AZ KD vs scramble and HP1β KD samples |

---

## Workflow

### Upstream Step — Compute condensability scores (`condensability.py`)

Raw condense-seq BEDGraph files (supernatant fractions and input) were converted to condensability scores at four resolutions. Nucleosome peaks were called from the input track first, and condensability was computed per peak and at sliding-window and binned resolutions:

```bash
python condensability.py \
    --pellet  data/Dennis8_supernatant.bedgraph \
    --input   data/Dennis14_input.bedgraph \
    --chrom_sizes  hs1.chrom.sizes \
    --out_prefix   condensability_output/20U
```

Condensability score formula:

```
S = ln((fraction_count + pseudocount) / (input_count + pseudocount))
```

Nucleosome peak calling used `scipy.signal.find_peaks` with `distance=107 bp` (NUC_SIZE − MAX_OVERLAP = 147 − 40), followed by a greedy overlap filter. Output resolutions:

| Output | Description |
| --- | --- |
| `*_per_nuc.bw` | Per-nucleosome condensability (171 bp half-window) |
| `*_slide171_25.bw` | Sliding window 171 bp / 25 bp step |
| `*_1kb.bw` | 1 kb bins |
| `*_10kb.bw` | 10 kb bins |

BigWig output used `bedGraphToBigWig` (UCSC kentUtils) with `pyBigWig` as fallback.

Sample mapping used in production:

| Dennis ID | Sample |
| --- | --- |
| Dennis8 | 20U supernatant rep1 |
| Dennis7 | 20U supernatant rep2 |
| Dennis14 | 20U input |

### Upstream Step — Invert bigWig files (`invert_bigwig.py`)

MNase sensitivity tracks are positively correlated with open/active chromatin (A compartment). Because the simulation's condensability convention encodes B-compartment enrichment as positive, sensitivity tracks were sign-inverted before use:

```bash
python invert_bigwig.py \
    --in   data/20U_supernatant_rep1_1kb.bw \
    --out  data/20U_supernatant_rep1_1kb_inverted.bw
```

This reads all intervals from the input bigWig, negates their values, and writes the result. The `*_inverted.bw` files are used as input to `correlation_analysis.py` and `subcompartment_analysis.py`.

### Upstream Step — Validate condensability sign (`check_condensability_sign.py`)

Before simulation, the orientation of each condensability track relative to the Hi-C compartment score was verified:

```bash
python check_condensability_sign.py \
    --numerator   data/Dennis8_supernatant.bedgraph \
    --denominator data/Dennis14_input.bedgraph \
    --compartment MCF10A_subcompartments_20000_T2T.bw \
    --chrom chr2 --start 10000000 --end 14000000
```

Computes `S = ln((numerator+1)/(denominator+1))` and reports its Pearson r with the compartment track. Prints an interpretation guide indicating whether the sign is consistent with expectations.

---

### Step 1 — Bin epigenomic tracks (`01_bin_bigwig.py`)

Condense-seq pellet and genomic MNase sensitivity bigWig files were binned over target chromosomal regions. Each bin value was averaged across the bin and output as a NumPy `.npz` file containing the φ vector, raw scores, bin start coordinates, bin size, chromosome, and region boundaries:

```bash
python 01_bin_bigwig.py \
    --bigwig  data/20U_supernatant_rep1_1kb.bw \
    --chrom   chr2 \
    --start   10000000 \
    --end     14000000 \
    --binsize 25000 \
    --out     data/phi_20U_rep1_chr2_4mb_25kb.npz
```

Use `--flip` to negate values (e.g. for MNase sensitivity tracks that must be inverted). The default bin size is 10 kb; production SLURM scripts used 25 kb.

The φ vector normalisation applied in the script:

```
q   = (p − min p) / (max p − min p)
φ   = q − mean(q) + φ̄          where φ̄ = 0.6
```

Three φ tracks were used as simulation inputs:

| Track label | Source |
| --- | --- |
| `rep1` | 20U condense-seq supernatant replicate 1 |
| `rep2` | 20U condense-seq supernatant replicate 2 |
| `mnase` | Genomic MNase sensitivity |

### Step 2 — Build energy matrices (`02_build_energy_matrix.py`)

The binned φ vectors were converted to pairwise bead interaction energies using an outer-product model with a genomic-separation term:

```bash
python 02_build_energy_matrix.py \
    --phi  data/phi_20U_rep1_chr2_4mb_25kb.npz \
    --out  data/energy_20U_rep1_chr2_4mb_25kb.npz
```

Energy terms:

| Term | Formula | Parameters |
| --- | --- | --- |
| ε_ij (condensability) | −β × (φ_i × φ_j − φ̄²) | β = 0.6, φ̄ = 0.6 |
| T_ij (separation) | ε₀ / \|j − i\|^α | ε₀ = −0.5, α = −0.5 |
| Prefactor | T_ij + ε_ij / 2 | — |

Optional flags: `--no_sep_term` (disable separation term), `--beta_scale` (scale β). Output `.npz` contains `prefactor`, `eps_ij`, `sep_term`, and run metadata.

### Step 3 — Run simulations (`03_run_simulation.py`)

Bead-spring polymer simulations were run using OpenMM (LangevinMiddleIntegrator) on NVIDIA A100 GPUs. Key physical constants and force terms:

| Physical constant | Value |
| --- | --- |
| SIGMA_NM | 10 nm (bead diameter) |
| Temperature | 300 K |
| kBT | 2.494 kJ/mol |
| Mass | 1 × 10⁶ g/mol |

| Force | Implementation | Parameters |
| --- | --- | --- |
| Bond | CustomBondForce (k2/k3/k4) | k2 = k3 = k4 = 20 kBT/σⁿ |
| Excluded volume | Softcore WCA | E_cut = 4 × kBT |
| Confinement | Harmonic sphere | κ = 1000 kBT/σ² |
| Compartment potential | Discrete Flory-Huggins contact | λ ramp 0→1 over 5,000 × 100 steps |

Contact cutoff: r_c = 1.42 × σ = 14.2 nm. Initialisation: Hilbert-like serpentine layout. Each trajectory writes `contact_traj_XX.npy`; the main function averages all trajectories to `contact_matrix_avg.npy` and `simulation_meta.npz`.

```bash
python 03_run_simulation.py \
    --energy  data/energy_20U_rep1_chr2_4mb_25kb.npz \
    --out_dir results_20U_rep1_chr2_4mb_25kb_gpu/ \
    --n_traj  16 \
    --strength 60.0 \
    --equil_steps 2000000 \
    --sim_steps   2000000 \
    --platform    CUDA
```

Optimal `--strength` was determined by parameter sweep (`sweep_strength.slurm`) over values of 5, 10, 20, 30, and 40 kBT with n_traj=4 prior to production runs.

### Step 4 — Analyse results (`04_analyse_results.py`)

```bash
python 04_analyse_results.py \
    --contact results_20U_rep1_chr2_4mb_25kb_gpu/contact_matrix_avg.npy \
    --phi     data/phi_20U_rep1_chr2_4mb_25kb.npz \
    --out_dir results_20U_rep1_chr2_4mb_25kb_gpu/
```

Processing steps:

1. O/E normalisation — each off-diagonal divided by its mean (`observed_over_expected()`)
2. A/B compartment assignment by PCA on the Pearson correlation matrix of the O/E contact map; PC1 is taken and sign-corrected so that positive values correspond to A compartment (low condensability)

Outputs per results directory:

```
results_<track>_<region>_gpu/
├── contact_map.png              # raw contact map heatmap
├── contact_map_OE.png           # O/E contact map heatmap
├── oe_matrix.npy                # O/E matrix array
├── oe_correlation_matrix.png    # Pearson correlation of O/E matrix
├── compartment_track.png        # PC1 compartment score track
├── compartment_pc1.npy          # PC1 array
├── compartment_scores.bedgraph  # compartment scores in BEDGraph format
├── condensability_binned.bedgraph
└── summary.txt                  # run parameters and convergence metrics
```

---

## SLURM Job Arrays

Production simulations were submitted as SLURM job arrays on the GPU partition (A100). Each array element ran one φ track through `03_run_simulation.py` followed immediately by `04_analyse_results.py`.

### Primary region (chr2:10–14 Mb, 25 kb bins)

```bash
sbatch run_simulations.slurm
```

Array: `0=20U_rep1`, `1=20U_rep2`, `2=mnase` — 3 jobs × 1 A100 GPU, 32 GB RAM, 4 CPUs each. strength=60.0, n_traj=16.

### Extended regions

```bash
sbatch run_30_40mb.slurm   # chr2:30–40 Mb; strength=60.0; tracks: 20U_rep1, 20U_rep2, mnase
sbatch run_35_60mb.slurm   # chr2:35–60 Mb; strength=50.0; tracks: 20U_rep1, 20U_rep2
```

### Parameter sweep (20U rep1 track only, chr2:10–14 Mb)

```bash
sbatch sweep_strength.slurm   # strength = 5, 10, 20, 30, 40 kBT; n_traj=4
```

---

## Simulation vs Hi-C Comparison

Simulated contact matrices were compared to MCF10A WT Hi-C data from GEO accession **GSE246947** (`GSE246947_hg38_MCF10A_WT_merged.mcool`, 20 kb resolution). The A/B compartment reference used for track validation and simulation calibration was `MCF10A_subcompartments_20000_T2T.bw` — a T2T/hs1-aligned subcompartment track (A0–B3) — distinct from the hg38-aligned Hi-C contact maps.

### Contact map visualization

Triangle contact maps were generated using `scipy.ndimage.rotate` (45° rotation, upper triangle). Simulation matrices at 25 kb were displayed alongside Hi-C at 20 kb on a shared genomic x-axis. Scripts generate several figure types:

- Stacked triangles: Hi-C (log scale, 95th pct) / Condense-seq supernatant / MNase sensitivity stacked vertically
- Split triangle: condense-seq supernatant (upper half) vs MNase sensitivity (lower half) on a shared matrix diagonal, in raw and O/E versions
- Difference maps: rank-normalised (simulation − Hi-C), diverging colormap

### Quantitative similarity metrics (`contact_comparison.py`, `fancy_contacts.py`)

Simulation matrices were resampled to Hi-C resolution (nearest-neighbour) for direct comparison:

- **P(s) curves** — mean contact probability vs genomic separation, linear and log-log scale
- **SCC** — stratum-adjusted correlation coefficient (HiCRep-style), computed per diagonal up to 2 Mb separation, on O/E-normalised matrices
- **Pearson r / Spearman ρ** — on O/E-normalised and rank-normalised upper triangle values

### TAD boundary analysis (`TADS_compare.py`, `TADS_cleaner.py`)

TAD boundaries were called as minima of the insulation score (Crane et al. 2015):

| Parameter | Hi-C | Simulation |
| --- | --- | --- |
| Window | 10 bins (\~200 kb at 20 kb) | 8 bins (\~200 kb at 25 kb) |
| Gaussian σ | 1.5 bins | 1.5 bins |
| Min. prominence | 0.05 | 0.05 |
| Min. distance | 3 bins | 3 bins |

Boundary overlap between Hi-C and simulation was assessed at ±200 kb tolerance; recall, precision, and F1 score were reported per condition.

---

## Track Orientation Validation

Prior to simulation, the sign convention of each input track was verified against the Hi-C compartment reference (`MCF10A_subcompartments_20000_T2T.bw`, positive = B compartment) over chr2:10–14 Mb using `compare_tracks3.py`. All tracks were used with `*_inverted.bw` files and `flip=False`:

- Condense-seq supernatant (`*_inverted.bw`): high values = B compartment enriched → positive r with B-positive compartment reference
- MNase sensitivity (`*_sensitivity.bw`): high sensitivity = A compartment → inverted for use

---

## Genome-wide Correlation Analyses

### Per-chromosome and pairwise correlations

`correlation_analysis.py` computed Pearson r between each track and the Hi-C compartment score at 1 Mb bins, genome-wide across all autosomes and chrX. A forest plot was produced showing Pearson r (filled circles) and Spearman ρ (open diamonds) with 95% confidence intervals via Fisher z-transform. Tracks were grouped by assay category (Condense-seq, MNase sensitivity, Chromatin marks).

Pre-computed hardcoded correlation values from this run were saved in `plot_correlation_forest.py` for rapid reproduction of the figure. Measured Pearson r values (genome-wide, 1 Mb bins, n ≈ 2,866):

| Track | Pearson r | Spearman ρ |
| --- | --- | --- |
| GM condense-seq (paper, GM12878) | −0.584 | −0.615 |
| 20U supernatant rep1 (MCF10A) | −0.542 | −0.625 |
| Genomic rep2 sensitivity | −0.430 | −0.458 |
| Genomic rep1 sensitivity | −0.262 | −0.307 |
| A_Inp rep1 sensitivity | −0.140 | −0.313 |

Sign convention: negative r = B-compartment enriched (as expected for condensability).

`pairwise_correlation.py` extended this to all track pairs, outputting a hierarchically clustered Pearson r heatmap and a pairwise forest plot comparing each track to the GM condense-seq and genomic sensitivity references.

`subcompartment_analysis.py` stratified 1 Mb bins by dominant subcompartment (A0–B3) and plotted violin distributions, a z-scored heatmap (tracks × subcompartments), and mean profiles across the A→B gradient. Spearman r and Kruskal-Wallis p-values vs subcompartment rank were reported per track.

### Multi-track correlation pipeline

`multi-track_correlation_clusters.py` computed genome-wide Spearman correlation matrices at 10 kb bins across 9 tracks (subcompartments, nucleosome distribution, nucleosome sensitivity, ATAC, GM NCP condensability, and 20U supernatant replicates) and produced TSS-centered profiles (±1 kb, 80 bins) stratified by promoter cluster. Clusters were derived from a deepTools k-means sort (`HP1_HPTM_2kb_sen_matrix_og_sort_k7.bed`, 7 clusters, cluster label in last column). Strand-aware averaging was applied (minus-strand TSSs reversed). Outputs include a clusters × tracks profile grid, per-track cluster overlay, and per-cluster Spearman correlation matrices. `replot_cluster_grid.py` regenerates the grid plot directly from bigWig files without recomputing intermediate matrices, using shared y-axis limits per track column.

---

## Nucleosome Sensitivity vs Compartment Analyses

### Genome-wide sensitivity stratified by subcompartment

`genomic_sensitivity_vs_compartment.py` assigned all genome-wide nucleosome calls to subcompartments A0–B3 and computed:

- Chi-squared test per subcompartment vs the genome-wide average Sensitive/Neutral/Resistant distribution
- Spearman ρ between subcompartment rank and sensitivity fraction
- Three figure types: grouped bars (counts), proportional stacked bars, fraction profiles across the A→B gradient

### TSS-level promoter classification

`nucleosome_sensitivity_vs_compartment.py` worked at the level of individual gene TSS regions. Each gene's nucleosomes were divided into upstream (first third) and downstream (last third) groups; an SR score \[(Sensitive − Resistant)/(Sensitive + Resistant)\] was computed for each group and compared to a threshold of 0 to assign one of four classes:

| Class | Upstream | Downstream |
| --- | --- | --- |
| Sens/Sens | Sensitive | Sensitive |
| Sens/Res | Sensitive | Resistant |
| Res/Res | Resistant | Resistant |
| Res/Sens | Resistant | Sensitive |

Results were plotted as stacked bar proportions per subcompartment with chi-squared significance stars.

`promoter_class_proportions.py` plots the same classification normalised to 100% within each subcompartment, comparing Scr (scramble), HP1α KD, and HP1β KD sample groups.

### Per-sample nucleosome sensitivity counts (R scripts)

Three R scripts (dplyr / ggplot2) read per-TSS sensitivity overview files, deduplicate to unique nucleosomes, and plot stacked bar figures:

| Script | Conditions | Output |
| --- | --- | --- |
| `nuc-summary_revised.R` | Scr (clone1/2), HP1α KD (clone5/6), HP1β KD (clone15/D) | Counts, proportions, and combined mean ± SD PDFs |
| `H2AZ_nucleosome_sensitivity.R` | Scr (clone1/2), H2AZ KD (rep1/rep4), HP1β KD (clone15/D) | Counts + proportions PDFs |

Input files: `TSS_10A_<condition>_<clone>_sensitivity_overview.txt` (tab-delimited, 10 columns: chr, start, end, category, gene_chr, gene_start, gene_end, gene_name, score, strand).

---

## Software Environment

| Package | Purpose |
| --- | --- |
| Python (conda env `chromatin_sim`) | Simulation and all Python analyses |
| OpenMM | GPU-accelerated molecular dynamics engine |
| NumPy | Array operations, matrix storage |
| SciPy | Peak finding, PCA, statistics, clustering |
| Matplotlib | All Python figure generation |
| pyBigWig | bigWig I/O for all Python scripts |
| cooler | Reading `.mcool` Hi-C files |
| UCSC kentUtils (`bedGraphToBigWig`) | BEDGraph → bigWig conversion |
| R (dplyr, tidyr, ggplot2) | Nucleosome sensitivity summary figures |
| CUDA (A100 GPU) | GPU acceleration via `--platform CUDA` |

---

## Reference Files

| File | Description |
| --- | --- |
| `MCF10A_subcompartments_20000_T2T.bw` | MCF10A Hi-C subcompartment calls (A0–B3) on T2T/hs1, 20 kb bins; positive = B |
| `GSE246947_hg38_MCF10A_WT_merged.mcool` | MCF10A WT Hi-C contact matrix (GEO: GSE246947, hg38, 20 kb) |
| `T2T_TSS_1kb_reduced_clean.bed` | TSS reference BED (T2T/hs1 coordinates) used for TSS profile analyses |
| `HP1_HPTM_2kb_sen_matrix_og_sort_k7.bed` | 7-cluster deepTools BED from k-means sorting on HP1/HPTM sensitivity matrix; cluster label in last column |
| `GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T.bw` | External GM12878 condense-seq condensability score (T2T coordinates) |

---

## Input / Output File Naming

```
data/
├── phi_<track>_<chrom>_<region>_<binsize>.npz         # binned φ vectors
└── energy_<track>_<chrom>_<region>_<binsize>.npz      # interaction energy matrices

results_<track>_<region>_gpu/
├── contact_traj_XX.npy                                 # per-trajectory contact matrices
├── contact_matrix_avg.npy                              # averaged contact matrix
├── oe_matrix.npy                                       # O/E normalised matrix
├── compartment_pc1.npy                                 # PC1 compartment scores
├── compartment_scores.bedgraph
├── condensability_binned.bedgraph
└── summary.txt

condensability_output/
├── <prefix>_per_nuc.bw
├── <prefix>_slide171_25.bw
├── <prefix>_1kb.bw
└── <prefix>_10kb.bw
```

Condensability tracks: `20U_rep1` / `20U_rep2` = condense-seq supernatant replicates; `mnase` = genomic MNase sensitivity.\
Regions used: `chr2_4mb_25kb` (10–14 Mb), `chr2_30_40mb_25kb` (30–40 Mb), `chr2_35_60mb_25kb` (35–60 Mb).