#!/bin/bash
#SBATCH --job-name=chromatin_sim
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-2
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jxb377@case.edu

# =============================================================================
# Chromatin polymer simulation — SLURM job array
# Array 0=supernatant_rep1, 1=rep2, 2=rep1
# Each job: 16 trajectories x 2M steps on A100 GPU
# =============================================================================

# Use explicit python path to avoid conda activation issues
PYTHON=/home/gdbishoplab/jxb377/.conda/envs/chromatin_sim/bin/python

WORKDIR=/home/gdbishoplab/jxb377/condense-seq
cd $WORKDIR
mkdir -p logs

TRACKS=(supernatant_rep1 rep2 rep1)
TRACK=${TRACKS[$SLURM_ARRAY_TASK_ID]}

ENERGY="data/energy_${TRACK}_chr2_4mb_25kb.npz"
PHI="data/phi_${TRACK}_chr2_4mb_25kb.npz"
OUTDIR="results_${TRACK}_chr2_25kb_gpu/"

echo "Track: $TRACK  Node: $SLURMD_NODENAME  GPU: $CUDA_VISIBLE_DEVICES"
echo "Start: $(date)"
$PYTHON --version

$PYTHON 03_run_simulation.py \
    --energy      $ENERGY \
    --out_dir     $OUTDIR \
    --n_traj      16 \
    --equil_steps 2000000 \
    --sim_steps   2000000 \
    --save_every  2000 \
    --strength    60.0 \
    --conc        2.0 \
    --platform    CUDA

$PYTHON 04_analyse_results.py \
    --contact $OUTDIR/contact_matrix_avg.npy \
    --phi     $PHI \
    --out_dir $OUTDIR/

echo "End: $(date)"
cat $OUTDIR/summary.txt
