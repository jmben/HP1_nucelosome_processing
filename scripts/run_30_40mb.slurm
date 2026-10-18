#!/bin/bash
#SBATCH --job-name=sim_3040mb
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-2
#SBATCH --output=logs/sim_3040mb_%A_%a.out
#SBATCH --error=logs/sim_3040mb_%A_%a.err

PYTHON=/home/gdbishoplab/jxb377/.conda/envs/chromatin_sim/bin/python
WORKDIR=/home/gdbishoplab/jxb377/condense-seq
cd $WORKDIR

TRACKS=(supernatant_rep1 rep2 rep1)
TRACK=${TRACKS[$SLURM_ARRAY_TASK_ID]}

echo "Track: $TRACK  Node: $SLURMD_NODENAME  Start: $(date)"

$PYTHON 03_run_simulation.py \
    --energy data/energy_${TRACK}_chr2_30_40mb_25kb.npz \
    --out_dir results_${TRACK}_chr2_30_40mb_gpu/ \
    --n_traj 16 \
    --equil_steps 2000000 \
    --sim_steps   2000000 \
    --save_every  2000 \
    --strength    60.0 \
    --conc        2.0 \
    --platform    CUDA

$PYTHON 04_analyse_results.py \
    --contact results_${TRACK}_chr2_30_40mb_gpu/contact_matrix_avg.npy \
    --phi     data/phi_${TRACK}_chr2_30_40mb_25kb.npz \
    --out_dir results_${TRACK}_chr2_30_40mb_gpu/

echo "End: $(date)"
cat results_${TRACK}_chr2_30_40mb_gpu/summary.txt
