#!/bin/bash
#SBATCH --job-name=chromatin_sweep
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-4
#SBATCH --output=logs/sweep_%A_%a.out
#SBATCH --error=logs/sweep_%A_%a.err
#SBATCH --time=48:00:00


# =============================================================================
# Parameter sweep: strength = 5, 10, 20, 30, 40 kBT
# Run supernatant rep1 track at each strength — find optimal before production
# =============================================================================

PYTHON=/home/gdbishoplab/jxb377/.conda/envs/chromatin_sim/bin/python

WORKDIR=/home/gdbishoplab/jxb377/condense-seq
cd $WORKDIR
mkdir -p logs

STRENGTHS=(5.0 10.0 20.0 30.0 40.0)
STRENGTH=${STRENGTHS[$SLURM_ARRAY_TASK_ID]}
LABEL=${STRENGTH/./_}
OUTDIR="results_sweep_s${LABEL}/"

echo "Strength: $STRENGTH kBT -> $OUTDIR"
echo "Start: $(date)"
$PYTHON --version

$PYTHON 03_run_simulation.py \
    --energy      data/energy_supernatant_rep1_chr2_4mb_25kb.npz \
    --out_dir     $OUTDIR \
    --n_traj      4 \
    --equil_steps 1000000 \
    --sim_steps   2000000 \
    --save_every  2000 \
    --strength    $STRENGTH \
    --conc        2.0 \
    --platform    CUDA

$PYTHON 04_analyse_results.py \
    --contact $OUTDIR/contact_matrix_avg.npy \
    --phi     data/phi_supernatant_rep1_chr2_4mb_25kb.npz \
    --out_dir $OUTDIR/

echo "End: $(date)"
cat $OUTDIR/summary.txt
