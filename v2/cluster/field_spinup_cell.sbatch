#!/bin/bash
#SBATCH --job-name=tc_spinup_cell
#SBATCH --partition=week
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/spinup_cell_%A_%a.out
#SBATCH --error=logs/spinup_cell_%A_%a.err
# One full (or RF-seeded) spin-up per SLURM array task.  Submit with
#   sbatch --array=$(cat "$TC_OUTPUT_DIR/array_indices_<all|subsample>.txt") field_spinup_cell.sbatch
# after `source field_spinup.env` (or submit_field_spinup.sh, which does all of it).
# Cost measured on Control 2 / Kantola case: ~48 s per Phase-1 year, ~3 min
# Phase 2, ~70 min per coupled Phase-3 year, 10-20 h per cell -> week partition.
set -euo pipefail
: "${TC_ROOT:?set TC_ROOT to the repository root (the folder holding v2)}"
module load MATLAB/2023b
cd "$TC_ROOT/v2/cluster"
mkdir -p logs
matlab -batch "Cluster_Field_Spinup_Task('cell', '${SLURM_ARRAY_TASK_ID}')"
