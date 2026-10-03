#!/bin/bash
#SBATCH --job-name=tc_spinup_cell
#SBATCH --partition=day
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/spinup_cell_%A_%a.out
#SBATCH --error=logs/spinup_cell_%A_%a.err
# One full (or RF-seeded) spin-up per SLURM array task.  Submit with
#   sbatch --array=$(cat "$TC_OUTPUT_DIR/array_indices_<all|subsample>.txt") field_spinup_cell.sbatch
# after `set -a; source field_spinup.env; set +a` -- a plain `source` does not
# export the variables, so sbatch would carry none of them to the job --
# or simply use submit_field_spinup.sh, which does all of it.
# Cost measured on Control 2 / Kantola case: ~48 s per Phase-1 year, ~6 min
# Phase 2, ~1 min per coupled (ERW-free) Phase-3 year, hours per cell -> day
# partition.  A task that hits the limit resumes its own checkpoints when the
# array is resubmitted.  MATLAB and node model as the v0.9 campaigns
# (runtime/*: MATLAB/2025a on cpumodel:9575f), so a state is produced on the
# platform the transient runs it on.
set -euo pipefail
: "${TC_ROOT:?set TC_ROOT to the repository root (the folder holding v2)}"
module load "${MATLAB_MODULE:-MATLAB/2025a}"
cd "$TC_ROOT/v2/cluster"
mkdir -p logs
matlab -batch "Cluster_Field_Spinup_Task('cell', '${SLURM_ARRAY_TASK_ID}')"
