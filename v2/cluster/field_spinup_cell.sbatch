#!/bin/bash
#SBATCH --job-name=tc_spinup_cell
#SBATCH --partition=day
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/spinup_%x_%A_%a.out
#SBATCH --error=logs/spinup_%x_%A_%a.err
# One plot-scale job per SLURM array task, TC_TASK naming which:
#   class_seed   the 1-D spin-up of class <id>     (flux_tracking)
#   replay       tracked cell <id>'s replay         (flux_tracking)
#   reconstruct  cell <id> from the forest          (flux_tracking)
#   cell         cell <id>'s own spin-up            (per_cell)
# Submit with
#   sbatch --export=ALL,TC_TASK=replay --array=$(cat "$TC_OUTPUT_DIR/array_indices_tracked.txt") field_spinup_cell.sbatch
# after `set -a; source field_spinup.env; set +a` -- a plain `source` does not
# export the variables, so sbatch would carry none of them to the job --
# or simply use submit_field_spinup.sh, which does all of it.
# Cost on one Control 2 cell (Kantola case, laptop): ~48 s per Phase-1
# year, a few minutes for Phase 2, ~1-2 min per coupled ERW-free Phase-3
# year over many rings -> hours per class seed or replay; a reconstruction
# is minutes (its ERW handoff only).  A task that hits the limit resumes its
# own checkpoints when the array is resubmitted.  MATLAB and node model as
# the v0.9 campaigns (runtime/*: MATLAB/2025a on cpumodel:9575f), so a state
# is produced on the platform the transient runs it on.
set -euo pipefail
: "${TC_ROOT:?set TC_ROOT to the repository root (the folder holding v2)}"
: "${TC_TASK:?set TC_TASK to class_seed | replay | reconstruct | cell}"
module load "${MATLAB_MODULE:-MATLAB/2025a}"
cd "$TC_ROOT/v2/cluster"
mkdir -p logs
matlab -batch "Cluster_Field_Spinup_Task('${TC_TASK}', '${SLURM_ARRAY_TASK_ID}')"
