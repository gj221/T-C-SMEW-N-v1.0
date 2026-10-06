#!/bin/bash
#SBATCH --job-name=tc_spinup_step
#SBATCH --partition=day
#SBATCH --time=06:00:00
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/spinup_%x_%j.out
#SBATCH --error=logs/spinup_%x_%j.err
#SBATCH --requeue
#SBATCH --signal=B:USR1@900
#SBATCH --open-mode=append
# The serial steps around the cell array: plan | assemble | validate.
# Submitted by submit_field_spinup.sh; by hand, first export the
# environment with `set -a; source field_spinup.env; set +a`, then:
#   sbatch --export=ALL,TC_TASK=plan     field_spinup_step.sbatch
#   sbatch --export=ALL,TC_TASK=assemble --dependency=afterok:<array job id> field_spinup_step.sbatch
#   sbatch --export=ALL,TC_TASK=validate --dependency=afterok:<assemble job id> field_spinup_step.sbatch
set -euo pipefail
: "${TC_ROOT:?set TC_ROOT to the repository root (the folder holding v2)}"
: "${TC_TASK:?set TC_TASK to plan | assemble | validate}"
module load "${MATLAB_MODULE:-MATLAB/2025a}"
cd "$TC_ROOT/v2/cluster"
mkdir -p logs
source ./requeue_on_timeout.sh     # near the time limit the job requeues itself
run_requeueable matlab -batch "Cluster_Field_Spinup_Task('${TC_TASK}')"
