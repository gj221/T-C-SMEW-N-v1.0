#!/bin/bash
#SBATCH --job-name=tc_spinup_track
#SBATCH --partition=day
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=logs/spinup_track_%j.out
#SBATCH --error=logs/spinup_track_%j.err
#SBATCH --requeue
#SBATCH --signal=B:USR1@900
#SBATCH --open-mode=append
# The 2-D tracking run of the flux-tracking spin-up (Lian et al. 2025,
# Fig. 3c): the whole field from its class seeds, ERW off, one segment per
# spin-up cycle year and then the handoff extension.  Every column steps
# every hour, so a segment costs about (cells x one ERW-free column-year):
# ~1 h for the 54-cell diagnostic grid, about a day for the 1461-cell
# production grid, where memory also grows with the cells (raise --mem).
# Finished segments are reused and the running one resumes from the 2-D
# driver's hourly rolling checkpoint, so near its limit the job requeues
# itself (requeue_on_timeout.sh) and carries on.  submit_field_spinup.sh
# also queues TC_TRACK_JOBS copies, each starting after the previous one
# ends -- the fallback where requeueing is refused -- and a copy that finds
# the run complete exits at once.
set -euo pipefail
: "${TC_ROOT:?set TC_ROOT to the repository root (the folder holding v2)}"
module load "${MATLAB_MODULE:-MATLAB/2025a}"
cd "$TC_ROOT/v2/cluster"
mkdir -p logs
source ./requeue_on_timeout.sh     # near the time limit the job requeues itself
run_requeueable matlab -batch "Cluster_Field_Spinup_Task('track')"
