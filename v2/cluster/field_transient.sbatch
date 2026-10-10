#!/bin/bash
#SBATCH --job-name=tc_transient
#SBATCH --partition=week
#SBATCH --time=4-00:00:00
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/transient_%j.out
#SBATCH --error=logs/transient_%j.err
#SBATCH --requeue
#SBATCH --signal=B:USR1@900
#SBATCH --open-mode=append
# The field transient (Run_Field_Transient): every cell from the assembled
# spin-up states (CellBGStates.mat) over the transient window, in 2-D with
# the ERW chemistry on, as segments of TC_SEGMENT_DAYS that each go on from
# the one before. A finished segment is reused and the running one resumes
# from the driver's rolling checkpoint (TC_TRANSIENT_CKPT_MINUTES), so near
# its limit the job requeues itself (requeue_on_timeout.sh) and carries on.
# MATLAB runs single-threaded, as the v0.9 2-D campaigns do; the daily ERW
# chemistry of the columns runs on TC_PARALLEL_WORKERS thread workers (one
# per CPU). Submit with submit_field_transient.sh, which exports the
# environment and chains TC_TRANSIENT_JOBS copies.
set -euo pipefail
: "${TC_ROOT:?set TC_ROOT to the repository root (the folder holding v2)}"
module load "${MATLAB_MODULE:-MATLAB/2025a}"
export TC_PARALLEL_WORKERS="${TC_PARALLEL_WORKERS:-${SLURM_CPUS_PER_TASK:-1}}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$TC_ROOT/v2/cluster"
mkdir -p logs
source ./requeue_on_timeout.sh     # near the time limit the job requeues itself
run_requeueable matlab -singleCompThread -batch "Cluster_Field_Spinup_Task('transient')"
