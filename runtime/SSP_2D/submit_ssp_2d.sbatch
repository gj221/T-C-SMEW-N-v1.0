#!/usr/bin/env bash
#SBATCH --job-name=SSP2D
#SBATCH --partition=week
#SBATCH --constraint=cpumodel:9575f
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=384G
#SBATCH --time=4-00:00:00
#SBATCH --requeue
#SBATCH --output=runtime/SSP_2D/slurm_%A_%a.out
#SBATCH --error=runtime/SSP_2D/slurm_%A_%a.err
#
# Runs 1-9 of one SSP campaign. TC_SSP and TC_HORIZON select the scenario and
# the run length; both enter the campaign tag, so this one file serves every
# scenario and both horizons and each campaign writes to its own tree.
#
# The directives above size the 35-year case. The controller overrides
# --partition, --time and --mem for 5yr, so those values are the ceiling
# rather than the only setting.
#
# Memory: the completed five-year 2-D runs peaked at 17.3 GB (control) and
# 24.2 GB (basalt) against a 48 G request. State arrays are indexed by day
# and hour, so seven times the span puts a 35-year run near 120-170 GB. The
# 384 G here is about 2.3x that, because the extrapolation is linear while
# memory need not be, and because only two of the nine 2-D treatments have
# ever been measured. At the old 48 G these jobs would be killed part-way
# through after many hours, which costs far more than the headroom does:
# 384 G is 17% of a 2251 GB week node. Confirm against the
# Worker_Scaling_Test MaxRSS figures before trusting the extrapolation.
#
# Twelve workers, one per column. See the worker-count comment in the task.

set -euo pipefail
umask 027

repo_root="${SLURM_SUBMIT_DIR:?Submit this job from the repository root}"
if [[ ! -f "$repo_root/v2/run/Run_One_Subfield.m" ]]; then
    echo "Submit from the T&C-SMEW-N v0.9 checkout root." >&2
    exit 2
fi
cd "$repo_root"

# Only the campaign controller may submit this file: it runs the preflight,
# pins the checkout and builds the dependency graph.
if [[ "${TC_SSP_CONTROLLER_TOKEN:-}" != "SSP_CAMPAIGN_V1" ]]; then
    echo "Submit through runtime/SSP_1D/submit_ssp.sh, not sbatch." >&2
    exit 7
fi
# The equivalence baseline is a fingerprinted dependency shared with earlier
# campaigns; rebuilding it would change every stamp.
if [[ "${TC_REBUILD_EQUIVALENCE:-0}" != "0" ]]; then
    echo "This campaign never rebuilds the equivalence baseline." >&2
    exit 7
fi

: "${TC_SSP:?TC_SSP must be ssp126, ssp245 or ssp585}"
: "${TC_HORIZON:?TC_HORIZON must be 5yr or 35yr}"
case "$TC_SSP" in ssp126|ssp245|ssp585) ;; *) echo "Bad TC_SSP: $TC_SSP" >&2; exit 2 ;; esac
case "$TC_HORIZON" in 5yr|35yr) ;; *) echo "Bad TC_HORIZON: $TC_HORIZON" >&2; exit 2 ;; esac

module purge
module load "${MATLAB_MODULE:-MATLAB/2025a}"

node_record="$(scontrol show node "$SLURMD_NODENAME")"
printf '%s\n' "$node_record"
if [[ "$node_record" != *"cpumodel:9575f"* ]]; then
    echo "The job is not running on the required cpumodel:9575f node." >&2
    exit 3
fi
hostname
lscpu
module list

export TC_PARALLEL_WORKERS=12
case_suffix=''
if [[ -n "${TC_CALCITE_CASE:-}" ]]; then
    case_suffix="_case${TC_CALCITE_CASE}"
fi
export TC_CAMPAIGN_TAG="SSP_2D_${TC_HORIZON}_${TC_SSP}${case_suffix}"
export TC_NUMERICAL_IDENTITY='cpumodel=9575f;matlab_release=2025a;single_comp_thread=1;parallel_workers=12'
export TC_CKPT_MINUTES=30
unset TC_CKPT_FORCE_RESUME || true
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

scratch_parent="${TMPDIR:-/tmp/$USER}"
mkdir -p "$scratch_parent"
scratch_parent="$(cd "$scratch_parent" && pwd -P)"
scratch_root="$scratch_parent/tc_smewn_ssp2d_${SLURM_JOB_ID}_${SLURM_RESTART_COUNT:-0}"
mkdir -p "$scratch_root/prefs" "$scratch_root/tmp"
cleanup_scratch() {
    case "$scratch_root" in
        "$scratch_parent"/tc_smewn_ssp2d_*) rm -rf -- "$scratch_root" ;;
        *) echo "Refusing to remove unexpected scratch path: $scratch_root" >&2 ;;
    esac
}
trap cleanup_scratch EXIT

export MATLAB_PREFDIR="$scratch_root/prefs"
export TMPDIR="$scratch_root/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

campaign='SSP_2D'
srun --cpu-bind=cores matlab -singleCompThread -batch \
    "addpath(fullfile(pwd,'runtime','$campaign')); Run_SSP_2D_Task"

echo "SSP 2-D ${TC_SSP} ${TC_HORIZON}${case_suffix} array element ${SLURM_ARRAY_TASK_ID} finished."
