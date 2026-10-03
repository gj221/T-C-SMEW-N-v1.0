#!/bin/bash
# Submit the whole field spin-up as dependent SLURM jobs on Bouchet.
#
#   cd <repo>/v2/cluster
#   cp field_spinup.env.example field_spinup.env   # edit paths and settings
#   bash submit_field_spinup.sh [field_spinup.env]
#
# Chain (Method flux_tracking, Lian et al. 2025):
#   plan -> array(class seeds) -> track (TC_TRACK_JOBS chained copies)
#        -> array(tracked cells: replays) -> array(other cells: forest)
#        -> assemble
# Chain (Method per_cell):  plan -> array(all cells) -> assemble -> validate
# The plan job runs first and blocks (sbatch --wait, a few minutes) because
# the array indices only exist once it has written them.
set -euo pipefail
ENV_FILE="${1:-field_spinup.env}"
if [[ ! -f "$ENV_FILE" ]]; then
  echo "environment file $ENV_FILE not found (copy field_spinup.env.example)" >&2
  exit 1
fi
set -a; source "$ENV_FILE"; set +a
: "${TC_ROOT:?}" "${TC_SHAPEFILE:?}" "${TC_AGRONOMIC_FILE:?}" "${TC_OUTPUT_DIR:?}"
METHOD="${TC_METHOD:-flux_tracking}"
cd "$TC_ROOT/v2/cluster"
mkdir -p logs "$TC_OUTPUT_DIR"

echo "plan ..."
sbatch --wait --export=ALL,TC_TASK=plan field_spinup_step.sbatch

# an index list is a SLURM --array range ("1-54") or comma list; strip whitespace
read_list() {
  [[ -f "$1" ]] || { echo "the plan did not write $1" >&2; exit 1; }
  local v; v=$(cat "$1"); echo "${v//[[:space:]]/}"
}
if [[ "$METHOD" == "flux_tracking" ]]; then
  classes=$(read_list "$TC_OUTPUT_DIR/array_indices_classes.txt")
  tracked=$(read_list "$TC_OUTPUT_DIR/array_indices_tracked.txt")
  rest=$(read_list "$TC_OUTPUT_DIR/array_indices_rf.txt")
  [[ -n "$classes" && -n "$tracked" ]] || { echo "the plan has no classes or no tracked cells" >&2; exit 1; }
  echo "class seeds: --array=$classes"
  dep=$(sbatch --parsable --export=ALL,TC_TASK=class_seed --array="$classes" field_spinup_cell.sbatch)
  n_track="${TC_TRACK_JOBS:-3}"
  for ((k = 1; k <= n_track; k++)); do
    # every copy after the first continues the previous one's work
    if (( k == 1 )); then cond="afterok:$dep"; else cond="afterany:$trk"; fi
    trk=$(sbatch --parsable --export=ALL --dependency="$cond" field_spinup_track.sbatch)
    echo "tracking run, job $k of $n_track: $trk ($cond)"
  done
  rep=$(sbatch --parsable --export=ALL,TC_TASK=replay --array="$tracked" \
        --dependency=afterok:"$trk" field_spinup_cell.sbatch)
  echo "replays of the tracked cells: --array=$tracked, job $rep"
  dep="$rep"
  if [[ -n "$rest" ]]; then
    rec=$(sbatch --parsable --export=ALL,TC_TASK=reconstruct --array="$rest" \
          --dependency=afterok:"$rep" field_spinup_cell.sbatch)
    echo "forest reconstruction of the other cells: --array=$rest, job $rec"
    dep="$rec"
  fi
else
  all=$(read_list "$TC_OUTPUT_DIR/array_indices_all.txt")
  [[ -n "$all" ]] || { echo "the plan has no cells" >&2; exit 1; }
  echo "full spin-ups for every cell: --array=$all"
  dep=$(sbatch --parsable --export=ALL,TC_TASK=cell --array="$all" field_spinup_cell.sbatch)
fi
asm=$(sbatch --parsable --export=ALL,TC_TASK=assemble --dependency=afterok:"$dep" field_spinup_step.sbatch)
echo "assemble job $asm (after $dep)"
if [[ "${TC_VALIDATE:-0}" == "1" ]]; then
  val=$(sbatch --parsable --export=ALL,TC_TASK=validate --dependency=afterok:"$asm" field_spinup_step.sbatch)
  echo "validation job $val (after $asm)"
fi
echo "monitor with: squeue -u \$USER ; logs in $TC_ROOT/v2/cluster/logs ; outputs in $TC_OUTPUT_DIR"
