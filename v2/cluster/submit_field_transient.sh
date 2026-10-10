#!/bin/bash
# Submit the field transient (Run_Field_Transient) on Bouchet, from the
# states the field spin-up assembled.
#
#   cd <repo>/v2/cluster
#   bash submit_field_transient.sh field_spinup_hamden.env            # states on disk
#   bash submit_field_transient.sh field_spinup_hamden.env <job id>   # after a job,
#                                   e.g. the spin-up's assemble job
#
# The same environment file as the spin-up: the grid, the column files and
# the DEM must be the ones the states were spun up with. The transient reads
# TC_STATES_FILE (default $TC_OUTPUT_DIR/CellBGStates.mat) and writes to
# TC_TRANSIENT_OUTPUT_DIR (default $TC_OUTPUT_DIR/transient). It queues
# TC_TRANSIENT_JOBS copies (default 2), each starting after the previous one
# ends; a copy continues where the last stopped and one that finds the run
# complete exits at once. Each job also requeues itself near its time limit,
# up to TC_TRANSIENT_MAX_REQUEUES times (default 40).
# TC_TRANSIENT_CEC="6.4 10.8 13.7" submits one such chain per CEC scenario
# (every cell at that CEC), each into <output dir>_CEC_<value>; unset, the
# cells keep their own CEC.
# TC_TRANSIENT_PARTITION, TC_TRANSIENT_TIME, TC_TRANSIENT_MEM and
# TC_TRANSIENT_CPUS override the batch script's week / 4-00:00:00 / 64G / 32.
set -euo pipefail
ENV_FILE="${1:-field_spinup.env}"
AFTER="${2:-}"
if [[ ! -f "$ENV_FILE" ]]; then
  echo "environment file $ENV_FILE not found (copy field_spinup_hamden.env.example)" >&2
  exit 1
fi
set -a; source "$ENV_FILE"; set +a
: "${TC_ROOT:?}" "${TC_SHAPEFILE:?}" "${TC_AGRONOMIC_FILE:?}" "${TC_OUTPUT_DIR:?}"
STATES="${TC_STATES_FILE:-$TC_OUTPUT_DIR/CellBGStates.mat}"
if [[ -z "$AFTER" && ! -f "$STATES" ]]; then
  echo "no spun-up states at $STATES: run the spin-up first, or pass the job id it ends with" >&2
  exit 1
fi
export TC_MAX_REQUEUES="${TC_TRANSIENT_MAX_REQUEUES:-40}"
cd "$TC_ROOT/v2/cluster"
mkdir -p logs
OUT="${TC_TRANSIENT_OUTPUT_DIR:-$TC_OUTPUT_DIR/transient}"

resources=()
[[ -n "${TC_TRANSIENT_PARTITION:-}" ]] && resources+=(--partition="$TC_TRANSIENT_PARTITION")
[[ -n "${TC_TRANSIENT_TIME:-}" ]] && resources+=(--time="$TC_TRANSIENT_TIME")
[[ -n "${TC_TRANSIENT_MEM:-}" ]] && resources+=(--mem="$TC_TRANSIENT_MEM")
[[ -n "${TC_TRANSIENT_CPUS:-}" ]] && resources+=(--cpus-per-task="$TC_TRANSIENT_CPUS")
n_jobs="${TC_TRANSIENT_JOBS:-2}"
read -r -a scenarios <<< "${TC_TRANSIENT_CEC:-}"
(( ${#scenarios[@]} > 0 )) || scenarios=("")
for cec in "${scenarios[@]}"; do
  export TC_TRANSIENT_CEC="$cec"
  out="$OUT"
  [[ -n "$cec" ]] && out="${OUT}_CEC_${cec//./p}"
  mkdir -p "$out"
  job=""
  for ((k = 1; k <= n_jobs; k++)); do
    dependency=()
    if (( k == 1 )); then
      [[ -n "$AFTER" ]] && dependency=(--dependency="afterok:$AFTER")
    else
      # every copy after the first continues the previous one's work
      dependency=(--dependency="afterany:$job")
    fi
    job=$(sbatch --parsable --export=ALL ${resources[@]+"${resources[@]}"} \
          ${dependency[@]+"${dependency[@]}"} field_transient.sbatch)
    echo "field transient${cec:+ (CEC $cec)}, job $k of $n_jobs: $job ${dependency[*]:-}"
  done
  echo "  output: $out"
done
echo "states: $STATES"
echo "monitor with: squeue -u \$USER ; logs in $TC_ROOT/v2/cluster/logs"
