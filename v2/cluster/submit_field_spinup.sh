#!/bin/bash
# Submit the whole field spin-up as dependent SLURM jobs on Bouchet.
#
#   cd <repo>/v2/cluster
#   cp field_spinup.env.example field_spinup.env   # edit paths and settings
#   bash submit_field_spinup.sh [field_spinup.env]
#
# Chain (Method per_cell):  plan -> array(all cells) -> assemble -> validate
# Chain (Method rf):        plan -> array(subsample, full spin-ups)
#                                -> array(remainder, RF-seeded) -> assemble
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
cd "$TC_ROOT/v2/cluster"
mkdir -p logs "$TC_OUTPUT_DIR"

echo "plan ..."
sbatch --wait --export=ALL,TC_TASK=plan field_spinup_step.sbatch
for f in all subsample rf; do
  [[ -f "$TC_OUTPUT_DIR/array_indices_$f.txt" ]] || { echo "plan did not write array_indices_$f.txt" >&2; exit 1; }
done

# an index list is a SLURM --array range ("1-54") or comma list; strip whitespace
read_list() { local v; v=$(cat "$1"); echo "${v//[[:space:]]/}"; }
if [[ "${TC_METHOD:-per_cell}" == "rf" ]]; then
  sub=$(read_list "$TC_OUTPUT_DIR/array_indices_subsample.txt")
  rest=$(read_list "$TC_OUTPUT_DIR/array_indices_rf.txt")
  [[ -n "$sub" ]] || { echo "the plan has no subsample cells" >&2; exit 1; }
  echo "subsample full spin-ups: --array=$sub"
  a1=$(sbatch --parsable --export=ALL --array="$sub" field_spinup_cell.sbatch)
  if [[ -n "$rest" ]]; then
    echo "RF-seeded cells: --array=$rest (after $a1)"
    a2=$(sbatch --parsable --export=ALL --array="$rest" --dependency=afterok:"$a1" field_spinup_cell.sbatch)
    dep="$a2"
  else
    echo "every cell is in the subsample: no RF-seeded array"
    dep="$a1"
  fi
else
  all=$(read_list "$TC_OUTPUT_DIR/array_indices_all.txt")
  [[ -n "$all" ]] || { echo "the plan has no cells" >&2; exit 1; }
  echo "full spin-ups for every cell: --array=$all"
  dep=$(sbatch --parsable --export=ALL --array="$all" field_spinup_cell.sbatch)
fi
asm=$(sbatch --parsable --export=ALL,TC_TASK=assemble --dependency=afterok:"$dep" field_spinup_step.sbatch)
echo "assemble job $asm (after $dep)"
if [[ "${TC_METHOD:-per_cell}" != "rf" && "${TC_VALIDATE:-1}" == "1" ]]; then
  val=$(sbatch --parsable --export=ALL,TC_TASK=validate --dependency=afterok:"$asm" field_spinup_step.sbatch)
  echo "validation job $val (after $asm)"
fi
echo "monitor with: squeue -u \$USER ; logs in $TC_ROOT/v2/cluster/logs ; outputs in $TC_OUTPUT_DIR"
