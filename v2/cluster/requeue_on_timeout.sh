# Sourced by the field spin-up batch scripts.  run_requeueable runs one
# step (the MATLAB call) so that a time limit does not end the chain: the
# script asks SLURM for a warning signal before the limit
# (#SBATCH --signal=B:USR1@900, with --requeue), and on it the job requeues
# itself under its own id -- an array task as <array id>_<task id> -- so
# every job that depends on it keeps waiting, and the next run resumes the
# step's own checkpoints (a plot-scale spin-up each Phase-3 year, the 2-D
# tracking run its hourly rolling checkpoint).  A step that fails is not
# requeued: only the time limit is.  At most TC_MAX_REQUEUES requeues
# (default 10) per job, so a step that cannot finish stops rather than
# cycling; the log says each time.
run_requeueable() {
  local job="${SLURM_JOB_ID:-}"
  if [[ -n "${SLURM_ARRAY_JOB_ID:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    job="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
  fi
  local restarts="${SLURM_RESTART_COUNT:-0}" cap="${TC_MAX_REQUEUES:-10}"
  echo "$(date '+%F %T') ${job:-local}: start (restart ${restarts} of at most ${cap})"
  requeue_near_limit() {
    if [[ -z "$job" ]]; then
      echo "$(date '+%F %T') time limit near, but not a SLURM job: nothing to requeue" >&2
    elif (( restarts < cap )); then
      echo "$(date '+%F %T') ${job}: time limit near, requeueing (restart $((restarts + 1)) of ${cap})"
      scontrol requeue "$job" || echo "$(date '+%F %T') ${job}: requeue refused; resubmit this step by hand" >&2
    else
      echo "$(date '+%F %T') ${job}: time limit near and ${cap} requeues used; not requeueing" >&2
    fi
  }
  trap requeue_near_limit USR1
  "$@" &
  local pid=$! rc=0
  # The trap interrupts wait; the step itself runs on until it ends or SLURM
  # stops the requeued job, so wait again while it is alive.
  while :; do
    rc=0
    wait "$pid" || rc=$?
    kill -0 "$pid" 2>/dev/null || break
  done
  trap - USR1
  return "$rc"
}
