# Field spin-up and transient on Bouchet

The field spin-up is the flux-tracking spin-up of T&C-BG-2D (Lian et al.
2025, *WRR*, section 2.3), adapted for ERW; the main README, *Flux-tracking
field spin-up*, describes the method. The field transient that follows it
is below (*The field transient*). The spin-up has four kinds of work:

- **class seeds** and **replays** are plot-scale spin-ups: about 48 s per
  forcing year in the T&C-only Phase 1, a few minutes in the stand-alone
  Phase 2 and about 1-2 min per coupled, ERW-free Phase-3 year, whose
  convergence test runs over complete rings of cycle years (one Control 2
  cell, this workstation). One is a job of many hours: the test cell's
  Phase 3 converged after 306 coupled years (about 5.4 h), its last drift
  held up by the slow inorganic phosphorus pools, so a task may need
  resubmitting to finish; it resumes from its checkpoint. A warm start
  (below) runs Phase 3 alone from a converged state: the same cell
  re-converged in 24 coupled years (17 min) after its pH, CEC and base
  saturation were revised.
- the **tracking run** is the whole field in 2-D for the ring and its
  handoff extension (ERW off): about one ERW-free column-year per cell
  per segment, so ~1 h per segment on the 54-cell diagnostic grid and
  about a day per segment on the 1461-cell production grid.
- a **reconstruction** is the forest's prediction plus the cell's ERW
  handoff gate: minutes.

Everything is resumable: finished class seeds, segments, replays and
reconstructions are reused under their provenance keys, and a running
plot-scale spin-up or tracking segment resumes its own checkpoints.

**Time limits are handled automatically.** Every job asks SLURM for a
warning 15 minutes before its limit (`--signal=B:USR1@900`, `--requeue`)
and, on it, requeues itself under its own job id (an array task as
`<array id>_<task id>`; `requeue_on_timeout.sh`), so every job that
depends on it keeps waiting and the requeued run resumes the step's
checkpoints: a plot-scale spin-up loses at most one Phase-3 year, the
tracking run at most an hour (its rolling checkpoint). A step that fails
is not requeued, only one that runs out of time, and at most
`TC_MAX_REQUEUES` times per job (default 10); the log records each
restart. Logs are appended across restarts. Should the cluster refuse a
requeue the log says so, the job ends at its limit, and that step is
resubmitted by hand as before (the chained tracking copies,
`TC_TRACK_JOBS`, remain the tracking run's fallback).

## What runs

| step | script | what it does |
|---|---|---|
| plan | `field_spinup_step.sbatch` (`TC_TASK=plan`) | `spinup_plan.csv`, the classes and the tracked subsample, and the array-index lists |
| class seeds | `field_spinup_cell.sbatch` (`TC_TASK=class_seed`, array over classes) | one plot-scale spin-up per land cover x soil class, on its most typical cell, with its ERW seed's inputs beside it (`classseed_<k>_erw_seed_inputs.mat`) |
| tracking | `field_spinup_track.sbatch` | the 2-D field run from the class seeds, recording the tracked cells' lateral exchange (`track/record_<id>.mat`, `track/ringend_<id>.mat`, `track/handoff_states.mat`) |
| replays | `field_spinup_cell.sbatch` (`TC_TASK=replay`, array over tracked cells) | each tracked cell's spin-up replaying its record |
| reconstruction | `field_spinup_cell.sbatch` (`TC_TASK=reconstruct`, array over the other cells) | the forest over the replayed cells, then the cell's own ERW handoff |
| assemble | `field_spinup_step.sbatch` (`TC_TASK=assemble`) | gathers `cellstate_<id>.mat` into the validated `CellBGStates.mat` the transient run reads |
| validate | `field_spinup_step.sbatch` (`TC_TASK=validate`) | `Run_RF_Validation_Experiment` on an assembled truth |
| transient | `field_transient.sbatch` (`submit_field_transient.sh`) | the field from `CellBGStates.mat` over the transient window, ERW on, in segments (below) |

`submit_field_spinup.sh` chains them with SLURM dependencies (with
`TC_METHOD=per_cell`: plan, one array of every cell's own spin-up,
assemble). All settings are environment variables in `field_spinup.env`
(copy and edit `field_spinup.env.example`); `Cluster_Field_Spinup_Task.m`
is the only MATLAB entry point. The tracking job may need more than one
day on a large grid: the submitter queues `TC_TRACK_JOBS` copies (default
3), each starting after the previous one ends; a copy continues where the
last stopped, and one that finds the run complete exits at once.

`TC_STOP_AFTER=class_seed` ends the chain after the class seeds. A site
whose exchanger emerges from the spin-up (`spinup_equilibrium`) refits its
exchange constants there on a class seed's converged soil solution
(`classseed_<k>_erw_seed_inputs.mat`, `v2/calibration/Fit_Exchange_Affinity.m`)
and submits again without the variable; the finished class seeds are
reused. Every replay and reconstruction hands off with the constants in
force when it runs, and a finished one is reused as it is, so constants
changed after them need a new `TC_OUTPUT_DIR`.

## Getting the repository and inputs onto the cluster

There is no SSH key on the Windows workstation, so use the web portal
or Globus:

1. Open OnDemand: https://ood-bouchet.ycrc.yale.edu -> Files. Upload a
   zip of the repository (exclude `v2/SIM`, `v2/logs`, `.git`) into
   `~/project_pi_<netid>/TC-SMEW-N` and unzip it there (Files -> Open in
   Terminal, `unzip`). For the DEM tiles (`v2/Spatial inputs/Bruce Farm
   Virginia/DEM`, gitignored) and the forcing `.mat` files in
   `v2/T&C_Code/Inputs`, Globus is faster for anything above a few
   hundred MB.
2. Put results on scratch (`~/scratch_pi_<netid>`, purged after 30
   days) and copy `CellBGStates.mat` and the CSVs back to project when
   a run finishes.
3. The per-cell agronomic CSV is produced locally (templates, kriging,
   gated merge; see the main README, *Parameterizing grid cells*) and
   uploaded with the repository. The Control 2 diagnostic file is
   `v2/T&C_Code/Inputs/Field_Spinup/diagnostic_Control2_15m/Cell_Agronomic_merged.csv`.

## Running

```bash
cd ~/project_pi_<netid>/TC-SMEW-N/v2/cluster
cp field_spinup.env.example field_spinup.env    # edit paths, method, grid
bash submit_field_spinup.sh
squeue -u $USER                                  # progress
tail logs/spinup_tc_spinup_cell_<jobid>_<task>.out   # one task's log
```

The plot-scale arrays use the `day` partition with a one-day limit per
task and 8 GB; the DEM loader crops each 1 m lidar tile to the field
window as it reads it, so a task holds one tile (about 0.8 GB) at its
peak. The tracking job holds every column's hourly state for one
segment; 64 GB is ample for the diagnostic grid, and the production grid
needs a large-memory node (raise `--mem` in `field_spinup_track.sbatch`).
If a task is preempted, or ends at its limit with no requeue left,
resubmit the same step: finished work is recognised by a provenance key (`Spinup_Provenance`: the cell,
its soil inputs, the workbook case, the MOD and ERW rows with their
forcing, meteo and transient windows, and what the run starts from --
cold, or the lateral record and tracking run it replays), so a changed
input starts that work afresh instead of resuming a stale run. A finished
state also records the instant it hands off at, and the field transient
(`Run_One_Subfield`) refuses a state for any other start. The scripts load
`${MATLAB_MODULE:-MATLAB/2025a}` and run on `cpumodel:9575f` nodes, the
platform of the v0.9 campaigns in `runtime/`, so a state is produced
where the transient will run it. The code needs the Statistics and
Machine Learning Toolbox (`TreeBagger`, `pdist2`) and the Mapping Toolbox
(`readgeoraster`, `projcrs`) -- check `ver` once on a compute node.

## Site settings and warm starts

A site other than the workbook's borrows a case (`TC_CASE`) and replaces
what differs through two one-row CSVs: `TC_MOD_COLUMNS_FILE` (location,
soil column, vegetation file, CO2 and deposition records, forcing and
windows) and `TC_ERW_COLUMNS_FILE` (the ERW settings). A number stays a
number, `NaN` included, anything else is text, as the workbook holds it;
the `TC_FORCING_FILE` ... `TC_TRANSIENT_END` variables still take
precedence. Each cell's own soil and land cover come only from
`TC_AGRONOMIC_FILE`: a column file that sets them is refused, and so is a
column nothing reads. `field_spinup_hamden.env.example` is the Hamden
site set up this way, on a 5 m grid under the majority rule
(`TC_CELL_M=5`, `TC_MIN_COVERAGE=majority`: each square kept when more
than half of it lies in the field, on the lattice of least misrepresented
area; main README, *The majority rule*). The per-cell file must be
written on the same grid (`Run_Hamden_Field_Inputs`); the plan refuses
one whose cell ids or centres differ.

Once site observations replace a field's provisional soil inputs, point
`TC_WARM_START_STATES` at the first run's `CellBGStates.mat` (copied out
of scratch), give a new `TC_OUTPUT_DIR` and submit again. The plan
refuses states of another grid or missing a cell. The class-seed array
then has nothing to do, the tracking starts from every cell's own
converged state, and each replay (or, with `per_cell`, each cell) runs
Phase 3 only; every handoff is built from the new inputs.

## The field transient

Once `CellBGStates.mat` is assembled, the same environment file runs the
field transient (`Run_Field_Transient`; main README, *Field transient*):
every cell from its handed-off state over the transient window, in 2-D
with the ERW chemistry on and no rock.

```bash
cd ~/project_pi_<netid>/TC-SMEW-N/v2/cluster
bash submit_field_transient.sh field_spinup_hamden.env            # states on disk
bash submit_field_transient.sh field_spinup_hamden.env <job id>   # or after the assemble job
```

It reads `TC_STATES_FILE` (default `$TC_OUTPUT_DIR/CellBGStates.mat`) on
the spin-up's grid, columns and DEM, and writes to
`TC_TRANSIENT_OUTPUT_DIR` (default `$TC_OUTPUT_DIR/transient`): a
`segment_<k>.mat` for every `TC_SEGMENT_DAYS` (30) with the next one's
`carry_<k>.mat`, then `Field_Transient_<from>_<to>.mat` and its field means
as a CSV. `field_transient.sbatch` asks for one node of the `week`
partition for four days, with 32 CPUs and 64 GB. MATLAB runs
single-threaded, as the v0.9 2-D campaigns do, and the columns' daily ERW
chemistry runs on `TC_PARALLEL_WORKERS` thread workers (default one per
CPU). Near its limit a job requeues itself, up to
`TC_TRANSIENT_MAX_REQUEUES` (40) times, and resumes: finished segments are
reused and the running one restarts from the driver's rolling checkpoint
(`TC_TRANSIENT_CKPT_MINUTES`, 60). `TC_TRANSIENT_JOBS` (2) chained copies
are the fallback where requeueing is refused.
`TC_TRANSIENT_CEC="6.4 10.8 13.7"` submits one chain per CEC scenario,
each into the output folder suffixed `_CEC_<value>`.
`TC_TRANSIENT_PARTITION`, `TC_TRANSIENT_TIME`, `TC_TRANSIENT_MEM` and
`TC_TRANSIENT_CPUS` override the batch script's resources.

Run a smoke test first, into its own folder: `TC_TRANSIENT_NDAYS=3`,
`TC_TRANSIENT_CELL_IDS="33 34 35 62"` and
`TC_TRANSIENT_OUTPUT_DIR=$TC_OUTPUT_DIR/transient_smoke`.

The run is long (main README, *Field transient*, *Cost*). Its daily ERW
chemistry is about 10,000–19,000 core-hours, which the workers share. Its
hourly stepping, roughly a week here, runs one column after another
whatever the worker count. With 32 workers expect two to four weeks, so
many requeues. Before committing to it, run two days of the whole field
into its own folder (`TC_TRANSIENT_NDAYS=2`). Its log gives these cores'
set-up time (936 parameter files), the driver's time per column-step and
the segment's wall time. `TC_TRANSIENT_CPUS=64` halves the ERW part. A
30-day segment of the 936 columns holds about 7.4 GB, which the rolling
checkpoint writes every `TC_TRANSIENT_CKPT_MINUTES`.

## The robust test

The forest exists so that only a subsample of the field's cells needs a
replay. Its justification is empirical, and the cluster makes the
experiment affordable:

1. `TC_SUBSAMPLE=1` on the diagnostic grid (54 cells at 15 m): every
   cell tracked and replayed. That is the truth.
2. `validate` (`TC_VALIDATE=1`): for subsample fractions 0.10-0.30 and
   five seeds, draw the stratified subsample from the truth, train the
   forest on it, and report the held-out relative error of every
   conserved total against Lian et al.'s 10 % target, beside the forest's
   own k-fold estimate and its skill against the training mean and the
   nearest tracked cell; check the exchangeable Ca/Mg/K of the handoff
   seed against the measured base saturation on the cells that carried
   it. `RF_Validation_<polygon>.csv/.png`,
   `Base_Saturation_Attainment_<polygon>.csv`.
3. A second run at the production fraction (`TC_SUBSAMPLE=0.10`) into
   another output directory, then `validate` with `TC_RF_STATES_DIR`
   pointing at it: every forest-reconstructed cell's conserved totals and
   handoff-gate result against its own replayed truth
   (`RF_Path_vs_Truth_<polygon>.csv`).

`TC_METHOD=per_cell` on the same grid gives the 1-D reference without
lateral exchange; the difference between it and the truth of step 1 is
the lateral redistribution the field imposes on the pools. The
production grid (1461 cells at 2.33 m on Control 2) runs through the same
scripts by dropping `TC_CELL_M` and `TC_MIN_COVERAGE`.

## Caveats carried from the workstation

- The workbook has no Bruce Farm case; the environment file supplies
  the Kantola forcing and windows as spin-up columns. Bruce Farm's own
  settings (forcing, soils, management) would go in column files, as
  Hamden's do.
- Base saturation reaches the ERW exchanger through the per-cell
  `f_CEC_in` seed (main README, *Base saturation*); Na and Al exchange
  start at zero. A cell without one on a `spinup_equilibrium` row (Hamden)
  takes the exchanger in equilibrium with its spun-up soil solution.
- Every spin-up stage is ERW-free (the v0.9 driver gates its chemistry on
  `OPT_ERW`); each cell's ERW seed is built and checked once, from the
  first two forcing days of the transient window, by its final handoff
  diagnostic.
- The 2-D driver steps its columns one after another every hour; the
  tracking run's wall time is therefore proportional to the cell count.
  Running the hourly column loop on parallel workers would divide it by
  the worker count and is the obvious next speed-up for the production
  grid (it would serve the field transient too).
