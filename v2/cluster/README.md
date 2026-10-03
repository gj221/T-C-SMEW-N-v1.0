# Field spin-up on Bouchet

The automatic spin-up costs about 48 s per forcing year in its T&C-only
Phase 1, about 6 min in the stand-alone Phase 2, and about 1 min per
coupled year in the ERW-free Phase 3, whose convergence test runs over
complete rings of cycle years (measured on one Control 2 cell, this
workstation). A cell is a job of hours, and a field is a SLURM array of
them. Everything here runs one cell per array task and reuses finished
cells, so a preempted or resubmitted array picks up where it stopped.

## What runs

| step | script | what it does |
|---|---|---|
| plan | `field_spinup_step.sbatch` (`TC_TASK=plan`) | `Run_Field_Cell_Spinup(...,'PlanOnly')`: `spinup_plan.csv` and the array-index lists |
| cells | `field_spinup_cell.sbatch` (array) | one cell per task: a full spin-up (`per_cell`, or an RF subsample cell), or the RF-seeded spin-up (`rf`, after the subsample array) |
| assemble | `field_spinup_step.sbatch` (`TC_TASK=assemble`) | gathers `cellstate_<id>.mat` into the validated `CellBGStates.mat` the transient run reads |
| validate | `field_spinup_step.sbatch` (`TC_TASK=validate`) | `Run_RF_Validation_Experiment` on the assembled truth |

`submit_field_spinup.sh` chains them with SLURM dependencies. All
settings are environment variables in `field_spinup.env` (copy and edit
`field_spinup.env.example`); `Cluster_Field_Spinup_Task.m` is the only
MATLAB entry point.

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
   an array finishes.
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
tail logs/spinup_cell_<jobid>_<task>.out         # one cell's log
```

The array uses the `day` partition with a one-day limit per task and
8 GB, which is ample: the DEM loader crops each 1 m lidar tile to the
field window as it reads it, so a task holds one tile (about 0.8 GB)
at its peak, not the six tiles that cover the farm. If a task is
preempted or times out, resubmit the same array: a finished cell is
reused, and a cell in progress resumes its own phase checkpoints. Both
are recognised by a provenance key (`Spinup_Provenance`: the cell, its
soil inputs, the workbook case, the MOD and ERW rows with their forcing,
meteo and transient windows, and the seed), so a changed input or
environment starts that cell afresh instead of resuming a stale run. A
finished state also records the instant it hands off at, and the field
transient (`Run_One_Subfield`) refuses a state for any other start.
The scripts load `${MATLAB_MODULE:-MATLAB/2025a}` and run on
`cpumodel:9575f` nodes, the platform of the v0.9 campaigns in
`runtime/`, so a state is produced where the transient will run it.
The code needs the Statistics and Machine Learning Toolbox
(`TreeBagger`, `pdist2`) and the Mapping Toolbox (`readgeoraster`,
`projcrs`) -- check `ver` once on a compute node.

## The robust test

The random-forest path exists so that a field of N cells needs only a
topographically stratified subsample of full spin-ups. Its
justification is empirical, and the cluster makes the experiment that
justifies it affordable:

1. `TC_METHOD=per_cell` on the diagnostic grid (54 cells at 15 m): every
   cell spun up in full, in parallel. That is the truth.
2. `validate` (automatic after assembly): for subsample fractions
   0.10-0.30 and five seeds, draw the stratified subsample from the
   truth, train the forest on it, and report the held-out relative
   error of every conserved total against Lian et al.'s 10 % target,
   beside the forest's own k-fold estimate; check the exchangeable
   Ca/Mg/K at handoff against the measured base saturation on the
   cells that carried it. `RF_Validation_<polygon>.csv/.png`,
   `Base_Saturation_Attainment_<polygon>.csv`.
3. A second run with `TC_METHOD=rf` into another output directory, then
   `validate` with `TC_RF_STATES_DIR` pointing at it: every RF-seeded
   cell's conserved totals, exchangeable pools and handoff-gate result
   against its own truth (`RF_Path_vs_Truth_<polygon>.csv`).

That tells you, per fraction, how much of the field's variance the
forest recovers and whether the seeded spin-up lands where the full one
does -- the two numbers the production grid (1461 cells at 2.33 m on
Control 2) needs before its subsample size is chosen. The production
grid runs through the same scripts by dropping `TC_CELL_M` and
`TC_MIN_COVERAGE`.

## Caveats carried from the workstation

- The workbook has no Bruce Farm case; the environment file supplies
  the Kantola forcing and windows as spin-up columns. A site row (Bruce
  forcing, soils, management) is the next step for a field product.
- Base saturation reaches the ERW exchanger through the per-cell
  `f_CEC_in` seed (main README, *Base saturation*); Na and Al exchange
  start at zero.
- Phase 3 is ERW-free (the v0.9 driver gates its chemistry on
  `OPT_ERW`); the ERW seed is built and checked once, from the first two
  forcing days of the transient window, by the final handoff diagnostic.
