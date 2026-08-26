# T&C-SMEW-N v1.0

## Repository layout

```
v2/run/          launchers: ERW_run, Run_Matrix_12, Run_One_Subfield,
                 Run_All_12_Parallel, Run_Subfield_Worker,
                 Terminate_Process_Tree
v2/plotting/     figure scripts: Plot_1D_Then_All_12, Plot_Grid_Setup,
                 Plot_Biogeo_Pools, Plot_Nflux_Diagnostics,
                 Plot_Kantola_Benchmark, and the Cflux_post_processing helper
v2/figures/      every PNG and FIG the plotting scripts write
v2/Spatial inputs/   field polygons used by full-grid simulations
v2/T&C_Code/     model source, Inputs, Parameters and Tests
v2/SIM/          simulation results (not tracked)
v2/logs/         worker logs and status files (not tracked)
v2/Experimental observations/   digitised field data used by the benchmark
```

Every plotting entry point is a runnable script: `cd` into `v2/plotting` and
run it by name, with no arguments. Each one puts the folders it needs on the
path, so it works from any working directory. `Cflux_post_processing` stays a
function because the benchmark script calls it as a helper.

## Simulations

The twelve simulations are defined in `v2/run/Run_Matrix_12.m`. Each group
contains control, basalt and calcite treatments.

| Runs | Scenario group | Domain and application |
|---:|---|---|
| 1 to 3 | whole field | 12 cells, amendment on the main field and 5 m margin |
| 4 to 6 | margin with crop | 12 cells, amendment on the cropped 5 m margin |
| 7 to 9 | margin with prairie | 12 cells, amendment on the prairie 5 m margin |
| 10 to 12 | vertical 1-D | one 3600 m2 column, amendment across the plot, no lateral transport |

The retained 12-cell geometry has four 25 m by 25 m field cells, four
50 m by 5 m side margins and four 5 m by 5 m corners. The prairie uses
big bluestem as the high C4 PFT and purple coneflower as the low C3 PFT.

### Irregular field grids

Runs 1 to 3 can use a polygon shapefile instead of the retained square
subfield. Selected records are unioned first, so duplicate or overlapping
features are not counted more than once. The union is then split into
connected polygon regions. Every region gets a separate grid, simulation,
checkpoint, result, status file and log. A large empty grid is therefore not
placed across the gap between fields. Geographic and projected source
coordinates are converted to the local WGS84 UTM zone before any distance or
area calculation. This requires MATLAB Mapping Toolbox.

The default cell width comes from the convex envelope of each connected
polygon. If its area is `A` and perimeter is `P`, the square width is
`2*A/P`. The orientation search checks the convex-envelope edge directions
and every angle where the padded row or column count changes. The retained
orientation has the smallest full square-grid area, with the minimum bounding
area used to break a tie. This avoids making cells small because of holes,
narrow gaps or noisy boundary detail. There is no maximum cell count.
`GridCellSizeM` can still prescribe a width when a fixed resolution is needed.

Every square in the rotated bounding grid is simulated, including squares
outside the polygon. The polygon edge is only an accounting boundary. It is
not a no-flow boundary, so laterally coupled cells outside the field can
affect cells inside it. The outer rectangular grid remains closed laterally,
while vertical drainage and leaching remain active.

For cell `i`, the stored accounting area is
`area(intersect(cell_i, field_polygon))`. An areal result `q_i` is converted
to a polygon total with `sum(q_i*cell_overlap_area_i)`. Polygon means use the
same overlap areas as weights. Site totals are the sum of the separate polygon
totals. pH is weighted through hydrogen activity, not by taking an arithmetic
mean of pH values. Full cell areas remain in `topo.cell_area_m2` for storage
and lateral transport. `topo.cell_overlap_area_m2` is used only for field
accounting.

The whole-site source contains three control polygons. Source feature 3,
`Control 3`, is associated with Field 7 and is excluded. Source features 1
and 2 are retained and run separately. Their areas are 8410.29 and 15282.73
m2. Geometry calibration gives a 46.29 m, 4 by 3 grid for `Control 2` and a
67.33 m, 3 by 3 grid for `Control 1`. This is 21 simulated cells across the
two independent setups. The retained polygon area is 23693.02 m2.

From the repository root:

```matlab
shape_file = fullfile('v2','Spatial inputs','Bruce Farm Virginia', ...
    'Whole farm','AdamWilbourne_Control_Fields.shp');
feature_indices = [1 2];
addpath(fullfile('v2','T&C_Code'),fullfile('v2','run'))
setups = Build_Shapefile_Grid_Setups(shape_file,[],feature_indices);
Run_One_Subfield(1,'FieldShapefile',shape_file, ...
    'FieldFeatureIndices',feature_indices,'FieldPolygonIndex',1,'NDays',1)
```

The parallel launcher expands each requested scenario into one job per
connected polygon:

```matlab
Run_All_12_Parallel('Runs',1:3,'FieldShapefile',shape_file, ...
    'FieldFeatureIndices',feature_indices, ...
    'MaxConcurrent',2,'MaxRuntimeHours',168)
```

To inspect the calibrated grid before a run:

```matlab
cd(fullfile('v2','plotting'))
Plot_Field_Grid
```

This writes one subplot per retained polygon, with the complete boundary over
the full grid. It also writes a CSV with separate area, cell-size, grid and
outside-area summaries for each polygon.

The result topology stores the UTM EPSG code, grid rotation, projected
origin, overlap areas, and projected cell centres. `Grid_Local_To_Projected`
maps any grid-local output coordinates back to the shapefile projection.

All runs use the control Kantola T&C parameters, forcing and initial
state. The active biogeochemical depth is 300 mm in every run, matching the
workbook value and the spin-up that produced Data_Kantola_BG.mat. Treatment rows only supply
the amendment properties. Untreated cells use the control chemistry row.
The configuration workbook contains only the five Kantola cases used by
the single-column comparison and the twelve-run matrix.

Field reactions and microbial rates use pore-water pH. Bulk-extract pH is
retained only as the 1:2.5 laboratory observation operator, with no additive
baseline offset. The Kantola seed explicitly uses the prescribed field-state
mode, so its initial bulk pH is a model prediction rather than a forced match.

The archived Kantola background file contains areal pool totals and no
depth metadata. It is reused unchanged in every case. This keeps the runs
comparable, and the 300 mm active depth is the depth it was spun up at.

The 12-cell runs exchange subsurface water with a finite-volume
Buckingham-Darcy calculation. Water advects the dissolved pools used by
`Biogeo_Leakage`: DOC, DON, DOP, NH4, NO3, mineral P, solution K, Ca, Mg
and Si. ERW-only Na, Al, residual anion and DIC use the same active-depth
water flux. Soil-gas CO2 diffuses across shared faces. Soil heat moves by
Fourier conduction and with liquid-water sensible heat. The outer plot
boundary is closed laterally. Existing vertical drainage and solute
leaching remain active. Frozen hydraulic conductivity uses the Hansson
et al. (2004) ice-impedance factor. Fluxes blocked by frozen pore space
are excluded when the adaptive lateral time step is calculated.

From PowerShell at the repository root:

```powershell
matlab -batch "cd('v2/run'); Run_All_12_Parallel('Runs',10:12,'MaxConcurrent',3,'MaxRuntimeHours',168)"
matlab -batch "cd('v2/run'); Run_All_12_Parallel('Runs',1:9,'MaxConcurrent',2,'MaxRuntimeHours',168)"
matlab -batch "cd('v2/plotting'); Plot_1D_Then_All_12"
```

The plotting command first draws the vertical 1-D weathering and nitrogen
figure (`vertical_1d_weathering_nitrogen.png` and `.fig`) and then writes
`Subfield_12_summary.png`, `Subfield_12_summary.fig`,
`Subfield_12_time_dynamics.png` and `Subfield_12_time_dynamics.fig` to
`v2/figures`, reading the results from `v2/SIM/Subfield_12`. Set
`FigureDirectory` beforehand to write the figures somewhere else. The
all-12 stage is skipped with a warning listing the missing runs until the
2-D results finish. The time plot
uses two columns per scenario group and no more than three stacked rows.
It shows control, basalt and calcite through time for cumulative N2O,
cumulative nitrate leaching, the soil nitrate pool, pore-water pH, NEE
and cumulative crop yield. NEE is respiration minus crown-weighted NPP;
positive values are carbon release to the atmosphere. Crop yield is
harvested fruit/grain carbon in g C m^-2, so prairie PFTs contribute zero
yield. Both carbon metrics use the full 3600 m^2 domain area as the
denominator.

The completed 1-D runs can be plotted before the 2-D runs finish, from
`v2/plotting`:

```matlab
Runs = 10:12; Plot_1D_Then_All_12
```

This writes updated six-panel PNG and FIG outputs for runs 10 to 12 only.
When runs 1 to 9 finish, the three 2-D groups can be plotted without
reloading the 1-D results:

```matlab
Runs = 1:9; Plot_1D_Then_All_12
```

This writes PNG and FIG copies of `Subfield_2D_summary` and
`Subfield_2D_time_dynamics`. Running the script without setting `Runs`
plots the 1-D group first and then all 12 scenarios together.

From an open MATLAB session at the repository root:

```matlab
cd v2/run
Run_All_12_Parallel('DryRun',true)
Run_All_12_Parallel('Runs',10:12,'MaxConcurrent',3,'MaxRuntimeHours',168)
Run_All_12_Parallel('Runs',1:9,'MaxConcurrent',2,'MaxRuntimeHours',168)
cd ../plotting
Plot_1D_Then_All_12
```

The first command completes the three vertical 1-D comparisons before the
2-D runs begin. Stop between the two commands if the second stage needs to
wait. `Ctrl+C` stops active child runs; completed result files are retained
and skipped when the same command is run again. On Windows the launcher
stops the complete worker process tree, so a MATLAB worker is not left behind.
Active runs restart from the beginning because transient checkpoints are not
written.

Results go to `v2/SIM/Subfield_12`. Logs and worker status files go to
`v2/logs/Subfield_12`. Existing results are skipped unless
`'Overwrite',true` is passed. Use overwrite after a solver change so all
comparisons use the same code. A one-day check uses:

```matlab
cd v2/run
Run_All_12_Parallel('MaxConcurrent',2,'NDays',1,'Overwrite',true)
addpath(fullfile('..','T&C_Code','Tests')); test_subfield_12_smoke_results(1)
cd ../plotting
NDays = 1; Plot_1D_Then_All_12
```

The September 2019 mineral-N regression uses a small saved fixture:

```matlab
addpath(fullfile('v2','T&C_Code','Tests')); test_no3_immobilization_replay
```

Two 2-D workers is the default for a 32 GB machine. The three 1-D runs can
run together. Full 12-column workers hold several GB of state, so increase
the second-stage limit only after checking available memory.
