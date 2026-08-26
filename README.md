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

### Bruce Field ERA5 forcing

`Preprocess_Bruce_Field_ERA5` reads the split ERA5-Land time-series NetCDF
files and writes `Data_Bruce_Field_Spinup_2014_2023.mat` and
`Data_Bruce_Field_Transient_2024_2025.mat` in `v2/T&C_Code/Inputs`.
The period labels use UTC boundaries. `Date` is fixed Eastern Standard Time
(GMT-5), as required by the radiation partitioning calculation.

From the repository root:

```matlab
addpath(fullfile('v2','T&C_Code','Inputs','Preprocess'))
Preprocess_Bruce_Field_ERA5
```

The ERA5 files contain the grid point at 36.7 N, 78.5 W. The preprocessing
uses a 127.4689 m USGS 3DEP elevation at that point. Negative precipitation
and shortwave values within the observed NetCDF packing noise are set to zero.

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
subfield. Each selected record is split into connected regions first. Only
regions that touch or overlap are unioned, so duplicate area is not counted
twice and distant polygons cannot alter one another's vertices. Every region
gets a separate grid, simulation, checkpoint, result, status file and log. A
large empty grid is therefore not placed across the gap between fields.
Geographic and projected source coordinates are converted to the local WGS84
UTM zone before any distance or area calculation. Polygon simplification is
disabled, collinear source points are retained, and the projection round trip
must close within 1e-6 m. This requires MATLAB Mapping Toolbox.

The automatic grid retains every full square with positive polygon overlap.
Squares in the bounding rectangle that do not touch the field are removed.
For polygon area `A`, retained cell count `N` and square width `h`, coverage is

```
field_fraction = A/(N*h^2)
```

and must be greater than 0.95. The grid is aligned to the minimum-area rotated
bounding rectangle. Cell-size candidates occur where a bounding width or
height fits another complete square. The calibration minimizes `N` first,
then minimizes `N*h^2-A`. There is no maximum cell count. `GridCellSizeM` can
still prescribe a width when a fixed resolution is needed, but the setup is
rejected if that width does not exceed 95% coverage.

Every retained square is simulated in full, including its area outside the
polygon. The polygon edge is only an accounting boundary. It is not a no-flow
boundary, so the outside parts of boundary cells remain part of lateral and
vertical transport. Missing neighbours around the retained-cell perimeter
form the closed lateral outer boundary. Vertical drainage and leaching remain
active.

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
m2. Geometry calibration gives 2.330 m cells for `Control 2` and 2.621 m cells
for `Control 1`. The independent setups contain 1630 and 2340 retained full
cells. Their in-polygon fractions are 95.012% and 95.105%, with 441.50 and
786.55 m2 outside. The retained polygon area is 23693.02 m2.

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
outside-area summaries for each polygon. The CSV includes the source and
projected vertex counts, projection round-trip error and source-to-polygon
vertex error.

The result topology stores the UTM EPSG code, grid rotation, projected
origin, overlap areas, and projected cell centres. `Grid_Local_To_Projected`
maps any grid-local output coordinates back to the shapefile projection.

### Cell soil inputs and automatic spin-up

Each disconnected polygon has its own agronomic input file and spin-up state
file. The files are tied to the shapefile-derived cell ids and projected cell
centres. This stops a state prepared for one polygon or grid orientation from
being used on another one.

For property `j` in cell `c`, the resolved input is

```
x(c,j) = use_measured(c,j)*measured(c,j) + ...
         (1-use_measured(c,j))*baseline(c,j)
```

The four toggles are `use_measured_bulk_pH`, `use_measured_texture`,
`use_measured_OM` and `use_measured_CEC`. A false toggle always uses the
interpolated baseline for that cell. A true toggle requires a finite measured
value. Texture uses numeric mass fractions:

```
silt_fraction = 1 - sand_fraction - clay_fraction - OM_fraction
```

All four fractions must be nonnegative. OM changes the physical soil organic
matter fraction used by the pedotransfer functions. It does not directly
repartition POC and MAOC. OM alone does not identify that split, so the carbon,
nitrogen and phosphorus pools are left to the automatic spin-up.

CEC is prescribed and is not calibrated. Its areal exchange inventory is

```
CEC_area = CEC_mass * 1e-5 * rho_bulk * active_depth * 1e6
```

where `1e-5` converts cmol charge kg-1 to mol charge g-1 and `1e6` converts
mol charge to micromol charge. Bulk pH is the laboratory observation. The
pore-water pH and exchange state are obtained by the existing strict inversion.

| Variable | Expected SI Units | Actual Units provided in the code |
|---|---|---|
| measurement toggle | 1 | logical or numeric 0/1 |
| bulk pH | 1 | pH index, dimensionless |
| hydrogen activity | 1 | `10^(-pH)` |
| sand, clay, silt, OM | kg kg-1 | decimal mass fraction |
| organic carbon fraction | kg C kg-1 soil | `OM_fraction/1.72` |
| CEC by soil mass | mol charge kg-1 | **Non-SI input:** cmol charge kg-1, multiplied by 0.01 |
| dry bulk density | kg soil m-3 | `rsd` in kg m-3, converted internally to g m-3 where required |
| active depth | m | **Non-SI input:** `Zbio` in mm, divided by 1000 |
| areal CEC inventory | mol charge m-2 | **Scaled unit:** micromol charge m-2 in the ERW arrays |
| biogeochemical pools | kg element m-2 | **Non-SI storage:** g element m-2 in the 79-pool state |
| projected cell centre | m | WGS84 UTM zone 17N metres |
| polygon accounting area | m2 | m2 |
| in-polygon grid fraction | 1 | `field_area_m2/simulated_area_m2`, required to be greater than 0.95 |

Templates are written directly from the shapefile geometry. Field 7 remains
excluded by selecting source features 1 and 2:

```matlab
shape_file = fullfile('v2','Spatial inputs','Bruce Farm Virginia', ...
    'Whole farm','AdamWilbourne_Control_Fields.shp');
addpath(fullfile('v2','T&C_Code'),fullfile('v2','spinup'))
files = Write_Field_Agronomic_Templates(shape_file, ...
    fullfile('v2','T&C_Code','Inputs','Field_Spinup','Bruce_Farm'),[],[1 2]);
```

Fill every baseline column in both files. Measurements may stay blank when
their toggle is false. Automatic spin-up is then run separately for each
polygon:

```matlab
[state_1,config_1] = Run_Field_Cell_Spinup(shape_file,files{1}, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',1);
[state_2,config_2] = Run_Field_Cell_Spinup(shape_file,files{2}, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',2);
```

The site workbook row must provide `meteo_file`, `meteo_start` and `meteo_end`
before spin-up can run. This is a hard preflight check. The imported machinery
is from `Model/TC-SMEW-N` commit `31e9e13`. It retains the three spin-up phases,
checkpoint rejection, convergence checks, 79-pool physicality screen, six
nutrient budgets and ERW handoff gate.

The transient launcher requires the matching agronomic and state file for
each polygon:

```matlab
Run_All_12_Parallel('Runs',1:3,'FieldShapefile',shape_file, ...
    'FieldFeatureIndices',[1 2], ...
    'CellAgronomicFiles',{files{1},files{2}}, ...
    'CellSpinupStateFiles',{state_1,state_2}, ...
    'MaxConcurrent',2,'MaxRuntimeHours',168)
```

Every retained full cell needs a baseline, including boundary cells that
extend outside the field polygon. Polygon overlap area is used only when
field totals and means are calculated.

Runs without cell soil files use the control Kantola T&C parameters, forcing
and archived initial state. The active biogeochemical depth is 300 mm in those
runs, matching the workbook value and the spin-up that produced
`Data_Kantola_BG.mat`. Runs with cell soil files use the matching validated
automatic spin-up state in each cell. Treatment rows only supply the amendment
properties. Untreated cells use the control chemistry row.
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
