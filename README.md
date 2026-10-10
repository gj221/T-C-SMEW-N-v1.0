# T&C-SMEW-N v1.0

v1.0 takes T&C-SMEW-N v0.9 from plot columns to field grids: irregular
polygon grids, kriged per-cell soil inputs, DEM topography with
gravity-driven lateral flow, and the flux-tracking field spin-up of
T&C-BG-2D (Lian et al. 2025), adapted for ERW, for the Bruce Farm
(Virginia) control fields. It carries the v0.9 model of
package w (2026-10-05) unchanged in the Energy Farm (Kantola) runs; the
section *Model lineage from v0.9* at the end records what each v0.9
package changed, and `RUNBOOK_2026-09-24.md` orders the Energy Farm
campaigns on Bouchet.

## Repository layout

```
v2/run/          launchers: ERW_run, Run_Matrix_12, Run_One_Subfield,
                 Run_All_12_Parallel, Run_Subfield_Worker,
                 Terminate_Process_Tree
v2/plotting/     figure scripts: Plot_1D_Then_All_12, Plot_Grid_Setup,
                 Plot_Biogeo_Pools, Plot_Uncertainty,
                 Plot_Topographic_Redistribution
v2/interpolation/  soil-sample kriging and the grid figure it feeds:
                 Run_Notebook_Interpolation, Run_AW_Bruce_Interpolation,
                 Run_Bruce_Interpolation, Run_Yield_Field_Interpolation,
                 Run_Cell_Size_Analysis, Plot_Field_Grid_Interpolation
v2/figures/      every PNG and FIG the plotting scripts write
v2/Spatial inputs/   field polygons used by full-grid simulations
v2/T&C_Code/     model source, Inputs, Parameters and Tests
v2/spinup/       the automatic spin-up and the field spin-up (flux tracking,
                 per cell) with its tracking run, replays and forest
v2/cluster/      SLURM scripts and entry point for the field spin-up on Bouchet
v2/SIM/          simulation results (not tracked)
v2/logs/         worker logs and status files (not tracked)
v2/Experimental observations/   digitised field data used by the benchmark
runtime/         Energy Farm campaign controllers for Bouchet (Base_1D/2D,
                 SSP_1D/2D, Calibration_Sweep_1D, Uncertainty_1D/2D,
                 Worker_Scaling_Test); each has its own README
*.md (root)      v0.9 change notes and RUNBOOK_2026-09-24.md
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
path, so it works from any working directory.

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

#### The majority rule

The coverage rule keeps every square the polygon touches, so it simulates a
band of outside land along the whole boundary, and its bound on that band,
`A/(N*h^2) > 0.95`, is also what sets the width. The majority rule
(`'majority'` wherever the coverage rule is passed: the fourth argument of
`Build_Shapefile_Grid_Setups`, `MinCoverage`, `FieldMinCoverage`,
`TC_MIN_COVERAGE`) settles the grid's three choices each on its own ground:

- **Which squares.** A square is simulated when more than half of it lies
  in the field. Over a given lattice that is the exact minimiser of the
  misrepresented area, the field land left out plus the outside land taken
  in: a kept square adds its outside part `h^2-a`, a dropped one loses its
  field part `a`, so `E = sum(min(a,h^2-a))`. There is no threshold to tune.
- **Where the lattice lies.** `E` still depends on the lattice, and
  `Place_Field_Polygon_Grid` chooses it by the same criterion. Along a
  boundary edge parallel to the grid the error can vanish, so the rotation
  is taken from the directions of the field's own edges at least one cell
  long; for each, the lattice offset comes from a 16 x 16 scan of the
  offset torus and a compass search. Every evaluation is exact: the
  overlap of each square comes from the polygon's boundary integral, each
  edge cut at the grid lines, and matches polygon clipping to 1e-6 m2.
- **How wide.** The width is prescribed from the terrain the lateral flow
  depends on (`Run_Topographic_Grid_Convergence`, `'CoverageRule'`), not
  from coverage; the rule refuses an empty width.

The rotation turns only the lattice. Every cell's terrain (elevation,
slope, aspect from true north, curvature, wetness index, sky view) comes
from the north-up DEM at the cell's true centre, so the vertical processes
and the sun geometry do not see the grid's angle, and the two-point face
fluxes of the lateral exchange resolve a uniform head gradient exactly at
any angle, so turning the grid re-aims the discretisation error without
biasing it. On the Hamden parcel at 5 m, the cells' mean slope, its spread,
the wetness index and its spread, the mean aspect and the gradient the
face fluxes see agree within 1 % at 0, 9.7, 19.3 and 45 degrees, and with
the DEM raster itself. What the rotation does change is the boundary: the
misrepresented area (548 m2 at 19.3 degrees, about 780 at 0 and 45) and
the number of perimeter cells where surface water leaves the grid (8 to
22), which moves where along the low edge the water exits, not how much.

The field land in the dropped squares, `topo.field_missed_area_m2`, is not
simulated. Field means are unchanged (overlap-weighted over the kept
squares); field totals are the mean times the whole field area `A`, so that
land enters at the field mean (`Field_Aggregation_Areas`,
`Domain.field_missed_area_m2`). The default rule is unchanged.
`test_field_grid_majority` checks a rectangle whose answer is known
exactly (two sides on grid lines, the 0.3 m column and 0.2 m row dropped),
the same rectangle turned 23 degrees, and a field with a hole, where the
boundary integral and the polygon clipping agree. On the Hamden parcel at
5 m it keeps 936 squares and misrepresents 2.3 % of the parcel, against
1,224 squares of 4.48 m and 5.0 % under the coverage rule (Hamden, below).

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
addpath(fullfile('v2','interpolation'))
Plot_Field_Grid_Interpolation(shape_file,[],feature_indices)
```

This writes one subplot per retained polygon, with the complete boundary over
the full grid. It also writes a CSV with separate area, cell-size, grid and
outside-area summaries for each polygon. The CSV includes the source and
projected vertex counts, projection round-trip error and source-to-polygon
vertex error. Given a sample table, an attribute and a sample stage the same
function draws the kriged surface under the grid; see the section below.

The result topology stores the UTM EPSG code, grid rotation, projected
origin, overlap areas, and projected cell centres. `Grid_Local_To_Projected`
maps any grid-local output coordinates back to the shapefile projection.

### Soil sample interpolation

`v2/interpolation` is the MATLAB form of `bruce_interpolation.ipynb`. Ordinary
kriging over the measured soil samples produces a continuous surface for one
attribute inside each polygon, on the same rotated grid the simulation uses.

The variogram keeps the pykrige parameterisation the notebook depended on, so
a MATLAB surface matches the Python one value for value. `range` is the
effective range: the exponential model decays over `range/3` and the gaussian
over `range*4/7`, and a `sill` given with a `nugget` is stored as the partial
sill `sill-nugget`. Left unset the parameters take the notebook recipe, a sill
equal to the population variance, an 80 m range and a nugget of 5% of the
variance. The notebook passed its `variogram_model` argument to the plot title
only and always fitted an exponential model; here the argument reaches the fit
and `exponential` remains the default.

A sampling stage is a timepoint, so stages are never pooled. Each polygon is
interpolated separately at each stage that leaves three or more finite values
inside its boundary, and the stages that do not are reported instead of being
dropped in silence.

There are two views, and they answer different questions.

`Run_Notebook_Interpolation` is the notebook figure: the selected fields are
unioned into one boundary, one kriging model is fitted to every sample inside
it, and one surface is drawn in Easting and Northing. Use it to look at the
site.

`Run_AW_Bruce_Interpolation`, `Run_Yield_Field_Interpolation` and
`Run_Bruce_Interpolation` are the simulation view: each disconnected polygon
gets its own grid, its own kriging model and its own surface, drawn under the
calibrated cells it will feed. Use it to prepare a run.

```matlab
addpath(fullfile('v2','interpolation'))
results = Run_Notebook_Interpolation();
report = Run_AW_Bruce_Interpolation();
report = Run_Yield_Field_Interpolation('Attribute','pH');
report = Run_Bruce_Interpolation();
```

`Run_Notebook_Interpolation` and `Run_AW_Bruce_Interpolation` default to
`AdamWilbourne_Fields.shp`, the notebook's own boundary. It is tracked in
`v2/Spatial inputs/Bruce Farm Virginia/Whole farm/` and mirrored in the survey
folder at `Whole farm/Map files/`. Records 65, 77, 115, 116, 117 and 118 are
AW Bruce 4, 2, 1, 6, 13 and 12, the notebook's six fields, selected by
`Farm == "AW Bruce"` and those six `Identifier` values:

```matlab
Run_Notebook_Interpolation('ShapeFile',field_shp, ...
    'AttributeFilters',struct('Farm','AW Bruce','Identifier', ...
    {{'AW Bruce 1','AW Bruce 2','AW Bruce 4','AW Bruce 6', ...
    'AW Bruce 12','AW Bruce 13'}}),'SampleSheet',1);
```

is the default and need not be passed; `ShapeFile` only matters for pointing
at a different boundary. If the tracked copy and the survey-folder copy are
both missing, `Run_Notebook_Interpolation` says so and falls back to the
Yield Field applied strips instead of erroring.

`Run_AW_Bruce_Interpolation` and `Run_Notebook_Interpolation` read sheet 1 of
the agronomic workbook, `Bruce_RJ_msft_yale_geochem - ag`, the notebook's
`bruce.csv`: 841 of its 1871 rows fall inside the six AW Bruce fields, with
every stage (BL, BLP, R1, R2) represented in every field.
`Run_Yield_Field_Interpolation` reads the three Yield Field applied strips and
the Yield Field soil workbook, which carries BL, BLP, R1, R2 and R3.
`Run_Bruce_Interpolation` reads the whole-site control polygons and sheet 2 of
the agronomic workbook, `Wilbourne_RJ_msft_yale_geochem`: sheet 1 covers the
Backside, Frontside and Yield Field block more than a kilometre from Control 1
and Control 2, so none of its rows fall inside them, but sheet 2 is the
whole-grower superset and does reach both. Sample workbooks are read from the
survey folder in OneDrive; set `BRUCE_FARM_ROOT` to read them from somewhere
else.

Each run writes one figure per timepoint, its grid summary CSV, and one
combined CSV of polygon by timepoint with sample counts, cell counts and
value ranges.

### Choosing the cell width

`Run_Cell_Size_Analysis` scores candidate widths for each polygon and writes
the table and the widths it recommends to `v2/figures`. Two rules get called
95% coverage and they do not bound the same error.

```
aggregate coverage   A/(N*h^2) over every cell with positive overlap. Bounds
                     the simulated area the grid adds in total. Says nothing
                     about any one cell: at the calibrated 14.077 m width,
                     Applied Strip 2 passes at 95.565% while its worst cell
                     is 0.92% inside the polygon.
per-cell fraction    overlap/h^2 for one cell on its own. A threshold on it
                     bounds how much of any single simulated column is
                     outside the field, and so bounds the area correction
                     that column needs. It is met by dropping the cells that
                     fail, which leaves a thin unsimulated band at the edge.
capture              retained overlap area over polygon area. The cost of
                     the per-cell rule, and the quantity that actually
                     decides the width.
```

Under the per-cell rule the retained band is what matters, and it narrows in
proportion to `h`: the missed area is close to `0.2*h` times the polygon
perimeter. The threshold value itself barely binds, because at any width
below about 2 m the overlap distribution is bimodal, and a kept cell is
typically more than 99.9% inside whether the threshold is 90%, 95% or 99%.
The width is the decision; the threshold is not.

Measured capture under a 95% per-cell rule:

| Polygon | Area | Calibrated | h = 2 m | h = 1 m | h = 0.5 m |
|---|---|---|---|---|---|
| Applied Strip 1 | 1697 m2 | 34 cells, 91.3% | 373, 87.9% | 1631, 96.1% | 6613, 97.4% |
| Applied Strip 2 | 8712 m2 | 43 cells, 95.9% | 1989, 91.3% | 8341, 95.7% | 34291, 98.4% |
| Applied Strip 3 | 2697 m2 | 58 cells, 95.0% | 606, 89.9% | 2638, 97.8% | 10590, 98.1% |
| Control 2 | 8410 m2 | 1461 cells, 94.3% | 1988, 94.5% | 8190, 97.4% | 33241, 98.8% |
| Control 1 | 15283 m2 | 2129 cells, 95.6% | 3693, 96.6% | 15032, 98.4% | 60651, 99.2% |
| AW Bruce 4 | 38887 m2 | 3785 cells, 96.0% | — | — | 154618, 99.4% |
| AW Bruce 2 | 37216 m2 | 2328 cells, 95.0% | — | — | 148038, 99.4% |
| AW Bruce 1 (part A) | 24374 m2 | 1985 cells, 95.0% | — | — | 96889, 99.4% |
| AW Bruce 1 (part B) | 13238 m2 | 1351 cells, 94.6% | — | — | 52520, 99.2% |
| AW Bruce 6 + 13 | 86648 m2 | 3550 cells, 94.9% | — | — | 344784, 99.5% |
| AW Bruce 12 | 27537 m2 | 2397 cells, 96.1% | — | — | 109540, 99.4% |

The three Yield Field strips are near-rectangles, so their coarse calibrated
widths pass the aggregate rule on alignment alone and the capture column is
the honest measure. Capture is not monotonic in `h` above about 1 m, because
the grid is re-oriented to the minimum-area rotated rectangle at every width
and a small change in angle moves many cells across the threshold at once.
Below 1 m the curve is smooth.

`AdamWilbourne_Fields.shp` is the least regular boundary here — 30 to 100
vertices per field against 9 to 30 for the strips and controls, and one
record, `AW Bruce 1`, is itself two disconnected parts under a single
`Identifier`, gridded separately. Their calibrated widths keep 94.6 to 96.1%
under the 95% per-cell rule, already close to the target because
`Calibrate_Field_Grid_Cell_Size` orients the grid once, at a width estimated
from the field's own extent, then searches width within that fixed
orientation; a width prescribed on its own, as 0.5 m is, re-orients fresh for
that width and lands closer to 99.2–99.5% captured, consistent with the other
two sites.

**The configured width is 0.5 m.** Half a metre sits well inside the soil
survey it is fed from: samples are 10 to 40 m apart and the variogram range
is 80 m, so nothing below about 2 m adds information about the soil, only
about the boundary. Going to 0.25 m quadruples the cell count for roughly one
more point of capture. Under the 95% per-cell rule every simulated column is
at least 95% inside its polygon, in practice above 99.9%, and 98 to 99% of
each polygon is simulated.

| Polygon | Cells at 0.5 m | Aggregate inside | Kept at >=95% | Captured |
|---|---|---|---|---|
| Applied Strip 1 | 6904 | 98.35% | 6613 | 97.38% |
| Applied Strip 2 | 35561 | 97.99% | 34291 | 98.40% |
| Applied Strip 3 | 11062 | 97.54% | 10590 | 98.15% |
| Control 2 | 34224 | 98.30% | 33241 | 98.81% |
| Control 1 | 61785 | 98.94% | 60651 | 99.21% |
| AW Bruce 4 | 156821 | 99.19% | 154618 | 99.40% |
| AW Bruce 2 | 149914 | 99.30% | 148038 | 99.44% |
| AW Bruce 1 (part A) | 98299 | 99.18% | 96889 | 99.38% |
| AW Bruce 1 (part B) | 53510 | 98.96% | 52520 | 99.18% |
| AW Bruce 6 + 13 | 349036 | 99.30% | 344784 | 99.48% |
| AW Bruce 12 | 111004 | 99.23% | 109540 | 99.44% |

53527 cells over the three Yield Field strips, 96009 over the two control
polygons, 918584 over the six AW Bruce fields. Pass the width to any driver:

```matlab
Run_AW_Bruce_Interpolation('CellSizeM',0.5);
Run_Yield_Field_Interpolation('CellSizeM',0.5);
Run_Bruce_Interpolation('CellSizeM',0.5);
```

The residual 1 to 2% band should be carried by scaling each polygon total by
polygon area over captured area, in the same way `cell_overlap_area_m2`
already scales a partial cell.

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
dem = fullfile('v2','Spatial inputs','Bruce Farm Virginia','DEM');
[state_1,config_1] = Run_Field_Cell_Spinup(shape_file,files{1}, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',1,'FieldDEM',dem);
[state_2,config_2] = Run_Field_Cell_Spinup(shape_file,files{2}, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',2,'FieldDEM',dem);
```

The site workbook row must provide `meteo_file`, `meteo_start` and `meteo_end`
before spin-up can run (an hourly forcing `.mat` is accepted as
`meteo_file`; a row that lacks them can be completed per cell with the
launcher's `MODParamColumns` option, see *Spin-up forcing* below). This
is a hard preflight check. The imported machinery
is from `Model/TC-SMEW-N` commit `31e9e13`. It retains the three spin-up phases,
checkpoint rejection, convergence checks, 79-pool physicality screen, six
nutrient budgets and ERW handoff gate. The field spin-up around it is the
flux-tracking spin-up of T&C-BG-2D (see *Flux-tracking field spin-up*);
`'Method','per_cell'` runs every cell's own spin-up instead, without lateral
exchange.

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

Aqueous ERW speciation and mineral reactions use pore-water pH. T&C microbial
rates, including nitrification and denitrification, use the reconstructed 1:2.5
DI-water bulk-extract pH, following the Beerling et al. preprint soil-pH method.
The Kantola seed uses the measured pretreatment bulk
pH as an inversion boundary while retaining the prescribed exchange inventory.
There is no additive pH offset.

The field charge balance includes dissolved Al hydroxide ions. `Alk_tot` is the
conserved strong-ion inventory; pore-water `Alk` is carbonate alkalinity.
Primary calcite contributes to the carbonate, not silicate, alkalinity source
diagnostic.

ERW cation partitioning follows the selected pH exchange scheme. `TC_GT` keeps
the legacy T&C fixed-pool rates and disables only its kinetic solution-exchange
rates. `appelo` disables both legacy rate pairs, leaving the conserved
SCEPTER/Appelo equilibrium as the only K, Ca and Mg partition operator.
Standalone T&C runs keep all legacy rates.

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

### Nitrogen transformation contract

The active ammonium and nitrate pools are `B(31)` and `B(32)`. They are
converted from g N m^-2 to g N per g dry soil inside
`BIOGEOCHEMISTRY_DYNAMIC3` and converted back to g N m^-2 d^-1 for output.
The corrected gross denitrification rate is:

```text
DOC_v = (B(12) + B(13)) * rho_soil                         [g C m^-3]
f_DOC = DOC_v / (DOC_v + Kc), Kc = 17.5 g C m^-3
DENIT = k_den * f_T * f_water * f_den_pH * f_DOC * B(32)  [g N gsoil^-1 d^-1]
```

This retains the first-order nitrate reservoir loss in Dickinson et al. (2002)
and legacy T&C-SMEW. The Li et al. (1992) nitrate half-saturation parameter is
volumetric, `Kn = 82.5 g N m^-3`, but remains inactive because TC-SMEW does not
carry the denitrifier biomass, yield and maintenance states required by the
DNDC double-Monod rate. The rejected `k*N^2/(N+Kn)` form is not a total
Michaelis-Menten rate.

Gross denitrification uses a clipped linear bulk-pH response from pH 4 to a
plateau at pH 7.5. The shape is adapted from Blanc-Betes et al. (2021), which
used a pH 7.0 breakpoint. 7.5 is the top of the range Olaya-Abril et al.
(2021) tested and inside the pH 7.0 to 8.0 optimum band of Simek and Cooper
(2002); it is not a universal optimum. The N2:N2O pH
factor is Blanc-Betes Eq. 9 with no floor of its own; the published 0.1 minimum
applies to the complete N2:N2O product, as in their Eq. 6. The nitrate,
respiration and WFPS terms are the legacy TC-SMEW DayCent forms and change only
the product split, so the composite ratio is a hybrid rather than the published
Blanc-Betes ratio. N2 plus denitrification-derived N2O
therefore equals the gross nitrate sink.

Ammonia volatilization remains the legacy empirical term
`VOL = kv*B(31)`, with `kv = 8.6e-5 d^-1`. It has no direct pH multiplier or
explicit NH3/TAN equilibrium. pH affects it indirectly because the pH response
of nitrification changes the competing NH4 demand and the NH4 pool available
on later days. Treat `VOL` as an NH4-based legacy loss diagnostic, not a full
surface NH3 transfer model.

Nitrate leakage remains an inventory-limited mass-flow term driven by drainage.
It has no direct pH multiplier. pH affects leakage indirectly by changing
nitrification, denitrification and the nitrate pool available to drain.

Every dissolved pool of the biogeochemical layer, and the ERW's, leaches
as in T&C-BG in every model iteration: with the leakage out of the bottom
of the soil column over the water of the whole column (`Biogeo_Leakage`,
`Lk/VSUM`; the ERW's `f_Lw_bio`). Leaching them with the layer's own
drainage over its own water was tested and not adopted;
[`LEACHING_WATER_2026-10-08.md`](LEACHING_WATER_2026-10-08.md) records why.

The detailed equations, dimensions and citations are in
[`v2/POST_LEGACY_EQUATION_AUDIT.md`](v2/POST_LEGACY_EQUATION_AUDIT.md).

### Phosphorus from weathering

Mineral P released by weathering returns once per ERW window to the T&C
mineral-P state `P(:,43)`, in g P m^-2, through
`ERW_Return_Weathered_Phosphorus`. The 12-cell and vertical 1-D layouts share
that one transfer. The window budget is checked at run time, so a replayed or
inconsistent window raises an error instead of creating P.

No feedstock row used by the twelve runs contains a P-bearing mineral, so
basalt-derived P supply is exactly zero in the current campaign. The Blue Ridge
metabasalt mineralogy follows Lewis et al. (2021),
doi:10.1016/j.apgeochem.2021.105023, which did not detect apatite or any
K-bearing mineral in that rock. Bulk P and K that XRD did not resolve into a
named phase are not released, because release is computed from named mineral
stoichiometry. Do not add an apatite fraction without a feedstock measurement
that assigns P to a phase.

### Field topography

Lian et al. (2025, *WRR* 10.1029/2025WR040260) show with T&C-BG-2D that
topography sets most of the spatial variance of soil carbon and nutrient
pools. The 2-D field driver resolves the same elevation-driven lateral and
vertical redistribution, but with v0.9's own lateral scheme rather than
T&C-BG-2D's kinematic-wave routing: the two-point finite-volume
Buckingham-Darcy exchange across every shared face, with horizontal
conductivity `K aR` and the frozen-soil impedance `10^(-5 F_ice)` of the
vertical T&C solver, each step kept inside `[Ohy, Osat-Oice]` by adaptive
sub-stepping. Topography adds one term, gravity. Attach a DEM to any
multi-cell topology and the driving potential switches from matric head
to **total hydraulic head** `H = psi + z`
(`Darcy_Lateral_Topology`; the flux runs along the slope distance and
crosses the slope-normal layer thickness `Zlayer cos(beta)`, likewise for
lateral heat conduction). The elevation term is the **nodal** drop
`z(a) - z(b)`: the two-point flux must be a difference of cell-centred
potentials, or a uniform total head is no longer an equilibrium and the
drops no longer sum to zero around graph cycles, leaving a closed domain
with a permanent circulation. The sub-field graph's buffer strips span
the whole domain, so their centres sit laterally offset from the cells
they touch; the **slope** of such a face is measured over the true
centre-to-centre distance (`topo.edge_centre_distance_m`), while the
resistance path keeps the graph's own face-separation convention scaled
by `1/cos(beta)`, so a zero-relief DEM leaves every flux bit-identical
to a run with no DEM. The first-order two-point error on an offset face
is thus the one the flat model has always carried, and it is
equilibrium-preserving. The review caught an earlier attempt to project
the drop onto the face normal instead: it broke both invariants above
and made the overland weights disagree with the routing sweep, so the
run aborted on the mass-balance guard the moment runoff occurred.
Overland runoff is routed down the DEM by
Quinn (1991) multiple-flow-direction weights (`(tan beta)^1.1` times the
face width) and leaves only at the topographic outlet
(`Surface_Runoff_Routing`), and the direct solar beam is projected onto
each cell's slope and aspect with a sky-view diffuse correction. Without
a DEM every run is flat and bit-identical to before
(`test_2d_driver_equivalence` still passes at tolerance 0).

Overland routing uses a **runoff-runon connectivity** rule rather than a
fixed number of hops per step: a cell that produced runoff itself this
step is at infiltration capacity (Hortonian) or saturated (Dunne), so
run-on cannot infiltrate there and cascades straight through, in the same
step, in descending-elevation order; a cell that produced no runoff has
spare capacity, so run-on stops there and is offered to its next vertical
step as `q_runon` (a rate, mm h-1). A hillslope generating runoff
everywhere therefore delivers to the outlet within the hour, as overland
travel times demand, while runoff that reaches dry ground is allowed to
re-infiltrate. A dry outlet cell absorbs like any other dry cell: it is a
full cell of field, not a channel.

Solute rides that water rather than being credited on arrival.
Saturation-excess runoff carries the active-zone soil-water
concentration; infiltration-excess runoff carries the **ponded mixture**
of rain and whatever run-on has not yet infiltrated, held per cell in a
surface store. The share of that store leaving again is `Rh / WIS`,
`WIS` being the water that actually reached the surface this step
(throughfall, drip, melt, run-on) and the quantity `Rh` is defined
against in `HYDROLOGIC_UNIT`, capped at 1. Crediting delivered mass
straight to the soil pool would strand it in the first absorbing cell,
because the water it arrived on has not infiltrated yet and the part that
runs off again must take its share. The **ERW aqueous stores** (Na, Al,
residual anion, DIC and hence alkalinity) get the identical two-limb,
pond-store treatment (`ERW_Surface_Export`), the surface counterpart of
`ERW_Lateral_Transport`; without it water leaving the field would carry
Ca and Mg but no DIC, the exported solution would not be charge balanced,
and the CDR ledger would miss the overland pathway. Field-origin
provenance is carried through the router as passive tracers **and
through the pond store**, because routed water mixes at confluences, the
outlet is a polygon-clipped perimeter cell whose own overlap fraction
says nothing about where its water came from, and a stop-over on such a
cell must not re-label the load.
`Tests/test_surface_runoff_routing.m` (21 checks) covers the cascade, the
absorbing stop and dry outlet, confluence mixing, conservation and pits;
`test_2d_topographic_smoke` adds a wet pass that drives the whole path
end to end. Two deliberate assumptions remain: the domain perimeter is
closed to overland flow except at the outlet cells (a side boundary on a
cross-slope does not leak), and cast shadows (`Shadow_Effect`) are not
applied -- both negligible on a smooth field with 0.4-16 % slopes and a
431 m terrain autocorrelation length.

**Solar time convention.** `SetSunVariables` averages the sun over
`[HR-t_bef, HR+t_aft]`, and the run scripts pass `t_bef=-0.67,
t_aft=1.67`, i.e. `[HR+0.67, HR+1.67]`. The ERA5-Land forcing labels each
record with the *end* of its accumulation hour, so record `HR` covers
`[HR-1, HR]`. Checked against the bundled forcing, the hour-ending window
tracks the observed shortwave far better (correlation 0.88 against 0.69;
diurnal centroid 12.35 h against 10.73 h for an observed 12.43 h), and
the two windows differ by 25 degrees of azimuth at midday -- enough to
reverse the east/west aspect contrast the projection exists to resolve.
The per-cell **beam projection** therefore uses the forcing's own
convention (`TOPO_SOLAR_WINDOW`, default `[1 0]`) and reports the
correlation at startup, while the global `t_bef`/`t_aft` feeding the 1-D
energy balance is left untouched, since that is a separate calibration
question and changing it would move every existing result.

The vertical column stays strictly 1-D in 2-D mode: `SOIL_MOISTURES_RICH_COMP`
drives its internal kinematic lateral loss from `sin(atan(Slo_top))`, so a
per-column slope would double-count the FV graph exchange. The topographic
effect on vertical drainage, infiltration and leaching instead comes
through the coupled state - lateral redistribution makes footslope columns
wetter and shoulder columns drier - plus per-cell radiation and
temperature. Over Bruce Farm slopes (1-8 %) the `Asur`/`cos(alpha)`
correction that is lost this way is below 0.2 %.

```matlab
%% ingest the DEM tiles and preview the terrain attributes
addpath(fullfile('v2','interpolation'));
Run_Field_DEM('AdamWilbourne_Control_Fields.shp', ...
    fullfile('v2','Spatial inputs','Bruce Farm Virginia','DEM'), ...
    'FeatureIndices',[1 2 3]);

%% compute the process-grid cell size from the site's own length scales
Run_Field_Length_Scales( ...
    fullfile('v2','Spatial inputs','Bruce Farm Virginia','Whole farm', ...
    'AdamWilbourne_Control_Fields.shp'), ...
    fullfile('v2','Spatial inputs','Bruce Farm Virginia','DEM'), ...
    'FeatureIndices',[1 2 3]);
```

`Load_Field_DEM` mosaics the tiles, reprojects NAD83 -> WGS84 UTM through
`projinv`/`projfwd`, and when the requested spacing is coarser than the
1 m lidar it returns the **area mean** of the source cells under each
output cell (an r x r box filter) rather than a point sample, so slope,
curvature and TWI statistics are not aliased by microtopography. Cell
aspects are sampled from the grid as a unit vector (sin/cos), never as
raw degrees across the 0/360 wrap.

`Run_Field_Length_Scales` fits an empirical variogram to each measured
Bruce Farm soil attribute (shortest range = `L_soil`) and to the
detrended-elevation semivariogram of the DEM (`L_terr`), and recommends
`h = max(DEM spacing, min(L_soil, L_terr) / 3)` capped by the 95 %
per-cell coverage rule. That fixed 80 m "notebook recipe" variogram range
in `Fit_Ordinary_Kriging` is a default, not a Bruce Farm measurement, and
is not used for grid sizing.

**Grid convergence.** `Run_Topographic_Grid_Convergence` confirms a width
by convergence rather than by assertion. Over a ladder of widths (20,
10, 5, 3.5, 2.5, 1.5 and 1 m) the lidar DEM is area-averaged to each
width and the terrain statistics the flow code depends on are derived on
the raster, masked to the polygon, so the question is put to the terrain
and is not blocked by the coverage rule (which admits no grid coarser
than 2.33 m here; cell count, outlet count and the share of the polygon
that drains to an outlet are reported wherever a grid can be built). A
statistic counts as converged from a width only if every finer step in
the ladder also changes it by less than the tolerance (5 %), and the
verdict is reported per statistic: one quiet step between two noisy ones
is not convergence. On Control 2
(`v2/figures/Grid_Convergence_AdamWilbourne_Control_Fields_polygon_01.csv`
and `.png`) mean slope and mean TWI are converged from 10 m and relief
from 20 m; the spatial std of slope only from 1.5 m and the spatial std
of TWI not within the ladder. |Plan curvature| is tabulated but excluded
from the verdict: a second derivative of a surface whose sub-metre
roughness (furrows, lidar noise) it picks up, it grows without limit as
the width shrinks (0.012 to 0.127 m-1 from 20 m to 1 m) and has no
convergence width to find. So the drivers of the lateral physics, slope
and head, are resolved by the 2.33 m process grid with margin; what is
left is the micro-roughness that the area-averaging removes at the
process width (the 3-5 % per halving in the slope std) and the known
grid dependence of TWI itself, whose specific catchment area scales with
the flow-tube width (Sorensen & Seibert 2007). With `RedistributionDays`
> 0 the ladder also makes a short topographic 2-D run at each buildable
width and adds lateral water moved per m2 of the simulated grid per day
and the spatial coefficient of variation of end-of-run water over the
whole soil profile (`profile_water_cv`).

**Redistribution diagnostic.** `Plot_Topographic_Redistribution` takes a
topographic run and its flat twin on the same grid and draws, per
variable, the 2-D field, the 2-D-minus-flat map and that difference
against TWI, plus a gradient-method feature-importance table
(elevation, TWI, aspect, curvature). It was exercised on Control 2 at
10 m (109 columns, coverage rule relaxed to 60 %, Kantola case 5
physics and forcing on the Bruce geometry) over the wettest warm week
of the forcing, 24-30 April 2017 (139 mm), after a January window had
given frozen soil and no signal at all. The flat run has zero spatial
variance by construction; the topographic run's departures from it are
small but coherent on this low-relief field (slopes 1.4-7.8 %, relief
5.6 m):

| variable | 2-D minus flat | relation to TWI |
|---|---|---|
| profile water, end of week | -0.52 to +0.55 mm | r = 0.30, wetter at high TWI |
| profile water, run mean | -0.17 to +0.17 mm | +0.033 mm per TWI unit |
| NO3 pool, run mean | -0.0008 to +0.0013 g N m-2 | aspect carries 61 % of the sensitivity |
| Ca in solution, run mean | -0.004 to +0.005 g Ca m-2 | no TWI trend |

No surface runoff was produced by 139 mm on this column (its
infiltration capacity exceeds the hourly intensities and it drains
freely at the base), so the whole signal is subsurface lateral flow and
aspect-dependent radiation. The 2-D driver's `Domain.lateral` ledger
now carries `water_edge_abs_m3`, the daily sum over faces, layers and
hours of |Q_edge|: the buffer/main ledgers are zero by construction on
a whole-field polygon grid, and this is the quantity that says how much
water the topography moved between cells (it is also what the
convergence ladder's `RedistributionDays` metric reads). That is a statement about the Kantola
column on Bruce terrain, not about Bruce Farm: the per-cell texture,
the perched-water-table lower boundary and the site forcing are the
deferred sweeps, and they are exactly what would turn lateral flow on.
`v2/figures/Topographic_Redistribution_Control2_10m_wetweek.png` (and
`_15d_winter.png` for the null January case).

A topographic run passes `FieldDEM` (and optionally `DEMResolutionM`) to
the launcher; the DEM is fetched with `Load_Field_DEM` (mosaic + reproject
+ resample) and attached with `Attach_Topography_To_Topology`:

```matlab
Run_One_Subfield(1,'FieldShapefile',shape_file, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',1, ...
    'FieldDEM',fullfile('v2','Spatial inputs','Bruce Farm Virginia','DEM'), ...
    'DEMResolutionM',5)
```

The polygon-grid builder's 95 % per-cell coverage rule admits nothing
coarser than the calibrated width (2.33 m on Control 2), so a coarse
diagnostic grid has to relax it explicitly: `FieldMinCoverage` on
`Run_One_Subfield`, `MinCoverage` on `Run_Field_Cell_Spinup` and
`Krige_Cell_Agronomic_Values`, and the fifth argument of
`Write_Field_Agronomic_Templates` are passed through to
`Build_Shapefile_Grid_Setups` / `Build_Projected_Field_Topology` and
become part of the setup cache key. Lowering the rule is a diagnostic
choice only: the retained grid then covers less of the polygon and the
field ledgers are for that interior core, not for the field. Checkpoint
and result tags carry `__dem<res>` for a topographic run and
`__cov<frac>` for a relaxed rule, so a topographic and a flat run of the
same grid, or two coverage rules, never share a file.

The `v2/Spatial inputs/**/DEM/` GeoTIFFs are fetched, not versioned
(gitignored). For the Adam Wilbourne control fields, the USGS 3DEP 1 m
tiles `USGS_1M_17_x71y405..x72y407` from the
`VA_FEMA-NRCS_SouthCentral_2017_D17` project cover both control polygons.

### Parameterizing grid cells from kriged soil measurements

`Write_Field_Agronomic_Templates` writes every `baseline_*`/`measured_*`
column as `NaN` and every `use_measured_*` toggle `false` -- filling them
from the kriged soil surfaces is a separate, deliberate step so an
unreviewed estimate can never become a model input by accident:

```matlab
addpath(fullfile('v2','interpolation'),fullfile('v2','spinup'));
shp = fullfile('v2','Spatial inputs','Bruce Farm Virginia','Whole farm', ...
    'AdamWilbourne_Control_Fields.shp');
sample_file = Bruce_Farm_Input_File(fullfile('Whole farm','Data', ...
    'RJ_msft_yale_geochem - agronomic.xlsx'));

%% 1. krige pH / Organic Matter / CEC onto every cell -> a REVIEW table
[~,kriged_files] = Krige_Cell_Agronomic_Values(shp,sample_file,[], ...
    'FeatureIndices',[1 2],'Stage','BLC');

%% 2. the deliberate, threshold-gated merge into the template
Apply_Kriged_Cell_Values(template_file,kriged_files{1}, ...
    'MaxRelativeStd',0.35,'OutputFile',parameterized_file);
```

A cell is promoted from baseline to measured only when the kriged
estimate is not an extrapolation (inside the convex hull of the fitting
samples), its relative standard deviation clears `MaxRelativeStd`, and
the value is physically admissible (fractions in (0, 1), pH in (0, 14),
CEC positive: ordinary-kriging weights can be negative, and a K
saturation of 1-3 % can be predicted below zero with a small std);
every other cell keeps exactly what the template already had. A cell
whose toggle is already on carries a value someone entered
deliberately and is never overwritten unless `OverwriteMeasured` is
set, and attributes that share a toggle (sand and clay; the three
base-saturation fractions, which must also sum to less than one) are
written together or not at all, so a cell's texture or exchanger never
mixes two sources. The log table counts, per attribute, the cells
accepted, written, rejected by each gate, and kept as already measured.

**Base saturation.** The workbook's `%Ca Saturation`, `%Mg Saturation`
and `%K Saturation` columns (`Ca[2+] Base Sat (%)` and so on in the
catalogue export) are kriged with the same ensemble as pH, OM and CEC,
as the fractions `bs_Ca_fraction`, `bs_Mg_fraction` and
`bs_K_fraction`, and enter the template as `measured_bs_*_fraction`
under one shared toggle, `use_measured_base_saturation`, so a cell is
promoted only when all three pass the gate. There is no pedotransfer
baseline for the exchanger, so the template carries no `baseline_bs_*`
column and `Apply_Kriged_Cell_Values` needs only a `measured_` column
and a toggle to recognise an attribute; a cell with neither a
measurement nor a baseline keeps the existing behaviour, in which the
handoff derives the exchanger composition from the pore water and the
bulk pH alone. A site may instead give its representative exchanger as
a baseline (`baseline_bs_Ca_fraction`, `baseline_bs_Mg_fraction`,
`baseline_bs_K_fraction`: all three, finite in every cell, like the
other baselines; Hamden, below), which a cell's measurement replaces
where its toggle is on. Both seed the exchanger the same way, and the
cell's identity records which (`exch=...:baseline` or `:measured`).
Where either is used, `Resolve_Cell_Agronomic_Inputs` checks it (each fraction in (0, 1), the
three summing to less than one; the remainder is the exchangeable
acidity), reports the exchangeable Ca, Mg and K it implies in
cmolc kg-1 against the cell's resolved CEC, and
`Configure_Cell_Agronomic_Inputs` writes the cell's ERW exchanger seed
`f_CEC_in = [Ca Mg K Na Al H]` as `[bs_Ca bs_Mg bs_K 0 0 1-sum]`: the
acidity share starts as exchangeable H+, and Na and Al exchange start
at zero and follow the handoff equilibrium. It also switches that row
to v0.9's `f_CEC_inventory` convention: the workbook rows use the
`Kantola` convention, whose seed prescribes an 80 % base saturation in
`ERW_pre_processing` and would overwrite the measurement (the anion
closure is full-charge under both). A row set to `spinup_equilibrium`
(v1.0) takes a measured base saturation only. Without one, no
composition is prescribed: the seed keeps the spun-up soil solution (the
biogeochemistry's dissolved Ca, Mg and K, the rain's Na), and
`Invert_Bulk_pH` puts the exchanger in equilibrium with it, by the run's
own exchange law, at the pore-water pH that reproduces the measured bulk
pH. The base saturation then emerges from the spin-up, the CEC and the
pH (`test_kantola_seed_charge`, `test_cell_agronomic_inputs`). That seed is the right
entry point mechanistically: `Invert_Bulk_pH` treats a supplied
exchanger seed as a conserved exchangeable inventory that it combines
with the dissolved cations and re-partitions by the same Gaines-Thomas
relations the daily ERW solver enforces, at the pore-water pH that
reproduces the measured bulk pH, so the handoff state is a fixed point
of the model's own exchange law with the measured stock, not a
prescribed composition that the first day relaxes away from. The
spin-up saves that handoff exchanger and the transient starts from it
(`Load_Cell_Parameters`), and the T&C exchangeable pools B(61), B(69),
B(53) are re-initialised from it at the first ERW step. (An earlier
version set those three pools as spin-up observation targets instead;
the ERW step overwrites them and the T&C exchange kinetics relax them
within days, so that route never reached the transient exchanger and
was removed.) `Tests/test_base_saturation_targets.m` covers the
conversion, its validation, the template columns, the merge and the
seed. Kriged base saturation is as well calibrated as the other
attributes on Control 2 (leave-one-out rms standardised error 1.10,
1.17 and 1.30 for Ca, Mg and K) and Control 3 (0.85, 1.08, 1.32).

**Uncertainty.** Ten samples per control plot cannot resolve a
semivariogram (`Fit_Empirical_Variogram_Range` reports `resolved =
false`), so the variogram range is neither fitted nor assumed: it is
integrated out by Monte Carlo (`Krige_Attribute_Range_Ensemble`, 200
draws) over a log-uniform prior whose support the sample geometry itself
sets (`Kriging_Range_Prior_Bounds`: median nearest-neighbour spacing to
sample-cloud diameter -- 23-100 m and 29-159 m on the two plots). The
**sill** is treated as uncertain too, since the kriging variance is
proportional to a sill that ten samples also only estimate: its sampling
distribution is the scaled inverse-chi-square with `n-1` degrees of
freedom, and because scaling sill and nugget together leaves the kriging
weights (hence the prediction) untouched and scales the variance
linearly, its expectation is applied in closed form -- `E[Var | range]`
times `nu/(nu-2)` -- rather than drawn, which would only add sampling
noise. The review table's `<stem>_std` is the total by the law of total
variance, `E[Var | range, sill] + Var(E | range)`, with its components
beside it; `RangeDraws=0` recovers the classical single fit at the 80 m
notebook default, flagged `range_assumed = true`. The leave-one-out
check re-derives sill, nugget and prior inside every fold from that
fold's own samples, so the held-out value can never leak into the
variance it is scored against.

Calibration was measured rather than assumed, on 300 synthetic Gaussian
fields whose range was itself unknown (drawn 10-150 m), 10 samples each,
scored at in-hull points:

| estimator | rms standardised error | 95 % interval covers |
|---|---|---|
| fixed 80 m range | 1.60 | 80.5 % |
| Monte Carlo range only | 1.27 | 88.3 % |
| Monte Carlo range + closed-form sill inflation | **1.13** | **92.2 %** |

So the classical single fit was substantially overconfident, and the
ensemble removes most but **not all** of that: at 1.13 it is still
mildly optimistic, because `n-1` degrees of freedom ignores that
spatially correlated samples carry a smaller effective sample size.
Treat `<stem>_std` as a good lower bound on the uncertainty, not an
exact one. Every attribute is also leave-one-out cross-validated
(`Krige_Attribute_Cross_Validate`), each fold re-deriving its own sill
and prior so nothing leaks. On Control 2 the rms standardised error is
1.05 (pH), 1.14 (OM) and 1.16 (CEC) with the ensemble against 1.39,
1.56 and 1.55 for the fixed 80 m fit -- calibrated. On Control 3 it is
0.88 for pH but 1.61 (OM) and 1.97 (CEC), against 0.99, 1.95 and 2.33
classical: better, but still overconfident for the two attributes whose
ten samples include one far-out value, and the review table says so
(`rms_z` per attribute in the printed report). A single sample changing
the verdict this much is itself the finding: for OM and CEC on that plot
the kriged surface is a screened estimate, not a measurement. `Tests/test_kriging_uncertainty.m` (25 checks)
covers the decomposition, the prior, the sill inflation factor,
reproducibility, exactness and the calibration check against a synthetic
field of known range.

The workbook's control-plot samples use `BLC`/`R1C`/`R2C` stage labels
(not the treated-field `BL`/`BLP`/`R1`/`R2`), and
`AdamWilbourne_Control_Fields.shp` carries **three** features: Control 2
(feature 1), Control 1 (feature 2) and Control 3 (feature 3);
`Krige_Cell_Agronomic_Values` defaults to all of them. **Control 2 and
Control 3** each hold 10 `BLC` samples (`Field ID` `Buck` in the
workbook, `Buck Control` in the complete sample-catalogue export -- the
same points) and krige cleanly: 580 of 1630 and 466 of 1357 cells fall
inside their convex hulls. At `MaxRelativeStd = 0.35` every in-hull cell
passes for pH on both plots and for everything on Control 2; on Control
3 the honest std rejects 175 (OM) and 146 (CEC) of the 466 in-hull cells
where the fixed-range fit would have rejected only 21 and 16. That plot
holds one sample far above the rest (OM 7.4 %, CEC 7.3 cmolc kg-1
against 2.4-3.7 % and 2.6-4.6 elsewhere), which is exactly the case
where ten samples understate the field's variability. **Control 1** (feature 2, 6.5 km south -- the HWY 58 Hay West
control) has no baseline-stage sample inside it at all, confirmed
against `v2/figures/Bruce_Control_Interpolation_Timepoints.csv`, only
`R1C`/`R2C`. **Decision: do not substitute `R1C` (or any other stage) as
a baseline proxy for Control 1.** It runs on pedotransfer baseline
values only, `use_measured_* = false` for every cell, until a genuine
baseline-stage sample exists there.

The complete sample-catalogue export (`MSFT_Lithos_Yale_Complete_Dataset.csv`
in the same `Data` folder; `sme_lng`/`sme_lat`, `pH Value`,
`Organic Matter (%)`, `CEC Value`) is accepted directly
(`'SampleSheet',1`; coordinate and attribute columns are auto-detected)
and reproduces the workbook predictions to within 0.005 pH and
0.02 cmolc kg-1. Soil **texture** (`Sand Mass (%)`, `Clay Mass (%)`) is
kriged when the file carries it, but no control-plot sample does (the
nearest textured sample is 3.3 km away), so SSURGO/SoilGrids remain its
source there. The bulk-density survey (`Bulk_Density_Full_Field.xlsx`,
fields F1-F8, 0-15 cm) lies 4 km from the control plots with no site
inside them, so per-cell bulk density is likewise not available; the
pedotransfer value from `Soil_parameters` stands.

### Flux-tracking field spin-up (T&C-BG-2D, adapted for ERW)

The field spin-up follows the T&C-BG-2D initialisation of Lian et al.
(2025, *WRR* 61, e2025WR040260, section 2.3 and Fig. 3): a plot-scale
flux-tracking spin-up combined with a random forest, so the slow soil pools
carry the lateral redistribution a 2-D field imposes without a multi-century
2-D run. `Run_Field_Cell_Spinup(..., 'Method','flux_tracking')` (the
default; it needs `FieldDEM`) runs four stages:

1. **Class seeds** (Fig. 3a-b). One plot-scale spin-up per land cover x
   soil class: `Run_BG_Spinup` (T&C-only ring, the stand-alone
   biogeochemistry under the long-term mean forcing to steady state, the
   coupled ring to dynamic steady state). Classes are the categorical
   columns of the agronomic table (`ClassColumns`, by default whichever of
   `land_cover`, `vegetation`, `pft`, `soil_class`, `texture_class` it
   carries; none -> one class), and each is spun up on its most typical
   cell, the one nearest the class median of the resolved soil
   properties. Every cell of a class starts from its seed.
2. **Tracking** (Fig. 3c). The whole field runs in 2-D from those seeds
   (`Run_Field_Lateral_Tracking`, `MAIN_FRAME_BG_ERW_2D` with the terrain
   attached) over the spin-up's meteorological ring and its handoff
   extension, split exactly as Phase 3 splits its ring: one segment per
   cycle year, then the partial year up to the transient start, every
   column restarting from the previous segment's end state through the
   same hand-over. Meanwhile `Lateral_Track_Record` records, for a
   topographically stratified ~10 % subsample (`Select_Topographic_Subsample`,
   `SubsampleFraction` 0.10 as Lian et al.), what entered and left each
   tracked cell: the hourly Darcy change of every layer's water content,
   the hourly overland run-on, and each day's exchange of every mobile
   pool as one exact affine map `T_end = retain .* T_start + input`,
   composed from the hourly upwind steps (subsurface transport,
   saturation-excess runoff and re-infiltrating run-on). `retain` is the
   share of the cell's own pool that stayed, `input` what arrived from its
   neighbours, so when the map is replayed the cell's export stays
   proportional to its own, evolving pool while what arrives from
   upslope is imposed as recorded: the upwind scheme with the
   neighbours' concentrations frozen at their tracked values.
3. **Replays** (Fig. 3d). Each tracked cell's own plot-scale spin-up with
   that exchange replayed in every phase (`LATERAL_REPLAY_FILE`): the
   water in Phase 1, the record's mean daily map in the stand-alone
   Phase 2, water and solute hour by hour and day by day in the Phase-3
   ring and the handoff extension, through inert hooks in the v0.9
   10-minute driver (`Lateral_Replay_Prepare`, `_Water`, `_Solute`), and
   with the topographic shortwave the cell received in the field run
   (`Topographic_Radiation_Projection`, one definition with the 2-D
   driver). The lateral gains and losses enter the six Phase-3 nutrient
   budgets as inputs and outputs, so they still close. The replay starts
   from the cell's own tracking-run state at the end of the ring and runs
   until Phase 3's ring-convergence test passes, the dynamic steady state
   Lian et al. reach with 1,000 years.
4. **Reconstruction** (Fig. 3e-f). A random forest
   (`Train_Pool_Random_Forest`, k-fold cross-validated; Lian et al. report
   < 10 %) maps terrain (elevation, slope, aspect, curvatures, TWI) and
   soil (sand, clay, OM, bulk pH, CEC) onto the replayed cells' conserved
   pool totals (POC, MAOC, DOC, microbial C, total N, NH4, NO3, total P),
   and every other cell gets its prediction mapped onto the 79-pool state
   of its nearest replayed cell, keeping that cell's partition and
   stoichiometry (`Initialize_BG_State_From_Observations`), with its own
   soil moisture from the 2-D run at the transient start
   (`Reconstruct_Cell_Pools_RF`). As in Lian et al. that is the cell's
   initial condition; no further spin-up runs.

**Adapted for ERW.** The ERW chemistry is off in every stage, as in every
spin-up phase of v0.9: the exchanger and the pore water are never spun
up, so neither the tracking nor the forest touches them. Every cell,
replayed or reconstructed, is handed on through the transient's own seed
path: the observation constraint at the transient start, then the ERW
seed built from the cell's own bulk pH, CEC and measured base saturation
(`Invert_Bulk_pH`, the `f_CEC_inventory` convention) or, under
`spinup_equilibrium` without one, from its spun-up soil solution, gated by
`ERW_Handoff_Readiness` (charge balance, pore and bulk pH, base
saturation). A reconstructed cell runs that gate alone
(`Run_Automatic_Site_Spinup` with `SPINUP_HANDOFF_ONLY`), keeping the
convergence record of the replayed cell it was built from
(`reconstruction_template_cell_id`, `reconstruction_rf_cv_rel_error`). The
forest predicts biological and nutrient pools only, never the exchanger,
pore-water pH, Al, DIC or the residual anion, so electroneutrality holds
by construction. What the tracking does carry for ERW is the lateral
redistribution of the solution K, Ca, Mg and Si pools and of NO3 and NH4,
the charged species the seed's anion budget closes against. A cell whose
reconstruction fails falls back to its own cold spin-up (flagged
`fallback`); one whose fallback fails too is flagged `failed` and stops
the assembly.

`Method 'per_cell'` remains the 1-D reference: every cell's own spin-up
with no lateral exchange.

**Warm start.** `'WarmStartStates'` (on the cluster
`TC_WARM_START_STATES`) re-spins a field from an earlier spin-up of the
same grid, its `CellBGStates.mat`, for when site observations replace the
provisional soil inputs the first spin-up used. Each earlier state is
matched by cell id and refused if it lies elsewhere (a cell's identity
records its centre), so states of another grid cannot be used, and the
warm run must write to another folder than the run it starts from. A warm
cell runs Phase 3 only (`Run_BG_Spinup`, `SPINUP_WARM_START`): Phases 1
and 2 exist to bring a cold soil near equilibrium and are skipped, the
ring converges under the new inputs from the earlier converged state,
and the observation constraint and the ERW seed are built at the handoff
from the new inputs, as always. Under `per_cell` each cell does this
alone; under `flux_tracking` no class seed is spun up, the tracking run
starts from every cell's own earlier state, each tracked cell replays its
new record in Phase 3, and the forest reconstructs the rest as before (a
cell whose reconstruction fails still falls back to a cold start). The
earlier state's provenance is part of the warm run's, so neither resumes
the other's checkpoints. `test_field_spinup_warm_start` (14 checks)
covers the refusals, the per-cell seed and provenance, and the
flux-tracking stages through the forest; the measured cost is under
*Cost* below.

`test_lateral_flux_tracking` checks the machinery: the radiation
projection, the exact composition of the hourly exchange into the daily
map (both directions of a face split layer by layer), the replay hooks'
coverage rules and exact ledgers, and, decisively, that a cell tracked in
a small topographic 2-D run and then run as a plot-scale column replaying
its record follows the 2-D cell far more closely than the same column
with the exchange left out. `test_spinup_handoff_only` runs a
reconstructed cell's handoff on Control 2 cell 1: no spin-up phase runs,
the seed built from the cell's own inputs passes the gate, and the pools,
soil moisture, convergence record and provenance are kept.


**Spin-up forcing.** The automatic spin-up cycles complete
meteorological years until its phases converge, then runs the handoff
window. It takes that forcing from the row's `meteo_file`: either a
daily station file (`Load_Daily_Meteo`, disaggregated to hourly) or,
through `Load_Spinup_Meteo`, an hourly T&C forcing `.mat` such as
`Data_Kantola.mat` or `Data_Bruce_Field_Spinup_2014_2023.mat`, used as
it is with its own diurnal cycle; the file's latitude, longitude and
clock offset are checked against the row's and a mismatch is reported.
The workbook's Kantola rows carry none of the spin-up columns
(`forcing_file`, `meteo_file`, `meteo_start`, `meteo_end`,
`transient_start`, `transient_end`), so
`Run_Field_Cell_Spinup(..., 'MODParamColumns', struct(...))` stamps them
onto every cell's parameter row; the spin-up accepts an override row
that carries extra optional columns, and both sides read the workbook
through `Read_Parameter_Workbook`, so an override is typed exactly like
the workbook row it replaces. `'ERWParamColumns'` does the same for the
ERW row, applied before each cell's own soil inputs, so a cell's CEC and
exchanger seed still win. Both refuse a column nothing reads and one the
cells set from the agronomic file (`PH`, `Psan`, `Pcla`, `Porg`,
`CEC_tot_in`, and the cover when the file carries it), and a row that
seeds the exchanger from base saturation (`f_CEC_inventory`) is refused
at once if any cell has none, rather than at its handoff hours later. On
the cluster both come from one-row CSVs (`TC_MOD_COLUMNS_FILE`,
`TC_ERW_COLUMNS_FILE`, read by `Read_Parameter_Columns`): a site's own
settings over the case rows it borrows. Three optional MOD columns serve
a site with records of its own: `co2_file` names a CO2 record in
`Inputs` of the form of v0.9's `Ca_Data.mat` (`Load_CO2_Record`,
interpolated linearly, refused outside its span); `pft_override` names a
vegetation script run after the base parameter file, in the 1-D and 2-D
drivers alike (the pattern of `PFT_Prairie_Margin`); and a non-empty
`rain_chem_file` hands the record `Rain_Chem_<Location>.mat`
(`Load_Rain_Chemistry`) to the ERW chemistry of both drivers: its rain
pH and its Na and strong-anion (SO4 + Cl) deposition. The Na seed is the
rain steady state, the Na flux over the long-term water loss of the ERW
layer (the transpiration drawn from it plus the leakage), which the
spin-up records for each ring year (`ERW_Mean_Water_Loss`) and hands off
with the state as the ring's mean (`seed_water_loss_l_m2_d`); a record
without it is refused. Two latent v0.9 flaws had kept the seed from that
steady state: the record's scalar flux was padded with zeros over the
forcing series, and the drivers seed on their second day, when the
series hold one day of water fluxes and zeros after it
(`test_rain_chemistry_seed`, `test_kantola_seed_charge`). N, K, Ca and
Mg enter through the biogeochemistry's deposition terms only, so each
element enters once. They are constant daily rates unless the vegetation
script gives `B_IO.Dep_events`, a year-by-year event record
(`Build_Site_Deposition`): then they arrive wet with each day's
precipitation (the year's precipitation-weighted mean concentration
times the day's rain or snow) and dry every day (the year's dry
deposition spread over its days), `Deposition_Day`. Both frames hand
`BIOGEO_UNIT` the day's rates and record them (`Dep_d`), Phase 2 takes
the ring's mean, and the Phase 3 budgets take what the frame applied;
every other site keeps its constant rates (`test_event_deposition`). The
record can also carry the ERW chemistry's Na and strong-acid anions (SO4
+ Cl; nitrate is N deposition), wet and dry the same way (`erw_*`
fields, `Deposition_Day`'s third output). Its long-term mean then seeds
the ERW pools in place of a rain record's wet-only fluxes, and with
`keyword_add` = 0 each day's own deposition is the ERW's Na and
strong-anion supply (`Dep_erw_d`), so deposition is the only source and
nothing returns what leaches. With `keyword_add` = 1 the steady closure
replaces it, as in v0.9. A vegetation script for a stand that is never
fertilised must also empty `B_IO.FertN_years`: the 10-minute driver
refills `FertN` with 202 kg N ha-1 in every year the borrowed maize-soy
row lists there, two years in three.
`Tests/test_spinup_meteo_loader.m` (10
checks) covers the slicing, complete-year detection, the cycle-year
validator and the refusals. The topographically stratified subsample
(`Select_Topographic_Subsample`) uses three quantile bins per terrain
feature and coarsens to two when three leave more non-empty strata than
the requested sample, so a small diagnostic grid is not made to spin up
two thirds of its cells.

Running the imported spin-up for the first time in this repository
surfaced four latent breaks, all fixed. Phase 1 (T&C only) now switches
ERW off together with biogeochemistry, because the ERW chemistry is
seeded from the biogeochemical pools that Phase 1 does not carry. The
10-minute driver records the daily biogeochemical environment (soil
temperature, saturation, matric potential, water sums, WFPS and the
anaerobic fraction of the active layer) that Phase 1 averages into the
mean forcing of the stand-alone Phase 2, and no longer clears the
field-capacity constant Phase 2 reads. A cold-start emergence had
tripped strict non-negativity guards on the plant nutrient reserve
(-4e-11 g N m-2); v0.9 package u restored the legacy T&C reserve, in
which a negative balance stands for diluted tissue, so that state is no
longer an error. And the spin-up and the field launcher read the
parameter workbook through one function, so a per-cell override row is
typed like the row it replaces.
Two more were interface drift between the two code lineages: the
stand-alone Phase 2 called `BIOGEO_UNIT` with the Model-repo signature
(it now passes the package-u tail, the Phase-1 mean anaerobic fraction
and WFPS of the active layer with the default exchange scheme and
`KNmaxd` that Phase 3 also uses, and reads v1.0's `P_inputs` output,
fertiliser, deposition and uplift phosphorus), and the Phase-3 nutrient
budgets read per-crown ledgers the Model-repo driver keeps, so the
10-minute driver now records the per-crown harvest reserve export and
planting input (`AgrHarNut_H/L`, `PlantingNut_H/L`, `PlantingP_H/L`)
and the daily phosphorus inputs, while the nitrogen gas loss is summed
from v1.0's own partition (volatilisation, N2 from denitrification and
total N2O; there is no NO pathway) and the fertiliser term is what
`BIOGEO_UNIT` added that day (`FertN_d`).

Running it on the v0.9 package-u driver surfaced the rest, all fixed. The
v0.9 driver gates its ERW chemistry on `OPT_ERW`, not on the Model
repo's `OPT_pH`, so the documented no-ERW Phase 3 had been solving the
daily ERW chemistry every year and reseeding the cation pools from the
workbook at each year boundary (+23 g Ca m-2 on the test cell, outside
every budget) at about 70 min per year; Phase 3 now clears `OPT_ERW`,
costs about a minute per year, and its fixed cation pools start on the
equilibrium Phase 2 reached under the same exchange model. Each Phase-3
year restarts from the previous year's end state, so a year boundary
moves no mass and the converged ring is a periodic orbit of the model;
the last year's mean is only the ring-convergence metric. Two pools are
not held to that orbit (`P3_UNGATED_POOLS`): occluded P (46), which has no
return flux, and secondary P (45), which turns over in 37 years toward
mineral P and some 82,000 years toward occluded P while a field's P
surplus accumulates in it ring after ring, so it cannot settle on a
spin-up's timescale. Phase 2 still converges secondary P under the mean
forcing, cheaply, and mineral P (43) stays gated in Phase 3. The final
handoff diagnostic waited for a readiness record the v0.9 driver never
writes (and ran the whole transient window doing so): it now runs the
two forcing days the ERW seed reads, with ERW on, and scores the seed the
driver itself built (`ERW_Handoff_Diagnostic`, `ERW_Handoff_Readiness`).
Its charge balance is taken over the species the seed holds: v0.9
speciates the pore water with activity-conditioned constants, and
re-speciating its Al with the thermodynamic constants had reported the
activity correction as an imbalance growing with dissolved Al (0.7
umolc/L at 25 umol/L on Control 2 cell 1, 75 umolc/L at 289 umol/L on
cell 17), failing every acidic Bruce Farm cell at the end of its spin-up.
The respiration and water-loss hand-offs are gone, because the v0.9
drivers seed the CO2 boundary and the Na rain seed from the run itself;
plant reserves start from the parameter file, as in v0.9. Phase 1 runs
the driver with biogeochemistry off, so the fertiliser schedule (which
lives in the biogeochemistry's inputs) is now read only with it on, and
the irrigation series the driver creates once is cleared between years
(a 2014 series broke the 2012 leap year). Phase 1 also refuses a ring in
which the vegetation never establishes (peak LAI below 0.5 in any cycle
year): before package u restored the legacy plant reserve the test cell's
crop never grew (annual-mean LAI 0.003; now 0.34 and 0.26 for the two
crowns). Every phase checkpoint and the final state carry a provenance
key (`Spinup_Provenance`: the cell, its soil inputs, the workbook case,
the MOD and ERW rows with their forcing, meteo and transient windows, and
the seed), so a task resumes only its own work and reuses only a state
of the same run, and a state records the instant it hands off at, which
`Run_One_Subfield` checks against its own first hour.

The first spin-up of a site with its own vegetation (Hamden: a mature
forest over grass) surfaced four more, all fixed; none changes what the
Kantola rows compute. Phase 2 ran only the base parameter file, not the
row's vegetation override that the driver runs in Phases 1 and 3, so it
equilibrated with the borrowed maize-soy roots, deposition and fertiliser
years; it now runs the override as the driver does. The Phase-3 nutrient
budgets took the plants' starting stock from row 1, the parameter file's
initial state, whose tissue nutrients (`TNIT` ... `TSIL`) the driver
never fills: harmless for a crop sown into bare soil, but a stand that
starts with wood (Hamden's 11 kg C m-2) gained its whole nutrient stock,
72 g N m-2, on the second day. The budgets now start from that stock
computed as the driver computes it (`Plant_Tissue_Nutrients`: each live
pool's carbon over its carbon:nutrient ratio, plus litter in transit).
And the driver hands each day's litter to the soil the next day
(`BIOGEO_UNIT` takes `ISOIL(j-1)`), so the litter of a loop's last day is
in neither the soil's nor the plants' end stock; the soil storage now
counts it (0.021 g N m-2, one December day's litter, on Hamden's first
loop, where the plant side then closed to 1e-14 g m-2 and the whole
budget to 3e-13). Like the year's plant growth, it is not carried into
the next loop, which restarts the vegetation, its reserves included, from
the parameter file and carries only the soil and its moisture. A budget
that fails reports every term. And the site's ERW row named its soil
"fine sandy loam", the SSURGO texture name, which the ERW chemistry (it
tabulates its exchange selectivities by USDA texture class) does not know,
so the final handoff check stopped; the row now gives Cheshire's class,
"sandy loam", and `test_hamden_field_configuration` checks the class
against the chemistry's own table.

A perennial stand -- trees, a sward -- lives from one year to the next. A
vegetation file declares one (`PFT_perennial = true`, then
`Perennial_Vegetation_Apply`). Its stand then starts only a cold start:
from Phase 3 on, every year, the handoff extension, the final handoff
check, the transient, a warm start and each segment of the field's 2-D
tracking start from the previous run's end state (`Veg_out`,
`Perennial_Vegetation_Capture`), carried through every checkpoint and the
state handed off. The 2-D driver hands each cell's stand to its vegetation
file (`Cell_BG_State_Context`). The state is every vegetation variable the
drivers read from the previous day (`Perennial_Vegetation_States`):
- the carbon pools, phenology and leaf age;
- the memories of growth, temperature, stress, uptake and light;
- the nutrient reserves and uptake suppression;
- the canopy the first day's hours see;
- the litter and exudates the soil takes the next day.

So a year boundary moves no mass, and a log line reports the stand each
year.

Phase 1 is excluded. It has no biogeochemistry, so the plants take up no
nutrients there, yet growth still draws on their reserves: carried
through it, Hamden's grass starved within a year and its trees within
four. Each Phase-1 year therefore starts from the file, as in v0.9.
Crops declare nothing and are unchanged (`test_perennial_vegetation`).

Restarting Hamden's forest every year had held it at the leaf mass its
initial fine roots allow. T&C stops leaf growth, and the spring flush
from the reserve, at LtR times the fine-root carbon (`VEGETATION_DYNAMIC`),
and the legacy Harvard Forest LtR is 0.3: peak LAI 2.5, every year.

Carrying the stand exposed T&C's sink limitation, which the v0.9
plot-scale driver switches on for every site (`OPT_EnvLimitGrowth = 1`;
T&C's own default, and the legacy T&C-BG-2D forest runs, leave it off).
With it on, NPP above the environmental growth limit goes to the
carbohydrate reserve past its cap of 2/3 of the sapwood. Carried from
year to year, the reserve grew from 0.8 to 2 kg C m-2 while the canopy
shrank to LAI 1.9. Hamden's vegetation file therefore switches it off,
and the 2-D driver now honours a cell's own setting (an inert extension
for v0.9 parameter files, which set none).

Carried, with the sink limitation off, the trees grow into a closed
canopy: Hamden's class seed reached peak LAI 4.9-5.1 within six coupled
years, near Harvard Forest's measured 4.5-5.

The site table's Harvard Forest row gives no wood turnover (`Wm` 0), a
value from water-carbon runs, where heartwood is inert. In a spin-up of
decades, heartwood then grows without limit and buries its sapwood's
nutrients. Hamden's doubled to 15.8 kg C m-2 in twenty coupled years,
burying 2 g N m-2 a year at T&C's heartwood C:N of 193 (three times the
deposition), and its base cations with it.

The legacy coupled forest runs of T&C-BG-2D turn heartwood over in 37.6
and 40.1 years, returning dead wood to the soil as woody litter.
Hamden's trees now do the same, at 40 years (kept by decision);
heartwood levels off where conversion from sapwood equals mortality.

With the turnover, the class seed's canopy holds through twenty-five
coupled years, cycling with the four-year forcing ring: peak LAI 4.8-5.0,
sapwood 820-950 gC m-2.

Without it, the same run held for fifteen years and then thinned, to LAI
1.4 by the forty-second, as the nitrogen buried in heartwood drained the
trees' reserve and the soil. By then the shrunken stand's reserves
sufficed again, so its uptake had switched off. Each year's stand is
logged (`Perennial_Stand_Summary`).

The grass started from seed, the state a sown crop starts from. Its
leaves were then capped at LtR times 20 gC m-2 of roots, LAI 0.36, and it
died out even in the open within two years. It now starts as an
established sward instead: roots of 140 gC m-2, which the legacy LtR and
specific leaf area turn into a mown lawn's peak LAI of about 2.5, and a
reserve at its cap. A 90 %-grass cell's lawn then thins over a decade to
a peak LAI of about 1.2 and holds there; it is unfertilised and its
clippings are removed each September. Beneath the closed canopy, the
sward is shaded out.

**Cost, measured on one Control 2 cell with the Kantola case (six cycle
years 2009-2014, a 2016-01-01 handoff after the 2015 extension), single
column, on this workstation while several other MATLAB jobs shared it:**
Phase 1 (T&C only) about 80 s per forcing year, 8 min for the ring;
Phase 2 (stand-alone biogeochemistry under the mean forcing) converged
after 787 mean-forcing years in 7 min; Phase 3, the coupled ERW-free
loop, about 1.3 min per forcing year, with its convergence test taken over
complete six-year rings. On the test cell (cell 17) the ring drift fell
from 0.67 to 0.015 per ring over the first 126 coupled years, an
e-folding time of about 60 years, and then flattened: at year 234 the
dynamic and slow drifts were still 0.011 and 0.012 against the 0.01
limit, held up by the inorganic phosphorus pools (pool 43, mineral P, and
pool 45, secondary P, which feeds the occluded pool), after about 3.5
hours of Phase 3. With secondary P out of the Phase-3 gate (above) the
resumed run converged at year 306 and passed its handoff gate, about 5.4
hours of Phase 3 in all. A plot-scale spin-up is therefore a job of many
hours, its length set by that slow phosphorus exchange; a cluster task
that reaches its time limit resumes from its checkpoint when resubmitted.
A warm start (below) is far cheaper: the same cell, its bulk pH, CEC and
base saturation revised (5.52 to 5.80, 3.47 to 5.0 cmolc kg-1, Ca 0.53 to
0.60), re-converged from its converged state in 24 coupled years (17 min)
and handed off at the new pH. In the flux-tracking spin-up that
cost is paid only by the class seeds and by the replays of the tracked
~10 % of cells, all array tasks; the 2-D tracking run costs about one
ERW-free column-year per cell per segment, because the 2-D driver steps
its columns one after another; and every other cell costs only its
two-day ERW handoff. The pipeline up to the spin-up call (templates,
kriging, gated merge, resolved per-cell inputs, topographic subsample,
spin-up configuration) was exercised on the real Control 2 grid (54 cells
at 15 m under a 50 % coverage rule; 14 in-hull cells took every kriged
attribute including base saturation) and is covered by the tests. The
flux-tracking stages were exercised on this workstation on a coarse
Control 2 grid: the 2-D tracking run with its records, a replay spin-up
of one cell with a lateral record through Phases 1 to 3 with its budgets
closing, and the handoff-only gate. The field spin-up itself is a cluster
job: `v2/cluster/README.md` describes the SLURM chain on Bouchet (plan,
class seeds, tracking, replays, reconstruction, assembly) and the
validation experiment that scores the random forest against a run in
which every cell is replayed.

### Field transient

`Run_Field_Transient` (`v2/run`) takes a field's spun-up states forward
over the site's transient window: every cell of the polygon in 2-D
(`MAIN_FRAME_BG_ERW_2D`, with the lateral water, solute and heat exchange
over the DEM's terrain that the tracking run uses), with the ERW chemistry
on and no rock. For Hamden it is the run the study asks for: bulk soil pH
under atmospheric deposition alone, 2020-01-01 to 2026-09-29.

- **Cells.** Each cell is configured as its spin-up configured it:
  `Run_Field_Cell_Spinup` with `DryRun` builds its rows (the workbook
  case, the site's MOD and ERW columns, the cell's own soil and cover),
  and the grid and terrain are built the same way. Each starts from the
  state it handed off at the window's first hour, with the exchanger its
  handoff seeded; `Load_Field_Spinup_States` refuses a state handed off at
  another instant or for other inputs. The forcing, CO2 and rain record
  are the ones the handoff check read (the rows' `meteo_file`, `co2_file`
  and `rain_chem_file`). A row with rock is refused.
- **Segments.** The 2-D driver keeps every column's hourly state, and each
  ERW column a ten-minute record, for the whole window it is given: about
  11 kB per column-hour, so some 570 GB for Hamden's 936 columns over
  2,464 days (a 30-day segment holds about 7.4 GB). The window
  therefore runs as segments of `SegmentDays` (30), each going on from the
  one before exactly as one run would (`Field_Transient_Carry`,
  `Field_Transient_Restore`):
  - a segment starts at 00:00 of the previous segment's last but one day,
    so the drivers' calendar days and daily stamps fall where one run puts
    them;
  - those two days' 48 hours become its first rows, so every hourly
    look-back reads what one run would read: the daily step's 24 hours,
    and the ERW day's water sums, which also read the hour before the day;
  - the last three daily rows become its first three, the columns'
    carried state (cover, rooting, pH, ponded run-on and its solute, the
    export ledgers) is handed on, and the time loop resumes at the hour
    after the carried days;
  - each ERW column is laid on the segment's ten-minute and daily axes,
    its last two days' rows first and every later row at the value its
    array held at the end: a segment that another follows runs no final
    ERW flush, so those rows are the array's preallocation. An array the
    ERW days grow row by row keeps its rows to the anchor, and the CO2
    record the first ERW day lays on its whole axis is laid on the new one;
  - every segment takes the whole window's longest day. The parameter
    files find it over the days a run is given (`Lmax_day = max(L_day)`),
    and it scales leaf photosynthesis with day length (`RELATIVE_PC`), so
    a winter segment would otherwise take its own short days for the
    longest.

  Four inert hooks in the driver serve this: the restore, the loop's
  first step, the window's longest day, and no final flush (nor
  particle-history assembly) in a segment another follows. Before it
  replaces the longest day, the driver checks the one the parameter files
  find over the segment against the runner's own finding. Segments are
  therefore at least two days long,
  and a column with particle classes (rock) is refused, since its
  application schedule counts from the start of the run. The ERW seed's
  rain Na had been the mean of the record's one long-term flux over the
  forcing's copies of it, whose last-bit rounding depended on the run's
  length; it is now that flux (`ERW_MOD_PARAM_excel`), so a run and its
  first segment seed alike. `test_field_transient` runs four Hamden cells
  for five days in one segment, in segments of 2, 2 and 1 days and in
  segments of 3 and 2: the daily records and the whole end state agree bit
  for bit.
- **CEC scenarios.** `'CEC'` (on the cluster `TC_TRANSIENT_CEC`, which
  takes a list) gives every cell one CEC at pH 7, the field-scale form of
  `Run_Plot_CEC_Scenarios`: the spin-up does not depend on the CEC, so the
  scenarios start from the one field spin-up, each with its handoff's
  exchanger fractions at its own capacity.
- **Checkpoint and resume.** A finished segment is written to
  `segment_<k>.mat` (its records) and `carry_<k>.mat` (the state the next
  starts from) under a key that chains the run's identity (the cells and
  their rows, the spun-up states, the terrain, the window, the segment
  plan and each segment's forcing) through every earlier segment, and a
  rerun reuses it. The code that wrote a segment is recorded; a segment of
  other code is refused unless `ForceResume` (`TC_CKPT_FORCE_RESUME=1`)
  splices them deliberately. Inside a segment the driver's rolling
  checkpoint resumes (`CheckpointMinutes`, 60). A job that runs out of
  time is therefore simply run again.
- **Output.** `Field_Transient_<from>_<to>.mat` (`result`) holds one record
  per ERW day for every cell, and the field's area-weighted means (pH by
  mean hydrogen activity), which are also written as a CSV. The model day
  runs from 01:00 to 01:00 local standard time, the drivers' daily step.
  Record 1 is the start (the handoff state as the ERW chemistry seeds it)
  and the last is the end of the window. Each record holds the state then:
  bulk and pore-water pH, base saturation, the exchanger's charge
  fractions (Ca, Mg, K, Na, Al, H) and its CEC, dissolved Ca, Mg, K, Na,
  free and total Al, the residual strong anion, DIC and the water content
  of the biogeochemical layer. It also holds the fluxes over the 24 hours
  before it: the soil biogeochemistry's leaching of Ca, Mg, K, NO3 and
  NH4, the ERW chemistry's leaching of Na, Al, strong anion and DIC,
  drainage and precipitation. `result.units` names the units.
- **Cost.** Measured on the workstation with other runs going:
  - The daily ERW chemistry is 98 % of a four-cell segment's time, about
    15–30 s per column-day, mostly the Appelo charge balance and the
    bulk-pH inversion. The columns run their ERW days in parallel on
    `ParallelWorkers` thread workers.
  - The hourly steps run one column after another, then the lateral
    exchange: about 12–25 ms per column-hour, and 7.3 s per hour for the
    lateral exchange of all 936 cells.
  - Each segment starts by loading 936 columns' parameter files (about
    1.3 s each).

  The whole Hamden window is therefore about 10,000–19,000 core-hours of
  ERW chemistry and, on the workstation, two to three weeks of serial
  hourly stepping that no number of workers shortens; on Bouchet's faster
  cores perhaps one (the spin-up's tracking run steps about 8 ms per
  cell-hour there). With 32 workers on Bouchet that is roughly two to four
  weeks in all; measure it with a two-day run of the whole field first
  (`TC_TRANSIENT_NDAYS=2`). Running the hourly column loop on workers too
  (`v2/cluster/README.md`) would cut the serial part.
- **On Bouchet.** `v2/cluster/submit_field_transient.sh` submits it with the
  spin-up's own environment file (`v2/cluster/README.md`, *The field
  transient*).

`Run_One_Subfield` is unchanged: it stays the Energy Farm's twelve-run
matrix. It is bound to the matrix (run ids, margin layouts, the rock
schedule), to the workbook's Kantola rows and window and to `Ca_Data.mat`,
and reads `Data_<Location>.mat` rather than the meteo file a spin-up's
handoff read. It also holds its whole window in one driver call, which the
twelve-cell Energy Farm runs afford and a 936-cell field does not. The
field transient shares its building blocks instead: the agronomic and
state loaders, the rain record, the grid, the DEM and the driver.

### Hamden (145 Blue Trail, CT)

A 2.34 ha residential-forest parcel in Hamden, Connecticut, set up for a
field spin-up on Bouchet with no rock: the runs track how differences in
CEC change bulk soil pH under atmospheric deposition alone, under
temperate deciduous trees over grass mown once a year, with no liming or
fertiliser. The source data, their provenance and the decisions are in
`v2/T&C_Code/Inputs/Hamden/README.md` (a documentation copy of every
input is kept outside the repository); the site polygon, the land-cover
window and, once downloaded, the 1 m DEM tile are in
`v2/Spatial inputs/Hamden Blue Trail`.

- **Grid**: 5 m squares in UTM 18N under the majority rule (above): 936
  squares, turned 19.31 degrees onto the parcel's 115 m and 71 m sides,
  taking in 274 m2 of outside land and leaving out 273 m2 of the parcel,
  2.3 % misrepresented (the coverage rule's 1,224 squares of 4.48 m took
  in 1,179 m2, 5.0 %). The width is the terrain's: 5 m is the coarsest
  width from which the mean slope, its spread, the wetness index's spread
  and the relief change by less than 5 % at every refinement down to 1 m
  (`Run_Topographic_Grid_Convergence`; the slope's spread is the last to
  settle); the mean wetness index never settles, falling as the 1 m
  micro-relief is resolved.
- **Land cover per cell** (`Cell_Land_Cover`): either from a land-cover
  raster -- each cell's footprint (`Cell_Footprints`, the square clipped to
  the field) intersected exactly with every pixel it touches, the pixels
  read through a window rather than the whole map, and each class split
  by `NLCD_Cover_Class_Table` among tree crown, open grass, urban, water,
  bare and rock -- or set by hand (`'manual'`, a table by cell id or one
  setting for every cell). The six `cover_*` columns of the agronomic
  file (all or none, each row summing to one) become each cell's
  `Ccrown = [tree grass]`, `Cbare`, `Curb`, `Cwat` and `Crock`; crown 1 is
  the trees with grass beneath them, crown 2 the open grass. NLCD 2020
  gives the parcel 80.3 % deciduous forest, 17.7 % woody wetland and
  2.0 % developed open space, so tree 0.980, grass 0.018 and urban 0.002;
  the kept squares' in-parcel parts hold 0.982, 0.016 and 0.002 (the 1 %
  of the parcel left in dropped boundary squares is grassier than the
  rest).
- **Soil per cell** (`Cell_SSURGO_Soil`): each cell's dominant SSURGO map
  unit by area, its dominant (or a named) major component, and the
  mineral horizons in 0-30 cm depth-weighted into sand, clay, OM, bulk
  density, pH (1:1 water) and CEC (CEC7, else ECEC). Texture is
  **Cheshire** in every cell, by decision; pH and OM are each cell's own
  map unit's: Cheshire (pH 5.3, OM 2.7 %) on the 646 cells of the
  Cheshire-Holyoke units, Wilbraham (pH 5.5, OM 4.5 %) on the 290 cells
  of the wet Wilbraham-Menlo unit (`Cell_Review_SSURGO_and_NLCD.csv`).
  Cheshire is very deep and well drained (bedrock below 1.8 m); the
  column is, by decision, the base parameter file's default, 20 layers to
  5 m with free drainage. The trees root to a ZR95 of 1 m and the grass to
  0.25 m, the legacy table's values for their rows.
- **Soil chemistry**: SSURGO cannot supply the exchanger (Cheshire has
  no CEC7, and an ECEC of 1.6 cmolc kg-1 below its own sum of bases, 10)
  and SoilGrids has CEC at pH 7 but no base saturation, so the
  representative exchanger comes from the NRCS Kellogg laboratory data
  of the eight characterised Cheshire pedons
  (`KSSL_Cheshire_pedons_0-30cm.csv`, six within 13 km).
  `Representative_Soil_Chemistry` pools the 0-30 cm layers whose pH lies
  within 0.3 of the cells' area-mean pH (5.36) and that carry a KCl Al:
  three layers at pH 5.2-5.3 (CEC-weighted 5.22). The model takes their
  CEC at pH 7 (NH4OAc), 10.8 cmolc kg-1, with Ca 0.444, Mg 0.135, K 0.012
  and Na 0.006 of it (base saturation 0.598; each fraction the layers'
  summed cation over their summed CEC), KCl Al 0.032 and the rest H. The
  Appelo exchanger has a fixed capacity and carries the pH-dependent
  charge as exchangeable H, so its CEC is a buffered one, as the Energy
  Farm's is (13 cmolc kg-1 by summation of Mehlich-3 bases and
  Sikora-buffer H, Pihlap et al. 2025); the effective CEC (bases + KCl
  Al, 6.8), the model's value until 8 October 2026, leaves that charge
  out. SoilGrids' CEC at pH 7 at the parcel (12.9 cmolc kg-1, 0-30 cm
  median) is the cross-check. The
  pH window matters, because the exchanger follows pH: the one forested
  pedon (7.9 km away, pH 4.6) holds Al on 30-60 % of its CEC at pH 7
  (67-83 % of its effective CEC). SoilGrids'
  pH (5.0-5.3) agrees with SSURGO; its organic carbon does not (269 g/kg
  at 0-5 cm, a forest-floor value), so OM is SSURGO's. All of this is
  every cell's baseline (`Run_Hamden_Field_Inputs`); sub-plot data fill
  the `measured_*` columns and a warm start re-spins the field. Until
  they do, the base saturation is left to emerge from the spin-up (below):
  imposed, the representative set fails the ERW handoff check.
- **Vegetation** (`PFT_Hamden_Deciduous_Grass`, the row's `pft_override`):
  trees from the T&C site table's Harvard Forest row, grass from its
  Fendt C3 meadow row; a mature, dormant stand re-prescribed each
  spin-up year (only soil, its water and the plant reserves are carried,
  so the wood increment is the stand's base-cation sink); the grass cut
  to LAI 0.5 on 1 September with 90 % of the cut removed; and the forcing
  reference height 2 m above the 24 m canopy (`zatm` 26 m, where the
  base file's 8 m would put it below the zero-plane displacement).
- **Deposition**, event by event and year by year in place of the global
  maps for N, K, Ca and Mg (`Inputs/Deposition_Hamden.mat`,
  `Build_Site_Deposition`, `Deposition_Day`). Wet: NADP NTN CT15's
  precipitation-weighted mean concentration of each year's valid samples
  times each day's NLDAS precipitation, so it comes with the storms and
  wet years deposit more; a year short of NADP's completeness criteria
  (2015, 2019, 2020, and 2014 and 2026 at the record's ends) pools the
  three nearest complete years. Dry: TDep's (MMF 2025.01) dry deposition
  of the parcel's 4 km cell, measured air concentrations times modelled
  deposition velocities, NH3 and the gaseous N included, spread evenly
  over each day; 2025-2026 take the 2022-2024 mean. Over 2015-2025
  (kg ha-1 yr-1): N 6.6 (wet 2.6, dry 4.0), S 1.6 (1.0, 0.6), Ca 0.85,
  Mg 0.42, K 0.36; dry S halves over the decade (0.97 to 0.45) and dry N
  falls a fifth. The wet-only record had given N 2.7. P keeps the global
  map. The ERW chemistry takes its rain pH (5.20) from CT15's rain record
  (`Rain_Chem_Hamden.mat`, `rain_chem_file`) and its Na and strong-acid
  anions (SO4 + Cl) from the event record, wet and dry like the rest:
  34.4 umol m-2 d-1 of Na and 62.1 umolc m-2 d-1 of SO4 + Cl over
  2015-2025 (wet alone 23.7 and 44.1). Deposition is their only supply
  (`keyword_add` = 0, by decision), and the pore water's Na is seeded at
  its steady state on that supply. No fertiliser or manure, no
  fertiliser year (`FertN_years` empty) and no biological N fixation:
  atmospheric deposition is the only external nutrient supply, beside the
  primary-mineral supply. That keeps the borrowed case-5 uplift's Ca + Mg
  charge (1.15 keq ha-1 yr-1) split for Hamden: `B_IO.Tup_Ca` 0.0038876
  and `Tup_Mg` 0.0014617 g m-2 d-1 (14.2 and 5.3 kg ha-1 yr-1, against the
  uplift's 11.0 and 7.3), so that with the deposition Ca and Mg arrive at
  1.58 to 1, the dissolved ratio in which the measured exchanger (Ca/Mg
  3.28) is at equilibrium under SCEPTER's Ca/Mg selectivity. The
  biogeochemistry treats Ca, Mg and K alike by mass, so its spun-up
  solution settles at the molar ratios of their supply; the uplift's
  split had left it Mg-rich (Mg/Ca 1.22). A cold start's primary-mineral
  pools begin at the steady state of their supply instead of the base
  file's placeholder 502 g m-2 each.
- **CO2**: NOAA Mauna Loa daily means as given
  (`Ca_Data_NOAA_MLO_daily.mat`, `co2_file`).
- **Forcing**: NLDAS-2 hourly primary forcing for the 0.125 degree cell
  holding the parcel (41.4375 N, 72.8125 W), 2014-12-31 00:00 to
  2026-09-30 12:00 UTC, fetched as one point time series per variable
  through NASA Harmony's Giovanni time-series service
  (`Inputs/Hamden/Meteo/download_nldas2_point_series.py`; its hour stamps
  equal the native files' to rounding) and converted by
  `Build_TC_Forcing_From_NLDAS2` into `Data_Hamden.mat` (102,996 hours,
  no gaps; Lat 41.437, Lon -72.866, Zbas 97 m from the USGS 1 m 3DEP
  point query, DeltaGMT -5): 10.9-12.1 C, 764-1,442 mm/yr (mean 1,111).
  The converter reads the hourly NetCDF files (any window around the
  site) or one text series per variable (Giovanni CSV, Data Rods), joins a
  variable's files, and follows NLDAS-2's own time convention (README,
  Table 8): temperature, humidity, pressure, wind and shortwave are
  instantaneous at the hour, precipitation the total of the hour ending
  there, so the state variables are taken at the end of each hour as the
  ERA5 route takes ERA5's, precipitation as the hour's total, and
  shortwave as the mean of its two instantaneous ends. Vapour pressure
  comes from specific humidity and pressure, capped at saturation, and
  the dew point is that vapour pressure's own; the shortwave is split by
  `Automatic_Radiation_Partition_I`, as the ERA5 route does. A missing
  hour, a site outside the grid window and other units are refused.
  `test_build_tc_forcing_from_nldas2` (14 checks) covers both formats.
- **Site rows**: `Hamden_MOD_columns.csv` and `Hamden_ERW_columns.csv`
  in `Inputs/Field_Spinup/Hamden` put the site over the
  borrowed case 5: Location `Hamden` (which names `Data_Hamden.mat`),
  the vegetation file, the CO2 and deposition records, the windows (ring 2015-2018, the 2019 handoff
  extension, the state handed off at 2020-01-01 for the simulation and
  validation runs, 2020-01-01 to 2026-09-29, the last whole day of the
  forcing), the base-saturation seed, pH by a 1:1 extract
  (`pH_DilutionFactor` 1, as SSURGO measures it) and no rock.
- **DEM**: the 1 m 3DEP tile (CT statewide 2016 lidar, NAD83 UTM 18N)
  cut to the parcel and 200 m around it
  (`USGS_1M_18_x67y459_CT_Statewide_C16_parcel_window.tif`, 1.4 MB,
  not versioned), which gives the cells the tile's terrain exactly: 82.7
  to 101.5 m, slope to 0.49 (median 0.11). The margin only lets the edge
  cells' slope, curvature and wetness index see real terrain; the
  simulation is the whole field, with no buffer region.
- **Bouchet**: `v2/cluster/field_spinup_hamden.env.example` runs the
  flux-tracking spin-up on the 5 m majority grid (`TC_CELL_M`,
  `TC_MIN_COVERAGE`): one class seed (every cell shares texture and
  vegetation), 94 tracked and replayed cells, 842 from the forest; its
  header lists what must exist first and the warm start that follows the
  site observations.

`test_cell_land_cover_and_soil` (19 checks) covers the cover and soil
mapping on the parcel's grid; `test_hamden_field_configuration` (31) the
site as the cluster configures it, from the grid down to one cell's
assembled parameters (the weathering split and the exchange constants
included) and, with the DEM, the cluster's plan; `test_fit_exchange_affinity`
the calibration tool and the Al(OH)3 phase option.

The exchanger: until sub-plot data arrive, Hamden's ERW row lets the
base saturation emerge from the spin-up (`spinup_equilibrium`), and a
cell given a measured one (`measured_bs_*`) is seeded from it. The
exchange constants are Hamden's own, fitted the way the Energy Farm's
were by `v2/calibration/Fit_Exchange_Affinity.m`, which works for any
site. SCEPTER's `logKH` are shifted by z·δ at alpha 3.4, and the pinned
Al(OH)3 phase is fitted too, so that on a converged class seed's spun-up
solution the emerging exchanger holds the representative layers' base
saturation (0.598 of the CEC at pH 7) and KCl Al (0.032) at their pH
(5.22). The fit gives δ −1.042 and log \*Kso 9.20 (crystalline gibbsite
8.04, microcrystalline 9.66). The Energy Farm's own shift is −1.18. Its
amorphous phase (10.8) cannot hold Hamden's bases at any shift.

The fit's input is `Inputs/Hamden/Calibration/ERW_seed_inputs_classseed_01.mat`,
from a cold class seed of this code converged at loop 216; every class
seed now leaves such a file. The exchanger then emerges with:

- at the layers' pH: Ca 0.406, Mg 0.179 and K 0.013 (measured 0.444,
  0.135 and 0.012);
- base saturation 0.625 at bulk pH 5.3 (Cheshire cells, pore pH 4.92)
  and 0.678 at 5.5 (Wilbraham, pore pH 5.12);
- dissolved Al 9 and 3.5 µmol/L.

Both pass the handoff check. The base file's Energy Farm constants (alpha
1.3, amorphous Al(OH)3) had held only 0.15–0.28, and imposing the
representative set under them failed the check (dissolved Ca 28.7 mM
against 20).

CEC scenarios: until sub-plot data give the cells their own CEC, every
cell has the representative 10.8 cmolc/kg. The CEC contrast therefore
comes from scenarios: 6.4, 10.8 and 13.7 cmolc/kg, the range of the
representative layers' CEC at pH 7 (`CEC_min_cmolc_kg`,
`CEC_max_cmolc_kg`).

The spin-up does not depend on the CEC (the ERW chemistry is off in
every spin-up stage), so the scenarios share one spin-up and differ in
the exchanger's capacity alone. Each keeps the handoff's exchanger
fractions at its own CEC; the seed re-partitions them, within 1.4e-4 of
a fresh seed at that CEC. `v2/run/Run_Plot_CEC_Scenarios.m` runs one
plot-scale column per scenario from one spun-up state (about 17 h each
for 2020-2026). `Summarize_CEC_Scenarios.m` tabulates their daily bulk
and pore pH, base saturation and exchangeable Al and H. A transient for
a site without a pH record starts from the cell's own bulk pH.

Deposition alone supplies Na and SO4 + Cl (`keyword_add` = 0, by
decision). The seeded strong anions are the measured pH's charge
balance, organic anions the model does not carry included. They are
40-60 times what deposition sustains (39 ueq/L at the ring's water
loss), and leave with the water the ERW layer loses, which is little
before the leaves are out. In a 200-day plot-scale test from 1 January
2020, under the Energy Farm constants, they stayed near 2 meq/L until
June, then went with the transpiration stream, to 84 ueq/L by 18 July.
Over that time the pore pH rose from 4.91 to 5.13 and base saturation
from 0.145 to 0.160. The layer's leakage is T&C's flux out of the bottom
of the 5 m column (`Lk`), which is small through the winter.

### Energy Farm Bouchet R2025a runs

`RUNBOOK_2026-09-24.md` holds the Bouchet deployment and the order of the
campaigns. The twelve-run validation goes through
`runtime/Base_1D/submit_all_12.sh`, which needs the R2025a equivalence
baseline in place (`regen_equiv_1d.sbatch`). It submits the 1-D and 2-D
validations (6 CPUs each, 40 GB and 32 GB) and both production jobs held,
then releases them; production waits for both validations, after which
runs 10 to 12 and the 2-D array of runs 1 to 9 become eligible together.
Every job pins `cpumodel:9575f` and records its hardware, and production
needs a fresh R2025a validation stamp, stopping if source, inputs, tests
or the numerical platform change after validation. Calling a campaign's
batch file directly bypasses these checks.

From PowerShell at the repository root:

```powershell
matlab -batch "cd('v2/run'); Run_All_12_Parallel('Runs',10:12,'MaxConcurrent',3,'MaxRuntimeHours',168)"
matlab -batch "cd('v2/run'); Run_All_12_Parallel('Runs',1:9,'MaxConcurrent',2,'ParallelWorkers',6,'MaxRuntimeHours',168)"
matlab -batch "cd('v2/plotting'); Plot_1D_Then_All_12"
```

The plotting command first draws the vertical 1-D weathering and nitrogen
figure (`vertical_1d_weathering_nitrogen.png` and `.fig`). It then writes one
PNG and FIG pair per summary metric and one pair per time-series metric to
`v2/figures`. The weathering panels report solution Ca, Mg and Si in g
element m^-2, and N2 plus total N2O as cumulative gaseous N loss. The
reported denitrification N2O share removes the saved nitrification N2O flux
before comparing denitrification N2O with N2. Filenames begin with
`Subfield_12_summary_` or
`Subfield_12_time_dynamics_`. Each figure has one axis and contains control,
basalt and calcite for the 1-D, 2-D whole-field, 2-D cropped-buffer and 2-D
native-vegetation-buffer configurations. Treatment controls colour and
configuration controls line style. Cumulative nitrate leaching and the soil
nitrate pool remain separate metrics. Set `FigureDirectory` beforehand to
write the figures somewhere else. The all-12 stage is skipped with a warning
listing the missing runs until the 2-D results finish. NEE is respiration
minus crown-weighted NPP; positive values are carbon release to the
atmosphere. Carbon is normalised by the full 3600 m^2 domain area.
Daily rows are dated by the BIOGEO day they hold: row j runs from 01:00 of
Date(1)+(j-3) to 00:00 of the next day, and row 2 is the lone first hour.
Daily WFPS averages the same hours, and the initial retained state has no
preceding interval and is plotted as `NaN`.

The completed 1-D runs can be plotted before the 2-D runs finish, from
`v2/plotting`:

```matlab
Runs = 10:12; Plot_1D_Then_All_12
```

This writes separate PNG and FIG files for each metric. With only runs 10 to
12 available, each comparison figure contains 1-D control, basalt and calcite.
When runs 1 to 9 finish, the three 2-D groups can be plotted without
reloading the 1-D results:

```matlab
Runs = 1:9; Plot_1D_Then_All_12
```

This writes separate `Subfield_2D_summary_*` and
`Subfield_2D_time_dynamics_*` PNG and FIG files. Each contains control,
basalt and calcite for the three available 2-D configurations. Running the
script without setting `Runs` plots the 1-D group first and then compares
all twelve treatment and configuration combinations.

From an open MATLAB session at the repository root:

```matlab
cd v2/run
Run_All_12_Parallel('DryRun',true)
Run_All_12_Parallel('Runs',10:12,'MaxConcurrent',3,'MaxRuntimeHours',168)
Run_All_12_Parallel('Runs',1:9,'MaxConcurrent',2,'ParallelWorkers',6,'MaxRuntimeHours',168)
cd ../plotting
Plot_1D_Then_All_12
```

The first command completes the three vertical 1-D comparisons before the
2-D runs begin. Stop between the two commands if the second stage needs to
wait. `Ctrl+C` stops active child runs; completed result files are retained
and skipped when the same command is run again. On Windows the launcher
stops the complete worker process tree, so a MATLAB worker is not left behind.
The chemistry reruns started with `run_phcharge(10)`, `run_phcharge(11)` or
`run_phcharge(12)` write a rolling checkpoint every 30 minutes to
`results/checkpoints`. A matching checkpoint resumes on the next launch. The
checkpoint is removed after the final result is saved.

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

`ParallelWorkers` only applies to the independent daily ERW chemistry in
each 2-D run. Keep `MaxConcurrent * ParallelWorkers` within the available
CPU allocation. Daily biology, hourly physics and lateral routing remain
serial. Full 12-column runs also hold several GB of state.

### Current campaigns (2026-09-25)

The calibrated set and every campaign that reads it live under `runtime/`;
`RUNBOOK_2026-09-24.md` at the checkout root orders them for Bouchet:

* `Base_1D` and `Base_2D`: the 5-year validation, all twelve runs,
  through `runtime/Base_1D/submit_all_12.sh` (optional
  `TC_CALCITE_CASE=8` adds the 5 t/ha calcite runs as a campaign of their
  own). `Campaign_Parameters.m` is the only file that spells the
  calibrated numbers, and `phobs.m` scores August pH against both
  observed depth intervals and `budget.m` prints the alkalinity budget
  and exchanger state of the control and basalt runs.
* `SSP_1D` and `SSP_2D`: the ssp126/245/585 projections to 2050 through
  `runtime/SSP_1D/submit_ssp.sh` (`TC_SSP`, `TC_HORIZON`, optional
  `TC_CALCITE_CASE=8`).
* `Calibration_Sweep_1D`: the joint sweep that calibrates alpha and
  KNmaxd, with the Kbot alternatives beside them.
* `Uncertainty_1D` and `Uncertainty_2D`: the parameter-uncertainty
  ensemble, the T&C-SMEW field paper's sensitivity analysis in the coupled
  model. Beside the original set, reused from the Base and SSP results,
  all sixteen parameters run at the lower and at the upper end of their
  ranges, in 1-D and 2-D, with the 5 t/ha arm, over both horizons
  (`runtime/Uncertainty_1D/submit_uncertainty.sh`). `Run_One_Subfield`
  carries its three knobs (`BiogeoOverrides`, `DissolutionFactor`,
  `SaturationExponent`), bit-identical when unset.
* `Worker_Scaling_Test`: sizes the 35-year 2-D jobs.

Figures: `v2/plotting/Plot_1D_Then_All_12.m` (`Scenario` for projection
results). Its `UncertaintyFile` bands and `v2/plotting/Plot_Uncertainty.m`
read the retired 36-member summary; the three-set ensemble's summary and
figures come with its post-processing package.

## Model lineage from v0.9

These are the v0.9 package notes, newest first, as shipped with package w.
v1.0 contains every change they describe; the files they cite are at the
repository root.

Package w is package v with one guard relaxed (`ROUNDOFF_NEGATIVE_2026-10-05.md`):
`ERW_Ionic_Strength` counts a concentration within 1e-10 umol-conv/L below
zero as zero instead of stopping the run. A 35-year 2-D uncertainty run
stopped on -4.5e-16, a cancellation residue. Every run that got past the
check gives the same numbers in v and w. One audit test file differs from
v's (`test_erw_activity_coefficient`), and the comparison utility at the
package root (`compare_u_v.m`) now also compares two results that both
keep daily particle rows, as v and w do, and compares the fsolve options
by their settings. `RUNBOOK_2026-09-24.md`, section 0w, is the procedure.
The paragraph below describes v, which w contains.

Package v is package u with two changes to how a run holds and saves its
state, and none to the model (`PARTICLE_RECORDS_DAILY_2026-10-04.md`).
The particle-class histories of an amended column (`psd`, `d`, `lamb`,
`SSA`, `delta_d`) keep their last three days whole and, of every earlier
day, its last row: the model reads only the previous row and each
application's first row, so every computed value is u's, and the results
hold those histories at the initial rows and the last row of every ERW
day (`weath.particle_rows`) instead of every ten-minute row. Held whole
they took 16.8 GiB per amended column over 35 years, about half of a
whole-field basalt run's memory. The SSP 2-D campaign checkpoints every
two hours instead of every thirty minutes, and the 35-year 2-D runs ask
256 G instead of 384 G, so six run together on the week partition. One
audit test is rewritten (`test_erw_particle_history`) and one pinned line
changes (`test_uncertainty_campaign`). `RUNBOOK_2026-09-24.md`, section
0v, is the procedure for this root. The paragraph below describes u,
which v contains.

Package u is package t with three changes. The soil-water paste of the
bulk-pH operator stands open to the laboratory air
(`PASTE_OPEN_TO_AIR_2026-10-02.md`): its dissolved CO2 is the air's at
25 C instead of the CO2 the field pore water carried, and the plot's solid
carbonates still dissolve in it, so a soil holding free calcite reads the
calcite-air equilibrium of its paste, about 8.1, and a soil without solid
carbonate reads within about 0.01 of the closed paste in its August mean
(up to 0.1 apart on single days). The
campaign's Appelo alpha is 3.4, the Appelo (1994) slope SCEPTER carries as
its default, in place of the 1.3 fitted on 26 September, and the lower and
upper sets of the uncertainty ensemble carry alpha 1.3 and 6 around it
(`runtime/Base_1D/Campaign_Parameters.m`,
`runtime/Uncertainty_1D/README.md`). Both campaign controllers can submit
their 1-D runs first (`TC_ALL12_DIMS=1d`, `TC_UC_DIMS=1d`), and
`readout_1d.m` at the package root prints the August read-out of those
runs. Every trajectory differs from t's: the seed is inverted through the
open paste, and the soil biology reads the open reading. Two audit test
files differ from t's (`test_pH_conversion`, which pins the open paste, and
`test_uncertainty_campaign`), with the all-12 configuration test.
`RUNBOOK_2026-09-24.md`, section 0u, is the procedure for this root. The
paragraph below describes t, which u contains.

Package t is package r with one change to the model and the 1-D sweep
returned to dissolution factor 1 (`AMPHIBOLE_KINETICS_2026-10-02.md`). The
basalt's amphibole (ferroactinolite, 11.6 % of the rock) ran on the
tremolite row of Palandri and Kharaka (2004), which is taken from Sverdrup
(1990) rather than fitted to tremolite experiments and is about 60 times
the measured rate at pH 6. That one mineral supplied 11.5 of the 12.5 eq
m-2 of silicate alkalinity by August 2020 and used up the exchanger's
10.3 eq m-2 of acidity in 2019. It now carries the amphibole rate law of
Hermanska et al. (2022). Re-evaluated offline from saved runs, the rock
then releases about 0.7, 1.5, 2.5 and 3.5 eq m-2 by August 2017 to 2020 at
dissolution factor 1, against a field exchanger gain of 1.2, 2.3, 3.0 and
3.0. Control and calcite runs hold no amphibole and are r's bit for bit.
The sweep is package o's again: seven alphas, four KNmaxd, four bottom
boundaries, 112 combinations, 336 trajectories, with `DissolutionFactor` 1
in every case. Two audit test files differ from r's
(`test_erw_mineral_kinetics`, which pins the new constants, and
`test_bouchet_1d_joint_sweep_configuration`). Packages p, q, r and s are
superseded for anything that holds basalt. The paragraph below describes
r, which t contains.

Package r is package q with the plot's solid Mg carbonate (the nesquehonite
pool) dissolving in the soil-water paste beside its calcite
(`PASTE_MG_CARBONATE_2026-10-01.md`). The pool is zero in every five-year
run, so the sweep, the controls, the calcite runs and the five-year basalt
runs are q's bit for bit; the 35-year basalt runs read a higher bulk pH
once the pool has formed. Two audit test files differ from q's:
`test_pH_conversion` (new checks) and
`test_bouchet_1d_joint_sweep_configuration` (the wording of one message).
The notes now say "dissolution factor" where p and q said "release
factor"; the option, the workbook column and every number are unchanged.
The paragraph below describes q, which r contains.

Package q is package p with one change to the model: the bulk-pH operator
lets the plot's solid calcite dissolve in the soil-water paste, as it does
in a real pH measurement, instead of freezing it with the silicates
(`PASTE_CALCITE_2026-10-01.md`). The pore-water chemistry is untouched; a
plot holding free calcite now reads 7.8 to 8.3 where it read 7.4 to 7.8,
and plots without solid calcite are unchanged bit for bit. Two audit tests
differ from o's (`test_bouchet_1d_joint_sweep_configuration`,
`test_pH_conversion`). The paragraph below describes p, which q contains.

Package p is package o with the dissolution factor (Run_One_Subfield's
`DissolutionFactor`, the workbook `diss_f`) added as the fourth dimension of
the 1-D calibration sweep and the other three grids trimmed to keep the
sweep's cost: 150 combinations, 450 trajectories
(`DISSOLUTION_FACTOR_SWEEP_2026-10-01.md`). The model, every campaign wrapper
and 93 of the 94 audit tests are o's unchanged; the inherited export
operator is kept as inherited by decision. Case 1 to 3 of the sweep remain
the reference tuple, so the one-year gate and the all-12 campaign label
(`a1p3_kn0p0048_kbfree_fd1`) carry the same physics as o. The paragraphs
below describe o, which p contains.

Package o is package n with its Bouchet validation repaired and made
diagnosable (`VALIDATION_FIXES_2026-10-01.md`): the equivalence baseline
now fingerprints only the four input files it reads, so the scenario
forcings copied into a root cannot stale it; a refused baseline names the
files that moved; both 1-D validators print a timestamped banner at
every stage; and the one audit test n's salt closure broke
(`test_erw_lateral_cation_daily_delta`, whose last-day expectation
predated a boundary Ca and Mg supply) now accounts for it. The model is
n's, unchanged, and the audit suite and smoke runs were executed on it
before packaging. Everything below describes
that model.

Package n lets the 0-30 cm box export what its drainage water carries.
Dissolved Ca and Mg now leach and are taken up at their concentration
(`aCa = aMg = 1`, as NO3 and K already were; the legacy 0.4 belonged to a
lumped pool), and the leached background anion returns as the salt it left
as instead of a bare anion. In package m the leachate carried 2.5 times
less Ca+Mg charge than anion charge, 84 % of the alkalinity the rock
released stayed in the box, the exchanger saturated by 2019 and the basalt
pH ran away; the field reached pH 6.9 at 83 to 90 % base saturation. The
sweep now scores that base saturation beside pH, over August, when the
soil was sampled, with the alpha grid extended to 6; the calibrated alpha
is to be re-chosen from it before production. Evidence, equations and
checks are in `LEACHING_CLOSURE_REVIEW_2026-10-01.md`. The Mg carbonate
sink of m (`MG_CARBONATE_REVIEW_2026-09-30.md`) and the particle-history
storage of l (`PARTICLE_HISTORY_REVIEW_2026-09-28.md`) are kept. Every run
changes, the control included, so this build needs its own Bouchet R2025a
equivalence baseline, validation and sweep before production; use the
runbook here.
