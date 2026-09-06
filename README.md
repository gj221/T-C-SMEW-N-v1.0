# T&C-SMEW-N v1.0

## Repository layout

```
v2/run/          launchers: ERW_run, Run_Matrix_12, Run_One_Subfield,
                 Run_All_12_Parallel, Run_Subfield_Worker,
                 Terminate_Process_Tree
v2/plotting/     figure scripts: Plot_1D_Then_All_12, Plot_Grid_Setup,
                 Plot_Biogeo_Pools, Plot_Nflux_Diagnostics,
                 Plot_Kantola_Benchmark, and the Cflux_post_processing helper
v2/interpolation/  soil-sample kriging and the grid figure it feeds:
                 Run_Notebook_Interpolation, Run_AW_Bruce_Interpolation,
                 Run_Bruce_Interpolation, Run_Yield_Field_Interpolation,
                 Run_Cell_Size_Analysis, Plot_Field_Grid_Interpolation
v2/figures/      every PNG and FIG the plotting scripts write
v2/Spatial inputs/   field polygons used by full-grid simulations
v2/T&C_Code/     model source, Inputs, Parameters and Tests
v2/spinup/       the automatic spin-up and the field (per-cell / RF) launchers
v2/cluster/      SLURM scripts and entry point for the field spin-up on Bouchet
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
[state_1,config_1] = Run_Field_Cell_Spinup(shape_file,files{1}, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',1);
[state_2,config_2] = Run_Field_Cell_Spinup(shape_file,files{2}, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',2);
```

The site workbook row must provide `meteo_file`, `meteo_start` and `meteo_end`
before spin-up can run (an hourly forcing `.mat` is accepted as
`meteo_file`; a row that lacks them can be completed per cell with the
launcher's `MODParamColumns` option, see *Spin-up forcing* below). This
is a hard preflight check. The imported machinery
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

### Field topography

Following Lian et al. (2025, *WRR* 10.1029/2025WR040260, T&C-BG-2D), the
2-D field driver can resolve elevation-driven lateral and vertical flow.
Attach a DEM to any multi-cell topology and the lateral solver switches
from matric-head-only to **total hydraulic head** `H = psi + z`
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
width and adds lateral water moved per m2 per day and the spatial
coefficient of variation of end-of-run root-zone water.

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
promoted only when all three pass the gate. They are measured-only
attributes: there is no pedotransfer baseline for the exchanger, so the
template carries no `baseline_bs_*` column, `Apply_Kriged_Cell_Values`
needs only a `measured_` column and a toggle to recognise an attribute,
and a cell without a measurement keeps the existing behaviour, in which
the handoff derives the exchanger composition from the pore water and
the bulk pH alone. Where the measurement is used,
`Resolve_Cell_Agronomic_Inputs` checks it (each fraction in (0, 1), the
three summing to less than one; the remainder is the exchangeable
acidity), reports the exchangeable Ca, Mg and K it implies in
cmolc kg-1 against the cell's resolved CEC, and
`Configure_Cell_Agronomic_Inputs` writes the cell's ERW exchanger seed
`f_CEC_in = [Ca Mg K Na Al H]` as `[bs_Ca bs_Mg bs_K 0 0 1-sum]`: the
acidity share starts as exchangeable H+, and Na and Al exchange start
at zero and follow the handoff equilibrium. That seed is the right
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

### Charge-balanced random-forest spin-up

`Run_Field_Cell_Spinup(..., 'Method','rf')` reconstructs the distributed
biogeochemistry pools from a topographically-stratified subsample of full
spin-ups plus a random forest on {elevation, slope, aspect, curvature,
TWI, soil class, vegetation} -> conserved pool totals (Lian et al.
Fig. 3e; needs the Statistics and Machine Learning Toolbox, falls back to
`'per_cell'` otherwise). Every reconstructed cell is seeded from its
nearest subsample neighbour and finished with a short
`Run_Automatic_Site_Spinup` that starts near equilibrium and runs the
standard ERW handoff gate (`ERW_Handoff_Readiness`), so electroneutrality
and pore-water pH are closed by the same code the transient run trusts;
cells that fail the gate fall back to a cold-start spin-up. The
charge-closure primitive `ERW_Initial_H_From_Anion_Budget` is ported from
the `Model/TC-SMEW-N` branch `lian-spinup-convergence`.

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
the workbook row it replaces. `Tests/test_spinup_meteo_loader.m` (10
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
temperature, saturation, matric potential, water sums and WFPS of the
active layer) that Phase 1 averages into the mean forcing of the
stand-alone Phase 2, and no longer clears the field-capacity constant
Phase 2 reads. The plant nutrient reserve is a physical inventory that
`Plant_Reserve_Growth_Limit` never overdraws, so `VEG_DYN_RES` clears
the integration round-off of a balance that lands on zero within the
same 1e-6 g m-2 resolved-mass tolerance the spin-up applies to the soil
pools (-4e-11 g N m-2 at a cold-start emergence had tripped the
non-negativity guards, which now report the offending value;
`Tests/test_plant_reserve_fraction.m`). And the spin-up and the field
launcher read the parameter workbook through one function, so a per-cell
override row is typed like the row it replaces.
Two more were interface drift between the two code lineages: the
stand-alone Phase 2 called `BIOGEO_UNIT` with the Model-repo signature
(it now passes v1.0's tail arguments and reads a new `P_inputs` output,
fertiliser, deposition and uplift phosphorus), and the Phase-3 nutrient
budgets read per-crown ledgers the Model-repo driver keeps, so the
10-minute driver now records the per-crown harvest reserve export and
planting input (`AgrHarNut_H/L`, `PlantingNut_H/L`, `PlantingP_H/L`)
and the daily phosphorus inputs, while the nitrogen gas loss is summed
from v1.0's own partition (volatilisation, N2 from denitrification and
total N2O; there is no NO pathway) and the fertiliser term is what
`BIOGEO_UNIT` added that day (`FertN_d`).

**Cost, measured on one Control 2 cell with the Kantola case (seven
cycle years, 2017 handoff), single column, this machine:** Phase 1
(T&C only) runs at 48 s per forcing year, 5.6 min for the ring; Phase 2
(stand-alone biogeochemistry under the mean forcing) converged after
641 mean-forcing years in about 3 min; Phase 3, the fully coupled loop
with the ERW chemistry solved daily, costs about 70 min per forcing
year, and its convergence test is taken over complete rings of seven
years. On the test cell the ring drift of the slow fixed-cation pools
(B(54), B(62), B(70)) was still 0.24 per ring against the 0.01 limit
after 27 coupled years, falling by about 4 % per year, so a cell is a
job of days, not hours.
That is the cost the random-forest path exists to avoid: the production
2.33 m Control 2 grid has 1461 cells, and even the ~10-15 % subsample
the RF trains on is weeks of single-threaded compute, so the field
spin-up has to be run cell-parallel (the `Run_All_12_Parallel` pattern)
and is not something a session completes end to end. The pipeline up to
the spin-up call (templates, kriging, gated merge, resolved per-cell
inputs, topographic subsample, spin-up configuration) was exercised on
the real Control 2 grid (54 cells at 15 m under a 50 % coverage rule;
14 in-hull cells took every kriged attribute including base
saturation) and is covered by the tests; the single-cell spin-up above
is the integration check of the spin-up itself: it completed Phases 1
and 2 and 27 coupled Phase-3 years with all six nutrient budgets closing
every year, and was stopped unconverged (above). The field spin-up is
therefore a cluster job: `v2/cluster/README.md` describes the SLURM
array on Bouchet (one cell per task, finished cells reused, plan /
assemble / validate steps) and the validation experiment that scores
the random-forest path against full spin-ups of every cell.

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
