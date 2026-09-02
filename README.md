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

### Field topography

Following Lian et al. (2025, *WRR* 10.1029/2025WR040260, T&C-BG-2D), the
2-D field driver can resolve elevation-driven lateral and vertical flow.
Attach a DEM to any multi-cell topology and the lateral solver switches
from matric-head-only to **total hydraulic head** `H = psi + z`
(`Darcy_Lateral_Topology`), overland runoff is routed down the DEM by
Quinn (1991) multiple-flow-direction weights and leaves only at the
topographic outlet (`Surface_Runoff_Routing`), and the direct solar beam
is projected onto each cell's slope and aspect with a sky-view diffuse
correction. Without a DEM every run is flat and bit-identical to before
(`test_2d_driver_equivalence` still passes at tolerance 0).

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
    'FeatureIndices',[1 2]);

%% compute the process-grid cell size from the site's own length scales
Run_Field_Length_Scales( ...
    fullfile('v2','Spatial inputs','Bruce Farm Virginia','Whole farm', ...
    'AdamWilbourne_Control_Fields.shp'), ...
    fullfile('v2','Spatial inputs','Bruce Farm Virginia','DEM'), ...
    'FeatureIndices',[1 2]);
```

`Run_Field_Length_Scales` fits an empirical variogram to each measured
Bruce Farm soil attribute (shortest range = `L_soil`) and to the
detrended-elevation semivariogram of the DEM (`L_terr`), and recommends
`h = max(DEM spacing, min(L_soil, L_terr) / 3)` capped by the 95 %
per-cell coverage rule. That fixed 80 m "notebook recipe" variogram range
in `Fit_Ordinary_Kriging` is a default, not a Bruce Farm measurement, and
is not used for grid sizing.

A topographic run passes `FieldDEM` (and optionally `DEMResolutionM`) to
the launcher; the DEM is fetched with `Load_Field_DEM` (mosaic + reproject
+ resample) and attached with `Attach_Topography_To_Topology`:

```matlab
Run_One_Subfield(1,'FieldShapefile',shape_file, ...
    'FieldFeatureIndices',[1 2],'FieldPolygonIndex',1, ...
    'FieldDEM',fullfile('v2','Spatial inputs','Bruce Farm Virginia','DEM'), ...
    'DEMResolutionM',5)
```

The `v2/Spatial inputs/**/DEM/` GeoTIFFs are fetched, not versioned
(gitignored). For the Adam Wilbourne control fields, the USGS 3DEP 1 m
tiles `USGS_1M_17_x71y405..x72y407` from the
`VA_FEMA-NRCS_SouthCentral_2017_D17` project cover both control polygons.

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
