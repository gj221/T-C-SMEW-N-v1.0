# Hamden (145 Blue Trail) — T&C-SMEW-N v1.0 model inputs

Site: Hamden, Connecticut, parcel 17852 (145 Blue Trail), 2.34 ha
(23,399 m², 666 m perimeter). Centre 41.437 N, 72.866 W; extent
41.436–41.438 N, 72.867–72.865 W (to 3 d.p.).

Study: no rock is added. The runs track how differences in CEC change bulk
soil pH under atmospheric deposition alone. Deciduous temperate broadleaf
trees over grass that is mown once a year (1 September); no liming,
fertiliser or other management. Field-scale run with the DEM, on Bouchet.

The same files are kept in two places:

| Where | Purpose |
|---|---|
| `T&C-SMEW-N v1.0/v2/T&C_Code/Inputs/Hamden/` | source data the model setup reads (deposition, CO₂, soils, vegetation; NLDAS-2 meteo once downloaded) |
| `T&C-SMEW-N v1.0/v2/T&C_Code/Inputs/` | model-ready records the loaders look for: `Data_Hamden.mat`, `Deposition_Hamden.mat`, `Ca_Data_NOAA_MLO_daily.mat`, and `Rain_Chem_Hamden.mat` (kept, not read; see below) |
| `T&C-SMEW-N v1.0/v2/Spatial inputs/Hamden Blue Trail/` | site polygon, land-cover clip, DEM (when downloaded) |
| `Documents/Model data/Hamden data/` | documentation copy of everything, by data type |

## Site boundary

- `Hamden_Blue_Trail_Parcel.shp/.shx/.dbf/.prj` — the parcel polygon from the
  Connecticut 2025 statewide parcel layer (`2025_Parcel_Layer.gdb.zip`,
  feature Location = "145 BLUE TRAIL", Town = Hamden, Parcel_ID 17852,
  Unit_Type "Single Fam M01"). Kept in the layer's own CRS, NAD83(2011) /
  Connecticut (ftUS). Attributes: Field, Parcel_ID, CAMA_Link, Town,
  Unit_Type, Area_m2, Perim_m. The layer carries no owner fields.
- `Hamden_Blue_Trail_Parcel_WGS84.zip` — the same polygon in geographic
  WGS84 for upload forms (USGS 3DEP); coordinates are NAD83(2011) latitude
  and longitude labelled WGS84 (≈1 m, irrelevant for selecting data).
- `Hamden_Blue_Trail_Parcel_StatePlane.zip` — the State Plane original, zipped.
- `Site_Boundary_Hamden_Blue_Trail.png` — the parcel among its Blue Trail
  neighbours.

## Meteorology

- **Source: NLDAS-2 hourly forcing** (`NLDAS_FORA0125_H` v2.0, 0.125°,
  about 11 × 14 km; the whole parcel lies in one grid cell), for the cell
  containing 41.437 N, 72.866 W. precip.ai is not used.
- **Downloaded** (6 October 2026) as one hourly point time series per
  variable -- Tair, Qair, PSurf, Wind_E, Wind_N, SWdown, Rainf -- for the
  0.125 degree cell holding the parcel (centre 41.4375 N, 72.8125 W),
  2014-12-31 00:00 to 2026-09-30 12:00 UTC (the latest available), through
  NASA Harmony's Giovanni time-series service:
  `Meteo/download_nldas2_point_series.py` (Earthdata Login from
  `~/.netrc`; the Giovanni application must be authorised once in the
  Earthdata profile). Fourteen CSVs (the record to 2026-01-01 and the 2026
  extension) are kept in the documentation copy (`Meteo/NLDAS-2/point_series`).
  Their hour stamps equal the native NLDAS files' to rounding, checked
  against five native `..._subsetted.nc4` hours (a one-hour shift would
  differ by 0.15-0.4 m/s in wind and 0.2 K in air temperature). The bulk
  alternatives -- "Get Original Files" (whole-continent granules, ~156 GB
  for these dates) and the per-hour subsetter (Transform Data, ~96,000
  one-cell files) -- are read by the converter too.
- **Converted**, in MATLAB from `v2/T&C_Code`:
  `Build_TC_Forcing_From_NLDAS2('<folder of the CSVs>','Inputs/Data_Hamden.mat',41.437,-72.866,97,-5)`
  -> `Inputs/Data_Hamden.mat`, 102,996 hours, 2014-12-30 20:00 to
  2026-09-30 07:00 local standard time, no gaps; 2015-2025 means
  10.9-12.1 C, 764-1,442 mm/yr precipitation (mean 1,111; 2015-16 the
  drought years), ~180 W m-2 shortwave, 4.1-4.7 m/s wind.
  It follows the repository's ERA5 route
  (`Build_TC_Forcing_From_ERA5_Land_Timeseries`) and NLDAS-2's own time
  convention (NLDAS-2 README, Table 8): temperature, humidity, pressure,
  wind and shortwave are instantaneous at the hour, precipitation the
  total of the hour ending there. So the state variables are taken at the
  end of each hour, precipitation as the hour's total and shortwave as the
  mean of its instantaneous values at the hour's two ends; vapour pressure
  comes from specific humidity and pressure (capped at saturation) and the
  dew point from that vapour pressure; T&C's
  `Automatic_Radiation_Partition_I` splits the shortwave into direct and
  diffuse, visible and near-infrared, with the cloudiness that sets
  longwave. Hours are converted to local standard time (UTC-5). Zbas 97 m
  is the USGS 3DEP 1 m elevation at the parcel centre. Each grid cell's
  own slope and aspect from the DEM then project the radiation onto the
  terrain. The spin-up reads 2015-01-01 to 2019-12-31 and the simulation
  and validation runs 2020-01-01 to 2026-09-29 (the last whole local
  day), both checked complete with the spin-up's own loader.
- **Against precip.ai** (PRISM 800 m daily, 2020-01-01 to 2026-09-30;
  not used as forcing). Compared over the same 24 hours, PRISM's day ends
  at 12 UTC on its date: daily r = 0.96 on that day, 0.34 on calendar days.
  - Rain falls on the same days in both, but NLDAS-2 carries 14 % less of
    it: lower in 69 of 81 months, with fewer heavy days (74 against 104
    of 25 mm or more).
  - Mean temperature agrees (+0.2 C). NLDAS-2's daily range is 3.2 C
    narrower and its dew point 2.5 C higher, so its summer maximum vapour
    pressure deficit is about 40 % lower.
  - Its shortwave is 15 % higher.
- **For the Bouchet run: CW3E**, the 1-km, 1-hourly forcing downscaled
  from NLDAS-2/HRRR with Stage IV and PRISM (Pan 2025,
  doi:10.5281/zenodo.14714512).
  - It is fetched for the parcel's 1 km cell from HydroFrame's
    `hf_hydrodata` service, 2015-01-01 to 2026-09-30, by
    `Inputs/Preprocess/Fetch_CW3E_Point.py`. That needs the account
    holder's HydroFrame API PIN, registered on the machine.
  - The fetch writes the converter's point-series format. The hour labels
    and units are checked against NLDAS-2 before `Data_Hamden.mat` and the
    deposition record are rebuilt from it.
  - Local tests keep NLDAS-2. NLDAS-3 (1 km) is set aside while it is in
    beta.

## Atmospheric deposition

- `Deposition/NTN-ct15-w-s-mg.csv` — NADP National Trends Network site CT15
  (Abington, CT), weekly wet-only samples, concentrations in mg/L
  (−9 = missing), 2014-12-30 to 2026-05-19, 594 weeks.
- `Deposition/CT15_wet_deposition_summary.csv` and `Rain_Chem_Hamden.mat` —
  derived with `Inputs/Preprocess/Extract_NTN_Rain_Chemistry.m`: the
  precipitation-weighted mean concentration of the 424 valid samples
  (validity code `w`, all ions and pH present), NADP's convention, times
  the site's own precipitation over the same window (30 December 2014 to
  19 May 2026): 1,098 mm/yr from the NLDAS-2 forcing, `Data_Hamden.mat`,
  where CT15's gauge, in north-eastern Connecticut, caught 1,415 mm/yr
  (82 % of it sampled validly). The concentrations are regional and the
  amount local, as in NADP's own deposition maps, and the deposition then
  rides on the rain the model receives. Rain pH 5.20. Wet deposition
  (kg ha⁻¹ yr⁻¹): N 2.68 (NH₄-N + NO₃-N), S 0.99, Ca 0.57, Mg 0.28,
  K 0.26, Na 1.99 (at the gauge's precipitation: N 3.45, S 1.28, Ca 0.74,
  Mg 0.36, K 0.33, Na 2.57). Strong-acid anions (SO₄ + Cl) 14.7 µeq/L =
  44.1 µmolc m⁻² d⁻¹. To rebuild, from `v2/T&C_Code`:

  ```matlab
  ntn = 'Inputs/Hamden/Deposition/NTN-ct15-w-s-mg.csv';
  w = Extract_NTN_Rain_Chemistry(ntn,'').window;
  F = load('Inputs/Data_Hamden.mat','Date','Pr');
  in = F.Date > w(1) & F.Date <= w(2);
  Extract_NTN_Rain_Chemistry(ntn,'Inputs/Rain_Chem_Hamden.mat',[],[], ...
      sum(F.Pr(in))/(w(2)-w(1))*365.25);
  ```
- `Deposition_Hamden.mat` (in `Inputs`) — what the biogeochemistry uses,
  event by event and year by year (`Inputs/Preprocess/Build_Site_Deposition.m`;
  the model applies it day by day with `Deposition_Day`):
  - **wet** — CT15's precipitation-weighted mean concentration of each
    calendar year's valid samples (`Deposition/CT15_wet_concentration_by_year.csv`)
    times each day's NLDAS precipitation, rain or snow, so it comes with
    the storms and wet years deposit more. A year counts when its samples
    span at least 90 % of it and its valid samples caught at least 75 % of
    its precipitation (NADP's completeness criteria); 2015 (74.5 %), 2019
    (69 %) and 2020 (50 %) fall short, and they and the record's ends,
    2014 and 2026, take the pooled concentration of the three nearest
    complete years;
  - **dry** — the NADP Total Deposition (TDep, MMF version 2025.01) dry
    deposition of the 4 km cell holding the parcel (NAD83 / Conus Albers,
    row 224, column 1066), each year spread evenly over its days
    (`Deposition/TDep_dry_deposition_parcel_cell_by_year.csv`): measured
    CASTNET air concentrations times CMAQ deposition velocities, N
    including HNO₃, particulate NO₃ and NH₄, NH₃ and the modelled gases
    (NO₂, PAN, N₂O₅, HONO, organic nitrate), S as SO₂ + SO₄. The grids,
    `{n,s,ca,mg,k,na,cl}_dw-2015.zip` to `-2024.zip` (70 zipped GeoTIFFs,
    95 MB, from gaftp.epa.gov/castnet/tdep/CURRENT_grids, downloaded
    6 October 2026), are in the documentation copy
    (`Atmospheric deposition/TDep dry deposition`); 2025 and 2026 take the
    2022-2024 mean.

  The year-by-year totals at the site's precipitation are in
  `Deposition/Hamden_deposition_by_year.csv`. Over 2015-2025
  (kg ha⁻¹ yr⁻¹):

  | | N | S | Ca | Mg | K | Na | Cl |
  |---|---|---|---|---|---|---|---|
  | wet | 2.56 | 0.99 | 0.57 | 0.30 | 0.27 | 2.14 | 3.81 |
  | dry | 4.03 | 0.61 | 0.28 | 0.13 | 0.09 | 0.74 | 0.70 |
  | total | 6.59 | 1.59 | 0.85 | 0.42 | 0.36 | 2.89 | 4.51 |

  Dry N falls from 4.6 (2015) to 3.7 (2024) and dry S halves (0.97 to
  0.45); wet N follows the year's rain, 1.6 in 2015 (764 mm) to 3.4 in
  2021. To rebuild, from `v2/T&C_Code`, with the TDep zips in `tdep`:

  ```matlab
  Build_Site_Deposition('Inputs/Hamden/Deposition/NTN-ct15-w-s-mg.csv', tdep, ...
      'Inputs/Data_Hamden.mat', 41.437, -72.866, 2014:2026, 'Inputs/Deposition_Hamden.mat');
  ```
- In the model each element enters once: N, K, Ca and Mg through the
  biogeochemistry's deposition terms (replacing the built-in global maps,
  which give 13.4 kg N ha⁻¹ yr⁻¹ at Hamden, wet plus dry, older and
  coarse); P is measured by neither network and stays on the global map;
  the ERW chemistry takes no rain record for now (the site row leaves
  `rain_chem_file` empty): the 2-D transient never reads one, so with it
  the spin-up's handoff check had seeded a rain Na the transient does not
  use. CT15's Na and strong anions (`Rain_Chem_Hamden.mat`, above) wait on
  the `keyword_add` decision; nitrate is not counted among the strong
  anions because it is already N deposition.

## Atmospheric CO₂

- `CO2/co2_daily_mlo_1974-2026.csv` — NOAA GML Mauna Loa daily means,
  1974-05-19 to 2026-10-05 (source file `co2_daily_mlo (1974-Oct 2026).csv`).
- `Ca_Data_NOAA_MLO_daily.mat` — the same daily values as given (Date_CO2 at
  00:00 of each day, Ca in ppm), in the format of v0.9's `Ca_Data.mat`;
  the model's loader (`Load_CO2_Record`) interpolates linearly between
  days. v0.9's own record ends in December 2022 and is left unchanged.

## Soils

- SSURGO via USDA NRCS Soil Data Access, queried with the parcel polygon
  (`Soils/SSURGO_query_parcel_polygon_wgs84.wkt`), 6 October 2026:
  - `SSURGO_parcel_mapunits.csv` — 77D Cheshire–Holyoke complex, 15–35 %
    slopes, very rocky (41 %); 6 Wilbraham and Menlo soils, 0–8 % slopes,
    extremely stony (31 %); 77C Cheshire–Holyoke complex, 3–15 % slopes,
    very rocky (28 %).
  - `SSURGO_parcel_horizons.csv` — horizons of the major components
    (texture, OM, bulk density, pH 1:1 water, CEC7/ECEC, sum of bases,
    Ksat, water retention, rock fragments).
  - `SSURGO_parcel_restrictions.csv` — Holyoke: basalt bedrock at 46 cm;
    Wilbraham and Menlo: dense till at 69 cm; Cheshire: none.
- Decision: **Cheshire** soil physics in every cell to begin with (well
  drained fine sandy loam, sand 55–62 %, clay 9.5 %, OM 3.5 % Ap / 1 % B,
  bulk density 1.30 / 1.53, pH 5.3, 12–22 % rock fragments, no restriction
  above 165 cm; the official series description puts bedrock below 6 ft
  and calls every horizon strongly acid). The soil column is, by
  decision, the base parameter file's default: 20 layers to 5 m with free
  drainage.
- **pH and OM per cell** from each cell's own map unit (dominant
  component, 0-30 cm): Cheshire pH 5.3, OM 2.7 % (646 cells); Wilbraham
  pH 5.5, OM 4.5 % (290 cells). pH is 1:1 water, so the model's
  extraction ratio is 1.
- **Exchanger (CEC and base saturation)**, one representative set for
  every cell, `Soils/Representative_soil_chemistry.csv`, derived by
  `Representative_Soil_Chemistry` (v2/interpolation): SSURGO gives
  Cheshire no CEC7 and an ECEC of 1.6 cmolc/kg below its own sum of bases
  (10), and SoilGrids has no base saturation, so it comes from the NRCS
  Kellogg laboratory data of the eight characterised Cheshire pedons
  (`Soils/KSSL_Cheshire_pedons_0-30cm.csv`, Soil Data Access lab tables;
  six within 13 km). The 0-30 cm layers within 0.3 of the cells' mean pH
  (5.36) that carry an effective CEC (bases + KCl Al) are pooled: three
  layers at pH 5.2-5.3 (85P0055, 97P0102). Effective CEC (at the soil's
  own pH, the model's `CEC_tot_in`) 6.8 cmolc/kg; Ca 0.709, Mg 0.216,
  K 0.019 of it.
- **SoilGrids** at the parcel (`Soils/SoilGrids_parcel_centre.csv`,
  v2.0 REST, 0-5/5-15/15-30 cm): pH 5.0/5.2/5.3, agreeing with SSURGO;
  CEC at pH 7 38.7/12.2/10.4 cmolc/kg (0-30 cm median 12.9, which at the
  pooled layers' effective/pH-7 ratio of 0.63 gives 8.1, the cross-check);
  organic carbon 269/50/22 g/kg, a forest-floor value at 0-5 cm (90 %
  range 27-528), so OM is SSURGO's, not SoilGrids'.
- The exchanger follows pH: the one forested Cheshire pedon (40A1239,
  7.9 km, pH 4.6) holds Al on 67-83 % of its sites and bases on 22 % of
  its 3.0 cmolc/kg. Were the site's pH revised to that, the set would
  change with it (`Representative_Soil_Chemistry` at the new pH).
- `Run_Hamden_Field_Inputs` (v2/interpolation) rebuilds the per-cell file
  from all of this; sub-plot data later fill its `measured_*` columns. The
  set passes the model's ERW handoff gate at both map units' pH.

## Vegetation and management

- `Vegetation/PFT_rows_Harvard_Forest_and_Fendt.csv` — the two rows taken
  from the T&C site-PFT lookup table of the legacy T&C-SMEW workflow
  (`Table_Sites_PFT_LOOKUP_latest.xlsx`, copied alongside):
  - trees (HIGH): Harvard Forest, MA (row 181, "BDT temperate":
    deciduous, aSE = 1, Vmax 75, rooting depth ZR95 1 m);
  - grass (LOW, beneath the trees and in the open): Fendt, DE (row 102,
    "C3 grass": aSE = 2, Vmax 50, ZR95 25 cm), a humid temperate meadow.
- Management: grass cut once a year on 1 September; no liming, fertiliser,
  irrigation or rock.
- Trees and grass are perennial. The stand in `PFT_Hamden_Deciduous_Grass`
  starts only a cold start: a dormant Harvard Forest stand over an
  established sward, its 140 gC m-2 of roots sized for a mown lawn's peak
  LAI of about 2.5 (a sward started from seed died within two years). From
  Phase 3 on, every year, the handoff and the transient start from the
  previous run's end state (`Perennial_Vegetation_*`).
  - The trees close their canopy within a decade and hold it: peak LAI
    4.8-5.0 through twenty-five coupled years, near Harvard Forest's
    measured 4.5-5. An open lawn settles at a peak LAI of about 1.2.
    Beneath the trees, the sward is shaded out.
  - Phase 1 has no biogeochemistry, so the plants take up no nutrients
    there, yet growth still draws on their reserves. Carried through it,
    the grass starved within a year and the trees within four, so each
    Phase-1 year starts from the file, as before.
- Heartwood turns over in 40 years (`Wm_H`), as in the legacy coupled
  T&C-BG-2D forest runs (37.6 and 40.1 years). The site table's 0, from
  water-carbon runs, let heartwood grow for ever and bury 2 g N m-2 a
  year (three times the deposition) and its base cations. The canopy then
  thinned after fifteen years, to LAI 1.4 by the forty-second.
- The trees grow without T&C's sink limitation (`OPT_EnvLimitGrowth = 0`),
  T&C's own default and the legacy T&C-BG-2D forest runs' setting.
  - The v0.9 plot-scale driver switches it on for every site. With it on,
    surplus NPP goes to the carbohydrate reserve past its cap of 2/3 of
    the sapwood.
  - Carried from year to year, the reserve then grew from 0.8 to 2 kg C
    m-2 while the canopy shrank.

## Land cover

- `Land cover/Annual_NLCD_LndCov_2020_CU_C1V2.zip` — USGS Annual NLCD 2020,
  CONUS, 30 m (the full 1.4 GB raster; documentation copy only).
- `Land cover/NLCD_2020_Hamden_Blue_Trail_window.tif` — the parcel window
  plus 300 m (28 × 29 pixels), in the source CRS (Albers equal area on
  WGS84), read without loading the CONUS raster; checked identical to the
  source pixels and projection.
- `Land cover/NLCD_2020_parcel_class_shares.csv` — exact area-weighted
  shares of the parcel: deciduous forest 80.3 %, woody wetlands 17.7 %,
  developed open space 2.0 %.
- `Land cover/NLCD_class_legend.csv` — NLCD class codes and names.
- Per-cell cover (tree crown, open grass, urban, water, bare, rock) comes
  from this map by exact area overlap of each grid cell with every pixel
  it touches (`Cell_Land_Cover`, classes split by `NLCD_Cover_Class_Table`),
  or is set by hand (`Cell_Land_Cover(...,'manual',...)`). The parcel:
  tree 0.980, grass 0.018, urban 0.002; the grid's kept squares, inside
  the parcel: 0.982, 0.016, 0.002 (the 1 % of the parcel in dropped
  boundary squares is grassier than the rest).

## DEM

- `DEM/USGS 1 Meter 18 x67y459 CT_Statewide_C16.csv` — the National Map
  product record for the 1 m 3DEP tile covering the parcel (UTM 18N NAD83,
  easting 670–680 km, northing 4580–4590 km; CT statewide 2016 lidar).
- `DEM/USGS_1M_18_x67y459_CT_Statewide_C16.tif` — the tile (406 MB as
  served now; NAD83 / UTM 18N, 1 m, 10012 x 10012), documentation copy.
- `DEM/USGS_1M_18_x67y459_CT_Statewide_C16_parcel_window.tif` — the
  parcel and 200 m around it (584 x 609 m, 1.4 MB), also in
  `v2/Spatial inputs/Hamden Blue Trail/DEM/` (not versioned in git). It
  gives the cells exactly the tile's terrain (elevation 82.7-101.5 m,
  slope to 0.49, median 0.11); the margin only lets the edge cells'
  slope and wetness index see real terrain. The simulation is the whole
  field, with no buffer region.

## Model configuration (Bouchet)

- Grid: 5 m squares, each kept when more than half of it lies in the
  parcel, on the lattice (rotation 19.31 degrees, along the parcel's
  115 m and 71 m sides, and offset) that minimises the misrepresented
  area: 936 squares, 274 m2 of outside land taken in and 273 m2 of the
  parcel left out, 2.3 % of the parcel (the 95 % coverage rule: 1,224
  squares of 4.48 m, 1,179 m2 taken in, 5.0 %). 5 m is the coarsest
  width from which the mean slope, its spread, the wetness index's spread
  and the relief change by less than 5 % at every refinement down to 1 m
  (`Run_Topographic_Grid_Convergence`). The rotation turns only the
  lattice: the cells' terrain comes from the north-up DEM at their true
  centres, and at 0, 9.7, 19.3 and 45 degrees their slope, wetness index,
  aspect and face gradients agree within 1 %. The cluster builds the same
  grid (`TC_CELL_M=5`, `TC_MIN_COVERAGE=majority`); the per-cell file is
  tied to it by cell id and centre.
- Per-cell file:
  `Inputs/Field_Spinup/Hamden/Cell_Agronomic_Hamden_Blue_Trail_Parcel_polygon_01.csv`
  (Cheshire texture, each map unit's pH and OM, the representative
  exchanger, NLCD cover), with
  `Cell_Review_SSURGO_and_NLCD.csv` beside it (every cell's own SSURGO map
  unit and NLCD class areas); both written by `Run_Hamden_Field_Inputs`.
- Site rows over the borrowed workbook case 5, in the same folder:
  `Hamden_MOD_columns.csv` (location, the vegetation file
  `PFT_Hamden_Deciduous_Grass`, the CO2 and deposition records, the
  forcing and windows: ring 2015–2018, the 2019 handoff extension,
  handoff at 2020-01-01 for the simulation and validation runs,
  2020-01-01 to 2026-09-29; the soil layers and rain chemistry are left
  to the borrowed row, its default column and none) and
  `Hamden_ERW_columns.csv` (no rock, the exchanger seeded from each cell's
  base saturation, pH by a 1:1 extract, the exchange selectivities of a
  sandy loam -- Cheshire's USDA texture class, the key the ERW chemistry
  tabulates them by).
- The vegetation file sets the forcing reference height 2 m above the
  24 m canopy (26 m), loads the event deposition and empties
  `FertN_years`: the borrowed maize-soy row lists two years in three
  there, and the plot-scale driver every spin-up stage runs would refill
  202 kg N ha⁻¹ of fertiliser in each of them.
- `v2/cluster/field_spinup_hamden.env.example` submits the flux-tracking
  spin-up (one class seed, 94 tracked and replayed cells, 842 from the
  forest); its header lists what must exist first (the forcing and the
  DEM window).
- **Warm start** once the site observations arrive: put them in the
  per-cell file's `measured_*` columns (or a new representative set), set
  `TC_WARM_START_STATES` to the first run's `CellBGStates.mat`, give a new
  output folder and submit again; every cell re-converges from its own
  converged state (Phase 3 only).
- Open decisions: whether the 2-D transient takes its Na and
  strong-anion inputs from CT15 rather than the steady closure
  (`keyword_add`); the base parameter file's exchange constants are the
  Energy Farm's calibration.
