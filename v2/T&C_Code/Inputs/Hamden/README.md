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
| `T&C-SMEW-N v1.0/v2/T&C_Code/Inputs/` | model-ready records the loaders look for: `Rain_Chem_Hamden.mat`, `Ca_Data_NOAA_MLO_daily.mat` |
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
- **Still to download** (Meteo/NLDAS-2): `Tair` (2 m air temperature, K),
  `Qair` (2 m specific humidity, kg/kg), `PSurf` (surface pressure, Pa),
  `Wind_E` and `Wind_N` (10 m wind, m/s), `Rainf` (hourly precipitation,
  kg m⁻² = mm) and `SWdown` (incoming shortwave, W m⁻²); `LWdown` optional,
  for checking. Hours in UTC. Suggested period: 2014-01-01 00 UTC to the
  latest available, giving ten complete years for the spin-up cycle plus
  the recent years.
- **Downloading.** Either route works with the converter:
  - GES DISC *Hydrology Data Rods* (one text time series per variable for
    the point 41.437 N, 72.866 W, the quickest for a single cell), or
    Giovanni's area-averaged time series (CSV); or
  - the hourly `NLDAS_FORA0125_H.A<date>.<HH>00.020.nc` files from the GES
    DISC subsetter or OPeNDAP, any small window around the site.
  Put the files of one route in `Meteo/NLDAS-2` (one format per folder).
- **Converting**, in MATLAB from `v2/T&C_Code`:
  `Build_TC_Forcing_From_NLDAS2('Inputs/Hamden/Meteo/NLDAS-2','Inputs/Data_Hamden.mat',41.437,-72.866,97,-5)`.
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
  terrain. The spin-up needs 2015-01-01 to 2025-12-31 (local).

## Atmospheric deposition

- `Deposition/NTN-ct15-w-s-mg.csv` — NADP National Trends Network site CT15
  (Abington, CT), weekly wet-only samples, concentrations in mg/L
  (−9 = missing), 2014-12-30 to 2026-05-19, 594 weeks.
- `Deposition/CT15_wet_deposition_summary.csv` and `Rain_Chem_Hamden.mat` —
  derived with `Inputs/Preprocess/Extract_NTN_Rain_Chemistry.m` by NADP's
  convention: precipitation-weighted mean concentration of the 424 valid
  samples (validity code `w`, all ions and pH present) × all gauge
  precipitation (1,415 mm/yr; 82 % of it sampled validly).
  Rain pH 5.20. Wet deposition (kg ha⁻¹ yr⁻¹): N 3.45 (NH₄-N + NO₃-N),
  S 1.28, Ca 0.74, Mg 0.36, K 0.33, Na 2.57. Strong-acid anions
  (SO₄ + Cl) 14.7 µeq/L = 56.9 µmolc m⁻² d⁻¹.
- In the model each element enters once: N, P, K, Ca, Mg through the
  biogeochemistry's deposition terms (these CT15 values replace the
  built-in global maps, which give 13.4 kg N ha⁻¹ yr⁻¹ at Hamden, wet plus
  dry, older and coarse); Na and the strong anions through the ERW
  chemistry; nitrate is not counted among the strong anions because it
  is already N deposition. P is not measured by NADP and stays on the
  global map.
- NTN is **wet deposition only**. Dry deposition (EPA TDep maps, or CASTNET
  site ABT147 at Abington) would add to it without overlap; not included.

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
  and calls every horizon strongly acid). The soil column is 2 m in 14
  layers with free drainage.
- CEC, base saturation and initial pH are every cell's baselines from one
  **representative set**, `Soils/Representative_soil_chemistry.csv`, to be
  replaced by the sub-plot data. pH is Cheshire's 5.3 (1:1 water, so the
  model's extraction ratio is 1). **CEC and base saturation are still to
  be chosen**: SSURGO gives Cheshire no CEC7 and an ECEC of 1.6 cmolc/kg
  that is below its own sum of bases (10), so neither can be taken from
  it. CEC is the exchanger's capacity at the soil's own pH (the model's
  `CEC_tot_in`); base saturation is Ca, Mg and K as fractions of it. Once
  set, `Run_Hamden_Field_Inputs` (v2/interpolation) rewrites the per-cell
  file. The spin-up refuses the grid until every cell has a base
  saturation: the site seeds its exchanger from it, and the borrowed
  Kantola row would otherwise impose its Mollisol's 80 %.

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
  or is set by hand (`Cell_Land_Cover(...,'manual',...)`). On the
  calibrated grid: tree 0.980, grass 0.018, urban 0.002 (area-weighted),
  matching the parcel's shares.

## DEM

- `DEM/USGS 1 Meter 18 x67y459 CT_Statewide_C16.csv` — the National Map
  product record for the 1 m 3DEP tile covering the parcel (UTM 18N NAD83,
  easting 670–680 km, northing 4580–4590 km; CT statewide 2016 lidar).
- **Still needed: the tile itself**, `USGS_1M_18_x67y459_CT_Statewide_C16.tif`
  (444.8 MB), from the URL in the record; it goes in
  `v2/Spatial inputs/Hamden Blue Trail/DEM/` (not versioned in git).

## Model configuration (Bouchet)

- Grid: the calibrated grid (95 % coverage rule), 4.48 m, 1,224 cells.
- Per-cell file:
  `Inputs/Field_Spinup/Hamden_calibrated_grid/Cell_Agronomic_Hamden_Blue_Trail_Parcel_polygon_01.csv`
  (Cheshire texture and OM, the representative chemistry, NLCD cover), with
  `Cell_Review_SSURGO_and_NLCD.csv` beside it (every cell's own SSURGO map
  unit and NLCD class areas); both written by `Run_Hamden_Field_Inputs`.
- Site rows over the borrowed workbook case 5, in the same folder:
  `Hamden_MOD_columns.csv` (location, the 2 m soil column, the vegetation
  file `PFT_Hamden_Deciduous_Grass`, the CO2 and deposition records, the
  forcing and windows: ring 2015–2023, the 2024 handoff extension, handoff
  at 2025-01-01, a placeholder until the study period is fixed) and
  `Hamden_ERW_columns.csv` (no rock, the exchanger seeded from each cell's
  base saturation, pH by a 1:1 extract).
- The vegetation file sets the forcing reference height 2 m above the
  24 m canopy (26 m).
- `v2/cluster/field_spinup_hamden.env.example` submits the flux-tracking
  spin-up; its header lists what must exist first (the forcing, the DEM
  tile, the representative chemistry).
- **Warm start** once the site observations arrive: put them in the
  per-cell file's `measured_*` columns (or a new representative set), set
  `TC_WARM_START_STATES` to the first run's `CellBGStates.mat`, give a new
  output folder and submit again; every cell re-converges from its own
  converged state (Phase 3 only).
- Open decisions: the representative CEC and base saturation; whether the
  2-D transient takes its Na and strong-anion inputs from CT15 rather than
  the steady closure (`keyword_add`); the base parameter file's exchange
  constants are the Energy Farm's calibration.
