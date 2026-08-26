# Post-legacy equation and mechanism audit

## Scope and acceptance rule

Legacy reference:

`C:\Users\Gregory\OneDrive\Yale\TC-SMEW-N\Other models\T&C-SMEW`

Audited model:

`C:\Users\Gregory\OneDrive\Yale\TC-SMEW-N\T&C-SMEW-N v0.9`

This audit covers equations, state transitions, and boundary couplings that are
absent from the legacy model or materially changed in the current model. A
mechanism is accepted only when all of the following hold:

1. The equation has a physically interpretable state and boundary.
2. Every term has consistent dimensions after the code's documented unit
   conversions.
3. Internal transfers cancel in the relevant element or energy ledger.
4. Source-pool availability is enforced before the state is updated.
5. A focused regression checks the equation, its limiting cases, or its
   conservation identity.
6. Empirical parameters are either traceable to literature or explicitly
   identified as model priors.

Passing a regression does not establish observational validity. That requires
the 92-day initialization gate and fresh full control, basalt, and calcite
trajectories.

## Highest-risk failure points

1. A diagnostic bulk-pH observation operator could be mistaken for pore-water
   chemistry and fed back into mineral kinetics. Bulk pH must remain diagnostic;
   pore-water pH drives speciation and reaction rates.
2. Hourly operator-split lateral transport can under-resolve a rapid reversal of
   water flux. It is conservative and refinement-convergent, but it is only
   first-order accurate and must not be used with an enlarged coupling step
   without repeating the refinement test.
3. Concentration-based saturation ratios assume `gamma_i = 1`. This is a stated
   dilute-solution approximation, not a general activity model. High ionic
   strength or concentrated fertilizer microsites require an activity model
   before saturation indices can be interpreted thermodynamically.

## Equation audit

### 1. ERW water volume and aqueous inventory

For an active soil depth `Zbio` in mm,

```text
W = Zbio * max(0, theta - theta_hyg)
M_i = C_i * W
```

`W` is L m^-2 because 1 mm water equals 1 L m^-2. `C_i` is umol L^-1
and `M_i` is umol m^-2. The previous layer-weighted expression multiplied
layer thickness twice and changed with mesh refinement. It was replaced by
`Biogeo_Mobile_Water_Volume`.

Status: corrected and tested.

### 2. Hydrogen concentration and carbonate speciation

```text
H = 10^(6-pH)                                      [umol L^-1]
h = H / 1e6                                        [mol L^-1]
CO2aq = DIC * h^2 / (h^2 + K1*h + K1*K2)
HCO3  = DIC * K1*h / (h^2 + K1*h + K1*K2)
CO3   = DIC * K1*K2 / (h^2 + K1*h + K1*K2)
```

The acid-base ratios are dimensionless when `h`, `K1`, and `K2` use molar
units. Species sum exactly to DIC. The retained code uses concentration ratios
with `gamma_i = 1`.

Status: mechanistically correct under the documented ideal-dilute assumption.

### 3. Al hydrolysis and finite Al(OH)3 inventory

For free Al3+ and hydrogen in the same concentration basis,

```text
Al_total_aq = Al3 * (H^4 + b1*H^3 + b2*H^2 + b3*H + b4) / H^4
Al3_gibbsite = 10^logK * (H/1e6)^3 * 1e6
```

The current initializer inverts these expressions rather than treating free
Al3+ as total dissolved Al. Secondary `AlOH3_solid` is retained as a finite
state, so precipitation cannot create an unlimited Al source and dissolution
cannot overdraw the solid. Crystalline gibbsite uses `logK = 8.04`; amorphous
Al(OH)3 uses `logK = 10.8` only when selected explicitly.

Status: corrected and tested for both exchange-controlled and solid-controlled
initialization.

### 4. Gaines-Thomas exchange and base saturation

The exchange state is solved on charge-equivalent fractions. The complete
state satisfies

```text
f_H + f_Al + f_Na + f_K + f_Ca + f_Mg = 1
base_saturation = f_Na + f_K + f_Ca + f_Mg
exchange_acidity = f_H + f_Al
```

Selectivity equations use one internally consistent concentration basis. The
Al/Ca inversion now receives `f_Al`, not `f_Ca`. Cation totals used by the
T&C/ERW daily coupling include solution plus exchange inventory, so internal
solution-exchange transfer cancels.

Status: corrected and tested for fraction closure, cation closure, pH response,
and finite Al.

### 5. Bulk-pH observation operator

Pore-water pH is the chemistry state. Bulk pH is reconstructed from the same
finite aqueous and exchange totals with a calibrated initial diagnostic offset:

```text
pH_bulk = Invert_Bulk_pH(finite totals, exchange state) + offset_initial
```

The offset aligns the initial measured bulk pH without changing pore-water H,
DIC, alkalinity, mineral saturation, or microbial pH. The operator is suitable
for comparison with soil-slurry observations but is not a mechanistic slurry
extraction model because field DIC transfer, extractant ratio, and atmospheric
CO2 equilibration are not represented explicitly.

Status: numerically correct diagnostic operator. It must not drive chemistry.

### 6. Mineral saturation index

For a reaction with stoichiometric exponents `nu_i`,

```text
IAP = product(a_i^nu_i)
Omega = IAP / Ksp
SI = log10(Omega)
```

With `gamma_i = 1`, the code uses concentration divided by the 1 mol L^-1
standard state as activity. Corrected stoichiometry includes:

```text
albite:       NaAlSi3O8
anorthite:    CaAl2Si2O8
zoisite:      Ca2Al3Si3O12(OH)
epidote:      Ca2Al2FeSi3O12(OH)
hydroxyapatite: Ca5(PO4)3OH
```

Unsupported Fe-bearing or custom phases return `NaN` saturation rather than a
fabricated number. Their dissolution can retain a kinetic affinity of one only
when the caller explicitly treats them as irreversible primary minerals.

Status: stoichiometric defects corrected and tested against direct IAP
reconstruction. Saturation values remain ideal-solution indices.

### 7. Primary-mineral shrinking-particle kinetics

For each particle class,

```text
r = r0 * (M/M0)^(1/3)
A = A0 * (M/M0)^(2/3)
dM/dt = -A * k(T,pH) * affinity
```

The particle-size distribution is formed from CDF mass increments, not CDF
slopes. Primary mineral mass is finite and dissolution is limited by remaining
inventory before update. Native mineral units, molar conversion, and g m^-2
are kept separate.

Status: corrected and tested for particle distribution, molar-mass shape,
finite inventory, and dissolution direction.

### 8. Calcite and carbonate solids

Primary calcite only dissolves. Dissolved carbonate can precipitate into a
separate secondary `CaCO3` pool. Dissolution follows the Plummer-Wigley-Parkhurst
surface-reaction structure. Secondary calcite precipitation uses the Kirk et al.
soil-volume rate coefficient:

```text
alpha = 5e-10 mol dm^-3 soil s^-1
      = 43200 umol m^-3 soil d^-1
Rprecip = alpha * f(Omega), Omega > 1
```

Direct anhydrous MgCO3 precipitation is disabled at ambient field temperature;
its dissolution pathway remains available. Precipitation and dissolution are
limited so `Omega` cannot be crossed by more material than is chemically or
physically available.

Status: corrected and tested for density, DIC stoichiometry, event routing,
primary inventory, secondary separation, flux sign, and Ca closure.

### 9. T&C to ERW daily element coupling

The first retained ERW state is already synchronized to the first T&C state.
It is not charged with the raw seed difference again. For later retained days,

```text
DeltaX_j = X_TC(j) - X_TC(j-1)
rate_X = DeltaX_j / Delta_t_window
sum(rate_X * dt_window) = DeltaX_j
```

Each daily delta is applied once. A repeated final partial flush carries zero
new T&C delta. K, Ca, and Mg totals include solution plus exchange pools;
therefore exchange is an internal redistribution rather than a boundary flux.

Status: duplicate-final-delta and partial-window defects corrected. Focused
coupling and 1-D/2-D equivalence tests cover the boundary.

### 10. Gross nitrification pH response

```text
f_nit(pH) = clamp(0.56 + atan(pi*0.45*(pH-5))/pi, 0, 1)
R_nit = R_nit_potential * f_nit(pH) * other_legacy_limits
```

The factor is dimensionless, bounded, monotonic, and uses the same-day pore or
microbial pH supplied to the N-cycle routine. It does not use bulk pH.

Status: extracted, bounded, and tested across acidic to neutral soil.

### 11. Nitrification N2O fraction

```text
f_N2O_nit = min(1, 721.86 * exp(-2.387*pH))
N2O_nit = R_nit * f_N2O_nit
```

The fraction is dimensionless and cannot exceed gross nitrification. This is
the Val Martin et al. implementation used for the Energy Farm comparison.

Status: tested for bounds, monotonic pH response, and donor limitation.

### 12. Gross denitrification and N2:N2O partition

```text
f_den_pH = 0.001                              pH <= 4
           (pH - 4)/3                         4 < pH < 7
           1                                  pH >= 7
f_ratio_pH = 0.0016 * exp(1.006*pH)
R_den = min(NO3_available/dt, R_den_potential * f_den_pH)
ratio_N2_to_N2O = ratio_base * f_ratio_pH
N2O_den = R_den / (1 + ratio_N2_to_N2O)
N2_den  = R_den - N2O_den
```

The published piecewise pH response controls gross denitrification. The
published exponential response controls `ratio_N2_to_N2O`. N2 and N2O
sum to gross denitrification. Gross denitrification remains limited by the
legacy moisture, carbon, temperature, and nitrate controls. Product-ratio
NO3, heterotrophic-respiration and WFPS modifiers use the explicitly selected
Parton formulation and do not alter the gross NO3 sink.

Status: tested for pH response, nitrate limitation, partition closure, and
nonnegative gas fluxes.

### 13. Biological N fixation carbon payment

```text
BfixN = Cfix_today / Costfix(T)
```

Both numerator and cost use g C m^-2 d^-1 and g C per g N, so the result is
g N m^-2 d^-1. The current model uses the same-day fixation allocation
`Rexmy(3)`, not its 365-day moving average. This prevents fixation after the
carbon payment stops. `Costfix(T)` is inherited from legacy T&C and is not a
new current-model equation. The legacy model does not retain a separate
symbiont-respiration pool, which remains a structural limitation.

Status: current boundary corrected and tested. Legacy cost function retained.

### 14. Plant reserve and crown-area ledgers

Reserve storage is physical mass, never a sentinel. Initialization uses zero
unless an explicit reserve is supplied. Per-crown fluxes convert to plot-area
mass exactly once:

```text
F_plot = Ccrown * F_crown
Delta reserve = inputs - remobilization - export
```

Harvest export is accumulated in plot-area units. The state update writes only
the current retained row; it does not overwrite the previous row. Standing
leaf nutrient closure uses one concentration for prior storage, same-day
senescence input, litter loss, and export.

Status: legacy sentinel, crown-zero division, row-overwrite, and harvest-unit
defects corrected and tested.

### 15. Lateral Darcy water flow

For edge half-distances `d_a` and `d_b`,

```text
Q = A_face * (h_a-h_b) / (d_a/K_a + d_b/K_b)
Delta theta_c = sum(signed Q * dt) / (A_cell * Zlayer)
```

`Q` is m^3 h^-1. The perimeter is closed. Every face is shared by two cells,
so water gain and loss cancel exactly. Adaptive substeps limit donor and
receiver water without clipping volume. On a flat field, elevation head
cancels and the lateral gradient is matric head only.

Status: tested for direction, unequal-cell resistance, closure, bounds,
heterogeneous soils, and applied-volume identity.

### 16. Frozen-pore impedance

```text
F_ice = theta_ice / (theta_liquid-mobile + theta_ice)
K_frozen = K_liquid * 10^(-6*F_ice)
```

The dimensionless empirical form and `Omega = 6` follow Swenson et al. (2012)
and CLM. The previous current-model value `5` was unsupported and was replaced.
Liquid plus ice cannot exceed mobile pore volume.

Status: corrected and tested.

### 17. Lateral solute advection

```text
Delta M_i = Q * dt * C_i,donor
```

`Q*dt` is m^3 and `C_i` is mass or amount per m^3, so `Delta M_i` is mass or
amount. Export is limited by donor inventory. The beginning-of-step donor
concentration makes this a first-order operator split. A fixed-total-time test
shows monotonic convergence at 1, 2, 4, and 64 substeps.

Status: conservative and refinement-convergent. No lateral dispersion or
molecular diffusion is represented.

### 18. Lateral heat and freeze-thaw enthalpy

```text
H = cv*T - rho_i*Lf*theta_ice
dtheta_liq/dT = Lf * (1000*Cw) / (g*Tref)
Capp = dH/dT
     = cv + (rho_i*Lf + T*(cv_w-cv_i))*dtheta_liq/dT
E_adv = rho_w*cp_w*Q*dt*T_donor
E_cond = G_edge*(T_a-T_b)*dt
```

`Cw` is the water-retention derivative in mm^-1. `Lf/(g*Tref)` is m K^-1;
the factor 1000 converts it to mm K^-1. Thus `dtheta_liq/dT` is K^-1 and
`Capp` is J m^-3 K^-1. The solver advances physical enthalpy and inverts it
with safeguarded Newton/bisection. `Capp` is only the derivative and stability
capacity; an invalid value now raises an error instead of silently reverting to
sensible heat capacity. A fixed offset carries T&C's post-cryosuction phase
departure without creating latent heat.

Status: equation, derivative, water-phase closure, energy closure, stiff-edge
bounds, and off-equilibrium state are tested.

### 19. ERW lateral aqueous and soil-gas transport

Aqueous Na, Al species, residual anion, and DIC use the same donor-volume
advection as other solutes. Soil-gas CO2 uses diffusive edge exchange and a
piston boundary when lateral water changes air-filled pore volume. The closed
lateral diffusion term cancels. Piston exchange is an explicit atmospheric IC
boundary and is reported separately.

Status: tested for species transport, saturation response, concentration
bounds, closed diffusion, and atmospheric IC closure.

### 20. 2-D topology and parameter adapters

The topology partitions one square domain into main and margin cells while
retaining exact cell areas, face widths, half-cell distances, and a closed
outer perimeter. Parameter adapters map legacy field names without changing
values. They do not introduce a new physical equation.

Status: geometry and adapter contracts tested. The 60 m square with a 5 m
margin is a model domain, not the literal Kantola experimental plot geometry.

### 21. Prairie-margin parameterization

`PFT_Prairie_Margin.m` uses legacy photosynthesis, hydraulics, phenology,
allocation, turnover, and nutrient equations. It introduces parameter priors,
not new governing equations. Big bluestem rooting beyond 2 m is literature
plausible, but exact root percentiles, coneflower rooting, phenology thresholds,
hydraulic vulnerability, turnover, initial biomass, and tissue stoichiometry
are not Kantola measurements.

Status: experimental parameterization. It is not used by 1-D runs 10-12 and
must not be presented as site-validated without margin observations and a
separate sensitivity analysis.

## Dimensional analysis

| Variable | Expected SI Units | Actual Units provided in code |
|---|---:|---:|
| `Zbio` | m | mm, converted through 1 mm = 1 L m^-2 |
| `theta`, `theta_hyg`, exchange fractions | m^3 m^-3 or 1 | dimensionless |
| `W` active water | m^3 m^-2 | L m^-2, numerically equal to mm |
| aqueous `C_i` | mol m^-3 | umol L^-1, numerically equal to mmol m^-3 |
| areal aqueous `M_i` | mol m^-2 | umol m^-2 |
| `H`, cations, anions | mol m^-3 | umol L^-1 |
| `h`, `K1`, `K2` in carbonate ratios | mol m^-3 on one basis | mol L^-1 after `/1e6` conversion |
| DIC species | mol m^-3 | umol L^-1 |
| alkalinity | eq m^-3 | ueq L^-1 or umol charge L^-1 |
| CEC inventory | eq m^-2 | umol charge m^-2 after soil-depth conversion |
| mineral `Ksp`, `IAP`, `Omega` | dimensionless standard-state ratio | concentration activities with `gamma_i=1` |
| primary mineral mass | kg m^-2 | native mineral amount plus g m^-2 diagnostics |
| mineral amount | mol m^-2 | umol m^-2 internally |
| mineral surface area | m^2 m^-2 | m^2 mineral m^-2 ground |
| mineral rate constant | mol m^-2 mineral s^-1 or source-specific equivalent | converted to umol m^-2 ground d^-1 |
| calcite precipitation `alpha` | mol m^-3 soil s^-1 | 5e-10 mol dm^-3 soil s^-1, converted to 43200 umol m^-3 d^-1 |
| T&C daily `DeltaX` | kg or mol m^-2 d^-1 over retained day | g m^-2 daily deltas converted to umol m^-2 |
| N pools | kg N m^-2 | g N m^-2 |
| N process rates | kg N m^-2 s^-1 | g N m^-2 d^-1 |
| reported N2O event flux | kg N m^-2 s^-1 | kg N ha^-1 d^-1 after diagnostic conversion |
| fixation carbon allocation | kg C m^-2 s^-1 | g C m^-2 d^-1 |
| fixation cost | kg C kg^-1 N | g C g^-1 N |
| crown fraction | 1 | dimensionless plot-area fraction |
| Darcy head `h` | m | mm, divided by 1000 in edge gradient |
| hydraulic conductivity `K` | m s^-1 | mm h^-1, divided by 1000 for m h^-1 |
| lateral face area | m^2 | m^2 |
| lateral water flux `Q` | m^3 s^-1 | m^3 h^-1 |
| lateral solute concentration | mol m^-3 or kg m^-3 | pool-specific amount m^-3 |
| lateral solute transfer | mol or kg per step | pool-specific amount per cell per step |
| temperature `T` | K difference | degC relative to 0 degC; increments equal K |
| `cv`, `Capp` | J m^-3 K^-1 | J m^-3 K^-1 |
| `Cw` retention derivative | m^-1 | mm^-1 |
| `dtheta_liq/dT` | K^-1 | K^-1 after explicit mm conversion |
| thermal conductivity `lambda` | W m^-1 K^-1 | W m^-1 K^-1 |
| edge thermal conductance `G` | W K^-1 | W K^-1 |
| cell enthalpy density `H` | J m^-3 | J m^-3 |
| advective/conductive heat | J | J |
| soil-gas CO2 | mol m^-3 air | umol L^-1 air |
| gas diffusion or piston transfer | mol per step | umol per cell per step |

No unresolved dimensional mismatch remains in the accepted post-legacy
equations. The concentration-as-activity approximation and empirical prairie
parameters are scientific scope limitations, not unit mismatches.

## Boundary conditions and assumptions

- Vertical T&C hydrology, energy, vegetation, and legacy soil C-N-P processes
  remain the base model.
- The 1-D domain has no lateral exchange. The 2-D perimeter is closed to
  lateral water, solute, and conductive heat unless an explicit atmosphere
  term is reported for soil gas.
- Pore-water chemistry uses `gamma_i = 1`.
- Bulk pH is diagnostic and must never replace pore pH in reaction equations.
- Primary mineral inventories are finite. Secondary carbonates and Al(OH)3 are
  distinct retained pools.
- Internal speciation, exchange, sorption, mineral transfer, biological
  immobilization, and plant-litter transfer cancel in element ledgers.
- Fertilizer, deposition, fixation, gas loss, leakage, lateral transfer,
  amendment, harvest, and export are explicit boundaries.
- A state clamp is not accepted as conservation. Source availability is
  enforced before update, and infeasible chemistry raises an error.

## Test-to-mechanism map

| Mechanism | Focused tests | Acceptance identity |
|---|---|---|
| Active water | `test_biogeo_mobile_water_volume` | depth scaling and nonnegative volume |
| Appelo exchange and Al | `test_appelo_*`, `test_kantola_seed_charge` | charge, fraction, and finite Al closure |
| Bulk/pore pH | `test_pH_conversion`, `test_kantola_bulk_pH_calibration` | initialization and monotonic observation operator |
| Saturation and mineral kinetics | `test_erw_mineral_kinetics`, `test_erw_particle_distribution` | IAP stoichiometry and finite inventory |
| Calcite | `test_calcite_*` | DIC/Ca stoichiometry, event routing, primary/secondary separation |
| ERW/T&C coupling | `test_erw_exchange_pool_coupling`, `test_2d_driver_equivalence` | daily delta exactly once and 1-D/2-D identity |
| N transformations | `test_nitrification_*`, `test_denitrification_*` | pH response, substrate limit, N2+N2O closure |
| Plant reserves | `test_crop_*`, `test_plant_*`, `test_biological_fixation_carbon_allocation` | plot-area reserve and export closure |
| Pool availability | `test_pool_overdraw_guards`, `test_no3_immobilization_replay` | no donor overdraw or mass-creating clamp |
| Lateral flow and solutes | `test_lateral_routing` | water/solute closure and timestep convergence |
| Lateral heat and ERW gas | `test_lateral_energy_erw` | physical enthalpy, phase, IC, and elemental closure |
| Initialization boundaries | `test_erw_initial_boundary_defaults`, `test_erw_equilibrium_initialization` | no arbitrary source or unconstrained seed |
| Adapters and resume safety | `test_model_adapter_contracts`, `test_checkpoint_*` | exact mapping and source fingerprint coverage |

## Primary literature used

- Appelo, 1994, cation exchange and transport:
  https://agupubs.onlinelibrary.wiley.com/doi/10.1029/94WR01048
- Palandri and Kharaka, 2004, mineral kinetic constants:
  https://pubs.usgs.gov/publication/ofr20041068
- Plummer, Wigley, and Parkhurst, 1978, calcite dissolution kinetics:
  https://ajsonline.org/article/59968-the-kinetics-of-calcite-dissolution-in-co-2-water-systems-at-5-degrees-to-60-degrees-c-and-0-0-to-1-0-atm-co-2.pdf
- Kirk et al., 2015, reactive-transport model of soil calcite precipitation:
  https://doi.org/10.1016/j.gca.2015.05.017
- May et al., 1979, gibbsite equilibrium:
  https://doi.org/10.1016/0016-7037(79)90224-2
- Kanzaki et al., 2024, SCEPTER soil weathering chemistry:
  https://gmd.copernicus.org/articles/17/4515/2024/gmd-17-4515-2024.pdf
- Parton et al., 1996, soil N gas process representation:
  https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/96GB01455
- Val Martin et al., 2023, Energy Farm N-cycle and N2O response:
  https://gmd.copernicus.org/articles/16/5783/2023/gmd-16-5783-2023.pdf
- Blanc-Betes et al., enhanced weathering and N2O field evidence:
  https://eprints.whiterose.ac.uk/id/eprint/169796/1/gcbb.12757.pdf
- Dall'Amico et al., 2011, conservative frozen-soil enthalpy:
  https://tc.copernicus.org/articles/5/469/2011/tc-5-469-2011.pdf
- Swenson et al., 2012, frozen-pore hydraulic impedance:
  https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2012MS000165
- Nippert et al., 2012, tallgrass root depth and water use:
  https://doi.org/10.1007/s11104-011-1112-4

## Evidence still required before production acceptance

1. Four clean executions of the post-legacy equation test set from one source
   state.
2. A cold-start 92-day control/basalt/calcite gate with treatment equality,
   all-pool finiteness, exchange/charge/element closure, and a complete N
   ledger.
3. Fresh full 1-D trajectories crossing application day 314, with matched
   audits at days 180, 314, 400, and 679.
4. Comparison of bulk pH, base saturation, and exchangeable acidity with
   Beerling et al.; N2O and mineral N with Val Martin et al.; and GPP, ecosystem
   respiration, NEE, and biomass with Kantola et al.
5. A separate calibration and validation exercise before the prairie-margin
   parameterization is treated as site-specific evidence.
