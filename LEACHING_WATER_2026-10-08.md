# The water the biogeochemical layer leaches with (2026-10-08)

Decision (Greg, 2026-10-09): use T&C-BG's legacy leaching, as below, in
every model iteration, with no option. This note records what the legacy
formulation does and the alternative that was tested and not adopted.

## The legacy formulation

`Biogeo_Leakage` leaches every dissolved pool of the biogeochemical layer
as

```text
LEAK_X = min(B_X, a_X * B_X * min(1, Lk/V))          [g m-2 d-1]
```

with `Lk` the day's leakage out of the bottom of the whole soil column
(`Leakage_Bottom`) and `V` = `VSUM`, the water of the whole column
(`Biogeo_environment`: `Vtot = sum(V)`). The ERW leaches its own pools
(Na, the residual anion, Al, DIC) the same way, `f_Lw_bio = Lk/VSUM` in
both drivers, and under `keyword_add = 0` caps that at the column's
drainage times the pore-water concentration.

This is T&C-BG's own formulation: the original code (`TeC_Source_Code`)
and T&C-BG-2D (Lian et al. 2025) call `Biogeo_Leakage(Ptm1,Lk,V,...)` with
`V = sum(V)`, under the comment "for simplicity leakage only at the base
of soil column". v0.9 pins it in `test_biogeo_mobile_water_volume`
("BIOGEO_UNIT must retain the legacy bottom-leakage water boundary").
`POST_LEGACY_EQUATION_AUDIT.md` section 9b and
`LEACHING_CLOSURE_REVIEW_2026-10-01.md` write the formula as `Lk/VT`; the
code uses `Lk/VSUM`.

The pools are held in the 300 mm layer, and passive uptake (`VT`), the
ERW's pore water (`cov_moll = Zbio*theta_bio`) and the 2-D lateral
transport (`Biogeo_Lateral_Concentration`) treat them so. Only the
vertical leaching spreads them over the whole column: they leave at their
concentration divided by `VSUM/(Zbio*theta_bio)`, carried by the column's
bottom drainage, and the steady-state solution is more concentrated by
that factor.

## Size at Hamden

The parent session's 200-day plot-scale test (1 January to 18 July 2020,
class-seed cell 1, default 20 layers to 5 m, `Zbio` 300 mm):

| | mean | range |
|---|---|---|
| `VSUM`, whole column | 1007 mm | 771-1172 |
| layer's water, `Zbio*theta_bio` | 62.6 mm | |
| layer's mobile water, `VT` | 47.4 mm | 0-82 |

Monthly water, mm:

| month | rain | `Lk`, column bottom | layer bottom, net | transpiration |
|---|---|---|---|---|
| Jan | 50 | 10 | 27 | 0 |
| Feb | 83 | 9 | 106 | 0 |
| Mar | 103 | 36 | 105 | 0 |
| Apr | 120 | 99 | 118 | 5 |
| May | 73 | 125 | 33 | 69 |
| Jun | 46 | 67 | -23 (up) | 117 |
| Jul (18 d) | 36 | 26 | -13 (up) | 62 |
| total | 512 | 373 | 353 | 254 |

The 5 m column delays its drainage by one to two months: it drains
little while the layer is flushed in winter and early spring, and keeps
draining in June and July while water rises into the drying layer. So
`Lk/VT` would be wrong too: it would leach the summer layer while none
leaves it.

Over the 200 days the pools as coded lost 0.083 g m-2 of Ca, 0.24 Mg,
0.50 K and 0.25 g N m-2 of NO3. At the same concentrations, the layer's
own drainage over its own water would have carried 1.75 Ca, 5.2 Mg,
10.4 K and 3.2 g N m-2 NO3. Over the spin-up's four meteorological years
the daily leaching fraction is 0.23 % d-1 as coded and 2.4 % d-1 at the
layer's own drainage over its own water.

## The alternative tested

Commit 91c7aee (reverted by 28fc00e) leached the pools with the drainage
out of the bottom of the biogeochemical layer over the layer's own water
`Zbio*theta_bio`, in both drivers, the ERW and the spin-up. T&C does not
return the fluxes between its layers, so the drainage was taken each hour
as the part of T&C's own step balance (`CK1`) that falls in the layers the
biogeochemistry spans; over every layer that sum returns `Lk` (to 5e-15
mm h-1 on the Hamden run). It was inert for v0.9 inputs (Kantola
equivalence baseline identical in all shared arrays, 26/26) and the two
drivers agreed bit for bit with it on.

### Hamden: steady state (probe class seeds) and 200 days from each

| | legacy (`Lk/VSUM`) | layer's drainage and water |
|---|---|---|
| dissolved Ca / Mg / K, mmol L-1 | 0.076 / 0.377 / 0.469 | 0.054 / 0.036 / 0.065 |
| NO3, mmol L-1; NH4, g N m-2 | 0.225; 2.07 | 0.00075; 0.025 |
| handoff base saturation, pore pH, Al | 0.145, 4.91, 420 umol L-1 | 0.054, 4.96, 273 umol L-1 |
| tree N reserve, g N m-2; N limitation FNC_H | 32.6; 1 | -1.3; 0.65 |
| tree sapwood, fine roots, reserves, g C m-2 | 1099, 608, 769 | 259, 140, 280 |
| seeded strong anions, ueq L-1, day 0/50/100/200 | 2193 / 2647 / 2112 / 84 | 712 / 263 / 54 / 35 |
| solution NO3 peak, mmol L-1 | 1.32 (mid-June) | 0.002 |
| NO3 leached, N2O, kg N ha-1 (200 d) | 2.44, 0.077 | 0.06, 0.0001 |

With the layer's own leaching the deposition-only solution and its strong
anions follow the layer's water balance (the seeded anions drain with the
winter's water to the 30-35 ueq L-1 deposition sustains), but the
unfertilised stand starves of N: T&C's uptake is first order in the
dissolved pool (`OPT_UPT = 2`) with leaching taking priority, and its
constants were set with the diluted leaching, so the dormant season's
drainage takes the N that mineralisation releases.

### Energy Farm: first year (control, Base_1D options)

| 2016 | legacy | layer's drainage and water |
|---|---|---|
| NO3 leached, kg N ha-1 | 12.3 | 62.1 |
| N2O, kg N ha-1 | 1.69 | 1.07 |
| Ca+Mg+K leached, eq m-2 | 0.53 | 6.33 |
| August bulk pH (field about 6.07) | 6.02 | 5.97 |

Maize-soybean NO3 leaching measured at 50 cm on the Energy Farm averaged
40.4 kg N ha-1 yr-1 (McIsaac et al. 2010, J. Environ. Qual.). The layer's
leaching brings the model to that order but cuts N2O by a third, which
would need KNmaxd refitted; the five-year runs were stopped when the
legacy formulation was chosen.

## Re-test with the vegetation crowns fixed (2026-10-09)

Greg asked for the layer's leaching back if the crown fix (ce82785: a
crown without cover carries no vegetation) cured the starvation. It does
not. A cold Hamden class seed with 91c7aee re-applied, `leaching_water =
layer` and the crown fix loses its stand as before.

By the tenth Phase 3 year (in brackets, the cold class seed of 2026-10-08
under the column's leaching, before the crown fix):

| | layer + crown fix | (column) |
|---|---|---|
| peak tree LAI | 1.7 | (4.9) |
| fine roots, g C m-2 | 198 | (525) |
| sapwood, g C m-2 | 282 | (894) |
| reserve, g C m-2 | 204 | (621) |

The trees' N reserve is below zero from the third year. The open-grass
crown is empty and the warnings are gone, so the empty crown was not the
cause.

The N does not leave by leaching. NO3, NH4 and DON together lose 0.01-0.02
g N m-2 yr-1, against 0.25-0.5 in the converged column seed. The soil stops
releasing it.

Phase 2 (408 mean-forcing years under the layer's leaching) ends with the
decomposers gone (in brackets, the 2026-10-08 column seed's Phase 2):

- bacteria (B18) and saprotrophic fungi (B19): ~1.5e-9 g C m-2 (11 and 41);
- their enzymes (B14-B17): ~1.5e-9 g C m-2 (0.05-0.15);
- the DOC the decomposers live on (B12, B13): 5 and 12 g C m-2 (58 each);
- no NH4 or NO3 (0.16 and 0.02 g N m-2);
- undecomposed POC (B10): 100 kg C m-2 (1);
- SOM N (B26): 5.0 kg N m-2 (0.7).

Phase 3 then starts with net mineralisation of -0.15 to +0.16 g N m-2 yr-1,
against 21-26 in the converged column seed. The trees live on the
deposition: uptake 0.5-0.75 against 19-27 g N m-2 yr-1. In the tenth year
mineralisation jumps to 10.9: the decomposers grow back from 1e-9 and
start on the organic matter that piled up.

The leaching water reaches the stand through T&C-BG's decomposers. Their
dissolved substrate and mineral N leave the 300 mm layer at its own
drainage about ten times faster than at the column's (2.4 % against
0.23 % d-1), and under Phase 2's constant forcing they die out. Making the
layer's leaching usable means revisiting that, not the vegetation: the
leaching of the decomposers' DOC and mineral N, or Phase 2's forcing.
