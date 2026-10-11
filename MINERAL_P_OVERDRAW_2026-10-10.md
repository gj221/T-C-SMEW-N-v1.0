# A single-precision overdraw of mineral P, 10 October 2026

One change to `BIOGEOCHEMISTRY_DYNAMIC3.m`. It comes with a test, `test_p_immobilization_exhaustion.m` (registered in `run_field_grid_tests`), and that test's fixture, `Tests/Fixtures/p_immobilization_hamden_20170723.mat`.

## What happened

Hamden's class seed under the layer's leaching (91c7aee re-applied, scratch only) stopped in Phase 3 loop 11, on 23 July 2017:

    Biogeochemistry pool 43 is -2.81845e-09 g m^-2 after the daily update.

A standalone daily-forced Phase 2 started from the same seed's Phase 2 state stopped the same way in its thirteenth year, at -2.72e-9 g m^-2. In both runs, fungi regrowing on a 100 kg C m^-2 POC pile wanted 0.061-0.064 g P m^-2 d^-1. The mineral pool and its day's inputs held 4e-5 g P m^-2.

The immobilisation limit was right: `Pmineral_avail` already holds uptake, leaching and the sorption exchange. What failed was its precision.

- `B_IO.DepP` comes from the global deposition map (`All_deposition_data.mat`), which holds P in single precision.
- `Deposition_Day` casts its record to double. `BIOGEO_UNIT` passes `B_IO.DepP` on unchanged.
- So `Incoming_P`, `P_external_base`, `Pmineral_avail` and the whole demand split after it are single.
- `P_mineral_transfer = PMIN_IN+Porg_imm+P_imm_deficit` equals `-Pmineral_imm` only to single precision. On that day the residue is 4e-8 of the demand: 6.5e-15 g g^-1, or 2.8e-9 g m^-2.
- When the microbes empty the pool, that residue is an overdraw. `Biogeo_Pool_Update` accepts 16 eps of the pool per gram of soil, 1.5e-9 g m^-2.

The fixture day reproduces the stop exactly. From that day, the old code overdraws for many states: 133 of 301 fungal biomasses from 150 to 900 g C m^-2, and 87 of 201 pools of 1e-6 to 0.1 g P m^-2 emptied in a day. The worst is 1.6e-6 g m^-2. With the deposition in double, the fixture day does not overdraw.

## The change

When the day's net microbial demand empties the mineral pool (`-PMIN_IN >= Pmineral_avail`), the pool's budget for that day is redone in double before the split. That budget is `Incoming_P`, the three exchange fluxes, `P_external_base`, `P_internal_solution` and `Pmineral_avail`. The split, `dB(43:46)` and `Min_P` are then double. The pool ends within a few 1e-18 g m^-2 of zero on the failing days (3e-15 at worst in the sweeps), which the update's band sets to zero. The day's phosphorus closes to 1e-13 g m^-2; before, it closed to 4e-9.

Every other day takes exactly the path it took before, so the numbers do not change:

- Days the demand does not empty the pool are untouched.
- Days with double inputs give the same numbers, because the same expressions are evaluated on the same values.

v0.9's Kantola never empties its mineral P. The pool holds 184.5 g m^-2 in the 10-day baseline and never falls below 181 g m^-2 in the 5-year Base_1D runs 10, 11 and 12. Its EQUIV_1D_TMP 10-day baseline, built from this code and from a16c1b0, is identical in every array (see Proof). `test_2d_driver_equivalence` passes 26 of 26.

## Not changed: the exchange margin

The same single-precision `base` reaches `limit_internal_fluxes`. There, the rounding margin `16*eps(max(1,abs(available)))` becomes 16 eps of a single, 1.9e-6 g g^-1, which is 0.82 g P m^-2 at Hamden.

- Sorption out of a mineral pool holding less than that is blocked, so the exchange then runs one way only. Hamden's column-leaching seed is below the margin on 1,865 of 2,914 days (its last eight Phase 3 loops; mineral P 0.27-1.2 g m^-2). Kantola and the Energy Farm never are.
- Desorption from secondary P stops at that level. In the layer runs, Hamden's pool 45 sits at exactly 0.822766 g m^-2.

On the fixture day, the deposition in double raises the microbes' mineral P supply from 4.0e-5 to 1.0e-4 g m^-2 d^-1. The cure for both effects is `DepP=double(B_IO.DepP)` in `BIOGEO_UNIT`, as `Deposition_Day` already does. It would change v0.9's numbers, because Kantola's deposition P is single too. So it is v0.9's decision and is not made here.

## Proof

All runs were local (R2026a, one computational thread), in scratch copies of a16c1b0 with and without this change.

- **Kantola EQUIV_1D_TMP, 10 days** (`TC_CASE=5`, `TC_NDAYS=10`). Of 803 variables, 13 differ: ten path strings (one is a `tempname`), `Computational_Time`, `equivalence_fingerprint`, and `ERW_data_ittn`, whose differing fields are function handles that capture paths. A deep comparison of 4,353 leaves, skipping 57 function handles, finds every array identical. This is the same set as the 2026-10-09 revert proof.
- **`test_2d_driver_equivalence`**: 26 passed, 0 failed.
- **Unit tests**: `test_p_immobilization_exhaustion` passes; on a16c1b0 it stops with the production error. `test_pool_overdraw_guards` passes 44 of 44, and `test_biogeo_pool_exhaustion`, `test_no3_immobilization_replay` and `test_invalid_state_guards` pass.
- **The standalone Phase 2** runs its 200 years. The pile burns off and the decomposers settle at 13.4 (bacteria) and 47.2 (fungi) g C m^-2. Before, it had 1,214 days of pool 43 below zero in 13 years, down to -2.7e-9 g m^-2. Now every such day is at most 2.2e-16 g m^-2 below zero: rounding in double.
- **The class seed** continues from its loop-10 checkpoint through loops 11 to 14. It tracks the scratch run that read small overdraws as empty, to within 1 % by loop 14: fungi average 2,833 g C m^-2 in loop 13, and by loop 14 the pile averages 1.2 kg C m^-2. Its days below zero are at most 3.1e-16 g m^-2.
- **Double inputs**: on 106 states, the two failing days and the sweeps with the deposition in double, the old and new code give bit-identical outputs.

## Using it

- Mineral P is emptied only in P-starved soils. Kantola, the Energy Farm and Hamden's column-leaching seeds never empty it (minimum 181, 180 and 0.27 g m^-2), so their results do not change. Only the layer-leaching prototypes do, and none of them is committed.
- The edit changes `ERW_Equivalence_Fingerprint`, which hashes `T&C_Code/*.m`. So the EQUIV_1D_TMP baseline must be rebuilt (`regen_equiv_1d.sbatch`, cpumodel:9575f) before the next campaign validation, as after any `T&C_Code` edit.
- The new test and fixture fall under the campaign stamp's `Tests` globs.
