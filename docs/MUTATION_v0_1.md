# Mutation testing v0.1: what the hank-limits suite catches and what it misses

Measured on branch `mutation-v0_1` (no tag, no release). Plan and predictions: `docs/MUTATION_PLAN_v0_1.md`, committed and pushed (`ccb7e38`) before any mutant was run. Runner: `scripts/mutation_v0_1.py`. Data: `out/mutation_v0_1.json` (every run's test outcomes and check summaries), `out/mutation_v0_1.csv` (one row per mutant × test group).

**Revision, 2026-10-22.** An independent review of commit `65f8f6d` found no discrepancy in the patches, the detection counts or the frozen checker, and asked for wording corrections. They are made here: the scope of the identity-test counts (sections 1 and 6); the split of the composite Ricardian test into its RA and HA assertions, and the description of the cohort and first-round checks (sections 1, 6 and 10); the classification of M21 (section 1); the execution record (sections 2 and 11). Section 8 adds a follow-up diagnosis of the C03 residual, run after and outside the registered measurement, and section 9 is updated to match. Patches, tolerances, test outcomes, the tables in sections 3–4 and `out/mutation_v0_1.{json,csv}` are unchanged. The unmodified suite still gives 41 passed, 2 xfailed (74 s, same environment as section 9).

## 1. Summary

- **Setup.** 25 mutants (small, plausible implementation errors in six categories) and 3 controls (changes that should not change any outcome). Each was applied to a copy of the model side (`models.py`, `core.py`), and the **full** pytest suite (43 tests) was run on the copy. Tolerances, check code, references and tests were not changed. The unmutated copy was run twice as the baseline (41 passed, 2 xfailed both times).
- **Suite-wide:** 21 of 25 mutants changed at least one test outcome. **Missed:** M12 (cached Jacobians passed to the T3 solve without matching options, F1 type), M21 (steady-state taxes omit debt service), M23 (TANK net-supply target missing the division by 1 − λ), M24 (Aiyagari borrowing limit −1 instead of 0). One of the 21 is incidental: M07 (quarterly δ) was flagged only by the top-of-grid mass check, because the larger capital stock pushed wealth towards the grid top. Without M07 the count is 20 of 25.
- **Controls:** no false alarms (0 of 3 changed any outcome).
- **Pre-registered predictions:** 18 mutants were predicted to be caught and 21 were. Of 170 firm group-level predictions, 163 held. The 7 that failed: M08 caught by T2, T5 and T6 (predicted missed); M09 caught by T1 (predicted missed); M12 missed by T3 (predicted caught); M13 caught by T5 (predicted missed); M15 caught by T3 (predicted missed). Section 7 explains each.
- **T5 and T6 as internal-consistency checks.** The review's concern holds for the identity parts:
  - The normally passing identity tests, that is the three IKC tests of T5 that pass on the baseline and the two decomposition tests of T6, caught three mutants between them and no others: M09 (finite-difference step mismatch), M11 (transposed Jacobians) and M18 (goods-market sign). This count excludes the two strict-xfail regression pins of T5 (F1, F2). Under the registered rule their XPASS also counts as a catch; with them, the union over all five IKC test functions is seven mutants: M05, M08, M09, M11, M18, M20, M25 (section 6).
  - The T6 matrix identity caught nothing.
  - Two household-block errors, Euler timing (M01) and the constraint at the wrong grid point (M16), and one fiscal steady-state error, steady-state taxes without debt service (M21), passed both T5 and T6 although their model executes the changed code.
  - Most T5/T6 catches came from checks other than the identities. T5's Walras test caught 9 mutants; 6 of them (M02, M05, M08, M13, M20, M25) no normally passing IKC test caught, and for 4 of those (M05, M08, M20, M25) the only other T5 signal was a strict-xfail pin turning XPASS. In T6, the test named `test_rank_benchmark_is_ricardian` caught 6 through two different assertions: the RA PV multiplier equal to 1, a closed form, caught M03, M05 and M20; the HA PV multiplier above 1 caught M08, M11 and M18. The nonlinear cohort simulation caught 5; it is a separate computation path on the same household code, policies and distribution, not an external reference.
- **T4:** the limit test itself (`r → 1/β − 1`, O(eps²)) and the monotonicity regressions caught nothing. All four T4 catches came from the v0.2 diagnostics: capital-market residual, Euler/KKT residuals, top-of-grid mass.
- **The baseline still passes:** the unmodified code gives 41 passed, 2 xfailed (75 s) on this machine, as in the README.

This is a measurement on a hand-picked sample of 25 errors. The rates below describe this sample only; they are not an estimate of how often the suite catches HANK implementation errors in general (section 9).

## 2. Setup in brief

Full details in the plan, sections 2–4.

- **Model side mutated, checker side untouched.** Mutations go into copies of `hank_limits/models.py` and `hank_limits/core.py`. `core.py` holds the household backward step, the borrowing constraint and the Young lottery, which the requested error types need; it is part of the implementation under test. Never changed: `checks.py`, `ra.py`, `metrics.py`, `diagnostics.py`, `suite.py`, `tests/`.
- **Detection.** A test catches a mutant if its outcome differs from the baseline: passed → failed or error, or xfailed → failed or XPASS(strict). A group catches a mutant if any of its tests does. Groups: T1–T6 = `tests/test_t1…t6`, D = `tests/test_diagnostics_v0_2.py`. This rule was pre-registered and was applied without change.
- **Reach.** "Reached" means the group's providers execute the changed lines. It was determined statically from the code map (plan, section 3) and confirmed dynamically: every group listed as not reached gave summaries identical to the baseline, apart from the one rounding-level number described below.
- **Runs.** 30 runs (2 baselines, 25 mutants, 3 controls), single-threaded, 4 in parallel, 300 s limit per fixture or test (never reached). About 75–150 s per run, 10.5 minutes in total.
- **Reproducibility.** The whole set was run twice: the first run, and a repeat with the runner as amended by A1, which only added storage of every run's check summaries. By the author's comparison, test outcomes and headline numbers are identical between the two; the first run's raw outputs were not kept, so this cannot be checked independently (section 11). Within the repeat, B0 and B1 are bit-identical except one number: T5's Walras residual in the deliberately mismatched-options solve, which differs by a relative 1.4e-11. That number lives on the ill-conditioned F2 path (LU element growth ~1e7), which amplifies rounding. The same 1.4e-11 shows up as "affected" for three mutants that do not reach T5 (M15, M19, C03); it is noise, not reach.

## 3. Mutant × test results

Cells: **n/N** = caught, with n of the group's N tests changing outcome; `·` = reached, not caught; `–` = not reached. `†` = the pre-registered prediction for this cell was wrong; `?` = no firm prediction was made.

| ID | Category | Change (short) | T1 | T2 | T3 | T4 | T5 | T6 | D | Suite |
|---|---|---|---|---|---|---|---|---|---|---|
| M01 | timing | Euler uses r_t, not r_{t+1} | **4/8** | **1/6** | **4/5** | · | · | · | · | caught |
| M02 | timing | debt service at r_{t-1} | – | – | – | – | **1/6** | · | – | caught |
| M03 | timing | tax rule on B_t, not B_{t-1} | – | – | – | – | · | **1/6** ? | – | caught |
| M04 | timing | T3 rate path shifted one period | – | – | **4/5** | – | – | – | – | caught |
| M05 | timing | gov. budget uses B_t, not B_{t-1} | – | – | – | – | **3/6** | **3/6** | – | caught |
| M06 | units | T2 rate path / 100 | – | **1/6** | – | – | – | – | – | caught |
| M07 | units | delta 0.02 (quarterly) in annual T4 | – | – | – | **1/7** ? | – | – | · | caught |
| M08 | units | income mean not normalised to 1 | · | **3/6** † | · | · | **2/6** † | **1/6** † | · | caught |
| M09 | jacobian | default FD step 1e-2 | **1/8** † | – | – | – | **4/6** | **2/6** | – | caught |
| M10 | jacobian | T2 saver J one-sided | – | **1/6** | – | – | – | – | – | caught |
| M11 | jacobian | Jacobians transposed | **3/8** | – | – | – | **4/6** | **3/6** | – | caught |
| M12 | jacobian | T3 cached Js, options dropped (F1 type) | – | – | · † | – | – | – | – | **missed** |
| M13 | distribution | Young lottery weights swapped | **1/8** | **1/6** | · ? | **1/7** | **1/6** † | **1/6** ? | **2/5** | caught |
| M14 | distribution | cohort distribution not renormalised | – | – | – | – | – | **1/6** | – | caught |
| M15 | distribution | T3 Markov matrix transposed | – | – | **2/5** † | – | – | – | – | caught |
| M16 | distribution | constraint at a_grid[1] | · | · | · | **1/7** | · | · | · | caught |
| M17 | distribution | T4 grid top 20, not 200 | – | – | – | **1/7** | – | – | **1/5** ? | caught |
| M18 | equilibrium | goods market Y - C + G | – | – | – | – | **4/6** | **2/6** | – | caught |
| M19 | equilibrium | T2 goods market Y = C_S | – | **2/6** | – | – | – | – | – | caught |
| M20 | equilibrium | gov. budget without interest | – | – | – | – | **2/6** | **2/6** | – | caught |
| M21 | equilibrium | Tax_ss omits r B | – | – | – | – | · | · | – | **missed** |
| M22 | equilibrium | T2 M_S from H block | – | **2/6** | – | – | – | – | – | caught |
| M23 | calibration | T2 net-supply target w/o /(1-lam) | – | · | – | – | – | – | – | **missed** |
| M24 | calibration | T4 borrowing limit -1 | – | – | – | · | – | – | · | **missed** |
| M25 | calibration | T5 beta target A = 4B | – | – | – | – | **2/6** | · | – | caught |
| C01 | control | T4 grid 1500 points | – | – | – | · | – | – | · | no alarm |
| C02 | control | Va = (1+r)/c^(1/eis) | · | · | · | · | · | · | · | no alarm |
| C03 | control | T3 target goods, not bonds | – | – | · | – | – | – | – | no alarm |

### Error sizes in the caught groups

The measured quantity of the failing test, its threshold, and the unmutated value.

| ID | Group: failing test | Measured | Threshold | Baseline |
|---|---|---|---|---|
| M01 | T1: dC/dr exact at eps = 0 | 1.00 | < 1e-7 | 2.9e-9 |
| | T2: GE IRF vs TANK (χ = 1), eps = 0.0125 | 0.25 (no decay with eps) | < 2e-3 | 1.0e-3 |
| | T3: proportional income vs RANK, b = 1e-4 | 0.249 | < 1e-3 | 4.3e-4 |
| M02 | T5: Walras residual, rate shock, T = 300 | 0.34 | < 1e-7 | 3.4e-10 |
| M03 | T6: \|RA PV multiplier − 1\|, T = 300 | 2.3e-10 | < 5·0.92^300 = 6.8e-11 | 2.5e-11 |
| M04 | T3: proportional income vs RANK, b = 1e-4 | 0.249 | < 1e-3 | 4.3e-4 |
| M05 | T5: Walras residual | 6e208 (explosive debt) | < 1e-7 | 2.0e-9 |
| | T6: RA PV multiplier, T = 150 | 4e103 | 1 ± 1.8e-5 | 1.000006 |
| M06 | T2: GE IRF vs TANK (χ = 1) | 0.99 at every eps | decreasing, < 2e-3 | 1.0e-3 |
| M07 | T4: mass in top 5 grid points | 1.7e-10 | < 1e-12 | 4.4e-29 |
| M08 | T2: H-type iMPC = I, eps ≤ 0.05 | 1.25e-3, 3.1e-4, 7.8e-5 (= eps²/2) | < 1e-10 | 3.9e-13 |
| | T5: Walras residual, T = 300 | 12.2 | < 1e-7 | 2.0e-9 |
| | T6: HA PV multiplier | −1.83 | > 1 | 3.38 |
| M09 | T1: J(T=50) = leading block of J(T=200) | NaN | < 1e-14 | 0 |
| | T5: IKC, rate shock (SSJ solve) | 8.4e-5 | < 1e-10 | ≤ 3.1e-13 |
| | T6: decomposition, rate shock | 7.0e-5 | < 1e-8 | 1.1e-11 |
| M10 | T2: GE goods-market residual | 5.9e-10 (1.2e-10 to 5.9e-10 over eps) | < 1e-10 | 6.9e-18 |
| M11 | T1: dC/dY exact at eps = 0 | 0.997 | < 1e-8 | 1.7e-9 |
| | T5: IKC, rate shock (SSJ solve) | 1.04 | < 1e-10 | ≤ 3.1e-13 |
| | T6: decomposition, rate shock | 11.7 | < 1e-8 | 1.1e-11 |
| M13 | T1: κ/eps² | ≈ 4e-4 at eps ≤ 0.05 | in (0.01, 0.1) | 0.027 |
| | T2: GE IRF vs TANK, eps = 0.0125 | 0.85, not monotone in eps | < 2e-3 | 1.0e-3 |
| | T4: capital-market residual | 1.1e-2 | < 1e-8 | 5.2e-10 |
| | T5: Walras residual, rate shocks | 2.0e-2 (G shock: 2.0e-9) | < 1e-7 | 3.4e-10 |
| | T6: cohort sum vs M x | 4.3e-4 | < 1e-4 | 3.2e-6 |
| | D: L1 distance to SSJ forward iteration | 0.42 (T4), 1.93 (T1), not shrinking with tolerance | must shrink 10× | 8e-9, 5e-6 |
| M14 | T6: cohort sum vs M x | 0.75 | < 1e-4 | 3.2e-6 |
| M15 | T3: countercyclical vs zero-liquidity solution, b = 1e-4 | 4.0e-2 | < 1e-3 | 6.0e-4 |
| M16 | T4: Euler error, mass-weighted mean | 1.3e-3 (max 8.2e-2) | < 1e-6 | 4.3e-8 |
| M17 | T4: mass in top 5 grid points | 2.6e-2 | < 1e-12 | 4.4e-29 |
| | D: L1 distance to SSJ forward iteration (T4) | 1.81, not shrinking | must shrink 10× | 8e-9 |
| M18 | T5: IKC, G shock | 632 | < 1e-10 | ≤ 1e-14 |
| | T6: decomposition, G shock, T = 150 | 1.0 | < 1e-8 | 9.8e-12 |
| M19 | T2: output amplification, χ = 1.5 relative to χ = 1 | 0.870 | 1.2727 ± 1e-3 | 1.2727 |
| M20 | T5: Walras residual, T = 300 | 2.9 | < 1e-7 | 2.0e-9 |
| | T6: RA PV multiplier, T = 150 | 4.08 | 1 ± 1.8e-5 | 1.000006 |
| M22 | T2: saver iMPC vs RA | 34 (M_S = I) | decreasing, < 2e-3 | 1.3e-3 |
| M25 | T5: Walras residual, rate shock, T = 300 | 51.9 | < 1e-7 | 3.4e-10 |

XPASS(strict) of the two xfail tests (F1, F2) happened under M05, M08, M20 and M25 (F2 test) and M05 (F1 test). These mutants change the ill-conditioned rule-based system. In each case another test in the same group also failed, so no group-level catch rests on an XPASS alone.

## 4. Detection rates

### By test group

| Group | Reference used | Caught / 25 | Caught / reached | Caught | Reached, missed |
|---|---|---|---|---|---|
| T1 | closed-form PIH Jacobians | 4 | 4 / 6 | M01, M09¹, M11, M13 | M08, M16 |
| T2 | closed-form TANK; GE residual | 7 | 7 / 9 | M01, M06, M08, M10, M13, M19, M22 | M16, M23 |
| T3 | closed-form RANK / zero-liquidity, built from the provider's reported Π, ζ | 3 | 3 / 7 | M01, M04, M15² | M08, M12, M13, M16 |
| T4 | 1/β − 1 (reported by the provider); v0.2 diagnostics | 4 | 4 / 7 | M07³, M13, M16, M17 | M01, M08, M24 |
| T5 | provider's own J and paths (IKC); Walras | 9 | 9 / 13 | M02, M05, M08, M09, M11, M13, M18, M20, M25 | M01, M03, M16, M21 |
| T6 | provider's own J (identities); RA closed form; nonlinear cohort simulation | 9 | 9 / 14 | M03⁴, M05, M08, M09, M11, M13, M14, M18, M20 | M01, M02, M16, M21, M25 |
| D | SSJ's own forward iteration | 2 | 2 / 7 | M13, M17 | M01, M07, M08, M16, M24 |
| **Suite** | | **21** | | | missed by all: M12, M21, M23, M24 |

¹ via NaN, not via the step size (section 7). ² by the size of a converging error, not by non-convergence (section 7). ³ incidental: grid adequacy, not the δ error. ⁴ the tail tolerance encodes the documented fiscal rule (a specification regression).

"Caught / reached" depends on what the mutants touch. T5 and T6 reach more mutants because six mutants target the fiscal model, which only they use. Compare the groups within the reach column, not across the "/ 25" column.

### By individual test (mutants only)

| Group | Test | Caught | Mutants |
|---|---|---|---|
| D | `test_inverse_iteration_report[audit_t4]` | 0 | - |
| D | `test_inverse_iteration_report[audit_t1_mid]` | 0 | - |
| D | `test_ssj_forward_iteration_converges_to_it[audit_t4]` | 2 | M13, M17 |
| D | `test_ssj_forward_iteration_converges_to_it[audit_t1_mid]` | 1 | M13 |
| D | `test_ssj_stopping_rule_is_not_an_accuracy_bound` | 0 | - |
| T1 | `test_exact_at_eps_zero` | 2 | M01, M11 |
| T1 | `test_errors_decrease_monotonically` | 1 | M01 |
| T1 | `test_second_order_convergence` | 2 | M01, M11 |
| T1 | `test_tightest_setting_accuracy` | 2 | M01, M11 |
| T1 | `test_impatience_vanishes_like_eps2` | 1 | M13 |
| T1 | `test_borrowing_constraint_slack` | 0 | - |
| T1 | `test_error_grows_with_observation_window` | 0 | - |
| T1 | `test_household_jacobian_independent_of_computing_horizon` | 1 | M09 |
| T2 | `test_htm_type_exact_below_threshold` | 1 | M08 |
| T2 | `test_saver_type_converges_to_RA` | 1 | M22 |
| T2 | `test_aggregate_impc_is_tank` | 1 | M22 |
| T2 | `test_ge_irf_converges_to_bilbiie` | 5 | M01, M06, M08, M13, M19 |
| T2 | `test_amplification_factor` | 2 | M08, M19 |
| T2 | `test_ge_solution_residual` | 1 | M10 |
| T3 | `test_acyclical_converges_to_rank` | 2 | M01, M04 |
| T3 | `test_countercyclical_converges_to_zero_liquidity_solution` | 3 | M01, M04, M15 |
| T3 | `test_countercyclical_deviates_from_rank` | 3 | M01, M04, M15 |
| T3 | `test_asymptotic_monotone` | 2 | M01, M04 |
| T3 | `test_steady_state_rate_converges` | 0 | - |
| T4 | `test_no_idiosyncratic_risk_limit` | 0 | - |
| T4 | `test_r_below_discount_rate_sanity_only` | 0 | - |
| T4 | `test_market_clearing_residual` | 1 | M13 |
| T4 | `test_euler_and_kkt_residuals` | 1 | M16 |
| T4 | `test_top_grid_mass_negligible` | 2 | M07, M17 |
| T4 | `test_monotone_in_risk_size_regression` | 0 | - |
| T4 | `test_monotone_in_persistence_regression` | 0 | - |
| T5 | `test_ikc_identity_refined_solve` | 3 | M09, M11, M18 |
| T5 | `test_ikc_identity_ssj_solver_fixed_rate` | 3 | M09, M11, M18 |
| T5 | `test_ikc_identity_ssj_solver_rule_short_horizon` | 3 | M09, M11, M18 |
| T5 | `test_ikc_identity_ssj_solver_rule_long_horizon` | 4 | M05, M08, M20, M25 |
| T5 | `test_ikc_identity_with_mismatched_options` | 1 | M05 |
| T5 | `test_walras_law` | 9 | M02, M05, M08, M09, M11, M13, M18, M20, M25 |
| T6 | `test_matrix_identity` | 0 | - |
| T6 | `test_decomposition_fiscal_shock` | 1 | M18 |
| T6 | `test_decomposition_rate_shock` | 2 | M09, M11 |
| T6 | `test_rank_benchmark_is_ricardian` | 6 | M03, M05, M08, M11, M18, M20 |
| T6 | `test_wealth_group_split_adds_up` | 5 | M05, M09, M11, M13, M14 |
| T6 | `test_first_round_split_independent_of_T` | 2 | M05, M20 |

11 of the 43 tests caught nothing in this sample: the two inverse-iteration reports and the stopping-rule test (D); constraint slack and observation window (T1); the steady-state rate (T3); the no-risk limit, `r < 1/β − 1` and both monotonicity regressions (T4); the matrix identity (T6). This says nothing against those tests' original purpose (most are not designed to detect implementation errors). It says what this sample of errors did not touch.

### By category

| Category | Mutants | Caught by the suite | Missed | Catches by group (T1 T2 T3 T4 T5 T6 D) |
|---|---|---|---|---|
| timing | 5 | 5 | – | 1 1 2 0 2 2 0 |
| units | 3 | 3 (2 without M07) | (M07 incidental) | 0 2 0 1 1 1 0 |
| Jacobian | 4 | 3 | M12 | 2 1 0 0 2 2 0 |
| distribution | 5 | 5 | – | 1 1 1 3 1 2 2 |
| equilibrium | 5 | 4 | M21 | 0 2 0 0 2 2 0 |
| calibration | 3 | 1 | M23, M24 | 0 0 0 0 1 0 0 |
| **all** | **25** | **21** | 4 | 4 7 3 4 9 9 2 |

Categories have 3–5 members each; a single mutant moves a category rate by 20–33 points.

## 5. Missed errors and why

**Missed by every test**

- **M12 — T3 solve with cached Jacobians but without matching options (F1 type).** SSJ re-differentiates only the propagation of the exogenous rate path (one-sided, h = 1e-4). The output response comes from H_U, which is built from the cached Jacobians. A rate step of 1e-4 is small relative to everything else, so dY moves by 5e-5 to 7e-5 (relative) at every b. The T3 errors at b = 1e-4 are 5.0e-4 / 6.7e-4 (baseline 4.3e-4 / 6.0e-4) and the fitted rates 0.93 / 0.95 (baseline 0.99), all within tolerance. The inconsistency is visible in the GE target: in a separate diagnostic run (not part of the registered measurement), the asset-market residual of the T3 solution is about 2e-7 (unmutated: ~1e-18). T3 has no test on this residual, and `check_werning` does not record it. T2 has one (threshold 1e-10) and caught the analogous M10 at 5.9e-10.
- **M21 — steady-state taxes omit debt service.** The steady-state government budget is violated (fiscal_res_ss = −rB). No test checks steady-state residuals. The linearised dynamics, the A = B calibration and Walras in deviations are all unaffected. The HA PV multiplier changes from 3.383 to 3.448, and nothing flags it.
- **M23 — TANK net-supply target without /(1 − λ).** As predicted, the error vanishes in the limit the test looks at: H types end at the zero limit, so A_H → 0. The T2 IRF error at eps = 0.2 jumps from 0.13 to 1.24, but it is identical to the baseline from eps = 0.05 down. The tests check monotone decay and the smallest-eps value, both of which pass.
- **M24 — Aiyagari borrowing limit −1 instead of 0.** The no-risk limit and the monotonicity regressions hold for any fixed borrowing limit. Gap at eps = 0.025: 5.4e-5 (baseline 6.2e-5), rate 1.94. The Euler/KKT diagnostics use the mutated grid consistently.

**Caught, but missed by groups that execute the code**

- **M01 (Euler timing) in T4, T5, T6, D.** The steady state is unchanged, so T4 and D see nothing. T5 and T6 build their identities from the same (wrong) household Jacobians that drive the GE paths, so the identities hold (T5 IKC 9.8e-12, T6 decomposition ≤ 3.4e-10, as in the baseline). Only the external closed forms in T1–T3 caught it.
- **M16 (constraint at a_grid[1]) in T1, T2, T3, T5, T6.** Only T4's Euler diagnostic flagged it: the states held at a_grid[1] are classified as unconstrained and show Euler gaps. In the fiscal model of T5/T6 (19.8% constrained) the HA PV multiplier moves from 3.383 to 3.391, and the T6 cohort check moves from 3.2e-6 to 3.3e-5, below its 1e-4 tolerance.
- **M08 (income mean not normalised) in T1, T3, T4.** T1 cannot see it: the error in dC/dY is eps²/2, the same order as the PIH limit itself (T1 errors 8.7e-4 vs 9.5e-4 at eps = 0.0125). T3 and T4 references are scale-free.
- **M13 (Young lottery swapped) in T3.** The T3 summaries change by about 1% only.
- **M02 and M25 in T6, M03 in T5.** These leave a model that is internally consistent for the paths T6 and T5 test (see section 6).
- **M07, M24 in D; M08, M16 in D.** The distribution audit compares two solvers for the same policy and is not designed to catch policy or calibration errors.

## 6. T5 and T6: what an internal-consistency test can and cannot see

The external review's point was that T5 and T6 compare quantities built from the same Jacobians, so an error that enters both sides passes. The sample confirms this for the identity parts, and shows that most T5/T6 catches came from the other checks in those files.

1. **The identity tests caught only inconsistencies between the cached Jacobians and the solver, and GE-structure errors.** The three normally passing IKC tests (T5) and the two decomposition tests (T6) caught exactly M09 (Jacobians at h = 1e-2, solver at 1e-4), M11 (transposed Jacobians handed to the solver) and M18 (goods-market sign). They caught none of the household-block errors that reach their model (M01, M08, M13, M16), none of the fiscal-block errors (M02, M03, M05, M20, M21), and not the calibration error M25.

   Scope: this excludes the two strict-xfail regression pins, `test_ikc_identity_ssj_solver_rule_long_horizon` (F2) and `test_ikc_identity_with_mismatched_options` (F1). The registered rule counts their XPASS as a catch, and the per-test table in section 4 does so: the F2 pin turned XPASS under M05, M08, M20 and M25, the F1 pin under M05. Over all five IKC test functions the union is therefore seven mutants: M05, M08, M09, M11, M18, M20, M25. An XPASS means that the known F1/F2 failure did not appear on the mutated system (section 3); it is not an identity violation caused by the error.
2. **The T6 matrix identity caught nothing.** (I − M_HA)^{-1} − (I − M_RA)^{-1} = (I − M_HA)^{-1}(M_HA − M_RA)(I − M_RA)^{-1} holds for any invertible matrices. It checks linear algebra, not the model.
3. **T5's Walras test is the most productive single test in the suite: 9 catches** (M02, M05, M08, M09, M11, M13, M18, M20, M25). Six of these no normally passing IKC test caught (excluding the two strict-xfail pins). M02 and M13 were caught by nothing else in T5; M05, M08, M20 and M25 otherwise only turned an xfail test into XPASS. It is not internal in the same sense: it compares the household budget (inside the HA block) with the government budget (a separately coded simple block). A timing, unit or calibration error on one side shows up as a residual (A_ss − B_ss) dr, B (dr_t − dr_{t−1}), r dB, … under rate shocks.
4. **T6's catches came from the checks other than the decomposition identities: one closed-form benchmark, one condition on the HA solution, and separate computation paths.**
   - `test_rank_benchmark_is_ricardian` (6 catches) is a composite test with two assertions, and they caught different mutants:

     | Assertion that failed | Mutants | Failing value (T) |
     |---|---|---|
     | RA PV multiplier = 1 within the tail bound (closed form) | M03, M05, M20 | 1.00000000023 (300), 4.4e103 (150), 4.08 (150) |
     | HA PV multiplier > 1 | M08, M11, M18 | −1.83, −33.2, −393 (150) |

     In the stored summaries of M03, M05 and M20 the HA multiplier is above 1 (3.46, 4.9e103, 35.3 at T = 150), so only the RA assertion fails there. In those of M08, M11 and M18 the RA multiplier equals the baseline (1.0000064 at T = 150), so only the HA assertion fails. The closed-form RA benchmark therefore caught three mutants; the other three were caught by a sign-and-size condition computed from the same HA solution as the rest of T6.
   - The nonlinear cohort simulation (M05, M09, M11, M13, M14) does not use the fake-news Jacobian, so it is a separate computation path. It runs on the same household code, policies and steady-state distribution, so it is not an external reference in the sense of T1–T3.
   - The T-independence of the first round (M05, M20) is a numerical-consistency condition on the same model.
5. **Errors that pass every T5/T6 test, although the code runs in their model:** M01 (Euler timing), M16 (borrowing constraint), M21 (steady-state taxes); and M03 for T5. These are the cases the review anticipated.

## 7. Pre-registered predictions versus results

**Firm predictions that failed (7 of 170 cells)**

| ID, group | Predicted | Observed | Why the prediction was wrong |
|---|---|---|---|
| M08, T2 | miss | caught | I reasoned that references use the reported steady-state levels and are scale-free. But Jacobians with respect to **aggregate** Y carry the factor E[e]: the hand-to-mouth iMPC becomes E[e]·I = (1 + eps²/2)·I, and T2's exactness test (1e-10) sees eps²/2 exactly. |
| M08, T5 | miss | caught | In the fiscal model sd = 0.6, so E[e] = exp(0.18) = 1.197. Household income no longer adds up to Y − Tax, so the household and government budgets disagree: Walras residual 12 at T = 300. |
| M08, T6 | miss | caught | Same cause: the effective MPC out of Y is 1.197 × the true one, and the HA PV multiplier turns negative (−1.83), failing "HA multiplier > 1". |
| M09, T1 | miss | caught | I argued that T1 passes h explicitly and that the horizon test compares two Jacobians with the same step. What I missed: a rate step of 1e-2 pushes households near the 0.99 × natural borrowing limit into negative consumption, the Jacobian contains NaN, and the horizon test fails on NaN. The test catches the NaN, not the step size. |
| M12, T3 | catch | missed | I expected the one-sided h = 1e-4 to move income by about b at b = 1e-4. But only the rate-path propagation is re-differentiated; Y's effect comes from H_U built from the cached Jacobians (section 5). |
| M13, T5 | miss | caught | I expected β recalibration to keep T5 consistent. But the swapped lottery does not preserve mean assets between end of period (A, calibrated to B) and the next beginning of period. Under rate shocks, households then earn interest on a different stock than the government pays on: Walras 2.0e-2 under rate shocks, 2.0e-9 (baseline) under the G shock. (Mechanism consistent with the pattern; not separately verified.) |
| M15, T3 | miss (L) | caught | The self-consistency argument was right: the reference uses the provider-reported Π, the proportional-income case passes, and the countercyclical error falls with b (0.99, 0.99, 0.94, 0.67, 0.33, 0.11, 0.040; fitted rate 0.91). But Π^T raises the compounding factor δ from 1.017 to 1.277 (zero-liquidity response 250 times RANK's), so the error at b = 1e-4 is 0.040 against a tolerance of 1e-3 set for the baseline calibration. The catch comes from the size of the constant, not from non-convergence; at small enough b this mutant would pass. |

**Cells without a firm prediction (`?`)**

- M03/T6: caught at T = 300 by the mechanism stated in the plan. The Ricardian tail bound 5·0.92^T encodes the documented rule's decay rate; the mutated rule decays at 0.927. Note that the baseline itself sits at 2.5e-11 against a bound of 6.8e-11 at T = 300, a margin of 2.7×.
- M07/T4: caught by the top-of-grid mass, as anticipated (incidental).
- M13/T3: missed.
- M13/T6: caught by the cohort check.
- M17/D: caught.

**Mechanisms right, magnitudes off**

- M10: predicted a GE residual of about 1e-6; observed 1.2e-10 to 5.9e-10, only 1.2–6 times the 1e-10 threshold. A slightly looser threshold would have missed it.
- M13 in T1: the wrong distribution made T1's Jacobian errors **smaller** (3.0e-4 vs 9.5e-4 at eps = 0.0125; fitted rates 2.29 / 2.19 vs 1.95 / 2.03). The convergence tests passed; only the κ/eps² band flagged it.
- M13 in T4: the limit test passed (gap rate 3.16 > 1.8); only the capital-market residual flagged it.

## 8. Controls

No control changed any test outcome.

- **C01 — T4 grid 1000 → 1500 points.** Gap at eps = 0.025: 6.05e-5 (baseline 6.17e-5); rate 1.937 (1.931).
- **C02 — Va = (1+r)/c^(1/eis).** Changes at rounding level only; for example the T3 errors change by 2.6e-6 relative.
- **C03 — T3 GE target goods market instead of bond market.** This is Walras-equivalent in exact arithmetic. At the b where the tests judge, dY changes by 1.2e-6 (proportional income) and 1.3e-4 (countercyclical, relative), with no outcome change.

  It is **not** numerically equivalent at large b with countercyclical incidence: dY differs by 35% at b = 0.1 and 11% at b = 1e-2. A separate diagnostic run found two contributing facts:
  - In the countercyclical model the linearised aggregate budget identity dC + dA = (1+r) dA_{−1} + dY holds in the household Jacobians only to 9.3e-6 at b = 0.1 and 9.5e-9 at b = 1e-4 (proportional income: 3.7e-12 and 5.4e-10). That is roughly proportional to the difference step h = 1e-3 b.
  - The goods-market system is badly conditioned: cond(J^{C,Y} − I) = 4.9e4.

  **Cause: not known at the time of the measurement.** The report then suggested that a kink, such as the borrowing constraint, is crossed. The follow-up below identifies the cause: SSJ differentiates the income hetinput one-sided, although it asks for two-sided differences. No kink is needed to explain the residual. Neither the residual nor its cause changes a test outcome.

### Follow-up diagnosis, 2026-10-22, outside the registered measurement

This was run after the registered measurement and after the review, on the unmodified branch. It does not change the registered patches, tolerances or outcomes, and it is not one of the 30 registered runs. The C03 result above stands: no outcome change.

Scripts: `scripts/c03_probe_v0_1.py` is the reviewer's probe (received file SHA-256 `c7926f8883bb2762f71e6dbacbb03de7f90c1a7a3ba9121910f6bb90be4956a9`); the only change is its output path, now `out/c03_probe_v0_1.json`. `scripts/repro_F3_v0_1.py` uses SSJ code only (received file, unchanged; SHA-256 `f220d369088c94c64a7701bbc617e2a87c9db958b353ea2945112e7f7bb9ebe2`). Environment as in section 9, one thread; the probe takes about 5 s.

**1. Cause: the hetinput derivative is one-sided.** In sequence-jacobian 1.0.0 (the PyPI wheel; SHA-256 `40715634…ef01` checked again here):

- `blocks/het_block.py`, `HetBlock.jac_backward_prelim`, lines 400–401: the comment "always use two-sided differentiation for hetinputs", then `self.hetinputs.differentiable(ss, h, True)`. This stores `default_twosided = True` (`utilities/function.py`, line 231).
- `utilities/function.py`, line 233: `DifferentiableCombinedExtendedFunction.diff(self, shock_dict, h=None, outputs=None, hide_zeros=False, twosided=False)`. The default is `False`, not `None`, so the branch `if twosided is None: twosided = self.default_twosided` (lines 234–235) is never taken when the caller omits the argument. The single-function class `DifferentiableExtendedFunction.diff` (line 97) has the default `None`, which that branch expects.
- `blocks/het_block.py`, `HetBlock.backward_fakenews`, line 302: `differentiable_hetinput.diff({input_shocked: 1})`, without `twosided`. The contemporaneous income response that enters the fake-news algorithm is therefore a forward difference.
- `blocks/stage_block.py` follows the same pattern: lines 422–423 (`differentiable(ss, h, True)`, same comment) and line 287 (`.diff(din_dict)`). The suite does not use `StageBlock`.
- GitHub master `eb4ac4c` has the same code according to the research session's check; this was not re-checked here.

This concerns the hetinput only. The household backward function is differentiated two-sided as requested (`het_block.py`, line 410). In this suite only T3's `income_cyclical` with ζ ≠ 0 is a nonlinear hetinput. `income_basic`, `income_wage`, `income_aftertax` and `income_cyclical` with ζ = 0 are linear in their inputs, so forward and central differences agree for them up to rounding.

`scripts/repro_F3_v0_1.py` builds the T3 income function with three income states and calls SSJ as `HetBlock` does. Output here:

```
h=0.01: |used-forward| 0.0e+00  |used-central| 2.4e-03  |used-analytic| 2.4e-03  |twosided=None - analytic| 8.4e-06
h=0.001: |used-forward| 0.0e+00  |used-central| 2.4e-04  |used-analytic| 2.4e-04  |twosided=None - analytic| 8.4e-08
h=0.0001: |used-forward| 0.0e+00  |used-central| 2.4e-05  |used-analytic| 2.4e-05  |twosided=None - analytic| 8.4e-10
```

The derivative SSJ uses is exactly the forward difference, with an O(h) error. Passing `twosided=None` applies the stored default and gives an O(h²) error.

**2. Size of the term.** With Yss = 1 the countercyclical income is y_i(Y) = e_i Y^(1+ζ_i), and the calibration gives Σ_i π_i e_i (1+ζ_i) = 1 (aggregate income equals Y). The forward difference of aggregate income is

    Σ_i π_i [y_i(1+h) − y_i(1)] / h = 1 + (h/2) Σ_i π_i e_i ζ_i (1+ζ_i) + O(h²) = 1 + 0.0927104 h + O(h²).

For proportional income (ζ = 0) the h term is zero. T3 uses h = 1e-3 b, so the term is 9.27e-5 b: 9.27e-6 at b = 0.1 and 9.27e-9 at b = 1e-4. These match the residuals reported above (9.3e-6 at b = 0.1; 9.5e-9 at b = 1e-4). At b = 1e-4 the term (9.27e-9) is of the same order as the rounding error of a forward difference with step 1e-7 (a few times 1e-9), and the residual no longer halves with h (table below).

The residual is R_Y = J^{C,Y} + J^{A,Y} − (1+r) L J^{A,Y} − I, with L the one-period lag matrix. At b = 0.1 and h = 1e-4 its diagonal is 9.27128e-6 and its largest off-diagonal element 4.9e-14. The aggregate income derivative that SSJ uses, measured directly, is 1 + 9.27129e-6; the formula gives 1 + 9.27104e-6, and the difference of 2.5e-10 is the O(h²) term. Subtracting (Σ_i π_i dy_i/dY − 1) I from R_Y leaves 1.81e-11.

**3. Results.** max |R_Y| and the dY gap = max |dY_goods − dY_assets| / max |dY_assets|, for the same steady state, T = 150 and dr_t = −0.01 · 0.8^t. Hetinput derivative: "original" = SSJ as installed; "two-sided" = the stored default applied; "analytic" = e_i(1+ζ_i). The two controls replace `DifferentiableCombinedExtendedFunction.diff` inside the probe process only. Review values: Windows, Python 3.13.3, the same package versions (review, section 6); "here": the environment of section 9. Rows without a review value were not in the review table. Countercyclical income unless stated otherwise.

| b | Hetinput derivative | h | max \|R_Y\|, review | max \|R_Y\|, here | dY gap, review | dY gap, here | cond(J^{C,Y} − I), here |
|---|---|---|---|---|---|---|---|
| 0.1, proportional | original | 1e-4 | – | 3.71e-12 | – | 0.000152% | 3.60e5 |
| 0.1 | original | 1e-4 | 9.27128e-6 | 9.27128e-6 | 35.1944% | 35.1944% | 4.86e4 |
| 0.1 | original | 5e-5 | 4.63557e-6 | 4.63557e-6 | 35.4089% | 35.4089% | 9.77e4 |
| 0.1 | original | 2.5e-5 | 2.31775e-6 | 2.31775e-6 | 35.8377% | 35.8377% | 1.98e5 |
| 0.1 | two-sided | 1e-4 | 2.33624e-10 | 2.32476e-10 | 0.145606% | 0.144877% | 8.21e6 |
| 0.1 | two-sided | 5e-5 | 4.49112e-11 | 4.48541e-11 | 0.0243933% | 0.0243017% | 8.19e6 |
| 0.1 | analytic | 1e-4 | 1.78646e-11 | 1.78848e-11 | 0.0157501% | 0.0157473% | 8.17e6 |
| 0.01 | original | 1e-5 | 9.27087e-7 | 9.27087e-7 | 10.6714% | 10.6715% | 1.82e5 |
| 0.01 | original | 5e-6 | – | 4.63542e-7 | – | 4.25136% | 1.45e5 |
| 0.01 | original | 2.5e-6 | – | 2.31687e-7 | – | 1.92890% | 1.31e5 |
| 0.01 | two-sided | 1e-5 | 2.33202e-11 | 2.32121e-11 | 0.000175713% | 0.000175170% | 1.20e5 |
| 0.01 | two-sided | 5e-6 | – | 2.78932e-12 | – | 0.0000172919% | 1.20e5 |
| 0.01 | analytic | 1e-5 | – | 1.76599e-11 | – | 0.000134687% | 1.20e5 |
| 1e-4 | original | 1e-7 | 9.41889e-9 | 9.54517e-9 | 0.0126776% | 0.0127509% | 3.05e4 |
| 1e-4 | original | 5e-8 | – | 3.86951e-9 | – | 0.00558774% | 3.05e4 |
| 1e-4 | original | 2.5e-8 | – | 6.97757e-9 | – | 0.00946485% | 3.05e4 |
| 1e-4 | two-sided | 1e-7 | – | 1.47269e-10 | – | 0.0000164992% | 3.05e4 |
| 1e-4 | two-sided | 5e-8 | – | 1.09690e-10 | – | 0.000402365% | 3.05e4 |
| 1e-4 | analytic | 1e-7 | – | 4.10195e-10 | – | 0.000733917% | 3.05e4 |

All 24 rows, including the other proportional-income rows (max |R_Y| 1.8e-12 to 6.8e-12 in every mode), are in `out/c03_probe_v0_1.json`. The original rows at b = 0.1 and 0.01 agree with the review to the printed digits, apart from the last digit of one dY gap. The two-sided and analytic rows agree to two to four digits; these are residuals of 1e-10 to 1e-11 in a system with condition number about 8e6. At b = 1e-4 the original row differs by about 1%, where the difference is dominated by rounding.

**4. What remains.**

- Two-sided differences alone do not make the two targets numerically equivalent at b = 0.1. The dY gap is 0.145% at h = 1e-4 and 0.024% at h = 5e-5, and 0.016% with the analytic derivative, where max |R_Y| = 1.8e-11 remains. cond(J^{C,Y} − I) rises from 4.9e4 to 8.2e6, while cond(J^{A,Y}) stays at about 174.
- In the countercyclical rows with b ≥ 0.01, the dY gap is 0.5 to 1.1 times cond(J^{C,Y} − I) · max |R_Y|. So the gap has the size of the budget-identity residual amplified by the goods-market condition number. With the original derivative at b = 0.1, cond(J^{C,Y} − I) doubles each time h is halved, so the product and the gap stay at 35–36% although R_Y halves. This describes the 24 rows; it is not a separate proof.
- A small residual in a solve's own target does not show that the two solves agree. With the original derivative at b = 0.1, the asset-market solution leaves a goods-market residual of 2.5e-7, and the goods-market solution leaves an asset-market residual of 1.1e-4.
- The sources of the 1.8e-11 floor and of the remaining 0.016% gap (other differentiation residuals, steady-state residuals, the sensitivity of the goods-market solve) have not been separated.

**5. Status of the registered result.** Unchanged. In the registered measurement C03 changed no test outcome; at b = 1e-4, where the T3 tests judge, the probe gives a dY gap of 1.3e-4 (relative), as reported above. The follow-up explains the O(h) residual at large b. It does not show that C03 is numerically equivalent at every b.

## 9. Limitations

- **A sample, not the error space.** The 25 mutants were chosen by hand, by the same author who wrote the predictions and who knew the suite. They are not drawn from any distribution of real errors. Detection rates are descriptive of this sample and must not be read as a general guarantee, for example "the suite catches 84% of HANK implementation errors".
- **One error at a time.** Interactions are not studied; two errors can offset each other.
- **"Caught" means an outcome changed.** It does not mean the failing test points to the cause: M07 was flagged as a grid problem, M09 in T1 as a NaN.
- **Catches depend on threshold margins.** The tolerances are regression thresholds for this calibration, and several catches have small margins: M10 at 1.2–6× its threshold, M03 at 3.3× (and the baseline at 2.7× below the same bound), M15 by the size of a converging error.
- **One environment.** Python 3.13.16, sequence-jacobian 1.0.0, numpy 2.5.3, scipy 1.18.1, numba 0.68.0 (llvmlite 0.50.0), scipy-openblas 0.3.34, pytest 9.1.1, Linux 6.18 x86_64, 4 CPUs, one thread per run. One rounding-level number (1.4e-11 relative) differed between runs.
- **Three controls only.** They show that the tests do not alarm on three equivalent changes. That is not a false-alarm rate.
- **Scope of "model side".** `core.py` was mutated as part of the model side (section 2). A reader who counts it as infrastructure should discount M01, M08, M13, M16 and C02.
- **Not known.** Whether M12-type mismatches are caught at other option settings; whether the M13 mechanism in T5 is exactly as described in section 7.
- **Identified after the measurement (follow-up diagnosis, section 8).** The cause of the O(h) budget-identity residual in the countercyclical T3 Jacobians: SSJ 1.0.0 differentiates the income hetinput one-sided. Still open: the 0.15% C03 dY gap at b = 0.1 that remains with the two-sided hetinput derivative.

## 10. What the misses suggest (not implemented here)

The checker is frozen on this branch, so these are notes for the author, not changes.

1. **T3: assert the GE target residual,** as T2 does. It would have caught M12 (2e-7 vs about 1e-18) and is cheap; `check_werning` already receives it in `diag['ge_target_resid']`.
2. **Steady-state residuals:** government budget, goods market and asset market at the steady state. These would catch M21 directly, and M08 more directly than now.
3. **An explicit income-normalisation check** (E[e] = 1, aggregate income = Y).
4. **T4:** the analytic limit caught none of the seven mutants that reach it, while the diagnostics caught four. The README could say that T4's error-detection value lies in the diagnostics.
5. **T6:** describe the matrix identity as a linear-algebra check. The parts that detected errors are the test named `test_rank_benchmark_is_ricardian`, whose RA assertion (a closed form) caught three mutants and whose HA-multiplier assertion caught three, and the cohort simulation and first-round T-independence, which are separate computation paths on the same household code, policies and distribution rather than external references. Splitting the Ricardian test into its two assertions would make this attribution visible in the test outcomes.
6. **Specification pins.** M03 (fiscal rule), M07 (δ) and M24 (borrowing limit) are deviations from the documented model that no limit test can see by construction. Only direct checks of the specification, such as asserting calibration constants or policy-rule coefficients, would catch them.

## 11. Reproduction

```bash
pip install sequence-jacobian==1.0.0 pytest          # numpy, scipy, numba as in section 9
python3 -m pytest                                    # unmodified suite: 41 passed, 2 xfailed
python3 scripts/mutation_v0_1.py --check             # apply and compile every patch, no tests
python3 scripts/mutation_v0_1.py --jobs 4            # baselines + 28 runs -> out/mutation_v0_1.{json,csv}
python3 scripts/c03_probe_v0_1.py .                  # follow-up C03 diagnosis (section 8) -> out/c03_probe_v0_1.json
python3 scripts/repro_F3_v0_1.py                     # one-sided hetinput derivative, SSJ code only
```

**Execution record.**

- Order (2026-10-21, Korea time): pre-registration `ccb7e38` 19:32:35; README correction `e51d7b8` 19:32:55; start of the repeat run 19:45:13 (`environment.utc` in the JSON, 10:45:13 UTC); amendment A1 `d311711` 19:48:26; results `65f8f6d` 20:00:44.
- The published outputs come from the repeat run. It started before A1 was committed, so it ran the A1 change as an uncommitted working copy of the runner. The JSON therefore records HEAD `e51d7b8`, and that record is left as it is. The A1 change only affects output storage.
- The runner as committed in `d311711` (unchanged in `65f8f6d`) is git blob `c78fcb88abc00b9663b52f2fc0d0e8713e6b89c6`. That the working copy used for the repeat had exactly this content is the author's statement. The runner's docstring was corrected later (amendment A2 in the plan; text only), so the current file is a different blob.
- `model_or_tests_dirty = false` in the JSON comes from `git status --porcelain hank_limits tests`. It covers `hank_limits/` and `tests/` only and says nothing about whether the runner itself was committed.
- The first full run's raw outputs, push times and run logs were not kept. Two statements are therefore the author's and have not been verified independently: that the pre-registration was pushed before the first run, and that the first run and the repeat gave identical test outcomes and headline numbers (section 2). What git shows is the commit order above.
- In the JSON, `settings.workdir` is the temporary directory used on the run machine. The `scratch` diagnostics cited in sections 5 and 8 (T3 provider calls under M12, C03 and M15, and the budget identity in the T3 Jacobians) were run outside the registered measurement and are not part of the outputs. The follow-up C03 diagnosis of section 8 is stored separately in `out/c03_probe_v0_1.json`.
