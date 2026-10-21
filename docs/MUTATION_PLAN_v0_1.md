# Mutation testing plan v0.1 (pre-registration)

Status: **pre-registered**. This file and `scripts/mutation_v0_1.py` (mutant list, predictions, runner) were committed and pushed to branch `mutation-v0_1` before any mutant was run. Results are in `docs/MUTATION_v0_1.md` and `out/mutation_v0_1.{json,csv}`. Changes made after results were seen are listed under "Amendments" with their reason and time.

## 1. Purpose

An external review noted that the evidence differs across tests: T1–T4 compare with an external reference (analytic solutions), while T5 and T6 are internal-consistency checks built from the same Jacobians, so an error that enters both sides can pass. The review's first recommendation was to inject timing, unit and similar errors on purpose and measure what the suite catches and what it misses. This is that measurement.

What it can show: for a sample of small, plausible implementation errors, which tests react. What it cannot show: a detection rate for HANK implementation errors in general. The mutant list is a hand-picked sample, not the error space.

## 2. Rules

1. Work only on branch `mutation-v0_1`; no merge into `main`, no tags, no GitHub Release (the Zenodo integration would mint a DOI).
2. Tolerances, check code and references are not changed. Mutations are applied to **copies** of the model side only:
   - `hank_limits/models.py`: providers and their model blocks;
   - `hank_limits/core.py`: the household block (backward step) and the stationary distribution (`FastHetBlock`) that every provider builds on.

   `core.py` is included because the timing, budget-constraint, borrowing-limit and Young-lottery code the review asks about lives there, not in `models.py`, and it is part of the implementation under test, not of the checker. The checker side, never changed: `checks.py`, `ra.py` (closed-form references), `metrics.py`, `diagnostics.py` (Euler/KKT residuals and distribution audit used by tests), `suite.py` (canonical settings), `tests/`.
3. Pinned dependencies: sequence-jacobian 1.0.0 (wheel SHA-256 `40715634…ef01`, as in the README), environment recorded in the result JSON.
4. Pre-registration: this file and the script docstring/data. Before the commit, two things were run: `--check` (each patch applied to a copy and byte-compiled; no tests) and a smoke test of the harness on the **unmutated** copy only (two baseline runs: 41 passed, 2 xfailed each, bit-identical summaries). No mutant was run before the commit.

## 3. Which test uses which provider code

| Group | Tests | Entry (suite.py) | Provider and model code executed (`models.py`, `core.py`) | Reference the result is compared with |
|---|---|---|---|---|
| T1 | `test_t1_rank_limit.py` (8) | `run_t1` → `check_rank_limit` | `ssj_rank_limit_provider` → `pih_household` (`rouwenhorst`, `make_household` + `income_basic`, `asset_grid` at eps = 0, β calibration by `brentq`, `_Dbeg_override`, `FastHetBlock` stationary distribution), `household_jacobians` (h passed explicitly). The horizon test calls `pih_household` and `household_jacobians` (default h) directly. | `ra.ra_jacobians` at the provider-reported (r, eis, C, A): external closed form. |
| T2 | `test_t2_tank_limit.py` (6) | `run_t2` → `check_tank_limit` | `ssj_tank_provider` → `tank_model` (`pih_household` for savers, H household, `tank_incidence`, `_remap_household`, `tank_aggregate`, `goods_clearing`), two-sided block Jacobians, `solve_impulse_linear` with `_het_opts`. | `ra_jacobians`, `tank_impc`, `tank_output_fixed_real_rate` (closed forms at reported r, eis, C_S, A_S, λ); GE residual of the solver. |
| T3 | `test_t3_werning.py` (5) | `run_t3` → `check_werning` | `ssj_werning_provider` → `werning_primitives` (`rouwenhorst`, β from r*, ζ), `werning_model` (`income_cyclical`, `huggett_clearing`, `brentq` on r), Jacobian h = 1e-3 b, `solve_impulse_linear` (target `asset_mkt`) with `_het_opts`. | `ra_output_fixed_real_rate`; `zero_liquidity_output` built from the **provider-reported** e_grid, Π, ζ, r_ss, eis. |
| T4 | `test_t4_aiyagari.py` (7) | `run_t4` → `check_aiyagari` | `ssj_aiyagari_r` (`AIYAGARI_DEFAULTS`, `income_wage`, firm FOCs, `brentq` on the capital market, `diagnostics.euler_kkt_residuals`). | 1/β − 1 as reported by the provider; Euler/KKT residuals; top-of-grid mass; monotonicity regressions. |
| T5 | `test_t5_ikc.py` (4 + 2 xfail) | `run_t5` → `check_ikc` | `ssj_fiscal_provider` → `fiscal_hank` (`income_aftertax`, `fiscal`, `fiscal_clearing`, `real_rate_rule`, β calibration), `household_jacobians` (default h), `solve_impulse_linear` or `solve_impulse_linear_refined`, `_het_opts`. | **Internal**: IKC built from the provider's own J and its own dG, dTax, dr paths; Walras residual of the asset market. |
| T6 | `test_t6_differencing.py` (6) | `run_t6` → `check_differencing` | `ssj_fiscal_provider` (refined solve), `cohort_impc_paths`, `suite._wealth_quartile_masks`. | **Internal** identities (resolvent, decomposition) on the provider's own J and paths; RA closed form at the reported (r, eis, C, A) for the Ricardian check; cohort sum vs M x. |
| D | `test_diagnostics_v0_2.py` (5) | fixtures in the test file | `ssj_aiyagari_r(0.2, 0.9)`, `pih_household(0.1)`, `core.stationary_distribution`. | SSJ's own forward iteration on the same policy (`diagnostics.distribution_audit`). |

## 4. Method

- **Run.** For each mutant: copy `hank_limits/`, `tests/`, `pytest.ini` to a temporary directory, apply the patch, run the full pytest suite there with the script as a plugin. Single-threaded BLAS per run; up to 4 runs in parallel. A fixture or test that runs longer than 300 s is stopped and counted as an error.
- **Baseline.** The unmutated copy, run twice (B0, B1) in the same harness. B0's outcome per test is the reference.
- **Detection.** A test catches a mutant if its outcome differs from B0 (passed → failed or error; xfailed → failed or xpassed). A group catches a mutant if any of its tests does. The mode is recorded: `assert` (an assertion failed), `exception`, `error` (fixture raised), `timeout`.
- **Error size.** Each check's result is summarised (errors along the limit sequence, fitted rates, residuals) and stored for every run; two headline numbers per group are tabulated against the baseline. `affected` = the largest relative change of any summarised number versus B0 (0 = bit-identical, i.e. the code was not reached or the change is exact).
- **Reach.** Static: the groups whose providers execute the changed lines (section 3). Detection rates are reported both over all 25 mutants and over the mutants that reach the group.
- **Rates.** Per group, per category, suite-wide (any test), and for T1–T4 versus T5–T6. Controls are reported separately as false alarms.

## 5. Mutants and predictions

Categories as in the request: timing, units, Jacobian, distribution, equilibrium, calibration, plus controls (equivalent changes). "Predicted catch" lists the groups expected to catch; all other groups are predicted to miss. "Unsure" marks groups without a firm prediction. Confidence H/M/L refers to the whole row.

| ID | Category | Target | Change | Reach | Predicted catch | Unsure | Conf. |
|---|---|---|---|---|---|---|---|
| M01 | timing | `core.py:_make_hh_function` | Euler equation discounts with the current rate r_t instead of r_{t+1}: uc = beta (1+r_t) E[c'^(-1/eis)], Va = c^(-1/eis) | T1 T2 T3 T4 T5 T6 D | T1 T2 T3 | - | H |
| M02 | timing | `models.py:fiscal` | government pays interest at r_{t-1} (ex-ante reading) instead of the ex-post r_t | T5 T6 | T5 | - | H |
| M03 | timing | `models.py:fiscal` | tax rule responds to current debt B_t instead of B_{t-1} | T5 T6 | none | T6 | H |
| M04 | timing | `models.py:ssj_werning_provider` | rate path shifted one period (ex-ante -> ex-post conversion applied to a path that is already ex-post) | T3 | T3 | - | H |
| M05 | timing | `models.py:fiscal` | government budget uses end-of-period debt B_t instead of B_{t-1}: B - (1+r) B - G + Tax | T5 T6 | T5 T6 | - | M |
| M06 | units | `models.py:ssj_tank_provider` | rate path divided by 100 (read as percent although the interface uses decimals) | T2 | T2 | - | H |
| M07 | units | `models.py:AIYAGARI_DEFAULTS` | quarterly depreciation delta = 0.02 in the annual calibration | T4 D | none | T4 | M |
| M08 | units | `core.py:rouwenhorst` | income levels e = exp(s) with E[log e] = 0 instead of E[e] = 1 (mean-income normalisation omitted) | T1 T2 T3 T4 T5 T6 D | none | - | M |
| M09 | jacobian | `models.py:household_jacobians` | default finite-difference step 1e-2 instead of 1e-4 | T1 T5 T6 | T5 T6 | - | H |
| M10 | jacobian | `models.py:ssj_tank_provider` | saver Jacobian one-sided while the GE options stay two-sided | T2 | T2 | - | H |
| M11 | jacobian | `models.py:household_jacobians` | Jacobians returned transposed (row = shock date, column = output date) | T1 T5 T6 | T1 T5 T6 | - | H |
| M12 | jacobian | `models.py:ssj_werning_provider` | cached Js (two-sided, h = 1e-3 b) passed without matching options: SSJ re-differentiates one-sided with h = 1e-4 (F1 type) | T3 | T3 | - | M |
| M13 | distribution | `core.py:transition_matrix` | Young lottery weights swapped between the lower and upper grid point | T1 T2 T3 T4 T5 T6 D | T1 T2 T4 D | T3 T6 | M |
| M14 | distribution | `models.py:cohort_impc_paths` | cohort distribution not divided by its mass | T6 | T6 | - | H |
| M15 | distribution | `models.py:werning_primitives` | Markov matrix transposed (column-stochastic convention) | T3 | none | - | L |
| M16 | distribution | `core.py:_make_hh_function` | borrowing constraint imposed at a_grid[1] instead of a_grid[0] (grid boundary) | T1 T2 T3 T4 T5 T6 D | T4 | - | H |
| M17 | distribution | `models.py:AIYAGARI_DEFAULTS` | asset grid upper bound 20 instead of 200 | T4 D | T4 | D | H |
| M18 | equilibrium | `models.py:fiscal_clearing` | goods-market clearing sign error: Y - C + G | T5 T6 | T5 T6 | - | H |
| M19 | equilibrium | `models.py:goods_clearing` | TANK goods market cleared with savers' consumption only (Y = C_S) | T2 | T2 | - | H |
| M20 | equilibrium | `models.py:fiscal` | government budget constraint without interest: B - B(-1) - G + Tax | T5 T6 | T5 T6 | - | H |
| M21 | equilibrium | `models.py:fiscal_hank` | steady-state taxes omit debt service: Tax_ss = G instead of G + r B | T5 T6 | none | - | M |
| M22 | equilibrium | `models.py:ssj_tank_provider` | saver iMPC matrix taken from the H-type block (copy-paste: J_H instead of J_S) | T2 | T2 | - | H |
| M23 | calibration | `models.py:tank_model` | zero-net-supply target for savers without the division by (1 - lam) | T2 | none | - | M |
| M24 | calibration | `models.py:ssj_aiyagari_r` | borrowing limit -1 instead of 0 | T4 D | none | - | M |
| M25 | calibration | `models.py:fiscal_hank` | beta calibrated to A = 4 B (debt target in quarterly-GDP units) while B = 1 | T5 T6 | T5 | - | H |
| C01 | control | `models.py:AIYAGARI_DEFAULTS` | asset grid refined from 1000 to 1500 points | T4 D | none | - | M |
| C02 | control | `core.py:_make_hh_function` | Va = (1 + r) / c ** (1 / eis) instead of (1 + r) * c ** (-1 / eis) | T1 T2 T3 T4 T5 T6 D | none | - | H |
| C03 | control | `models.py:ssj_werning_provider` | GE target goods market (C = Y) instead of bond market (A = 0) | T3 | none | - | M |

### Mechanism behind each prediction

- **M01** (core.py:_make_hh_function): Steady state unchanged (r constant), so T4/D see nothing. dC/dr is shifted by one date: T1 (closed-form dC/dr, exact at eps=0), T2/T3 (fixed-rate GE IRF vs closed form). T5/T6 build their references from the same household Jacobians.
- **M02** (models.py:fiscal): IKC uses the solver's own dTax path, so the IKC tests pass; households earn r_t on A_{t-1} while the government pays r_{t-1} on B_{t-1}, so the asset market misses by B (dr_t - dr_{t-1}) under rate shocks: T5 Walras test. T6 uses dr = 0 for its Ricardian check and identities that hold for any path.
- **M03** (models.py:fiscal): A different but internally consistent and stable fiscal rule; no test pins the rule. Only T6's Ricardian tail bound 5 (1 + r - phi_B)^T encodes the original decay rate (new rate (1+r)/(1+phi_B) = 0.927 vs 0.92), so T6 may fail at T = 300.
- **M04** (models.py:ssj_werning_provider): dY_t depends on sum_{s>t} dr_s; a one-period shift changes it by O(1) relative.
- **M05** (models.py:fiscal): Steady state still satisfied; dynamics of B explosive (factor phi_B / r = 5). Asset market misses (T5 Walras) and the RA benchmark is no longer Ricardian (T6).
- **M06** (models.py:ssj_tank_provider): GE IRF 100 times too small against the closed form.
- **M07** (models.py:AIYAGARI_DEFAULTS): The no-risk limit r -> 1/beta - 1 does not depend on delta, and monotonicity in risk is generic. Capital is about 3.5 times larger, so the top-of-grid mass check (T4) might trigger for a reason unrelated to the error.
- **M08** (core.py:rouwenhorst): Scales income levels by about exp(sd^2/2). References are built from the reported steady state (C, A) and are scale-free; the factor goes to 1 in the eps -> 0 limits. No test checks the steady-state goods-market residual.
- **M09** (models.py:household_jacobians): T1's provider passes h explicitly; T1's horizon test compares two Jacobians with the same step. T5/T6: Js (h=1e-2) and SSJ's re-differentiation (h=1e-4) disagree, so the GE solution no longer satisfies the IKC built from Js (F1 mechanism).
- **M10** (models.py:ssj_tank_provider): H_U from the one-sided Js, propagation re-differentiated two-sided: goods-market residual ~1e-6 > 1e-10 (test_ge_solution_residual).
- **M11** (models.py:household_jacobians): M_RA has identical rows, its transpose does not (T1). The transposed Js enter the SSJ solve while impulses are re-differentiated (T5, T6).
- **M12** (models.py:ssj_werning_provider): At b = 1e-4 the default step moves income by about b, so propagation and H_U disagree at O(1) in the b-scaled problem.
- **M13** (core.py:transition_matrix): Only FastHetBlock's stationary distribution changes (SSJ's own operator is intact): D compares with SSJ forward iteration; mean assets and so beta, r change (T1, T2, T4). T5 recalibrates beta to A = B under the wrong distribution and stays consistent. Calibration loops may not terminate (timeout counts as caught by error).
- **M14** (models.py:cohort_impc_paths): Group paths scaled by mu_k: sum_k mu_k dC_k != M x (test_wealth_group_split_adds_up).
- **M15** (models.py:werning_primitives): The zero-liquidity reference is built from the provider-reported Pi, and the marginal-saver Euler argument holds for any positive Pi used consistently; the transposed chain keeps unit spectral radius (uniform income-state marginal). Low confidence: Pi.T rows do not sum to one.
- **M16** (core.py:_make_hh_function): States held at a_grid[1] are classified as unconstrained by the Euler diagnostic and show large Euler gaps (T4). T1 has no constrained mass; T2 H types still have dC/dY = 1; T3 shift is ~1e-4 b; T5/T6 consistent.
- **M17** (models.py:AIYAGARI_DEFAULTS): Mass at the top of the grid (test_top_grid_mass_negligible). D compares clipped and unclipped forward operators and may also move.
- **M18** (models.py:fiscal_clearing): IKC with +dG fails (T5), decomposition and multipliers fail (T6).
- **M19** (models.py:goods_clearing): chi = 1 gives the same answer (Y_S = Y), chi = 1.5 does not: T2 IRF and amplification.
- **M20** (models.py:fiscal): IKC holds on the solver's paths; asset market misses by r dB (T5 Walras); taxes cover phi_B/(r+phi_B) = 83% of spending in PV, so the RA multiplier is not 1 (T6).
- **M21** (models.py:fiscal_hank): Steady-state government budget violated (fiscal_res_ss = -r B), but no test checks steady-state residuals; dynamics, A = B and Walras in deviations are unaffected.
- **M22** (models.py:ssj_tank_provider): M_S = I instead of M_RA.
- **M23** (models.py:tank_model): A_H -> 0 as eps -> 0 (H types end at the zero limit), so both targets go to 0; the error vanishes in the limit the test looks at.
- **M24** (models.py:ssj_aiyagari_r): The no-risk limit and the monotonicity regressions hold for any fixed borrowing limit; Euler/KKT diagnostics use the mutated grid consistently.
- **M25** (models.py:fiscal_hank): A_ss != B_ss: in deviations dA - dB picks up (A_ss - B_ss) dr under rate shocks (T5 Walras). T6 compares with an RA at the same (C, A) and is consistent.
- **C01** (models.py:AIYAGARI_DEFAULTS): Refinement within the range the README reports; thresholds are regression thresholds.
- **C02** (core.py:_make_hh_function): Same expression; rounding-level differences only.
- **C03** (models.py:ssj_werning_provider): Walras: dC - dY = (1+r) dA_{t-1} - dA_t with A_ss = 0, so C = Y for all t with dA_{-1} = 0 implies dA = 0; equivalent up to O(h^2) and rounding.

### Summary of the predictions

- Suite-wide: 18 of 25 mutants predicted caught (M01, M02, M04, M05, M06, M09, M10, M11, M12, M13, M14, M16, M17, M18, M19, M20, M22, M25); 7 predicted missed (M03, M07, M08, M15, M21, M23, M24), of which M03 (T6) and M07 (T4) are uncertain.
- Per group (predicted catches): T1 3 (M01, M11, M13); T2 6 (M01, M06, M10, M13, M19, M22); T3 3 (M01, M04, M12); T4 3 (M13, M16, M17); T5 7 (M02, M05, M09, M11, M18, M20, M25); T6 6 (M05, M09, M11, M14, M18, M20); D 1 (M13).
- Main hypothesis on T5/T6: they miss errors inside the household block that change the Jacobians and the GE paths together (M01, M08, M13), and errors that leave a consistent but different model (M03, M21). They catch errors in GE structure and in the consistency between cached Jacobians and the solver (M09, M11, M18), and T5's Walras test is the one check there that compares two separately coded budget constraints (M02, M05, M20, M25).
- Expected misses across the suite: calibration and unit errors whose effect vanishes in the limit being tested or is scale-free (M07, M08, M23, M24), steady-state budget errors (M21), and a specification change that is internally consistent (M03).
- Controls: no outcome changes.

## 6. Amendments

None at the time of the pre-registration commit. Any later change to the mutant list, the predictions, the detection rule or the harness will be listed here with date, reason and whether results had been seen.

- **A1 (2026-10-21, after the first full run; results seen).** The output JSON stored the check summaries of the two baseline runs only; for mutant runs it kept test outcomes and two headline numbers per group. Changed `write_outputs` to store every run's outcomes and check summaries (key `runs`, which replaces `baselines`), so that error sizes can be reported in full. Mutant list, patches, predictions, detection rule and plugin code are unchanged. All 30 runs were repeated with the amended script; the published `out/mutation_v0_1.{json,csv}` are from the repeat, and test outcomes and headline numbers were compared with the first run (result stated in `docs/MUTATION_v0_1.md`).

## 7. Environment at pre-registration

Linux 6.18 x86_64, Python 3.13.16, sequence-jacobian 1.0.0, numpy 2.5.3, scipy 1.18.1, numba 0.68.0 (llvmlite 0.50.0), scipy-openblas 0.3.34, pytest 9.1.1; 4 CPUs, one thread per run. Unmodified suite on this machine: 41 passed, 2 xfailed in 75 s.
