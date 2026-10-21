# v0.3 checks and mutation testing v0.2: plan (pre-registration)

Status: **pre-registered** in commit P on branch `v0_3`. This file, the v0.2 mutant list with its predictions and the runner `scripts/mutation_v0_2.py` were committed and pushed before any v0.3 check code was written and before any mutant was run. The threshold values (section 4) are filled in by commit T, before any mutant is run on either suite. Changes after P are listed under "Amendments" (section 10) with date, reason and whether results had been seen.

Starting point: PR #1 (mutation testing v0.1) was not merged when this branch was created, so `v0_3` starts from `mutation-v0_1` at `aa85365`. Its `hank_limits/`, `tests/` and `pytest.ini` are identical to tag `v0.2.1` (`75544da`).

## 1. Purpose

Mutation testing v0.1 (`docs/MUTATION_v0_1.md`) found that the v0.2.1 suite caught 21 of 25 hand-made errors and missed M12 (cached Jacobians with mismatched options in T3), M21 (steady-state taxes without debt service), M23 (TANK net-supply target without 1 − λ) and M24 (Aiyagari borrowing limit −1); M07 (quarterly δ) was caught only incidentally. Section 10 of that report listed six remedies. This plan does two things:

1. Implement those remedies in the provider style of the suite, as checks V1–V6 (version 0.3, not released).
2. Measure, on new mutants that were generated without knowledge of V1–V6, what the v0.2.1 suite and the v0.3 suite each catch.

What it can show: for a second hand-picked sample of plausible errors, whether the v0.3 checks catch errors that the v0.2.1 suite misses, and at what cost in false alarms on controls. What it cannot show: a detection rate for HANK implementation errors in general (section 9).

## 2. Rules

1. **Scope.** Work only on branch `v0_3`. No merge, no tag, no GitHub Release (the Zenodo integration would mint a DOI). `pyproject.toml` version, `CITATION.cff` and `.zenodo.json` are not changed.
2. **Existing tests.** Tolerances and references are not changed. The one exception is splitting `test_rank_benchmark_is_ricardian` into its RA and HA assertions (V5), with the assertions and tolerances as they are. On the unmodified code, every number in the v0.1 summaries of the existing checks must be **bit-identical** to v0.2.1: the baseline runs of the v0.3 suite are compared with runs B0 and B1 of `out/mutation_v0_1.json` (`--compare-v01` in the runner), and every number must equal its value in one of them (section 8 explains why one number has two values).
3. **Model side** (`models.py`, `core.py`): only **new fields** in the dicts the providers return. No computation is changed, so the v0.1 patches still apply and rule 2 holds. New fields are inserted as new lines wherever possible, so that existing lines (and patch anchors) stay as they are.
4. **Checker side.**
   - `checks.py` does not depend on SSJ.
   - The new checks compute residuals **directly** from primitive quantities the provider reports (Tax, G, r, B, Y, C, A, K, w, e_grid, Π, the income grid, policies, paths). They do not judge on a residual that the provider computed, except V1(a), which records the solver residual the T3 provider already reports, as T2 does; V1(b) computes the same residual directly.
   - Documented specification values live in a new checker module `hank_limits/spec.py`. They are **not** read from `models.py` (`PIH_DEFAULTS`, `AIYAGARI_DEFAULTS`, `FISCAL_DEFAULTS`, ...), so that a mutated default is compared with the documentation, not with itself.
   - Where possible the specification is checked on reported equilibrium quantities rather than on self-reported parameters: firm first-order conditions r + δ = α (K/L)^(α−1) and w = (1 − α)(K/L)^α, the Euler equation at the documented β and EIS, the fiscal rule and the government budget along the path.
   - New fields are **optional**. If a provider does not report a field, the checks that need it are skipped, and the result records which checks were skipped and why.
5. **Thresholds** follow the rule in section 4, which uses only the unmodified baseline values and solver tolerances. Values are committed in T before any mutant is run on either suite. If a threshold has to change after results are seen, the change is recorded as an amendment "after results", and the results under the original value are reported as well.
6. **Order.** P (this plan, the mutants, the predictions, the runner) pushed → v0.3 check code written → T (checks and thresholds) pushed → runs → report. Before P only two things were run: `--check` (patches applied and compiled; no tests) and a smoke test of the unmutated v0.2.1 tree (section 8).
7. **v0.1 material** (`docs/MUTATION_PLAN_v0_1.md`, `docs/MUTATION_v0_1.md`, `scripts/mutation_v0_1.py`, `out/mutation_v0_1.*`) is not changed, with one exception: one dated sentence appended to section 8 of `docs/MUTATION_v0_1.md` on the hetoutput part of F3.

## 3. New checks V1–V6

All quantities are computed in `hank_limits/checks.py` (numpy only) from fields the providers report, and stored under the key `v0_3` of the check result, so the v0.1 part of every result is unchanged. Documented values come from `hank_limits/spec.py`. The tests are in a new file `tests/test_v0_3_checks.py` (group **V** in the runner), except V5, which splits an existing T6 test. Reason for a separate file: the v0.3 groups T1–T6 and D then contain exactly the v0.2.1 tests, so the contribution of the new checks is the group V column, and the comparison of the two suites is not mixed with a change of the old groups. V1 is therefore in group V and not in the T3 file. This departs from "recorded and judged in `check_werning`, as T2 does" only in the location of the test; the quantity is recorded in `check_werning` as asked.

One test function per model and condition, so that a failing test names one condition (lesson of the composite Ricardian test, v0.1 report, section 6). Tests with several quantities (V3, V4(a)) list each quantity in the failure message.

| ID | Test(s) in `tests/test_v0_3_checks.py` | Quantity (max over the check's settings) | Fields used (new fields in *italics*) |
|---|---|---|---|
| V1(a) | `test_v1_werning_ge_target_resid_reported` | provider-reported GE target residual of the T3 solve, max \|asset_mkt\| / b, both incidences | `diag['ge_target_resid']` (exists) |
| V1(b) | `test_v1_werning_ge_target_resid_direct` | max \|J^{A,Y} dY + J^{A,r} dr\| / b, computed in the check | dY (exists), *`J_A`* = {J^{A,Y}, J^{A,r}} |
| V2 | `test_v2_fiscal_ss_government_budget` | \|Tax − G − r B\| / Y | *`ss_detail`*: Tax, G, r, B, Y |
| V2 | `test_v2_fiscal_ss_goods_market` | \|Y − C − G\| / Y | *`ss_detail`*: Y, C, G |
| V2 | `test_v2_fiscal_ss_asset_market` | \|A − B\| / B | *`ss_detail`*: A, B |
| V2 | `test_v2_tank_ss_goods_market` | \|Y − λ C_H − (1 − λ) C_S\| / Y | *`ss_detail`*: Y, λ, C_H, C_S |
| V2 | `test_v2_tank_ss_asset_market` | \|λ A_H + (1 − λ) A_S\| / Y | *`ss_detail`*: λ, A_H, A_S, Y |
| V2 | `test_v2_werning_ss_asset_market` | \|A\| / b | *`ss_detail`*: A, b |
| V2 | `test_v2_werning_ss_goods_market` | \|C − Y\| / Y | *`ss_detail`*: C, Y |
| V2 | `test_v2_aiyagari_ss_capital_market` | \|A − K\| / K | *`ss_detail`*: A, K |
| V2 | `test_v2_aiyagari_ss_firm_focs` | \|r + δ − α (K/L)^(α−1)\| and \|w − (1 − α)(K/L)^α\| / w, documented α, δ; L = Σ π e from the reported chain | *`ss_detail`*: r, w, K, e_grid, Π |
| V2 | `test_v2_aiyagari_ss_goods_market` | \|K^α L^(1−α) − C − δ K\| / Y, documented α, δ | *`ss_detail`*: K, C, e_grid, Π |
| V2 | `test_v2_household_budget[t1, t2, t3, t4, t5]` | aggregate household budget at the steady state, \|C − r A − Σ_e π_e y_e\| / C (per type in T2), π the stationary distribution of the reported Π | *`ss_detail`*: C, A, r, y, Π |
| V3 | `test_v3_income_normalisation[t1, t2, t3, t4, t5]` | (i) Π row-stochastic: max\|Σ_j Π_ij − 1\| and max(0, −min Π); (ii) \|E_π[e] − 1\|; (iii) aggregate income \|Σ π y − Y_ref\| / Y_ref with Y_ref = Y (T1, T3), Y_H and Y_S per type and λ Y_H + (1 − λ) Y_S = Y (T2), w × documented labour 1 (T4), Y − Tax (T5); (iv) T3 countercyclical: \|Σ π e ζ\| | e_grid, Π (T3: exist), *`ss_detail`*: y, Y, Tax, w, Y_H, Y_S |
| V4(a) | `test_v4_spec_parameters[t1, t2, t3, t4, t5]` | reported parameters against `spec.py`: T1 r, EIS, mean assets A = 5, borrowing limit 0.99 × natural limit; T2 λ, β_H, r, EIS, H limit 0, S limit as T1; T3 EIS, limit −b, ζ (acyclical 0; countercyclical ζ_low = 1), r* = 1% from the reported β, Π, e (`ra.zero_liquidity_rstar`); T4 β, EIS, limit 0; T5 r, EIS, B, G, φ_B, limit 0 | *`ss_detail`* |
| V4(b) | `test_v4_aiyagari_euler_at_spec_beta` | Euler errors of the reported policy at the documented β and EIS: mass-weighted mean and max over unconstrained states with mass, KKT violation at binding states (same definitions as `diagnostics.euler_kkt_residuals`, reimplemented in numpy) | *`policy`*: a_grid, a, c, D; *`ss_detail`*: r, Π |
| V4(b) | `test_v4_fiscal_rule_on_path` | max_t \|dTax_t − φ_B dB_{t−1}\| / s, documented φ_B | irfs: dTax, dB (exist) |
| V4(b) | `test_v4_government_budget_on_path` | max_t \|dB_t − (1 + r) dB_{t−1} − B dr_t − dG_t + dTax_t\| / s, documented r, B | irfs: dB, dr, dG, dTax (exist) |
| V4(b) | `test_v4_real_rate_rule_on_path` | rule-based shocks: max_t \|dr_t − φ_y dY_t − drshock_t\| / s, documented φ_y | irfs: dr, dY (exist), *`drshock`* |
| V4(b) | `test_v4_goods_market_on_path` | max_t \|dY_t − dC_t − dG_t\| / s | irfs: dY, dC, dG (exist) |
| V5 | `tests/test_t6_differencing.py`: `test_rank_benchmark_is_ricardian` (RA assertion only), `test_ha_multiplier_exceeds_one` (HA assertion only) | as in v0.2.1 | as in v0.2.1 |
| V6 | `test_v6_f3_hetinput_hetoutput_two_sided` (strict xfail, SSJ code only) | (i) the hetinput derivative that `HetBlock` uses, `CombinedExtendedFunction.differentiable(ss, h, True).diff({x: 1})`, equals the central difference; (ii) `HetBlock.jacobian(..., twosided=True)` of a hetoutput equals the one obtained when the stored `twosided` default is honoured | none (SSJ `hh_sim` example block) |

s = max \|exogenous input path\| of the shock (0.01 here). The path checks V4(b) use the T5 solve with matched options and iterative refinement (`ssj_consistent_refined`), as the Walras test does, so that F1 and F2 do not enter them. V2, V3 and V4(a) on the fiscal model also use that setting (the steady state is the same in all three).

Section-10 remedies of the v0.1 report and where they go: 1 → V1; 2 → V2; 3 → V3; 4 → README (T4's error detection comes from its diagnostics) and V2/V4 for Aiyagari; 5 → V5 and README (the T6 matrix identity is a linear-algebra check); 6 → V4. V6 pins finding F3.

Documented values in `spec.py` (from the README "Exact setting" paragraphs): T1/T2 savers r = 0.03, EIS 0.5, mean assets 5, borrowing limit 0.99 × (−y_min/r), Y = 1; T2 λ = 0.3, β_H = 0.95, H borrowing limit 0; T3 EIS 0.5, r* = 0.01, borrowing limit −b, ζ_low = 1 in the countercyclical case; T4 β = 0.96, α = 0.36, δ = 0.08, EIS 0.5, borrowing limit 0, labour 1; T5 r = 0.02, EIS 0.5, B = 1, G = 0.2, φ_B = 0.1, φ_y = 0.25, borrowing limit 0. Grid sizes, grid tops and the income-process parameters are numerical choices or enter through the reported chain, and are not pinned.

## 4. Thresholds

**Rule.** For each quantity q of a new test:

    threshold(q) = roundup_1sig( max(100 × b_max(q), 10 × τ(q)) )

- b_max(q): the largest value of q over all settings of the check (all eps, b, T, labels, incidences) in **two** runs of the unmodified v0.3 suite at the T candidate, before T is committed;
- τ(q): the tolerance of the numerical step that determines q, from the table below;
- roundup_1sig: round up to one significant digit (3.2e-11 → 4e-11).

| Class | Quantities | τ |
|---|---|---|
| exact up to rounding | Π row sums and sign, E_π[e], Σ π e ζ, parameter pins, r* from β, firm FOCs (K is computed from them), GE target residuals and path residuals (linear solves) | 1e-14 |
| steady-state aggregates | household budget, aggregate income, goods markets | 1e-14 (distribution sweep tolerance 1e-15; the identities hold per state) |
| calibrated market clearing | T5 A − B (brentq on β, xtol 1e-14), T4 A − K (brentq on r, xtol 1e-14), T2 λ A_H + (1 − λ) A_S (brentq on log κ, xtol 1e-13), T3 A / b (brentq on r, xtol 1e-15), T1 mean assets (xtol 1e-13) | the root finder's xtol |
| Euler errors | V4(b) Euler and KKT | backward tolerance 1e-12 (`SS_OPTS`) |

The values of b_max, τ and the thresholds are written into the table below and into `tests/test_v0_3_checks.py` in commit T. Existing tests keep their tolerances.

| Quantity | b_max (run 1, run 2) | τ | Threshold |
|---|---|---|---|
| (filled in by T) | | | |

## 5. New provider fields (all optional)

Documented in `hank_limits/interface.py` in T.

- `RankLimitProvider`: `ss_detail` = {Y, C, A, r, eis, beta, a_min, e_grid, Pi, y}; y = income by income state at the steady state.
- `TankProvider`: `ss_detail` = {lam, r, eis, Y, Y_H, Y_S, C_H, C_S, A_H, A_S, beta_H, beta_S, a_min_H, a_min_S, e_grid_H, Pi_H, y_H, e_grid_S, Pi_S, y_S}.
- `WerningProvider`: `ss_detail` = {b, Y, C, A, r, eis, beta, a_min, y}; `J_A` = {'Y': J^{A,Y}, 'r': J^{A,r}} (T × T, the Jacobians passed to the solve).
- `AiyagariProvider`: `ss_detail` = {r, w, K, A, C, beta, eis, a_min, e_grid, Pi, y}; `policy` = {a_grid, a, c, D}.
- `IKCProvider`: `ss_detail` = {Y, C, A, B, G, Tax, r, eis, beta, phi_B, a_min, e_grid, Pi, y}; each `irfs[label]` gains `drshock` (the exogenous rate-shock path; zeros when there is none).

All of these are read from objects the providers already compute (the SSJ steady-state dict, the calibration dict, the Jacobians and paths passed to or returned by the solve). No computation is changed.

## 6. New mutants v0.2

### Generation

A separate agent (a fresh general-purpose subagent with no access to this session's context) generated the candidates. It received copies of `hank_limits/models.py`, `core.py` and `interface.py` (byte-identical to `v0.2.1`; SHA-256 in appendix A) and a list of the v0.1 categories and the 28 v0.1 changes with their exact patches (no predictions, no results). It did not receive the v0.1 report, its section 10, the tests, `checks.py` or this plan, and it was told not to open any other file. Whether it kept to that instruction was not audited. The prompt and the agent's complete output are in appendix A.

Timing (UTC, 2026-10-21): prompt written 21:37:46; the check design of sections 3 and 4 written 21:48–21:51; candidates delivered 21:58:36. Sections 3 and 4 were not changed after the candidates were received.

### Selection rule and selection

Conditions: (a) the patch applies cleanly (each `old` occurs exactly once) and compiles; (b) a plausible implementation error (mutants) or an equivalent change (controls); (c) not a duplicate of a v0.1 change (same line and same change); (d) the changed code is executed by at least one test group. Number: mutants are taken in the agent's order until **20** qualify (the upper end of the allowed 16–20), controls until **3** qualify. The cut at 20 and 3 was fixed when the list was received; the list had been read at that point. Candidates are not chosen by whether a new check would catch them.

Result: N01–N20 and K01–K03 all meet (a)–(d), so the selection is exactly the first 20 mutants and the first 3 controls. Not used (beyond the cut, not for failing a condition): N21–N26, K04, K05. Notes on (c): N18 changes M23's line but flips a sign (M23 removed a factor); N19 changes the `Va` line that M01 and C02 change, with a different change; N10 is an options error like M12 but in `_het_opts`, and drops only the step size. `--check` before P: all 23 patches apply and compile on the v0.2.1 tree, and on the current `v0_3` tree (identical model code), as do the 28 v0.1 patches.

Categories of the 20 mutants (as labelled by the agent): timing 3 (N04, N08, N12), units 3 (N01, N07, N19), jacobian 3 (N03, N10, N15), distribution 4 (N02, N06, N09, N20), equilibrium 3 (N11, N14, N16), calibration 4 (N05, N13, N17, N18).

### Mutants, reach and predictions

"Predicted catch" lists the groups expected to catch; all other groups are predicted to miss, except those under "Unsure" (no firm prediction). V exists only in the v0.3 suite. Confidence H/M/L refers to the whole row. The same data are in `MUTANTS` in `scripts/mutation_v0_2.py`.

| ID | Category | Target | Change | Reach | Predicted catch, v0.2.1 | Predicted catch, v0.3 | Unsure (v0.2.1 / v0.3) | Conf. |
|---|---|---|---|---|---|---|---|---|
| N01 | units | `models.py:ssj_aiyagari_r` | The log-income sd passed to the provider is multiplied by sqrt(1 - rho^2) before discretisation, so the chain has unconditional sd sd*sqrt(1 - rho^2) instead of sd. | T4 D | none | none | T4 / T4 | M |
| N02 | distribution | `core.py:transition_matrix` | The exogenous step of the Young-lottery chain uses kron(Pi, I) instead of kron(Pi', I), so mass moves with the rows of Pi and the stationary income marginal becomes uniform instead of the Markov chain's stationary distribution. | T1 T2 T3 T4 T5 T6 D | T5 D | T5 D V | T1 T2 T3 T6 / T1 T2 T3 T6 | M |
| N03 | jacobian | `models.py:cohort_impc_paths` | The central finite difference of the cohort consumption paths is divided by h instead of 2h, doubling every cohort iMPC path. | T6 | T6 | T6 | - / - | H |
| N04 | timing | `models.py:real_rate_rule` | The real-rate rule responds to last period's output gap Y_{t-1} instead of the current gap Y_t. | T5 T6 | none | V | T5 / T5 | H |
| N05 | calibration | `models.py:werning_primitives` | The countercyclical incidence is normalised so that sum_i pi_i zeta_i = 0 instead of sum_i pi_i e_i zeta_i = 0, so aggregate household income no longer equals Y off the steady state. | T3 | none | V | - / - | M |
| N06 | distribution | `core.py:FastHetBlock.forward_steady_state` | The post-shock distribution D is set equal to the beginning-of-period distribution Dbeg, so steady-state aggregates weight policies indexed by e_t with the e_{t-1} distribution. | T1 T2 T3 T4 T5 T6 D | T5 | T5 V | T1 T2 T3 T6 / T1 T2 T3 T6 | M |
| N07 | units | `models.py:FISCAL_DEFAULTS` | The fiscal-HANK calibration uses the quarterly real rate r = 0.005 inside an otherwise annual calibration (B = 1 annual GDP, annual persistence). | T5 T6 | none | V | T5 / T5 | M |
| N08 | timing | `models.py:fiscal_clearing` | Asset-market clearing compares households' end-of-period assets A_t with last period's debt B_{t-1} instead of B_t. | T5 T6 | T5 | T5 | - / - | H |
| N09 | distribution | `models.py:cohort_impc_paths` | Wealth cohorts are built from the post-shock distribution D instead of the beginning-of-period distribution Dbeg, so the cohort weights and initial conditions are wrong and the cohort paths no longer sum to M x. | T6 | T6 | T6 | - / - | M |
| N10 | jacobian | `models.py:_het_opts` | The options passed to SSJ's impulse_linear drop the step size, so the HetBlock is re-differentiated with SSJ's default h = 1e-4 while the supplied Js use the caller's h (1e-3 b in the Werning provider). | T2 T3 T5 T6 | none | V | - / - | H |
| N11 | equilibrium | `models.py:tank_incidence` | H-type income is Y (Y/Yss)^chi, so chi acts as an excess elasticity (total elasticity 1 + chi) instead of the elasticity of Y_H with respect to Y. | T2 | T2 | T2 | - / - | H |
| N12 | timing | `core.py:transition_matrix` | The one-period operator is composed as X @ lot (savings lottery first, then income shock) instead of lot @ X, so the vector it returns as Dbeg is really the post-shock distribution. | T1 T2 T3 T4 T5 T6 D | T5 D | T5 D V | T1 T2 T3 T6 / T1 T2 T3 T6 | M |
| N13 | calibration | `models.py:pih_household` | The final discount factor is set to 1/(1+r) - kappa instead of (1 - kappa)/(1+r), which does not match the parameterisation used in the calibration root-find, so mean assets miss A_target and kappa no longer equals 1 - beta (1+r). | T1 T2 D | none | V | - / - | M |
| N14 | equilibrium | `models.py:tank_aggregate` | Aggregate consumption weights the H type by 1 - lam and savers by lam, the reverse of the population shares. | T2 | T2 | T2 | - / - | H |
| N15 | jacobian | `models.py:household_jacobians` | The household Jacobians are computed with SSJ's default one-sided differences instead of two-sided ones, while the GE options still re-differentiate two-sided. | T1 T5 T6 | T1 T5 T6 | T1 T5 T6 V | - / - | H |
| N16 | equilibrium | `models.py:solve_impulse_linear_refined` | After solving for the unknowns, outputs are propagated from dU alone, so the direct effects of the exogenous shocks (G, r, rshock) are missing from all reported paths. | T5 T6 | T5 | T5 V | - / - | H |
| N17 | calibration | `models.py:tank_model` | The saver block in the TANK model is given the H type's asset grid (limit 0, top 50) instead of its own near-natural-limit grid, while beta_S was calibrated on the latter. | T2 | T2 | T2 V | - / - | H |
| N18 | calibration | `models.py:tank_model` | The savers' asset target has the wrong sign (+lam A_H/(1-lam)), so both types hold positive bonds and the bond is not in zero net supply. | T2 | none | V | - / - | M |
| N19 | units | `core.py:_make_hh_function` | The parameter eis is used as the CRRA coefficient in both the Euler inversion and the envelope condition, so the household's actual EIS is 1/eis. | T1 T2 T3 T4 T5 T6 D | T1 T2 T3 T4 | T1 T2 T3 T4 V | T5 / T5 | H |
| N20 | distribution | `core.py:mass_at_constraint` | The constrained share weights the policy mask (indexed by e_t) with the beginning-of-period distribution Dbeg instead of D. | T1 T2 T3 T4 T5 T6 D | none | none | - / - | H |
| K01 | control | `models.py:ssj_aiyagari_r` | The wage is computed from zero profits, w = K^alpha - (r + delta) K, instead of the marginal product (1 - alpha) K^alpha. | T4 D | none | none | - / - | H |
| K02 | control | `core.py:transition_matrix` | The redundant np.clip of the policy before interpolate_coord_robust is removed, so the out-of-grid handling relies on the existing clip of the lottery weights. | T1 T2 T3 T4 T5 T6 D | none | none | - / - | H |
| K03 | control | `core.py:_make_hh_function` | The endogenous-grid step interpolates on beginning-of-period assets (c + a' - y)/(1+r) evaluated at a_grid instead of on cash on hand c + a' evaluated at coh. | T1 T2 T3 T4 T5 T6 D | none | none | - / - | M |

### Mechanism behind each prediction

- **N01** (`models.py:ssj_aiyagari_r`): The chain gets unconditional sd sd*sqrt(1-rho^2) (0.44 sd at rho = 0.9). The no-risk limit, its O(eps^2) rate and the monotonicity in sd hold for any rescaling of risk. Monotonicity in rho at sd = 0.4 now compares effective sd 0.40, 0.38, 0.32, 0.17 for rho = 0, 0.3, 0.6, 0.9, so persistence and risk size move r in opposite directions (T4 unsure). Steady-state identities, E[e] = 1 and the pins hold; the income-process sd is not pinned (V miss).
- **N02** (`core.py:transition_matrix`): FastHetBlock's chain moves mass with Pi instead of Pi': its stationary Dbeg has a uniform income marginal, and SSJ's D = Pi' Dbeg is not stationary under SSJ's own law of motion. D compares with SSJ's forward iteration (catch). End-of-period and next beginning-of-period mean assets differ, so households earn interest on a different stock than the government pays on: T5 Walras under rate shocks (as M13). T4: brentq still clears the market on the mis-aggregated A, and the Euler errors use the policy (miss). V2: household budget C - r A - sum(pi y) != 0.
- **N03** (`models.py:cohort_impc_paths`): Cohort paths doubled: sum_k mu_k dC_k = 2 M x (test_wealth_group_split_adds_up). No V test reads t6.
- **N04** (`models.py:real_rate_rule`): A different but internally consistent rule: IKC and Walras hold on the solver's paths, the fixed-rate shocks and T6 (fixed r) are untouched. The F2 pin may turn XPASS because the rule-based system changes (T5 unsure). V4(b) checks the documented contemporaneous rule on the path: residual phi_y (dY_{t-1} - dY_t).
- **N05** (`models.py:werning_primitives`): zeta = (1, -1/3, -1/3), so sum pi e zeta != 0. The zero-liquidity reference is built from the reported zeta and needs only the marginal saver's Euler equation, so T3 still converges to it; delta stays > 1 and the deviation from RANK > 0.25. V3 checks sum pi e zeta = 0.
- **N06** (`core.py:FastHetBlock.forward_steady_state`): SSJ aggregates with Dbeg in place of D (a_{t-1} paired with the wrong income state). Dbeg itself is correct, so the D audits pass, and eps = 0 is unchanged. A = sum Dbeg a' differs from the mean beginning-of-period assets: Walras residual under rate shocks (T5). V2: C - r A - sum(pi y) = (1+r)(A_beg - A) != 0.
- **N07** (`models.py:FISCAL_DEFAULTS`): A consistent model at r = 0.005: IKC and Walras hold; T6's Ricardian bound 5 * 0.92^T is looser than the new tail 0.905^T, and the HA multiplier stays > 1. The F2 pin may turn XPASS (T5 unsure). V4(a) pins r = 0.02; V4(b) checks the budget at the documented r.
- **N08** (`models.py:fiscal_clearing`): Only the reported asset_mkt changes (it is not a target): dA_t - dB_{t-1} = dB_t - dB_{t-1} != 0 (T5 Walras). V2 uses A and B at the steady state, where B(-1) = B, and no V test reads asset_mkt (V miss).
- **N09** (`models.py:cohort_impc_paths`): Cohorts start from the post-shock D instead of Dbeg, so sum_k mu_k dC_k != M x by more than 1e-4 (expected, not computed). No V test reads t6.
- **N10** (`models.py:_het_opts`): Only T3 uses h != 1e-4: propagation re-differentiated two-sided at h = 1e-4, H_U from Js at 1e-3 b. Like M12 (one-sided; missed by T3) but smaller, so T3 misses. T2, T5, T6 are bit-identical. The T3 asset-market residual rises from ~1e-18 to well above its threshold (V1(a), V1(b)).
- **N11** (`models.py:tank_incidence`): Elasticity of Y_H is 1 + chi: GE IRF vs TANK and the amplification factor (1.6 instead of 1.273) fail. The steady state is unchanged and lam Y_H + (1 - lam) Y_S = Y holds there (V miss).
- **N12** (`core.py:transition_matrix`): The reversed chain's stationary vector is the post-shock distribution, used as Dbeg; SSJ then applies Pi' again. D: SSJ's forward iteration converges to a different Dbeg (catch). Mean-asset mismatch between end and beginning of period: T5 Walras; V2 household budget.
- **N13** (`models.py:pih_household`): Actual impatience 1 - beta (1+r) = (1+r) kappa, 3% above the reported kappa, so mean assets miss the target by a few percent. The closed forms use the reported A and the kappa band the reported kappa (T1 miss); the TANK net-supply error vanishes with A_H (T2 miss, as M23). V4(a): T1 mean assets != 5; V2: TANK lam A_H + (1 - lam) A_S != 0.
- **N14** (`models.py:tank_aggregate`): The goods market clears with weights (1 - lam, lam): GE IRF and amplification fail. The V2 TANK goods market is computed directly from C_H, C_S with the reported lam and holds (V miss).
- **N15** (`models.py:household_jacobians`): One-sided Js: O(h) ~ 5e-5 error in dC/dr (T1 exactness at eps = 0 requires < 1e-7). T5/T6: H_U from one-sided Js, propagation two-sided (F1 mechanism): IKC and decomposition errors ~1e-5. V4(b): goods market on the path.
- **N16** (`models.py:solve_impulse_linear_refined`): Refined solve: the unknowns are right, the outputs are propagated without the shocks. IKC uses dY, dTax and the input paths (holds); Walras under the r and rshock shocks fails (dA misses the direct rate effect). T6 uses only dY, dG, dTax, dr (miss). V4(b): dY - dC - dG != 0 under rate shocks.
- **N17** (`models.py:tank_model`): Savers on the H grid (limit 0): many savers constrained, M_S no longer tends to M_RA. V4(a): saver borrowing limit; V2: TANK asset market.
- **N18** (`models.py:tank_model`): Net supply lam A_H + (1 - lam) A_S = 2 lam A_H, which vanishes with A_H as eps -> 0 (as M23): T2 miss. V2 catches it at the eps where A_H > 0.
- **N19** (`core.py:_make_hh_function`): Actual EIS 2 instead of 0.5. T1: dC/dr at eps = 0 depends on the EIS (closed form at the reported 0.5); T2: GE IRF; T3: the zero-liquidity rate leaves the root bracket (error, per the generating agent's note); T4: Euler errors of the diagnostic at the reported EIS. T5/T6 internal checks hold; F2 pin unsure; D audits are distribution-only (miss). V: Euler at the documented beta and EIS (V4(b)); the t3 fixture error also fails the V tests on t3.
- **N20** (`core.py:mass_at_constraint`): Diagnostic only: the constrained share is weighted with Dbeg. The T1 slack test (< 1e-5; baseline 2.9e-6) still passes; no other test and no V check uses it.
- **K01** (`models.py:ssj_aiyagari_r`): w = K^alpha - (r + delta) K equals (1 - alpha) K^alpha up to rounding, because K is computed from the same FOC.
- **K02** (`core.py:transition_matrix`): Weights outside [0, 1] are clipped afterwards, which gives the same lottery as clipping the policy.
- **K03** (`core.py:_make_hh_function`): Linear interpolation is invariant to the per-state affine map coh -> (coh - y)/(1+r); rounding-level differences only. The new thresholds have a floor of 10 tau.

### Summary of the predictions

- v0.2.1 suite: 12 of 20 mutants predicted caught (N02, N03, N06, N08, N09, N11, N12, N14, N15, N16, N17, N19); predicted missed: N01, N04, N05, N07, N10, N13, N18, N20, of which N01, N04, N07 have an unsure group.
- v0.3 suite: 18 of 20 predicted caught; predicted missed: N01, N20. Gain predicted from V alone: N04, N05, N07, N10, N13, N18.
- Controls: no outcome change in either suite.
- Main hypothesis: the v0.3 checks add the errors that leave a consistent model with the wrong specification or the wrong steady state (N04, N05, N07, N13, N18), and an options mismatch in T3 of the M12 type (N10). They do not add errors that only affect the T2 GE solve, the T6 cohort simulation, a reported path that is not a target (N08), or a diagnostic (N20).

## 7. Predictions for the v0.1 mutants on the v0.3 suite (design iii)

For T1–T6 and D the v0.3 groups contain the v0.2.1 tests, so the prediction is the v0.1 outcome (`out/mutation_v0_1.json`); the T6 split does not change a group-level outcome. Predictions for group V (also in `PRED_V01_V` in the runner):

| ID | V | Conf. | Mechanism |
|---|---|---|---|
| M01 | miss | H | Steady state unchanged; the GE paths satisfy the (unchanged) fiscal block and market clearing; T3 solve consistent. |
| M02 | catch | H | Government budget on the path at the documented timing misses by B (dr_{t−1} − dr_t) under rate shocks (V4(b)). |
| M03 | catch | H | Fiscal rule on the path: dTax_t − φ_B dB_{t−1} = φ_B (dB_t − dB_{t−1}) ≠ 0 (V4(b)). |
| M04 | catch | H | V1(b) uses the rate path the check passed; the solve used the shifted path, so J^{A,Y} dY + J^{A,r} dr ≠ 0. V1(a) (solver's own residual) passes. |
| M05 | catch | H | Government budget on the path (V4(b)); explosive debt. |
| M06 | miss | H | The T2 V checks are steady-state checks; the steady state is unchanged. |
| M07 | catch | H | Firm FOC and goods market at the documented δ (V2 Aiyagari). |
| M08 | catch | H | E_π[e] ≠ 1 (V3) wherever risk is positive. |
| M09 | catch | H | Js (h = 1e-2) and re-differentiated propagation disagree, so the refined solve leaves goods-market and budget residuals on the path of order 1e-4 (V4(b)). |
| M10 | miss | H | T2 V checks are steady-state checks. |
| M11 | catch | H | Transposed Js in H_U: goods-market residual on the T5 path (V4(b)). |
| M12 | catch | H | GE target residual of the T3 solve ~2e-7 (v0.1 scratch diagnostic) against a baseline near 1e-18 (V1(a), V1(b)). |
| M13 | catch | H | The swapped lottery does not preserve mean assets, so C − r A − Σ π y = (1 + r)(A_beg − A) ≠ 0 at the steady state (V2 household budget), in every model. |
| M14 | miss | H | Cohort paths are T6 only; no V test reads t6. |
| M15 | catch | H | Π^T is not row-stochastic (V3, T3). |
| M16 | catch | H | Euler check at the documented β classifies states held at a_grid[1] as unconstrained (V4(b) Aiyagari), as T4's own diagnostic did. The reported limit a_grid[0] is unchanged, so V4(a) passes. |
| M17 | ? | L | Steady-state identities hold unless policies leave the grid top with mass; the clipped lottery then breaks the household budget. Not predicted. |
| M18 | catch | H | Goods market on the T5 path at the documented sign: residual 2 dG (V4(b)). The steady-state goods market Y − C − G is computed directly and holds. |
| M19 | miss | H | T2 V checks are steady-state checks; the steady state does not use goods clearing. |
| M20 | catch | H | Government budget on the path: residual r dB_{t−1} + B dr_t (V4(b)). |
| M21 | catch | H | Steady-state government budget: Tax − G − r B = −r B = −0.02 (V2). |
| M22 | miss | H | T2 V checks are steady-state checks. |
| M23 | catch | H | TANK asset market λ A_H + (1 − λ) A_S = λ² A_H ≠ 0 at the eps where A_H > 0 (V2). |
| M24 | catch | H | Reported borrowing limit −1 against the documented 0 (V4(a)). |
| M25 | catch | H | Fiscal asset market A − B = 3 B (V2). |
| C01 | miss | M | Grid refinement within the reported range; thresholds are 100 × baseline. |
| C02 | miss | H | Rounding-level change. |
| C03 | catch (false alarm) | H | The asset-market residual of the goods-market solve is 1.1e-4 at b = 0.1 because of F3 (v0.1 report, section 8), far above a baseline near 1e-18 (V1(a), V1(b)). Reported as a false alarm; the threshold is not adjusted. |

Summary: V is predicted to catch 18 of the 25 v0.1 mutants (M17 unsure), to miss M01, M06, M10, M14, M19, M22, and to give one false alarm (C03).

## 8. Run design

- **Designs.** (i) v0.2 mutants and controls × v0.2.1 suite (23 runs); (ii) the same × v0.3 suite (23 runs); (iii) the 28 v0.1 mutants and controls × v0.3 suite (28 runs). With two baselines per suite, 78 runs in all.
- **Suites.** v0.2.1 = the tree of tag `v0.2.1` (`75544da`); v0.3 = the tree of commit T. The runner extracts `hank_limits/`, `tests/`, `pytest.ini` with `git archive`, so an uncommitted change cannot enter a run, and it records the commit and tree hashes.
- **Baselines.** Two unmutated runs per suite (B0, B1). Detection is relative to the same suite's B0.
- **Detection rule, time limits, parallelism:** as in v0.1. A test catches a mutant when its outcome differs from B0 (passed → failed or error; xfailed → failed or XPASS(strict)); a group catches when any of its tests does; 300 s per fixture or test; 4 runs in parallel, one BLAS thread each.
- **Groups.** T1–T6 and D as in v0.1; V = `tests/test_v0_3_checks.py` (v0.3 only). In v0.3, T6 has 7 tests (V5 split). "Reached" for T1–T6 and D as in v0.1 (static, from the code map); V is reached whenever any of T1–T5 is, because the V tests read the t1–t5 check results.
- **Bit identity (rule 2).** The runner compares the v0.1 part of each baseline summary with runs B0 and B1 of `out/mutation_v0_1.json`. In the smoke test before P, the v0.2.1 tree reproduced all 276 summarised numbers bit for bit except one, and that one (`t5.ssj_inconsistent_options.walras_max.rshock_r_rule`, on the ill-conditioned F2 path) took the value of v0.1 B1 in one run and of v0.1 B0 in the other (relative difference 1.4e-11). That number is known to take two values (v0.1 report, section 2). Rule 2 is therefore applied as: every summarised number equals its value in v0.1 B0 or in v0.1 B1; the comparison with each is reported.
- **Re-anchoring.** If a v0.1 or v0.2 patch does not apply to the v0.3 tree because a line it anchors on was changed by a new field, the same change is applied at the new location (`REANCHOR` in the runner) and listed in the report. No other change to a patch is allowed.
- **Runner and outputs.** `scripts/mutation_v0_2.py` (copied from the v0.1 runner; plugin, detection rule and v0.1 summaries unchanged). Outputs: `out/mutation_v0_2.json` (environment, suite commits and trees, runner blob and whether it equals the committed one, mutants with patches, every run's outcomes and summaries, analyses) and `out/mutation_v0_2.csv` (one row per design × mutant × group).

## 9. What can and cannot be claimed

- **(iii) is not an independent evaluation.** The v0.3 checks were designed after seeing which v0.1 mutants were missed. Catching M12, M21, M23, M24 there shows only that the remedies do what they were built for.
- **The fair comparison is (i) against (ii)**: the same independently generated mutants, run on the suite without and with the new checks. Even so, the mutants are a hand-picked sample (by an agent, from a prompt written by the author of the checks), not draws from a distribution of real errors. Detection counts describe this sample only.
- The agent that generated the mutants saw the model code and the v0.1 mutant list, not the tests, the v0.1 report or this plan. It did know the six categories. The selection (section 6) uses its order and objective conditions only.
- Controls show whether the new checks alarm on three equivalent changes. That is not a false-alarm rate.
- The new checks pin the documented calibration. A code under test with a different calibration needs its own `spec` values; the pins are regression checks for this suite, like the existing tolerances.

## 10. Amendments

None at the time of the pre-registration commit.

## Appendix A. Mutant generation: prompt, inputs and complete output

### A.1 Inputs given to the agent

Directory with four files (SHA-256):

```
03b6bcb393bd3a7feeae96b0c3437058a06db7abaf4714f5f58098b143b29b08  models.py     (= hank_limits/models.py at v0.2.1)
39db2e80d4b093a0f97c756edf63348a4673a226eba90dc5a01f3c1c942fb0c5  core.py       (= hank_limits/core.py at v0.2.1)
04e8815c59d92cdaa8a507bf5f46f07a0d287ac5c3d6d5bd479022e257046c0a  interface.py  (= hank_limits/interface.py at v0.2.1)
4abc5aca6ab65a7fd01426b9c1892ea754e02a2d956e196cd581f56e164a453f  v0_1_mutants.md (A.3)
```

### A.2 Prompt (verbatim)

The agent received this text as its whole task (it started without this session's context).

```text
You are helping with a mutation-testing study of a small research code base: heterogeneous-agent
New Keynesian (HANK) models built on the Python package sequence-jacobian (SSJ) 1.0.0, which is
installed in this environment. Your job is to propose new mutants: small, plausible
implementation errors that a careful economist-programmer could realistically make in this code.
You work independently: you are not told how the code is tested, and you must not try to find out.

INPUTS. Read only these files, all in the directory INPUT_DIR (given below):
  - models.py  : model builders and "providers" (functions ssj_*), annual calibrations
  - core.py    : household block (backward step), stationary distribution, grids, income
  - interface.py : what each provider returns
  - v0_1_mutants.md : the six categories and the 28 changes used in an earlier round
Do not open, list or search any other file or directory on this machine (in particular nothing
under /home/user), and do not use the network. You may read the installed sequence_jacobian
package source if you need to understand an SSJ call. Scratch work goes in INPUT_DIR/work only.

TASK.
1. Propose 26 mutant candidates and 5 control candidates.
   - Mutants: each is one small, plausible implementation error in models.py or core.py. Cover
     all six categories (timing, units, jacobian, distribution, equilibrium, calibration), at
     least three candidates per category. Prefer errors that run and give wrong numbers over
     errors that raise at once. An error may touch more than one line if a real programmer would
     make it as one mistake (e.g. a consistent sign convention change), but keep it small.
   - Not a duplicate of v0.1: not the same change on the same line as any entry in
     v0_1_mutants.md. A different error in the same function is fine.
   - Controls: changes that a reviewer would accept as not changing the model or its results
     beyond rounding: algebraically equivalent rewrites, harmless refactors, or a modest
     numerical refinement. Also not duplicates of the v0.1 controls.
2. Order: list the mutant candidates in your order of preference, most realistic first, and the
   controls likewise. The study takes them in your order and cuts at a fixed number, so the
   order matters. Do not reorder by how easy or hard you think they are to detect.
3. Each candidate must be a mechanical text patch: a list of (file, old, new) where `old` occurs
   EXACTLY ONCE in that file (exact whitespace, including indentation and newlines) and the
   patched file must byte-compile. Verify this yourself with a short Python script that copies
   the two files into INPUT_DIR/work, applies each candidate separately, checks the occurrence
   count and runs py_compile. You do not need to run the models (that is slow); do not run them.

OUTPUT. Write INPUT_DIR/work/candidates.json and also return it in full as your final message,
followed by the output of your verification script. JSON format:
{"mutants": [ {"id": "N01", "category": "...", "target": "models.py:function_name",
               "change": "one sentence: what is wrong",
               "why_plausible": "one sentence: how a programmer would make this mistake",
               "patches": [{"file": "models.py", "old": "...", "new": "..."}]}, ... ],
 "controls": [ {"id": "K01", "category": "control", "target": "...", "change": "...",
                "why_equivalent": "one sentence", "patches": [...]}, ... ]}
Use ids N01..N26 for mutants and K01..K05 for controls in your order of preference. Do not say
which tests or checks would catch a candidate (you do not have the tests).

INPUT_DIR = /tmp/claude-0/-home-user-hank-limits/dac78a62-da54-5c22-9fc7-31ee064ed8c5/scratchpad/mutgen
```

### A.3 `v0_1_mutants.md` (verbatim)

````markdown
# Mutants used in the earlier round (v0.1)

Categories: timing, units, jacobian, distribution, equilibrium, calibration; plus controls (changes that should not change any result).

| ID | Category | Target (file:function) | Change |
|---|---|---|---|
| M01 | timing | core.py:_make_hh_function | Euler equation discounts with the current rate r_t instead of r_{t+1}: uc = beta (1+r_t) E[c'^(-1/eis)], Va = c^(-1/eis) |
| M02 | timing | models.py:fiscal | government pays interest at r_{t-1} (ex-ante reading) instead of the ex-post r_t |
| M03 | timing | models.py:fiscal | tax rule responds to current debt B_t instead of B_{t-1} |
| M04 | timing | models.py:ssj_werning_provider | rate path shifted one period (ex-ante -> ex-post conversion applied to a path that is already ex-post) |
| M05 | timing | models.py:fiscal | government budget uses end-of-period debt B_t instead of B_{t-1}: B - (1+r) B - G + Tax |
| M06 | units | models.py:ssj_tank_provider | rate path divided by 100 (read as percent although the interface uses decimals) |
| M07 | units | models.py:AIYAGARI_DEFAULTS | quarterly depreciation delta = 0.02 in the annual calibration |
| M08 | units | core.py:rouwenhorst | income levels e = exp(s) with E[log e] = 0 instead of E[e] = 1 (mean-income normalisation omitted) |
| M09 | jacobian | models.py:household_jacobians | default finite-difference step 1e-2 instead of 1e-4 |
| M10 | jacobian | models.py:ssj_tank_provider | saver Jacobian one-sided while the GE options stay two-sided |
| M11 | jacobian | models.py:household_jacobians | Jacobians returned transposed (row = shock date, column = output date) |
| M12 | jacobian | models.py:ssj_werning_provider | cached Js (two-sided, h = 1e-3 b) passed without matching options: SSJ re-differentiates one-sided with h = 1e-4 (F1 type) |
| M13 | distribution | core.py:transition_matrix | Young lottery weights swapped between the lower and upper grid point |
| M14 | distribution | models.py:cohort_impc_paths | cohort distribution not divided by its mass |
| M15 | distribution | models.py:werning_primitives | Markov matrix transposed (column-stochastic convention) |
| M16 | distribution | core.py:_make_hh_function | borrowing constraint imposed at a_grid[1] instead of a_grid[0] (grid boundary) |
| M17 | distribution | models.py:AIYAGARI_DEFAULTS | asset grid upper bound 20 instead of 200 |
| M18 | equilibrium | models.py:fiscal_clearing | goods-market clearing sign error: Y - C + G |
| M19 | equilibrium | models.py:goods_clearing | TANK goods market cleared with savers' consumption only (Y = C_S) |
| M20 | equilibrium | models.py:fiscal | government budget constraint without interest: B - B(-1) - G + Tax |
| M21 | equilibrium | models.py:fiscal_hank | steady-state taxes omit debt service: Tax_ss = G instead of G + r B |
| M22 | equilibrium | models.py:ssj_tank_provider | saver iMPC matrix taken from the H-type block (copy-paste: J_H instead of J_S) |
| M23 | calibration | models.py:tank_model | zero-net-supply target for savers without the division by (1 - lam) |
| M24 | calibration | models.py:ssj_aiyagari_r | borrowing limit -1 instead of 0 |
| M25 | calibration | models.py:fiscal_hank | beta calibrated to A = 4 B (debt target in quarterly-GDP units) while B = 1 |
| C01 | control | models.py:AIYAGARI_DEFAULTS | asset grid refined from 1000 to 1500 points |
| C02 | control | core.py:_make_hh_function | Va = (1 + r) / c ** (1 / eis) instead of (1 + r) * c ** (-1 / eis) |
| C03 | control | models.py:ssj_werning_provider | GE target goods market (C = Y) instead of bond market (A = 0) |

Exact patches (file: old -> new; a newline inside a string is shown as \n):

- M01 core.py: `uc_nextgrid = beta * Va_p` -> `uc_nextgrid = beta * (1 + r) * Va_p`
- M01 core.py: `Va = (1 + r) * c ** (-1 / eis)\n        return Va, a, c` -> `Va = c ** (-1 / eis)\n        return Va, a, c`
- M02 models.py: `fiscal_res = B - (1 + r) * B(-1) - G + Tax` -> `fiscal_res = B - (1 + r(-1)) * B(-1) - G + Tax`
- M03 models.py: `Tax = Tax_ss + phi_B * (B(-1) - B_ss)` -> `Tax = Tax_ss + phi_B * (B - B_ss)`
- M04 models.py: `irf = m['model'].solve_impulse_linear(ss, ['Y'], ['asset_mkt'], {'r': dr}, Js={'hh': J},` -> `irf = m['model'].solve_impulse_linear(ss, ['Y'], ['asset_mkt'], {'r': np.r_[0.0, dr[:-1]]}, Js={'hh': J},`
- M05 models.py: `fiscal_res = B - (1 + r) * B(-1) - G + Tax` -> `fiscal_res = B - (1 + r) * B - G + Tax`
- M06 models.py: `irf = model.solve_impulse_linear(ss_chi, ['Y'], ['goods_mkt'], {'r': dr}, Js=Js,` -> `irf = model.solve_impulse_linear(ss_chi, ['Y'], ['goods_mkt'], {'r': dr / 100}, Js=Js,`
- M07 models.py: `AIYAGARI_DEFAULTS = dict(beta=0.96, alpha=0.36, delta=0.08, eis=0.5, n_e=7, n_a=1000,` -> `AIYAGARI_DEFAULTS = dict(beta=0.96, alpha=0.36, delta=0.02, eis=0.5, n_e=7, n_a=1000,`
- M08 core.py: `return sj.grids.markov_rouwenhorst(rho, sd, n)` -> `e, pi, Pi = sj.grids.markov_rouwenhorst(rho, sd, n)\n    return np.exp(np.log(e) - pi @ np.log(e)), pi, Pi`
- M09 models.py: `def household_jacobians(block, ss, T, inputs, outputs, h=1e-4):` -> `def household_jacobians(block, ss, T, inputs, outputs, h=1e-2):`
- M10 models.py: `J_S = m['hh_S'].jacobian(ss, ['Y_S', 'r'], ['C_S', 'A_S'], T=T, twosided=True)` -> `J_S = m['hh_S'].jacobian(ss, ['Y_S', 'r'], ['C_S', 'A_S'], T=T, twosided=False)`
- M11 models.py: `return {o: {i: np.asarray(J[o][i]) for i in inputs} for o in outputs}` -> `return {o: {i: np.asarray(J[o][i]).T for i in inputs} for o in outputs}`
- M12 models.py: `options=_het_opts(['hh'], h=1e-3 * b))` -> `options={})`
- M13 core.py: `lot = sp.csr_matrix((np.concatenate([pi.ravel(), 1 - pi.ravel()]),` -> `lot = sp.csr_matrix((np.concatenate([1 - pi.ravel(), pi.ravel()]),`
- M14 models.py: `ss_k.internals['hh']['Dbeg'] = Dk / mu` -> `ss_k.internals['hh']['Dbeg'] = Dk`
- M15 models.py: `e, pi, Pi = rouwenhorst(p['rho_e'], p['sd_e'], p['n_e'])\n    gam = 1 / p['eis']` -> `e, pi, Pi = rouwenhorst(p['rho_e'], p['sd_e'], p['n_e'])\n    Pi = Pi.T\n    gam = 1 / p['eis']`
- M16 core.py: `misc.setmin(a, a_grid[0])` -> `misc.setmin(a, a_grid[1])`
- M17 models.py: `a_max=200.0)` -> `a_max=20.0)`
- M18 models.py: `goods_mkt = Y - C - G` -> `goods_mkt = Y - C + G`
- M19 models.py: `@sj.simple\ndef goods_clearing(Y, C):\n    goods_mkt = Y - C` -> `@sj.simple\ndef goods_clearing(Y, C_S):\n    goods_mkt = Y - C_S`
- M20 models.py: `fiscal_res = B - (1 + r) * B(-1) - G + Tax` -> `fiscal_res = B - B(-1) - G + Tax`
- M21 models.py: `Tax_ss = p['G'] + p['r'] * p['B']` -> `Tax_ss = p['G']`
- M22 models.py: `return dict(M_H=np.asarray(J_H['C_H']['Y_H']), M_S=np.asarray(J_S['C_S']['Y_S']),` -> `return dict(M_H=np.asarray(J_H['C_H']['Y_H']), M_S=np.asarray(J_H['C_H']['Y_H']),`
- M23 models.py: `A_S_target = -lam * A_H / (1 - lam)` -> `A_S_target = -lam * A_H`
- M24 models.py: `a_grid = p['a_max'] * u ** 2` -> `a_grid = -1.0 + (p['a_max'] + 1.0) * u ** 2`
- M25 models.py: `f = lambda beta: hh.steady_state({**calib, 'beta': beta}, **SS_OPTS)['A'] - p['B']` -> `f = lambda beta: hh.steady_state({**calib, 'beta': beta}, **SS_OPTS)['A'] - 4 * p['B']`
- C01 models.py: `AIYAGARI_DEFAULTS = dict(beta=0.96, alpha=0.36, delta=0.08, eis=0.5, n_e=7, n_a=1000,` -> `AIYAGARI_DEFAULTS = dict(beta=0.96, alpha=0.36, delta=0.08, eis=0.5, n_e=7, n_a=1500,`
- C02 core.py: `Va = (1 + r) * c ** (-1 / eis)\n        return Va, a, c` -> `Va = (1 + r) / c ** (1 / eis)\n        return Va, a, c`
- C03 models.py: `irf = m['model'].solve_impulse_linear(ss, ['Y'], ['asset_mkt'], {'r': dr}, Js={'hh': J},` -> `irf = m['model'].solve_impulse_linear(ss, ['Y'], ['goods_mkt'], {'r': dr}, Js={'hh': J},`
````

### A.4 The agent's final report (verbatim)

The report text, with the JSON block it contained replaced by a pointer to A.5 (the JSON it returned is byte-identical to the file it wrote).

````markdown
I wrote 26 mutant candidates and 5 control candidates to `/tmp/claude-0/-home-user-hank-limits/dac78a62-da54-5c22-9fc7-31ee064ed8c5/scratchpad/mutgen/work/candidates.json`. Every patch applies on its own, each `old` string occurs exactly once, and every patched file byte-compiles. I did not run any models.

**Coverage:** timing 4, units 4, jacobian 4, distribution 4, equilibrium 5, calibration 5.

**Not duplicates of v0.1:** a script compared every patch with the 29 exact v0.1 patches and found no repeats. Two candidates touch a v0.1 line but make a different change, which your rules allow:
- N18 flips the sign on M23's line (M23 dropped the division).
- N19 changes the exponent on the `Va` line that M01 and C02 also change, plus the `c_nextgrid` line.

**Order:** most realistic first. The last few (N21–N26) are the ones I find least convincing.

**Behaviour, checked by reading the code and the SSJ source rather than running it:**
- **One candidate raises in one provider.** N19 (EIS/CRRA swap) should give wrong numbers everywhere else, but in the Werning provider it moves the zero-liquidity rate to about 12%, outside the root bracket. The root-find there should then fail with an error.
- **Some mutants only affect diagnostics:** N20 and N24 change reported diagnostic values only, and N08 changes only the reported `asset_mkt` path, not the GE solution.
- **Some mutants are only reached on certain code paths:**
  - N16 and N22 only matter when the fiscal provider uses the refined solver (`solver != 'ssj'`).
  - N10 only matters in the Werning provider, where the step is 1e-3·b rather than the default 1e-4.
  - N01 is identical to the original when rho = 0.

**Controls:** K02 should give bitwise-identical results. K01, K03, K04 and K05 should differ only at rounding level. K03 rounding differences could be amplified to about 1e-10 relative near the PIH limit, because of the slow backward-iteration stopping rule.

Files are in `/tmp/claude-0/-home-user-hank-limits/dac78a62-da54-5c22-9fc7-31ee064ed8c5/scratchpad/mutgen/work/`:
- `candidates.json`
- `make_candidates.py`
- `verify_candidates.py`
- `orig/` and `patched/<id>/` (scratch copies made by the verification script)

## candidates.json (in full)

(identical to `candidates.json`, reproduced below)

## Verification output (`python -I verify_candidates.py`)

```
N01  units        models.py:ssj_aiyagari_r                      patches=1  files=models.py OK
N02  distribution core.py:transition_matrix                     patches=1  files=core.py   OK
N03  jacobian     models.py:cohort_impc_paths                   patches=1  files=models.py OK
N04  timing       models.py:real_rate_rule                      patches=1  files=models.py OK
N05  calibration  models.py:werning_primitives                  patches=1  files=models.py OK
N06  distribution core.py:FastHetBlock.forward_steady_state     patches=1  files=core.py   OK
N07  units        models.py:FISCAL_DEFAULTS                     patches=1  files=models.py OK
N08  timing       models.py:fiscal_clearing                     patches=1  files=models.py OK
N09  distribution models.py:cohort_impc_paths                   patches=1  files=models.py OK
N10  jacobian     models.py:_het_opts                           patches=1  files=models.py OK
N11  equilibrium  models.py:tank_incidence                      patches=1  files=models.py OK
N12  timing       core.py:transition_matrix                     patches=1  files=core.py   OK
N13  calibration  models.py:pih_household                       patches=1  files=models.py OK
N14  equilibrium  models.py:tank_aggregate                      patches=1  files=models.py OK
N15  jacobian     models.py:household_jacobians                 patches=1  files=models.py OK
N16  equilibrium  models.py:solve_impulse_linear_refined        patches=1  files=models.py OK
N17  calibration  models.py:tank_model                          patches=1  files=models.py OK
N18  calibration  models.py:tank_model                          patches=1  files=models.py OK
N19  units        core.py:_make_hh_function                     patches=2  files=core.py   OK
N20  distribution core.py:mass_at_constraint                    patches=1  files=core.py   OK
N21  jacobian     models.py:ssj_werning_provider                patches=1  files=models.py OK
N22  equilibrium  models.py:solve_impulse_linear_refined        patches=1  files=models.py OK
N23  timing       core.py:_make_hh_function                     patches=1  files=core.py   OK
N24  units        models.py:werning_model                       patches=1  files=models.py OK
N25  equilibrium  models.py:ssj_aiyagari_r                      patches=2  files=models.py OK
N26  calibration  models.py:pih_household                       patches=1  files=models.py OK
K01  control      models.py:ssj_aiyagari_r                      patches=1  files=models.py OK
K02  control      core.py:transition_matrix                     patches=1  files=core.py   OK
K03  control      core.py:_make_hh_function                     patches=1  files=core.py   OK
K04  control      models.py:solve_impulse_linear_refined        patches=1  files=models.py OK
K05  control      models.py:tank_incidence                      patches=1  files=models.py OK

mutants: 26 controls: 5
mutants per category: {'calibration': 5, 'distribution': 4, 'equilibrium': 5, 'jacobian': 4, 'timing': 4, 'units': 4}
every category >= 3: True | ids N01..N26, K01..K05: True
no identical patch sets: True
ALL CHECKS PASSED
```

v0.1 duplicate check (exact old/new comparison against the 29 parsed v0.1 patches):
```
v0.1 patches parsed: 29
same line as M23 (different change): N18
same line as M01 (different change): N19
same line as C02 (different change): N19
exact duplicates of v0.1: 0
```
````

### A.5 `candidates.json` (verbatim)

```json
{
 "mutants": [
  {
   "id": "N01",
   "category": "units",
   "target": "models.py:ssj_aiyagari_r",
   "change": "The log-income sd passed to the provider is multiplied by sqrt(1 - rho^2) before discretisation, so the chain has unconditional sd sd*sqrt(1 - rho^2) instead of sd.",
   "why_plausible": "Aiyagari (1994) writes the innovation as sigma*(1 - rho^2)^(1/2), and a programmer who believes the Rouwenhorst routine takes the innovation sd passes that product instead of the unconditional sd.",
   "patches": [
    {
     "file": "models.py",
     "old": "        e, pi, Pi = rouwenhorst(rho, sd, p['n_e'])",
     "new": "        e, pi, Pi = rouwenhorst(rho, sd * np.sqrt(1 - rho ** 2), p['n_e'])"
    }
   ]
  },
  {
   "id": "N02",
   "category": "distribution",
   "target": "core.py:transition_matrix",
   "change": "The exogenous step of the Young-lottery chain uses kron(Pi, I) instead of kron(Pi', I), so mass moves with the rows of Pi and the stationary income marginal becomes uniform instead of the Markov chain's stationary distribution.",
   "why_plausible": "Row- versus column-stochastic conventions are easy to mix up when building a sparse forward operator by hand, and kron(Pi, I) still yields a normalisable distribution without any error.",
   "patches": [
    {
     "file": "core.py",
     "old": "    X = sp.kron(sp.csr_matrix(Pi.T), sp.identity(n_a), format='csr')",
     "new": "    X = sp.kron(sp.csr_matrix(Pi), sp.identity(n_a), format='csr')"
    }
   ]
  },
  {
   "id": "N03",
   "category": "jacobian",
   "target": "models.py:cohort_impc_paths",
   "change": "The central finite difference of the cohort consumption paths is divided by h instead of 2h, doubling every cohort iMPC path.",
   "why_plausible": "Forgetting the factor 2 when turning a one-sided difference into a two-sided one is a classic slip.",
   "patches": [
    {
     "file": "models.py",
     "old": "        paths.append((np.asarray(up) - np.asarray(dn)) / (2 * h))",
     "new": "        paths.append((np.asarray(up) - np.asarray(dn)) / h)"
    }
   ]
  },
  {
   "id": "N04",
   "category": "timing",
   "target": "models.py:real_rate_rule",
   "change": "The real-rate rule responds to last period's output gap Y_{t-1} instead of the current gap Y_t.",
   "why_plausible": "Since r_t is the ex-post return on assets bought at t-1, a programmer can reason that the rate was set at t-1 and so should respond to date t-1 output.",
   "patches": [
    {
     "file": "models.py",
     "old": "    r = rss + phi_y * (Y - 1) + rshock",
     "new": "    r = rss + phi_y * (Y(-1) - 1) + rshock"
    }
   ]
  },
  {
   "id": "N05",
   "category": "calibration",
   "target": "models.py:werning_primitives",
   "change": "The countercyclical incidence is normalised so that sum_i pi_i zeta_i = 0 instead of sum_i pi_i e_i zeta_i = 0, so aggregate household income no longer equals Y off the steady state.",
   "why_plausible": "Making the incidence 'average to zero' with the stationary probabilities alone, and forgetting the income weights e_i, is a natural shortcut.",
   "patches": [
    {
     "file": "models.py",
     "old": "        zeta[1:] = -pi[0] * e[0] * zeta[0] / np.sum(pi[1:] * e[1:])",
     "new": "        zeta[1:] = -pi[0] * zeta[0] / np.sum(pi[1:])"
    }
   ]
  },
  {
   "id": "N06",
   "category": "distribution",
   "target": "core.py:FastHetBlock.forward_steady_state",
   "change": "The post-shock distribution D is set equal to the beginning-of-period distribution Dbeg, so steady-state aggregates weight policies indexed by e_t with the e_{t-1} distribution.",
   "why_plausible": "Someone who thinks 'D and Dbeg coincide in a stationary state' would skip the exogenous step when returning the pair.",
   "patches": [
    {
     "file": "core.py",
     "old": "        D = exog.forward(Dbeg)\n",
     "new": "        D = Dbeg\n"
    }
   ]
  },
  {
   "id": "N07",
   "category": "units",
   "target": "models.py:FISCAL_DEFAULTS",
   "change": "The fiscal-HANK calibration uses the quarterly real rate r = 0.005 inside an otherwise annual calibration (B = 1 annual GDP, annual persistence).",
   "why_plausible": "r = 0.005 is the standard quarterly value in the SSJ example models and is easy to copy over unchanged.",
   "patches": [
    {
     "file": "models.py",
     "old": "FISCAL_DEFAULTS = dict(r=0.02, ",
     "new": "FISCAL_DEFAULTS = dict(r=0.005, "
    }
   ]
  },
  {
   "id": "N08",
   "category": "timing",
   "target": "models.py:fiscal_clearing",
   "change": "Asset-market clearing compares households' end-of-period assets A_t with last period's debt B_{t-1} instead of B_t.",
   "why_plausible": "SSJ's aggregate A is the end-of-period choice but is often read as wealth brought into the period, which suggests matching it to B(-1).",
   "patches": [
    {
     "file": "models.py",
     "old": "    asset_mkt = A - B\n",
     "new": "    asset_mkt = A - B(-1)\n"
    }
   ]
  },
  {
   "id": "N09",
   "category": "distribution",
   "target": "models.py:cohort_impc_paths",
   "change": "Wealth cohorts are built from the post-shock distribution D instead of the beginning-of-period distribution Dbeg, so the cohort weights and initial conditions are wrong and the cohort paths no longer sum to M x.",
   "why_plausible": "Both arrays are 'the distribution' in ss.internals, and 'D' looks like the canonical one.",
   "patches": [
    {
     "file": "models.py",
     "old": "    Dbeg = ss.internals['hh']['Dbeg']",
     "new": "    Dbeg = ss.internals['hh']['D']"
    }
   ]
  },
  {
   "id": "N10",
   "category": "jacobian",
   "target": "models.py:_het_opts",
   "change": "The options passed to SSJ's impulse_linear drop the step size, so the HetBlock is re-differentiated with SSJ's default h = 1e-4 while the supplied Js use the caller's h (1e-3 b in the Werning provider).",
   "why_plausible": "Forgetting to forward one keyword argument in a small helper is an everyday mistake, and nothing changes visibly when h equals the default.",
   "patches": [
    {
     "file": "models.py",
     "old": "    return {n: {'twosided': True, 'h': h} for n in names}",
     "new": "    return {n: {'twosided': True} for n in names}"
    }
   ]
  },
  {
   "id": "N11",
   "category": "equilibrium",
   "target": "models.py:tank_incidence",
   "change": "H-type income is Y (Y/Yss)^chi, so chi acts as an excess elasticity (total elasticity 1 + chi) instead of the elasticity of Y_H with respect to Y.",
   "why_plausible": "The same file writes cyclical income as Y e (Y/Yss)^zeta with zeta the excess elasticity, and copying that form into the TANK incidence is easy.",
   "patches": [
    {
     "file": "models.py",
     "old": "    Y_H = Yss * (Y / Yss) ** chi",
     "new": "    Y_H = Y * (Y / Yss) ** chi"
    }
   ]
  },
  {
   "id": "N12",
   "category": "timing",
   "target": "core.py:transition_matrix",
   "change": "The one-period operator is composed as X @ lot (savings lottery first, then income shock) instead of lot @ X, so the vector it returns as Dbeg is really the post-shock distribution.",
   "why_plausible": "Reading 'shock, then lottery' left to right and writing the matrix product in that order is a common composition-order slip.",
   "patches": [
    {
     "file": "core.py",
     "old": "    return (lot @ X).tocsc()",
     "new": "    return (X @ lot).tocsc()"
    }
   ]
  },
  {
   "id": "N13",
   "category": "calibration",
   "target": "models.py:pih_household",
   "change": "The final discount factor is set to 1/(1+r) - kappa instead of (1 - kappa)/(1+r), which does not match the parameterisation used in the calibration root-find, so mean assets miss A_target and kappa no longer equals 1 - beta (1+r).",
   "why_plausible": "Writing the patience gap additively in one place and multiplicatively in the other is a typical inconsistency between a residual function and the final assignment.",
   "patches": [
    {
     "file": "models.py",
     "old": "    calib['beta'] = beta_ra * (1 - kappa)",
     "new": "    calib['beta'] = beta_ra - kappa"
    }
   ]
  },
  {
   "id": "N14",
   "category": "equilibrium",
   "target": "models.py:tank_aggregate",
   "change": "Aggregate consumption weights the H type by 1 - lam and savers by lam, the reverse of the population shares.",
   "why_plausible": "Whether lambda denotes the hand-to-mouth or the saver share differs across TANK papers, so the weights are easily swapped in one block.",
   "patches": [
    {
     "file": "models.py",
     "old": "    C = lam * C_H + (1 - lam) * C_S",
     "new": "    C = (1 - lam) * C_H + lam * C_S"
    }
   ]
  },
  {
   "id": "N15",
   "category": "jacobian",
   "target": "models.py:household_jacobians",
   "change": "The household Jacobians are computed with SSJ's default one-sided differences instead of two-sided ones, while the GE options still re-differentiate two-sided.",
   "why_plausible": "Dropping a keyword to halve the cost or tidy the call looks harmless.",
   "patches": [
    {
     "file": "models.py",
     "old": "    J = block.jacobian(ss, inputs, outputs, T=T, twosided=True, h=h)",
     "new": "    J = block.jacobian(ss, inputs, outputs, T=T, h=h)"
    }
   ]
  },
  {
   "id": "N16",
   "category": "equilibrium",
   "target": "models.py:solve_impulse_linear_refined",
   "change": "After solving for the unknowns, outputs are propagated from dU alone, so the direct effects of the exogenous shocks (G, r, rshock) are missing from all reported paths.",
   "why_plausible": "Thinking of the final step as 'propagate the solved unknowns' and leaving out the shocks is an easy omission.",
   "patches": [
    {
     "file": "models.py",
     "old": "    return (inputs | dU) | model.impulse_linear(ss, dU | inputs, None, Js, options)",
     "new": "    return (inputs | dU) | model.impulse_linear(ss, dU, None, Js, options)"
    }
   ]
  },
  {
   "id": "N17",
   "category": "calibration",
   "target": "models.py:tank_model",
   "change": "The saver block in the TANK model is given the H type's asset grid (limit 0, top 50) instead of its own near-natural-limit grid, while beta_S was calibrated on the latter.",
   "why_plausible": "This is a copy-paste slip inside a long calibration dict where _H and _S entries sit side by side.",
   "patches": [
    {
     "file": "models.py",
     "old": "                 a_grid_S=S['a_grid'])",
     "new": "                 a_grid_S=a_grid_H)"
    }
   ]
  },
  {
   "id": "N18",
   "category": "calibration",
   "target": "models.py:tank_model",
   "change": "The savers' asset target has the wrong sign (+lam A_H/(1-lam)), so both types hold positive bonds and the bond is not in zero net supply.",
   "why_plausible": "Dropping a minus sign when solving lam A_H + (1-lam) A_S = 0 for A_S is a common algebra slip.",
   "patches": [
    {
     "file": "models.py",
     "old": "    A_S_target = -lam * A_H / (1 - lam)",
     "new": "    A_S_target = lam * A_H / (1 - lam)"
    }
   ]
  },
  {
   "id": "N19",
   "category": "units",
   "target": "core.py:_make_hh_function",
   "change": "The parameter eis is used as the CRRA coefficient in both the Euler inversion and the envelope condition, so the household's actual EIS is 1/eis.",
   "why_plausible": "Many codes take sigma (risk aversion) rather than the EIS, and a consistent reciprocal swap in both lines still converges without warning.",
   "patches": [
    {
     "file": "core.py",
     "old": "        c_nextgrid = uc_nextgrid ** (-eis)",
     "new": "        c_nextgrid = uc_nextgrid ** (-1 / eis)"
    },
    {
     "file": "core.py",
     "old": "        Va = (1 + r) * c ** (-1 / eis)\n        return Va, a, c",
     "new": "        Va = (1 + r) * c ** (-eis)\n        return Va, a, c"
    }
   ]
  },
  {
   "id": "N20",
   "category": "distribution",
   "target": "core.py:mass_at_constraint",
   "change": "The constrained share weights the policy mask (indexed by e_t) with the beginning-of-period distribution Dbeg instead of D.",
   "why_plausible": "Dbeg is the distribution most codes call 'the' stationary distribution, so it is easy to reach for it here.",
   "patches": [
    {
     "file": "core.py",
     "old": "    return float(np.sum(inter['D'][inter['a'] <= a_grid[0] + tol]))",
     "new": "    return float(np.sum(inter['Dbeg'][inter['a'] <= a_grid[0] + tol]))"
    }
   ]
  },
  {
   "id": "N21",
   "category": "jacobian",
   "target": "models.py:ssj_werning_provider",
   "change": "The Werning household Jacobian uses a fixed step h = 1e-3 instead of one scaled with the liquidity b, so the step is of the order of the whole wealth distribution for small b and no longer matches the b-scaled step in the GE options.",
   "why_plausible": "Hard-coding the step in one call while the neighbouring line keeps the b-scaling is an easy oversight.",
   "patches": [
    {
     "file": "models.py",
     "old": "    J = m['block'].jacobian(ss, ['Y', 'r'], ['C', 'A'], T=T, twosided=True, h=1e-3 * b)",
     "new": "    J = m['block'].jacobian(ss, ['Y', 'r'], ['C', 'A'], T=T, twosided=True, h=1e-3)"
    }
   ]
  },
  {
   "id": "N22",
   "category": "equilibrium",
   "target": "models.py:solve_impulse_linear_refined",
   "change": "The iterative-refinement residual is written as dH - H_U x instead of -dH - H_U x, so the refined solver converges to the solution of H_U dU = +dH (sign-flipped unknowns).",
   "why_plausible": "Writing the residual as 'dH minus H_U x' and forgetting that the right-hand side is -dH is a typical sign slip.",
   "patches": [
    {
     "file": "models.py",
     "old": "        x = x + sla.lu_solve(lu, -dH - H_U @ x)",
     "new": "        x = x + sla.lu_solve(lu, dH - H_U @ x)"
    }
   ]
  },
  {
   "id": "N23",
   "category": "timing",
   "target": "core.py:_make_hh_function",
   "change": "Cash on hand is (1+r)(a + y), so current income earns interest within the period, instead of (1+r) a + y.",
   "why_plausible": "Some textbooks use the timing a' = (1+r)(a + y - c), and half-adopting it in the cash-on-hand line is a plausible mix-up.",
   "patches": [
    {
     "file": "core.py",
     "old": "        c_nextgrid = uc_nextgrid ** (-eis)\n        coh = (1 + r) * a_grid[np.newaxis, :] + y[:, np.newaxis]",
     "new": "        c_nextgrid = uc_nextgrid ** (-eis)\n        coh = (1 + r) * (a_grid[np.newaxis, :] + y[:, np.newaxis])"
    }
   ]
  },
  {
   "id": "N24",
   "category": "units",
   "target": "models.py:werning_model",
   "change": "max_wealth_over_b measures the top occupied grid point from zero instead of from the borrowing limit -b, so it is understated by exactly 1 and is inconsistent with grid_top_over_b.",
   "why_plausible": "Normalising wealth by b and forgetting the shift by the limit is an easy slip when the grid starts at -b.",
   "patches": [
    {
     "file": "models.py",
     "old": "                max_wealth_over_b=float((a_grid[used] + b) / b),",
     "new": "                max_wealth_over_b=float(a_grid[used] / b),"
    }
   ]
  },
  {
   "id": "N25",
   "category": "equilibrium",
   "target": "models.py:ssj_aiyagari_r",
   "change": "The capital stock is set to the capital-output ratio alpha/(r+delta) instead of (alpha/(r+delta))^(1/(1-alpha)), so the capital-market condition A = K and the wage use a K that is inconsistent with the firm's first-order condition.",
   "why_plausible": "K/Y = alpha/(r+delta) is the familiar steady-state formula, and treating it as K (implicitly Y = 1) is a natural shortcut.",
   "patches": [
    {
     "file": "models.py",
     "old": "    def resid(r):\n        K = (al / (r + de)) ** (1 / (1 - al))",
     "new": "    def resid(r):\n        K = al / (r + de)"
    },
    {
     "file": "models.py",
     "old": "    gap, ss = resid(r)\n    K = (al / (r + de)) ** (1 / (1 - al))",
     "new": "    gap, ss = resid(r)\n    K = al / (r + de)"
    }
   ]
  },
  {
   "id": "N26",
   "category": "calibration",
   "target": "models.py:pih_household",
   "change": "The near-natural borrowing limit is computed as y_min/(1+r) instead of y_min/r, which is far tighter, so the saver is no longer effectively unconstrained.",
   "why_plausible": "Confusing the perpetual-income bound y_min/r with a one-period present value y_min/(1+r) is a plausible slip.",
   "patches": [
    {
     "file": "models.py",
     "old": "    a_min = -p['nbl_frac'] * y_min / r",
     "new": "    a_min = -p['nbl_frac'] * y_min / (1 + r)"
    }
   ]
  }
 ],
 "controls": [
  {
   "id": "K01",
   "category": "control",
   "target": "models.py:ssj_aiyagari_r",
   "change": "The wage is computed from zero profits, w = K^alpha - (r + delta) K, instead of the marginal product (1 - alpha) K^alpha.",
   "why_equivalent": "Under Cobb-Douglas with L = 1 and r + delta = alpha K^(alpha-1), the two expressions are identical, so only rounding differs.",
   "patches": [
    {
     "file": "models.py",
     "old": "        w = (1 - al) * K ** al",
     "new": "        w = K ** al - (r + de) * K"
    }
   ]
  },
  {
   "id": "K02",
   "category": "control",
   "target": "core.py:transition_matrix",
   "change": "The redundant np.clip of the policy before interpolate_coord_robust is removed, so the out-of-grid handling relies on the existing clip of the lottery weights.",
   "why_equivalent": "interpolate_coord_robust extrapolates with weight < 0 (above) or > 1 (below), and clipping those weights to [0, 1] reproduces exactly the indices and weights obtained from the clipped policy.",
   "patches": [
    {
     "file": "core.py",
     "old": "    i, pi = interpolate_coord_robust(grid, np.clip(a_pol, grid[0], grid[-1]))",
     "new": "    i, pi = interpolate_coord_robust(grid, a_pol)"
    }
   ]
  },
  {
   "id": "K03",
   "category": "control",
   "target": "core.py:_make_hh_function",
   "change": "The endogenous-grid step interpolates on beginning-of-period assets (c + a' - y)/(1+r) evaluated at a_grid instead of on cash on hand c + a' evaluated at coh.",
   "why_equivalent": "Both the data and query points are mapped by the same increasing affine map within each income state, and linear interpolation and extrapolation are invariant to that map.",
   "patches": [
    {
     "file": "core.py",
     "old": "        a = interpolate.interpolate_y(c_nextgrid + a_grid, coh, a_grid)",
     "new": "        a = interpolate.interpolate_y((c_nextgrid + a_grid - y[:, np.newaxis]) / (1 + r), a_grid, a_grid)"
    }
   ]
  },
  {
   "id": "K04",
   "category": "control",
   "target": "models.py:solve_impulse_linear_refined",
   "change": "One more step of iterative refinement (n_refine 3 instead of 2).",
   "why_equivalent": "Iterative refinement converges to the same solution of H_U dU = -dH, and an extra step only removes residual rounding error.",
   "patches": [
    {
     "file": "models.py",
     "old": "                                 n_refine=2):",
     "new": "                                 n_refine=3):"
    }
   ]
  },
  {
   "id": "K05",
   "category": "control",
   "target": "models.py:tank_incidence",
   "change": "Saver income is written as Y minus the per-saver transfer implied by H-type incidence, Y_S = Y - lam (Y_H - Y)/(1 - lam), instead of (Y - lam Y_H)/(1 - lam).",
   "why_equivalent": "The two expressions are algebraically identical, so the steady state (Y_S = 1) and the derivative (1 - lam chi)/(1 - lam) are unchanged.",
   "patches": [
    {
     "file": "models.py",
     "old": "    Y_S = (Y - lam * Y_H) / (1 - lam)",
     "new": "    Y_S = Y - lam * (Y_H - Y) / (1 - lam)"
    }
   ]
  }
 ]
}
```
