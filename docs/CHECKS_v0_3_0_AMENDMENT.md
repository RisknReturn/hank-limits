# v0.3.0 release candidate: review-driven changes to the v0.3 checks, and their re-run

Status: **design and predictions** committed as A (`471fc47`) on branch `v0_3_0-rc`, before any code of this amendment and before any re-run (sections 1–9 as committed then). The thresholds of the new quantities (section 5) were filled in by commit B (`a290c75`, code, tests and thresholds), before any mutant was run. The re-run outputs are commit C (`ad1857d`), the results (section 10) commit D. Release metadata (version 0.3.0) is commit E.

Starting point: `main` at `281de19` (v0.3, unreleased; the last release is v0.2.1, concept DOI 10.5281/zenodo.23489176). Pre-registration of v0.3: `docs/CHECKS_PLAN_v0_3.md` ("the plan"); its results: `docs/MUTATION_v0_2.md`.

## 1. Why

A pre-release review of v0.3 (GPT, 2026-10-22) found no P0 or P1 defect. It reproduced **74 passed, 3 xfailed** in a separate Windows environment, confirmed that the 88 thresholds follow the rule of the plan (section 4), and that the pre-registration order and the reported numbers (14/20, 19/20, 25/25) agree with the record. It raised four P2 items, fixed here before the release, and P3 documentation items (section 8):

- **R1. Missing, non-finite and skipped values can hide a verdict** (`tests/test_v0_3_checks.py`).
  - `_assert_below({'hh_budget': [0.0, NaN]}, ...)` passes: `max` keeps the leading 0.0. With the order reversed it fails.
  - `[0.0, None]` passes: None is dropped.
  - If one quantity of a test has no value, the whole test is skipped, even when another quantity of the same test is above its threshold. Example: V1 direct with the proportional-income value absent and 10.0 for the countercyclical one gives a skip.
  - `_assert_paths` uses the same `max`, and `_paths` collects only the paths the provider returned, so a missing part of a path set is not visible.
- **R2. The optional-field convention of `hank_limits/interface.py` and the code disagree.**
  - The T5 path checks read `dB` and `dC` without checking that they exist (`checks.py`, lines 252 and 257). A provider without them gets a `KeyError` under `v0_3` → `errors`, which fails every V test on T5.
  - T3: `interface.py` puts `beta` in `ss_detail`; the code reads `out['beta']`.
- **R3. The checks do not compare the shock the provider solved with the shock the caller requested.** A linear provider that applies twice the requested `G`, `r` and `rshock`, and returns paths consistent with that, passes the IKC test and every path check (largest residual 6.07e-16).
- **R4. V1 assumes that T3 is solved on the asset-market target.** An external code solved on the goods-market target, equivalent in exact arithmetic, alarms (C03).

The V6 pin of F3 is also split in two (section 3.5).

## 2. Rules

1. **Scope.** Branch `v0_3_0-rc` from `main` `281de19`. A pull request to `main` is opened and not merged here. No tag and no GitHub Release: the user creates them, and Zenodo then mints the version DOI.
2. **Historical material is not changed:** `docs/CHECKS_PLAN_v0_3.md`, `docs/MUTATION_PLAN_v0_1.md`, `docs/MUTATION_v0_1.md`, `scripts/mutation_v0_1.py`, `scripts/mutation_v0_2.py`, `out/mutation_v0_1.*`, `out/mutation_v0_2.*`. Exceptions: two dated corrections in `docs/MUTATION_v0_2.md` in the form of its existing "Correction, 2026-10-22" paragraph (section 8), and one line in section 10 (Amendments) of the plan that points to this file.
3. **Model side** (`models.py`, `core.py`). No computation is changed. Allowed: new fields in the dicts the providers return, optional arguments whose default is the current behaviour, docstrings. The 276 numbers of the v0.1 check summaries must equal run B0 of `out/mutation_v0_1.json` bit for bit in the re-run baselines.
4. **Thresholds.** The 88 existing thresholds are not changed. The new quantities get thresholds by the rule of the plan, section 4, `roundup_1sig(max(100 × b_max, 10 × τ))`, with b_max from two unmodified runs of the code of commit B. They are committed in B, before any mutant is run, in a table of the plan's form (section 5).
5. **Checker principles, as in v0.3.** `checks.py` does not depend on SSJ. A new check computes its quantity directly from reported primitive quantities and from **the inputs the caller passed**. Documented values are read from `spec.py`.
6. **Order.** A (this file) → B (code, tests, thresholds) → C (re-run outputs) → D (results) → E (release metadata). Each commit is pushed on its own.

## 3. Design

### 3.1 R1: a verdict that missing or non-finite values cannot hide

**Coverage is declared by the check, not inferred from what the provider returned.** In each check function, at each setting, the quantities that are expected there are declared first, with a placeholder `None` that is replaced when the quantity is computed.

| Model | Settings | Quantities expected at each setting |
|---|---|---|
| T1 | each eps | V2, V3, V4(a) of T1 |
| T2 | each eps | per income type H and S, and the aggregates |
| T3 | each b | per incidence (acyclical, countercyclical); V1 by the reported GE target (section 3.4) |
| T4 | each eps | V2, V3, V4(a), V4(b) Euler |
| T5 | each T | steady-state quantities; path quantities **for every shock label the caller passed** (rate rule only in rule-based experiments, `exog_r` only in fixed-rate experiments) |

Records in `res['v0_3']` (new keys; the value lists keep their v0.3 form, one entry per setting):

- `settings`: the setting labels in order (e.g. `eps=0.2`, `b=0.001`, `T=200`);
- `at[q]`: the setting of each value in `res['v0_3'][q]` (paths: `paths[label]['at'][q]`);
- `status[q]`: coverage per quantity: number expected, number reported, the settings with a missing value, those with a non-finite value, and a state `complete`, `incomplete`, `not_reported` or `nonfinite`. Thresholds are applied by the tests, not here.

**Verdict of one V test** (`checks.v03_verdict(v, quantities, thresholds)`, numpy only; every V test uses it through `_assert_below` and `_assert_paths`):

1. **failed** if a v0.3 computation raised an exception (as in v0.3), or if any recorded value of any of the test's quantities is non-finite (NaN, ±inf), negative, or not below its threshold. This holds whatever else is missing. Every V quantity is an absolute value or a maximum of absolute values, so a negative number can only be a malformed report (V1(a) is a value the provider reports);
2. otherwise **skipped** ("not reported", with the reasons recorded under `skipped`) if none of the test's quantities has a value at any setting: the provider does not support the check, as in v0.3;
3. otherwise **incomplete** if any expected value is missing: pytest skip with a message that starts with "incomplete" and lists the missing settings. It is never a pass;
4. otherwise **passed**.

A quantity that the check declared at no setting does not apply and is ignored, e.g. the goods-market V1 quantities when every T3 solve is on the asset market. For a hand-made `v0_3` dict without coverage records (as in the regression tests), every listed quantity is expected and each list entry is one setting.

**Regression tests** use the real helpers `_assert_below` and `_assert_paths` of `tests/test_v0_3_checks.py` and `checks.v03_verdict`:

- `{'hh_budget': [0.0, NaN]}` → failed (also with the order reversed); `[0.0, inf]` → failed;
- `{'hh_budget': [0.0, None]}` → incomplete;
- V1 direct with the proportional-income quantity absent and `[10.0]` for the countercyclical one → failed;
- T5 paths: `[0.0, NaN]` → failed; a declared path value that is missing, and a shock label the provider did not return (through `check_ikc` and a fake provider) → incomplete.

They are in a new file `tests/test_v0_3_0_checker.py`. It uses no model code, so no mutant can change it; the mutation runner maps it to no group.

### 3.2 R2: one field convention, in the documentation and in the code

- **T5 paths: fields per equation.** Fiscal rule (dB, dTax); government budget (dB, dTax, dr, dG); goods market (dY, dC, dG); rate rule (dr, dY, and the caller's rshock, section 3.3). A missing field skips only that equation, recorded under `skipped` with the path label and the missing field. It is no longer an error that fails every V test.
- **Steady-state quantities: fields per equation as well.** A missing income vector y no longer skips the market-clearing residuals that do not use it, and so on.
- **T3 `beta`** is read from `ss_detail['beta']`, as documented. The SSJ provider already reports it there; the top-level `beta` it also returns stays, so existing callers do not break.
- **`interface.py` lists every field that the checks read**, v0.1 part and V part, as required or optional, per provider, with the V quantities each optional field feeds. The code reads no other field (enforced by the fake providers below).
- **Two fake providers per check** in the regression tests:
  - one that returns exactly the documented required and optional fields: every V test reaches a verdict (passed or failed) and no error is recorded;
  - one that returns only the required fields: no error is recorded. Every V test that needs an optional field is skipped or incomplete, and `skipped` lists the reasons. The checks that need only required fields and the caller's inputs (the T5 rate rule, `exog_G`, `exog_r`; section 3.3) reach a verdict.

### 3.3 R3: the solved shock against the requested shock

- **New path quantities** in `_v03_fiscal`, normalised by s (the largest absolute requested input) like the other path checks:
  - `exog_G` = max_t |dG_t − G_t| / s, with G the path the caller requested, or zero when the experiment has no G shock;
  - `exog_r` = max_t |dr_t − r_t| / s in the fixed-rate experiments (requested r, or zero). In rule-based experiments dr is endogenous and the rate rule checks it;
  - `exog_rshock` = max_t |drshock_t − rshock_t| / s (requested rshock, or zero; `drshock` stays an optional field).
- **The rate-rule residual uses the requested rshock**, not the reported drshock: max_t |dr_t − φ_y dY_t − rshock_t| / s.
- **The fiscal rule, the government budget and the goods market keep using the reported dG and dr.** The exog quantities compare those with the request. So a wrong exogenous path is attributed to `exog_*`, and a violated equation to that equation; together they are equivalent to evaluating the equations at the requested paths.
- **Tests** (group V): `test_v4_exogenous_G_as_requested`, `test_v4_exogenous_r_as_requested`, `test_v4_exogenous_rshock_as_requested`.
- **Regression test:** a small linear fake fiscal provider that applies twice the requested G, r and rshock, and returns paths consistent with that. It passes the IKC check (v0.1 part) and the fiscal-rule, budget and goods-market path checks, and fails the three exog tests and the rate rule.
- **Threshold** by the rule of section 2. The SSJ provider takes the exogenous paths from the inputs (`get` in `ssj_fiscal_provider`), so b_max is expected to be exactly 0, and τ = 1e-14 ("exact up to rounding", plan section 4) gives 1e-13.

### 3.4 R4: the GE target of T3 in the provider convention

**Provider.**

- `WerningProvider` reports `ge_target` ∈ {`'asset_mkt'`, `'goods_mkt'`}: the market whose clearing the solve imposed.
  - Absent: `'asset_mkt'` is assumed and recorded as assumed (backward compatible).
  - Any other value: V1 is skipped at that setting as "applicability condition not met". The placeholder stays missing, so a mix with settings whose target is known is incomplete, not passed.
- New optional field `J_C` = {'Y': J^{C,Y}, 'r': J^{C,r}} (T × T) for the goods-market residual. `J_A` stays as it is.
- SSJ reference provider: optional argument `ge_target='asset_mkt'` (the default is the current behaviour). `targets = [ge_target]` is passed to the solve, and **the same variable** gives the reported `ge_target` and `diag['ge_target_resid']`, the residual of the solved target. A change of the solve's target, like C03, therefore shows in what is reported.

**Checker.** V1 computes the residual of the reported target directly and divides it by b. dr is the rate path the check passed, as before.

| Quantity | Target | Definition | Threshold |
|---|---|---|---|
| `ge_target_resid_direct_{acyc,cyc}` | asset market | max \|J_A['Y'] dY + J_A['r'] dr\| / b | 2e-12 (unchanged) |
| `ge_target_resid_direct_goods_{acyc,cyc}` | goods market | max \|dY − J_C['Y'] dY − J_C['r'] dr\| / b | new (section 5) |
| `ge_target_resid_reported_{acyc,cyc}` | asset market | V1(a): `diag['ge_target_resid']` / b, must be finite and non-negative | 2e-12 (unchanged) |
| `ge_target_resid_reported_goods_{acyc,cyc}` | goods market | V1(a) as above | new (section 5) |

- **The other market's residual** is recorded without a threshold, where its Jacobians are reported: `v0_3['diagnostics']`, `goods_mkt_resid_direct_*` for an asset-market solve, `asset_mkt_resid_direct_*` for a goods-market solve. The target of each setting is recorded in `v0_3['ge_target']`.
- **Thresholds of the goods-market quantities:** two unmodified runs of the reference provider with `ge_target='goods_mkt'`, both incidences and all seven b, by the rule of section 2, τ = 1e-14 (linear solve).
- **Documentation:** `interface.py` and the V1 rows of both READMEs: the target convention, and that V1 applies to a T3 solve on one of the two markets.

### 3.5 V6 split

`test_v6_f3_hetinput_hetoutput_two_sided` becomes two strict xfails, each with SSJ code only:

- `test_v6a_f3_hetinput_two_sided`: the hetinput derivative that `HetBlock` uses equals the central difference (fails: always one-sided);
- `test_v6b_f3_hetoutput_two_sided`: the hetoutput Jacobian with `twosided=True` equals the one obtained when the stored default is honoured (fails: one-sided even when two-sided differences are requested).

If only one part is fixed upstream, only that test turns XPASS.

### 3.6 Model side, in full, and patch anchors

Only `ssj_werning_provider` changes:

- signature (`ge_target='asset_mkt'`);
- a new line `targets = [ge_target]`, and the solve line passes `targets`;
- new return fields `ge_target` and `J_C`;
- `diag['ge_target_resid']` reads `irf[targets[0]]`;
- the docstring: it said "(goods market)" while the target is the asset market (P3).

No other line of `models.py` or `core.py` changes.

Patch anchors (checked against the plan's patch lists):

- **M04 and C03** (v0.1) anchor on the T3 solve line, which now passes `targets`. They are re-anchored, the same change at the new location:
  - M04: the same one-period shift of the rate path, on the new solve line;
  - C03: `targets = [ge_target]` → `targets = ['goods_mkt']`. The target the solve imposes changes, as in the original patch; by design (R4) the reported target follows it.
- **M12** anchors on the `options=` line, which does not change.
- **No v0.2 patch** anchors on a changed line.

`scripts/mutation_v0_2.py` has an empty `REANCHOR` and is not changed (rule 2). The re-run therefore uses a copy, `scripts/mutation_v0_3_0.py`, with only these changes: the two `REANCHOR` entries, its own module name as the pytest plugin and in the environment record, the default output path, and the docstring.

## 4. New quantities and tests

| Test (group V) | Quantity | Model | Fields used |
|---|---|---|---|
| `test_v1_werning_ge_target_resid_reported` | `ge_target_resid_reported_goods_{acyc,cyc}` (new), besides the asset-market ones | T3 | `ge_target`, `diag['ge_target_resid']` |
| `test_v1_werning_ge_target_resid_direct` | `ge_target_resid_direct_goods_{acyc,cyc}` (new), besides the asset-market ones | T3 | `ge_target`, dY, *`J_C`* (new, optional); dr from the caller |
| `test_v4_exogenous_G_as_requested` | `exog_G` | T5 paths | irfs dG; G from the caller |
| `test_v4_exogenous_r_as_requested` | `exog_r` (fixed-rate experiments) | T5 paths | irfs dr; r from the caller |
| `test_v4_exogenous_rshock_as_requested` | `exog_rshock` | T5 paths | irfs drshock; rshock from the caller |
| `test_v6a_f3_hetinput_two_sided` | as V6 (i) | – | SSJ only |
| `test_v6b_f3_hetoutput_two_sided` | as V6 (ii) | – | SSJ only |

- **Changed meaning:** `rate_rule` uses the requested rshock (section 3.3).
- **Diagnostics without threshold:** `goods_mkt_resid_direct_*`, `asset_mkt_resid_direct_*` (section 3.4).
- **Group V** after B has 37 tests: 33, minus the V6 test, plus V6a, V6b and the three exog tests. Every other test name is unchanged.
- **The checker regression tests** are in `tests/test_v0_3_0_checker.py`, outside group V.

## 5. Thresholds of the new quantities

Rule and τ classes: plan, section 4. The table is filled in by commit B from two unmodified runs of the B candidate: the full suite for the T5 quantities, and `suite.run_t3(ge_target='goods_mkt')` for the goods-market quantities.

**Filled in by commit B.** b_max from two runs of `scripts/thresholds_v0_3_0.py` on the B candidate (2026-10-22, 14:46 and 14:47 UTC). The candidate is a local commit, `e9fd7e1`, with the same `hank_limits/`, `tests/` and `pytest.ini` trees as B; B adds this table and the two run files `out/thresholds_v0_3_0_run{1,2}.json`. The script calls the check functions exactly as the pytest fixtures do. Environment as in v0.2: Python 3.13.16, sequence-jacobian 1.0.0, numpy 2.5.3 (scipy-openblas 0.3.34), scipy 1.18.1, numba 0.68.0. The two runs gave identical numbers. Path quantities are the largest over the four shocks and three horizons; T3 quantities the largest over the seven b.

| Model | Quantity | b_max (run 1, run 2) | τ (source) | Threshold |
|---|---|---|---|---|
| T3, goods-market target | `ge_target_resid_reported_goods_acyc` | 8.67e-14, 8.67e-14 | 1e-14 (exact up to rounding (linear solve)) | 9e-12 |
| T3, goods-market target | `ge_target_resid_reported_goods_cyc` | 1.39e-13, 1.39e-13 | 1e-14 (exact up to rounding (linear solve)) | 2e-11 |
| T3, goods-market target | `ge_target_resid_direct_goods_acyc` | 8.75e-14, 8.75e-14 | 1e-14 (exact up to rounding (linear solve)) | 9e-12 |
| T3, goods-market target | `ge_target_resid_direct_goods_cyc` | 1.27e-13, 1.27e-13 | 1e-14 (exact up to rounding (linear solve)) | 2e-11 |
| T5 | `exog_G` (path, max over labels) | 0, 0 | 1e-14 (exact up to rounding) | 1e-13 |
| T5 | `exog_r` (path, fixed-rate labels) | 0, 0 | 1e-14 (exact up to rounding) | 1e-13 |
| T5 | `exog_rshock` (path, max over labels) | 0, 0 | 1e-14 (exact up to rounding) | 1e-13 |

**Notes on the table**

- **The goods-market quantities are about ten times the asset-market ones** (b_max 1.1e-14 to 1.5e-14 there). Both are rounding-level residuals of a linear solve. Both are largest at b = 1e-4, where the division by b is largest. Before the division, the goods-market residual is about 1e-17 and the asset-market one about 1e-18.
- **Residuals of the other market,** recorded without threshold. Divided by b and maximised over the seven b:
  - the asset-market solve of the reference provider leaves a goods-market residual of 1.36e-7 (proportional income) and 2.52e-6 (countercyclical income);
  - the goods-market solve leaves an asset-market residual of 2.96e-6 and 2.44e-3. These equal, bit for bit, the V1 values with which C03 alarmed in v0.2.
- **The 88 existing quantities** have, in these runs, the b_max of the plan's table to three significant digits, except the T5 path `goods_market`: 3.11e-13 here, 3.12e-13 in the plan. This is the run-to-run rounding variation described below. Their thresholds are unchanged.

**Observations at B, before any mutant was run**

- **Full pytest on the candidate:** 107 passed, 4 xfailed. That is 77 tests in the suite's groups and the 30 checker tests of `tests/test_v0_3_0_checker.py`; the xfails are F1, F2, V6a and V6b.
- **The 276 v0.1 numbers.** Two baseline runs of the re-run runner (`scripts/mutation_v0_3_0.py --baseline-only v0.3`; output not kept) gave, in both runs, all 276 v0.1 summary numbers equal to v0.1 run B0. Against v0.1 run B1 only the known two-valued number differs. Test outcomes were identical.
- **Rounding-level differences between B0 and B1.** The two baselines differed in five T5 goods-market path residuals (dY − dC − dG; values between 1e-15 and 4e-6):
  - one is in the refined setting that the V tests use: `G_debt_r_rule` at T = 300, 3.11556e-13 in one run and 3.11383e-13 in the other, against a threshold of 4e-11;
  - the other four are in the two settings that no V test uses; two of them are the quantities that the plan (section 4) already saw vary between runs.
  - The refined-setting quantity also took the value 3.11383e-13 in the threshold runs, and 3.11556e-13 (the v0.2 value) in a pytest run in the repository. This is the run-to-run rounding variation of dC that the plan describes. No number of the v0.1 summaries and no test outcome is affected.
- **Effect on the predictions.** Prediction P1 ("B0 and B1 identical in every summarised number") and prediction P4 ("equal to the v0.2 value", where the refined setting is not exempted) are therefore expected to fail at rounding level for these quantities. The predictions stay as registered, and the outcome is reported in section 10.

## 6. Re-run design

- **Suite.** The v0.3.0 candidate is the tree of commit B (`hank_limits/`, `tests/`, `pytest.ini`, extracted with `git archive` as before).
- **Designs.**
  - (ii) the 20 v0.2 mutants and 3 controls × the candidate (23 runs);
  - (iii) the 28 v0.1 mutants and controls × the candidate (28 runs);
  - two unmutated baselines (B0, B1): 53 runs in all.
  - Design (i), the v0.2.1 suite, does not change and is not re-run; comparisons with it use its v0.2 results.
- **Runner.** `python3 scripts/mutation_v0_3_0.py --v03-ref <B> --design ii iii --out out/mutation_v0_3_0_rerun` (the copy of section 3.6). Detection rule, time limits, parallelism (4 runs, one BLAS thread each) and summaries are those of v0.2. The patches are those of the plan and of v0.1, re-anchored only as in section 3.6, and recorded in the output.
- **Outputs.** `out/mutation_v0_3_0_rerun.{json,csv}`.
- **Comparison with v0.2,** per run and group: verdicts, the tests that changed outcome, and the values of the v0.3 quantities.
- **The original results stay as the record of the checks as they were** (19/20, 25/25, the C03 alarm). This is a re-run of the suite as revised after review.

## 7. Predictions

**Basis.** The code changes above. For R1, a scan of `out/mutation_v0_2.json` made while writing this section: in the 53 runs of the v0.3 suite (both baselines, designs (ii) and (iii)), no v0.3 quantity is None, NaN or ±inf, no v0.3 computation recorded an error or a skip reason, and no V test was skipped. The one fixture error is t3 under N19 (root bracket, ValueError). For R3: the SSJ provider returns the requested exogenous paths because it takes them from the inputs, and no v0.1 or v0.2 patch touches `ssj_fiscal_provider` or the inputs. For the same reason the reported `drshock` equals the requested rshock in every run.

**P1. Baselines.**
- All tests pass except four strict xfails (F1, F2, V6a, V6b).
- B0 and B1 are identical in every outcome and every summarised number.
- All 276 numbers of the v0.1 summaries equal v0.1 run B0. Against v0.1 run B1 only the known two-valued number differs (plan, section 8).

**P2. Groups T1–T6 and D.** In all 51 runs, the same verdict and the same changed tests as in the v0.2 runs.

**P3. Existing V tests.** Every test of v0.2's group V except the V6 test has the same outcome as in v0.2, in every run, **except C03** (P6).

**P4. Values.** Every v0.3 quantity that exists in the v0.2 data equals its v0.2 value in the same run, except:
- C03's V1 quantities, now goods-market ones;
- possibly rounding-level differences in the two T5 settings that no V test uses (as seen between threshold runs in the plan, section 4).

**P5. New tests.**
- `exog_G`, `exog_r` and `exog_rshock` are exactly 0 in every run, and the three exog tests pass in every run: they catch nothing.
- V6a and V6b xfail in every run.

**P6. C03** (GE target goods market).
- The provider reports `goods_mkt`, and V1 is judged on the goods-market residual. Its values equal those of the goods-market threshold runs, so V1 passes, and C03 no longer alarms.
- **This is not independent evidence.** Its thresholds come from unmodified runs with the same setting as C03.
- The asset-market residual is now a diagnostic without threshold, and equals the values that alarmed in v0.2. Divided by b, the maximum over the seven b is 2.96e-6 with proportional income and 2.44e-3 with countercyclical income.

**P7. Suite level.**
- (ii): 19 of 20 mutants caught, N20 missed, no alarm on K01–K03.
- (iii): 25 of 25 caught, no alarm on C01–C03.
- Against the original pre-registered predictions: (ii) unchanged (167 of 168 firm cells held). (iii): the C03/V cell now differs from the pre-registered "catch (false alarm)", by design of R4.

**P8. Checker self-tests** (`tests/test_v0_3_0_checker.py`) pass in every run.

**Cells expected to differ from v0.2** (cells as in `docs/MUTATION_v0_2.md`, sections 3 and 9):

| Design | ID | Group | v0.2 | Predicted in the re-run | Why |
|---|---|---|---|---|---|
| (iii) | C03 | V | **2/33** (alarm: V1 reported and direct) | · (no alarm) | R4: judged on the reported goods-market target, with thresholds from the same setting |
| (ii), (iii) | every run | V | n/33 | the same n, out of 37 | V6 split, three exog tests; none expected to change outcome |

No other cell is expected to change. R1 to R3 change how missing, non-finite and wrongly sized exogenous values are judged. No run of the v0.2 data produced such values, so their evidence is the regression tests with fake providers (sections 3.1–3.3), not this re-run.

## 8. Documentation (P3)

| Item | Where | Commit |
|---|---|---|
| `docs/MUTATION_v0_2.md` line 17: "Every V catch … factor of 1e4 or more" is wrong. N18's largest margin is about 7,577 (asset market) and 9,742 (goods market). By failing test, N12's TANK goods-market test is about 526. Recompute from the data, separating the largest margin per mutant from the smallest margin per failing test. Dated correction; no effect on 19/20 | `docs/MUTATION_v0_2.md` | D |
| Same file, line 242: "any goods-market solve" → what was observed with this reference provider, calibration and solver setting. Dated correction | `docs/MUTATION_v0_2.md` | D |
| `models.py` (~line 310): the `ssj_werning_provider` docstring says "(goods market)"; the target is `asset_mkt` | `models.py` | B |
| V1 rows: target convention and scope (R4) | both READMEs | B |
| V3 rows: aggregate income = Y, with the reference per model: T2 income per type and the weighted sum, T4 w L, T5 Y − Tax | both READMEs | E |
| "bit-identical": limited to "in the report's baseline runs, the 276 existing summary numbers were exactly equal" | both READMEs | E |
| "41 passed, 2 xfailed" (README.md ~line 69, README_ko.md ~line 68) and section 11: marked as the v0.2.1 result, with a link to the current result | both READMEs | E |

## 9. What can and cannot be claimed

- The re-run shows whether the revision changed detection on the same two samples (20 + 3 v0.2 runs, 28 v0.1 runs). It is not new independent evidence beyond v0.2.
- If C03 no longer alarms, that holds by construction (section 3.4); it shows only that the target convention works as built.
- R1–R3 address failure modes (missing, non-finite and wrongly sized exogenous values) that neither the reference provider nor the mutants produced. Their evidence is the regression tests with fake providers.
- The thresholds of the new quantities are regression thresholds for this calibration, like the others.

## 10. Results (commit D)

Data: `out/mutation_v0_3_0_rerun.json` (every run's outcomes and check summaries) and `out/mutation_v0_3_0_rerun.csv` (one row per design × mutant × group). Comparison with v0.2: `python3 scripts/rerun_v0_3_0_compare.py --json out/rerun_v0_3_0_compare.json`.

**These are results of a re-run of the suite as revised after the review.** The v0.2 results in `docs/MUTATION_v0_2.md` (19 of 20, 25 of 25, the C03 alarm) remain the record for the checks as they were then.

### 10.1 Execution record (2026-10-22, UTC)

| Step | Commit | Time | Content |
|---|---|---|---|
| Design and predictions A | `471fc47` | 14:28:46 | this file, sections 1–9; one line in the plan, section 10 |
| Threshold runs | (candidate `e9fd7e1`) | 14:46, 14:47 | `scripts/thresholds_v0_3_0.py` twice; one baseline-only run of the runner (section 5) |
| Code, tests, thresholds B | `a290c75` | 14:51:57 | R1–R4, V6 split, thresholds, runner copy |
| Re-run | – | 14:52:08 start, 1071 s | 53 runs, 4 in parallel, 68–152 s each |
| Outputs C | `ad1857d` | 15:10:12 | `out/mutation_v0_3_0_rerun.{json,csv}` |
| Results D | this commit | | this section; corrections in `docs/MUTATION_v0_2.md` |

- The runner record shows HEAD `a290c75` with a clean worktree, and `runner_committed = true` (the copy committed in B).
- The environment is that of v0.2 (section 5).
- Re-anchored patches, as registered in section 3.6: M04 and C03, recorded under `reanchored` in the JSON. No other patch changed. Before B, the unmodified v0.2 runner failed on exactly these two patches, and on no other (`--check`).
- No mutant was run before B was pushed.

### 10.2 Baselines

- Both baselines: 107 passed, 4 xfailed (F1, F2, V6a, V6b).
- B0 and B1 are identical in every outcome and every summarised number; the runner reports `deterministic = true`.
- **All 276 numbers of the v0.1 summaries equal v0.1 run B0 bit for bit, in both baselines.** Against v0.1 run B1 only the known two-valued number differs (`t5.ssj_inconsistent_options.walras_max.rshock_r_rule`, relative 1.4e-11).

### 10.3 Changes against v0.2

Cells as in `docs/MUTATION_v0_2.md`, sections 3 and 9: per group, whether the mutant was caught and which tests changed outcome.

| Design | Runs | Cells that differ from v0.2 (detection or the set of changed tests) |
|---|---|---|
| (ii) | 23 | **none**. Every mutant and control has the same verdict in every group, and the same changed tests, as in v0.2. |
| (iii) | 28 | **one**: C03, group V. v0.2: 2/33 (`test_v1_werning_ge_target_resid_reported`, `..._direct`); re-run: 0/37. |

**Suite level.**

- (ii): 19 of 20 mutants caught, N20 missed. Five of them are caught only by V (N04, N05, N07, N10, N18). No control alarms (K01–K03).
- (iii): 25 of 25 caught, four only by V (M12, M21, M23, M24). No control alarms (C01, C02, C03).
- Group-V counts are those of v0.2 with 37 in the denominator. For example, in (ii): N02 11/37, N06 9/37, N12 9/37, N19 8/37 (N19: the t3 fixture error, as before), N17 4/37. In (iii): M13 14/37, M08 13/37.

**Per V test** (the table of `docs/MUTATION_v0_2.md`, section 5):

- Each existing V test catches the same mutants as in v0.2, except that C03 no longer appears under the two V1 tests in (iii).
- The new tests caught nothing: `test_v4_exogenous_G_as_requested`, `..._r_...`, `..._rshock_...` passed and `test_v6a_...`, `test_v6b_...` xfailed in all 53 runs.

**Pre-registered predictions of v0.2, scored by the runner.**

- (ii): 167 of 168 firm cells held. The one miss is N13/T2, as in v0.2.
- (iii): 222 of 223 held. The one miss is C03/V, predicted "catch (false alarm)" in v0.2, which R4 was designed to change.

### 10.4 Predictions of this amendment (section 7)

| Prediction | Result | Detail |
|---|---|---|
| P1 baselines | **held** | 107 passed, 4 xfailed; B0 = B1 in every outcome and summarised number; 276 = v0.1 B0 in both. (The two baselines of the runner check at B differed at rounding level, section 5; these two did not.) |
| P2 T1–T6, D | **held** | same verdict and same changed tests as v0.2 in all 51 runs |
| P3 existing V tests | **held** | same outcomes as v0.2 in every run except C03 |
| P4 values | **failed at rounding level** | see below |
| P5 new tests | **held** | `exog_G`, `exog_r`, `exog_rshock` exactly 0 in every run; the three tests pass, V6a and V6b xfail, in every run |
| P6 C03 | **held** | the provider reports `goods_mkt` for both incidences; V1 passes. The four goods-market values equal the threshold runs' b_max bit for bit (8.67e-14, 1.39e-13, 8.75e-14, 1.27e-13 against 9e-12 and 2e-11). The asset-market diagnostic equals the values that alarmed in v0.2 (2.96e-6 and 2.44e-3, divided by b) |
| P7 suite level | **held** | (ii) 19/20, N20 missed, no alarm on K01–K03; (iii) 25/25, no alarm on C01–C03; original predictions 167/168 and 222/223 (C03/V) |
| P8 checker self-tests | **held** | the 30 tests of `tests/test_v0_3_0_checker.py` passed in all 53 runs |

**P4 in detail.** Every v0.3 quantity in the v0.2 data was compared with the same run's value in the re-run. Only one kind of quantity differs, the T5 goods-market path residual max_t |dY_t − dC_t − dG_t| / s; there are no differences in any other quantity, and none in the v0.1 summaries apart from the known two-valued number.

- **Runs affected:** 31 of the 53. The baselines differ from v0.2's in five such values, which are the ones seen in section 5.
- **Settings affected:** all three T5 settings.
  - In the refined setting that the V tests use: 27 runs, 43 values, absolute differences at most 2.8e-15. The threshold is 4e-11, and the baseline values are near 3e-13.
  - In the two settings that no V test uses: up to 5.8e-9 in absolute terms, but at most 4.5e-12 relative (M08, whose residuals there are of order 1e3).
- **Cause.** This is the run-to-run rounding of dC described in the plan (section 4), which also shows in the v0.1 two-valued number. The prediction exempted only the two unused settings, and was wrong for the refined one.
- **No test outcome changed because of it.** Section 5 had flagged this risk before any mutant was run.

**C03 is not independent evidence.** Its V1 thresholds come from unmodified runs of the same setting (section 3.4), so its pass holds by construction. It shows that the target convention works as built: a goods-market solve is judged on the goods market, and its asset-market residual is kept as a diagnostic.

### 10.5 Contribution of the new quantities

| Quantity | (ii) | (iii) | Note |
|---|---|---|---|
| `exog_G`, `exog_r`, `exog_rshock` | 0 catches; value 0 in every run | 0 catches; value 0 in every run | the reference provider takes the exogenous paths from the inputs, and no mutant changes that; the evidence for R3 is the doubled-shock regression test |
| `ge_target_resid_*_goods_*` | not applicable (asset-market target in every run) | C03 only: passes | the four values equal the threshold runs' |
| `rate_rule` with the requested rshock | same catches as v0.2 (N04, N16) | same as v0.2 (M05, M20, M25) | the reported drshock equals the requested rshock in every run |
| R1 verdict (non-finite, missing, incomplete) | no run had a missing or non-finite V value | as (ii) | evidence: regression tests (section 3.1) |
| R2 per-equation fields | every field reported in every run | as (ii) | evidence: fake providers (section 3.2) |
| V6a, V6b | xfail in every run | xfail in every run | SSJ code only; no model-side mutant reaches them |

### 10.6 Controls

- No control changed any test outcome in the re-run: K01–K03 in (ii), C01–C03 in (iii).
- C03 no longer alarms, for the reason in section 10.4. Six controls do not make a false-alarm rate.

### 10.7 Margins

- The margins of the V catches are the same as in v0.2. The values are bit-identical, apart from the rounding-level goods-market path residuals of section 10.4, none of which is a catching value.
- In (ii), every caught mutant has a failing V test above its threshold by a factor of 9.7e3 or more (N18). By individual failing test, the smallest factor is 526 (N12, TANK goods market).
- In (iii), the smallest of the largest factors per mutant is 274 (M16). By individual failing test, the smallest factor is 1.3 (M15, V1).
- These figures correct a statement in `docs/MUTATION_v0_2.md`, section 1; its dated correction is part of this commit.

### 10.8 What this shows

- The revision changed one detection outcome on the two samples, the one it was designed to change (C03). The 19 of 20 and 25 of 25 of v0.2 hold for the revised suite as well.
- The evidence for R1–R3 is in the regression tests with fake providers: missing, non-finite and wrongly sized exogenous values. Neither the reference provider nor any mutant produced such a case.
- The statements of `docs/CHECKS_PLAN_v0_3.md`, section 9, and of `docs/MUTATION_v0_2.md`, section 10, hold unchanged: two hand-picked samples, design (iii) not independent, one environment.
