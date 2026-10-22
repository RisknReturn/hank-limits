# v0.3.0 release candidate: review-driven changes to the v0.3 checks, and their re-run

Status: **design and predictions** (commit A on branch `v0_3_0-rc`). Written and pushed before any code of this amendment and before any re-run. Code, tests and the thresholds of the new quantities follow in commit B (section 5 is filled in there, before any mutant is run); the re-run outputs in commit C; the results in commit D (section 10). Release metadata (version 0.3.0) is commit E.

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

*(filled in by commit B)*

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

## 10. Results

*(commit D)*
