# hank-limits v0.2.1: analytic-limit and aggregate-identity tests for HANK codes

[![DOI](https://zenodo.org/badge/1426447298.svg)](https://doi.org/10.5281/zenodo.23489176)

Cite all versions with the concept DOI [10.5281/zenodo.23489176](https://doi.org/10.5281/zenodo.23489176); this release (v0.2.1) is [10.5281/zenodo.23489177](https://doi.org/10.5281/zenodo.23489177). GitHub's "Cite this repository" button uses `CITATION.cff`.

This repository collects tests that a heterogeneous-agent New Keynesian (HANK) code should pass at limits with known answers and at exact aggregate identities. It is implemented on the public package `sequence-jacobian` (SSJ) 1.0.0. Each test reports not only pass/fail at one point but **the error as the limit parameter goes to its limit and the observed order of convergence**.

- Licence: MIT (`LICENSE`). Archived on Zenodo (concept DOI [10.5281/zenodo.23489176](https://doi.org/10.5281/zenodo.23489176)). **Not peer reviewed.**
- Model-independent check functions (`hank_limits/checks.py`) are separated from the SSJ reference implementation (`hank_limits/models.py`), so the checks can be reused as a validation harness for other HANK codes.
- **What this suite does and does not show.** A pass is numerical evidence that a specific, limited statement holds in this implementation. It is not a statistical significance test, and it does not certify the predictive performance of HANK models or the accuracy of every HANK code.
- Korean version of this README: `README_ko.md`.

### Changes v0.1 → v0.2 (review R4 by GPT 6.1 Sol; see section 17)

| Item | v0.2 change |
|---|---|
| Names and scope | T1 = PIH **household-block** limit (not a full NK general-equilibrium RANK limit). T2 = **TANK benchmark with imposed Bilbiie-style income incidence** (χ is not derived structurally). T4 = **no-idiosyncratic-risk limit** (not a complete-markets limit). |
| T1 horizon | The T = 50 and 100 errors are sub-blocks (observation windows) of the same T = 200 Jacobian. Jacobians computed with those horizons are **exactly equal** to the sub-blocks (difference 0). The wording "truncation horizon" was corrected. |
| T3 | Aggregate-income preservation under countercyclical incidence holds **only to first order** at the steady state; its size is reported. The transition matrix Π is fixed, so cyclical employment probabilities are not replicated. |
| T4 | `r < 1/β − 1` is forced by the root-search bracket and is no longer counted as independent evidence. Added capital-market residual, Euler/KKT residuals, top-of-grid mass, and asset-limit and resolution sensitivity. Added a comparison with two cells of Aiyagari (1994) Table II. |
| T6 | States that the divergence of group-wise general-equilibrium attribution is a property of forcing that violates the zero-present-value condition, not an "impossibility theorem". |
| FastHetBlock | Diagnostics added for negative mass, termination of inverse iteration, sensitivity to the shift and the iteration count, and comparison with SSJ's own forward iteration (`hank_limits/diagnostics.py`). |
| F1, F2 | Minimal reproduction scripts that use only SSJ code; for F2, reordering comparisons and normalised backward error; versions and wheel hash. |
| Tables | The numbers of the main tests and grid studies (239 rows) are exported as CSV with columns for limit parameter, grid, difference step, observation window and solver (`out/tables_v0_2/`). Diagnostics, difference steps, Aiyagari, F1 and F2 results are in their JSON files. |

**v0.2 → v0.2.1 (R4 follow-up review, same reviewer):** (1) removed the diagnostic field `normalisation_err_last_sweep`, which was always zero by construction; (2) the second Aiyagari cell is now reported as "not reproduced, cause unresolved", with a comparison against the level CV and autocorrelation of the original Table I and a grid audit; (3) the practical impact of F2 is limited to this example; (4) horizon independence in T1 is limited to the linear household Jacobian at a fixed steady state; (5) revised wording of the SSJ issue drafts. Test results are unchanged.

---

## 1. Summary

The error metric is the **relative max-entry error** `max|x − ref| / max|ref|`. It is not an induced matrix norm and does not guarantee relative accuracy of small tail entries. For references near zero, `metrics.check_close` combines absolute and relative tolerances. Convergence orders are the neighbour-pair order `log(e_i/e_{i+1}) / log(h_i/h_{i+1})` and the least-squares slope on a log-log scale (fitted slope).

| Test | Statement (scope) | Limit parameter | Error at the tightest setting | Observed order (theory) | Status |
|---|---|---|---|---|---|
| T1 PIH household-block limit | HA household Jacobians dC/dY, dC/dr → permanent-income closed form (partial equilibrium) | income risk eps → 0 (β recalibrated along the path) | dC/dY 9.5e-4, dC/dr 4.0e-4 (eps = 0.0125); 1.7e-9 at eps = 0 | 1.95 / 2.03 (2) | pass |
| T2 TANK benchmark, imposed incidence | aggregate iMPC = λχ I + (1−λχ) M_RA; GE output = (1−λ)/(1−λχ) × RANK | eps → 0 (both types) | aggregate iMPC 8.0e-5, IRF (χ = 1.5) 1.05e-3 | savers 1.95, IRF 1.79 (2) | pass |
| T3 Werning vanishing liquidity | with proportional income, the output response to a real-rate path = RANK; generalised Euler relation when income elasticities vary | borrowing limit b → 0 | proportional 4.3e-4; countercyclical 6.0e-4 (deviation 0.319, closed form 0.320) | 0.99 / 0.99 (1) | pass |
| T4 Aiyagari steady state | r → 1/β − 1 (idiosyncratic risk → 0); monotone in risk size and persistence in this calibration (regression check) | log-income sd eps → 0 | gap 6.2e-5 (eps = 0.025) | 1.93 (2) | pass |
| T5 IKC consistency | SSJ general equilibrium dY = (I − M)^{-1} dC_partial (consistency test built from the same block Jacobians) | identity (T = 100, 200, 300) | refined solve 9.8e-12; SSJ default solve 5.9e-7 | n/a | pass + 2 expected failures (F1, F2) |
| T6 HANK − RANK decomposition | finite-T resolvent identity and differencing decomposition along the same paths | identity (T = 150, 300) | identity 2.2e-14, decomposition 1.9e-10 (G) / 2.2e-10 (rate) | n/a | pass |

**Interpreting the orders.** The O(eps²) in T1 and T2 is the order along a model path on which risk shrinks while β is recalibrated; in T4, β is fixed at 0.96 and O(eps²) is the order of the gap 1/β − 1 − r as risk shrinks. Neither is evidence that the spatial discretisation is second-order accurate. The O(b) in T3 is also the approximation order of the liquidity limit. Do not read T1 as saying that every model with vanishing risk and a fixed β converges to the PIH solution.

**Status of the tolerances.** Test tolerances are **regression-check thresholds for this calibration**. They are not a conformance specification requiring other implementations to show the same slopes or the same 1e-10 threshold.

Full pytest run: **41 passed, 2 xfailed, about 65–75 s** (section 11). Report script about 110 s, diagnostics script about 80 s.

### Key findings

1. **T1: observation-window dependence.** The difference between the HA and PIH Jacobians concentrates in long-lag entries (t = 199, s = 0): PIH wealth has a unit root while HA wealth mean-reverts at rate O(eps²). At eps = 0.0125 the dC/dY errors for observation windows 50, 100 and 200 are 2.0e-4, 3.8e-4 and 9.5e-4. Jacobians computed with horizons T = 50 and 100 are **entry-by-entry identical** to the leading blocks of the T = 200 Jacobian (max difference 0.0). The entries of SSJ household Jacobians do not depend on the computation horizon, so this is not a truncation effect: the error grows as longer lags enter the observation window. The eps → 0 limit holds for a fixed window but not when the window grows without bound at fixed eps. This horizon independence is a result for the **linear household Jacobian at a fixed steady state**; it does not mean that the inverse of the GE H_U, endogenous price paths, nonlinear transitions or terminal conditions are independent of T.
2. **Source of the discretisation error.** The order dropping slightly below 2 at the smallest eps happens because the numerical diffusion of Young's lottery method becomes comparable to the true diffusion when eps is small. Refining the grid raises the order towards 2 (T2 IRF: 1.47, 1.72, 1.85 at 2400, 4800, 9600 points; T1 dC/dr 2.006 at 9600 points).
3. **F1: consistency between cached Jacobians and output propagation (SSJ 1.0.0, master eb4ac4c).** `HetBlock._impulse_linear` does not use the `Js` passed to it and re-differentiates with its own options. In a minimal reproduction that uses only SSJ code, passing `Js` multiplied by 2 **does not change** the output (difference 0). The GE solve (H_U) and the propagation of direct shocks and outputs then use different linear operators, so the IKC identity for a rate shock that enters the household directly holds only to about 1e-4. With matched options it holds to 3e-13. This is a problem with the reuse convention, not with SSJ's differentiation theory; matching the options is a workaround, not the same as reusing the cache. Reported upstream as [#39](https://github.com/shade-econ/sequence-jacobian/issues/39).
4. **F2: LU stability that depends on equation ordering.** In a fiscal HANK with an interest-rate rule, partial-pivoting LU on H_U as SSJ stacks it (blocks by variable) shows element growth of 19, 3.7e4 and 1.1e7 at T = 100, 200, 300 (condition numbers 3.6e2, 2.9e3, 2.1e4). **Reordering the same system by date gives growth 1.02; reversing the block order gives 1.00**, and the normalised backward error falls from 3e-8 to about 1e-16. The condition number alone does not explain the problem: LU element growth in this ordering causes the loss of accuracy. So the statement is not "SSJ solves are inaccurate at long T" but **"element growth arises for this block ordering combined with partial-pivoting LU"**. We do not claim that reordering, iterative refinement or QR fixes every model. In this example the IKC error at T = 200 and 300 is 2.2e-9 to 5.9e-7, and the relative difference from a pivoted-QR solution is about 1e-7 (1.0e-7 here, 1.2e-7 in the reviewing model's environment). This matters for strict identity checks (1e-10). The practical importance depends on the model, target and tolerance and is not compared here. Reported upstream as [#40](https://github.com/shade-econ/sequence-jacobian/issues/40).
5. **T6: scope of group-wise GE attribution.** First-round gaps by wealth quartile are well defined and independent of T. But the group-wise forcing has non-zero present value (income is (Y − Tax) × e and the composition of e differs by group), so passing each group's forcing through the infinite-horizon consumption multiplier excites a mode that is inconsistent with the budget constraint and diverges with T. The mode cancels in the total. This is not a theorem that group-wise attribution is impossible. Group-wise GE attribution requires additional choices of terminal conditions, fiscal incidence and present-value compensation paths, and the attribution depends on those choices. Zero present value is a necessary compatibility condition only; it does not guarantee stability or uniqueness.
6. **A small stationarity residual does not guarantee distributional accuracy (an observation, not an SSJ defect).** In a slowly mixing chain the inverse operator on the zero-mass subspace is large, so a small residual does not imply a small distribution error. For the medium-risk PIH household (eps = 0.1), SSJ's original forward iteration with its default stopping rule (every 10 iterations, max|ΔDbeg| < 1e-10) stops after 37,910 iterations, and the mean assets of that distribution differ from the inverse-iteration solution by **1.7e-3 (relative)**. Tightening the tolerance to 1e-12 stops after 54,590 iterations with a difference of 1.7e-5, converging to the same solution. At eps = 0.0125 it does not stop within 600,000 iterations. In the fast-mixing T4 and T5 models the differences are 1e-6 and 5e-9 (section 12).

---

## 2. How to run

```bash
pip install sequence-jacobian==1.0.0 pytest matplotlib   # needs numpy, scipy, numba
cd hank-limits
python3 -m pytest                          # about 74 s; summary output with python3 -m pytest -s
python3 scripts/convergence_report.py      # about 110 s: out/convergence_v0_2.json, 6 PNGs, out/tables_v0_2/*.csv
python3 scripts/diagnostics_v0_2.py        # about 80 s: out/diagnostics_v0_2.json (section 12)
python3 scripts/aiyagari_table2_v0_2.py    # about 20 s: out/aiyagari_table2_v0_2.json (section 7)
python3 scripts/repro_F1_v0_2.py           # SSJ code only: out/repro_F1_v0_2.json
python3 scripts/repro_F2_v0_2.py           # SSJ code only: out/repro_F2_v0_2.json, out/F2_HU_T300_v0_2.npz
```

- Every entry point limits threads to one (`OMP_NUM_THREADS=1` etc.).
- Versions used: **sequence-jacobian 1.0.0** (PyPI wheel SHA-256 `40715634f6fcf8f695d79687c471bedcecacdeff94522557f5186b630f46ef01`), Python 3.13.16, numpy 2.5.3, scipy 1.18.1, numba 0.68.0, OpenBLAS 0.3.34 (scipy-openblas), pytest 9.1.1. GitHub master (commit `eb4ac4c3157357a4cacbea4b5a7a281f5b48cab9`) has the same code as F1.
- The reviewing model (GPT 6.1 Sol, section 17) re-ran v0.1 in its own environment (Python 3.12.14, numpy 2.5.4, scipy 1.17.0) and reproduced **32 passed, 2 xfailed**. Last digits and timings may vary with BLAS and environment.

### Layout

```
hank_limits/
  core.py         household block, FastHetBlock (direct stationary-distribution solve, with diagnostics), grids
  diagnostics.py  v0.2: Euler/KKT residuals, distribution audit (shift and iteration sensitivity, SSJ forward iteration)
  ra.py           RA, TANK and zero-liquidity closed forms
  models.py       SSJ reference implementation (providers)
  checks.py       model-independent check functions
  interface.py    provider conventions
  metrics.py      error metrics, convergence orders, mixed absolute/relative tolerance
  suite.py        standard settings (shared by tests and scripts)
tests/            test_t1 ... test_t6, test_diagnostics_v0_2
scripts/          convergence_report.py, diagnostics_v0_2.py, aiyagari_table2_v0_2.py,
                  repro_F1_v0_2.py, repro_F2_v0_2.py, a1_markov_checks_v0_2.py (for the design note, section 13)
out/              JSON, PNG, tables_v0_2/*.csv, pytest output
```

---

## 3. Conventions and numerics

- **Timing.** SSJ one-asset household: budget c_t + a_t = (1 + r_t) a_{t−1} + y_t, Euler equation c_t^{−1/σ} = β (1 + r_{t+1}) E c_{t+1}^{−1/σ}. r_t is the ex-post real return paid in period t on assets held at t−1.
- **Annual calibration.** The suite uses an annual calibration with steady-state output Y = 1. Do not mix it with the quarterly UK calibration note.
- **Jacobians.** SSJ fake-news algorithm with **two-sided differences**, h = 1e-4. SSJ's default one-sided differences leave an O(h) error of about 5e-5 in dC/dr. Step-size sensitivity is in section 4.
- **`FastHetBlock`.** A subclass of SSJ's `HetBlock` that changes only how the stationary distribution is found. The unit-eigenvalue problem of Young's lottery transition matrix is solved by sparse inverse iteration (one LU of (1+τ)I − L, τ = 1e-10, asset-major ordering to keep the matrix banded), and backward iteration starts from the previous solution. Everything else (fake-news algorithm, GE solves) is SSJ as is. Because this is our own numerical component, v0.2 records the following (`FastHetBlock._ss_diag`, section 12):
  - negative mass before taking absolute values (relative to the total absolute mass), number of inverse iterations and whether the stopping rule was met, size of the last update, mass drift in the polishing steps. (The v0.2 field `normalisation_err_last_sweep` was always zero by construction and was removed in v0.2.1; negative mass is reported by `neg_mass_rel_max`.)
  - L1 and max-entry residuals of the resulting distribution, for the clipped chain and for SSJ's own forward operator.
  - changes in the distribution, mean assets and consumption when the shift τ (1e-12 to 1e-6) and the iteration count vary.
  - comparison with SSJ's original forward iteration (`HetBlock.forward_steady_state`) at the same policy.
- **Policies off the grid.** The sparse transition gathers policies beyond the top of the grid at the end point (reflecting barrier), while SSJ's forward operator is evaluated as is. The mass affected (`mass_policy_off_grid`) and the SSJ residual are always recorded in `FastHetBlock._ss_diag`. At the T4 points recorded in the v0.2 diagnostics it is at most 3.1e-30 at the default upper limit (200) and at most 8.7e-13 in the sensitivity check with limit 100 (`out/diagnostics_v0_2.json` B).

---

## 4. T1: household-block permanent-income limit

**Scope.** An exact Jacobian comparison for the riskless PIH household block, and a convergence test along a specific eps path on which β is recalibrated to an asset target. **This is a partial-equilibrium household-block test, not a test of the RANK limit of a full NK general equilibrium.**

**Statement.** A household with no income risk, β(1+r) = 1 and a non-binding borrowing constraint is a permanent-income consumer, and with q = 1/(1+r)
- dC_t/dY_s = (1−q) q^s
- dC_t/dr_s = (1−q) A q^s + σ C q (1{t ≥ s} − q^s)

(derived from the Euler equation and the intertemporal budget constraint, `hank_limits/ra.py`). Source: linearisation of the standard permanent-income hypothesis. SSJ methodology: Auclert, Bardóczy, Rognlie, Straub (2021).

**Exact setting.** r = 0.03, σ (EIS) = 0.5, y = Y e, Rouwenhorst 3 states, persistence 0.9, unconditional log-income sd = eps. Borrowing limit = 0.99 × the natural borrowing limit −y_min/r (constrained mass at most 3e-6). β is calibrated so that mean assets A = 5, and κ = 1 − β(1+r) goes to 0 like eps² (κ/eps² = 0.028). eps = 0 is the exact permanent-income case (β(1+r) = 1, all mass at A = 5; the stationary wealth distribution is not unique in this case, so a point mass is chosen). T = 200, asset grid 4800 points (linear), upper limit 250.

**Results.**

| eps | 0 | 0.2 | 0.1 | 0.05 | 0.025 | 0.0125 |
|---|---|---|---|---|---|---|
| dC/dY error | 1.7e-9 | 2.08e-1 | 5.52e-2 | 1.40e-2 | 3.60e-3 | **9.53e-4** |
| dC/dr error | 2.9e-9 | 1.12e-1 | 2.63e-2 | 6.44e-3 | 1.60e-3 | **3.96e-4** |
| κ | 0 | 1.14e-3 | 2.74e-4 | 6.79e-5 | 1.71e-5 | 4.35e-6 |

- Neighbour orders: dC/dY 1.92, 1.98, 1.96, 1.92 (least squares 1.95); dC/dr 2.09, 2.03, 2.01, 2.01 (2.03).
- The 1e-9 error at eps = 0 is only the O(h²) error of the finite differences; it shows that SSJ's timing convention and the closed form agree.

**Observation window and computation horizon (v0.2).** eps = 0.0125:

| Observation window (leading block of the same T = 200 Jacobian) | 50 | 100 | 200 |
|---|---|---|---|
| dC/dY error | 1.96e-4 | 3.79e-4 | 9.53e-4 |
| dC/dr error | 9.25e-5 | 9.60e-5 | 3.96e-4 |
| Max entry difference from a Jacobian computed with that horizon | 0.0 | 0.0 | — |

The same holds at eps = 0.1 (`out/diagnostics_v0_2.json` C). Terminal-horizon effects need a separate audit for the GE solves (T5, T6) (for example, the T6 RANK multiplier of 1 + 6.4e-6 comes from the truncated tail of the debt path at T = 150).

**Difference step h sensitivity (two-sided differences, T = 200, v0.2).**

| eps | h = 1e-3 | 1e-4 (default) | 1e-5 | 1e-6 |
|---|---|---|---|---|
| 0: dC/dY / dC/dr error | 2.8e-9 / undefined* | 1.7e-9 / 2.9e-9 | 1.2e-8 / 9.7e-10 | 2.3e-7 / 6.7e-9 |
| 0.0125: dC/dY / dC/dr error | 9.527886e-4 / undefined* | 9.527901e-4 / 3.9569e-4 | 9.527900e-4 / 3.9608e-4 | 9.527900e-4 / 3.9611e-4 |

\* Raising r by 1e-3 moves the natural borrowing limit −y_min/(r+h) above the grid's lower bound (0.99 × −y_min/r), so consumption becomes negative. With a near-natural-limit grid the difference step in r must be below about 3e-4. At eps = 0.0125 the dC/dY error is the same to within 1.1e-7 (relative) for h ≤ 1e-4, so it is independent of the step; at eps = 0, rounding error grows as h shrinks.

**Grid convergence.** eps = 0.0125 dC/dY error: 1200 points 1.41e-3, 2400 points 1.08e-3, 4800 points 9.53e-4, 9600 points 9.01e-4. At 9600 points the eps 0.025 → 0.0125 order is 1.96 for dC/dY and 2.006 for dC/dr.

**Known limitations (v0.1 results).** The standard deviation of the limiting distribution stays wide (about 14, in units of annual income) regardless of eps, so a near-natural borrowing limit and a wide grid are needed. With a limit of 0.95 × the natural limit (600 points), a constrained mass of 1.8e-4 creates an error floor (about 6e-3). For eps below 0.0125, numerical diffusion exceeds true diffusion and finer grids are needed.

---

## 5. T2: TANK benchmark with imposed Bilbiie-style incidence

**Scope.** A test of TANK aggregation and real-rate transmission for a given income incidence χ. χ is not derived from the labour market or profit distribution, so this is not a "structural replication of the Bilbiie model". Imposing χ exogenously keeps its value as a benchmark.

**Statement.** Hand-to-mouth households (H) with share λ and savers (S) with share 1−λ. If H income has elasticity χ with respect to aggregate income, then (i) the aggregate iMPC matrix entering general equilibrium is λχ I + (1 − λχ) M_RA (λ I + (1−λ) M_RA when χ = 1); (ii) with zero liquidity (zero net bond supply) the aggregate Euler-IS equation is c_t = c_{t+1} − σ (1−λ)/(1−λχ) r_{t+1}, and the GE output response to a fixed real-rate path is (1−λ)/(1−λχ) times RANK.

**Conditions.** Both types have the same steady-state consumption and EIS, type H is exactly hand-to-mouth in the limit, the real-rate path is fixed, and **1 − λχ > 0** (here λ = 0.3, χ ∈ {1, 1.5}, so 1 − λχ = 0.55). The amplification factor 0.7/0.55 = 1.27273 is the value at λ = 0.3, χ = 1.5.

**Sources.** Bilbiie (2020), "The New Keynesian cross", JME 114, 90-108 (the aggregate Euler-IS equation and "amplification only when χ > 1" checked in the 2018 ECB conference version; Bilbiie's r_t is r_{t+1} in SSJ timing). Bilbiie (2008), JET 140(1), 162-196.

**Exact setting.** λ = 0.3. Type H: borrowing limit 0, β_H = 0.95, income risk eps, 500 points. Type S: the same household as T1 (4800 points), β_S calibrated so that λ A_H + (1−λ) A_S = 0. Incidence: Y_H = Y^χ, Y_S = (Y − λ Y_H)/(1−λ). General equilibrium: dr_t = −0.01 × 0.8^t, Y = C, unknown Y, T = 200. The same eps → 0 for both types.

**Results.**

| eps | 0.2 | 0.1 | 0.05 | 0.025 | 0.0125 |
|---|---|---|---|---|---|
| type H: M_H vs I | 0.596 | 0.342 | 3.9e-13 | 1.7e-13 | 3.9e-13 |
| type S: M_S vs M_RA | 0.275 | 7.29e-2 | 1.85e-2 | 4.76e-3 | **1.26e-3** |
| aggregate iMPC (χ = 1) vs λ I + (1−λ) M_RA | 0.556 | 0.320 | 1.18e-3 | 3.03e-4 | **8.02e-5** |
| GE IRF vs (1−λ)/(1−λχ) × RANK (χ = 1.5) | 0.131 | 5.44e-2 | 1.24e-2 | 3.46e-3 | **1.05e-3** |
| HA amplification (χ = 1.5 / χ = 1) | 1.150 | 1.209 | 1.27273 | 1.27273 | 1.27273 |

At χ = 1.5 the difference between λχ M_H + (1−λχ) M_S and λχ I + (1−λχ) M_RA is 4.3e-5 at eps = 0.0125 (`err_Magg_lamchi_vs_TANK_chi1.5` in `out/tables_v0_2/`). Type H becomes exactly hand-to-mouth once eps falls below a finite threshold. Grid convergence (eps = 0.0125 IRF): 2400 points 1.50e-3, 4800 points 1.05e-3, 9600 points 8.70e-4.

**Known limitations.** χ is a reduced form. With zero net bond supply the wealth-effect term in dC/dr is not tested here (it is tested in T1). Fitting a slope to the aggregate iMPC alone is not meaningful, so only neighbour orders are reported.

---

## 6. T3: Werning (2015) vanishing-liquidity benchmark

**Scope.** Tests the vanishing-liquidity refinement, Werning's neutrality under proportional income, and the generalised Euler relation when income elasticities vary. The countercyclical case is the **variation in income elasticities** of Werning's section 3.3; the model of section 3.4 with varying employment probabilities is not implemented. The transition matrix Π is fixed, so cyclical employment probabilities are not replicated.

**Statement.** With no borrowing and no outside assets, everyone consumes their income in equilibrium and the real rate is set by the Euler equation of the household most willing to save (Werning eq. 12). At exactly zero liquidity the rate can be indeterminate within an interval, but the vanishing-liquidity refinement selects the marginal rate. With power utility (eq. 14) and household income proportional to aggregate income (eq. 15), the interest sensitivity of the aggregate Euler equation equals that of the RA (Proposition 2). Source: Werning (2015), NBER WP 21448, sections 3.1-3.4. Implemented as a Huggett (1993) economy.

**Cyclicality of income risk.** y_it = Y_t e_it (Y_t/Y)^{ζ_i}, with e a Markov chain independent of Y (Rouwenhorst 3 states, persistence 0.9, log sd 0.3).
- Proportional income: ζ = 0.
- Countercyclical risk: ζ = (1, −0.185, −0.185). Σ_i π_i e_i ζ_i = 0 preserves aggregate income at the steady state **only to first order**, which suffices for the linear SSJ test. For nonlinear paths, normalise by Σ π_i e_i (Y/Y_ss)^{ζ_i} or report the aggregate-income deviation separately. Actual deviation (v0.2): 9.3e-8 at Y = 1 ± 0.001, 9.2e-6 to 9.3e-6 at ± 0.01, 2.2e-4 to 2.4e-4 at ± 0.05 (O(dY²)).

Zero-liquidity closed form: linearising the Euler equation of the marginal household h = argmax_e E[(e'/e)^{−1/σ} | e] gives (1+ζ_h) ŷ_t = (1+ζ̄) ŷ_{t+1} − σ dr_{t+1}/(1+r), where ζ̄ is a marginal-utility-weighted average. In the countercyclical case δ = 1.0174 and the elasticity σ/(1+ζ_h) = 0.614. If the argmax over marginal saving states comes close to a tie, differentiability needs a separate check (in this calibration the top state is the unique marginal saver).

**Exact setting.** σ = 0.5, β set so that r* = 1% (β = 0.86997, v0.1 result). Borrowing limit −b, zero net bond supply, grid −b + 20 b u² (400 points), difference step h = 1e-3 b. dr_t = −0.01 × 0.8^t, T = 150.

**Results.**

| b | 1e-1 | 3e-2 | 1e-2 | 3e-3 | 1e-3 | 3e-4 | 1e-4 |
|---|---|---|---|---|---|---|---|
| proportional income: error vs RANK | 6.0e-2 | 6.6e-3 | 1.29e-2 | 1.23e-2 | 4.23e-3 | 1.29e-3 | **4.30e-4** |
| countercyclical: error vs zero-liquidity solution | 0.100 | 2.3e-2 | 9.7e-3 | 1.70e-2 | 5.89e-3 | 1.79e-3 | **5.99e-4** |
| countercyclical: deviation from RANK | 0.452 | 0.350 | 0.333 | 0.297 | 0.312 | 0.317 | **0.319** |
| (r(b) − r*)/b | 0.58 | 1.15 | 1.85 | 2.07 | 2.081 | 2.087 | 2.088 |

Order in the asymptotic range (b ≤ 1e-3): 0.993 (both cases). With grids of 200, 400 and 800 points the errors agree to four digits. The numerical Jacobian is not differentiated directly at b = 0 because of the kink; the result is verified only through the same b → 0 limit as the refinement.

---

## 7. T4: Aiyagari steady state (no-idiosyncratic-risk limit)

**Scope.** (i) The no-idiosyncratic-risk limit: as income risk → 0, r → 1/β − 1. The code has no completion of insurance markets, so this is not a "complete-markets limit". (ii) A **regression check** that, in this calibration, r < 1/β − 1 with risk and r falls with the size and persistence of risk. Not used as a general theorem.

**Mechanical and non-mechanical checks.** `ssj_aiyagari_r` finds the root in (−δ + 1e-3, 1/β − 1 − 1e-10), so r < 1/β − 1 holds by construction and is not counted as independent accuracy evidence (kept only as a sanity check). Non-mechanical checks in v0.2:

| Check | Max over all T4 points |
|---|---|
| capital-market residual \|A − K\|/K | 5.2e-10 |
| Euler error (unconstrained states, distribution-weighted mean) | 1.1e-7 |
| Euler error (unconstrained states, max over points with mass) | 1.6e-5 |
| KKT violation (constrained states: only where u'(c) < βR E u'(c')) | 0 |
| mass on the top 5 grid points / mass with off-grid policy | 4.4e-29 / 3.1e-30 |

The Euler error is unit-free, c_impl/c − 1, with c(e', ·) linearly interpolated on the asset grid. A positive Euler gap in constrained states is not counted as numerical error (`hank_limits/diagnostics.euler_kkt_residuals`).

**Asset-limit and resolution sensitivity (sd = 0.4, persistence 0.9, the widest wealth distribution).** r (%):

| a_max \ points | 500 | 1000 | 2000 |
|---|---|---|---|
| 100 | 2.8818 | 2.8822 | 2.8823 |
| 200 (default) | 2.8813 | 2.8821 | 2.8823 |
| 400 | 2.8801 | 2.8818 | 2.8822 |

**Source.** Aiyagari (1994), QJE 109(3), 659-684.

**Exact setting.** β = 0.96, α = 0.36, δ = 0.08, σ (EIS) = 0.5 (risk aversion 2), borrowing limit 0, labour 1, Rouwenhorst 7 states, 1000 points on [0, 200] (u² spacing).

**Results.**

| eps | 0.4 | 0.2 | 0.1 | 0.05 | 0.025 |
|---|---|---|---|---|---|
| r (%) | 2.8821 | 3.8119 | 4.0742 | 4.1429 | 4.1605 |
| 1/β − 1 − r | 1.28e-2 | 3.55e-3 | 9.24e-4 | 2.38e-4 | **6.17e-5** |
| gap/eps² | 0.080 | 0.089 | 0.092 | 0.095 | 0.099 |

r by standard deviation (persistence 0.9): 4.074, 3.812, 3.404, 2.882%. r by persistence (sd 0.4): 3.956, 3.792, 3.475, 2.882%. gap/eps² is still rising slightly, partly a grid effect (0.096 at 2000 points).

**Comparison with Aiyagari (1994) Table II (v0.2, `scripts/aiyagari_table2_v0_2.py`).** In the original Table II, σ is described as the coefficient of variation of earnings, and the actual moments of the 7-state Markov approximation differ from the inputs (original Table I, footnote 38; checked by the reviewing model). Table II also has columns for risk aversion μ = 1, 3, 5, none of which equals the suite's EIS 0.5 (μ = 2). So two cells were computed with μ matched, varying the chain (Tauchen 1986 with width 3 sd; Rouwenhorst) and the reading of σ (log sd; coefficient of variation → log sd = √log(1+σ²)). β = 0.96, α = 0.36, δ = 0.08, borrowing limit 0, mean labour normalised to 1 (with a zero borrowing limit the household problem is homogeneous in the labour level, so the normalisation does not change r), 2000 points, upper limit 400.

| Cell (σ, ρ, μ) | Original* r / saving rate (%) | Tauchen, log sd | Tauchen, CV | Rouwenhorst, log sd | Rouwenhorst, CV |
|---|---|---|---|---|---|
| 0.2, 0.6, 3 | 3.8767 / 24.25 | **3.8782** / 24.25 | 3.8842 | 3.9127 | 3.9179 |
| 0.4, 0.9, 5 | −0.3456 / 37.63 | −0.0856 / 36.39 | 0.1041 | 0.7265 | 0.9007 |

\* The original values were read from a scan of the original (printed p. 678) by the reviewing model (GPT 6.1 Sol). We could not read the scan's text and have not checked them directly. The saving rate s = δα/(r+δ) is a function of r, so a matching saving rate at the same r is not an independent second check.

- **The first cell** matches to within 0.0016 percentage points with the Tauchen chain (3.8782 vs 3.8767).
- **The second cell is not reproduced and the cause is unresolved** (−0.0856 vs −0.3456, about 0.26 percentage points). The original Table I (printed p. 676) reports the **level income CV and autocorrelation** of the chain, and the values of the Tauchen 3-sd chain already match the original two-digit rounded values.

  | Cell (σ, ρ) | Original Table I: CV / corr (read by the reviewing model) | Tauchen 3-sd chain: level CV / level corr |
  |---|---|---|
  | 0.2, 0.6 | 0.21 / 0.59 | 0.2142 / 0.5935 |
  | 0.4, 0.9 | 0.49 / 0.89 | 0.4885 / 0.8913 |

  We therefore withdraw the v0.2 interpretation that the gap "is explained by differences in the actual risk of the chain". Equal moments do not guarantee equal transition matrices, so a chain difference is not ruled out either. Changing the asset grid does not close the gap (r: −0.0862 at 1000 points and limit 400, −0.0856 at 2000/400, −0.0855 at 4000/400, −0.0858% at 2000/800; top-of-grid and off-grid mass 0). Candidate causes (computation from simulated series in the original, the transition matrix, other implementation differences) are all **[unverified]**.
- For reference, Tauchen widths of 2, 2.5 and 3 give chain log sds of 0.399, 0.440 and 0.468 and r of 0.915, 0.334 and −0.086%: the second cell is very sensitive to the chain. We did not tune the chain to match the original value.
- Public Econ-ARK replications of Aiyagari (HARK, DolARK) exist; we did not audit the accuracy of that code.

---

## 8. T5: intertemporal Keynesian cross consistency

**Scope.** A GE-IKC **consistency test** built from the same block Jacobians. It is distinct from an independent certification of the economic model's accuracy. It is practical in that it caught F1 and F2.

**Statement.** In a model whose only HA block is the household, GE output is dY = (I − M)^{-1} [dG + J^{C,Tax} dTax + J^{C,r} dr] (M = J^{C,Y}; dTax and dr are the GE paths). Source: Auclert, Rognlie, Straub (2024), JPE 132(12), 4068-4121 (NBER WP 25020).

**Exact setting.** Household: borrowing limit 0, income (Y − Tax) e, Rouwenhorst 7 states (persistence 0.92, sd 0.6), 500 points, σ 0.5. r = 0.02, B = 1 = A (β = 0.9105, constrained mass 19.8%), G = 0.2, Tax_t = Tax + 0.1 (B_{t−1} − B). A fixed-real-rate model and a model with the real-rate rule r_t = r + 0.25 (Y_t − 1) + shock. Unknowns Y, B. Four shocks, T = 100, 200, 300. Three solves: (a) SSJ with matched options, (b) the same algorithm plus iterative refinement, (c) mismatched options (F1).

**Results (relative error, T = 100 / 200 / 300).**

| Shock | (a) SSJ default | (b) refined | (c) mismatched options |
|---|---|---|---|
| G, fixed rate | 1.6e-15 / 1.5e-14 / 9.5e-14 | 5.3e-16 / 9.7e-15 / 1.9e-15 | 1.6e-15 / 1.5e-14 / 9.5e-14 |
| rate, fixed rate | 3.9e-15 / 4.8e-14 / 3.1e-13 | 8.5e-16 / 1.7e-14 / 1.4e-13 | **9.9e-5** / 9.9e-5 / 9.9e-5 |
| G, rate rule | 5.0e-13 / 2.2e-9 / **5.9e-7** | 2.0e-15 / 3.2e-15 / 3.0e-15 | 5.0e-13 / 2.2e-9 / 5.9e-7 |
| rate shock, rate rule | 4.8e-13 / 2.3e-9 / **3.3e-7** | 3.1e-13 / 1.3e-12 / 9.8e-12 | **9.0e-5** / 9.0e-5 / 9.0e-5 |

- Tolerance 1e-10. The rate-rule cases of (a) at T ≥ 200 and the rate shocks of (c) are `xfail(strict=True)`. **These strict expected failures apply only to the regression test pinned to SSJ 1.0.0.** If SSJ is fixed they become XPASS and signal an update. Do not carry a rule that other solvers "must fail" into a general specification.
- Causes and minimal reproductions of F1 and F2: section 14.
- Walras-law residual: asset market at most 9e-10 (refined solve, T = 300).

---

## 9. T6: HANK − RANK differencing identity and decomposition

**Statement.** With A = I − M_HA and B = I − M_RA, A^{-1} − B^{-1} = A^{-1}(B − A)B^{-1} (finite-T resolvent identity). For the same steady state and the same policy paths (dG, dTax, dr), dY_HA − dY_RA = (I−M_HA)^{-1} [ (M_HA − M_RA)(dY_RA − dTax) + (J^r_HA − J^r_RA) dr ].

**Scope.** Verifies the finite-T resolvent identity and the differencing decomposition along the same paths. It does not verify a unique group-wise GE attribution at an infinite horizon (key finding 5 in section 1).

**Exact setting.** The fixed-real-rate fiscal HANK of T5. The RA is the closed form with the same r, σ, C, A (= B) (β_RA = 1/(1+r)). Shocks: debt-financed G (0.01 × 0.8^t) and rate −0.01 × 0.8^t. T = 150, 300.

**Results.**
- Identity error: 3.8e-15 at T = 150, 2.2e-14 at T = 300. Condition numbers: cond(I − M_RA) 31 / 1131, cond(I − M_HA) 294 / 5933.
- Decomposition error: G shock 9.8e-12 / 1.9e-10, rate shock 1.1e-11 / 2.2e-10.
- The RANK multiplier is 1 by Ricardian equivalence; the 1 + 6.4e-6 at T = 150 comes from the truncated tail of the debt path (0.92^T) (1 + 2.5e-11 at T = 300). HANK impact multiplier 1.979, present-value cumulative multiplier 3.383.

| t | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| dY_HA − dY_RA | 0.00979 | 0.00932 | 0.00880 | 0.00825 | 0.00770 | 0.00716 |
| first round | 0.00316 | 0.00274 | 0.00231 | 0.00190 | 0.00151 | 0.00115 |
| GE amplification | 0.00664 | 0.00658 | 0.00648 | 0.00635 | 0.00619 | 0.00601 |

- First-round gaps by wealth quartile (t = 0): Q1 1.05e-3, Q2 8.7e-4, Q3 6.5e-4, Q4 5.9e-4. Their sum is within 3.2e-6 of the Jacobian total and agrees between T = 150 and T = 300 to within 1e-8.
- Present value of the first-round gaps: −2.09e-3, −8.5e-4, 3.5e-4, 2.58e-3 for Q1-Q4 (sum zero). Applying (I−M_HA)^{-1} group by group gives values that grow with T (lowest quartile, impact, T = 150: −2.0e-3; T = 300: −5.8e-3; T = 450: −8.7e-3). So by group we report only the first round, and GE amplification only in total.

---

## 10. Plugging in your own HANK code

The check functions (`hank_limits/checks.py`) do not depend on SSJ; they take provider functions only. The conventions are in `hank_limits/interface.py`. Jacobians are T × T matrices (row = outcome date t, column = shock date s) in levels, dr in absolute units (0.01 = 1 percentage point), and **the timing convention must match SSJ's.**

| Check | Provider signature | Returns |
|---|---|---|
| `check_rank_limit` | `f(eps, T)` | `J={'C': {'Y', 'r'}}`, `ss={'C','A','r','eis'}` |
| `check_tank_limit` | `f(eps, T, dr, chis)` | `M_H, M_S`, `dY={chi: path}`, `lam`, `ss` |
| `check_werning` | `f(b, T, dr, cyclical)` | `dY, r_ss, eis, e_grid, Pi, zeta`, `diag` |
| `check_aiyagari` | `f(sd, rho)` | `r, rho_disc`, optional `diag` (residuals) |
| `check_ikc`, `check_differencing` | `f(T, shocks)` | `J`, `ss`, `irfs` |

Recommendations: use T1 and T4 as single-type household-block unit tests, T2 in a form generalised to any number of types, and T5 and T6 as regression tests of your own GE solver. The numerical tolerances are values for this calibration, so reset them for your own calibration. Differences in timing convention or units usually show up at once in the eps = 0 exactness test (T1, at the 1e-9 level).

---

## 11. pytest output (actual)

The full output is in `out/pytest_output_v0_2_1.txt`.

```
tests/test_t5_ikc.py::test_ikc_identity_ssj_solver_rule_long_horizon XFAIL
tests/test_t5_ikc.py::test_ikc_identity_with_mismatched_options XFAIL
tests/test_t5_ikc.py::test_walras_law PASSED
...
tests/test_diagnostics_v0_2.py::test_inverse_iteration_report[audit_t4] PASSED
tests/test_diagnostics_v0_2.py::test_inverse_iteration_report[audit_t1_mid] PASSED
tests/test_diagnostics_v0_2.py::test_ssj_forward_iteration_converges_to_it[audit_t4] PASSED
tests/test_diagnostics_v0_2.py::test_ssj_forward_iteration_converges_to_it[audit_t1_mid] PASSED
tests/test_diagnostics_v0_2.py::test_ssj_stopping_rule_is_not_an_accuracy_bound PASSED
...
=========================== short test summary info ============================
XFAIL tests/test_t5_ikc.py::test_ikc_identity_ssj_solver_rule_long_horizon - Finding F2 (SSJ 1.0.0): partial-pivoting LU on H_U packed in variable-block order shows element growth ~1e7 at T=300 (cond ~2e4); date-interleaved ordering or iterative refinement removes it.
XFAIL tests/test_t5_ikc.py::test_ikc_identity_with_mismatched_options - Finding F1 (SSJ 1.0.0): HetBlock linear impulses recompute Jacobians despite supplied Js (own options, default one-sided h=1e-4); with two-sided Js the r-shock IKC holds only to ~1e-4.
=================== 41 passed, 2 xfailed in 65.36s (0:01:05) ===================
```

---

## 12. FastHetBlock distribution audit (v0.2)

`scripts/diagnostics_v0_2.py` part A, `out/diagnostics_v0_2.json`. Baseline is the default setting (τ = 1e-10, at most 20 sweeps, 5 polishing steps).

| Model (mixing speed) | inverse iterations | negative mass | L1 residual | change in mean assets for τ 1e-12 to 1e-6 | SSJ forward iteration, tol 1e-10: iterations / mean-asset difference | tol 1e-12 |
|---|---|---|---|---|---|---|
| T1 PIH eps = 0.1 (slow) | 3 | 0 | 1.2e-16 | ≤ 1.1e-14 | 37,910 / **1.7e-3** | 54,590 / 1.7e-5 |
| T1 PIH eps = 0.0125 (very slow) | 4 | 0 | 1.0e-16 | ≤ 1.1e-12 | does not stop within 600,000 | — |
| T4 sd = 0.2 (fast) | 3 | 0 | 1.8e-16 | ≤ 1.1e-15 | 3,010 / 1.1e-6 | 3,730 / 1.1e-8 |
| T5 fiscal HANK (fast) | 3 | 0 | 2.0e-16 | ≤ 2.2e-16 | 210 / 4.9e-9 | 250 / 7.5e-11 |

- With a single inverse iteration and no polishing, mean assets at T1 eps = 0.0125 differ by 6.8e-4; with 2 iterations by 2.0e-8; with 3 by 5.4e-13. The default iterates until the stopping rule (1e-15) is met.
- SSJ's forward iteration converges to the inverse-iteration solution as the tolerance is tightened. SSJ's stopping rule (difference between successive iterates) is not an accuracy bound under slow mixing, so care is needed when β is calibrated to an asset target in near-RA calibrations. This is an observation about usage, not a report of an SSJ defect.
- A replica loop that follows SSJ's stopping rule and counts iterations was checked to agree with SSJ entry by entry.

---

## 13. Related material, not part of the test suite

- `scripts/a1_markov_checks_v0_2.py`: first checks for a design note on adjoint-based error estimation (A1): the identity with the transfer term on a small chain, the order (2.00) of the first-order prediction remainder under conservative perturbations of a real household chain, and the identity and approximate-adjoint effectivity for coarse → fine stationary-target errors. **It does not implement or validate an adaptive-grid algorithm.** Results are in a separate design note.
- If you plan to compare with Dynare's heterogeneous-household features, use **Dynare 7.1 or later**. The 7.1 release (2026-05-28) fixed an error in 7.0 where complementarity conditions in `model(heterogeneity=…)` blocks were interpreted with the opposite sign, and an error where aggregate variables appearing only in heterogeneous blocks were not recognised. Pinning to 7.0 could contaminate borrowing-constraint tests.

---

## 14. Minimal reproductions of F1 and F2

**F1** (`scripts/repro_F1_v0_2.py`, SSJ code only): one block `sequence_jacobian.hetblocks.hh_sim.hh`, T = 300.

| Comparison | Max relative difference |
|---|---:|
| `impulse_linear(..., Js={'hh': J2})` vs `J2 @ dr` (J2 two-sided) | 9.3e-5 |
| the same call vs `J1 @ dr` (J1 default one-sided) | 0 |
| call with `Js={'hh': 2 × J2}` vs call with `Js={'hh': J2}` | **0** |
| `impulse_linear(..., Js={'hh': J2}, twosided=True)` vs `J2 @ dr` | 0 |

**F2** (`scripts/repro_F2_v0_2.py`, SSJ code only, stationary distribution from SSJ's own forward iteration). `np.linalg.solve` follows the same path as `Block.solve_impulse_linear`.

| T | cond₂(H_U) | LU element growth | normalised backward error η (SSJ path) | η after 2 refinement steps |
|---:|---:|---:|---:|---:|
| 100 | 3.6e2 | 1.9e1 | 3.7e-14 | 5.7e-17 |
| 200 | 2.9e3 | 3.7e4 | 1.6e-10 | 7.6e-17 |
| 300 | 2.1e4 | 1.1e7 | 3.1e-8 | 1.0e-16 |

η = ‖H_U x − b‖_∞ / (‖H_U‖_∞ ‖x‖_∞ + ‖b‖_∞). At T = 300, reordering by date gives growth 1.02 and η 4.3e-16; reversing the block order gives growth 1.00 and η 9.3e-17. H_U and b are in `out/F2_HU_T300_v0_2.npz`.

Reports to the SSJ maintainers (2026-10-20): F1 [shade-econ/sequence-jacobian#39](https://github.com/shade-econ/sequence-jacobian/issues/39), F2 [shade-econ/sequence-jacobian#40](https://github.com/shade-econ/sequence-jacobian/issues/40). Titles describe the observed behaviour ("HetBlock linear impulses recompute Jacobians despite supplied Js").

---

## 15. Figures and tables

- Figures: `out/t1_rank_limit.png` … `out/t6_decomposition.png` (log axes, limit to the right).
- Machine-readable tables: `out/tables_v0_2/suite_results_long_v0_2.csv` (columns: test, limit_param, limit_value, n_a, a_max, fd_h, fd_twosided, window, horizon_T, solver, metric, value), `out/convergence_v0_2.json`, `out/diagnostics_v0_2.json`, `out/aiyagari_table2_v0_2.json`, `out/repro_F1_v0_2.json`, `out/repro_F2_v0_2.json`.

## 16. References

- Aiyagari, S. R. (1994). Uninsured Idiosyncratic Risk and Aggregate Saving. QJE 109(3), 659-684. https://ideas.repec.org/a/oup/qjecon/v109y1994i3p659-684..html
- Auclert, A., Bardóczy, B., Rognlie, M., Straub, L. (2021). Using the Sequence-Space Jacobian to Solve and Estimate Heterogeneous-Agent Models. Econometrica 89(5), 2375-2408. https://doi.org/10.3982/ECTA17434 (Definition 1 and eq. (23) state the expectation vectors explicitly, and the paper itself reports Jacobian accuracy checks against direct computation and automatic differentiation. This suite complements those checks with economic limits.)
- Auclert, A., Rognlie, M., Straub, L. (2024). The Intertemporal Keynesian Cross. JPE 132(12), 4068-4121. https://doi.org/10.1086/732531 (NBER WP 25020)
- Bayer, C., Luetticke, R., Weiss, M., Winkelmann, Y. (2026). An Endogenous Gridpoint Method for Distributional Dynamics. JME 158, 103895. https://doi.org/10.1016/j.jmoneco.2026.103895 (a distribution update different from the lottery method; not compared here)
- Bilbiie, F. O. (2008). JET 140(1), 162-196. https://ideas.repec.org/a/eee/jetheo/v140y2008i1p162-196.html
- Bilbiie, F. O. (2020). The New Keynesian cross. JME 114, 90-108. https://doi.org/10.1016/j.jmoneco.2019.03.003
- Dynare Team (2026). Release of Dynare 7.1. https://forum.dynare.org/t/release-of-dynare-7-1/28376
- Huggett, M. (1993). JEDC 17(5-6), 953-969. https://ideas.repec.org/a/eee/dyncon/v17y1993i5-6p953-969.html
- Santos, M. S. (2000). Accuracy of Numerical Solutions Using the Euler Equation Residuals. Econometrica 68(6), 1377-1402. https://doi.org/10.1111/1468-0262.00165 (background for the T4 Euler residuals; not cited as a guarantee for constrained households, distributions or GE targets)
- Tauchen, G. (1986). Finite State Markov-Chain Approximations to Univariate and Vector Autoregressions. Economics Letters 20(2), 177-181.
- Werning, I. (2015). Incomplete Markets and Aggregate Demand. NBER WP 21448. https://www.nber.org/papers/w21448
- Young, E. R. (2010). Solving the Incomplete Markets Model with Aggregate Uncertainty Using the Krusell-Smith Algorithm and Non-Stochastic Simulations. JEDC 34(1), 36-41. https://doi.org/10.1016/j.jedc.2008.11.010
- Econ-ARK, Aiyagari (1994) replication. https://econ-ark.org/materials/aiyagari1994qje/
- sequence-jacobian: https://github.com/shade-econ/sequence-jacobian , https://pypi.org/project/sequence-jacobian/1.0.0/

## 17. AI usage disclosure

Generative AI tools were used extensively in this project. Claude (Anthropic; Claude Opus 5.5, October 2026) wrote most of the code, tests, numerical experiments and documentation, and derived and checked the closed-form benchmarks, working from the author's instructions. GPT 6.1 Sol (OpenAI) reviewed the code and documents in two rounds; the changes made in response are listed at the top of this README. The author chose the research questions and the tests to include, directed and reviewed the work, decided which findings to report and how to word them, and takes full responsibility for the accuracy, originality and licensing of all submitted material.

## 18. Licence

MIT. See `LICENSE`.
