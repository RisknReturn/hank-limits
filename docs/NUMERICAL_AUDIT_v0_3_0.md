# Numerical audit before v0.3.0 (2026-10-28)

An independent audit (GPT, run on the author's computer) of the numerical safeguards in `FastHetBlock` (`hank_limits/core.py`) and of the calibration paths that use them, carried out before the v0.3.0 release. The judgement criteria were written and hashed before any numerical output was read (SHA-256 `a5699008e2a6637cd9307f57174fc5b4bdb9af5a7d685a6c6aa43cfb8f3bb2d5`, fixed 2026-10-28T03:37:32Z). The audit worked on an archive of commit `35ea89c`; the full audit record (criteria, 8,381 diagnostic rows, sensitivity copies and diffs, file hashes) is kept by the author and is not part of this repository.

## What was audited

| Lines of `core.py` (at `35ea89c`) | Treatment |
|---|---|
| 138–140 | Policies clipped to the ends of the asset grid; lottery weights clipped to [0, 1] |
| 178–183 | Negative mass of the inverse iterate recorded, then the absolute value normalised |
| 189–202 | Residual of the clipped chain recorded after 5 polishing steps |
| 80–101 | SSJ's own unclipped forward residual and off-grid policy mass recorded separately |

Variants run on copies only: (a) zero truncation instead of the absolute value, (b) signed normalisation, (c) no weight clip, and (c2) neither policy nor weight clip (an extrapolating, non-stochastic transition, used as a probe and not as an alternative solver).

## Findings

1. **Absolute value and weight clip: no effect.** In 967 inverse-iteration solves the iterate never had negative mass before the absolute value. Variants (a), (b) and (c) gave 114 passed, 4 xfailed, like the baseline, and left 2,064 summary values unchanged (the only differences, in (b), were 8 T5 auxiliary residuals of at most 1.7e-16, which the uninstrumented original showed too). A further 66,860 values from extra settings (T1/T2 at 9,600 points, T3 at 800 points and on the goods-market target, T4 at 2,000 points and upper limit 100) were also unchanged.
2. **Final distributions: boundary mass is negligible.** Over 68 final distributions the mass above the grid was at most 6.1e-28, below it 0.
3. **Root-search trial states are different.** A T3 trial state had 0.58 of its mass at clipped policies and an SSJ forward residual of 19.6 while its clipped-chain residual was 5.6e-17. Trial-state and final-state diagnostics, and clipped-chain and SSJ residuals, must therefore be kept apart.
4. **Small upper-limit study (T4, upper limit 100, 1,000 points).** Clipped mass 7.6e-13: below the test's 1e-12 mass threshold, above the audit's noise screen.
5. **PIH calibration path dependence (the reason for the code change).** The calibration of β (T1, T2 savers) warm-started each root-search trial from the previous trial's backward-iteration solution. SSJ stops the backward iteration on the change between iterates, which is not an error bound, so the calibrated state depended on the search path. With the clipped transition unchanged, moving the root bracket changed T2's displayed IRF error from 1.05e-3 to 1.03e-3 and its fitted order from 1.79 to 1.80. All tests still passed, so the suite did not detect this. Tightening the backward tolerance to 1e-13 or 1e-14 did not converge within 100,000 iterations.

The overall verdict was "minimal code change needed before release", for finding 5, not for the absolute value or the clips.

## Change in v0.3.0

`models.pih_household` clears the cached backward seed (`_va_cache`) before every calibration trial and before the final recomputation, so each solve starts from the same cold seed. In the audit's copy the two brackets then gave the same T2 values (1.04e-3, order 1.79) and the full suite gave 114 passed, 4 xfailed.

Effect on the outputs of this repository, regenerated with the same package versions as recorded in `out/convergence_v0_2.json` (Python 3.13.16, sequence-jacobian 1.0.0, numpy 2.5.3, scipy 1.18.1, numba 0.68.0):

| Quantity | Before | v0.3.0 |
|---|---|---|
| T2 IRF error vs (1−λ)/(1−λχ) × RANK, eps = 0.0125 | 1.0500e-3 | 1.0415e-3 |
| T2 IRF error, eps = 0.025 | 3.456e-3 | 3.455e-3 |
| T2 grid study, eps = 0.0125: 2,400 / 4,800 / 9,600 points | 1.498e-3 / 1.050e-3 / 8.698e-4 | 1.494e-3 / 1.041e-3 / 8.647e-4 |
| T2 IRF fitted order (χ = 1.5) | 1.7908 | 1.7932 |
| T1 asset target, max relative residual (V4 `spec_A`) | 5.7e-8 | 7.1e-11 |
| T2 asset / goods market, max residual (V2) | 2.2e-7 / 6.6e-9 | 2.0e-8 / 6.1e-10 |
| T1 errors, κ, orders, windows, difference steps | — | changed by less than 1e-7 relative |
| T3–T6 | — | unchanged |

The smaller calibration residuals indicate that the warm start had also left the calibrated state further from its target than the root-search tolerance implies.

`out/convergence_v0_2.json` had last been generated on 2026-10-19 with the v0.2 code. Regenerated now, it also carries the `v0_3` block that the v0.3 checks attach to each result (T1–T5); the rerun of `35ea89c` described next produces the same blocks, so they are not an effect of the change.

Before the change, a rerun of `35ea89c` in the environment above reproduced the committed T1–T4 and T6 numbers exactly and the T5 residuals to within 3e-9 (relative), so the differences in the table are due to the change.

## Not changed, and limits

- The Werning (T3), Aiyagari (T4) and fiscal HANK (T5) calibrations still warm-start their root searches. The audit re-solved T3 on a narrow bracket and its displayed values were unchanged; T4 and T5 were not bracket-tested.
- Figures in `out/*.png` were not redrawn (changes below plot resolution).
- The audit also recommended a regression test that compares two calibration brackets for T2; it is not part of v0.3.0. The comment in `models.werning_model` that said trial savers stay inside the grid was corrected.
- The mutation tests (`docs/MUTATION_v0_1.md`, `docs/MUTATION_v0_2.md`, `docs/CHECKS_v0_3_0_AMENDMENT.md`) were run before this change and were not re-run.
- Deferred to a later version: a general convergence guarantee for warm starts, state-wise stopping rules and error bounds, a common standard for provider diagnostics, wider grid and initial-value stability checks, and a separate solver for non-stochastic extrapolated transitions.
- None of these checks is a general error bound.
