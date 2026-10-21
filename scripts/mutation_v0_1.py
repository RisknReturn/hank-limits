"""Mutation testing v0.1: inject small, plausible implementation errors into the model side and
record which tests of the suite catch them (docs/MUTATION_PLAN_v0_1.md, docs/MUTATION_v0_1.md).

PRE-REGISTRATION. The mutant list and the predictions below (and in `MUTANTS`) were committed
before any mutant was run. Two things were done before that commit: `--check`, which applies
each patch to a copy of the package and byte-compiles it and runs no tests, and a smoke test of
the harness on the unmutated copy only (two baseline runs; plan, section 2, rule 4). Changes made
after seeing results are logged in docs/MUTATION_PLAN_v0_1.md, section "Amendments".

Rules. Mutations touch only the model side: copies of `hank_limits/models.py` (providers and
their model blocks) and `hank_limits/core.py` (the household block and stationary distribution
the providers build on). The checker side is never changed: `checks.py`, `ra.py`, `metrics.py`,
`diagnostics.py`, `suite.py` (canonical settings) and `tests/` (tolerances, references).

Method. Each run copies `hank_limits/`, `tests/` and `pytest.ini` to a temporary directory,
applies one mutant's patch there and runs the full pytest suite on the copy, with this file
loaded as a pytest plugin (`-p mutation_v0_1`). The plugin records every test outcome and a
summary of each check's result (errors along the limit sequence, rates, residuals), and stops a
fixture or test that runs longer than `--fixture-timeout`. The unmutated copy is run twice
(B0, B1) as the baseline; a test "catches" a mutant when its outcome differs from the baseline
outcome (passed -> failed/error, or xfailed -> failed/xpassed). A test group (T1..T6, and D =
tests/test_diagnostics_v0_2.py) catches a mutant when at least one of its tests does.

Mutants (category; file:function; change) -> PREDICTED catching groups, confidence H/M/L.
Groups not listed are predicted to miss. '?' = no firm prediction.

  M01 timing  core:hh             Euler discounts with r_t instead of r_{t+1}  -> T1 T2 T3 (H)
  M02 timing  models:fiscal       debt service at r_{t-1} (ex-ante) not r_t    -> T5 via Walras (H)
  M03 timing  models:fiscal       tax rule on B_t instead of B_{t-1}           -> none (H); T6 ?
  M04 timing  models:werning prov rate path shifted one period (ex-ante fix)   -> T3 (H)
  M05 timing  models:fiscal       budget uses B_t instead of B_{t-1}           -> T5 T6 (M)
  M06 units   models:tank prov    rate path divided by 100 (percent)           -> T2 (H)
  M07 units   models:aiyagari     quarterly depreciation 0.02 in annual model  -> none (M); T4 ?
  M08 units   core:rouwenhorst    income levels not normalised to mean 1       -> none (M)
  M09 jacob.  models:hh_jacobians default FD step 1e-2 instead of 1e-4         -> T5 T6 (H)
  M10 jacob.  models:tank prov    one-sided saver Jacobian, two-sided options  -> T2 (H)
  M11 jacob.  models:hh_jacobians Jacobians returned transposed                -> T1 T5 T6 (H)
  M12 jacob.  models:werning prov cached Js passed, options dropped (F1 type)  -> T3 (M)
  M13 distr.  core:transition     Young lottery weights swapped                -> T1 T2 T4 D (M); T3 T6 ?
  M14 distr.  models:cohort paths cohort distribution not renormalised         -> T6 (H)
  M15 distr.  models:werning prim Markov matrix transposed                     -> none (L)
  M16 distr.  core:hh             borrowing constraint at a_grid[1]            -> T4 (H)
  M17 distr.  models:aiyagari     asset grid top 20 instead of 200             -> T4 (H); D ?
  M18 equil.  models:clearing     goods market Y - C + G                       -> T5 T6 (H)
  M19 equil.  models:tank         goods market cleared with savers' C only     -> T2 (H)
  M20 equil.  models:fiscal       government budget without interest           -> T5 T6 (H)
  M21 equil.  models:fiscal_hank  steady-state taxes omit debt service r B     -> none (M)
  M22 equil.  models:tank prov    saver iMPC taken from the H-type block       -> T2 (H)
  M23 calib.  models:tank_model   zero-net-supply target without /(1-lam)      -> none (M)
  M24 calib.  models:aiyagari     borrowing limit -1 instead of 0              -> none (M)
  M25 calib.  models:fiscal_hank  beta target A = 4B (quarterly-GDP units)     -> T5 via Walras (H)
  C01 control models:aiyagari     asset grid 1500 instead of 1000 points       -> none (M)
  C02 control core:hh             Va = (1+r) / c^(1/eis) (same expression)     -> none (H)
  C03 control models:werning prov target goods market instead of bond market   -> none (M)
                                  (Walras-equivalent: A_ss = 0, dA_{-1} = 0)

Usage:
  python3 scripts/mutation_v0_1.py --check            # apply + compile every patch, run no tests
  python3 scripts/mutation_v0_1.py --table            # print the mutant table (markdown)
  python3 scripts/mutation_v0_1.py --jobs 4           # baselines + all mutants -> out/mutation_v0_1.{json,csv}
  python3 scripts/mutation_v0_1.py --only M01 C02     # a subset (still runs the baselines)
Every run is single-threaded (OMP_NUM_THREADS=1 etc.); `--jobs` runs that many in parallel.
About 75-90 s per run; 30 runs on 4 workers take roughly 10-15 minutes.
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import math
import os
import platform
import py_compile
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GROUPS = ('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'D')
TEST_FILE_GROUP = {'test_t1_rank_limit.py': 'T1', 'test_t2_tank_limit.py': 'T2',
                   'test_t3_werning.py': 'T3', 'test_t4_aiyagari.py': 'T4',
                   'test_t5_ikc.py': 'T5', 'test_t6_differencing.py': 'T6',
                   'test_diagnostics_v0_2.py': 'D'}
ALL = ' '.join(GROUPS)

# ---------------------------------------------------------------------------
# Mutants. `patches`: (file, old, new), `old` must occur exactly once in the file.
# `reach`: groups whose providers execute the changed code (static, from the code map in
# docs/MUTATION_PLAN_v0_1.md). `expect`: groups predicted to catch; `unsure`: no prediction.
# ---------------------------------------------------------------------------

MUTANTS = [
    dict(id='M01', category='timing', target='core.py:_make_hh_function',
         change='Euler equation discounts with the current rate r_t instead of r_{t+1}: '
                'uc = beta (1+r_t) E[c\'^(-1/eis)], Va = c^(-1/eis)',
         patches=[('core.py', '        uc_nextgrid = beta * Va_p',
                   '        uc_nextgrid = beta * (1 + r) * Va_p'),
                  ('core.py', '        Va = (1 + r) * c ** (-1 / eis)\n        return Va, a, c',
                   '        Va = c ** (-1 / eis)\n        return Va, a, c')],
         reach=ALL, expect='T1 T2 T3', conf='H',
         why='Steady state unchanged (r constant), so T4/D see nothing. dC/dr is shifted by one '
             'date: T1 (closed-form dC/dr, exact at eps=0), T2/T3 (fixed-rate GE IRF vs closed '
             'form). T5/T6 build their references from the same household Jacobians.'),
    dict(id='M02', category='timing', target='models.py:fiscal',
         change='government pays interest at r_{t-1} (ex-ante reading) instead of the ex-post r_t',
         patches=[('models.py', '    fiscal_res = B - (1 + r) * B(-1) - G + Tax',
                   '    fiscal_res = B - (1 + r(-1)) * B(-1) - G + Tax')],
         reach='T5 T6', expect='T5', conf='H',
         why='IKC uses the solver\'s own dTax path, so the IKC tests pass; households earn r_t on '
             'A_{t-1} while the government pays r_{t-1} on B_{t-1}, so the asset market misses '
             'by B (dr_t - dr_{t-1}) under rate shocks: T5 Walras test. T6 uses dr = 0 for its '
             'Ricardian check and identities that hold for any path.'),
    dict(id='M03', category='timing', target='models.py:fiscal',
         change='tax rule responds to current debt B_t instead of B_{t-1}',
         patches=[('models.py', '    Tax = Tax_ss + phi_B * (B(-1) - B_ss)',
                   '    Tax = Tax_ss + phi_B * (B - B_ss)')],
         reach='T5 T6', expect='', unsure='T6', conf='H',
         why='A different but internally consistent and stable fiscal rule; no test pins the rule. '
             'Only T6\'s Ricardian tail bound 5 (1 + r - phi_B)^T encodes the original decay rate '
             '(new rate (1+r)/(1+phi_B) = 0.927 vs 0.92), so T6 may fail at T = 300.'),
    dict(id='M04', category='timing', target='models.py:ssj_werning_provider',
         change='rate path shifted one period (ex-ante -> ex-post conversion applied to a path '
                'that is already ex-post)',
         patches=[('models.py', "    irf = m['model'].solve_impulse_linear(ss, ['Y'], ['asset_mkt'], {'r': dr}, Js={'hh': J},",
                   "    irf = m['model'].solve_impulse_linear(ss, ['Y'], ['asset_mkt'], {'r': np.r_[0.0, dr[:-1]]}, Js={'hh': J},")],
         reach='T3', expect='T3', conf='H',
         why='dY_t depends on sum_{s>t} dr_s; a one-period shift changes it by O(1) relative.'),
    dict(id='M05', category='timing', target='models.py:fiscal',
         change='government budget uses end-of-period debt B_t instead of B_{t-1}: '
                'B - (1+r) B - G + Tax',
         patches=[('models.py', '    fiscal_res = B - (1 + r) * B(-1) - G + Tax',
                   '    fiscal_res = B - (1 + r) * B - G + Tax')],
         reach='T5 T6', expect='T5 T6', conf='M',
         why='Steady state still satisfied; dynamics of B explosive (factor phi_B / r = 5). Asset '
             'market misses (T5 Walras) and the RA benchmark is no longer Ricardian (T6).'),
    dict(id='M06', category='units', target='models.py:ssj_tank_provider',
         change='rate path divided by 100 (read as percent although the interface uses decimals)',
         patches=[('models.py', "        irf = model.solve_impulse_linear(ss_chi, ['Y'], ['goods_mkt'], {'r': dr}, Js=Js,",
                   "        irf = model.solve_impulse_linear(ss_chi, ['Y'], ['goods_mkt'], {'r': dr / 100}, Js=Js,")],
         reach='T2', expect='T2', conf='H',
         why='GE IRF 100 times too small against the closed form.'),
    dict(id='M07', category='units', target='models.py:AIYAGARI_DEFAULTS',
         change='quarterly depreciation delta = 0.02 in the annual calibration',
         patches=[('models.py', 'AIYAGARI_DEFAULTS = dict(beta=0.96, alpha=0.36, delta=0.08, eis=0.5, n_e=7, n_a=1000,',
                   'AIYAGARI_DEFAULTS = dict(beta=0.96, alpha=0.36, delta=0.02, eis=0.5, n_e=7, n_a=1000,')],
         reach='T4 D', expect='', unsure='T4', conf='M',
         why='The no-risk limit r -> 1/beta - 1 does not depend on delta, and monotonicity in risk '
             'is generic. Capital is about 3.5 times larger, so the top-of-grid mass check (T4) '
             'might trigger for a reason unrelated to the error.'),
    dict(id='M08', category='units', target='core.py:rouwenhorst',
         change='income levels e = exp(s) with E[log e] = 0 instead of E[e] = 1 '
                '(mean-income normalisation omitted)',
         patches=[('core.py', '    return sj.grids.markov_rouwenhorst(rho, sd, n)',
                   '    e, pi, Pi = sj.grids.markov_rouwenhorst(rho, sd, n)\n'
                   '    return np.exp(np.log(e) - pi @ np.log(e)), pi, Pi')],
         reach=ALL, expect='', conf='M',
         why='Scales income levels by about exp(sd^2/2). References are built from the reported '
             'steady state (C, A) and are scale-free; the factor goes to 1 in the eps -> 0 limits. '
             'No test checks the steady-state goods-market residual.'),
    dict(id='M09', category='jacobian', target='models.py:household_jacobians',
         change='default finite-difference step 1e-2 instead of 1e-4',
         patches=[('models.py', 'def household_jacobians(block, ss, T, inputs, outputs, h=1e-4):',
                   'def household_jacobians(block, ss, T, inputs, outputs, h=1e-2):')],
         reach='T1 T5 T6', expect='T5 T6', conf='H',
         why='T1\'s provider passes h explicitly; T1\'s horizon test compares two Jacobians with the '
             'same step. T5/T6: Js (h=1e-2) and SSJ\'s re-differentiation (h=1e-4) disagree, so the '
             'GE solution no longer satisfies the IKC built from Js (F1 mechanism).'),
    dict(id='M10', category='jacobian', target='models.py:ssj_tank_provider',
         change='saver Jacobian one-sided while the GE options stay two-sided',
         patches=[('models.py', "    J_S = m['hh_S'].jacobian(ss, ['Y_S', 'r'], ['C_S', 'A_S'], T=T, twosided=True)",
                   "    J_S = m['hh_S'].jacobian(ss, ['Y_S', 'r'], ['C_S', 'A_S'], T=T, twosided=False)")],
         reach='T2', expect='T2', conf='H',
         why='H_U from the one-sided Js, propagation re-differentiated two-sided: goods-market '
             'residual ~1e-6 > 1e-10 (test_ge_solution_residual).'),
    dict(id='M11', category='jacobian', target='models.py:household_jacobians',
         change='Jacobians returned transposed (row = shock date, column = output date)',
         patches=[('models.py', '    return {o: {i: np.asarray(J[o][i]) for i in inputs} for o in outputs}',
                   '    return {o: {i: np.asarray(J[o][i]).T for i in inputs} for o in outputs}')],
         reach='T1 T5 T6', expect='T1 T5 T6', conf='H',
         why='M_RA has identical rows, its transpose does not (T1). The transposed Js enter the SSJ '
             'solve while impulses are re-differentiated (T5, T6).'),
    dict(id='M12', category='jacobian', target='models.py:ssj_werning_provider',
         change='cached Js (two-sided, h = 1e-3 b) passed without matching options: SSJ '
                're-differentiates one-sided with h = 1e-4 (F1 type)',
         patches=[('models.py', "                                          options=_het_opts(['hh'], h=1e-3 * b))",
                   "                                          options={})")],
         reach='T3', expect='T3', conf='M',
         why='At b = 1e-4 the default step moves income by about b, so propagation and H_U '
             'disagree at O(1) in the b-scaled problem.'),
    dict(id='M13', category='distribution', target='core.py:transition_matrix',
         change='Young lottery weights swapped between the lower and upper grid point',
         patches=[('core.py', '    lot = sp.csr_matrix((np.concatenate([pi.ravel(), 1 - pi.ravel()]),',
                   '    lot = sp.csr_matrix((np.concatenate([1 - pi.ravel(), pi.ravel()]),')],
         reach=ALL, expect='T1 T2 T4 D', unsure='T3 T6', conf='M',
         why='Only FastHetBlock\'s stationary distribution changes (SSJ\'s own operator is intact): '
             'D compares with SSJ forward iteration; mean assets and so beta, r change (T1, T2, T4). '
             'T5 recalibrates beta to A = B under the wrong distribution and stays consistent. '
             'Calibration loops may not terminate (timeout counts as caught by error).'),
    dict(id='M14', category='distribution', target='models.py:cohort_impc_paths',
         change='cohort distribution not divided by its mass',
         patches=[('models.py', "        ss_k.internals['hh']['Dbeg'] = Dk / mu",
                   "        ss_k.internals['hh']['Dbeg'] = Dk")],
         reach='T6', expect='T6', conf='H',
         why='Group paths scaled by mu_k: sum_k mu_k dC_k != M x (test_wealth_group_split_adds_up).'),
    dict(id='M15', category='distribution', target='models.py:werning_primitives',
         change='Markov matrix transposed (column-stochastic convention)',
         patches=[('models.py', "    e, pi, Pi = rouwenhorst(p['rho_e'], p['sd_e'], p['n_e'])\n    gam = 1 / p['eis']",
                   "    e, pi, Pi = rouwenhorst(p['rho_e'], p['sd_e'], p['n_e'])\n    Pi = Pi.T\n    gam = 1 / p['eis']")],
         reach='T3', expect='', conf='L',
         why='The zero-liquidity reference is built from the provider-reported Pi, and the '
             'marginal-saver Euler argument holds for any positive Pi used consistently; the '
             'transposed chain keeps unit spectral radius (uniform income-state marginal). '
             'Low confidence: Pi.T rows do not sum to one.'),
    dict(id='M16', category='distribution', target='core.py:_make_hh_function',
         change='borrowing constraint imposed at a_grid[1] instead of a_grid[0] (grid boundary)',
         patches=[('core.py', '        misc.setmin(a, a_grid[0])',
                   '        misc.setmin(a, a_grid[1])')],
         reach=ALL, expect='T4', conf='H',
         why='States held at a_grid[1] are classified as unconstrained by the Euler diagnostic and '
             'show large Euler gaps (T4). T1 has no constrained mass; T2 H types still have '
             'dC/dY = 1; T3 shift is ~1e-4 b; T5/T6 consistent.'),
    dict(id='M17', category='distribution', target='models.py:AIYAGARI_DEFAULTS',
         change='asset grid upper bound 20 instead of 200',
         patches=[('models.py', '                         a_max=200.0)',
                   '                         a_max=20.0)')],
         reach='T4 D', expect='T4', unsure='D', conf='H',
         why='Mass at the top of the grid (test_top_grid_mass_negligible). D compares clipped and '
             'unclipped forward operators and may also move.'),
    dict(id='M18', category='equilibrium', target='models.py:fiscal_clearing',
         change='goods-market clearing sign error: Y - C + G',
         patches=[('models.py', '    goods_mkt = Y - C - G', '    goods_mkt = Y - C + G')],
         reach='T5 T6', expect='T5 T6', conf='H',
         why='IKC with +dG fails (T5), decomposition and multipliers fail (T6).'),
    dict(id='M19', category='equilibrium', target='models.py:goods_clearing',
         change='TANK goods market cleared with savers\' consumption only (Y = C_S)',
         patches=[('models.py', '@sj.simple\ndef goods_clearing(Y, C):\n    goods_mkt = Y - C\n',
                   '@sj.simple\ndef goods_clearing(Y, C_S):\n    goods_mkt = Y - C_S\n')],
         reach='T2', expect='T2', conf='H',
         why='chi = 1 gives the same answer (Y_S = Y), chi = 1.5 does not: T2 IRF and amplification.'),
    dict(id='M20', category='equilibrium', target='models.py:fiscal',
         change='government budget constraint without interest: B - B(-1) - G + Tax',
         patches=[('models.py', '    fiscal_res = B - (1 + r) * B(-1) - G + Tax',
                   '    fiscal_res = B - B(-1) - G + Tax')],
         reach='T5 T6', expect='T5 T6', conf='H',
         why='IKC holds on the solver\'s paths; asset market misses by r dB (T5 Walras); taxes cover '
             'phi_B/(r+phi_B) = 83% of spending in PV, so the RA multiplier is not 1 (T6).'),
    dict(id='M21', category='equilibrium', target='models.py:fiscal_hank',
         change='steady-state taxes omit debt service: Tax_ss = G instead of G + r B',
         patches=[('models.py', "    Tax_ss = p['G'] + p['r'] * p['B']", "    Tax_ss = p['G']")],
         reach='T5 T6', expect='', conf='M',
         why='Steady-state government budget violated (fiscal_res_ss = -r B), but no test checks '
             'steady-state residuals; dynamics, A = B and Walras in deviations are unaffected.'),
    dict(id='M22', category='equilibrium', target='models.py:ssj_tank_provider',
         change='saver iMPC matrix taken from the H-type block (copy-paste: J_H instead of J_S)',
         patches=[('models.py', "    return dict(M_H=np.asarray(J_H['C_H']['Y_H']), M_S=np.asarray(J_S['C_S']['Y_S']),",
                   "    return dict(M_H=np.asarray(J_H['C_H']['Y_H']), M_S=np.asarray(J_H['C_H']['Y_H']),")],
         reach='T2', expect='T2', conf='H',
         why='M_S = I instead of M_RA.'),
    dict(id='M23', category='calibration', target='models.py:tank_model',
         change='zero-net-supply target for savers without the division by (1 - lam)',
         patches=[('models.py', '    A_S_target = -lam * A_H / (1 - lam)', '    A_S_target = -lam * A_H')],
         reach='T2', expect='', conf='M',
         why='A_H -> 0 as eps -> 0 (H types end at the zero limit), so both targets go to 0; the '
             'error vanishes in the limit the test looks at.'),
    dict(id='M24', category='calibration', target='models.py:ssj_aiyagari_r',
         change='borrowing limit -1 instead of 0',
         patches=[('models.py', "    a_grid = p['a_max'] * u ** 2",
                   "    a_grid = -1.0 + (p['a_max'] + 1.0) * u ** 2")],
         reach='T4 D', expect='', conf='M',
         why='The no-risk limit and the monotonicity regressions hold for any fixed borrowing '
             'limit; Euler/KKT diagnostics use the mutated grid consistently.'),
    dict(id='M25', category='calibration', target='models.py:fiscal_hank',
         change='beta calibrated to A = 4 B (debt target in quarterly-GDP units) while B = 1',
         patches=[('models.py', "    f = lambda beta: hh.steady_state({**calib, 'beta': beta}, **SS_OPTS)['A'] - p['B']",
                   "    f = lambda beta: hh.steady_state({**calib, 'beta': beta}, **SS_OPTS)['A'] - 4 * p['B']")],
         reach='T5 T6', expect='T5', conf='H',
         why='A_ss != B_ss: in deviations dA - dB picks up (A_ss - B_ss) dr under rate shocks (T5 '
             'Walras). T6 compares with an RA at the same (C, A) and is consistent.'),
    # controls: changes that should not change any outcome
    dict(id='C01', category='control', target='models.py:AIYAGARI_DEFAULTS',
         change='asset grid refined from 1000 to 1500 points',
         patches=[('models.py', 'AIYAGARI_DEFAULTS = dict(beta=0.96, alpha=0.36, delta=0.08, eis=0.5, n_e=7, n_a=1000,',
                   'AIYAGARI_DEFAULTS = dict(beta=0.96, alpha=0.36, delta=0.08, eis=0.5, n_e=7, n_a=1500,')],
         reach='T4 D', expect='', conf='M',
         why='Refinement within the range the README reports; thresholds are regression thresholds.'),
    dict(id='C02', category='control', target='core.py:_make_hh_function',
         change='Va = (1 + r) / c ** (1 / eis) instead of (1 + r) * c ** (-1 / eis)',
         patches=[('core.py', '        Va = (1 + r) * c ** (-1 / eis)\n        return Va, a, c',
                   '        Va = (1 + r) / c ** (1 / eis)\n        return Va, a, c')],
         reach=ALL, expect='', conf='H',
         why='Same expression; rounding-level differences only.'),
    dict(id='C03', category='control', target='models.py:ssj_werning_provider',
         change='GE target goods market (C = Y) instead of bond market (A = 0)',
         patches=[('models.py', "    irf = m['model'].solve_impulse_linear(ss, ['Y'], ['asset_mkt'], {'r': dr}, Js={'hh': J},",
                   "    irf = m['model'].solve_impulse_linear(ss, ['Y'], ['goods_mkt'], {'r': dr}, Js={'hh': J},")],
         reach='T3', expect='', conf='M',
         why='Walras: dC - dY = (1+r) dA_{t-1} - dA_t with A_ss = 0, so C = Y for all t with '
             'dA_{-1} = 0 implies dA = 0; equivalent up to O(h^2) and rounding.'),
]
MUTANT_IDS = [m['id'] for m in MUTANTS]


def _groups(s):
    return [g for g in s.split() if g]


# ---------------------------------------------------------------------------
# pytest plugin part (active only when MUTATION_V0_1_OUT is set)
# ---------------------------------------------------------------------------

_OUT = os.environ.get('MUTATION_V0_1_OUT')
_TIMEOUT = int(os.environ.get('MUTATION_V0_1_TIMEOUT', '300'))
_WATCH = ('t1', 't2', 't3', 't4', 't5', 't6', 'audit_t4', 'audit_t1_mid')
_STATE = dict(outcomes={}, fixtures={}, t0=time.time())


class MutationTimeout(Exception):
    pass


def _on_alarm(signum, frame):
    raise MutationTimeout(f'stopped after {_TIMEOUT} s')


def _alarm(on):
    signal.signal(signal.SIGALRM, _on_alarm)
    signal.alarm(_TIMEOUT if on else 0)


def _mx(x):
    return max(x) if len(x) else None


def _sum_t1(r):
    return dict(eps=r['eps'], err_M=r['err_M'], err_Jr=r['err_Jr'],
                fitted_rate_err_M=r['fitted_rate_err_M'], fitted_rate_err_Jr=r['fitted_rate_err_Jr'],
                min_rate_err_M=min(r['rates_err_M']), min_rate_err_Jr=min(r['rates_err_Jr']),
                kappa=r['kappa'], mass_constrained_max=_mx(r['mass_constrained']),
                err_M_by_window=r['err_M_by_window'])


def _sum_t2(r):
    return dict(eps=r['eps'], err_MH_vs_I=r['err_MH_vs_I'], err_MS_vs_MRA=r['err_MS_vs_MRA'],
                err_Magg_vs_TANK=r['err_Magg_vs_TANK'],
                err_Magg_chi_vs_TANK=r['err_Magg_chi_vs_TANK'], err_IRF_vs_TANK=r['err_IRF_vs_TANK'],
                fitted_rate_err_MS_vs_MRA=r['fitted_rate_err_MS_vs_MRA'],
                fitted_rate_err_IRF={c: r[f'fitted_rate_err_IRF_chi{c}'] for c in r['err_IRF_vs_TANK']},
                amplification_HA_last={c: v[-1] for c, v in r['amplification_HA'].items()},
                amplification_factor=r['amplification_factor'],
                ge_target_resid=r['ge_target_resid'])


def _sum_t3(r):
    return dict(b=r['b'], err_acyc_vs_RANK=r['err_acyc_vs_RANK'], err_cyc_vs_ZL=r['err_cyc_vs_ZL'],
                dev_cyc_vs_RANK=r['dev_cyc_vs_RANK'], zl_cyc_vs_RANK=r['zl_cyc_vs_RANK'],
                fitted_rate_acyc=r['fitted_rate_err_acyc_vs_RANK_asymptotic'],
                fitted_rate_cyc=r['fitted_rate_err_cyc_vs_ZL_asymptotic'],
                r_gap_over_b_acyc=r['r_gap_over_b_acyc'], delta=r['zero_liquidity_cyc']['delta'])


def _sum_t4(r):
    return dict(eps=r['eps'], r=r['r'], gap=r['gap'], fitted_rate_gap=r['fitted_rate_gap'],
                capital_market_resid_rel_max=_mx(r['capital_market_resid_rel']),
                euler_err_mean_weighted_max=_mx(r['euler_err_mean_weighted']),
                euler_err_max_with_mass_max=_mx(r['euler_err_max_with_mass']),
                kkt_violation_max_with_mass_max=_mx(r['kkt_violation_max_with_mass']),
                mass_top_5pts_max=_mx(r['mass_top_5pts']),
                mass_policy_off_grid_max=_mx(r['mass_policy_off_grid']),
                r_by_sd=r['r_by_sd']['r'], r_by_rho=r['r_by_rho']['r'],
                monotone_in_sd=r['monotone_in_sd'], monotone_in_rho=r['monotone_in_rho'],
                all_below_rho=r['all_below_rho'])


def _sum_t5(r):
    return {k: dict(max_err=v['max_err'], err=v['err'],
                    walras_max={lab: _mx(x) for lab, x in v['walras_asset_mkt'].items()})
            for k, v in r.items()}


def _sum_t6(r):
    out = {}
    for T, v in r['by_T'].items():
        out[T] = dict(identity_err=v['G']['identity_err'],
                      decomposition_err_G=v['G']['decomposition_err'],
                      decomposition_err_r=v['r']['decomposition_err'],
                      cum_multiplier_RA=v['G']['cum_multiplier_RA'],
                      cum_multiplier_HA=v['G']['cum_multiplier_HA'],
                      group_sum_vs_M_x=v['group_sum_vs_M_x'],
                      group_first_round_sum_vs_total=v['group_first_round_sum_vs_total'],
                      first_round_0=[g['first_round'][0] for g in v['groups']])
    return out


def _sum_audit(a):
    ref = a['reference']
    return dict(sweeps=ref['sweeps'], sweep_tol_met=ref['sweep_tol_met'],
                neg_mass_rel_max=ref['neg_mass_rel_max'], resid_L1=ref['resid_L1'],
                shift_L1_diff_max=_mx([v['L1_diff'] for v in a['shift'].values()]),
                shift_dA_rel_max=_mx([v['dA_rel'] for v in a['shift'].values()]),
                ssj_forward={k: {f: v.get(f) for f in ('L1_diff', 'dA_rel', 'mirror_identical',
                                                       'iterations', 'error')}
                             for k, v in a['ssj_forward'].items()})


_SUMMARISE = dict(t1=_sum_t1, t2=_sum_t2, t3=_sum_t3, t4=_sum_t4, t5=_sum_t5, t6=_sum_t6,
                  audit_t4=_sum_audit, audit_t1_mid=_sum_audit)


@pytest.hookimpl(wrapper=True)
def pytest_fixture_setup(fixturedef, request):
    name = fixturedef.argname
    if not _OUT or name not in _WATCH:
        return (yield)
    t0 = time.time()
    _alarm(True)
    try:
        result = yield
    except BaseException as exc:
        _STATE['fixtures'][name] = dict(error=f'{type(exc).__name__}: {exc}'[:500],
                                        runtime_s=time.time() - t0)
        raise
    finally:
        _alarm(False)
    rec = dict(runtime_s=time.time() - t0)
    try:
        rec['summary'] = _clean(_SUMMARISE[name](result))
    except Exception as exc:              # a malformed result must not stop the run
        rec['summary_error'] = f'{type(exc).__name__}: {exc}'[:500]
    _STATE['fixtures'][name] = rec
    return result


@pytest.hookimpl(wrapper=True)
def pytest_runtest_call(item):
    if not _OUT:
        return (yield)
    _alarm(True)
    try:
        return (yield)
    finally:
        _alarm(False)


def _message(report):
    text = getattr(report, 'longreprtext', '') or ''
    lines = [ln[1:].strip() for ln in text.splitlines() if ln.startswith('E ')]
    msg = ' | '.join(lines[:4]) if lines else (text.strip().splitlines() or [''])[-1]
    return msg[:400]


def pytest_runtest_logreport(report):
    if not _OUT:
        return
    rec = _STATE['outcomes'].setdefault(report.nodeid, dict(outcome=None, message=''))
    if report.when == 'setup':
        if report.failed:
            rec.update(outcome='error', message=_message(report))
        elif report.skipped:
            rec.update(outcome='xfailed' if hasattr(report, 'wasxfail') else 'skipped')
    elif report.when == 'call':
        if hasattr(report, 'wasxfail'):
            rec['outcome'] = 'xfailed' if report.skipped else 'xpassed'
        else:
            rec['outcome'] = report.outcome
        if report.failed:
            rec['message'] = _message(report)
    elif report.when == 'teardown' and report.failed and rec['outcome'] == 'passed':
        rec.update(outcome='error', message=_message(report))


def pytest_sessionfinish(session, exitstatus):
    if not _OUT:
        return
    with open(_OUT, 'w') as f:
        json.dump(dict(outcomes=_STATE['outcomes'], fixtures=_STATE['fixtures'],
                       exitstatus=int(exitstatus), runtime_s=time.time() - _STATE['t0']), f)


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def _clean(x):
    """JSON-safe copy: numpy scalars -> Python, non-finite floats -> strings."""
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if hasattr(x, 'item') and not isinstance(x, (str, bytes)):
        try:
            x = x.item()
        except (ValueError, AttributeError):
            return str(x)
    if isinstance(x, float) and not math.isfinite(x):
        return str(x)
    return x


def _copy_tree(dst):
    shutil.copytree(os.path.join(ROOT, 'hank_limits'), os.path.join(dst, 'hank_limits'),
                    ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(os.path.join(ROOT, 'tests'), os.path.join(dst, 'tests'),
                    ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy(os.path.join(ROOT, 'pytest.ini'), dst)


def _apply(dst, patches):
    for fname, old, new in patches:
        path = os.path.join(dst, 'hank_limits', fname)
        with open(path) as f:
            src = f.read()
        n = src.count(old)
        if n != 1:
            raise ValueError(f'{fname}: patch target occurs {n} times: {old[:60]!r}')
        with open(path, 'w') as f:
            f.write(src.replace(old, new))
    for fname in sorted({p[0] for p in patches}):
        py_compile.compile(os.path.join(dst, 'hank_limits', fname), doraise=True)


def _thread_env():
    env = dict(os.environ)
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMBA_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[v] = '1'
    return env


def run_one(run_id, patches, workdir, fixture_timeout, run_timeout):
    tmp = tempfile.mkdtemp(prefix=f'mut_{run_id}_', dir=workdir)
    try:
        _copy_tree(tmp)
        _apply(tmp, patches)
        out = os.path.join(tmp, 'result.json')
        env = _thread_env()
        env.update(MUTATION_V0_1_OUT=out, MUTATION_V0_1_TIMEOUT=str(fixture_timeout),
                   PYTHONPATH=os.path.join(ROOT, 'scripts'), PYTHONDONTWRITEBYTECODE='1')
        t0 = time.time()
        try:
            p = subprocess.run([sys.executable, '-m', 'pytest', '-q', '-p', 'mutation_v0_1',
                                '-p', 'no:cacheprovider'], cwd=tmp, env=env,
                               capture_output=True, text=True, timeout=run_timeout)
            tail = p.stdout.strip().splitlines()[-1:] if p.stdout else []
            status = dict(returncode=p.returncode, pytest_summary=tail[0] if tail else '')
        except subprocess.TimeoutExpired:
            status = dict(returncode=None, pytest_summary=f'run stopped after {run_timeout} s')
        wall = time.time() - t0
        res = json.load(open(out)) if os.path.exists(out) else dict(outcomes={}, fixtures={})
        res.update(status, wall_s=wall)
        return res
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _flatten(x, prefix=''):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from _flatten(v, f'{prefix}{k}.')
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from _flatten(v, f'{prefix}{i}.')
    else:
        yield prefix.rstrip('.'), x


def _group_fixtures(group):
    return {'T1': ['t1'], 'T2': ['t2'], 'T3': ['t3'], 'T4': ['t4'], 'T5': ['t5'], 'T6': ['t6'],
            'D': ['audit_t4', 'audit_t1_mid']}[group]


def _max_rel_change(run, base, group):
    """Largest relative change of any summarised number of the group's checks vs baseline
    (0.0 = bit-identical; inf = the check failed to produce a result)."""
    worst = 0.0
    for fx in _group_fixtures(group):
        a, b = run['fixtures'].get(fx, {}), base['fixtures'].get(fx, {})
        if 'summary' not in a or 'summary' not in b:
            return math.inf if 'summary' in b else None
        fa, fb = dict(_flatten(a['summary'])), dict(_flatten(b['summary']))
        if fa.keys() != fb.keys():
            return math.inf
        for k, vb in fb.items():
            va = fa[k]
            if isinstance(vb, bool) or not isinstance(vb, (int, float)) or \
                    not isinstance(va, (int, float)):
                if va != vb:
                    worst = math.inf
                continue
            if va != vb:
                worst = max(worst, abs(va - vb) / max(abs(vb), 1e-300))
    return worst


def _headline(group, summ):
    """Two indicative error measures per group (the tightest setting of each check)."""
    try:
        if group == 'T1':
            s = summ['t1']['summary']
            return [('err_Jr_at_eps0', s['err_Jr'][0]),
                    ('err_max_at_eps_min', max(s['err_M'][-1], s['err_Jr'][-1]))]
        if group == 'T2':
            s = summ['t2']['summary']
            return [('err_IRF_chi1.5_at_eps_min', s['err_IRF_vs_TANK']['1.5'][-1]),
                    ('ge_target_resid_max', max(s['ge_target_resid']))]
        if group == 'T3':
            s = summ['t3']['summary']
            return [('err_vs_ref_at_b_min', max(s['err_acyc_vs_RANK'][-1], s['err_cyc_vs_ZL'][-1])),
                    ('fitted_rate_acyc', s['fitted_rate_acyc'])]
        if group == 'T4':
            s = summ['t4']['summary']
            return [('gap_at_eps_min', s['gap'][-1]), ('fitted_rate_gap', s['fitted_rate_gap'])]
        if group == 'T5':
            s = summ['t5']['summary']['ssj_consistent_refined']
            return [('ikc_max_err_refined', s['max_err']),
                    ('walras_max', max(v for v in s['walras_max'].values()))]
        if group == 'T6':
            s = summ['t6']['summary']
            return [('decomposition_err_max', max(max(v['decomposition_err_G'], v['decomposition_err_r'])
                                                  for v in s.values())),
                    ('abs_cum_mult_RA_minus_1_max', max(abs(v['cum_multiplier_RA'] - 1)
                                                        for v in s.values()))]
        if group == 'D':
            v = [summ[k]['summary']['ssj_forward']['1e-12'].get('dA_rel') for k in ('audit_t4', 'audit_t1_mid')]
            return [('ssj_fwd_dA_rel_1e-12_max', max(x for x in v if x is not None)),
                    ('neg_mass_rel_max', max(summ[k]['summary']['neg_mass_rel_max']
                                             for k in ('audit_t4', 'audit_t1_mid')))]
    except (KeyError, TypeError, ValueError):
        pass
    return [('n/a', None), ('n/a', None)]


def _mode(rec):
    if rec['outcome'] == 'error':
        return 'timeout' if 'MutationTimeout' in rec['message'] else 'error'
    if rec['outcome'] == 'failed':
        m = rec['message']
        if 'XPASS' in m:
            return 'xpass'
        return 'assert' if ('assert' in m or 'AssertionError' in m) else 'exception'
    return rec['outcome']


def analyse(runs, base_id='B0'):
    base = runs[base_id]
    base_out = {k: v['outcome'] for k, v in base['outcomes'].items()}
    table = {}
    for m in MUTANTS:
        r = runs.get(m['id'])
        if r is None:
            continue
        reach, expect, unsure = set(_groups(m['reach'])), set(_groups(m['expect'])), set(_groups(m.get('unsure', '')))
        per = {}
        for g in GROUPS:
            ids = [k for k in base_out if TEST_FILE_GROUP.get(k.split('::')[0].split('/')[-1]) == g]
            changed = []
            for k in ids:
                rec = r['outcomes'].get(k, dict(outcome='missing', message=''))
                if rec['outcome'] != base_out[k]:
                    changed.append(dict(test=k.split('::')[-1], baseline=base_out[k],
                                        outcome=rec['outcome'], mode=_mode(rec),
                                        message=rec['message']))
            detected = bool(changed)
            pred = 'catch' if g in expect else ('?' if g in unsure else 'miss')
            modes = sorted({c['mode'] for c in changed})
            per[g] = dict(reached=g in reach, affected=_max_rel_change(r, base, g),
                          predicted=pred, detected=detected,
                          prediction_ok=None if pred == '?' else (detected == (pred == 'catch')),
                          modes=modes, n_tests=len(ids), changed=changed,
                          headline=_headline(g, r['fixtures']),
                          headline_baseline=_headline(g, base['fixtures']))
        table[m['id']] = dict(category=m['category'], target=m['target'], change=m['change'],
                              conf=m['conf'], detected_any=any(v['detected'] for v in per.values()),
                              groups=per, pytest_summary=r.get('pytest_summary'),
                              wall_s=r.get('wall_s'))
    # rates (mutants only; controls reported separately)
    mut = [m for m in MUTANTS if m['category'] != 'control' and m['id'] in table]
    rates = {}
    for g in GROUPS:
        reached = [m['id'] for m in mut if table[m['id']]['groups'][g]['reached']]
        caught = [i for i in reached if table[i]['groups'][g]['detected']]
        caught_all = [m['id'] for m in mut if table[m['id']]['groups'][g]['detected']]
        rates[g] = dict(caught=len(caught_all), of_all=len(mut), caught_reached=len(caught),
                        reached=len(reached), caught_ids=caught_all)
    cats = {}
    for m in mut:
        c = cats.setdefault(m['category'], dict(n=0, caught=0, ids_missed=[], by_group={g: 0 for g in GROUPS}))
        c['n'] += 1
        c['caught'] += table[m['id']]['detected_any']
        if not table[m['id']]['detected_any']:
            c['ids_missed'].append(m['id'])
        for g in GROUPS:
            c['by_group'][g] += table[m['id']]['groups'][g]['detected']
    suite = dict(caught=sum(table[m['id']]['detected_any'] for m in mut), of=len(mut),
                 caught_by_T1_T4=sum(any(table[m['id']]['groups'][g]['detected'] for g in ('T1', 'T2', 'T3', 'T4'))
                                     for m in mut),
                 caught_by_T5_T6=sum(any(table[m['id']]['groups'][g]['detected'] for g in ('T5', 'T6'))
                                     for m in mut))
    controls = {m['id']: dict(false_alarm=table[m['id']]['detected_any'],
                              groups=[g for g in GROUPS if table[m['id']]['groups'][g]['detected']])
                for m in MUTANTS if m['category'] == 'control' and m['id'] in table}
    mismatches = [dict(id=i, group=g, predicted=v['predicted'], detected=v['detected'])
                  for i, t in table.items() for g, v in t['groups'].items()
                  if v['prediction_ok'] is False]
    return dict(table=table, rates_by_group=rates, rates_by_category=cats, suite=suite,
                controls=controls, prediction_mismatches=mismatches)


def environment():
    from importlib import metadata
    pk = {}
    for name in ('sequence-jacobian', 'numpy', 'scipy', 'numba', 'llvmlite', 'pytest'):
        try:
            pk[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            pk[name] = None
    blas = None
    try:
        import numpy as np
        cfg = np.show_config(mode='dicts')
        blas = cfg.get('Build Dependencies', {}).get('blas')
    except Exception:
        pass
    try:
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True,
                                text=True).stdout.strip()
        dirty = bool(subprocess.run(['git', 'status', '--porcelain', 'hank_limits', 'tests'],
                                    cwd=ROOT, capture_output=True, text=True).stdout.strip())
    except OSError:
        commit, dirty = None, None
    return dict(platform=platform.platform(), machine=platform.machine(),
                python=sys.version.split()[0], packages=pk, numpy_blas=blas,
                cpu_count=os.cpu_count(), git_commit=commit, model_or_tests_dirty=dirty,
                threads_per_run=1, utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec='seconds'))


def write_outputs(prefix, env, runs, res, args):
    rows = []
    for i, t in res['table'].items():
        for g, v in t['groups'].items():
            (h1, x1), (h2, x2) = v['headline']
            (_, b1), (_, b2) = v['headline_baseline']
            rows.append(dict(id=i, category=t['category'], target=t['target'], group=g,
                             reached=v['reached'], affected_max_rel_change=v['affected'],
                             predicted=v['predicted'], detected=v['detected'],
                             prediction_ok=v['prediction_ok'], modes=';'.join(v['modes']),
                             n_tests=v['n_tests'], n_changed=len(v['changed']),
                             changed_tests=';'.join(c['test'] for c in v['changed']),
                             metric1=h1, value1=x1, baseline1=b1, metric2=h2, value2=x2, baseline2=b2,
                             first_message=v['changed'][0]['message'] if v['changed'] else ''))
    if rows:
        with open(prefix + '.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    # every run's outcomes and check summaries (v0.1 amendment A1: mutant runs were not stored)
    base = {k: dict(pytest_summary=runs[k].get('pytest_summary'), wall_s=runs[k].get('wall_s'),
                    fixtures=runs[k]['fixtures'], outcomes=runs[k]['outcomes'])
            for k in runs}
    with open(prefix + '.json', 'w') as f:
        json.dump(_clean(dict(environment=env, settings=vars(args),
                              mutants=[{k: v for k, v in m.items() if k != 'patches'} | dict(
                                  patches=[dict(file=p[0], old=p[1], new=p[2]) for p in m['patches']])
                                       for m in MUTANTS],
                              runs=base, baseline_deterministic=res.get('baseline_deterministic'),
                              results=res)), f, indent=1)


def print_table():
    print('| ID | Category | Target | Change | Reach | Predicted catch | Unsure | Conf. |')
    print('|---|---|---|---|---|---|---|---|')
    for m in MUTANTS:
        print(f"| {m['id']} | {m['category']} | `{m['target']}` | {m['change']} | {m['reach']} | "
              f"{m['expect'] or 'none'} | {m.get('unsure', '') or '-'} | {m['conf']} |")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--check', action='store_true', help='apply and compile every patch; run no tests')
    ap.add_argument('--table', action='store_true', help='print the mutant table and exit')
    ap.add_argument('--only', nargs='*', default=None, help='mutant ids to run')
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--fixture-timeout', type=int, default=300)
    ap.add_argument('--run-timeout', type=int, default=3600)
    ap.add_argument('--workdir', default=None, help='parent directory for the temporary copies')
    ap.add_argument('--out', default=os.path.join(ROOT, 'out', 'mutation_v0_1'))
    args = ap.parse_args()
    if args.table:
        print_table()
        return
    sel = [m for m in MUTANTS if args.only is None or m['id'] in args.only]
    if args.check:
        for m in sel:
            tmp = tempfile.mkdtemp(prefix=f'chk_{m["id"]}_', dir=args.workdir)
            try:
                _copy_tree(tmp)
                _apply(tmp, m['patches'])
                print(f"{m['id']}: patch applies and compiles")
            finally:
                shutil.rmtree(tmp, ignore_errors=True)
        return
    env = environment()
    jobs = [('B0', []), ('B1', [])] + [(m['id'], [(p[0], p[1], p[2]) for p in m['patches']]) for m in sel]
    runs = {}
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_one, i, p, args.workdir, args.fixture_timeout, args.run_timeout): i
                for i, p in jobs}
        for fut in futs:
            i = futs[fut]
            runs[i] = fut.result()
            print(f"[{time.time() - t0:7.1f}s] {i}: {runs[i].get('pytest_summary')}", flush=True)
    res = analyse(runs)
    res['baseline_deterministic'] = all(
        _max_rel_change(runs['B1'], runs['B0'], g) == 0.0 for g in GROUPS) and \
        {k: v['outcome'] for k, v in runs['B0']['outcomes'].items()} == \
        {k: v['outcome'] for k, v in runs['B1']['outcomes'].items()}
    write_outputs(args.out, env, runs, res, args)
    print(f"baselines: {runs['B0'].get('pytest_summary')} / {runs['B1'].get('pytest_summary')}; "
          f"deterministic: {res['baseline_deterministic']}")
    print(f"suite caught {res['suite']['caught']} of {res['suite']['of']} mutants")
    for g, v in res['rates_by_group'].items():
        print(f"  {g}: caught {v['caught']}/{v['of_all']} (of reached: {v['caught_reached']}/{v['reached']})")
    print('controls:', res['controls'])
    print('prediction mismatches:', [(x['id'], x['group'], x['predicted']) for x in res['prediction_mismatches']])


if __name__ == '__main__':
    main()
