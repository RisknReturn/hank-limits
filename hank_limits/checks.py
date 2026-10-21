"""Model-agnostic analytic-limit checks.

Every check takes a *provider* (a function implemented by the HANK code under test; see
`hank_limits.interface`) and returns a JSON-serialisable dict with the errors along the
limit sequence and the observed convergence rates. The tests in `tests/` and the report
script call these with the SSJ providers from `hank_limits.models`.
"""
from __future__ import annotations

import time
import numpy as np

from .metrics import rel_max_err, rel_fro_err, pairwise_rates, fitted_rate
from .ra import (ra_jacobians, ra_output_fixed_real_rate, tank_output_fixed_real_rate,
                 tank_impc, zero_liquidity_output, zero_liquidity_rstar)


def _tolist(x):
    return [float(v) for v in np.asarray(x).ravel()]


# ---------------------------------------------------------------------------
# v0.3 checks (docs/CHECKS_PLAN_v0_3.md): steady-state residuals (V2), income normalisation
# (V3) and specification pins (V4), computed here from primitive quantities the provider
# reports in the optional fields `ss_detail`, `J_A`, `policy` and the paths. They are stored
# under res['v0_3'] (one list entry per setting of the check), so the v0.1 part of every
# result is unchanged. A missing field skips the quantities that need it; the reason is
# recorded in res['v0_3']['skipped'].
# ---------------------------------------------------------------------------

def stationary_marginal(Pi):
    """Stationary distribution pi of a reported transition matrix (pi Pi = pi, sum pi = 1)."""
    Pi = np.asarray(Pi, float)
    n = Pi.shape[0]
    A = Pi.T - np.eye(n)
    A[-1] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    return np.linalg.solve(A, b)


def income_report(e, Pi, y=None, Y_ref=None):
    """V3: Pi row-stochastic, E_pi[e] = 1 and, if y is given, sum_e pi_e y_e = Y_ref."""
    Pi = np.asarray(Pi, float)
    pi = stationary_marginal(Pi)
    out = dict(pi_row_sum_err=float(np.max(np.abs(Pi.sum(1) - 1))),
               pi_negative=float(max(0.0, -Pi.min())),
               mean_e_err=float(abs(pi @ np.asarray(e, float) - 1)))
    if y is not None and Y_ref is not None:
        out['agg_income_err'] = float(abs(pi @ np.asarray(y, float) - Y_ref) / abs(Y_ref))
    return out, pi


def household_budget_resid(C, A, r, y, pi):
    """V2: aggregate household budget at the steady state, |C - r A - sum_e pi_e y_e| / C."""
    return float(abs(C - r * A - pi @ np.asarray(y, float)) / abs(C))


def euler_errors(Pi, a_grid, a_pol, c_pol, D, r, beta, eis, bind_tol=1e-10, mass_tol=1e-12):
    """V4(b): Euler errors of a reported steady-state policy at a given (beta, eis). Same
    definitions as `hank_limits.diagnostics.euler_kkt_residuals` (numpy only): c_impl =
    (beta (1+r) E[c(e', a')^(-1/eis)])^(-eis), err = c_impl / c - 1; unit-free |err| on
    unconstrained states, max(0, -err) (KKT violation) on binding states."""
    Pi, a_grid = np.asarray(Pi, float), np.asarray(a_grid, float)
    a_pol, c_pol, D = (np.asarray(x, float) for x in (a_pol, c_pol, D))
    n_e, n_a = a_pol.shape
    up = np.empty_like(c_pol)
    for e in range(n_e):
        x = a_pol[e]
        j = np.clip(np.searchsorted(a_grid, x) - 1, 0, n_a - 2)
        w = (x - a_grid[j]) / (a_grid[j + 1] - a_grid[j])
        cn = (1 - w)[None, :] * c_pol[:, j] + w[None, :] * c_pol[:, j + 1]
        up[e] = Pi[e] @ np.maximum(cn, 1e-300) ** (-1 / eis)
    err = (beta * (1 + r) * up) ** (-eis) / c_pol - 1
    binding = a_pol <= a_grid[0] + bind_tol
    unc = ~binding & ~(a_pol > a_grid[-1])
    has_mass = D > mass_tol
    w = D[unc]
    return dict(
        euler_err_mean_weighted=float((np.abs(err[unc]) * w).sum() / max(w.sum(), 1e-300)),
        euler_err_max_with_mass=float(np.abs(err[unc & has_mass]).max()) if np.any(unc & has_mass) else 0.0,
        kkt_violation_max_with_mass=float(np.maximum(0.0, -err[binding & has_mass]).max())
        if np.any(binding & has_mass) else 0.0)


def _new_v03():
    return dict(skipped=[])


def _guard(v, fn, *args):
    """Run one v0.3 computation; an exception is recorded (the V tests then fail) and never
    stops the v0.1 part of the check, so groups T1-T6 behave as in v0.2.1."""
    try:
        fn(v, *args)
    except Exception as exc:          # noqa: BLE001
        v.setdefault('errors', []).append(f'{fn.__name__}: {type(exc).__name__}: {exc}'[:300])


def _put(v, key, value):
    v.setdefault(key, []).append(None if value is None else float(value))


def _skip(v, what, missing):
    msg = f"{what}: provider does not report {', '.join(missing)}"
    if msg not in v['skipped']:
        v['skipped'].append(msg)


def _need(d, keys):
    return [k for k in keys if d is None or k not in d]


def _pin(v, key, value, ref, rel=False):
    _put(v, key, abs(value - ref) / (abs(ref) if rel else 1.0))


def _income(v, suffix, e, Pi, y, Y_ref):
    inc, pi = income_report(e, Pi, y, Y_ref)
    for k, x in inc.items():
        _put(v, k + suffix, x)
    return pi


def _v03_rank(v, out, spec):
    d = out.get('ss_detail')
    miss = _need(d, ('Y', 'C', 'A', 'r', 'eis', 'a_min', 'e_grid', 'Pi', 'y'))
    if miss:
        return _skip(v, 'V2/V3/V4 (T1)', miss)
    pi = _income(v, '', d['e_grid'], d['Pi'], d['y'], d['Y'])
    _put(v, 'hh_budget', household_budget_resid(d['C'], d['A'], d['r'], d['y'], pi))
    if spec:
        _pin(v, 'spec_r', d['r'], spec['r'])
        _pin(v, 'spec_eis', d['eis'], spec['eis'])
        _pin(v, 'spec_A', d['A'], spec['A'], rel=True)
        nbl = -spec['nbl_frac'] * spec['Y'] * np.min(d['e_grid']) / spec['r']
        _pin(v, 'spec_a_min', d['a_min'], nbl, rel=True)


def _v03_tank(v, out, spec):
    d = out.get('ss_detail')
    keys = ['lam', 'r', 'eis', 'Y', 'Y_H', 'Y_S'] + [f'{k}_{t}' for t in 'HS' for k in
                                                    ('C', 'A', 'beta', 'a_min', 'e_grid', 'Pi', 'y')]
    miss = _need(d, keys)
    if miss:
        return _skip(v, 'V2/V3/V4 (T2)', miss)
    lam = d['lam']
    _put(v, 'ss_goods', abs(d['Y'] - lam * d['C_H'] - (1 - lam) * d['C_S']) / d['Y'])
    _put(v, 'ss_assets', abs(lam * d['A_H'] + (1 - lam) * d['A_S']) / d['Y'])
    for t in 'HS':
        pi = _income(v, f'_{t}', d[f'e_grid_{t}'], d[f'Pi_{t}'], d[f'y_{t}'], d[f'Y_{t}'])
        _put(v, f'hh_budget_{t}', household_budget_resid(d[f'C_{t}'], d[f'A_{t}'], d['r'],
                                                         d[f'y_{t}'], pi))
    _put(v, 'incidence_err', abs(lam * d['Y_H'] + (1 - lam) * d['Y_S'] - d['Y']) / d['Y'])
    if spec:
        _pin(v, 'spec_lam', lam, spec['lam'])
        _pin(v, 'spec_beta_H', d['beta_H'], spec['beta_H'])
        _pin(v, 'spec_r', d['r'], spec['r'])
        _pin(v, 'spec_eis', d['eis'], spec['eis'])
        _pin(v, 'spec_a_min_H', d['a_min_H'], spec['a_min_H'])
        nbl = -spec['nbl_frac'] * spec['Y'] * np.min(d['e_grid_S']) / spec['r']
        _pin(v, 'spec_a_min_S', d['a_min_S'], nbl, rel=True)


def _v03_werning(v, out, b, dr, cyc, spec):
    tag = 'cyc' if cyc else 'acyc'
    gt = out.get('diag', {}).get('ge_target_resid')
    _put(v, f'ge_target_resid_reported_{tag}', None if gt is None else gt / b)
    JA = out.get('J_A')
    if JA is None:
        _skip(v, 'V1(b)', ['J_A'])
    else:
        _put(v, f'ge_target_resid_direct_{tag}',
             np.max(np.abs(np.asarray(JA['Y']) @ out['dY'] + np.asarray(JA['r']) @ dr)) / b)
    d = out.get('ss_detail')
    miss = _need(d, ('Y', 'C', 'A', 'r', 'a_min', 'y'))
    if miss:
        return _skip(v, 'V2/V3/V4 (T3)', miss)
    e, Pi, zeta = (np.asarray(out[k], float) for k in ('e_grid', 'Pi', 'zeta'))
    _put(v, f'ss_assets_{tag}', abs(d['A']) / b)
    _put(v, f'ss_goods_{tag}', abs(d['C'] - d['Y']) / d['Y'])
    pi = _income(v, f'_{tag}', e, Pi, d['y'], d['Y'])
    _put(v, f'hh_budget_{tag}', household_budget_resid(d['C'], d['A'], d['r'], d['y'], pi))
    if cyc:
        _put(v, 'incidence_err_cyc', abs(pi @ (e * zeta)))
    if spec:
        _pin(v, f'spec_eis_{tag}', out['eis'], spec['eis'])
        _pin(v, f'spec_a_min_{tag}', d['a_min'], -b, rel=True)
        if cyc:
            _pin(v, 'spec_zeta_low_cyc', zeta[0], spec['zeta_low_cyc'])
        else:
            _put(v, 'spec_zeta_acyc', np.max(np.abs(zeta)))
        _pin(v, f'spec_rstar_{tag}', zero_liquidity_rstar(out['beta'], out['eis'], e, Pi),
             spec['r_star'])


def _v03_aiyagari(v, out, spec):
    d = out.get('ss_detail')
    miss = _need(d, ('r', 'w', 'K', 'A', 'C', 'beta', 'eis', 'a_min', 'e_grid', 'Pi', 'y'))
    if miss:
        return _skip(v, 'V2/V3/V4 (T4)', miss)
    L_doc = spec['L'] if spec else 1.0
    pi = _income(v, '', d['e_grid'], d['Pi'], d['y'], d['w'] * L_doc)
    _put(v, 'ss_capital', abs(d['A'] - d['K']) / d['K'])
    _put(v, 'hh_budget', household_budget_resid(d['C'], d['A'], d['r'], d['y'], pi))
    if not spec:
        return
    al, de = spec['alpha'], spec['delta']
    L = float(pi @ np.asarray(d['e_grid'], float))
    KL = d['K'] / L
    _put(v, 'foc_r', abs(d['r'] + de - al * KL ** (al - 1)))
    _put(v, 'foc_w', abs(d['w'] - (1 - al) * KL ** al) / d['w'])
    Y = d['K'] ** al * L ** (1 - al)
    _put(v, 'ss_goods', abs(Y - d['C'] - de * d['K']) / Y)
    _pin(v, 'spec_beta', d['beta'], spec['beta'])
    _pin(v, 'spec_eis', d['eis'], spec['eis'])
    _pin(v, 'spec_a_min', d['a_min'], spec['a_min'])
    pol = out.get('policy')
    miss = _need(pol, ('a_grid', 'a', 'c', 'D'))
    if miss:
        return _skip(v, 'V4(b) Euler (T4)', miss)
    for k, x in euler_errors(d['Pi'], pol['a_grid'], pol['a'], pol['c'], pol['D'], d['r'],
                             spec['beta'], spec['eis']).items():
        _put(v, 'spec_' + k, x)


def _lag(x):
    return np.concatenate([[0.0], np.asarray(x, float)[:-1]])


def _v03_fiscal(v, out, shocks, spec):
    d = out.get('ss_detail')
    miss = _need(d, ('Y', 'C', 'A', 'B', 'G', 'Tax', 'r', 'eis', 'phi_B', 'a_min', 'e_grid',
                     'Pi', 'y'))
    if miss:
        _skip(v, 'V2/V3/V4(a) (T5)', miss)
    else:
        _put(v, 'ss_government_budget', abs(d['Tax'] - d['G'] - d['r'] * d['B']) / d['Y'])
        _put(v, 'ss_goods', abs(d['Y'] - d['C'] - d['G']) / d['Y'])
        _put(v, 'ss_assets', abs(d['A'] - d['B']) / d['B'])
        pi = _income(v, '', d['e_grid'], d['Pi'], d['y'], d['Y'] - d['Tax'])
        _put(v, 'hh_budget', household_budget_resid(d['C'], d['A'], d['r'], d['y'], pi))
        if spec:
            for k in ('r', 'eis', 'B', 'G', 'phi_B', 'a_min'):
                _pin(v, f'spec_{k}', d[k], spec[k])
    if not spec:
        return
    paths = v.setdefault('paths', {})
    for label, irf in out['irfs'].items():
        kind, inputs = shocks[label]
        s = max(float(np.max(np.abs(p))) for p in inputs.values())
        p = paths.setdefault(label, {})
        dB, dTax = np.asarray(irf['dB'], float), np.asarray(irf['dTax'], float)
        _put(p, 'fiscal_rule', np.max(np.abs(dTax - spec['phi_B'] * _lag(dB))) / s)
        _put(p, 'government_budget',
             np.max(np.abs(dB - (1 + spec['r']) * _lag(dB) - spec['B'] * irf['dr'] - irf['dG']
                           + dTax)) / s)
        _put(p, 'goods_market', np.max(np.abs(irf['dY'] - irf['dC'] - irf['dG'])) / s)
        if kind == 'r_rule':
            if 'drshock' not in irf:
                _skip(v, 'V4(b) rate rule (T5)', ['irfs[label][drshock]'])
            else:
                _put(p, 'rate_rule',
                     np.max(np.abs(irf['dr'] - spec['phi_y'] * irf['dY'] - irf['drshock'])) / s)


# ---------------------------------------------------------------------------
# T1: RANK limit
# ---------------------------------------------------------------------------

def check_rank_limit(provider, eps_list, T, T_sub=(50, 100), spec=None, **kw):
    """PIH household-block limit (partial equilibrium; not a test of the full NK GE RANK
    limit): ||J_HA(eps) - J_PIH|| for J in {dC/dY, dC/dr}, J_PIH the permanent-income closed
    form at the provider's steady state (r, eis, C, A), along the path eps -> 0 on which beta is
    recalibrated to keep mean assets at the target.

    Also reports the error on the leading T_sub x T_sub block of the SAME T x T Jacobian
    (`err_*_by_window`; `err_*_by_T` is kept as an alias). This is an observation-window
    diagnostic: it shows that the gap sits at long lags. It is not a terminal-horizon audit;
    SSJ's household Jacobian entries do not depend on the computing horizon (checked in
    scripts/diagnostics_v0_2.py, section C: T=50 and T=100 Jacobians computed separately equal
    the leading blocks of the T=200 one exactly)."""
    t0 = time.time()
    res = dict(eps=[], err_M=[], err_Jr=[], err_M_fro=[], err_Jr_fro=[], kappa=[],
               mass_constrained=[], argmax_err_M=[], err_M_by_T={str(t): [] for t in T_sub},
               err_Jr_by_T={str(t): [] for t in T_sub}, ssj_forward_resid=[])
    v03 = _new_v03()
    for eps in eps_list:
        out = provider(eps, T, **kw)
        _guard(v03, _v03_rank, out, spec)
        s = out['ss']
        M_ra, Jr_ra = ra_jacobians(T, s['r'], s['eis'], s['C'], s['A'])
        M, Jr = out['J']['C']['Y'], out['J']['C']['r']
        res['eps'].append(float(eps))
        res['err_M'].append(rel_max_err(M, M_ra))
        res['err_Jr'].append(rel_max_err(Jr, Jr_ra))
        res['err_M_fro'].append(rel_fro_err(M, M_ra))
        res['err_Jr_fro'].append(rel_fro_err(Jr, Jr_ra))
        i, j = np.unravel_index(np.argmax(np.abs(M - M_ra)), M.shape)
        res['argmax_err_M'].append([int(i), int(j)])
        for t in T_sub:
            res['err_M_by_T'][str(t)].append(rel_max_err(M[:t, :t], M_ra[:t, :t]))
            res['err_Jr_by_T'][str(t)].append(rel_max_err(Jr[:t, :t], Jr_ra[:t, :t]))
        d = out.get('diag', {})
        res['kappa'].append(float(d.get('kappa', np.nan)))
        res['mass_constrained'].append(float(d.get('mass_constrained', np.nan)))
        res['ssj_forward_resid'].append(float(d.get('ssj_forward_resid', np.nan)))
    pos = [k for k, e in enumerate(res['eps']) if e > 0]
    eps_pos = [res['eps'][k] for k in pos]
    for key in ('err_M', 'err_Jr'):
        vals = [res[key][k] for k in pos]
        res[f'rates_{key}'] = pairwise_rates(eps_pos, vals)
        res[f'fitted_rate_{key}'] = fitted_rate(eps_pos, vals)
    res['err_M_by_window'] = res['err_M_by_T']
    res['err_Jr_by_window'] = res['err_Jr_by_T']
    res['v0_3'] = v03
    res['T'] = T
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T2: TANK limit
# ---------------------------------------------------------------------------

def check_tank_limit(provider, eps_list, T, dr, chis=(1.0, 1.5), spec=None, **kw):
    """TANK benchmark with an imposed Bilbiie-style income incidence chi (reduced form; chi is
    not derived from labour-market or profit-sharing primitives). Checks, as eps -> 0 for both
    types: (i) chi = 1: aggregate iMPC = lam I + (1-lam) M_RA; (ii) any chi: the GE-relevant
    matrix lam chi M_H + (1 - lam chi) M_S -> lam chi I + (1 - lam chi) M_RA; (iii) GE output
    for a fixed real-rate path = (1-lam)/(1-lam chi) x RANK. Conditions: same steady-state
    consumption and EIS, H types exactly hand-to-mouth in the limit, 1 - lam chi > 0."""
    t0 = time.time()
    res = dict(eps=[], err_MH_vs_I=[], err_MS_vs_MRA=[], err_Magg_vs_TANK=[],
               err_Magg_chi_vs_TANK={str(c): [] for c in chis},
               err_IRF_vs_TANK={str(c): [] for c in chis}, unconstrained_H=[], A_H=[],
               kappa_S=[], mass_constrained_S=[], ge_target_resid=[],
               amplification_HA={str(c): [] for c in chis})
    I = np.eye(T)
    v03 = _new_v03()
    for eps in eps_list:
        out = provider(eps, T, dr, chis=chis, **kw)
        _guard(v03, _v03_tank, out, spec)
        lam, s = out['lam'], out['ss']
        M_ra, _ = ra_jacobians(T, s['r'], s['eis'], s['C_S'], s['A_S'])
        M_agg = lam * out['M_H'] + (1 - lam) * out['M_S']
        res['eps'].append(float(eps))
        res['err_MH_vs_I'].append(rel_max_err(out['M_H'], I))
        res['err_MS_vs_MRA'].append(rel_max_err(out['M_S'], M_ra))
        res['err_Magg_vs_TANK'].append(rel_max_err(M_agg, tank_impc(T, s['r'], lam, 1.0)))
        for c in chis:
            # GE-relevant aggregate iMPC with income incidence chi: lam*chi*M_H + (1-lam*chi)*M_S
            M_chi = lam * c * out['M_H'] + (1 - lam * c) * out['M_S']
            res['err_Magg_chi_vs_TANK'][str(c)].append(
                rel_max_err(M_chi, tank_impc(T, s['r'], lam, c)))
            ref = tank_output_fixed_real_rate(dr, s['r'], s['eis'], lam, c, s['Y'])
            res['err_IRF_vs_TANK'][str(c)].append(rel_max_err(out['dY'][float(c)], ref))
        for c in chis:
            # impact output relative to the chi=1 (equal incidence) economy
            res['amplification_HA'][str(c)].append(
                float(out['dY'][float(c)][0] / out['dY'][float(chis[0])][0]))
        d = out.get('diag', {})
        for k in ('unconstrained_H', 'A_H', 'kappa_S', 'mass_constrained_S', 'ge_target_resid'):
            res[k].append(float(d.get(k, np.nan)))
    res['rates_err_MS_vs_MRA'] = pairwise_rates(res['eps'], res['err_MS_vs_MRA'])
    res['fitted_rate_err_MS_vs_MRA'] = fitted_rate(res['eps'], res['err_MS_vs_MRA'])
    # the H-type error vanishes exactly below a finite eps, so a single fitted slope for the
    # aggregate would be meaningless; report pairwise rates only
    res['rates_err_Magg_vs_TANK'] = pairwise_rates(res['eps'], res['err_Magg_vs_TANK'])
    for c in chis:
        e = res['err_IRF_vs_TANK'][str(c)]
        res[f'rates_err_IRF_chi{c}'] = pairwise_rates(res['eps'], e)
        res[f'fitted_rate_err_IRF_chi{c}'] = fitted_rate(res['eps'], e)
    res['amplification_factor'] = {str(c): (1 - lam) / (1 - lam * c) / ((1 - lam) / (1 - lam * chis[0]))
                                   for c in chis}
    res['lam'] = lam
    res['v0_3'] = v03
    res['T'] = T
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T3: Werning (2015) zero / vanishing liquidity
# ---------------------------------------------------------------------------

def check_werning(provider, b_list, T, dr, asymptotic_from=1e-3, spec=None, **kw):
    """For each b (borrowing limit -b, zero net supply, b -> 0):
      acyclical (income proportional to Y): error of dY_HA vs RANK  (Werning Prop. 2);
      countercyclical incidence: error vs the generalized zero-liquidity solution, and
      the (non-vanishing) deviation from RANK."""
    t0 = time.time()
    res = dict(b=[], err_acyc_vs_RANK=[], err_cyc_vs_ZL=[], dev_cyc_vs_RANK=[],
               zl_cyc_vs_RANK=[], r_gap_over_b_acyc=[], r_gap_over_b_cyc=[],
               mass_constrained=[], max_wealth_over_b=[])
    info_cyc = None
    v03 = _new_v03()
    for b in b_list:
        a = provider(b, T, dr, cyclical=False, **kw)
        c = provider(b, T, dr, cyclical=True, **kw)
        _guard(v03, _v03_werning, a, b, dr, False, spec)
        _guard(v03, _v03_werning, c, b, dr, True, spec)
        rank_a = ra_output_fixed_real_rate(dr, a['r_ss'], a['eis'])
        rank_c = ra_output_fixed_real_rate(dr, c['r_ss'], c['eis'])
        zl_c, info_cyc = zero_liquidity_output(dr, c['r_ss'], c['eis'], c['e_grid'], c['Pi'],
                                               c['zeta'])
        res['b'].append(float(b))
        res['err_acyc_vs_RANK'].append(rel_max_err(a['dY'], rank_a))
        res['err_cyc_vs_ZL'].append(rel_max_err(c['dY'], zl_c))
        res['dev_cyc_vs_RANK'].append(rel_max_err(c['dY'], rank_c))
        res['zl_cyc_vs_RANK'].append(rel_max_err(zl_c, rank_c))
        res['r_gap_over_b_acyc'].append(float(a['diag']['r_gap_over_b']))
        res['r_gap_over_b_cyc'].append(float(c['diag']['r_gap_over_b']))
        res['mass_constrained'].append(float(a['diag']['mass_constrained']))
        res['max_wealth_over_b'].append(float(a['diag']['max_wealth_over_b']))
    asym = [k for k, b in enumerate(res['b']) if b <= asymptotic_from]
    bs = [res['b'][k] for k in asym]
    for key in ('err_acyc_vs_RANK', 'err_cyc_vs_ZL'):
        res[f'rates_{key}'] = pairwise_rates(res['b'], res[key])
        res[f'fitted_rate_{key}_asymptotic'] = fitted_rate(bs, [res[key][k] for k in asym])
    res['asymptotic_from'] = asymptotic_from
    res['v0_3'] = v03
    res['zero_liquidity_cyc'] = {k: (float(v) if not isinstance(v, int) else v)
                                 for k, v in info_cyc.items()}
    res['T'] = T
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T4: Aiyagari steady state
# ---------------------------------------------------------------------------

def check_aiyagari(provider, eps_list, rho_limit, sd_list, rho_list, sd_for_rho, rho_for_sd,
                   spec=None, **kw):
    """No-idiosyncratic-risk limit (r -> 1/beta - 1 as risk -> 0; the code has no insurance
    markets, so this is not a complete-markets limit) and comparative statics in risk size and
    persistence for THIS calibration (regression checks, not general theorems).
    `all_below_rho` is mechanical when the provider searches r below 1/beta - 1 (the SSJ
    provider does) and is not counted as accuracy evidence. Non-mechanical diagnostics, if the
    provider reports them in diag: capital-market residual, Euler/KKT residuals, top-grid mass."""
    t0 = time.time()
    res = dict(eps=[], r=[], gap=[], mass_constrained=[], capital_market_resid_rel=[],
               euler_err_mean_weighted=[], euler_err_max_with_mass=[], kkt_violation_max_with_mass=[],
               mass_top_5pts=[], mass_policy_off_grid=[])
    v03 = _new_v03()
    for eps in eps_list:
        out = provider(eps, rho_limit, **kw)
        _guard(v03, _v03_aiyagari, out, spec)
        res['eps'].append(float(eps))
        res['r'].append(out['r'])
        res['gap'].append(out['rho_disc'] - out['r'])
        d = out['diag']
        res['mass_constrained'].append(d['mass_constrained'])
        ee = d.get('euler', {})
        res['capital_market_resid_rel'].append(float(d.get('capital_market_resid_rel', np.nan)))
        res['euler_err_mean_weighted'].append(float(ee.get('euler_err_mean_weighted', np.nan)))
        res['euler_err_max_with_mass'].append(float(ee.get('euler_err_max_with_mass', np.nan)))
        res['kkt_violation_max_with_mass'].append(float(ee.get('kkt_violation_max_with_mass', np.nan)))
        res['mass_top_5pts'].append(float(d.get('mass_top_5pts', np.nan)))
        res['mass_policy_off_grid'].append(float(d.get('mass_policy_off_grid', np.nan)))
    res['rho_disc'] = out['rho_disc']
    res['gap_over_eps2'] = [g / e ** 2 for g, e in zip(res['gap'], res['eps'])]
    res['rates_gap'] = pairwise_rates(res['eps'], res['gap'])
    res['fitted_rate_gap'] = fitted_rate(res['eps'], res['gap'])
    res['r_by_sd'] = dict(sd=list(map(float, sd_list)), rho=rho_for_sd,
                          r=[provider(s, rho_for_sd, **kw)['r'] for s in sd_list])
    res['r_by_rho'] = dict(rho=list(map(float, rho_list)), sd=sd_for_rho,
                           r=[provider(sd_for_rho, rh, **kw)['r'] for rh in rho_list])
    allr = res['r'] + res['r_by_sd']['r'] + res['r_by_rho']['r']
    res['all_below_rho'] = bool(all(r < res['rho_disc'] for r in allr))
    res['all_below_rho_is_mechanical'] = True     # SSJ provider: root bracket below 1/beta - 1
    res['monotone_in_sd'] = bool(np.all(np.diff(res['r_by_sd']['r']) < 0))
    res['monotone_in_rho'] = bool(np.all(np.diff(res['r_by_rho']['r']) < 0))
    res['v0_3'] = v03
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T5: intertemporal Keynesian cross
# ---------------------------------------------------------------------------

def ikc_rhs(J, dG, dTax, dr):
    """(I - M)^{-1} [dG + J^{C,Tax} dTax + J^{C,r} dr] with M = J^{C,Y}."""
    M = J['C']['Y']
    T = M.shape[0]
    partial = dG + J['C']['Tax'] @ dTax + J['C']['r'] @ dr
    return np.linalg.solve(np.eye(T) - M, partial), partial


def check_ikc(provider, T_list, shocks_fn, spec=None, **kw):
    """shocks_fn(T) -> {label: (kind, {input: path})}. For each T and shock, compare the
    provider's GE output path with the IKC formula built from its own block Jacobians."""
    t0 = time.time()
    res = dict(T=[], err={}, walras_asset_mkt={})
    v03 = _new_v03()
    for T in T_list:
        shocks = shocks_fn(T)
        out = provider(T, shocks, **kw)
        _guard(v03, _v03_fiscal, out, shocks, spec)
        res['T'].append(int(T))
        for label, irf in out['irfs'].items():
            rhs, _ = ikc_rhs(out['J'], irf['dG'], irf['dTax'], irf['dr'])
            res['err'].setdefault(label, []).append(rel_max_err(irf['dY'], rhs))
            res['walras_asset_mkt'].setdefault(label, []).append(
                float(np.max(np.abs(irf['asset_mkt']))))
    res['max_err'] = float(max(max(v) for v in res['err'].values()))
    res['v0_3'] = v03
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T6: HANK - RANK differencing identity and decomposition
# ---------------------------------------------------------------------------

def check_differencing(J, ss, irf, n_dates=12):
    """Identity (I-M_HA)^{-1} - (I-M_RA)^{-1} = (I-M_HA)^{-1}(M_HA-M_RA)(I-M_RA)^{-1}, and
    decomposition dY_HA - dY_RA = (I-M_HA)^{-1}[(M_HA-M_RA)(dY_RA - dTax) + (Jr_HA-Jr_RA) dr]
    for the same steady state (r, eis, C, A) and the same policy paths (dG, dTax, dr)."""
    M = J['C']['Y']
    T = M.shape[0]
    I = np.eye(T)
    M_ra, Jr_ra = ra_jacobians(T, ss['r'], ss['eis'], ss['C'], ss['A'])
    inv_ha = np.linalg.inv(I - M)
    inv_ra = np.linalg.inv(I - M_ra)
    lhs = inv_ha - inv_ra
    rhs = inv_ha @ (M - M_ra) @ inv_ra
    dG, dTax, dr = irf['dG'], irf['dTax'], irf['dr']
    dY_ha = irf['dY']
    dY_ra = inv_ra @ (dG - M_ra @ dTax + Jr_ra @ dr)
    first_round = (M - M_ra) @ (dY_ra - dTax) + (J['C']['r'] - Jr_ra) @ dr
    dY_gap = inv_ha @ first_round
    q = 1 / (1 + ss['r'])
    disc = q ** np.arange(T)
    out = dict(
        identity_err=rel_max_err(rhs, lhs),
        cond_I_minus_M_RA=float(np.linalg.cond(I - M_ra)),
        cond_I_minus_M_HA=float(np.linalg.cond(I - M)),
        decomposition_err=rel_max_err(dY_gap, dY_ha - dY_ra),
        dates=list(range(n_dates)),
        dY_HA=_tolist(dY_ha[:n_dates]), dY_RA=_tolist(dY_ra[:n_dates]),
        dY_gap=_tolist(dY_gap[:n_dates]),
        first_round=_tolist(first_round[:n_dates]),
        ge_propagation=_tolist((dY_gap - first_round)[:n_dates]),
        cum_multiplier_HA=float(disc @ dY_ha / (disc @ dG)) if np.any(dG) else None,
        cum_multiplier_RA=float(disc @ dY_ra / (disc @ dG)) if np.any(dG) else None,
        impact_multiplier_HA=float(dY_ha[0] / dG[0]) if dG[0] != 0 else None,
        impact_multiplier_RA=float(dY_ra[0] / dG[0]) if dG[0] != 0 else None,
    )
    out['_arrays'] = dict(inv_ha=inv_ha, M_ra=M_ra, dY_ra=dY_ra, x=dY_ra - dTax)
    return out
