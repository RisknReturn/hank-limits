"""Canonical settings of the suite (v0.2; unchanged from v0.1), shared by tests/ and scripts/convergence_report.py.
v0.3: the checks also receive the documented specification (`hank_limits.spec`) for the V checks."""
from __future__ import annotations

import time
import numpy as np

from . import models as M
from . import checks as C
from . import spec as SPEC

T1 = dict(eps=[0.0, 0.2, 0.1, 0.05, 0.025, 0.0125], T=200, T_sub=(50, 100))
T2 = dict(eps=[0.2, 0.1, 0.05, 0.025, 0.0125], T=200, chis=(1.0, 1.5), rho_r=0.8, dr0=-0.01)
T3 = dict(b=[1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4], T=150, rho_r=0.8, dr0=-0.01,
          asymptotic_from=1e-3)
T4 = dict(eps=[0.4, 0.2, 0.1, 0.05, 0.025], rho_limit=0.9, sd_list=[0.1, 0.2, 0.3, 0.4],
          rho_list=[0.0, 0.3, 0.6, 0.9], sd_for_rho=0.4, rho_for_sd=0.9)
T5 = dict(T_list=[100, 200, 300], rho_shock=0.8, size=0.01)
T6 = dict(T_list=[150, 300], T=300, rho_shock=0.8, size=0.01, n_dates=12)


def real_rate_path(T, dr0, rho):
    return dr0 * rho ** np.arange(T)


def run_t1(**kw):
    return C.check_rank_limit(M.ssj_rank_limit_provider, T1['eps'], T1['T'],
                              T_sub=T1['T_sub'], spec=SPEC.T1, **kw)


def run_t2(**kw):
    dr = real_rate_path(T2['T'], T2['dr0'], T2['rho_r'])
    return C.check_tank_limit(M.ssj_tank_provider, T2['eps'], T2['T'], dr,
                              chis=T2['chis'], spec=SPEC.T2, **kw)


def run_t3(**kw):
    dr = real_rate_path(T3['T'], T3['dr0'], T3['rho_r'])
    return C.check_werning(M.ssj_werning_provider, T3['b'], T3['T'], dr,
                           asymptotic_from=T3['asymptotic_from'], spec=SPEC.T3, **kw)


def run_t4(**kw):
    return C.check_aiyagari(M.ssj_aiyagari_r, T4['eps'], T4['rho_limit'], T4['sd_list'],
                            T4['rho_list'], T4['sd_for_rho'], T4['rho_for_sd'], spec=SPEC.T4,
                            **kw)


def fiscal_shocks(T):
    t = np.arange(T)
    path = T5['size'] * T5['rho_shock'] ** t
    return {'G_debt_fixed_r': ('fixed_r', {'G': path}),
            'r_fixed_r': ('fixed_r', {'r': -path}),
            'G_debt_r_rule': ('r_rule', {'G': path}),
            'rshock_r_rule': ('r_rule', {'rshock': -path})}


def run_t5():
    """Three solver settings: the default-consistent SSJ solve, SSJ with refinement, and
    SSJ with Js/options mismatch (documents finding F1)."""
    out = {}
    for label, kw in [('ssj_consistent', dict()),
                      ('ssj_consistent_refined', dict(solver='refined')),
                      ('ssj_inconsistent_options', dict(ssj_consistent=False))]:
        out[label] = C.check_ikc(M.ssj_fiscal_provider, T5['T_list'], fiscal_shocks,
                                 spec=SPEC.T5, **kw)
    return out


def _wealth_quartile_masks(ss):
    da = ss.internals['hh']['Dbeg'].sum(0)
    cdf = np.cumsum(da)
    grp = np.searchsorted([0.25, 0.5, 0.75], cdf - da / 2)
    return [(grp == k).astype(float) for k in range(4)]


def run_t6():
    t0 = time.time()
    res = dict(by_T={})
    for T in T6['T_list']:
        t = np.arange(T)
        G = T6['size'] * T6['rho_shock'] ** t
        out = M.ssj_fiscal_provider(T, {'G': ('fixed_r', {'G': G}),
                                        'r': ('fixed_r', {'r': -G})}, solver='refined')
        dG = C.check_differencing(out['J'], out['ss'], out['irfs']['G'], T6['n_dates'])
        dR = C.check_differencing(out['J'], out['ss'], out['irfs']['r'], T6['n_dates'])
        # wealth-quartile split of the first-round gap (G shock, dr = 0)
        m = out['model']
        masks = _wealth_quartile_masks(m['ss'])
        arr = dG.pop('_arrays')
        dR.pop('_arrays')
        x, M_ra, inv_ha = arr['x'], arr['M_ra'], arr['inv_ha']
        mus, paths = M.cohort_impc_paths(m, x, masks)
        Mx = out['J']['C']['Y'] @ x
        g = [mu * (p - M_ra @ x) for mu, p in zip(mus, paths)]
        q = (1 / (1 + out['ss']['r'])) ** t
        a_grid = m['a_grid']
        groups = []
        for k in range(4):
            sel = masks[k] > 0
            groups.append(dict(
                mass=mus[k], wealth_min=float(a_grid[sel].min()),
                wealth_max=float(a_grid[sel].max()),
                first_round=[float(v) for v in g[k][:T6['n_dates']]],
                pv_first_round=float(q @ g[k]),
                ge_image_impact=float((inv_ha @ g[k])[0])))
        res['by_T'][str(T)] = dict(
            G=dG, r=dR, groups=groups,
            group_sum_vs_M_x=float(np.max(np.abs(sum(mu * p for mu, p in zip(mus, paths)) - Mx))
                                   / np.max(np.abs(Mx))),
            group_first_round_sum_vs_total=float(
                np.max(np.abs(sum(g)[:T6['n_dates']] - np.array(dG['first_round'])))
                / np.max(np.abs(dG['first_round']))),
            pv_x=float(q @ x), fiscal_diag=out['model']['diag'])
    res['runtime_s'] = time.time() - t0
    return res
