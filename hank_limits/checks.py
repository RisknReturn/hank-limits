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
                 tank_impc, zero_liquidity_output)


def _tolist(x):
    return [float(v) for v in np.asarray(x).ravel()]


# ---------------------------------------------------------------------------
# T1: RANK limit
# ---------------------------------------------------------------------------

def check_rank_limit(provider, eps_list, T, T_sub=(50, 100), **kw):
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
    for eps in eps_list:
        out = provider(eps, T, **kw)
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
    res['T'] = T
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T2: TANK limit
# ---------------------------------------------------------------------------

def check_tank_limit(provider, eps_list, T, dr, chis=(1.0, 1.5), **kw):
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
    for eps in eps_list:
        out = provider(eps, T, dr, chis=chis, **kw)
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
    res['T'] = T
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T3: Werning (2015) zero / vanishing liquidity
# ---------------------------------------------------------------------------

def check_werning(provider, b_list, T, dr, asymptotic_from=1e-3, **kw):
    """For each b (borrowing limit -b, zero net supply, b -> 0):
      acyclical (income proportional to Y): error of dY_HA vs RANK  (Werning Prop. 2);
      countercyclical incidence: error vs the generalized zero-liquidity solution, and
      the (non-vanishing) deviation from RANK."""
    t0 = time.time()
    res = dict(b=[], err_acyc_vs_RANK=[], err_cyc_vs_ZL=[], dev_cyc_vs_RANK=[],
               zl_cyc_vs_RANK=[], r_gap_over_b_acyc=[], r_gap_over_b_cyc=[],
               mass_constrained=[], max_wealth_over_b=[])
    info_cyc = None
    for b in b_list:
        a = provider(b, T, dr, cyclical=False, **kw)
        c = provider(b, T, dr, cyclical=True, **kw)
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
    res['zero_liquidity_cyc'] = {k: (float(v) if not isinstance(v, int) else v)
                                 for k, v in info_cyc.items()}
    res['T'] = T
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T4: Aiyagari steady state
# ---------------------------------------------------------------------------

def check_aiyagari(provider, eps_list, rho_limit, sd_list, rho_list, sd_for_rho, rho_for_sd,
                   **kw):
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
    for eps in eps_list:
        out = provider(eps, rho_limit, **kw)
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


def check_ikc(provider, T_list, shocks_fn, **kw):
    """shocks_fn(T) -> {label: (kind, {input: path})}. For each T and shock, compare the
    provider's GE output path with the IKC formula built from its own block Jacobians."""
    t0 = time.time()
    res = dict(T=[], err={}, walras_asset_mkt={})
    for T in T_list:
        out = provider(T, shocks_fn(T), **kw)
        res['T'].append(int(T))
        for label, irf in out['irfs'].items():
            rhs, _ = ikc_rhs(out['J'], irf['dG'], irf['dTax'], irf['dr'])
            res['err'].setdefault(label, []).append(rel_max_err(irf['dY'], rhs))
            res['walras_asset_mkt'].setdefault(label, []).append(
                float(np.max(np.abs(irf['asset_mkt']))))
    res['max_err'] = float(max(max(v) for v in res['err'].values()))
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
