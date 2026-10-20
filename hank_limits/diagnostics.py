"""Numerical diagnostics added in v0.2 (R4 review).

* `euler_kkt_residuals`: unit-free Euler errors on unconstrained states and KKT violations on
  constrained states of a one-asset household policy (SSJ timing).
* `distribution_audit`: audit of FastHetBlock's stationary distribution: inverse-iteration
  report, sensitivity to the shift tau and to the number of sweeps, and a comparison with
  SSJ's own forward iteration (`HetBlock.forward_steady_state`) on the same policy.

A small stationarity residual alone does not bound the distribution error when the chain
mixes slowly; these diagnostics report what else was checked.
"""
from __future__ import annotations

import time
import numpy as np
from sequence_jacobian.blocks.het_block import HetBlock
from sequence_jacobian.utilities.optimized_routines import within_tolerance

from .core import stationary_distribution


def euler_kkt_residuals(Pi, a_grid, a_pol, c_pol, D, r, beta, eis, bind_tol=1e-10,
                        mass_tol=1e-12):
    """Euler/KKT residuals of a steady-state policy (state = (e, beginning-of-period a)).

    Euler: c^{-1/eis} = beta (1+r) E[c(e', a')^{-1/eis} | e], with c(e', .) linearly
    interpolated on the asset grid. c_impl = (beta (1+r) E[.])^{-eis}; err = c_impl/c - 1.
    Unconstrained states (a' > a_min + bind_tol): |err| is the unit-free Euler error.
    Constrained states: KKT requires u'(c) >= beta (1+r) E u'(c'), i.e. err >= 0; only
    max(0, -err) is a violation (a positive Euler gap at a binding state is NOT an error).
    Policies above the top grid point are extrapolated linearly and counted separately.
    """
    n_e, n_a = a_pol.shape
    up = np.empty_like(c_pol)
    for e in range(n_e):
        # c(e', a'(e, .)) for every e', linear inter/extrapolation in a
        cn = np.empty((n_e, n_a))
        for e2 in range(n_e):
            x = a_pol[e]
            j = np.clip(np.searchsorted(a_grid, x) - 1, 0, n_a - 2)
            w = (x - a_grid[j]) / (a_grid[j + 1] - a_grid[j])
            cn[e2] = (1 - w) * c_pol[e2, j] + w * c_pol[e2, j + 1]
        up[e] = Pi[e] @ np.maximum(cn, 1e-300) ** (-1 / eis)
    c_impl = (beta * (1 + r) * up) ** (-eis)
    err = c_impl / c_pol - 1
    binding = a_pol <= a_grid[0] + bind_tol
    off_top = a_pol > a_grid[-1]
    has_mass = D > mass_tol
    unc = ~binding & ~off_top
    absu = np.abs(err[unc])
    w = D[unc]
    viol = np.maximum(0.0, -err[binding])
    out = dict(
        euler_err_mean_weighted=float((absu * w).sum() / max(w.sum(), 1e-300)),
        euler_err_max_with_mass=float(np.abs(err[unc & has_mass]).max()) if np.any(unc & has_mass) else 0.0,
        euler_err_max_all_grid=float(absu.max()) if absu.size else 0.0,
        log10_euler_err_mean_weighted=float(np.log10(max((absu * w).sum() / max(w.sum(), 1e-300), 1e-300))),
        kkt_violation_max=float(viol.max()) if viol.size else 0.0,
        kkt_violation_max_with_mass=float(np.maximum(0.0, -err[binding & has_mass]).max())
        if np.any(binding & has_mass) else 0.0,
        share_binding=float(D[binding].sum()),
        mass_policy_above_grid=float(D[off_top].sum()),
        n_unconstrained=int(unc.sum()), n_binding=int(binding.sum()))
    return out


def _moments(Dbeg, st):
    D = st['Pi'].T @ Dbeg
    return dict(A=float(np.vdot(D, st['a'])), C=float(np.vdot(D, st['c'])))


def _ssj_forward(block, st, tol, maxit):
    """Mirror of SSJ 1.0.0 HetBlock.forward_steady_state with an iteration counter (the
    result is checked against a direct call of SSJ's method)."""
    st = {k: v for k, v in st.items() if k not in ('Dbeg', 'D')}
    t0 = time.time()
    Dbeg_ssj, _ = HetBlock.forward_steady_state(block, st, tol=tol, maxit=maxit)
    sec = time.time() - t0
    exog = block.make_exog_law_of_motion(st)
    endog = block.make_endog_law_of_motion(st)
    pis = [exog[0].stationary(None)]
    Dbeg = np.outer(pis[0], np.full(len(st['a_grid']), 1 / len(st['a_grid'])))
    D = exog.forward(Dbeg)
    n_it = None
    for it in range(maxit):
        Dbeg_new = endog.forward(D)
        D_new = exog.forward(Dbeg_new)
        if it % 10 == 0 and within_tolerance(Dbeg, Dbeg_new, tol):
            n_it = it
            break
        Dbeg, D = Dbeg_new, D_new
    same = bool(np.array_equal(Dbeg, Dbeg_ssj))
    return Dbeg_ssj, dict(seconds=sec, iterations=n_it, mirror_identical=same, tol=tol)


def distribution_audit(block, ss, name='hh', taus=(1e-12, 1e-10, 1e-8, 1e-6),
                       sweeps=(1, 2, 3), ssj_tols=(1e-10, 1e-12), ssj_maxit=2_000_000,
                       threshold=None):
    """Audit FastHetBlock's stationary distribution at steady state `ss` of `block`.

    Reference = default settings (tau=1e-10, up to 20 sweeps, 5 polishing steps).
    Reports, relative to the reference: max |D - D_ref|, L1 |D - D_ref|, and the change in
    mean assets A and consumption C, for (i) other shifts tau, (ii) fewer sweeps without
    polishing, (iii) SSJ's own forward iteration from its default seed at several
    tolerances. `threshold`: also report the mass with beginning-of-period assets below it.
    """
    st = {**ss.toplevel, **ss.internals[name]}
    Pi, a_pol, grid = st['Pi'], st['a'], st['a_grid']
    t0 = time.time()
    D0, info0 = stationary_distribution(Pi, a_pol, grid, return_info=True)
    sec0 = time.time() - t0
    m0 = _moments(D0, st)

    def cmp(D):
        m = _moments(D, st)
        out = dict(max_abs_diff=float(np.abs(D - D0).max()), L1_diff=float(np.abs(D - D0).sum()),
                   dA_rel=float(abs(m['A'] - m0['A']) / max(abs(m0['A']), 1e-12)),
                   dC_rel=float(abs(m['C'] - m0['C']) / abs(m0['C'])))
        if threshold is not None:
            out['d_mass_below_threshold'] = float(abs(D[:, grid < threshold].sum()
                                                      - D0[:, grid < threshold].sum()))
        return out

    res = dict(reference=dict(info0, seconds=sec0, **m0), shift={}, sweeps={}, ssj_forward={})
    if threshold is not None:
        res['reference']['mass_below_threshold'] = float(D0[:, grid < threshold].sum())
        res['threshold'] = float(threshold)
    for tau in taus:
        D, info = stationary_distribution(Pi, a_pol, grid, tau=tau, return_info=True)
        res['shift'][f'{tau:.0e}'] = dict(cmp(D), sweeps=info['sweeps'],
                                          sweep_tol_met=info['sweep_tol_met'],
                                          neg_mass_rel_max=info['neg_mass_rel_max'])
    for n in sweeps:
        D, info = stationary_distribution(Pi, a_pol, grid, maxit=n, n_polish=0, return_info=True)
        res['sweeps'][str(n)] = dict(cmp(D), last_step_max=info['last_step_max'],
                                     resid_max=info['resid_max'])
    for tol in ssj_tols:
        try:
            D, info = _ssj_forward(block, st, tol, ssj_maxit)
            res['ssj_forward'][f'{tol:.0e}'] = dict(cmp(D), **info)
        except ValueError as exc:            # SSJ raises if maxit is reached
            res['ssj_forward'][f'{tol:.0e}'] = dict(error=str(exc), tol=tol)
    return res
