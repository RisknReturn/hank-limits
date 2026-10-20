"""Core building blocks shared by all limit tests.

* `make_household`: the standard SSJ one-asset household (identical backward step to
  `sequence_jacobian.hetblocks.hh_sim.hh`, with an initial guess that is valid for
  negative asset grids, e.g. near-natural borrowing limits).
* `FastHetBlock`: a thin subclass of SSJ's `HetBlock` that only changes *how the
  stationary distribution is found* (sparse inverse iteration on the eigenvalue-one
  problem of the Young-lottery chain instead of forward iteration) and warm-starts
  the backward iteration.
  The Jacobian algorithm (fake-news), the GE solvers and everything else are SSJ's own.

  Why: near the representative-agent limit the stationary distribution mixes at a
  rate 1 - O(eps^2) per period, so forward iteration would need >10^5 sweeps.
"""
from __future__ import annotations

import copy
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import sequence_jacobian as sj
from sequence_jacobian.blocks.het_block import HetBlock
from sequence_jacobian import interpolate, misc
from sequence_jacobian.utilities.interpolate import interpolate_coord_robust


# ---------------------------------------------------------------------------
# Household block
# ---------------------------------------------------------------------------

def _hh_init(a_grid, y, r, eis):
    """Initial guess for Va valid also when a_grid[0] < 0 (cash on hand above the limit)."""
    coh = (1 + r) * a_grid[np.newaxis, :] + y[:, np.newaxis]
    c0 = 0.2 * (coh - a_grid[0]) + 1e-3
    Va = (1 + r) * c0 ** (-1 / eis)
    return Va


# SSJ parses the *source* of the backward function for its return list, so we need a
# function whose final line is literally `return Va, a, c`.
def _make_hh_function():
    def hh(Va_p, a_grid, y, r, beta, eis):
        uc_nextgrid = beta * Va_p
        c_nextgrid = uc_nextgrid ** (-eis)
        coh = (1 + r) * a_grid[np.newaxis, :] + y[:, np.newaxis]
        a = interpolate.interpolate_y(c_nextgrid + a_grid, coh, a_grid)
        misc.setmin(a, a_grid[0])
        c = coh - a
        Va = (1 + r) * c ** (-1 / eis)
        return Va, a, c
    return hh


class FastHetBlock(HetBlock):
    """SSJ HetBlock with (i) stationary distribution by sparse inverse iteration,
    (ii) warm-started backward iteration. Diagnostics of the last steady state are
    stored in `_ss_diag` (SSJ forward-operator residual, mass whose policy leaves the
    grid, mass in the top 5 grid points).

    Optional instance attribute `_Dbeg_override`: if set (array of shape (n_e, n_a)),
    it is used as the 'stationary' beginning-of-period distribution. This is only
    legitimate when the steady-state policy leaves it invariant (checked).
    """

    def initialize_backward(self, ss):
        key = self.backward[0]
        cache = self.__dict__.get('_va_cache')
        shape = (len(ss[self.exogenous[0]]), len(ss[self.policy[0] + '_grid']))
        if key not in ss and cache is not None and cache.shape == shape:
            ss[key] = cache.copy()
        else:
            super().initialize_backward(ss)

    def backward_steady_state(self, ss, tol=1e-8, maxit=5000):
        out = super().backward_steady_state(ss, tol=tol, maxit=maxit)
        self._va_cache = np.array(out[self.backward[0]], copy=True)
        return out

    def forward_steady_state(self, ss, tol=1e-10, maxit=100_000):
        exog = self.make_exog_law_of_motion(ss)
        endog = self.make_endog_law_of_motion(ss)
        override = self.__dict__.get('_Dbeg_override')
        grid = ss[self.policy[0] + '_grid']
        a_pol = ss[self.policy[0]]
        info = {}
        if override is not None:
            Dbeg = np.array(override, dtype=float)
        else:
            Pi = ss[self.exogenous[0]]
            Dbeg, info = stationary_distribution(Pi, a_pol, grid, return_info=True)
        # diagnostics: residual under SSJ's own (unclipped) forward operator, and the
        # mass whose savings policy leaves the grid (where the solver clips the lottery)
        fwd = endog.forward(exog.forward(Dbeg)) - Dbeg
        resid = float(np.max(np.abs(fwd)))
        D = exog.forward(Dbeg)
        outside = float(np.sum(D[(a_pol > grid[-1]) | (a_pol < grid[0])]))
        self._ss_diag = dict(ssj_forward_resid=resid, ssj_forward_resid_L1=float(np.abs(fwd).sum()),
                             mass_policy_off_grid=outside,
                             mass_top_5pts=float(Dbeg[:, -5:].sum()),
                             **{f'inv_iter_{k}': v for k, v in info.items()})
        if override is not None and resid > 1e-8:
            raise ValueError(f'_Dbeg_override is not invariant (err={resid:.2e})')
        return Dbeg, D


def make_fast(block: HetBlock) -> FastHetBlock:
    """Return a copy of an SSJ HetBlock whose class is FastHetBlock."""
    b = copy.copy(block)
    b.__class__ = FastHetBlock
    return b


def make_household(hetinputs=None, name='hh') -> FastHetBlock:
    """One-asset household HetBlock with optional hetinputs (e.g. income)."""
    f = _make_hh_function()
    blk = sj.het(exogenous='Pi', policy='a', backward='Va', backward_init=_hh_init)(f)
    if hetinputs:
        blk = blk.add_hetinputs(hetinputs)
    blk = make_fast(blk)
    if name != 'hh':
        blk = blk.rename(name)
    return blk


# ---------------------------------------------------------------------------
# Stationary distribution by sparse direct solve
# ---------------------------------------------------------------------------

def transition_matrix(Pi, a_pol, grid):
    """Sparse matrix L such that vec(Dbeg') = L vec(Dbeg).

    SSJ timing: D = Pi' Dbeg (exogenous shock), then Dbeg' = lottery(a_pol) D.
    Policies outside the grid are clipped (reflecting barriers) so that L is a proper
    stochastic matrix; the mass affected is reported in FastHetBlock._ss_diag.
    """
    n_e, n_a = a_pol.shape
    i, pi = interpolate_coord_robust(grid, np.clip(a_pol, grid[0], grid[-1]))
    i = i.astype(np.int64)
    pi = np.clip(pi, 0.0, 1.0)
    rows_lo = (np.arange(n_e)[:, None] * n_a + i).ravel()
    rows_hi = rows_lo + 1
    cols = np.arange(n_e * n_a)
    lot = sp.csr_matrix((np.concatenate([pi.ravel(), 1 - pi.ravel()]),
                         (np.concatenate([rows_lo, rows_hi]), np.concatenate([cols, cols]))),
                        shape=(n_e * n_a, n_e * n_a))
    X = sp.kron(sp.csr_matrix(Pi.T), sp.identity(n_a), format='csr')
    return (lot @ X).tocsc()


def stationary_distribution(Pi, a_pol, grid, tau=1e-10, maxit=20, tol=1e-15, n_polish=5,
                            return_info=False):
    """Stationary distribution of the (clipped) Young-lottery chain by inverse iteration.

    Solves ((1 + tau) I - L) x_{n+1} = x_n with one sparse LU factorisation. Each sweep
    shrinks every non-stationary component by tau / (1 - lambda_2 + tau), so a handful
    of sweeps reach machine precision even when the chain mixes extremely slowly
    (1 - lambda_2 ~ 1e-6 near the RA limit). States are ordered asset-major
    (j * n_e + e) so the matrix is banded and the LU is cheap.

    A small stationarity residual does NOT by itself bound the distribution error when the
    chain mixes slowly (the inverse of I - L on the mass-zero subspace grows like
    1/(1 - lambda_2)); `return_info=True` therefore also reports what the iteration did:
    negative mass before the absolute value is taken (relative to total absolute mass), the
    number of sweeps, whether the sweep tolerance was met, the size of the last update, the
    mass drift during polishing, and L1 / max residuals of the result under the (clipped)
    chain. (v0.2.1: the v0.2 field `normalisation_err_last_sweep` was identically zero by
    construction and has been removed; negative mass is reported by `neg_mass_rel_max`.) See README section 3 and `hank_limits.diagnostics`.
    """
    n_e, n_a = a_pol.shape
    N = n_e * n_a
    L = transition_matrix(Pi, a_pol, grid)
    perm = (np.arange(n_e)[None, :] * n_a + np.arange(n_a)[:, None]).ravel()  # new -> old
    Lp = L[perm][:, perm].tocsc()
    lu = spla.splu(((1.0 + tau) * sp.identity(N, format='csc') - Lp).tocsc(), permc_spec='NATURAL')
    dp = np.full(N, 1.0 / N)
    neg_mass_max, step, done, it = 0.0, np.nan, False, 0
    for it in range(1, maxit + 1):
        new = lu.solve(dp)
        tot = np.abs(new).sum()
        neg_mass_max = max(neg_mass_max, float(-new[new < 0].sum() / tot))
        new = np.abs(new)
        new /= new.sum()
        step = float(np.max(np.abs(new - dp)))
        done = step < tol
        dp = new
        if done:
            break
    for _ in range(n_polish):   # polish with the chain itself
        dp = Lp @ dp
    pre = float(abs(dp.sum() - 1.0))
    dp /= dp.sum()
    d = np.empty(N)
    d[perm] = dp
    d = d.reshape(n_e, n_a)
    if not return_info:
        return d
    res = Lp @ dp - dp
    info = dict(tau=tau, sweeps=it, sweep_tol_met=bool(done), last_step_max=step,
                neg_mass_rel_max=neg_mass_max,
                mass_drift_polish=pre, resid_L1=float(np.abs(res).sum()),
                resid_max=float(np.abs(res).max()), min_mass=float(d.min()))
    return d, info


# ---------------------------------------------------------------------------
# Grids / income processes
# ---------------------------------------------------------------------------

def rouwenhorst(rho, sd, n):
    """SSJ Rouwenhorst discretisation; sd = unconditional sd of log income.

    For sd == 0 returns a degenerate (all-ones) income grid with the same Pi.
    """
    if sd == 0:
        e, pi, Pi = sj.grids.markov_rouwenhorst(rho, 1.0, n)
        return np.ones(n), pi, Pi
    return sj.grids.markov_rouwenhorst(rho, sd, n)


def asset_grid(a_min, a_max, n, include=None):
    """SSJ double-exponential grid on [a_min, a_max]; optionally force a point in the grid."""
    g = sj.grids.asset_grid(a_min, a_max, n)
    if include is not None:
        j = np.argmin(np.abs(g - include))
        g[j] = include
        assert np.all(np.diff(g) > 0)
    return g


def income_basic(Y, e_grid):
    y = Y * e_grid
    return y


def mass_at_constraint(ss, a_grid, block_name='hh', tol=1e-10):
    """Share of households whose policy a' is at the borrowing limit a_grid[0]."""
    inter = ss.internals[block_name]
    return float(np.sum(inter['D'][inter['a'] <= a_grid[0] + tol]))
