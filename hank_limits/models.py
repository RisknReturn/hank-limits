"""SSJ model builders: the reference implementation ("adapter") of every limit test.

Each `ssj_*` function returns plain numpy objects in the format documented in
`hank_limits.interface`, so the model-agnostic checks in `hank_limits.checks` can be run
on any other HANK code that produces the same objects.

Calibrations are annual. All SSJ household Jacobians are computed with two-sided finite
differences (SSJ's default one-sided step h=1e-4 leaves an O(h) error of ~5e-5 in dC/dr,
which would put a floor under the convergence curves).
"""
from __future__ import annotations

import copy
import numpy as np
import scipy.linalg as sla
from scipy.optimize import brentq

import sequence_jacobian as sj

from .core import (make_household, income_basic, rouwenhorst, asset_grid,
                   mass_at_constraint)

SS_OPTS = dict(backward_tol=1e-12, backward_maxit=100_000)


def _het_opts(names, consistent=True, h=1e-4):
    """Options making SSJ's HetBlock.impulse_linear use the same differentiation as the
    supplied Js. (SSJ 1.0.0: HetBlock._impulse_linear ignores Js and re-differentiates
    with its own defaults h=1e-4, one-sided; see README, finding F1.)"""
    if not consistent:
        return {}
    return {n: {'twosided': True, 'h': h} for n in names}


def solve_impulse_linear_refined(model, ss, unknowns, targets, inputs, Js, options,
                                 n_refine=2):
    """Same algorithm as SSJ's Block.solve_impulse_linear (H_U dU = -dH, then propagate),
    but with iterative refinement of the LU solve. Needed because Gaussian elimination
    with partial pivoting can show large element growth on stacked sequence-space
    systems (README, finding F2)."""
    inputs = sj.ImpulseDict(inputs)
    T = inputs.T
    H_U = model.jacobian(ss, unknowns, targets, T, Js, options).pack(T)
    dH = model.impulse_linear(ss, inputs, targets, Js, options).get(targets).pack()
    lu = sla.lu_factor(H_U)
    x = sla.lu_solve(lu, -dH)
    for _ in range(n_refine):
        x = x + sla.lu_solve(lu, -dH - H_U @ x)
    dU = sj.ImpulseDict.unpack(x, unknowns, T)
    return (inputs | dU) | model.impulse_linear(ss, dU | inputs, None, Js, options)


# v0.3: steady-state primitives for the checks in `hank_limits.checks` (optional provider
# fields `ss_detail`; read from SSJ's steady-state dict, nothing is recomputed).

def _household_detail(ss, name='hh', suffix=''):
    """C, A, beta, borrowing limit, income chain and income by state of one household block."""
    inter = ss.internals[name]
    d = dict(C=float(ss['C' + suffix]), A=float(ss['A' + suffix]),
             beta=float(ss['beta' + suffix]), a_min=float(np.asarray(ss['a_grid' + suffix])[0]),
             e_grid=np.asarray(ss['e_grid' + suffix]), Pi=np.asarray(inter['Pi']),
             y=np.asarray(inter['y']))
    return {k + suffix: v for k, v in d.items()}


# ---------------------------------------------------------------------------
# Income hetinputs
# ---------------------------------------------------------------------------

def income_cyclical(Y, Yss, e_grid, zeta):
    y = Y * e_grid * (Y / Yss) ** zeta
    return y


def income_wage(w, e_grid):
    y = w * e_grid
    return y


def income_aftertax(Y, Tax, e_grid):
    y = (Y - Tax) * e_grid
    return y


# ---------------------------------------------------------------------------
# T1 / T2 savers: permanent-income limit household
# ---------------------------------------------------------------------------

PIH_DEFAULTS = dict(r=0.03, eis=0.5, A_target=5.0, n_e=3, rho_e=0.9,
                    n_a=4800, a_max=250.0, nbl_frac=0.99, Y=1.0)


def pih_household(eps, **kw):
    """HA household with idiosyncratic log-income sd `eps` and no binding constraint.

    * borrowing limit a_min = nbl_frac * (natural borrowing limit -y_min/r);
    * beta calibrated so that mean assets equal A_target (a stationary distribution
      with risk requires beta (1+r) < 1; we report kappa = 1 - beta (1+r) -> 0);
    * eps == 0: exact permanent-income consumer, beta (1+r) = 1, all mass at A_target
      (any distribution is stationary then; the point mass is a legitimate choice).
    Returns dict(block, ss, calib, a_grid, kappa, diag).
    """
    p = {**PIH_DEFAULTS, **kw}
    r, eis, A_t = p['r'], p['eis'], p['A_target']
    e, pi, Pi = rouwenhorst(p['rho_e'], eps, p['n_e'])
    y_min = p['Y'] * e.min()
    a_min = -p['nbl_frac'] * y_min / r
    a_grid = np.linspace(a_min, p['a_max'], p['n_a'])
    if eps == 0:
        a_grid = asset_grid(a_min, p['a_max'], p['n_a'], include=A_t)
    blk = make_household([income_basic], name=kw.get('name', 'hh'))
    calib = dict(Y=p['Y'], e_grid=e, Pi=Pi, a_grid=a_grid, r=r, eis=eis)
    beta_ra = 1.0 / (1.0 + r)
    if eps == 0:
        j = int(np.argmin(np.abs(a_grid - A_t)))
        D = np.zeros((p['n_e'], len(a_grid)))
        D[:, j] = pi
        blk._Dbeg_override = D
        kappa = 0.0
    else:
        def resid(lk):
            ss = blk.steady_state({**calib, 'beta': beta_ra * (1 - np.exp(lk))}, **SS_OPTS)
            return ss['A'] - A_t
        lo, hi = np.log(1e-3 * eps ** 2), np.log(10 * eps ** 2)
        while resid(lo) < 0:
            lo -= 2.0
        while resid(hi) > 0:
            hi += 1.0
        lk = brentq(resid, lo, hi, xtol=1e-13, rtol=1e-13)
        kappa = float(np.exp(lk))
    calib['beta'] = beta_ra * (1 - kappa)
    ss = blk.steady_state(calib, **SS_OPTS)
    diag = dict(blk._ss_diag)
    diag['mass_constrained'] = mass_at_constraint(ss, a_grid, blk.name)
    return dict(block=blk, ss=ss, calib=calib, a_grid=a_grid, kappa=kappa, diag=diag,
                e_grid=e, Pi=Pi, pi=pi)


def household_jacobians(block, ss, T, inputs, outputs, h=1e-4):
    J = block.jacobian(ss, inputs, outputs, T=T, twosided=True, h=h)
    return {o: {i: np.asarray(J[o][i]) for i in inputs} for o in outputs}


def ssj_rank_limit_provider(eps, T, h=1e-4, **kw):
    """Interface `RankLimitProvider`: household Jacobians at idiosyncratic risk eps
    (two-sided finite differences with step h)."""
    m = pih_household(eps, **kw)
    ss = m['ss']
    J = household_jacobians(m['block'], ss, T, ['Y', 'r'], ['C'], h=h)
    return dict(J=J, ss=dict(C=float(ss['C']), A=float(ss['A']), r=float(ss['r']),
                             eis=float(ss['eis']), beta=float(ss['beta'])),
                ss_detail=dict(_household_detail(ss, m['block'].name), Y=float(ss['Y']),
                               r=float(ss['r']), eis=float(ss['eis'])),
                diag=dict(kappa=m['kappa'], **m['diag']))


# ---------------------------------------------------------------------------
# T2: TANK limit
# ---------------------------------------------------------------------------

TANK_DEFAULTS = dict(lam=0.3, beta_H=0.95, n_a_H=500, a_max_H=50.0)


@sj.simple
def tank_incidence(Y, Yss, lam, chi):
    Y_H = Yss * (Y / Yss) ** chi
    Y_S = (Y - lam * Y_H) / (1 - lam)
    return Y_H, Y_S


@sj.simple
def tank_aggregate(C_H, C_S, A_H, A_S, lam):
    C = lam * C_H + (1 - lam) * C_S
    A = lam * A_H + (1 - lam) * A_S
    return C, A


@sj.simple
def goods_clearing(Y, C):
    goods_mkt = Y - C
    return goods_mkt


def _remap_household(blk, suffix):
    names = ['Y', 'beta', 'e_grid', 'Pi', 'a_grid', 'C', 'A']
    return blk.remap({k: f'{k}_{suffix}' for k in names}).rename(f'hh_{suffix}')


def tank_model(eps, **kw):
    """Two permanent types sharing r and eis:
    H (share lam): borrowing limit 0, impatient beta_H, income risk eps -> hand-to-mouth;
    S (share 1-lam): near-natural limit, income risk eps, beta_S calibrated so that the
    bond is in zero net supply (lam A_H + (1-lam) A_S = 0).
    """
    p = {**PIH_DEFAULTS, **TANK_DEFAULTS, **kw}
    lam = p['lam']
    e, pi, Pi = rouwenhorst(p['rho_e'], eps, p['n_e'])
    hhH = make_household([income_basic], name='hh')
    a_grid_H = asset_grid(0.0, p['a_max_H'], p['n_a_H'])
    calib_H = dict(Y=1.0, e_grid=e, Pi=Pi, a_grid=a_grid_H, r=p['r'], eis=p['eis'],
                   beta=p['beta_H'])
    ssH = hhH.steady_state(calib_H, **SS_OPTS)
    A_H = float(ssH['A'])
    unconstrained_H = 1.0 - mass_at_constraint(ssH, a_grid_H, 'hh')
    A_S_target = -lam * A_H / (1 - lam)
    S = pih_household(eps, **{**{k: p[k] for k in PIH_DEFAULTS}, 'A_target': A_S_target})
    hh_H = _remap_household(hhH, 'H')
    hh_S = _remap_household(S['block'], 'S')
    model = sj.create_model([tank_incidence, hh_H, hh_S, tank_aggregate, goods_clearing],
                            name='TANK-limit')
    calib = dict(Y=1.0, Yss=1.0, lam=lam, chi=1.0, r=p['r'], eis=p['eis'],
                 beta_H=p['beta_H'], e_grid_H=e, Pi_H=Pi, a_grid_H=a_grid_H,
                 beta_S=S['calib']['beta'], e_grid_S=S['e_grid'], Pi_S=S['Pi'],
                 a_grid_S=S['a_grid'])
    ss = model.steady_state(calib, **SS_OPTS)
    return dict(model=model, ss=ss, hh_H=hh_H, hh_S=hh_S, lam=lam,
                diag=dict(A_H=A_H, unconstrained_H=unconstrained_H, kappa_S=S['kappa'],
                          mass_constrained_S=S['diag']['mass_constrained'],
                          goods_mkt_ss=float(ss['goods_mkt'])))


def ssj_tank_provider(eps, T, dr, chis=(1.0, 1.5), **kw):
    """Interface `TankProvider`: per-type iMPC matrices and GE output paths."""
    m = tank_model(eps, **kw)
    ss, model = m['ss'], m['model']
    J_H = m['hh_H'].jacobian(ss, ['Y_H', 'r'], ['C_H', 'A_H'], T=T, twosided=True)
    J_S = m['hh_S'].jacobian(ss, ['Y_S', 'r'], ['C_S', 'A_S'], T=T, twosided=True)
    Js = {'hh_H': J_H, 'hh_S': J_S}
    dY, ge_resid = {}, 0.0
    for chi in chis:
        ss_chi = copy.deepcopy(ss)
        ss_chi['chi'] = chi
        irf = model.solve_impulse_linear(ss_chi, ['Y'], ['goods_mkt'], {'r': dr}, Js=Js,
                                         options=_het_opts(['hh_H', 'hh_S']))
        dY[float(chi)] = np.asarray(irf['Y'])
        ge_resid = max(ge_resid, float(np.max(np.abs(irf['goods_mkt']))))
    return dict(M_H=np.asarray(J_H['C_H']['Y_H']), M_S=np.asarray(J_S['C_S']['Y_S']),
                Jr_S=np.asarray(J_S['C_S']['r']), dY=dY, lam=m['lam'],
                ss=dict(r=float(ss['r']), eis=float(ss['eis']), C_S=float(ss['C_S']),
                        A_S=float(ss['A_S']), Y=1.0),
                ss_detail=dict(lam=float(ss['lam']), r=float(ss['r']), eis=float(ss['eis']),
                               Y=float(ss['Y']), Y_H=float(ss['Y_H']), Y_S=float(ss['Y_S']),
                               **_household_detail(ss, 'hh_H', '_H'),
                               **_household_detail(ss, 'hh_S', '_S')),
                diag={**m['diag'], 'ge_target_resid': ge_resid})


# ---------------------------------------------------------------------------
# T3: Werning (2015) vanishing liquidity
# ---------------------------------------------------------------------------

WERNING_DEFAULTS = dict(eis=0.5, n_e=3, rho_e=0.9, sd_e=0.3, r_star=0.01, n_a=400,
                        grid_span=20.0)


@sj.simple
def huggett_clearing(Y, C, A):
    goods_mkt = C - Y
    asset_mkt = A
    return goods_mkt, asset_mkt


def werning_primitives(cyclical=False, **kw):
    """Income process, incidence zeta and beta such that the zero-liquidity rate is r_star.

    Income y_i = Y e_i (Y/Yss)^{zeta_i}. Acyclical benchmark (Werning eq. 15): zeta = 0.
    'Countercyclical risk' variant: zeta_low = 1 (lowest state's income moves 2-for-1 with
    Y), other states' zeta set so that sum_i pi_i e_i zeta_i = 0 (aggregate income = Y).
    """
    p = {**WERNING_DEFAULTS, **kw}
    e, pi, Pi = rouwenhorst(p['rho_e'], p['sd_e'], p['n_e'])
    gam = 1 / p['eis']
    ratio = (Pi @ e ** (-gam)) / e ** (-gam)
    beta = 1.0 / ((1 + p['r_star']) * ratio.max())
    zeta = np.zeros(p['n_e'])
    if cyclical:
        zeta[0] = 1.0
        zeta[1:] = -pi[0] * e[0] * zeta[0] / np.sum(pi[1:] * e[1:])
    return dict(e=e, pi=pi, Pi=Pi, beta=beta, zeta=zeta, eis=p['eis'], p=p)


def werning_model(b, cyclical=False, **kw):
    """Huggett economy: borrowing limit -b, bonds in zero net supply, r solves A(r) = 0."""
    prim = werning_primitives(cyclical, **kw)
    p = prim['p']
    blk = make_household([income_cyclical], name='hh')
    u = np.linspace(0, 1, p['n_a'])
    a_grid = -b + b * p['grid_span'] * u ** 2     # grid scaled with b
    calib = dict(Y=1.0, Yss=1.0, e_grid=prim['e'], zeta=prim['zeta'], Pi=prim['Pi'],
                 a_grid=a_grid, beta=prim['beta'], eis=p['eis'])
    opts = dict(backward_tol=1e-9 * b, backward_maxit=100_000)
    f = lambda r: blk.steady_state({**calib, 'r': r}, **opts)['A'] / b
    r_star = p['r_star']
    # A(r) < 0 for r < r* (nobody wants to save); r(b) - r* is O(b) (about 2.1 b here), so
    # [r* - b/2, r* + 5b] brackets the root while keeping trial savers inside the b-scaled grid
    lo, hi = r_star - 0.5 * b, min(r_star + 5 * b, 1 / prim['beta'] - 1 - 1e-6)
    r_ss = brentq(f, lo, hi, xtol=1e-15)
    model = sj.create_model([blk, huggett_clearing], name='Huggett-vanishing-liquidity')
    ss = model.steady_state({**calib, 'r': r_ss}, **opts)
    D = ss.internals['hh']['D']
    used = np.where(D.sum(0) > 1e-12)[0].max()
    diag = dict(r_ss=r_ss, r_star=r_star, r_gap_over_b=(r_ss - r_star) / b,
                mass_constrained=mass_at_constraint(ss, a_grid, 'hh'),
                max_wealth_over_b=float((a_grid[used] + b) / b),
                grid_top_over_b=float(p['grid_span']))
    return dict(model=model, block=blk, ss=ss, prim=prim, diag=diag)


def ssj_werning_provider(b, T, dr, cyclical=False, ge_target='asset_mkt', **kw):
    """Interface `WerningProvider`: GE output path for a real-rate path dr. The GE target is the
    asset market (A = 0, bonds in zero net supply) unless `ge_target='goods_mkt'`; the target
    passed to the solve is the one reported (`ge_target`, `diag['ge_target_resid']`)."""
    m = werning_model(b, cyclical, **kw)
    ss = m['ss']
    J = m['block'].jacobian(ss, ['Y', 'r'], ['C', 'A'], T=T, twosided=True, h=1e-3 * b)
    targets = [ge_target]
    irf = m['model'].solve_impulse_linear(ss, ['Y'], targets, {'r': dr}, Js={'hh': J},
                                          options=_het_opts(['hh'], h=1e-3 * b))
    prim = m['prim']
    return dict(dY=np.asarray(irf['Y']), r_ss=float(ss['r']), eis=prim['eis'],
                e_grid=prim['e'], Pi=prim['Pi'], zeta=prim['zeta'], beta=prim['beta'],
                ss_detail=dict(_household_detail(ss), b=float(b), Y=float(ss['Y']),
                               r=float(ss['r']), eis=float(ss['eis'])),
                J_A={'Y': np.asarray(J['A']['Y']), 'r': np.asarray(J['A']['r'])},
                J_C={'Y': np.asarray(J['C']['Y']), 'r': np.asarray(J['C']['r'])},
                ge_target=targets[0],
                diag={**m['diag'],
                      'ge_target_resid': float(np.max(np.abs(irf[targets[0]])))})


# ---------------------------------------------------------------------------
# T4: Aiyagari steady state
# ---------------------------------------------------------------------------

AIYAGARI_DEFAULTS = dict(beta=0.96, alpha=0.36, delta=0.08, eis=0.5, n_e=7, n_a=1000,
                         a_max=200.0)


def ssj_aiyagari_r(sd, rho, chain=None, **kw):
    """Interface `AiyagariProvider`: steady-state r for log-income sd and persistence.

    The root is searched on (-delta + 1e-3, 1/beta - 1 - 1e-10), so `r < 1/beta - 1` holds
    by construction and is NOT an independent accuracy check (README section 7). The
    non-mechanical diagnostics returned in `diag` are: the capital-market residual at the
    root, unit-free Euler errors on unconstrained states and KKT violations on constrained
    states (`hank_limits.diagnostics.euler_kkt_residuals`), the mass in the top grid points
    and the mass whose policy leaves the grid.

    `chain=(e, Pi)` replaces the Rouwenhorst chain (used for the Aiyagari Table II cells,
    scripts/aiyagari_table2_v0_2.py); `eis` = 1/risk aversion.
    """
    from .diagnostics import euler_kkt_residuals
    p = {**AIYAGARI_DEFAULTS, **kw}
    al, de, be = p['alpha'], p['delta'], p['beta']
    if chain is None:
        e, pi, Pi = rouwenhorst(rho, sd, p['n_e'])
    else:
        e, Pi = chain
    u = np.linspace(0, 1, p['n_a'])
    a_grid = p['a_max'] * u ** 2
    blk = make_household([income_wage], name='hh')

    def resid(r):
        K = (al / (r + de)) ** (1 / (1 - al))
        w = (1 - al) * K ** al
        ss = blk.steady_state(dict(w=w, e_grid=e, Pi=Pi, a_grid=a_grid, r=r, beta=be,
                                   eis=p['eis']), **SS_OPTS)
        return ss['A'] - K, ss

    rho_disc = 1 / be - 1
    lo, hi = -de + 1e-3, rho_disc - 1e-10
    r = brentq(lambda r: resid(r)[0], lo, hi, xtol=1e-14)
    gap, ss = resid(r)
    K = (al / (r + de)) ** (1 / (1 - al))
    st = {**ss.toplevel, **ss.internals['hh']}
    ee = euler_kkt_residuals(st['Pi'], a_grid, st['a'], st['c'], st['D'], r, be, p['eis'])
    diag = dict(mass_constrained=mass_at_constraint(ss, a_grid, 'hh'),
                mass_top_5pts=blk._ss_diag['mass_top_5pts'],
                mass_policy_off_grid=blk._ss_diag['mass_policy_off_grid'],
                ssj_forward_resid=blk._ss_diag['ssj_forward_resid'],
                capital_market_resid_rel=float(abs(gap) / K),
                root_distance_to_upper_bracket=float(hi - r),
                euler=ee, K=float(K), saving_rate=float(de * al / (r + de)),
                a_max=float(p['a_max']), n_a=int(p['n_a']))
    return dict(r=float(r), rho_disc=rho_disc, diag=diag, _model=dict(block=blk, ss=ss),
                ss_detail=dict(_household_detail(ss), r=float(ss['r']), w=float(ss['w']),
                               K=float(K), eis=float(ss['eis'])),
                policy=dict(a_grid=a_grid, a=st['a'], c=st['c'], D=st['D']))


# ---------------------------------------------------------------------------
# T5 / T6: fiscal HANK with a fixed (or rule-based) real rate
# ---------------------------------------------------------------------------

FISCAL_DEFAULTS = dict(r=0.02, eis=0.5, B=1.0, G=0.2, phi_B=0.1, n_e=7, rho_e=0.92,
                       sd_e=0.6, n_a=500, a_max=60.0, phi_y=0.25)


@sj.simple
def fiscal(B, G, r, Tax_ss, B_ss, phi_B):
    Tax = Tax_ss + phi_B * (B(-1) - B_ss)
    fiscal_res = B - (1 + r) * B(-1) - G + Tax
    return Tax, fiscal_res


@sj.simple
def fiscal_clearing(Y, C, G, A, B):
    goods_mkt = Y - C - G
    asset_mkt = A - B
    return goods_mkt, asset_mkt


@sj.simple
def real_rate_rule(Y, rss, phi_y, rshock):
    r = rss + phi_y * (Y - 1) + rshock
    return r


def fiscal_hank(**kw):
    """One-asset HANK (limit 0), income (Y - Tax) e, government debt B = A, beta calibrated."""
    p = {**FISCAL_DEFAULTS, **kw}
    e, pi, Pi = rouwenhorst(p['rho_e'], p['sd_e'], p['n_e'])
    a_grid = asset_grid(0.0, p['a_max'], p['n_a'])
    hh = make_household([income_aftertax], name='hh')
    Tax_ss = p['G'] + p['r'] * p['B']
    calib = dict(Y=1.0, Tax=Tax_ss, e_grid=e, Pi=Pi, a_grid=a_grid, r=p['r'], eis=p['eis'])
    f = lambda beta: hh.steady_state({**calib, 'beta': beta}, **SS_OPTS)['A'] - p['B']
    beta = brentq(f, 0.80, 1 / (1 + p['r']) - 1e-6, xtol=1e-14)
    calib.update(beta=beta, B=p['B'], B_ss=p['B'], G=p['G'], Tax_ss=Tax_ss,
                 phi_B=p['phi_B'], rss=p['r'], phi_y=p['phi_y'], rshock=0.0)
    model_fixed = sj.create_model([hh, fiscal, fiscal_clearing], name='fiscal-HANK-fixed-r')
    model_rule = sj.create_model([real_rate_rule, hh, fiscal, fiscal_clearing],
                                 name='fiscal-HANK-r-rule')
    ss = model_fixed.steady_state(calib, **SS_OPTS)
    ss_rule = model_rule.steady_state(calib, **SS_OPTS)
    diag = dict(beta=beta, mass_constrained=mass_at_constraint(ss, a_grid, 'hh'),
                goods_mkt_ss=float(ss['goods_mkt']), asset_mkt_ss=float(ss['asset_mkt']),
                fiscal_res_ss=float(ss['fiscal_res']))
    return dict(hh=hh, ss=ss, ss_rule=ss_rule, model_fixed=model_fixed,
                model_rule=model_rule, a_grid=a_grid, calib=calib, p=p, diag=diag)


def ssj_fiscal_provider(T, shocks, ssj_consistent=True, solver='ssj', **kw):
    """Interface `IKCProvider`: household Jacobians + SSJ GE paths for each shock.

    `shocks` maps a label to (model_kind, {input: path}) with model_kind in
    {'fixed_r', 'r_rule'}; exogenous inputs are G and r (fixed_r) or G and rshock (r_rule).
    """
    m = fiscal_hank(**kw)
    ss = m['ss']
    J = household_jacobians(m['hh'], ss, T, ['Y', 'Tax', 'r'], ['C', 'A'])
    Jd = {'hh': sj.JacobianDict({o: {i: J[o][i] for i in J[o]} for o in J}, name='hh', T=T)}
    out = {}
    for label, (kind, inputs) in shocks.items():
        model = m['model_fixed'] if kind == 'fixed_r' else m['model_rule']
        ssk = ss if kind == 'fixed_r' else m['ss_rule']
        opts = _het_opts(['hh'], ssj_consistent)
        if solver == 'ssj':
            irf = model.solve_impulse_linear(ssk, ['Y', 'B'], ['goods_mkt', 'fiscal_res'],
                                             inputs, Js=Jd, options=opts)
        else:
            irf = solve_impulse_linear_refined(model, ssk, ['Y', 'B'],
                                               ['goods_mkt', 'fiscal_res'], inputs, Jd, opts)
        T_ = len(next(iter(inputs.values())))

        def get(k):
            if k in inputs:                      # exogenous paths: take them as given
                return np.asarray(inputs[k], float)
            return np.asarray(irf.toplevel[k]) if k in irf.toplevel else np.zeros(T_)
        out[label] = dict(dY=get('Y'), dC=get('C'), dG=get('G'), dTax=get('Tax'),
                          dr=get('r'), dB=get('B'), dA=get('A'),
                          drshock=get('rshock'),
                          asset_mkt=get('asset_mkt'), goods_mkt=get('goods_mkt'))
    return dict(J=J, ss=dict(C=float(ss['C']), A=float(ss['A']), r=float(ss['r']),
                             eis=float(ss['eis']), Y=float(ss['Y']), G=float(ss['G'])),
                ss_detail=dict(_household_detail(ss), Y=float(ss['Y']), B=float(ss['B']),
                               G=float(ss['G']), Tax=float(ss['Tax']), r=float(ss['r']),
                               eis=float(ss['eis']), phi_B=float(ss['phi_B'])),
                irfs=out, model=m)


def cohort_impc_paths(m, x, groups, h=1e-4):
    """Consumption response of cohorts that start in wealth group k to a disposable-income
    path x (fixed r, fixed taxes), by two-sided nonlinear simulation in SSJ.

    groups: list of boolean masks over the asset grid (beginning-of-period wealth).
    Returns (mu_k, list of dC_k paths) with sum_k mu_k dC_k = M x.
    """
    hh, ss = m['hh'], m['ss']
    Dbeg = ss.internals['hh']['Dbeg']
    mus, paths = [], []
    for mask in groups:
        Dk = Dbeg * mask[None, :]
        mu = float(Dk.sum())
        ss_k = copy.deepcopy(ss)
        ss_k.internals['hh']['Dbeg'] = Dk / mu
        up = hh.impulse_nonlinear(ss, {'Y': h * x}, outputs=['C'], ss_initial=ss_k)['C']
        dn = hh.impulse_nonlinear(ss, {'Y': -h * x}, outputs=['C'], ss_initial=ss_k)['C']
        mus.append(mu)
        paths.append((np.asarray(up) - np.asarray(dn)) / (2 * h))
    return mus, paths
