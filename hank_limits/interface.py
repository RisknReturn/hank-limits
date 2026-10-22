"""Interface for plugging another HANK code into the analytic-limit checks.

A code under test supplies *providers*: plain Python callables that return numpy arrays.
Nothing in `hank_limits.checks` depends on SSJ. The SSJ reference providers are in
`hank_limits.models` (functions `ssj_*`). Shapes: T x T for Jacobians (rows = output date
t, columns = shock date s), length-T vectors for impulse responses. Units: levels
(dC, dY in units of steady-state output; dr in absolute rate units, e.g. 0.01 = 1pp).

Timing convention (must match, or convert): r_t is the ex-post real return paid at
date t on assets held from t-1 (SSJ convention), so the household Euler equation between
t and t+1 involves r_{t+1}, and a date-0 change in r_0 only has a wealth effect.

Fields. Every field that `hank_limits.checks` reads is listed below, and no other is read. A
REQUIRED field is used by the v0.1 part of a check (T1-T6); without it the check fails. An
OPTIONAL field is used only by the v0.3 checks V1-V4 (docs/CHECKS_PLAN_v0_3.md; revised for
v0.3.0, docs/CHECKS_v0_3_0_AMENDMENT.md) or by a diagnostic. If a provider does not report an
optional field, the V quantities that need it are skipped and the reason is recorded under
res['v0_3']['skipped']; each quantity is computed from the fields its own equation needs, so the
others are still computed. y = income by income state at the steady state, Pi = the income
transition matrix (rows = current state), a_min = borrowing limit.

RankLimitProvider(eps, T) -> dict(
    J  = {'C': {'Y': M (T x T), 'r': Jr (T x T)}},   # REQUIRED: household Jacobians at risk eps
    ss = {'C', 'A', 'r', 'eis'},                    # REQUIRED: steady-state aggregates (floats)
    diag = {'kappa', 'mass_constrained', 'ssj_forward_resid'},   # optional diagnostics
    ss_detail = {...})                              # optional (V checks, below)

TankProvider(eps, T, dr, chis) -> dict(
    M_H, M_S (T x T)  : REQUIRED, iMPC matrices of the hand-to-mouth and saver types,
    dY = {chi: path}  : REQUIRED, GE output response to real-rate path dr with goods-market
                        clearing and H-income elasticity chi,
    lam, ss = {'r', 'eis', 'C_S', 'A_S', 'Y'} : REQUIRED,
    diag = {'unconstrained_H', 'A_H', 'kappa_S', 'mass_constrained_S', 'ge_target_resid'}: optional,
    ss_detail = {...})                          : optional (V checks)

WerningProvider(b, T, dr, cyclical) -> dict(
    dY        : REQUIRED, GE output response to dr (zero net supply, borrowing limit -b),
    r_ss, eis : REQUIRED, steady-state rate and EIS,
    e_grid, Pi, zeta : REQUIRED, the idiosyncratic income chain and incidence
                (y_i = Y e_i (Y/Yss)^zeta_i),
    diag = {'r_gap_over_b', 'mass_constrained', 'max_wealth_over_b'}: REQUIRED;
           diag['ge_target_resid']: optional, max |residual of the GE target the solve imposed|,
    ge_target : optional, 'asset_mkt' or 'goods_mkt': the market whose clearing the solve
                imposed (v0.3.0). Absent: 'asset_mkt' is assumed and recorded as assumed. Any
                other value: V1 does not apply ("applicability condition not met", skipped).
    J_A, J_C, ss_detail : optional (V checks))

AiyagariProvider(sd, rho) -> dict(
    r, rho_disc = 1/beta - 1 : REQUIRED,
    diag = {'mass_constrained'}: REQUIRED; diag 'capital_market_resid_rel', 'mass_top_5pts',
           'mass_policy_off_grid', 'euler' = {'euler_err_mean_weighted',
           'euler_err_max_with_mass', 'kkt_violation_max_with_mass'}: optional,
    ss_detail, policy : optional (V checks))

IKCProvider(T, shocks) -> dict(
    J = {'C': {'Y', 'Tax', 'r'}},   # REQUIRED: household Jacobians (income is (Y - Tax) e)
    ss = {'C', 'A', 'r', 'eis'},    # REQUIRED for T6 (`check_differencing`)
    irfs = {label: {'dY', 'dG', 'dTax', 'dr', 'asset_mkt'}},   # REQUIRED: GE paths from the
                                    # code's solver, one entry per experiment the caller requested
    ss_detail, irfs[label]['dC', 'dB', 'drshock'] : optional (V checks))
`shocks` maps a label to (kind, {input: path}), kind 'fixed_r' (exogenous G and/or r) or
'r_rule' (exogenous G and/or rshock; r_t = r + phi_y (Y_t - 1) + rshock_t). The checks take the
requested paths from `shocks` and compare the returned dG, dr, drshock with them (v0.3.0, R3);
an input that an experiment does not request is a zero path.

The T6 decomposition (`checks.check_differencing`) needs only the IKC provider output.

Optional fields of the V checks, and the V quantities each one feeds (v0.3.0). Fields in
ss_detail are floats unless noted (e_grid, Pi, y: arrays of the income chain).

RankLimitProvider  ss_detail: Y, C, A, r, eis, a_min, e_grid, Pi, y
    V2 hh_budget (C, A, r, y, Pi); V3 Pi checks (Pi), mean_e_err (e_grid, Pi), agg_income_err
    (y, Pi, Y); V4(a) spec_r (r), spec_eis (eis), spec_A (A), spec_a_min (a_min, e_grid).
TankProvider       ss_detail: lam, r, eis, Y, Y_H, Y_S and, for t in (H, S), C_t, A_t, a_min_t,
                   e_grid_t, Pi_t, y_t, and beta_H
    V2 ss_goods (Y, lam, C_H, C_S), ss_assets (lam, A_H, A_S, Y), hh_budget_t (C_t, A_t, r, y_t,
    Pi_t); V3 per type (e_grid_t, Pi_t, y_t, Y_t), incidence_err (lam, Y_H, Y_S, Y); V4(a) lam,
    beta_H, r, eis, a_min_H, a_min_S (with e_grid_S).
WerningProvider    ss_detail: Y, C, A, r, a_min, y, beta;  J_A = {'Y': J^{A,Y}, 'r': J^{A,r}},
                   J_C = {'Y': J^{C,Y}, 'r': J^{C,r}} (T x T, the Jacobians used in the solve)
    V1 (by `ge_target`; residual / b, rate path dr as passed by the check): reported =
    diag['ge_target_resid']; direct = max |J^{A,Y} dY + J^{A,r} dr| (asset market, J_A) or
    max |dY - J^{C,Y} dY - J^{C,r} dr| (goods market, J_C). The other market's residual, if its
    Jacobians are reported, is recorded without threshold under res['v0_3']['diagnostics'].
    V1 applies to a T3 solve on one of these two markets.
    V2 ss_assets (A), ss_goods (C, Y), hh_budget (C, A, r, y, Pi); V3 (e_grid, Pi; y, Y),
    incidence_err_cyc (e_grid, Pi, zeta); V4(a) spec_eis (eis), spec_a_min (a_min),
    spec_zeta (zeta), spec_rstar (ss_detail beta; eis, e_grid, Pi).
AiyagariProvider   ss_detail: r, w, K, A, C, beta, eis, a_min, e_grid, Pi, y;
                   policy = {'a_grid', 'a', 'c', 'D'} (policies and distribution, n_e x n_a)
    V2 ss_capital (A, K), hh_budget (C, A, r, y, Pi), foc_r (r, K, e_grid, Pi), foc_w (w, K,
    e_grid, Pi), ss_goods (K, C, e_grid, Pi); V3 (e_grid, Pi; y, w); V4(a) beta, eis, a_min;
    V4(b) Euler (policy, Pi, r).
IKCProvider        ss_detail: Y, C, A, B, G, Tax, r, eis, phi_B, a_min, e_grid, Pi, y;
                   irfs[label]: dC, dB, drshock (exogenous rate-shock path; zeros if none)
    V2 ss_government_budget (Tax, G, r, B, Y), ss_goods (Y, C, G), ss_assets (A, B), hh_budget
    (C, A, r, y, Pi); V3 (e_grid, Pi; y, Y, Tax); V4(a) r, eis, B, G, phi_B, a_min.
    Paths (each experiment the caller requested, each T): fiscal_rule (dB, dTax),
    government_budget (dB, dTax, dr, dG), goods_market (dY, dC, dG), rate_rule (dr, dY and the
    requested rshock; rule-based experiments), exog_G (dG), exog_r (dr; fixed-rate
    experiments), exog_rshock (drshock).

The SSJ reference providers report more than this (e.g. beta in several ss_detail dicts, and a
top-level beta in the Werning provider, kept for existing callers); the checks do not read it.

The documented values the checks compare with are in `hank_limits.spec` (passed as `spec=`);
another code with another calibration passes its own.

Minimal example for a user code `mycode` (pseudo-code):

    from hank_limits.checks import check_rank_limit
    def my_provider(eps, T):
        model = mycode.build(income_sd=eps, borrowing='natural', target_assets=5.0)
        Jc = mycode.household_jacobian(model, inputs=['Y', 'r'], T=T)
        return dict(J={'C': {'Y': Jc['Y'], 'r': Jc['r']}},
                    ss=dict(C=model.C, A=model.A, r=model.r, eis=model.eis))
    res = check_rank_limit(my_provider, eps_list=[0.2, 0.1, 0.05, 0.025], T=200)
    assert res['fitted_rate_err_M'] > 1.5
"""
from typing import Callable, Dict, Any

RankLimitProvider = Callable[..., Dict[str, Any]]
TankProvider = Callable[..., Dict[str, Any]]
WerningProvider = Callable[..., Dict[str, Any]]
AiyagariProvider = Callable[..., Dict[str, Any]]
IKCProvider = Callable[..., Dict[str, Any]]
