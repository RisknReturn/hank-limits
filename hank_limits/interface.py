"""Interface for plugging another HANK code into the analytic-limit checks.

A code under test supplies *providers*: plain Python callables that return numpy arrays.
Nothing in `hank_limits.checks` depends on SSJ. The SSJ reference providers are in
`hank_limits.models` (functions `ssj_*`). Shapes: T x T for Jacobians (rows = output date
t, columns = shock date s), length-T vectors for impulse responses. Units: levels
(dC, dY in units of steady-state output; dr in absolute rate units, e.g. 0.01 = 1pp).

Timing convention (must match, or convert): r_t is the ex-post real return paid at
date t on assets held from t-1 (SSJ convention), so the household Euler equation between
t and t+1 involves r_{t+1}, and a date-0 change in r_0 only has a wealth effect.

RankLimitProvider(eps, T) -> dict(
    J  = {'C': {'Y': M (T x T), 'r': Jr (T x T)}},   # household Jacobians at risk eps
    ss = {'C', 'A', 'r', 'eis'},                    # steady-state aggregates (floats)
    diag = {...optional diagnostics: 'kappa', 'mass_constrained'})

TankProvider(eps, T, dr, chis) -> dict(
    M_H, M_S (T x T)  : iMPC matrices of the hand-to-mouth and saver types,
    dY = {chi: path}  : GE output response to real-rate path dr with goods-market
                        clearing and H-income elasticity chi,
    lam, ss = {'r', 'eis', 'C_S', 'A_S', 'Y'}, diag)

WerningProvider(b, T, dr, cyclical) -> dict(
    dY        : GE output response to dr (zero net supply, borrowing limit -b),
    r_ss, eis : steady-state rate and EIS,
    e_grid, Pi, zeta : the idiosyncratic income chain and incidence (y_i = Y e_i (Y/Yss)^zeta_i),
    diag = {'r_gap_over_b', 'mass_constrained', 'max_wealth_over_b'})

AiyagariProvider(sd, rho) -> dict(r, rho_disc = 1/beta - 1, diag={'mass_constrained'})

IKCProvider(T, shocks) -> dict(
    J = {'C': {'Y', 'Tax', 'r'}},   # household Jacobians (income is (Y - Tax) e)
    ss = {'C', 'A', 'r', 'eis', 'Y', 'G'},
    irfs = {label: {'dY','dG','dTax','dr','asset_mkt'}})   # GE paths from the code's solver

The T6 decomposition (`checks.check_differencing`) needs only the IKC provider output.

v0.3: optional fields for the steady-state, normalisation and specification checks (V1-V4,
docs/CHECKS_PLAN_v0_3.md). The checks compute their residuals from these primitive quantities
themselves; a provider that does not report a field gets those checks skipped (the reason is
recorded in the result under 'v0_3' -> 'skipped'). y = income by income state at the steady
state, Pi = the income transition matrix (rows = current state), a_min = borrowing limit.

RankLimitProvider: ss_detail = {'Y', 'C', 'A', 'r', 'eis', 'beta', 'a_min', 'e_grid', 'Pi', 'y'}
TankProvider:      ss_detail = {'lam', 'r', 'eis', 'Y', 'Y_H', 'Y_S',
                                'C_t', 'A_t', 'beta_t', 'a_min_t', 'e_grid_t', 'Pi_t', 'y_t'
                                for t in (H, S)}
WerningProvider:   ss_detail = {'b', 'Y', 'C', 'A', 'r', 'eis', 'beta', 'a_min', 'e_grid', 'Pi', 'y'},
                   J_A = {'Y': J^{A,Y}, 'r': J^{A,r}}   (T x T, the Jacobians used in the solve)
AiyagariProvider:  ss_detail = {'r', 'w', 'K', 'A', 'C', 'beta', 'eis', 'a_min', 'e_grid', 'Pi', 'y'},
                   policy = {'a_grid', 'a', 'c', 'D'}   (policies and distribution, n_e x n_a)
IKCProvider:       ss_detail = {'Y', 'C', 'A', 'B', 'G', 'Tax', 'r', 'eis', 'beta', 'phi_B',
                                'a_min', 'e_grid', 'Pi', 'y'},
                   irfs[label] also 'dC', 'dB', 'drshock' (exogenous rate-shock path, zeros if none)

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
