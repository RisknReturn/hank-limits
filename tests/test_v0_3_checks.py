"""v0.3 checks V1-V4 and V6 (docs/CHECKS_PLAN_v0_3.md, section 3; V5 is in test_t6_differencing.py).

The quantities are computed in `hank_limits.checks` (numpy only) from primitive quantities the
providers report in the optional fields `ss_detail`, `J_A`, `policy` and the T5 paths, against
the documented values in `hank_limits.spec`, and stored under the key 'v0_3' of each check
result. A provider that does not report a field gets a skip, with the reason.

Thresholds (plan, section 4): threshold(q) = roundup_1sig(max(100 * b_max(q), 10 * tau(q))),
b_max = largest value over all settings in two runs of the unmodified suite before commit T,
tau = tolerance of the numerical step that determines q. These are regression thresholds for
this calibration, like the tolerances of T1-T6.
"""
import numpy as np
import pytest

from sequence_jacobian.hetblocks.hh_sim import hh as _hh_sim, make_grids as _make_grids
from sequence_jacobian.utilities import function as _ssj_function

# fixture -> quantity -> threshold (plan, section 4, table; values from two runs of the unmodified
# suite at commit T). 't5' and 't5_paths' are the T5 setting with matched options and refinement.
THRESH = {
    't1': {
        'hh_budget': 4e-10, 'pi_row_sum_err': 1e-13, 'pi_negative': 1e-13, 'mean_e_err': 2e-10,
        'agg_income_err': 2e-10, 'spec_r': 1e-13, 'spec_eis': 1e-13, 'spec_A': 6e-06,
        'spec_a_min': 1e-13,
    },
    't2': {
        'ss_goods': 7e-07, 'ss_assets': 3e-05, 'hh_budget_H': 1e-13, 'hh_budget_S': 1e-13,
        'pi_row_sum_err_H': 1e-13, 'pi_negative_H': 1e-13, 'mean_e_err_H': 2e-10,
        'agg_income_err_H': 2e-10, 'pi_row_sum_err_S': 1e-13, 'pi_negative_S': 1e-13,
        'mean_e_err_S': 2e-10, 'agg_income_err_S': 2e-10, 'incidence_err': 1e-13,
        'spec_lam': 1e-13, 'spec_beta_H': 1e-13, 'spec_r': 1e-13, 'spec_eis': 1e-13,
        'spec_a_min_H': 1e-13, 'spec_a_min_S': 1e-13,
    },
    't3': {
        'ge_target_resid_reported_acyc': 2e-12, 'ge_target_resid_reported_cyc': 2e-12,
        'ge_target_resid_direct_acyc': 2e-12, 'ge_target_resid_direct_cyc': 2e-12,
        'ss_assets_acyc': 2e-08, 'ss_assets_cyc': 2e-08, 'ss_goods_acyc': 4e-10,
        'ss_goods_cyc': 4e-10, 'hh_budget_acyc': 1e-13, 'hh_budget_cyc': 1e-13,
        'pi_row_sum_err_acyc': 1e-13, 'pi_negative_acyc': 1e-13, 'mean_e_err_acyc': 4e-10,
        'agg_income_err_acyc': 4e-10, 'pi_row_sum_err_cyc': 1e-13, 'pi_negative_cyc': 1e-13,
        'mean_e_err_cyc': 4e-10, 'agg_income_err_cyc': 4e-10, 'incidence_err_cyc': 2e-09,
        'spec_eis_acyc': 1e-13, 'spec_eis_cyc': 1e-13, 'spec_a_min_acyc': 1e-13,
        'spec_a_min_cyc': 1e-13, 'spec_zeta_acyc': 1e-13, 'spec_zeta_low_cyc': 1e-13,
        'spec_rstar_acyc': 1e-13, 'spec_rstar_cyc': 1e-13,
    },
    't4': {
        'ss_capital': 6e-08, 'foc_r': 2e-10, 'foc_w': 7e-10, 'ss_goods': 7e-09,
        'hh_budget': 2e-13, 'pi_row_sum_err': 1e-13, 'pi_negative': 1e-13, 'mean_e_err': 2e-09,
        'agg_income_err': 2e-09, 'spec_beta': 1e-13, 'spec_eis': 1e-13, 'spec_a_min': 1e-13,
        'spec_euler_err_mean_weighted': 5e-06, 'spec_euler_err_max_with_mass': 3e-04,
        'spec_kkt_violation_max_with_mass': 1e-11,
    },
    't5': {
        'ss_government_budget': 1e-13, 'ss_goods': 4e-09, 'ss_assets': 2e-12, 'hh_budget': 6e-13,
        'pi_row_sum_err': 1e-13, 'pi_negative': 1e-13, 'mean_e_err': 5e-09,
        'agg_income_err': 5e-09, 'spec_r': 1e-13, 'spec_eis': 1e-13, 'spec_B': 1e-13,
        'spec_G': 1e-13, 'spec_phi_B': 1e-13, 'spec_a_min': 1e-13,
    },
    't5_paths': {
        'fiscal_rule': 1e-13, 'government_budget': 1e-13, 'rate_rule': 1e-13,
        'goods_market': 4e-11,
    },
}

FISCAL = 'ssj_consistent_refined'     # T5 setting used by the V checks (as the Walras test)


def _v03(res):
    if 'v0_3' not in res:
        pytest.skip('check result has no v0.3 quantities')
    if res['v0_3'].get('errors'):
        pytest.fail(f"v0.3 computation failed: {res['v0_3']['errors'][:3]}")
    return res['v0_3']


def _assert_below(v, keys, fx):
    """Every quantity in `keys` (lists over the check's settings) below its threshold."""
    missing = [k for k in keys if not [x for x in v.get(k, []) if x is not None]]
    if missing:
        pytest.skip(f"not reported: {', '.join(missing)}; skipped: {v.get('skipped')}")
    bad = {}
    for k in keys:
        worst = max(x for x in v[k] if x is not None)
        if not worst < THRESH[fx][k]:
            bad[k] = (worst, THRESH[fx][k])
    assert not bad, bad


def _fixture_v03(request, name):
    res = request.getfixturevalue(name)
    return _v03(res[FISCAL] if name == 't5' else res)


# ---------------------------------------------------------------------------
# V1: T3 GE target residual (asset market, zero net supply), divided by b
# ---------------------------------------------------------------------------

def test_v1_werning_ge_target_resid_reported(t3):
    # the solver residual the T3 provider reports, recorded as in T2
    _assert_below(_v03(t3), ['ge_target_resid_reported_acyc', 'ge_target_resid_reported_cyc'], 't3')


def test_v1_werning_ge_target_resid_direct(t3):
    # max |J^{A,Y} dY + J^{A,r} dr| / b from the reported Jacobians and the rate path passed
    _assert_below(_v03(t3), ['ge_target_resid_direct_acyc', 'ge_target_resid_direct_cyc'], 't3')


# ---------------------------------------------------------------------------
# V2: steady-state residuals
# ---------------------------------------------------------------------------

def test_v2_fiscal_ss_government_budget(t5):
    _assert_below(_v03(t5[FISCAL]), ['ss_government_budget'], 't5')


def test_v2_fiscal_ss_goods_market(t5):
    _assert_below(_v03(t5[FISCAL]), ['ss_goods'], 't5')


def test_v2_fiscal_ss_asset_market(t5):
    _assert_below(_v03(t5[FISCAL]), ['ss_assets'], 't5')


def test_v2_tank_ss_goods_market(t2):
    _assert_below(_v03(t2), ['ss_goods'], 't2')


def test_v2_tank_ss_asset_market(t2):
    _assert_below(_v03(t2), ['ss_assets'], 't2')


def test_v2_werning_ss_asset_market(t3):
    _assert_below(_v03(t3), ['ss_assets_acyc', 'ss_assets_cyc'], 't3')


def test_v2_werning_ss_goods_market(t3):
    _assert_below(_v03(t3), ['ss_goods_acyc', 'ss_goods_cyc'], 't3')


def test_v2_aiyagari_ss_capital_market(t4):
    _assert_below(_v03(t4), ['ss_capital'], 't4')


def test_v2_aiyagari_ss_firm_focs(t4):
    # documented alpha, delta; L = mean of the reported income states
    _assert_below(_v03(t4), ['foc_r', 'foc_w'], 't4')


def test_v2_aiyagari_ss_goods_market(t4):
    _assert_below(_v03(t4), ['ss_goods'], 't4')


_HH_BUDGET = {'t1': ['hh_budget'], 't2': ['hh_budget_H', 'hh_budget_S'],
              't3': ['hh_budget_acyc', 'hh_budget_cyc'], 't4': ['hh_budget'], 't5': ['hh_budget']}


@pytest.mark.parametrize('name', sorted(_HH_BUDGET))
def test_v2_household_budget(name, request):
    # |C - r A - sum_e pi_e y_e| / C, pi the stationary distribution of the reported Pi
    _assert_below(_fixture_v03(request, name), _HH_BUDGET[name], name)


# ---------------------------------------------------------------------------
# V3: income normalisation
# ---------------------------------------------------------------------------

_INC = ('pi_row_sum_err', 'pi_negative', 'mean_e_err', 'agg_income_err')
_INCOME = {'t1': list(_INC), 't2': [f'{k}_{t}' for t in 'HS' for k in _INC] + ['incidence_err'],
           't3': [f'{k}_{t}' for t in ('acyc', 'cyc') for k in _INC] + ['incidence_err_cyc'],
           't4': list(_INC), 't5': list(_INC)}


@pytest.mark.parametrize('name', sorted(_INCOME))
def test_v3_income_normalisation(name, request):
    # Pi row-stochastic, E_pi[e] = 1, aggregate income = Y (T2 per type and incidence; T4
    # w x documented labour 1; T5 Y - Tax); T3 countercyclical: sum pi e zeta = 0
    _assert_below(_fixture_v03(request, name), _INCOME[name], name)


# ---------------------------------------------------------------------------
# V4: specification pins
# ---------------------------------------------------------------------------

_SPEC = {'t1': ['spec_r', 'spec_eis', 'spec_A', 'spec_a_min'],
         't2': ['spec_lam', 'spec_beta_H', 'spec_r', 'spec_eis', 'spec_a_min_H', 'spec_a_min_S'],
         't3': ['spec_eis_acyc', 'spec_eis_cyc', 'spec_a_min_acyc', 'spec_a_min_cyc',
                'spec_zeta_acyc', 'spec_zeta_low_cyc', 'spec_rstar_acyc', 'spec_rstar_cyc'],
         't4': ['spec_beta', 'spec_eis', 'spec_a_min'],
         't5': ['spec_r', 'spec_eis', 'spec_B', 'spec_G', 'spec_phi_B', 'spec_a_min']}


@pytest.mark.parametrize('name', sorted(_SPEC))
def test_v4_spec_parameters(name, request):
    # reported parameters against hank_limits.spec (T3: r* from the reported beta, Pi, e)
    _assert_below(_fixture_v03(request, name), _SPEC[name], name)


def test_v4_aiyagari_euler_at_spec_beta(t4):
    # Euler errors of the reported policy at the documented beta and EIS
    _assert_below(_v03(t4), ['spec_euler_err_mean_weighted', 'spec_euler_err_max_with_mass',
                             'spec_kkt_violation_max_with_mass'], 't4')


def _paths(t5, key, labels=None):
    paths = _v03(t5[FISCAL]).get('paths', {})
    vals = {lab: p[key] for lab, p in paths.items() if key in p and (labels is None or lab in labels)}
    if not vals:
        pytest.skip(f'no path quantity {key}')
    return {f'{key}[{lab}]': x for lab, x in vals.items()}, key


def _assert_paths(t5, key, labels=None):
    vals, key = _paths(t5, key, labels)
    th = THRESH['t5_paths'][key]
    bad = {k: (max(x), th) for k, x in vals.items() if not max(x) < th}
    assert not bad, bad


def test_v4_fiscal_rule_on_path(t5):
    # dTax_t = phi_B dB_{t-1}, documented phi_B
    _assert_paths(t5, 'fiscal_rule')


def test_v4_government_budget_on_path(t5):
    # dB_t = (1 + r) dB_{t-1} + B dr_t + dG_t - dTax_t, documented r and B
    _assert_paths(t5, 'government_budget')


def test_v4_real_rate_rule_on_path(t5):
    # dr_t = phi_y dY_t + drshock_t (rule-based shocks), documented phi_y
    _assert_paths(t5, 'rate_rule')


def test_v4_goods_market_on_path(t5):
    # dY_t = dC_t + dG_t
    _assert_paths(t5, 'goods_market')


# ---------------------------------------------------------------------------
# V6: finding F3 (SSJ code only)
# ---------------------------------------------------------------------------

def _income_nl(w, e_grid):
    zeta = np.linspace(0.6, -0.4, e_grid.size)
    y = w ** (1 + zeta) * e_grid
    return y


def _income_lin(w, e_grid):
    y = w * e_grid
    return y


def _c2(c):
    c2 = c ** 2
    return c2


# strict xfail = regression pin for sequence-jacobian 1.0.0 only (as F1, F2)
@pytest.mark.xfail(strict=True, reason="Finding F3 (SSJ 1.0.0): DifferentiableCombinedExtendedFunction"
                   ".diff defaults to twosided=False, so HetBlock differentiates hetinputs one-sided "
                   "although it asks for two-sided differences, and hetoutputs one-sided even with "
                   "jacobian(..., twosided=True).")
def test_v6_f3_hetinput_hetoutput_two_sided():
    F = _ssj_function
    # (i) the hetinput derivative HetBlock uses (jac_backward_prelim + backward_fakenews)
    e, h = np.array([0.5, 1.0, 1.5]), 1e-4
    d = F.CombinedExtendedFunction([_income_nl]).differentiable(dict(w=1.0, e_grid=e), h, True)
    used = d.diff({'w': 1.0})['y']
    central = (_income_nl(1 + h, e) - _income_nl(1 - h, e)) / (2 * h)
    hetinput_err = float(np.max(np.abs(used - central)))
    # (ii) a hetoutput Jacobian with twosided=True, against the one with the stored default
    cal = dict(min_a=0, max_a=1000, rho_e=0.975, sd_e=0.7, n_a=200, n_e=7, w=1.0, r=0.01 / 4,
               beta=1 - 0.08 / 4, eis=1)
    blk = _hh_sim.add_hetinputs([_income_lin, _make_grids]).add_hetoutputs([_c2])
    ss = blk.steady_state(cal)
    orig = F.DifferentiableCombinedExtendedFunction.diff

    def stored_default(self, shock_dict, h=None, outputs=None, hide_zeros=False, twosided=None):
        return orig(self, shock_dict, h, outputs, hide_zeros, twosided)

    J0 = np.asarray(blk.jacobian(ss, ['r'], ['C2'], T=30, twosided=True, h=1e-4)['C2']['r'])
    F.DifferentiableCombinedExtendedFunction.diff = stored_default
    try:
        J1 = np.asarray(blk.jacobian(ss, ['r'], ['C2'], T=30, twosided=True, h=1e-4)['C2']['r'])
    finally:
        F.DifferentiableCombinedExtendedFunction.diff = orig
    hetoutput_err = float(np.max(np.abs(J0 - J1)) / np.max(np.abs(J1)))
    assert hetinput_err < 1e-9 and hetoutput_err < 1e-9, (hetinput_err, hetoutput_err)
