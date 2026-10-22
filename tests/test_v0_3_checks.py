"""v0.3 checks V1-V4 and V6 (docs/CHECKS_PLAN_v0_3.md, section 3; V5 is in test_t6_differencing.py),
as revised for v0.3.0 (docs/CHECKS_v0_3_0_AMENDMENT.md).

The quantities are computed in `hank_limits.checks` (numpy only) from primitive quantities the
providers report in the optional fields `ss_detail`, `J_A`, `J_C`, `policy` and the T5 paths, and
from the inputs the check passed, against the documented values in `hank_limits.spec`, and stored
under the key 'v0_3' of each check result. Every V test judges its quantities with
`checks.v03_verdict`: a non-finite, negative or too large value fails the test whatever else is
missing; a test with no failure but a missing value is skipped as "incomplete" (never passed); a
test none of whose quantities the provider reports is skipped, with the reason.

Thresholds (plan, section 4): threshold(q) = roundup_1sig(max(100 * b_max(q), 10 * tau(q))),
b_max = largest value over all settings in two runs of the unmodified suite before commit T
(quantities added in v0.3.0: before commit B of the amendment, its section 5), tau = tolerance of
the numerical step that determines q. These are regression thresholds for this calibration, like
the tolerances of T1-T6.
"""
import numpy as np
import pytest

from sequence_jacobian.hetblocks.hh_sim import hh as _hh_sim, make_grids as _make_grids
from sequence_jacobian.utilities import function as _ssj_function

from hank_limits import checks as C

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

# v0.3.0 quantities (amendment, section 5; two unmodified runs of the commit-B candidate):
# V1 on a goods-market GE target (`ssj_werning_provider(ge_target='goods_mkt')`), and the T5
# exogenous paths against the request (R3). The 88 thresholds above are unchanged.
THRESH['t3'].update({
    'ge_target_resid_reported_goods_acyc': 9e-12, 'ge_target_resid_reported_goods_cyc': 2e-11,
    'ge_target_resid_direct_goods_acyc': 9e-12, 'ge_target_resid_direct_goods_cyc': 2e-11,
})
THRESH['t5_paths'].update({'exog_G': 1e-13, 'exog_r': 1e-13, 'exog_rshock': 1e-13})

FISCAL = 'ssj_consistent_refined'     # T5 setting used by the V checks (as the Walras test)


def _t3(*names):
    return [f'{n}_{t}' for n in names for t in ('acyc', 'cyc')]


_HH_BUDGET = {'t1': ['hh_budget'], 't2': ['hh_budget_H', 'hh_budget_S'],
              't3': _t3('hh_budget'), 't4': ['hh_budget'], 't5': ['hh_budget']}

_INC = ('pi_row_sum_err', 'pi_negative', 'mean_e_err', 'agg_income_err')
_INCOME = {'t1': list(_INC), 't2': [f'{k}_{t}' for t in 'HS' for k in _INC] + ['incidence_err'],
           't3': [f'{k}_{t}' for t in ('acyc', 'cyc') for k in _INC] + ['incidence_err_cyc'],
           't4': list(_INC), 't5': list(_INC)}

_SPEC = {'t1': ['spec_r', 'spec_eis', 'spec_A', 'spec_a_min'],
         't2': ['spec_lam', 'spec_beta_H', 'spec_r', 'spec_eis', 'spec_a_min_H', 'spec_a_min_S'],
         't3': ['spec_eis_acyc', 'spec_eis_cyc', 'spec_a_min_acyc', 'spec_a_min_cyc',
                'spec_zeta_acyc', 'spec_zeta_low_cyc', 'spec_rstar_acyc', 'spec_rstar_cyc'],
         't4': ['spec_beta', 'spec_eis', 'spec_a_min'],
         't5': ['spec_r', 'spec_eis', 'spec_B', 'spec_G', 'spec_phi_B', 'spec_a_min']}

# Every V test (group V of the mutation runner) -> (check result, quantities, threshold table).
# The tests below judge exactly these; tests/test_v0_3_0_checker.py runs the same list on fake
# providers. Threshold table 't5_paths': T5 path quantities, over all shock labels.
V_TESTS = {
    'test_v1_werning_ge_target_resid_reported':
        ('t3', _t3('ge_target_resid_reported', 'ge_target_resid_reported_goods'), 't3'),
    'test_v1_werning_ge_target_resid_direct':
        ('t3', _t3('ge_target_resid_direct', 'ge_target_resid_direct_goods'), 't3'),
    'test_v2_fiscal_ss_government_budget': ('t5', ['ss_government_budget'], 't5'),
    'test_v2_fiscal_ss_goods_market': ('t5', ['ss_goods'], 't5'),
    'test_v2_fiscal_ss_asset_market': ('t5', ['ss_assets'], 't5'),
    'test_v2_tank_ss_goods_market': ('t2', ['ss_goods'], 't2'),
    'test_v2_tank_ss_asset_market': ('t2', ['ss_assets'], 't2'),
    'test_v2_werning_ss_asset_market': ('t3', _t3('ss_assets'), 't3'),
    'test_v2_werning_ss_goods_market': ('t3', _t3('ss_goods'), 't3'),
    'test_v2_aiyagari_ss_capital_market': ('t4', ['ss_capital'], 't4'),
    'test_v2_aiyagari_ss_firm_focs': ('t4', ['foc_r', 'foc_w'], 't4'),
    'test_v2_aiyagari_ss_goods_market': ('t4', ['ss_goods'], 't4'),
    **{f'test_v2_household_budget[{n}]': (n, k, n) for n, k in _HH_BUDGET.items()},
    **{f'test_v3_income_normalisation[{n}]': (n, k, n) for n, k in _INCOME.items()},
    **{f'test_v4_spec_parameters[{n}]': (n, k, n) for n, k in _SPEC.items()},
    'test_v4_aiyagari_euler_at_spec_beta':
        ('t4', ['spec_euler_err_mean_weighted', 'spec_euler_err_max_with_mass',
                'spec_kkt_violation_max_with_mass'], 't4'),
    'test_v4_fiscal_rule_on_path': ('t5', ['fiscal_rule'], 't5_paths'),
    'test_v4_government_budget_on_path': ('t5', ['government_budget'], 't5_paths'),
    'test_v4_real_rate_rule_on_path': ('t5', ['rate_rule'], 't5_paths'),
    'test_v4_goods_market_on_path': ('t5', ['goods_market'], 't5_paths'),
    'test_v4_exogenous_G_as_requested': ('t5', ['exog_G'], 't5_paths'),
    'test_v4_exogenous_r_as_requested': ('t5', ['exog_r'], 't5_paths'),
    'test_v4_exogenous_rshock_as_requested': ('t5', ['exog_rshock'], 't5_paths'),
}


def _v03(res):
    if 'v0_3' not in res:
        pytest.skip('check result has no v0.3 quantities')
    if res['v0_3'].get('errors'):
        pytest.fail(f"v0.3 computation failed: {res['v0_3']['errors'][:3]}")
    return res['v0_3']


def _assert_below(v, keys, fx, paths=False):
    """One V test on the quantities `keys` of the v0.3 record `v`, thresholds THRESH[fx]
    (`checks.v03_verdict`): fails on any non-finite, negative or too large value, whatever else
    is missing; skips as "incomplete" if values are missing without failure; skips with the
    recorded reasons if the provider reports none of the quantities."""
    r = C.v03_verdict(v, keys, {k: THRESH[fx][k] for k in keys}, paths=paths)
    q = {k: p for k, p in r['quantities'].items() if p['status'] != 'n/a'}
    if r['errors']:
        pytest.fail(f"v0.3 computation failed: {r['errors'][:3]}")
    bad = {k: [(w, x, p['threshold']) for w, x in p['failed']][:4]
           for k, p in q.items() if p['status'] == 'failed'}
    missing = {k: p['missing'][:6] for k, p in q.items() if p['missing']}
    assert not bad, (bad, f'missing: {missing}' if missing else '')
    if r['status'] == 'incomplete':
        pytest.skip(f"incomplete (no failure, values missing): {missing}; "
                    f"not reported: {[k for k, p in q.items() if p['status'] == 'skipped']}; "
                    f"skipped: {r['skipped'][:6]}")
    if r['status'] == 'skipped':
        pytest.skip(f"not reported: {', '.join(keys)}; skipped: {r['skipped'][:6]}")


def _assert_paths(t5, key):
    """V4(b) path quantity `key` of the T5 setting FISCAL, over every shock label the check
    expected (every experiment the caller requested) and every T."""
    _assert_below(_v03(t5[FISCAL]), [key], 't5_paths', paths=True)


def _judge(name, res):
    fx, keys, table = V_TESTS[name]
    if table == 't5_paths':
        return _assert_paths(res, keys[0])
    _assert_below(_v03(res[FISCAL] if fx == 't5' else res), keys, table)


# ---------------------------------------------------------------------------
# V1: T3 GE target residual, divided by b. The target is the market the provider reports in
# `ge_target` ('asset_mkt' assumed if absent): asset market J^{A,Y} dY + J^{A,r} dr, goods market
# dY - J^{C,Y} dY - J^{C,r} dr (v0.3.0, R4); the other market is recorded as a diagnostic.
# ---------------------------------------------------------------------------

def test_v1_werning_ge_target_resid_reported(t3):
    # the residual of its target that the T3 provider reports, recorded as in T2
    _judge('test_v1_werning_ge_target_resid_reported', t3)


def test_v1_werning_ge_target_resid_direct(t3):
    # the residual of the reported target, from the reported Jacobians and the rate path passed
    _judge('test_v1_werning_ge_target_resid_direct', t3)


# ---------------------------------------------------------------------------
# V2: steady-state residuals
# ---------------------------------------------------------------------------

def test_v2_fiscal_ss_government_budget(t5):
    _judge('test_v2_fiscal_ss_government_budget', t5)


def test_v2_fiscal_ss_goods_market(t5):
    _judge('test_v2_fiscal_ss_goods_market', t5)


def test_v2_fiscal_ss_asset_market(t5):
    _judge('test_v2_fiscal_ss_asset_market', t5)


def test_v2_tank_ss_goods_market(t2):
    _judge('test_v2_tank_ss_goods_market', t2)


def test_v2_tank_ss_asset_market(t2):
    _judge('test_v2_tank_ss_asset_market', t2)


def test_v2_werning_ss_asset_market(t3):
    _judge('test_v2_werning_ss_asset_market', t3)


def test_v2_werning_ss_goods_market(t3):
    _judge('test_v2_werning_ss_goods_market', t3)


def test_v2_aiyagari_ss_capital_market(t4):
    _judge('test_v2_aiyagari_ss_capital_market', t4)


def test_v2_aiyagari_ss_firm_focs(t4):
    # documented alpha, delta; L = mean of the reported income states
    _judge('test_v2_aiyagari_ss_firm_focs', t4)


def test_v2_aiyagari_ss_goods_market(t4):
    _judge('test_v2_aiyagari_ss_goods_market', t4)


@pytest.mark.parametrize('name', sorted(_HH_BUDGET))
def test_v2_household_budget(name, request):
    # |C - r A - sum_e pi_e y_e| / C, pi the stationary distribution of the reported Pi
    _judge(f'test_v2_household_budget[{name}]', request.getfixturevalue(name))


# ---------------------------------------------------------------------------
# V3: income normalisation
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('name', sorted(_INCOME))
def test_v3_income_normalisation(name, request):
    # Pi row-stochastic, E_pi[e] = 1, aggregate income = Y (T2 per type and incidence; T4
    # w x documented labour 1; T5 Y - Tax); T3 countercyclical: sum pi e zeta = 0
    _judge(f'test_v3_income_normalisation[{name}]', request.getfixturevalue(name))


# ---------------------------------------------------------------------------
# V4: specification pins
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('name', sorted(_SPEC))
def test_v4_spec_parameters(name, request):
    # reported parameters against hank_limits.spec (T3: r* from the reported beta, Pi, e)
    _judge(f'test_v4_spec_parameters[{name}]', request.getfixturevalue(name))


def test_v4_aiyagari_euler_at_spec_beta(t4):
    # Euler errors of the reported policy at the documented beta and EIS
    _judge('test_v4_aiyagari_euler_at_spec_beta', t4)


def test_v4_fiscal_rule_on_path(t5):
    # dTax_t = phi_B dB_{t-1}, documented phi_B
    _judge('test_v4_fiscal_rule_on_path', t5)


def test_v4_government_budget_on_path(t5):
    # dB_t = (1 + r) dB_{t-1} + B dr_t + dG_t - dTax_t, documented r and B
    _judge('test_v4_government_budget_on_path', t5)


def test_v4_real_rate_rule_on_path(t5):
    # dr_t = phi_y dY_t + rshock_t (rule-based shocks), documented phi_y; rshock as requested
    _judge('test_v4_real_rate_rule_on_path', t5)


def test_v4_goods_market_on_path(t5):
    # dY_t = dC_t + dG_t
    _judge('test_v4_goods_market_on_path', t5)


# v0.3.0 (R3): the exogenous paths the provider solved, against the paths the check requested
# (zero where an experiment requests none)

def test_v4_exogenous_G_as_requested(t5):
    _judge('test_v4_exogenous_G_as_requested', t5)


def test_v4_exogenous_r_as_requested(t5):
    # fixed-rate experiments; in rule-based ones dr is endogenous (rate rule)
    _judge('test_v4_exogenous_r_as_requested', t5)


def test_v4_exogenous_rshock_as_requested(t5):
    _judge('test_v4_exogenous_rshock_as_requested', t5)


# ---------------------------------------------------------------------------
# V6: finding F3 (SSJ code only), split in two (v0.3.0): if only one part is fixed upstream,
# only that test turns XPASS
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


_F3 = ("Finding F3 (SSJ 1.0.0): DifferentiableCombinedExtendedFunction.diff defaults to "
       "twosided=False instead of None, so the stored setting is ignored when the caller omits it. ")


# strict xfail = regression pin for sequence-jacobian 1.0.0 only (as F1, F2)
@pytest.mark.xfail(strict=True, reason=_F3 + "V6a: HetBlock differentiates hetinputs one-sided "
                   "although it asks for two-sided differences.")
def test_v6a_f3_hetinput_two_sided():
    # the hetinput derivative HetBlock uses (jac_backward_prelim + backward_fakenews)
    F = _ssj_function
    e, h = np.array([0.5, 1.0, 1.5]), 1e-4
    d = F.CombinedExtendedFunction([_income_nl]).differentiable(dict(w=1.0, e_grid=e), h, True)
    used = d.diff({'w': 1.0})['y']
    central = (_income_nl(1 + h, e) - _income_nl(1 - h, e)) / (2 * h)
    hetinput_err = float(np.max(np.abs(used - central)))
    assert hetinput_err < 1e-9, hetinput_err


@pytest.mark.xfail(strict=True, reason=_F3 + "V6b: hetoutputs are differentiated one-sided even "
                   "with jacobian(..., twosided=True).")
def test_v6b_f3_hetoutput_two_sided():
    # a hetoutput Jacobian with twosided=True, against the one with the stored default
    F = _ssj_function
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
    assert hetoutput_err < 1e-9, hetoutput_err
