"""Regression tests of the v0.3 checker itself, with fake providers (v0.3.0 amendment,
docs/CHECKS_v0_3_0_AMENDMENT.md, sections 3.1-3.4). No model code is used, so the mutation runner
maps this file to no group.

R1  a missing or non-finite value cannot hide a verdict (the real helpers of
    tests/test_v0_3_checks.py and `checks.v03_verdict`);
R2  a provider with exactly the documented fields (`hank_limits.interface`) gets a verdict on
    every V test without error; one with the required fields only gets skips with reasons;
R3  the solved exogenous paths are compared with the requested ones;
R4  V1 is judged on the GE target the provider reports.
"""
import numpy as np
import pytest

from hank_limits import checks as C
from hank_limits import spec as SPEC
from hank_limits.ra import ra_jacobians, tank_output_fixed_real_rate, zero_liquidity_output

import test_v0_3_checks as V

E = np.array([0.5, 1.5])                       # income states, E_pi[e] = 1
PI = np.array([[0.9, 0.1], [0.1, 0.9]])        # stationary distribution (0.5, 0.5)
DR = -0.01 * 0.8 ** np.arange(8)


def _fiscal_shocks(T):
    path = 0.01 * 0.8 ** np.arange(T)
    return {'G_debt_fixed_r': ('fixed_r', {'G': path}), 'r_fixed_r': ('fixed_r', {'r': -path}),
            'G_debt_r_rule': ('r_rule', {'G': path}), 'rshock_r_rule': ('r_rule', {'rshock': -path})}


# ---------------------------------------------------------------------------
# Fake providers: `full=True` returns exactly the documented required and optional fields,
# `full=False` the required ones only (hank_limits.interface). Steady states are consistent, so
# the full fakes pass most V checks; what matters here is that every check reaches a verdict.
# ---------------------------------------------------------------------------

def fake_rank(eps, T, full=True):
    r, eis, A, Y = 0.03, 0.5, 5.0, 1.0
    C_ = Y + r * A
    M, Jr = ra_jacobians(T, r, eis, C_, A)
    d = eps ** 2 + 1e-6
    out = dict(J={'C': {'Y': M + d, 'r': Jr + d}}, ss=dict(C=C_, A=A, r=r, eis=eis))
    if full:
        out['diag'] = dict(kappa=eps ** 2, mass_constrained=0.0, ssj_forward_resid=0.0)
        out['ss_detail'] = dict(Y=Y, C=C_, A=A, r=r, eis=eis, a_min=-0.99 * Y * E.min() / r,
                                e_grid=E, Pi=PI, y=Y * E)
    return out


def fake_tank(eps, T, dr, chis=(1.0, 1.5), full=True):
    lam, r, eis = 0.3, 0.03, 0.5
    M_ra, _ = ra_jacobians(T, r, eis, 1.0, 0.0)
    d = eps ** 2
    out = dict(M_H=np.eye(T) + d, M_S=M_ra + d, lam=lam,
               dY={float(c): tank_output_fixed_real_rate(dr, r, eis, lam, c) * (1 + d) for c in chis},
               ss=dict(r=r, eis=eis, C_S=1.0, A_S=0.0, Y=1.0))
    if full:
        out['diag'] = dict(unconstrained_H=0.0, A_H=0.0, kappa_S=d, mass_constrained_S=0.0,
                           ge_target_resid=0.0)
        det = dict(lam=lam, r=r, eis=eis, Y=1.0, Y_H=1.0, Y_S=1.0, beta_H=0.95)
        for t, a_min in (('H', 0.0), ('S', -0.99 * E.min() / r)):
            det.update({f'C_{t}': 1.0, f'A_{t}': 0.0, f'a_min_{t}': a_min, f'e_grid_{t}': E,
                        f'Pi_{t}': PI, f'y_{t}': E})
        out['ss_detail'] = det
    return out


def fake_werning(b, T, dr, cyclical=False, full=True, target=('asset_mkt', 'goods_mkt'),
                 solve_on=None):
    """target: the GE target reported (acyclical, countercyclical); None = not reported.
    solve_on: the market whose residual the returned dY satisfies (default: both)."""
    eis, gam = 0.5, 2.0
    zeta = np.zeros(2)
    if cyclical:
        zeta = np.array([1.0, -0.5 * E[0] / (0.5 * E[1])])
    ratio = (PI @ E ** (-gam)) / E ** (-gam)
    beta = 1.0 / ((1 + 0.01) * ratio.max())
    r_ss = 0.01 + 2 * b
    dY = zero_liquidity_output(dr, r_ss, eis, E, PI, zeta)[0] * (1 + b)
    out = dict(dY=dY, r_ss=r_ss, eis=eis, e_grid=E, Pi=PI, zeta=zeta,
               diag=dict(r_gap_over_b=2.0, mass_constrained=0.0, max_wealth_over_b=1.0))
    if full:
        # Jacobians constructed so that dY clears the asset market and the goods market
        k = dr / (dr @ dr)
        JA_Y, JC_Y = b * np.eye(T), 0.5 * np.eye(T)
        JA_r = np.outer(-JA_Y @ dY, k) if solve_on in (None, 'asset_mkt') else np.zeros((T, T))
        JC_r = np.outer(dY - JC_Y @ dY, k) if solve_on in (None, 'goods_mkt') else np.zeros((T, T))
        out.update(J_A={'Y': JA_Y, 'r': JA_r}, J_C={'Y': JC_Y, 'r': JC_r},
                   ss_detail=dict(Y=1.0, C=1.0, A=0.0, r=r_ss, a_min=-b, y=E, beta=beta))
        out['diag']['ge_target_resid'] = 0.0
        t = target[int(cyclical)]
        if t is not None:
            out['ge_target'] = t
    return out


def fake_aiyagari(sd, rho, full=True):
    beta, al, de = 0.96, 0.36, 0.08
    r = 1 / beta - 1 - 0.5 * sd ** 2 * (1 + rho)
    out = dict(r=r, rho_disc=1 / beta - 1, diag=dict(mass_constrained=0.1))
    if full:
        K = (al / (r + de)) ** (1 / (1 - al))
        w = (1 - al) * K ** al
        out['diag'].update(capital_market_resid_rel=0.0, mass_top_5pts=0.0, mass_policy_off_grid=0.0,
                           euler=dict(euler_err_mean_weighted=0.0, euler_err_max_with_mass=0.0,
                                      kkt_violation_max_with_mass=0.0))
        out['ss_detail'] = dict(r=r, w=w, K=K, A=K, C=K ** al - de * K, beta=beta, eis=0.5,
                                a_min=0.0, e_grid=E, Pi=PI, y=w * E)
        out['policy'] = dict(a_grid=np.array([0.0, 1.0, 2.0]), a=np.tile([0.5, 1.0, 1.5], (2, 1)),
                             c=np.ones((2, 3)), D=np.full((2, 3), 1 / 6))
    return out


def fake_fiscal(T, shocks, full=True, factor=1.0, drop=()):
    """Linear fiscal HANK (income (Y - Tax) e, fixed r or r_t = r + phi_y dY_t + rshock_t) solved
    exactly, with every exogenous input multiplied by `factor` (R3: factor 2 = a provider that
    applies twice the requested shock and returns consistent paths). `drop`: labels not returned."""
    r, B, G, phi_B, phi_y, eis, Y = 0.02, 1.0, 0.2, 0.1, 0.25, 0.5, 1.0
    Tax, C_, A = G + r * B, Y - G, B
    M, Jr = ra_jacobians(T, r, eis, C_, A)
    JT, I, L = -M, np.eye(T), np.eye(T, k=-1)
    debt = I - (1 + r) * L + phi_B * L
    irfs = {}
    for label, (kind, inputs) in shocks.items():
        if label in drop:
            continue
        z = np.zeros(T)
        dG, rx, rs = (factor * np.asarray(inputs.get(k, z), float) for k in ('G', 'r', 'rshock'))
        if kind == 'fixed_r':
            dr = rx
            dB = np.linalg.solve(debt, B * dr + dG)
            dTax = phi_B * L @ dB
            dY = np.linalg.solve(I - M, dG + JT @ dTax + Jr @ dr)
        else:
            H = np.block([[I - M - phi_y * Jr, -phi_B * JT @ L], [-B * phi_y * I, debt]])
            x = np.linalg.solve(H, np.concatenate([dG + Jr @ rs, B * rs + dG]))
            dY, dB = x[:T], x[T:]
            dr, dTax = phi_y * dY + rs, phi_B * L @ dB
        irfs[label] = dict(dY=dY, dG=dG, dTax=dTax, dr=dr, asset_mkt=np.zeros(T))
        if full:
            irfs[label].update(dC=dY - dG, dB=dB, drshock=rs)
    out = dict(J={'C': {'Y': M, 'Tax': JT, 'r': Jr}}, ss=dict(C=C_, A=A, r=r, eis=eis), irfs=irfs)
    if full:
        out['ss_detail'] = dict(Y=Y, C=C_, A=A, B=B, G=G, Tax=Tax, r=r, eis=eis, phi_B=phi_B,
                                a_min=0.0, e_grid=E, Pi=PI, y=(Y - Tax) * E)
    return out


def run_fake(model, full=True, **kw):
    """The check of `model` on its fake provider; returns the result in the form of the fixture."""
    if model == 't1':
        return C.check_rank_limit(lambda eps, T: fake_rank(eps, T, full), [0.0, 0.2, 0.1], 8,
                                  T_sub=(4, 6), spec=SPEC.T1)
    if model == 't2':
        return C.check_tank_limit(lambda eps, T, dr, chis: fake_tank(eps, T, dr, chis, full),
                                  [0.2, 0.1], 8, DR, spec=SPEC.T2)
    if model == 't3':
        return C.check_werning(lambda b, T, dr, cyclical: fake_werning(b, T, dr, cyclical, full, **kw),
                               [1e-2, 1e-3], 8, DR, asymptotic_from=1e-2, spec=SPEC.T3)
    if model == 't4':
        return C.check_aiyagari(lambda sd, rho: fake_aiyagari(sd, rho, full), [0.2, 0.1], 0.9,
                                [0.1, 0.2], [0.0, 0.5], 0.2, 0.9, spec=SPEC.T4)
    return {V.FISCAL: C.check_ikc(lambda T, shocks: fake_fiscal(T, shocks, full, **kw), [8, 12],
                                  _fiscal_shocks, spec=SPEC.T5)}


def verdicts(model, res):
    """{V test: v03_verdict} for every V test that reads `model` (tests/test_v0_3_checks.V_TESTS)."""
    v = res[V.FISCAL]['v0_3'] if model == 't5' else res['v0_3']
    return {name: C.v03_verdict(v, keys, {k: V.THRESH[table][k] for k in keys},
                                paths=(table == 't5_paths'))
            for name, (fx, keys, table) in V.V_TESTS.items() if fx == model}


MODELS = ['t1', 't2', 't3', 't4', 't5']


def test_registry_lists_every_v_test():
    funcs = {n for n in dir(V) if n.startswith('test_v') and not n.startswith('test_v6')}
    assert {n.split('[')[0] for n in V.V_TESTS} == funcs


# ---------------------------------------------------------------------------
# R1
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('vals', [[0.0, np.nan], [np.nan, 0.0], [0.0, np.inf], [0.0, -1e-3]])
def test_r1_nonfinite_or_negative_value_fails(vals):
    # v0.3: [0.0, NaN] passed because max() kept the leading 0.0
    with pytest.raises(AssertionError):
        V._assert_below({'hh_budget': vals}, ['hh_budget'], 't1')
    assert C.v03_verdict({'hh_budget': vals}, ['hh_budget'], {'hh_budget': 1.0})['status'] == 'failed'


def test_r1_missing_value_is_incomplete_not_passed():
    # v0.3: [0.0, None] passed because None was dropped
    with pytest.raises(pytest.skip.Exception, match='^incomplete'):
        V._assert_below({'hh_budget': [0.0, None]}, ['hh_budget'], 't1')
    assert C.v03_verdict({'hh_budget': [0.0, None]}, ['hh_budget'],
                         {'hh_budget': 1.0})['status'] == 'incomplete'


def test_r1_missing_quantity_does_not_hide_a_failure():
    # v0.3: the proportional-income value absent -> the whole V1 test skipped, 10.0 unseen
    v = {'ge_target_resid_direct_cyc': [10.0]}
    with pytest.raises(AssertionError):
        V._assert_below(v, ['ge_target_resid_direct_acyc', 'ge_target_resid_direct_cyc'], 't3')
    with pytest.raises(AssertionError):
        V._judge('test_v1_werning_ge_target_resid_direct', {'v0_3': v})
    r = C.v03_verdict(v, ['ge_target_resid_direct_acyc', 'ge_target_resid_direct_cyc'],
                      {'ge_target_resid_direct_acyc': 2e-12, 'ge_target_resid_direct_cyc': 2e-12})
    assert r['status'] == 'failed'
    assert r['quantities']['ge_target_resid_direct_acyc']['status'] == 'skipped'


def test_r1_missing_quantity_without_failure_is_incomplete():
    v = {'ge_target_resid_direct_cyc': [0.0]}
    with pytest.raises(pytest.skip.Exception, match='^incomplete'):
        V._assert_below(v, ['ge_target_resid_direct_acyc', 'ge_target_resid_direct_cyc'], 't3')


def test_r1_nothing_reported_is_skipped_with_reason():
    with pytest.raises(pytest.skip.Exception, match='^not reported'):
        V._assert_below({'skipped': ['hh_budget (T1): provider does not report ss_detail']},
                        ['hh_budget'], 't1')


def test_r1_path_nonfinite_fails():
    # v0.3: _assert_paths took max() per label as well
    t5 = {V.FISCAL: {'v0_3': {'paths': {'G': {'goods_market': [0.0, np.nan]},
                                        'r': {'goods_market': [0.0]}}}}}
    with pytest.raises(AssertionError):
        V._assert_paths(t5, 'goods_market')


def test_r1_path_of_requested_experiment_missing_is_incomplete():
    # v0.3: _paths collected only the labels the provider returned
    res = run_fake('t5', drop=('rshock_r_rule',))
    for key in ('goods_market', 'rate_rule', 'exog_G'):
        with pytest.raises(pytest.skip.Exception, match='^incomplete'):
            V._assert_paths(res, key)
    st = res[V.FISCAL]['v0_3']['status']
    assert st['paths.rshock_r_rule.goods_market']['state'] == 'not_reported'
    assert st['paths.G_debt_fixed_r.goods_market']['state'] == 'complete'


def test_r1_coverage_records():
    v = run_fake('t1')['v0_3']
    assert v['settings'] == ['eps=0', 'eps=0.2', 'eps=0.1']
    assert v['at']['hh_budget'] == v['settings']
    assert v['status']['hh_budget'] == dict(state='complete', n_expected=3, n_reported=3,
                                            missing=[], nonfinite=[])


# ---------------------------------------------------------------------------
# R2
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('model', MODELS)
def test_r2_documented_fields_reach_a_verdict(model):
    res = run_fake(model, full=True)
    v = res[V.FISCAL]['v0_3'] if model == 't5' else res['v0_3']
    assert not v.get('errors') and not v['skipped'], (v.get('errors'), v['skipped'])
    out = verdicts(model, res)
    assert out and all(r['status'] in ('passed', 'failed') for r in out.values()), \
        {n: (r['status'], r['quantities']) for n, r in out.items() if r['status'] not in ('passed', 'failed')}


# with the required fields only, these V tests are still computable (T5 paths: dr, dY, dG and
# the caller's inputs)
REQUIRED_ONLY_JUDGED = {'test_v4_real_rate_rule_on_path', 'test_v4_exogenous_G_as_requested',
                        'test_v4_exogenous_r_as_requested'}


@pytest.mark.parametrize('model', MODELS)
def test_r2_required_fields_only_skip_with_reason(model):
    res = run_fake(model, full=False)
    v = res[V.FISCAL]['v0_3'] if model == 't5' else res['v0_3']
    assert not v.get('errors'), v.get('errors')
    for name, r in verdicts(model, res).items():
        if name in REQUIRED_ONLY_JUDGED:
            assert r['status'] in ('passed', 'failed'), (name, r)
        else:
            assert r['status'] in ('skipped', 'incomplete') and r['skipped'], (name, r)


def test_r2_t5_path_equation_skipped_alone():
    # a provider without dB and dC: the fiscal rule, budget and goods market are skipped with the
    # field named; no KeyError, and the other path quantities are still computed
    res = run_fake('t5', full=False)
    v = res[V.FISCAL]['v0_3']
    assert any("goods_market (T5 path G_debt_fixed_r): provider does not report "
               "irfs['G_debt_fixed_r']['dC']" in s for s in v['skipped'])
    assert v['status']['paths.G_debt_fixed_r.exog_G']['state'] == 'complete'


def test_r2_werning_beta_from_ss_detail():
    # interface: beta is in ss_detail; the fake has no top-level beta
    res = run_fake('t3')
    assert res['v0_3']['status']['spec_rstar_acyc']['state'] == 'complete'
    assert max(res['v0_3']['spec_rstar_acyc']) < 1e-13


# ---------------------------------------------------------------------------
# R3
# ---------------------------------------------------------------------------

def test_r3_doubled_shock_provider_fails_exogenous_checks():
    """A linear provider that applies twice the requested G, r and rshock and returns paths
    consistent with that passes the IKC check and the fiscal-rule, budget and goods-market path
    checks; the exog checks and the rate rule (with the requested rshock) fail."""
    res = run_fake('t5', factor=2.0)
    assert res[V.FISCAL]['max_err'] < 1e-10
    for key in ('fiscal_rule', 'government_budget', 'goods_market'):
        V._assert_paths(res, key)
    for key in ('exog_G', 'exog_r', 'exog_rshock', 'rate_rule'):
        with pytest.raises(AssertionError):
            V._assert_paths(res, key)


def test_r3_consistent_provider_passes():
    res = run_fake('t5')
    for key in ('fiscal_rule', 'government_budget', 'goods_market', 'rate_rule', 'exog_G',
                'exog_r', 'exog_rshock'):
        V._assert_paths(res, key)


# ---------------------------------------------------------------------------
# R4
# ---------------------------------------------------------------------------

def test_r4_target_judged_by_report():
    # acyclical solve reported on the asset market, countercyclical on the goods market
    res = run_fake('t3')
    v = res['v0_3']
    assert set(v['ge_target']['acyc']) == {'asset_mkt'} and set(v['ge_target']['cyc']) == {'goods_mkt'}
    assert 'ge_target_resid_direct_acyc' in v and 'ge_target_resid_direct_goods_cyc' in v
    assert 'ge_target_resid_direct_cyc' not in v and 'ge_target_resid_direct_goods_acyc' not in v
    assert set(v['diagnostics']) == {'at', 'goods_mkt_resid_direct_acyc', 'asset_mkt_resid_direct_cyc'}
    for name in ('test_v1_werning_ge_target_resid_reported', 'test_v1_werning_ge_target_resid_direct'):
        V._judge(name, res)


def test_r4_reported_target_is_the_one_checked():
    # dY clears the goods market only: passes when reported as a goods-market solve, fails when
    # reported as an asset-market solve
    good = run_fake('t3', target=('goods_mkt', 'goods_mkt'), solve_on='goods_mkt')
    V._judge('test_v1_werning_ge_target_resid_direct', good)
    bad = run_fake('t3', target=('asset_mkt', 'asset_mkt'), solve_on='goods_mkt')
    with pytest.raises(AssertionError):
        V._judge('test_v1_werning_ge_target_resid_direct', bad)


def test_r4_absent_target_assumed_asset_market():
    res = run_fake('t3', target=(None, None))
    assert set(res['v0_3']['ge_target']['acyc']) == {'asset_mkt (assumed)'}
    V._judge('test_v1_werning_ge_target_resid_direct', res)


def test_r4_unknown_target_not_applicable():
    res = run_fake('t3', target=('bond_mkt', 'bond_mkt'))
    assert any('applicability condition not met' in s for s in res['v0_3']['skipped'])
    with pytest.raises(pytest.skip.Exception, match='^not reported'):
        V._judge('test_v1_werning_ge_target_resid_direct', res)
    mixed = run_fake('t3', target=('asset_mkt', 'bond_mkt'))
    with pytest.raises(pytest.skip.Exception, match='^incomplete'):
        V._judge('test_v1_werning_ge_target_resid_direct', mixed)
