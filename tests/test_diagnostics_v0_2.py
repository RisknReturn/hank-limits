"""v0.2: audit of FastHetBlock's stationary distribution (R4 review).

A small stationarity residual does not by itself bound the distribution error when the chain
mixes slowly. These tests check what else is reported: no negative mass before the absolute
value is taken, the sweep tolerance is met, the result does not move with the shift tau, and
SSJ's own forward iteration converges towards it as its tolerance is tightened.
"""
import numpy as np
import pytest

from hank_limits import models as M, diagnostics as DG


@pytest.fixture(scope='module')
def audit_t4():
    out = M.ssj_aiyagari_r(0.2, 0.9)
    return DG.distribution_audit(out['_model']['block'], out['_model']['ss'],
                                 ssj_tols=(1e-10, 1e-12), threshold=1.0)


@pytest.fixture(scope='module')
def audit_t1_mid():
    m = M.pih_household(0.1)
    return DG.distribution_audit(m['block'], m['ss'], ssj_tols=(1e-10, 1e-12), threshold=0.0)


@pytest.mark.parametrize('name', ['audit_t4', 'audit_t1_mid'])
def test_inverse_iteration_report(name, request):
    a = request.getfixturevalue(name)
    ref = a['reference']
    assert ref['sweep_tol_met'] and ref['sweeps'] <= 6
    assert ref['neg_mass_rel_max'] < 1e-12
    assert ref['resid_L1'] < 1e-12
    for row in a['shift'].values():                 # tau from 1e-12 to 1e-6
        assert row['L1_diff'] < 1e-11 and row['dA_rel'] < 1e-11


@pytest.mark.parametrize('name', ['audit_t4', 'audit_t1_mid'])
def test_ssj_forward_iteration_converges_to_it(name, request):
    a = request.getfixturevalue(name)
    f10, f12 = a['ssj_forward']['1e-10'], a['ssj_forward']['1e-12']
    assert f10['mirror_identical'] and f12['mirror_identical']
    assert f12['L1_diff'] < f10['L1_diff'] / 10       # tighter SSJ tolerance -> closer to ours
    assert f12['dA_rel'] < 1e-4


def test_ssj_stopping_rule_is_not_an_accuracy_bound(audit_t1_mid):
    # mid-risk PIH household (eps = 0.1, slow mixing): SSJ's default stopping rule
    # (max |Dbeg_{n+1} - Dbeg_n| < 1e-10, checked every 10 sweeps) stops ~1e-3 away in mean assets
    assert audit_t1_mid['ssj_forward']['1e-10']['dA_rel'] > 1e-4
