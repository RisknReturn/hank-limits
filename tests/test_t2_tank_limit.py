"""T2: TANK limit (Bilbiie). Share lam hand-to-mouth, share 1-lam savers without risk.

Theory: aggregate iMPC = lam I + (1-lam) M_RA (equal income incidence); with H-income
elasticity chi and zero liquidity the fixed-real-rate GE output response is
(1-lam)/(1-lam chi) times the RANK response (Bilbiie 2020, aggregate Euler-IS curve).
"""
import numpy as np


def test_htm_type_exact_below_threshold(t2):
    # impatient agent with limit 0 is strictly constrained once risk is small: exact identity
    eps, err = t2['eps'], t2['err_MH_vs_I']
    small = [e for x, e in zip(eps, err) if x <= 0.05]
    assert max(small) < 1e-10
    assert all(a >= b - 1e-12 for a, b in zip(err, err[1:]))


def test_saver_type_converges_to_RA(t2):
    e = t2['err_MS_vs_MRA']
    assert all(a > b for a, b in zip(e, e[1:]))
    assert t2['fitted_rate_err_MS_vs_MRA'] > 1.8
    assert e[-1] < 2e-3


def test_aggregate_impc_is_tank(t2):
    e = t2['err_Magg_vs_TANK']
    assert all(a > b for a, b in zip(e, e[1:]))
    assert e[-1] < 2e-4
    for chi, errs in t2['err_Magg_chi_vs_TANK'].items():
        assert errs[-1] < 2e-4, chi


def test_ge_irf_converges_to_bilbiie(t2):
    for chi, errs in t2['err_IRF_vs_TANK'].items():
        assert all(a > b for a, b in zip(errs, errs[1:])), chi
        assert errs[-1] < 2e-3, chi
    # observed ~1.7-1.8, limited by the asset grid at the smallest eps (see README)
    assert t2['fitted_rate_err_IRF_chi1.5'] > 1.5
    assert t2['fitted_rate_err_IRF_chi1.0'] > 1.5


def test_amplification_factor(t2):
    target = t2['amplification_factor']['1.5']
    assert abs(t2['amplification_HA']['1.5'][-1] - target) < 1e-3


def test_ge_solution_residual(t2):
    assert max(t2['ge_target_resid']) < 1e-10
