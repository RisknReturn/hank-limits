"""T1: permanent-income (PIH) limit of a one-asset HA household block (partial equilibrium;
income risk eps -> 0 along a path on which beta is recalibrated to the asset target; constraint
slack). This is a household-block test, not a test of the full NK general-equilibrium RANK limit.

Theory: with no income risk and beta (1+r) = 1 the household is a permanent-income
consumer whose consumption is linear in wealth, so aggregation is exact and
  dC_t/dY_s = (1-q) q^s,  dC_t/dr_s = (1-q) A q^s + eis C q (1{t>=s} - q^s),  q = 1/(1+r).
With eps > 0 a stationary distribution requires kappa = 1 - beta(1+r) > 0; both kappa and
the precautionary distortion of MPCs are O(eps^2).
"""
import numpy as np


def test_exact_at_eps_zero(t1):
    # eps = 0: SSJ fake-news Jacobian vs closed form; only O(h^2) finite-difference error
    assert t1['eps'][0] == 0.0
    assert t1['err_M'][0] < 1e-8
    assert t1['err_Jr'][0] < 1e-7


def test_errors_decrease_monotonically(t1):
    for key in ('err_M', 'err_Jr'):
        e = t1[key][1:]
        assert all(a > b for a, b in zip(e, e[1:])), (key, e)


def test_second_order_convergence(t1):
    # theory: O(eps^2); observed slopes are ~1.9-2.1
    assert t1['fitted_rate_err_M'] > 1.8
    assert t1['fitted_rate_err_Jr'] > 1.8
    assert min(t1['rates_err_Jr']) > 1.8
    assert min(t1['rates_err_M']) > 1.8


def test_tightest_setting_accuracy(t1):
    assert t1['err_M'][-1] < 2e-3
    assert t1['err_Jr'][-1] < 1e-3


def test_impatience_vanishes_like_eps2(t1):
    eps = np.array(t1['eps'][1:])
    kap = np.array(t1['kappa'][1:])
    assert np.all(np.diff(kap) < 0)
    ratio = kap / eps ** 2
    assert np.all((ratio > 0.01) & (ratio < 0.1)), ratio


def test_borrowing_constraint_slack(t1):
    # near-natural limit: (numerically) nobody at the constraint
    assert max(t1['mass_constrained']) < 1e-5


def test_error_grows_with_observation_window(t1):
    # HA wealth mean-reverts at rate O(eps^2), PIH wealth has a unit root: the gap is
    # concentrated at long lags. err_M_by_window are leading blocks of the SAME T=200 Jacobian
    # (observation window), not Jacobians computed with a shorter horizon.
    e50, e100, e200 = (t1['err_M_by_window']['50'][-1], t1['err_M_by_window']['100'][-1],
                       t1['err_M'][-1])
    assert e50 < e100 < e200


def test_household_jacobian_independent_of_computing_horizon():
    # SSJ fake-news entries J[t, s] do not depend on the computing horizon T: a T=50 Jacobian
    # computed on its own equals the leading 50 x 50 block of the T=200 Jacobian
    from hank_limits import models as M
    m = M.pih_household(0.05)
    J50 = M.household_jacobians(m['block'], m['ss'], 50, ['Y', 'r'], ['C'])
    J200 = M.household_jacobians(m['block'], m['ss'], 200, ['Y', 'r'], ['C'])
    for k in ('Y', 'r'):
        assert np.max(np.abs(J50['C'][k] - J200['C'][k][:50, :50])) < 1e-14
