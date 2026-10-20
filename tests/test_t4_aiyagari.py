"""T4: Aiyagari (1994) steady state.

No-idiosyncratic-risk limit r -> 1/beta - 1 as income risk -> 0 (there are no insurance markets
in the code, so this is not called a complete-markets limit). For THIS calibration: r below
1/beta - 1 and decreasing in risk size and persistence at fixed unconditional sd; these are
regression checks, not general theorems. `r < 1/beta - 1` is enforced by the root bracket and
is kept only as a sanity check. Non-mechanical checks: capital-market residual, Euler errors on
unconstrained states, KKT violations on constrained states, mass in the top grid points.
Published Table II cells: scripts/aiyagari_table2_v0_2.py (README section 7).
"""
import numpy as np


def test_no_idiosyncratic_risk_limit(t4):
    gap = np.array(t4['gap'])
    assert np.all(gap > 0)
    assert np.all(np.diff(gap) < 0)
    assert gap[-1] < 1e-4
    assert t4['fitted_rate_gap'] > 1.8          # observed ~1.9: gap = O(sd^2) on this path


def test_r_below_discount_rate_sanity_only(t4):
    # mechanical for the SSJ provider (root searched below 1/beta - 1); not accuracy evidence
    assert t4['all_below_rho']


def test_market_clearing_residual(t4):
    assert max(t4['capital_market_resid_rel']) < 1e-8


def test_euler_and_kkt_residuals(t4):
    # unit-free Euler errors on unconstrained states; positive Euler gaps at binding states are
    # not errors, only KKT violations (u'(c) < beta R E u'(c') at a binding state) are
    assert max(t4['euler_err_mean_weighted']) < 1e-6
    assert max(t4['euler_err_max_with_mass']) < 1e-4
    assert max(t4['kkt_violation_max_with_mass']) < 1e-10


def test_top_grid_mass_negligible(t4):
    assert max(t4['mass_top_5pts']) < 1e-12
    assert max(t4['mass_policy_off_grid']) < 1e-12


def test_monotone_in_risk_size_regression(t4):
    assert t4['monotone_in_sd']


def test_monotone_in_persistence_regression(t4):
    assert t4['monotone_in_rho']
