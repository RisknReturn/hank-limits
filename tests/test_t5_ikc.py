"""T5: intertemporal Keynesian cross (Auclert-Rognlie-Straub, NBER w25020).

Model: one-asset HANK, income (Y - Tax) e, debt-financed G, fixed real rate or a real-rate
rule r = r_ss + phi_y (Y - 1) + shock. The only HA block is the household, so the SSJ GE
output path must satisfy dY = (I - M)^{-1} [dG + J^{C,Tax} dTax + J^{C,r} dr].
"""
import pytest

TOL = 1e-10


def test_ikc_identity_refined_solve(t5):
    # SSJ blocks + SSJ's algorithm with iterative refinement of the linear solve
    res = t5['ssj_consistent_refined']
    assert res['max_err'] < TOL, res['err']


def test_ikc_identity_ssj_solver_fixed_rate(t5):
    res = t5['ssj_consistent']
    for lab in ('G_debt_fixed_r', 'r_fixed_r'):
        assert max(res['err'][lab]) < TOL, (lab, res['err'][lab])


def test_ikc_identity_ssj_solver_rule_short_horizon(t5):
    res = t5['ssj_consistent']
    for lab in ('G_debt_r_rule', 'rshock_r_rule'):
        assert res['err'][lab][0] < TOL, (lab, res['err'][lab])   # T = 100


# strict xfail = regression pin for sequence-jacobian 1.0.0 only; a general conformance suite
# should not require other solvers to fail here.
@pytest.mark.xfail(strict=True, reason="Finding F2 (SSJ 1.0.0): partial-pivoting LU on H_U packed "
                   "in variable-block order shows element growth ~1e7 at T=300 (cond ~2e4); "
                   "date-interleaved ordering or iterative refinement removes it.")
def test_ikc_identity_ssj_solver_rule_long_horizon(t5):
    res = t5['ssj_consistent']
    for lab in ('G_debt_r_rule', 'rshock_r_rule'):
        assert max(res['err'][lab]) < TOL, (lab, res['err'][lab])


@pytest.mark.xfail(strict=True, reason="Finding F1 (SSJ 1.0.0): HetBlock linear impulses recompute "
                   "Jacobians despite supplied Js (own options, default one-sided h=1e-4); with "
                   "two-sided Js the r-shock IKC holds only to ~1e-4.")
def test_ikc_identity_with_mismatched_options(t5):
    res = t5['ssj_inconsistent_options']
    assert res['max_err'] < TOL, res['err']


def test_walras_law(t5):
    res = t5['ssj_consistent_refined']
    for lab, v in res['walras_asset_mkt'].items():
        assert max(v) < 1e-7, (lab, v)
