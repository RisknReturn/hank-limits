"""T3: Werning (2015) zero / vanishing liquidity.

Huggett economy, bonds in zero net supply, borrowing limit -b, b -> 0 (Werning's
"vanishing liquidity" refinement). Income y_i = Y e_i (Y/Yss)^{zeta_i}.
  zeta = 0 (income proportional to aggregate income, Werning eq. 15 / Prop. 2):
      GE output response to a real-rate path -> RANK's.
  zeta_low = 1 (low earners' income more cyclical = countercyclical risk):
      response -> the generalized zero-liquidity solution (discounting delta != 1),
      and stays away from RANK.
"""


def test_acyclical_converges_to_rank(t3):
    assert t3['err_acyc_vs_RANK'][-1] < 1e-3
    r = t3['fitted_rate_err_acyc_vs_RANK_asymptotic']
    assert 0.8 < r < 1.2, r


def test_countercyclical_converges_to_zero_liquidity_solution(t3):
    assert t3['err_cyc_vs_ZL'][-1] < 1e-3
    r = t3['fitted_rate_err_cyc_vs_ZL_asymptotic']
    assert 0.8 < r < 1.2, r


def test_countercyclical_deviates_from_rank(t3):
    # deviation does not vanish as b -> 0: it converges to the closed-form deviation
    zl = t3['zl_cyc_vs_RANK'][-1]
    assert zl > 0.25
    assert abs(t3['dev_cyc_vs_RANK'][-1] - zl) < 0.01
    assert t3['zero_liquidity_cyc']['delta'] > 1.0   # compounding: amplification


def test_asymptotic_monotone(t3):
    b, e = t3['b'], t3['err_acyc_vs_RANK']
    asym = [x for bb, x in zip(b, e) if bb <= t3['asymptotic_from']]
    assert all(x > y for x, y in zip(asym, asym[1:]))


def test_steady_state_rate_converges(t3):
    g = t3['r_gap_over_b_acyc']
    assert all(x > 0 for x in g)                 # r(b) > r* and r(b) - r* = O(b)
    assert abs(g[-1] - g[-2]) / g[-1] < 0.01
