"""T6: HANK - RANK differencing identity and the induced decomposition of dY_HA - dY_RA.

(I-M_HA)^{-1} - (I-M_RA)^{-1} = (I-M_HA)^{-1} (M_HA - M_RA) (I-M_RA)^{-1}, hence, for the
same steady state and the same policy paths (dG, dTax, dr),
dY_HA - dY_RA = (I-M_HA)^{-1} [(M_HA - M_RA)(dY_RA - dTax) + (Jr_HA - Jr_RA) dr].
"""


def test_matrix_identity(t6):
    for T, v in t6['by_T'].items():
        assert v['G']['identity_err'] < 1e-10, (T, v['G']['identity_err'])


def test_decomposition_fiscal_shock(t6):
    for T, v in t6['by_T'].items():
        assert v['G']['decomposition_err'] < 1e-8, (T, v['G']['decomposition_err'])


def test_decomposition_rate_shock(t6):
    for T, v in t6['by_T'].items():
        assert v['r']['decomposition_err'] < 1e-8, (T, v['r']['decomposition_err'])


def test_rank_benchmark_is_ricardian(t6):
    # same policy path in RANK: deficit-financed G has (PV) multiplier exactly 1, up to the
    # part of the debt path cut off at T, which decays like (1 + r - phi_B)^T = 0.92^T
    # (about 4e-6 at T = 150 and 1e-11 at T = 300)
    # v0.3: the closed-form RA assertion only; the HA assertion is test_ha_multiplier_exceeds_one
    # (until v0.2.1 both were one test; docs/MUTATION_v0_1.md, section 6)
    for T, v in t6['by_T'].items():
        tail = (1 + 0.02 - 0.1) ** int(T)
        assert abs(v['G']['cum_multiplier_RA'] - 1) < 5 * tail, (T, v['G']['cum_multiplier_RA'])


def test_ha_multiplier_exceeds_one(t6):
    # a property of the HA solution in this calibration, not a closed form (split from
    # test_rank_benchmark_is_ricardian in v0.3)
    for T, v in t6['by_T'].items():
        assert v['G']['cum_multiplier_HA'] > 1, (T, v['G']['cum_multiplier_HA'])


def test_wealth_group_split_adds_up(t6):
    for T, v in t6['by_T'].items():
        assert v['group_sum_vs_M_x'] < 1e-4                 # nonlinear simulation, h = 1e-4
        assert v['group_first_round_sum_vs_total'] < 1e-4


def test_first_round_split_independent_of_T(t6):
    a, b = (t6['by_T'][k]['groups'] for k in sorted(t6['by_T'], key=int))
    for ga, gb in zip(a, b):
        assert abs(ga['first_round'][0] - gb['first_round'][0]) < 1e-8
