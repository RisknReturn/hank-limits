import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import pytest  # noqa: E402

from hank_limits import suite  # noqa: E402


def _summary(name, res, keys):
    lines = [f"[{name}] runtime {res.get('runtime_s', float('nan')):.1f}s"]
    for k in keys:
        v = res
        for part in (k if isinstance(k, tuple) else (k,)):
            v = v[part]
        if isinstance(v, list):
            v = ', '.join(f'{x:.3e}' if isinstance(x, float) else str(x) for x in v)
        elif isinstance(v, float):
            v = f'{v:.4g}'
        lines.append(f"  {'.'.join(k) if isinstance(k, tuple) else k}: {v}")
    return '\n'.join(lines)


@pytest.fixture(scope='session')
def t1():
    r = suite.run_t1()
    print('\n' + _summary('T1', r, ['eps', 'err_M', 'err_Jr', 'fitted_rate_err_M',
                                    'fitted_rate_err_Jr', 'mass_constrained']))
    return r


@pytest.fixture(scope='session')
def t2():
    r = suite.run_t2()
    print('\n' + _summary('T2', r, ['eps', 'err_MH_vs_I', 'err_MS_vs_MRA', 'err_Magg_vs_TANK',
                                    ('err_IRF_vs_TANK', '1.5'), 'fitted_rate_err_MS_vs_MRA',
                                    'fitted_rate_err_IRF_chi1.5']))
    return r


@pytest.fixture(scope='session')
def t3():
    r = suite.run_t3()
    print('\n' + _summary('T3', r, ['b', 'err_acyc_vs_RANK', 'err_cyc_vs_ZL',
                                    'dev_cyc_vs_RANK', 'fitted_rate_err_acyc_vs_RANK_asymptotic',
                                    'fitted_rate_err_cyc_vs_ZL_asymptotic']))
    return r


@pytest.fixture(scope='session')
def t4():
    r = suite.run_t4()
    print('\n' + _summary('T4', r, ['eps', 'gap', 'fitted_rate_gap', ('r_by_sd', 'r'),
                                    ('r_by_rho', 'r')]))
    return r


@pytest.fixture(scope='session')
def t5():
    r = suite.run_t5()
    for k, v in r.items():
        print(f"\n[T5:{k}] max err {v['max_err']:.2e}; " +
              '; '.join(f"{lab}: {', '.join(f'{e:.1e}' for e in errs)}"
                        for lab, errs in v['err'].items()))
    return r


@pytest.fixture(scope='session')
def t6():
    r = suite.run_t6()
    for T, v in r['by_T'].items():
        print(f"\n[T6:T={T}] identity err {v['G']['identity_err']:.2e}, decomposition err "
              f"{v['G']['decomposition_err']:.2e}, cum multiplier HA {v['G']['cum_multiplier_HA']:.3f}"
              f" RA {v['G']['cum_multiplier_RA']:.3f}")
    return r
