"""Aiyagari (1994) Table II: two cells, with the input definitions made explicit.

Reference values (net return to capital, %, and saving rate, %) as transcribed by the R4
reviewer from the original Table II (QJE 109(3), printed p. 678):
    sigma=0.2, rho=0.6, mu=3:  3.8767 / 24.25
    sigma=0.4, rho=0.9, mu=5: -0.3456 / 37.63
We could not open the scanned original ourselves (no machine-readable text); the R4 follow-up
reviewer checked Tables I and II visually. s = delta*alpha/(r + delta) is a function of r, so a matching
saving rate is not a second, independent check.
Status (v0.2.1): cell 1 is close (Tauchen 3 sd: 3.8782 vs 3.8767). Cell 2 is NOT reproduced (-0.086 vs -0.3456);
the Tauchen 3 sd chain's level CV and level autocorrelation already match Table I to two decimals, so the gap is
not explained by the chain's risk; the cause is unresolved (candidates: Aiyagari's simulation-based solution,
the exact transition matrix, other implementation details - none verified).

What differs from T4 of the suite: risk aversion mu = 3 and 5 (T4 uses EIS 0.5, i.e. mu = 2);
Aiyagari's sigma is described as the coefficient of variation of earnings and his inputs are
those of the AR(1) for log labour endowment, approximated by a 7-state Markov chain whose
actual moments differ from the inputs (his Table I, footnote 38, per the reviewer).
We therefore report a small grid of readings:
    chain  in {Tauchen (1986), m = 3 sd; Rouwenhorst},
    sigma  read as {sd of log endowment; coefficient of variation => sd_log = sqrt(log(1+sigma^2))},
with beta=0.96, alpha=0.36, delta=0.08, zero borrowing limit, mean endowment normalised to 1
(with a zero limit the household problem is homogeneous in the endowment level, so the
normalisation does not move r). The point is to explain the gap, not to tune to it.
Usage: python3 scripts/aiyagari_table2_v0_2.py   Writes out/aiyagari_table2_v0_2.json
"""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import numpy as np
from scipy.stats import norm
import sequence_jacobian as sj
from hank_limits import models as M

REF = {(0.2, 0.6, 3): (3.8767, 24.25), (0.4, 0.9, 5): (-0.3456, 37.63)}
al, de = 0.36, 0.08


def tauchen(rho, sd_log, n=7, m=3.0):
    se = sd_log * np.sqrt(1 - rho ** 2)
    x = np.linspace(-m * sd_log, m * sd_log, n)
    step = x[1] - x[0]
    P = np.empty((n, n))
    for i in range(n):
        z = (x - rho * x[i]) / se
        lo = norm.cdf(z - step / (2 * se))
        hi = norm.cdf(z + step / (2 * se))
        P[i] = hi - lo
        P[i, 0] = hi[0]
        P[i, -1] = 1 - lo[-1]
    return x, P


def chain(kind, rho, sd_log):
    if kind == 'tauchen':
        x, P = tauchen(rho, sd_log)
    else:
        x, _, P = sj.grids.markov_rouwenhorst(rho, sd_log, 7)
        x = np.log(x)
    w, v = np.linalg.eig(P.T)
    pi = np.real(v[:, np.argmin(np.abs(w - 1))]); pi /= pi.sum()
    e = np.exp(x); e /= pi @ e
    le = np.log(e); mle = pi @ le
    sdl = np.sqrt(pi @ (le - mle) ** 2)
    ac = (pi * (le - mle)) @ P @ (le - mle) / sdl ** 2
    cv = np.sqrt(pi @ (e - 1) ** 2)                 # mean of e is 1, so this is the level CV
    ac_lev = (pi * (e - 1)) @ P @ (e - 1) / cv ** 2  # level autocorrelation (Aiyagari Table I reports level CV and corr)
    return e, P, dict(sd_log=float(sdl), autocorr_log=float(ac), cv_level=float(cv), autocorr_level=float(ac_lev))


rows = []
for (sig, rho, mu), (r_ref, s_ref) in REF.items():
    s_check = 100 * de * al / (r_ref / 100 + de)
    for kind in ('tauchen', 'rouwenhorst'):
        for reading, sd_log in (('sd_log', sig), ('cv', float(np.sqrt(np.log(1 + sig ** 2))))):
            e, P, mom = chain(kind, rho, sd_log)
            t0 = time.time()
            out = M.ssj_aiyagari_r(None, None, chain=(e, P), eis=1 / mu, n_a=2000, a_max=400.0)
            rows.append(dict(sigma=sig, rho=rho, mu=mu, chain=kind, sigma_reading=reading,
                             sd_log_input=sd_log, chain_moments=mom, r_pct=100 * out['r'],
                             saving_rate_pct=100 * out['diag']['saving_rate'],
                             ref_r_pct=r_ref, ref_saving_pct=s_ref, ref_saving_from_r=s_check,
                             diff_r_pp=100 * out['r'] - r_ref,
                             euler_err_mean=out['diag']['euler']['euler_err_mean_weighted'],
                             mass_top_5pts=out['diag']['mass_top_5pts'], seconds=time.time() - t0))
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()
                   if k != 'chain_moments'}, mom, flush=True)
# sensitivity of the hard cell to the Tauchen grid width (explains the gap; not used to tune)
sens = []
for m in (2.0, 2.5, 3.0, 3.5):
    x, P = tauchen(0.9, 0.4, m=m)
    w_, v_ = np.linalg.eig(P.T)
    pi = np.real(v_[:, np.argmin(np.abs(w_ - 1))]); pi /= pi.sum()
    e = np.exp(x); e /= pi @ e
    le = np.log(e)
    row = dict(m=m, sd_log_chain=float(np.sqrt(pi @ (le - pi @ le) ** 2)), cv_level=float(np.sqrt(pi @ (e - 1) ** 2)))
    try:
        row['r_pct'] = 100 * M.ssj_aiyagari_r(None, None, chain=(e, P), eis=0.2, n_a=2000, a_max=400.0)['r']
    except ValueError as exc:          # backward iteration does not converge at the top of the bracket
        row['error'] = str(exc)[:80]
    sens.append(row)
    print(row, flush=True)
# asset-grid audit of the second cell with the Tauchen 3 sd chain fixed
grid_audit = []
e2, P2, _ = chain('tauchen', 0.9, 0.4)
for n_a, a_max in ((1000, 400.0), (2000, 400.0), (4000, 400.0), (2000, 800.0)):
    o = M.ssj_aiyagari_r(None, None, chain=(e2, P2), eis=0.2, n_a=n_a, a_max=a_max)
    grid_audit.append(dict(n_a=n_a, a_max=a_max, r_pct=100 * o['r'], mass_top_5pts=o['diag']['mass_top_5pts'],
                           mass_policy_off_grid=o['diag']['mass_policy_off_grid']))
    print(grid_audit[-1], flush=True)
# Aiyagari (1994) Table I (printed p. 676, as checked visually by the R4 follow-up reviewer): level CV / level corr
TABLE_I = {(0.2, 0.6): (0.21, 0.59), (0.4, 0.9): (0.49, 0.89)}
json.dump(dict(rows=rows, note=__doc__, tauchen_width_sensitivity_sigma04_rho09_mu5=sens,
               grid_audit_sigma04_rho09_mu5=grid_audit,
               table_I_level_cv_corr={str(k): v for k, v in TABLE_I.items()},
               status_cell2='not reproduced; cause unresolved (chain level CV/corr match Table I to two decimals)', settings=dict(n_a=2000, a_max=400.0, grid='a_max u^2', borrowing_limit=0.0,
                beta=0.96, alpha=0.36, delta=0.08, n_e=7, tauchen_width_sd=3.0)), open(os.path.join(ROOT, 'out', 'aiyagari_table2_v0_2.json'), 'w'), indent=1)
