"""Minimal reproduction of finding F2 (SSJ only; no code from this repository).

Observed behaviour: for a small one-asset fiscal HANK with a real-rate rule, the stacked
H_U (unknowns Y, B; targets goods_mkt, fiscal_res) that SSJ builds in variable-block order
shows large element growth under Gaussian elimination with partial pivoting (LAPACK getrf,
used by np.linalg.solve in Block.solve_impulse_linear), although cond(H_U) is moderate.
Reordering the same system date by date removes the growth. This is a solver-robustness
observation for this system, not a claim about every model or about SSJ's method.

Reports, for T = 100, 200, 300: cond_2(H_U), growth = max|U| / max|H_U| of scipy's LU,
residual and normalised backward error eta = ||H_U x - b||_inf / (||H_U||_inf ||x||_inf + ||b||_inf)
for (i) np.linalg.solve as in SSJ, (ii) LU + 2 steps of iterative refinement, (iii) pivoted
QR (gelsy), and at T=300 (iv) the date-interleaved and (v) reversed-block orderings.
Saves H_U and b at T=300 to out/F2_HU_T300_v0_2.npz.
Usage: python3 scripts/repro_F2_v0_2.py      Writes out/repro_F2_v0_2.json
"""
import os, sys, json, platform
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
from importlib.metadata import version
import numpy as np
import scipy
import scipy.linalg as la
from scipy.optimize import brentq
import sequence_jacobian as sj
from sequence_jacobian.hetblocks.hh_sim import hh as hh_base

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out')
os.makedirs(OUTDIR, exist_ok=True)


def income(Y, Tax, e_grid):
    y = (Y - Tax) * e_grid
    return y


@sj.simple
def fiscal(B, G, r, Tax_ss, B_ss, phi_B):
    Tax = Tax_ss + phi_B * (B(-1) - B_ss)
    fiscal_res = B - (1 + r) * B(-1) - G + Tax
    return Tax, fiscal_res


@sj.simple
def clearing(Y, C, G):
    goods_mkt = Y - C - G
    return goods_mkt


@sj.simple
def rate_rule(Y, rss, phi_y, rshock):
    r = rss + phi_y * (Y - 1) + rshock
    return r


hh = hh_base.add_hetinputs([income])
e, pi, Pi = sj.grids.markov_rouwenhorst(rho=0.92, sigma=0.6, N=7)
r, B, G = 0.02, 1.0, 0.2
calib = dict(Y=1.0, Tax=G + r * B, e_grid=e, Pi=Pi, a_grid=sj.grids.asset_grid(0.0, 60.0, 500),
             r=r, eis=0.5)
beta = brentq(lambda b: hh.steady_state({**calib, 'beta': b})['A'] - B, 0.80, 1 / (1 + r) - 1e-6,
              xtol=1e-14)
calib.update(beta=beta, B=B, B_ss=B, G=G, Tax_ss=G + r * B, phi_B=0.1, rss=r, phi_y=0.25, rshock=0.0)
model = sj.create_model([rate_rule, hh, fiscal, clearing], name='fiscal-HANK-r-rule')
ss = model.steady_state(calib)


def measures(A, b, x, xref):
    res = la.norm(A @ x - b, np.inf)
    return dict(residual_inf=float(res),
                backward_eta=float(res / (la.norm(A, np.inf) * la.norm(x, np.inf) + la.norm(b, np.inf))),
                rel_diff_from_qr=float(la.norm(x - xref, np.inf) / la.norm(xref, np.inf)))


def growth(A):
    lu, _ = la.lu_factor(A)
    return float(np.abs(np.triu(lu)).max() / np.abs(A).max()), lu


out = dict(rows=[], versions=dict(python=platform.python_version(), numpy=np.__version__,
                                  scipy=scipy.__version__, sequence_jacobian=version('sequence-jacobian')),
           setup=dict(beta=float(beta), n_e=7, n_a=500, rho_e=0.92, sd_e=0.6, r=r, B=B, G=G, phi_B=0.1,
                      phi_y=0.25, shock='G_t = 0.01 * 0.8^t (debt financed)', h=1e-4, twosided=True,
                      ss_goods_mkt=float(ss['goods_mkt']), ss_fiscal_res=float(ss['fiscal_res'])))
try:
    out['versions']['blas'] = np.show_config(mode='dicts')['Build Dependencies']['blas']
except Exception:
    pass
unknowns, targets = ['Y', 'B'], ['goods_mkt', 'fiscal_res']
opts = {'hh': {'twosided': True, 'h': 1e-4}}
for T in (100, 200, 300):
    Js = {'hh': hh.jacobian(ss, ['Y', 'Tax', 'r'], ['C', 'A'], T=T, twosided=True, h=1e-4)}
    A = model.jacobian(ss, unknowns, targets, T, Js, opts).pack(T)
    b = -model.impulse_linear(ss, {'G': 0.01 * 0.8 ** np.arange(T)}, targets, Js, opts).get(targets).pack()
    xq, *_ = la.lstsq(A, b, lapack_driver='gelsy')
    g, lu = growth(A)
    piv = la.lu_factor(A)
    x_np = np.linalg.solve(A, b)                      # what Block.solve_impulse_linear does
    x_ref = la.lu_solve(piv, b)
    for _ in range(2):
        x_ref = x_ref + la.lu_solve(piv, b - A @ x_ref)
    row = dict(T=T, cond2=float(np.linalg.cond(A)), lu_growth=g,
               np_linalg_solve=measures(A, b, x_np, xq), lu_refined_2=measures(A, b, x_ref, xq),
               pivoted_qr=measures(A, b, xq, xq))
    if T == 300:
        np.savez_compressed(os.path.join(OUTDIR, 'F2_HU_T300_v0_2.npz'), H_U=A, b=b)
        row['orderings'] = {}
        for name, perm in {'interleaved_by_date': np.arange(2 * T).reshape(2, T).T.ravel(),
                           'reversed_variable_blocks': np.r_[np.arange(T, 2 * T), np.arange(T)]}.items():
            Ap = A[np.ix_(perm, perm)]
            gp, _ = growth(Ap)
            xp = np.empty(2 * T)
            xp[perm] = np.linalg.solve(Ap, b[perm])
            row['orderings'][name] = dict(lu_growth=gp, **measures(A, b, xp, xq))
    out['rows'].append(row)
    print(json.dumps({k: v for k, v in row.items()}, indent=None)[:700], flush=True)
json.dump(out, open(os.path.join(OUTDIR, 'repro_F2_v0_2.json'), 'w'), indent=1)
