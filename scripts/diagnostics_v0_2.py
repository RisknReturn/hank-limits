"""v0.2 numerical diagnostics requested in the R4 review. Writes out/diagnostics_v0_2.json.

A. FastHetBlock distribution audit (inverse-iteration report, shift / sweep sensitivity,
   SSJ forward-iteration comparison) at a mid-risk point (T1 eps=0.1, T4 sd=0.2, T5 model)
   and near the limit (T1 eps=0.0125, SSJ forward iteration with a capped budget).
B. T4 non-mechanical residuals at every T4 point + asset-bound / resolution sensitivity.
C. T1 observation window vs computing horizon: Jacobians computed separately at T=50,100,200.
D. T1 finite-difference step sensitivity (h) at eps = 0 and eps = 0.0125.
E. T3 countercyclical incidence: aggregate income is preserved to first order only.
Usage: python3 scripts/diagnostics_v0_2.py   (about 4-6 minutes on one core)
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import json, time
import numpy as np
from hank_limits import models as M, suite, diagnostics as DG
from hank_limits.ra import ra_jacobians
from hank_limits.metrics import rel_max_err

OUT = os.path.join(ROOT, 'out', 'diagnostics_v0_2.json')
res = {}
t_all = time.time()

def save():
    json.dump(res, open(OUT, 'w'), indent=1)

# ---------------- A. distribution audit -------------------------------------------
print('A. distribution audit', flush=True)
A = {}
for eps, tols in [(0.1, (1e-10, 1e-12)), (0.0125, (1e-10,))]:
    m = M.pih_household(eps)
    t0 = time.time()
    A[f'T1_eps{eps}'] = DG.distribution_audit(m['block'], m['ss'], ssj_tols=tols,
                                              ssj_maxit=600_000, threshold=0.0)
    A[f'T1_eps{eps}']['kappa'] = m['kappa']
    print(f'  T1 eps={eps}: {time.time()-t0:.0f}s', json.dumps(A[f'T1_eps{eps}']['ssj_forward'])[:400], flush=True)
    res['A_distribution_audit'] = A; save()
out = M.ssj_aiyagari_r(0.2, 0.9)
A['T4_sd0.2_rho0.9'] = DG.distribution_audit(out['_model']['block'], out['_model']['ss'],
                                             threshold=1.0)
fh = M.fiscal_hank()
A['T5_fiscal'] = DG.distribution_audit(fh['hh'], fh['ss'], threshold=0.1)
res['A_distribution_audit'] = A; save()
print('  T4/T5 done', flush=True)

# ---------------- B. T4 residuals ---------------------------------------------------
print('B. T4 residuals', flush=True)
B = dict(points=[], bounds=[])
pts = [(e, suite.T4['rho_limit']) for e in suite.T4['eps']]
pts += [(s, suite.T4['rho_for_sd']) for s in suite.T4['sd_list']]
pts += [(suite.T4['sd_for_rho'], rh) for rh in suite.T4['rho_list']]
for sd, rho in dict.fromkeys(pts):
    o = M.ssj_aiyagari_r(sd, rho)
    d = {k: v for k, v in o['diag'].items()}
    B['points'].append(dict(sd=sd, rho=rho, r=o['r'], **d))
for a_max in (100.0, 200.0, 400.0):
    for n_a in (500, 1000, 2000):
        o = M.ssj_aiyagari_r(0.4, 0.9, a_max=a_max, n_a=n_a)
        B['bounds'].append(dict(sd=0.4, rho=0.9, a_max=a_max, n_a=n_a, r=o['r'],
                                mass_top_5pts=o['diag']['mass_top_5pts'],
                                mass_policy_off_grid=o['diag']['mass_policy_off_grid'],
                                euler_err_mean_weighted=o['diag']['euler']['euler_err_mean_weighted'],
                                capital_market_resid_rel=o['diag']['capital_market_resid_rel']))
res['B_T4_residuals'] = B; save()
print('  done', flush=True)

# ---------------- C. T1 window vs computing horizon --------------------------------
print('C. T1 horizon', flush=True)
C = {}
for eps in (0.0125, 0.1):
    Js = {T: M.ssj_rank_limit_provider(eps, T) for T in (50, 100, 200)}
    s = Js[200]['ss']
    row = dict(eps=eps)
    for T in (50, 100):
        for k in ('Y', 'r'):
            a, b = Js[T]['J']['C'][k], Js[200]['J']['C'][k][:T, :T]
            row[f'max_abs_diff_T{T}_vs_T200_window_{k}'] = float(np.abs(a - b).max())
    for W in (50, 100, 200):
        M_ra, Jr_ra = ra_jacobians(W, s['r'], s['eis'], s['C'], s['A'])
        row[f'err_M_window{W}_computedT{W}'] = rel_max_err(Js[W]['J']['C']['Y'], M_ra)
        row[f'err_Jr_window{W}_computedT{W}'] = rel_max_err(Js[W]['J']['C']['r'], Jr_ra)
    C[str(eps)] = row
res['C_T1_horizon'] = C; save()
print('  ', json.dumps(C)[:600], flush=True)

# ---------------- D. finite-difference step -----------------------------------------
print('D. h sensitivity', flush=True)
D = []
for eps in (0.0, 0.0125):
    for h in (1e-3, 1e-4, 1e-5, 1e-6):
        o = M.ssj_rank_limit_provider(eps, 200, h=h)
        s = o['ss']
        M_ra, Jr_ra = ra_jacobians(200, s['r'], s['eis'], s['C'], s['A'])
        D.append(dict(eps=eps, h=h, twosided=True, T=200,
                      err_M=rel_max_err(o['J']['C']['Y'], M_ra),
                      err_Jr=rel_max_err(o['J']['C']['r'], Jr_ra)))
res['D_T1_fd_step'] = D; save()
print('  ', D, flush=True)

# ---------------- E. T3 income preservation -----------------------------------------
prim = M.werning_primitives(cyclical=True)
e, pi, z = prim['e'], prim['pi'], prim['zeta']
E = []
for dY in (1e-3, 1e-2, 5e-2):
    for sgn in (1, -1):
        Y = 1 + sgn * dY
        agg = Y * np.sum(pi * e * Y ** z)
        E.append(dict(Y=Y, aggregate_income=float(agg), rel_dev=float(agg / Y - 1)))
res['E_T3_income_preservation'] = dict(zeta=z.tolist(), e=e.tolist(), pi=pi.tolist(), rows=E,
    note='sum_i pi_i e_i zeta_i = 0 preserves aggregate income to first order; the deviation is O(dY^2)')
res['settings'] = dict(A=dict(reference=dict(tau=1e-10, maxit=20, tol=1e-15, n_polish=5),
                              shift_taus=[1e-12, 1e-10, 1e-8, 1e-6], sweeps_variant=dict(maxit=[1, 2, 3], n_polish=0),
                              ssj_forward=dict(check_every=10, maxit=600000, tols=[1e-10, 1e-12]),
                              T1=dict(n_a=4800, a_max=250.0, n_e=3), T4=dict(n_a=1000, a_max=200.0, n_e=7),
                              T5=dict(n_a=500, a_max=60.0, n_e=7)),
                       B=dict(default_grid=dict(n_a=1000, a_max=200.0)), D=dict(T=200, twosided=True))
res['runtime_s'] = time.time() - t_all
save()
print('done', res['runtime_s'], flush=True)
