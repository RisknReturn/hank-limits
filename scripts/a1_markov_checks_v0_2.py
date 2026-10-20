"""A1 (adjoint design note) first checks on real household chains, v0.2. Partial equilibrium,
stationary targets, fixed primitives and fixed prices. This is NOT the adaptive-grid algorithm.

1. Small chain: identity Q_f - Q_h = v_f' r + (q_f' I d_h - q_h' d_h), own code (independent
   of the reviewer's script, same toy numbers).
2. Real household chain (T4 Aiyagari household, sd=0.2, rho=0.9, at its GE prices; and the
   T1 household at eps=0.1): conservative perturbation of the lottery weights,
   pi -> pi + eps * s * pi (1 - pi) (s = +-1 random), rows keep summing to one.
   Linear prediction dQ = eps v' B' d vs exact Q(eps) - Q(0); remainder order should be 2.
   Also: partial sums of sum_t (Lambda^t q - Q 1) vs the Poisson solution (how many SSJ-type
   expectation vectors would be needed), for a fast-mixing and a slow-mixing chain.
3. Coarse -> fine on the same household (nested grids u^2, fine n=2001; coarse 126/251/501,
   plus a non-nested coarse grid 300 with lottery prolongation). Targets: mean beginning-of-
   period assets, mass below a fixed threshold a < 1, mean consumption.
   Variant P: fine chain built from the coarse policy interpolated to the fine grid (isolates
   the distribution / lottery discretisation given the policy).
   Variant F: fine chain from the fine policy (full stationary error, policy error included).
   For each: exact error, identity check with the exact fine adjoint (machine precision
   expected; a code check), and the effectivity of a cheap adjoint (the coarse adjoint
   interpolated to the fine grid, linearly and by a cubic spline):
   eta = (v_h->f' r + transfer) / (exact error).
Usage: python3 scripts/a1_markov_checks_v0_2.py   Writes out/a1_markov_checks_v0_2.json
"""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as la
from sequence_jacobian.utilities.interpolate import interpolate_coord_robust
from hank_limits import models as M
from hank_limits.core import make_household, transition_matrix, stationary_distribution

OUT = os.path.join(ROOT, 'out', 'a1_markov_checks_v0_2.json')
res = {}


def save():
    json.dump(res, open(OUT, 'w'), indent=1)


# ---------------------------------------------------------------- helpers (row-stochastic Lam)
def stat_dense(P):
    A = np.eye(len(P)) - P.T
    A[-1] = 1
    b = np.zeros(len(P)); b[-1] = 1
    return la.solve(A, b)


def poisson_sparse(Lam, q, d):
    """(I - Lam) v = q - Q 1, d'v = 0; Q returned as the extra unknown (= d'q)."""
    N = Lam.shape[0]
    K = sp.bmat([[sp.identity(N) - Lam, sp.csr_matrix(np.ones((N, 1)))],
                 [sp.csr_matrix(d[None, :]), None]], format='csc')
    x = spla.spsolve(K, np.r_[q, 0.0])
    return x[:-1], x[-1]


def stationary_sparse(Lam):
    N = Lam.shape[0]
    K = (sp.identity(N) - Lam.T).tolil()
    K[N - 1, :] = np.ones(N)
    b = np.zeros(N); b[-1] = 1
    return spla.spsolve(K.tocsc(), b)


# ---------------------------------------------------------------- 1. small chain
P = np.array([[.7, .25, .05], [.15, .7, .15], [.05, .25, .7]])
Pc = np.array([[.8, .2], [.3, .7]])
I_hf = np.array([[.5, 0], [.5, .5], [0, .5]])
qc, q = np.array([.5, 1.5]) ** 2, np.arange(3.) ** 2
d, dc = stat_dense(P), stat_dense(Pc)
dt = I_hf @ dc
A = np.block([[np.eye(3) - P, np.ones((3, 1))], [d[None, :], np.zeros((1, 1))]])
v = la.solve(A, np.r_[q, 0.])[:-1]
r = P.T @ dt - dt
res['1_small_chain'] = dict(actual=float(q @ d - qc @ dc), residual_pairing=float(v @ r),
                            transfer=float(q @ dt - qc @ dc),
                            identity_err=float(abs(q @ d - qc @ dc - v @ r - (q @ dt - qc @ dc))),
                            residual_mass=float(r.sum()))
print('1', res['1_small_chain'], flush=True)
save()


# ---------------------------------------------------------------- household utilities
def solve_hh(a_grid, e, Pi, r, w, beta, eis):
    blk = make_household([M.income_wage], name='hh')
    ss = blk.steady_state(dict(w=w, e_grid=e, Pi=Pi, a_grid=a_grid, r=r, beta=beta, eis=eis),
                          **M.SS_OPTS)
    st = ss.internals['hh']
    return dict(a=np.asarray(st['a']), c=np.asarray(st['c']), Dbeg=np.asarray(st['Dbeg']),
                grid=a_grid)


def lam_from_policy(Pi, a_pol, grid):
    return transition_matrix(Pi, a_pol, grid).T.tocsr()          # row-stochastic


def targets(sol, Pi, thr):
    n_e = Pi.shape[0]
    g = sol['grid']
    qA = np.tile(g, (n_e, 1)).ravel()
    qM = np.tile((g < thr).astype(float), (n_e, 1)).ravel()
    qC = (Pi @ sol['c']).ravel()                                   # C = E over e' of c, timed on Dbeg
    # NOTE (v0.2.1): at fixed r, w with stationarity and the budget constraint, C = w E[e] + r A, so C_f - C_h =
    # r (A_f - A_h): A_beg and C are the SAME check here, not two independent targets. A2 and logA are added as
    # independent smooth targets.
    qA2 = np.tile(g ** 2, (n_e, 1)).ravel()
    qlog = np.tile(np.log1p(g), (n_e, 1)).ravel()
    return dict(A_beg=qA, mass_below_thr=qM, C=qC, A2_beg=qA2, logA_beg=qlog)


# ---------------------------------------------------------------- 2. perturbation + series
aiy = M.ssj_aiyagari_r(0.2, 0.9)
r_eq = aiy['r']
K = (0.36 / (r_eq + 0.08)) ** (1 / 0.64)
w_eq = 0.64 * K ** 0.36
e7, _, Pi7 = M.rouwenhorst(0.9, 0.2, 7)
PRIM = dict(e=e7, Pi=Pi7, r=r_eq, w=w_eq, beta=0.96, eis=0.5, thr=1.0)
res['primitives'] = dict(r=r_eq, w=w_eq, beta=0.96, eis=0.5, sd=0.2, rho=0.9, n_e=7,
                         grid='a_max u^2, a_max=200', threshold=1.0)


def perturbation_check(Pi, sol, q, seed=0):
    n_e, n_a = sol['a'].shape
    grid = sol['grid']
    i, pi = interpolate_coord_robust(grid, np.clip(sol['a'], grid[0], grid[-1]))
    pi = np.clip(pi, 0, 1)
    s = np.random.default_rng(seed).choice([-1.0, 1.0], size=pi.shape)
    bdir = s * pi * (1 - pi)

    def lam(eps):
        pe = pi + eps * bdir
        rows_lo = (np.arange(n_e)[:, None] * n_a + i).ravel()
        cols = np.arange(n_e * n_a)
        lot = sp.csr_matrix((np.r_[pe.ravel(), 1 - pe.ravel()],
                             (np.r_[rows_lo, rows_lo + 1], np.r_[cols, cols])), shape=(n_e * n_a,) * 2)
        X = sp.kron(sp.csr_matrix(Pi.T), sp.identity(n_a), format='csr')
        return (lot @ X).T.tocsr()
    L0 = lam(0.0)
    d0 = stationary_sparse(L0)
    v0, Q0 = poisson_sparse(L0, q, d0)
    B = (lam(1.0) - L0)
    lin = float(v0 @ (B.T @ d0))
    rows = []
    for eps in (0.04, 0.02, 0.01, 0.005, 0.0025):
        de = stationary_sparse(lam(eps))
        dQ = float(q @ de - q @ d0)
        rows.append(dict(eps=eps, actual=dQ, linear=eps * lin, remainder=abs(dQ - eps * lin)))
    rem = [x['remainder'] for x in rows]
    order = float(np.polyfit(np.log([x['eps'] for x in rows]), np.log(rem), 1)[0])
    return dict(rows=rows, remainder_order=order, Q0=float(Q0), dQ_dir=lin,
                row_sums_B_max=float(np.abs(B @ np.ones(B.shape[0])).max()),
                v_normalisation=float(d0 @ v0))


def series_check(Lam, q, d, T_list=(10, 100, 1000, 10000, 50000)):
    v, Q = poisson_sparse(Lam, q, d)
    v = v - d @ v
    out, x, acc, t = [], q.copy(), np.zeros_like(q), 0
    for T in T_list:
        while t < T:
            acc += x - Q
            x = Lam @ x
            t += 1
        vt = acc - d @ acc
        out.append(dict(T=T, rel_err_inf=float(np.abs(vt - v).max() / np.abs(v).max())))
    return out


t0 = time.time()
grid1000 = 200.0 * np.linspace(0, 1, 1000) ** 2
sol = solve_hh(grid1000, PRIM['e'], PRIM['Pi'], PRIM['r'], PRIM['w'], 0.96, 0.5)
tq = targets(sol, PRIM['Pi'], 1.0)
res['2_perturbation'] = {k: perturbation_check(PRIM['Pi'], sol, qv) for k, qv in tq.items()
                         if k in ('A_beg', 'mass_below_thr', 'C')}
Lam = lam_from_policy(PRIM['Pi'], sol['a'], sol['grid'])
dd = stationary_sparse(Lam)
res['2_series_T4'] = series_check(Lam, tq['A_beg'], dd)
pih = M.pih_household(0.1)
st = pih['ss'].internals['hh']
Lam1 = lam_from_policy(pih['Pi'], np.asarray(st['a']), pih['a_grid'])
d1 = stationary_sparse(Lam1)
q1 = np.tile(pih['a_grid'], (len(pih['e_grid']), 1)).ravel()
res['2_series_T1_eps0.1'] = series_check(Lam1, q1, d1, T_list=(10, 100, 1000, 10000, 100000))
print('2', {k: v['remainder_order'] for k, v in res['2_perturbation'].items()},
      res['2_series_T4'][-1], res['2_series_T1_eps0.1'][-1], f'{time.time()-t0:.0f}s', flush=True)
save()


# ---------------------------------------------------------------- 3. coarse -> fine
def prolong_mass(g_h, g_f):
    """Mass- and mean-preserving prolongation (lottery split of each coarse node)."""
    i, pi = interpolate_coord_robust(g_f, g_h)
    n_h, n_f = len(g_h), len(g_f)
    return sp.csr_matrix((np.r_[pi, 1 - pi], (np.r_[i, i + 1], np.r_[np.arange(n_h), np.arange(n_h)])),
                         shape=(n_f, n_h))


def interp_fun(g_h, g_f, vals):            # adjoint (a function) -> linear interpolation
    return np.vstack([np.interp(g_f, g_h, row) for row in vals])


def interp_cubic(g_h, g_f, vals):          # enriched reconstruction of the coarse adjoint
    from scipy.interpolate import CubicSpline
    return np.vstack([CubicSpline(g_h, row)(g_f) for row in vals])


n_e = len(PRIM['e'])
rows3 = []


def coarse_to_fine(nf, cases):
  gf = 200.0 * np.linspace(0, 1, nf) ** 2
  fine = solve_hh(gf, PRIM['e'], PRIM['Pi'], PRIM['r'], PRIM['w'], 0.96, 0.5)
  Lam_f = lam_from_policy(PRIM['Pi'], fine['a'], gf)
  d_f = stationary_sparse(Lam_f)
  q_f = targets(fine, PRIM['Pi'], 1.0)
  v_f = {k: poisson_sparse(Lam_f, qv, d_f)[0] for k, qv in q_f.items()}
  for nh, nested in cases:
      gh = 200.0 * np.linspace(0, 1, nh) ** 2
      co = solve_hh(gh, PRIM['e'], PRIM['Pi'], PRIM['r'], PRIM['w'], 0.96, 0.5)
      Lam_h = lam_from_policy(PRIM['Pi'], co['a'], gh)
      d_h = stationary_sparse(Lam_h)
      q_h = targets(co, PRIM['Pi'], 1.0)
      Ih = sp.kron(sp.identity(n_e), prolong_mass(gh, gf), format='csr')
      dt = Ih @ d_h
      # variant P: fine chain from the coarse policy interpolated to the fine grid
      a_interp = interp_fun(gh, gf, co['a'])
      c_interp = interp_fun(gh, gf, co['c'])
      Lam_P = lam_from_policy(PRIM['Pi'], a_interp, gf)
      d_P = stationary_sparse(Lam_P)
      q_P = targets(dict(c=c_interp, grid=gf, a=a_interp), PRIM['Pi'], 1.0)
      for variant, Lam_x, d_x, q_x in (('P_coarse_policy', Lam_P, d_P, q_P), ('F_fine_policy', Lam_f, d_f, q_f)):
          rr = Lam_x.T @ dt - dt
          for k in q_x:
              Qx, Qh = float(q_x[k] @ d_x), float(q_h[k] @ d_h)
              err = Qx - Qh
              vx = v_f[k] if variant.startswith('F') else poisson_sparse(Lam_x, q_x[k], d_x)[0]
              transfer = float(q_x[k] @ dt - Qh)
              v_hc = poisson_sparse(Lam_h, q_h[k], d_h)[0].reshape(n_e, nh)
              v_cheap = interp_fun(gh, gf, v_hc).ravel()
              v_cub = interp_cubic(gh, gf, v_hc).ravel()
              est_exact = float(vx @ rr) + transfer
              est_cheap = float(v_cheap @ rr) + transfer
              est_cub = float(v_cub @ rr) + transfer
              rows3.append(dict(n_f=nf, n_h=nh, nested=nested, variant=variant, target=k, Q_h=Qh, Q_x=Qx,
                                error=err, rel_error=err / Qx if Qx else None, transfer=transfer,
                                residual_pairing_exact=float(vx @ rr),
                                identity_err=abs(est_exact - err), est_cheap=est_cheap,
                                effectivity_cheap=est_cheap / err if err else None,
                                effectivity_cheap_without_transfer=float(v_cheap @ rr) / err if err else None,
                                effectivity_cubic=est_cub / err if err else None,
                                residual_L1=float(np.abs(rr).sum()), residual_mass=float(rr.sum())))
              x = rows3[-1]
              print(f"3 n_f={nf} n_h={nh:4d} {variant:16s} {k:15s} err={err:+.3e} transfer={transfer:+.2e} "
                    f"id_err={x['identity_err']:.1e} eff_cheap={x['effectivity_cheap']:.3f} "
                    f"(w/o transfer {x['effectivity_cheap_without_transfer']:.3f}) eff_cubic={x['effectivity_cubic']:.3f}", flush=True)
      if nested:
          # projection compatibility (v0.2.1): I = injection of coarse masses, E = piecewise-linear interpolation of
          # coarse function values; then I' Lam_P E = Lam_h and E' I = Id, so E' r = 0 for the stationary d_h
          Ew = sp.csr_matrix(interp_fun(gh, gf, np.eye(nh)).T)               # n_f x n_h
          E = sp.kron(sp.identity(n_e), Ew, format='csr')
          lhs = (Ih.T @ Lam_P @ E).toarray()
          rr_P = Lam_P.T @ dt - dt
          res.setdefault('3_projection_check', []).append(dict(
              n_f=nf, n_h=nh, max_abs_ItLamPE_minus_Lamh=float(np.abs(lhs - Lam_h.toarray()).max()),
              max_abs_EtI_minus_Id=float(np.abs((E.T @ Ih).toarray() - np.eye(n_e * nh)).max()),
              max_abs_Et_r=float(np.abs(E.T @ rr_P).max())))
          print('  projection', res['3_projection_check'][-1], flush=True)
      res['3_coarse_to_fine'] = dict(rows=rows3)
      res['settings'] = dict(perturbation_chain=dict(n_e=7, n_a=1000, a_max=200.0, eps=[0.04, 0.02, 0.01, 0.005, 0.0025],
                                                     direction='s*pi*(1-pi), s=+-1, seed 0'),
                             coarse_to_fine=dict(runs=[dict(n_f=2001, n_h=[126, 251, 501, 300], nested=[True, True, True, False]),
                                                       dict(n_f=4001, n_h=[251], nested=[True])], threshold=1.0,
                                                 targets=['A_beg', 'mass_below_thr', 'C (= w E[e] + r A here)', 'A2_beg', 'logA_beg']))
      save()


coarse_to_fine(2001, ((126, True), (251, True), (501, True), (300, False)))
coarse_to_fine(4001, ((251, True),))
print('done', flush=True)
