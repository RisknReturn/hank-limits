"""Run all analytic-limit checks plus discretisation studies; write JSON and PNG plots.

Usage:  python3 scripts/convergence_report.py            (from the project root)
Output: out/convergence_v0_2.json, out/*.png, out/tables_v0_2/suite_results_long_v0_2.csv
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import json  # noqa: E402
import platform  # noqa: E402
import time  # noqa: E402
from datetime import datetime, timezone  # noqa: E402
from importlib.metadata import version  # noqa: E402

import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import sequence_jacobian as sj  # noqa: E402

from hank_limits import suite, checks as C, models as M  # noqa: E402

OUT = os.path.join(ROOT, 'out')
os.makedirs(OUT, exist_ok=True)

# --- chart styling (reference palette slots 1-3, light surface) -------------------------
SURF, INK, INK2, GRID, REF = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df', '#8a8984'
S1, S2, S3 = '#2a78d6', '#eb6834', '#1baf7a'
plt.rcParams.update({
    'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
    'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2,
    'ytick.color': INK2, 'text.color': INK, 'axes.grid': True, 'grid.color': GRID,
    'grid.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False,
    'font.size': 10, 'axes.titlesize': 11, 'legend.frameon': False, 'lines.linewidth': 2,
    'lines.markersize': 6})


def ref_line(ax, x, y_anchor, order, label):
    x = np.asarray(x, float)
    ax.plot(x, y_anchor * (x / x[-1]) ** order, ls='--', lw=1.2, color=REF, label=label)


def finish(fig, ax, name, xlabel, ylabel, title, legend_kw=None, invert=False):
    if invert:                      # the limit lies to the right
        ax.invert_xaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc='left')
    ax.legend(fontsize=9, **(legend_kw or dict(loc='best')))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=150)
    plt.close(fig)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items() if not str(k).startswith('_')}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    return o


def timed(f, *a, **k):
    t0 = time.time()
    r = f(*a, **k)
    return r, time.time() - t0


def main():
    t_start = time.time()
    report = dict(meta=dict(
        version='0.2.0', generated_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'),
        python=platform.python_version(), sequence_jacobian=version('sequence-jacobian'),
        numpy=np.__version__, scipy=version('scipy'), numba=version('numba'),
        settings=dict(T1=suite.T1, T2=suite.T2, T3=suite.T3, T4=suite.T4, T5=suite.T5,
                      T6=suite.T6),
        error_metric='relative max-abs: max|x - ref| / max|ref| over all entries'))

    print('T1 ...', flush=True)
    report['T1'] = suite.run_t1()
    print('T2 ...', flush=True)
    report['T2'] = suite.run_t2()
    print('T3 ...', flush=True)
    report['T3'] = suite.run_t3()
    print('T4 ...', flush=True)
    report['T4'] = suite.run_t4()
    print('T5 ...', flush=True)
    report['T5'] = suite.run_t5()
    print('T6 ...', flush=True)
    report['T6'] = suite.run_t6()

    # ---------------- discretisation studies -----------------------------------------
    print('grid studies ...', flush=True)
    studies = {}
    # T1: asset-grid refinement at the two smallest eps
    g = dict(n_a=[1200, 2400, 4800, 9600], eps=[0.025, 0.0125], err_M={}, err_Jr={})
    for na in g['n_a']:
        r = C.check_rank_limit(M.ssj_rank_limit_provider, g['eps'], suite.T1['T'], n_a=na)
        g['err_M'][str(na)] = r['err_M']
        g['err_Jr'][str(na)] = r['err_Jr']
    studies['T1_grid'] = g
    # T2: asset-grid refinement of the saver block, GE IRF error (chi = 1.5)
    dr2 = suite.real_rate_path(suite.T2['T'], suite.T2['dr0'], suite.T2['rho_r'])
    g = dict(n_a=[2400, 4800, 9600], eps=[0.025, 0.0125], err_IRF={}, rate={})
    for na in g['n_a']:
        r = C.check_tank_limit(M.ssj_tank_provider, g['eps'], suite.T2['T'], dr2,
                               chis=(1.5,), n_a=na)
        g['err_IRF'][str(na)] = r['err_IRF_vs_TANK']['1.5']
        g['rate'][str(na)] = r['rates_err_IRF_chi1.5'][0]
    studies['T2_grid'] = g
    # T3: grid refinement at b = 1e-3 and 1e-4
    dr3 = suite.real_rate_path(suite.T3['T'], suite.T3['dr0'], suite.T3['rho_r'])
    g = dict(n_a=[200, 400, 800], b=[1e-3, 1e-4], err_acyc={}, err_cyc={})
    for na in g['n_a']:
        r = C.check_werning(M.ssj_werning_provider, g['b'], suite.T3['T'], dr3,
                            asymptotic_from=1.0, n_a=na)
        g['err_acyc'][str(na)] = r['err_acyc_vs_RANK']
        g['err_cyc'][str(na)] = r['err_cyc_vs_ZL']
    studies['T3_grid'] = g
    # T4: grid refinement at the smallest eps
    g = dict(n_a=[500, 1000, 2000], eps=0.025, gap={})
    for na in g['n_a']:
        out = M.ssj_aiyagari_r(0.025, 0.9, n_a=na)
        g['gap'][str(na)] = out['rho_disc'] - out['r']
    studies['T4_grid'] = g
    # T5: LU element growth on SSJ's stacked H_U (finding F2)
    m = M.fiscal_hank()
    g = dict(T=[100, 200, 300], growth_rule=[], cond_rule=[], growth_fixed=[], cond_fixed=[])
    for T in g['T']:
        J = M.household_jacobians(m['hh'], m['ss'], T, ['Y', 'Tax', 'r'], ['C', 'A'])
        Jd = {'hh': sj.JacobianDict(J, name='hh', T=T)}
        for kind, model, ss in (('rule', m['model_rule'], m['ss_rule']),
                                ('fixed', m['model_fixed'], m['ss'])):
            H = model.jacobian(ss, ['Y', 'B'], ['goods_mkt', 'fiscal_res'], T=T, Js=Jd,
                               options=M._het_opts(['hh'])).pack(T)
            _, _, U = sla.lu(H)
            g[f'growth_{kind}'].append(float(np.abs(U).max() / np.abs(H).max()))
            g[f'cond_{kind}'].append(float(np.linalg.cond(H)))
    studies['T5_lu_growth'] = g
    report['discretisation_studies'] = studies
    report['meta']['total_runtime_s'] = time.time() - t_start

    with open(os.path.join(OUT, 'convergence_v0_2.json'), 'w') as f:
        json.dump(jsonable(report), f, indent=1)

    make_plots(report)
    export_tables(report)
    print(f"done in {report['meta']['total_runtime_s']:.0f}s -> {OUT}")


def export_tables(rep):
    """Machine-readable long-format tables (out/tables_v0_2/*.csv), one row per number, with
    the limit parameter, state grid, finite-difference step, observation window / computing
    horizon and linear solver kept in separate columns."""
    import csv
    d = os.path.join(OUT, 'tables_v0_2')
    os.makedirs(d, exist_ok=True)
    cols = ['test', 'limit_param', 'limit_value', 'n_a', 'a_max', 'fd_h', 'fd_twosided',
            'window', 'horizon_T', 'solver', 'metric', 'value']
    S = suite
    rows = []

    def add(**k):
        rows.append({c: k.get(c, '') for c in cols})
    t1 = rep['T1']
    for i, e in enumerate(t1['eps']):
        for met in ('err_M', 'err_Jr', 'kappa', 'mass_constrained'):
            add(test='T1', limit_param='eps', limit_value=e, n_a=4800, a_max=250.0, fd_h=1e-4,
                fd_twosided=True, window=S.T1['T'], horizon_T=S.T1['T'], solver='none (PE)',
                metric=met, value=t1[met][i])
        for w in S.T1['T_sub']:
            for met in ('err_M', 'err_Jr'):
                add(test='T1', limit_param='eps', limit_value=e, n_a=4800, a_max=250.0, fd_h=1e-4,
                    fd_twosided=True, window=w, horizon_T=S.T1['T'], solver='none (PE)',
                    metric=met, value=t1[f'{met}_by_window'][str(w)][i])
    g = rep['discretisation_studies']['T1_grid']
    for na in g['n_a']:
        for i, e in enumerate(g['eps']):
            for met in ('err_M', 'err_Jr'):
                add(test='T1_grid', limit_param='eps', limit_value=e, n_a=na, a_max=250.0, fd_h=1e-4,
                    fd_twosided=True, window=S.T1['T'], horizon_T=S.T1['T'], solver='none (PE)',
                    metric=met, value=g[met][str(na)][i])
    t2 = rep['T2']
    for i, e in enumerate(t2['eps']):
        for met in ('err_MH_vs_I', 'err_MS_vs_MRA', 'err_Magg_vs_TANK'):
            add(test='T2', limit_param='eps', limit_value=e, n_a=4800, a_max=250.0, fd_h=1e-4,
                fd_twosided=True, window=S.T2['T'], horizon_T=S.T2['T'], solver='none (PE)',
                metric=met, value=t2[met][i])
        for c in t2['err_IRF_vs_TANK']:
            add(test='T2', limit_param='eps', limit_value=e, n_a=4800, a_max=250.0, fd_h=1e-4,
                fd_twosided=True, window=S.T2['T'], horizon_T=S.T2['T'], solver='SSJ solve_impulse_linear',
                metric=f'err_IRF_vs_TANK_chi{c}', value=t2['err_IRF_vs_TANK'][c][i])
            add(test='T2', limit_param='eps', limit_value=e, n_a=4800, a_max=250.0, fd_h=1e-4,
                fd_twosided=True, window=S.T2['T'], horizon_T=S.T2['T'], solver='none (PE)',
                metric=f'err_Magg_lamchi_vs_TANK_chi{c}', value=t2['err_Magg_chi_vs_TANK'][c][i])
    g = rep['discretisation_studies']['T2_grid']
    for na in g['n_a']:
        for i, e in enumerate(g['eps']):
            add(test='T2_grid', limit_param='eps', limit_value=e, n_a=na, a_max=250.0, fd_h=1e-4,
                fd_twosided=True, window=S.T2['T'], horizon_T=S.T2['T'], solver='SSJ solve_impulse_linear',
                metric='err_IRF_vs_TANK_chi1.5', value=g['err_IRF'][str(na)][i])
    t3 = rep['T3']
    for i, b in enumerate(t3['b']):
        for met in ('err_acyc_vs_RANK', 'err_cyc_vs_ZL', 'dev_cyc_vs_RANK', 'r_gap_over_b_acyc'):
            add(test='T3', limit_param='b', limit_value=b, n_a=400, a_max='-b + 20 b', fd_h='1e-3 b',
                fd_twosided=True, window=S.T3['T'], horizon_T=S.T3['T'], solver='SSJ solve_impulse_linear',
                metric=met, value=t3[met][i])
    g = rep['discretisation_studies']['T3_grid']
    for na in g['n_a']:
        for i, b in enumerate(g['b']):
            for met in ('err_acyc', 'err_cyc'):
                add(test='T3_grid', limit_param='b', limit_value=b, n_a=na, a_max='-b + 20 b', fd_h='1e-3 b',
                    fd_twosided=True, window=S.T3['T'], horizon_T=S.T3['T'],
                    solver='SSJ solve_impulse_linear', metric=met, value=g[met][str(na)][i])
    t4 = rep['T4']
    for i, e in enumerate(t4['eps']):
        for met in ('r', 'gap', 'capital_market_resid_rel', 'euler_err_mean_weighted',
                    'euler_err_max_with_mass', 'kkt_violation_max_with_mass', 'mass_top_5pts'):
            add(test='T4', limit_param='sd_log_income', limit_value=e, n_a=1000, a_max=200.0,
                solver='brentq on A(r) - K(r)', metric=met, value=t4[met][i])
    g = rep['discretisation_studies']['T4_grid']
    for na in g['n_a']:
        add(test='T4_grid', limit_param='sd_log_income', limit_value=g['eps'], n_a=na, a_max=200.0,
            solver='brentq on A(r) - K(r)', metric='gap', value=g['gap'][str(na)])
    for label, v in rep['T5'].items():
        for shock, errs in v['err'].items():
            for T, e in zip(v['T'], errs):
                add(test='T5', limit_param='identity', n_a=500, a_max=60.0, fd_h=1e-4,
                    fd_twosided=True, window=T, horizon_T=T, solver=label, metric=f'ikc_err_{shock}', value=e)
    g = rep['discretisation_studies']['T5_lu_growth']
    for i, T in enumerate(g['T']):
        for met in ('growth_rule', 'cond_rule', 'growth_fixed', 'cond_fixed'):
            add(test='T5_F2', limit_param='identity', n_a=500, a_max=60.0, fd_h=1e-4, fd_twosided=True,
                window=T, horizon_T=T, solver='scipy.linalg.lu (partial pivoting)', metric=met, value=g[met][i])
    for T, v in rep['T6']['by_T'].items():
        for met in ('identity_err', 'decomposition_err', 'cum_multiplier_HA', 'cum_multiplier_RA'):
            add(test='T6', limit_param='identity', n_a=500, a_max=60.0, fd_h=1e-4, fd_twosided=True,
                window=T, horizon_T=T, solver='refined LU', metric=f'G_{met}', value=v['G'][met])
    with open(os.path.join(d, 'suite_results_long_v0_2.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f'tables: {len(rows)} rows -> {d}')


def make_plots(rep):
    # T1
    r = rep['T1']
    eps = r['eps'][1:]
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.loglog(eps, r['err_M'][1:], marker='o', color=S1, label='dC/dY (iMPC matrix)')
    ax.loglog(eps, r['err_Jr'][1:], marker='s', ls='-.', color=S2, label='dC/dr')
    ref_line(ax, eps, r['err_Jr'][-1], 2, 'slope 2 reference')
    finish(fig, ax, 't1_rank_limit.png', 'income risk eps (sd of log income), decreasing ->',
           'relative max error vs RA closed form',
           f'T1 RANK limit: HA Jacobians -> RA (T={r["T"]})', invert=True)

    # T2
    r = rep['T2']
    eps = r['eps']
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.loglog(eps, r['err_MS_vs_MRA'], marker='o', color=S1, label='saver iMPC vs M_RA')
    ax.loglog(eps, r['err_Magg_vs_TANK'], marker='s', ls='-.', color=S2,
              label='aggregate iMPC vs lam I + (1-lam) M_RA')
    ax.loglog(eps, r['err_IRF_vs_TANK']['1.5'], marker='^', ls=':', color=S3,
              label='GE output IRF vs Bilbiie (chi = 1.5)')
    ref_line(ax, eps, r['err_MS_vs_MRA'][-1], 2, 'slope 2 reference')
    finish(fig, ax, 't2_tank_limit.png', 'income risk eps (both types), decreasing ->',
           'relative max error', 'T2 TANK limit (H type exact for eps <= 0.05)', invert=True)

    # T3
    r = rep['T3']
    b = r['b']
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.loglog(b, r['err_acyc_vs_RANK'], marker='o', color=S1,
              label='income proportional to Y: error vs RANK')
    ax.loglog(b, r['err_cyc_vs_ZL'], marker='s', ls='-.', color=S2,
              label='countercyclical: error vs zero-liquidity solution')
    ax.loglog(b, r['dev_cyc_vs_RANK'], marker='^', ls=':', color=S3,
              label='countercyclical: deviation from RANK')
    ref_line(ax, [x for x in b if x <= 1e-3], r['err_acyc_vs_RANK'][-1], 1, 'slope 1 reference')
    finish(fig, ax, 't3_werning.png', 'borrowing limit b (liquidity), decreasing ->',
           'relative max error / deviation', 'T3 Werning (2015) vanishing liquidity',
           invert=True)

    # T4
    r = rep['T4']
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.loglog(r['eps'], r['gap'], marker='o', color=S1, label='1/beta - 1 - r')
    ref_line(ax, r['eps'], r['gap'][-1], 2, 'slope 2 reference')
    finish(fig, ax, 't4_aiyagari.png', 'income risk eps (sd of log income), decreasing ->',
           'gap to complete-markets rate', 'T4 Aiyagari: r -> 1/beta - 1 as risk -> 0',
           invert=True)

    # T5
    r = rep['T5']
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for (key, lab), col, mk, ls in zip(
            [('ssj_consistent_refined', 'SSJ blocks, refined solve'),
             ('ssj_consistent', 'SSJ solve_impulse_linear'),
             ('ssj_inconsistent_options', 'SSJ, Js / options mismatch')],
            (S1, S2, S3), ('o', 's', '^'), ('-', '-.', ':')):
        worst = np.max(np.array(list(r[key]['err'].values())), axis=0)
        ax.semilogy(r[key]['T'], worst, marker=mk, ls=ls, color=col, label=lab)
    ax.axhline(1e-10, color=REF, ls='--', lw=1.2, label='test tolerance 1e-10')
    finish(fig, ax, 't5_ikc.png', 'truncation horizon T', 'worst relative IKC error (4 shocks)',
           'T5 Intertemporal Keynesian cross identity',
           legend_kw=dict(loc='upper left', bbox_to_anchor=(0.0, 0.9)))

    # T6: decomposition by date (levels, not log)
    v = rep['T6']['by_T'][str(suite.T6['T'])]['G']
    d = np.array(v['dates'])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.bar(d - 0.2, v['first_round'], width=0.38, color=S1,
           label='first round (M_HA - M_RA)(dY_RA - dTax)')
    ax.bar(d + 0.2, v['ge_propagation'], width=0.38, color=S2,
           label='GE propagation via (I - M_HA)^-1')
    ax.plot(d, v['dY_gap'], marker='o', color=INK2, lw=1.5, label='dY_HA - dY_RA (total)')
    ax.axhline(0, color=GRID, lw=1)
    finish(fig, ax, 't6_decomposition.png', 'year after deficit-financed G shock',
           'output gap HANK - RANK (share of SS output)',
           'T6 HANK - RANK decomposition by date')


if __name__ == '__main__':
    if '--plots-only' in sys.argv:          # redraw from an existing JSON
        with open(os.path.join(OUT, 'convergence_v0_2.json')) as f:
            make_plots(json.load(f))
    else:
        main()
