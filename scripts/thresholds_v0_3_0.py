"""Thresholds of the quantities added in v0.3.0 (docs/CHECKS_v0_3_0_AMENDMENT.md, section 5).

One run of the unmodified suite: the check functions T1-T5 exactly as the pytest fixtures call
them (`hank_limits.suite`), plus T3 with the reference provider solved on the goods-market target
(`suite.run_t3(ge_target='goods_mkt')`). It writes b_max = the largest value over all settings
of every V quantity (T5: the setting with matched options and refinement; path quantities: the
largest over all shock labels and horizons), and the other market's V1 residual (diagnostic).

    python3 scripts/thresholds_v0_3_0.py --out out/thresholds_v0_3_0_run1.json
    python3 scripts/thresholds_v0_3_0.py --out out/thresholds_v0_3_0_run2.json
    python3 scripts/thresholds_v0_3_0.py --table out/thresholds_v0_3_0_run1.json out/thresholds_v0_3_0_run2.json

Rule (plan, section 4): threshold = roundup_1sig(max(100 b_max, 10 tau)), b_max over both runs.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import platform
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

FISCAL = 'ssj_consistent_refined'
NEW = {   # quantity -> (model, tau, tau class); all are linear-solve residuals or exact identities
    'ge_target_resid_reported_goods_acyc': ('T3 goods', 1e-14, 'exact up to rounding (linear solve)'),
    'ge_target_resid_reported_goods_cyc': ('T3 goods', 1e-14, 'exact up to rounding (linear solve)'),
    'ge_target_resid_direct_goods_acyc': ('T3 goods', 1e-14, 'exact up to rounding (linear solve)'),
    'ge_target_resid_direct_goods_cyc': ('T3 goods', 1e-14, 'exact up to rounding (linear solve)'),
    'exog_G': ('T5 paths', 1e-14, 'exact up to rounding'),
    'exog_r': ('T5 paths', 1e-14, 'exact up to rounding'),
    'exog_rshock': ('T5 paths', 1e-14, 'exact up to rounding'),
}


def roundup_1sig(x):
    if x == 0:
        return 0.0
    e = math.floor(math.log10(x))
    m = round(x / 10 ** e, 9)
    if m >= 10:
        m, e = m / 10, e + 1
    return float(f'{math.ceil(m)}e{e}')


def _bmax(vals):
    xs = [x for x in vals if x is not None]
    return max(xs) if xs else None


def _quantities(v):
    """{quantity: b_max} over the settings of one v0.3 record (paths: over labels and T)."""
    out = {k: _bmax(v[k]) for k in v.get('at', {})}
    for p in v.get('paths', {}).values():
        for k in p.get('at', {}):
            b = _bmax(p[k])
            if b is not None:
                out[k] = b if out.get(k) is None else max(out[k], b)
    return out


def run():
    from hank_limits import suite
    t0 = time.time()
    res, times = {}, {}
    for name, fn in [('t1', suite.run_t1), ('t2', suite.run_t2), ('t3', suite.run_t3),
                     ('t4', suite.run_t4), ('t5', lambda: suite.run_t5()[FISCAL]),
                     ('t3_goods', lambda: suite.run_t3(ge_target='goods_mkt'))]:
        t = time.time()
        r = fn()
        v = r['v0_3']
        times[name] = time.time() - t
        res[name] = dict(b_max=_quantities(v), errors=v.get('errors', []), skipped=v['skipped'],
                         ge_target=v.get('ge_target'),
                         diagnostics_max={k: _bmax(x) for k, x in v.get('diagnostics', {}).items()
                                          if k != 'at'})
    from importlib import metadata
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()
    dirty = subprocess.run(['git', 'status', '--porcelain'], cwd=ROOT, capture_output=True,
                           text=True).stdout.splitlines()
    env = dict(utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec='seconds'),
               python=sys.version.split()[0], platform=platform.platform(), git_head=head,
               worktree_status=dirty, runtime_s=time.time() - t0, runtime_by_check_s=times,
               packages={p: metadata.version(p) for p in
                         ('sequence-jacobian', 'numpy', 'scipy', 'numba', 'llvmlite')})
    return dict(environment=env, results=res)


def table(paths):
    runs = [json.load(open(p))['results'] for p in paths]
    rows = []
    for q, (model, tau, cls) in NEW.items():
        key = 't3_goods' if model == 'T3 goods' else 't5'
        b = [r[key]['b_max'].get(q) for r in runs]
        bm = max(x for x in b if x is not None)
        th = roundup_1sig(max(100 * bm, 10 * tau))
        rows.append((model, q, b, tau, cls, th))
        print(f"| {model} | `{q}` | {', '.join(f'{x:.3g}' for x in b)} | {tau:g} ({cls}) | {th:g} |")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', help='write one run to this JSON file')
    ap.add_argument('--table', nargs='+', help='threshold table from run files')
    args = ap.parse_args()
    if args.table:
        table(args.table)
        return
    r = run()
    with open(args.out, 'w') as f:
        json.dump(r, f, indent=1)
    print(json.dumps({k: v['b_max'] for k, v in r['results'].items()
                      if k in ('t3_goods', 't5')}, indent=1))


if __name__ == '__main__':
    main()
