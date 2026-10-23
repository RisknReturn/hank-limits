"""Compare the v0.3.0 re-run (out/mutation_v0_3_0_rerun.json) with the v0.2 runs of the same
designs (out/mutation_v0_2.json), and compute the margins of the V catches
(docs/CHECKS_v0_3_0_AMENDMENT.md, section 10; docs/MUTATION_v0_2.md, correction of 2026-10-22).

    python3 scripts/rerun_v0_3_0_compare.py            # report (markdown fragments)
    python3 scripts/rerun_v0_3_0_compare.py --json F   # also write the numbers to F

Margin of a failing V test = the largest value of one of its quantities over its threshold
(over all settings; T5 paths over all shock labels and horizons). A mutant's margin in group V
is reported both as its largest test margin and its smallest one. Tests that failed without a
measured value (fixture error) have no margin.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'tests'))
sys.path.insert(0, ROOT)

import test_v0_3_checks as V        # noqa: E402  (thresholds and the V test registry)

V02 = os.path.join(ROOT, 'out', 'mutation_v0_2.json')
RERUN = os.path.join(ROOT, 'out', 'mutation_v0_3_0_rerun.json')
GROUPS = ('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'D', 'V')


def _num(x):
    if isinstance(x, str):
        try:
            return float(x)
        except ValueError:
            return math.nan
    return x


def _values(summ, fx, key, table):
    """All recorded values of one V quantity in a run's fixture summaries."""
    s = summ.get(fx, {}).get('summary')
    if s is None:
        return None
    v = s.get('v0_3', {})
    if fx == 't5':
        v = v.get(V.FISCAL, {})
    if table == 't5_paths':
        out = [x for p in v.get('paths', {}).values() for x in p.get(key, [])]
        return out if any(key in p for p in v.get('paths', {}).values()) else None
    return v.get(key)


def test_margin(fixtures, name):
    """(margin, quantity) of a V test in one run: largest value / threshold over its quantities."""
    base = name.split('::')[-1]
    if base not in V.V_TESTS:
        return None
    fx, keys, table = V.V_TESTS[base]
    best = None
    for k in keys:
        vals = _values(fixtures, fx, k, table)
        if not vals:
            continue
        xs = [_num(x) for x in vals if x is not None]
        if not xs:
            continue
        m = max(xs) / V.THRESH[table][k]
        if best is None or m > best[0]:
            best = (m, k)
    return best


def margins(data, design):
    """{mutant: {test: (margin, quantity, outcome)}} for the V tests that changed outcome."""
    res = data['results']['designs'][design]['table']
    out = {}
    for i, t in res.items():
        g = t['groups'].get('V')
        if not g:
            continue
        run = data['runs'][f'v0.3/{i}']['fixtures']
        rows = {}
        for c in g['changed']:
            m = test_margin(run, c['test'])
            rows[c['test']] = (m[0] if m else None, m[1] if m else None, c['outcome'])
        if rows:
            out[i] = rows
    return out


def margin_summary(mg):
    per_mutant, per_test = {}, []
    for i, rows in mg.items():
        ms = [(m, t, q) for t, (m, q, o) in rows.items() if m is not None and o == 'failed']
        if not ms:
            continue
        per_mutant[i] = dict(largest=max(ms), smallest=min(ms), n=len(ms))
        per_test += [(m, i, t, q) for m, t, q in ms]
    return per_mutant, sorted(per_test)


def compare(v02, rr, design):
    a, b = v02['results']['designs'][design]['table'], rr['results']['designs'][design]['table']
    rows = []
    for i in b:
        for g in GROUPS:
            ga, gb = a[i]['groups'].get(g), b[i]['groups'].get(g)
            if ga is None and gb is None:
                continue
            ta = sorted(c['test'] for c in ga['changed']) if ga else []
            tb = sorted(c['test'] for c in gb['changed']) if gb else []
            if (ga and ga['detected']) != (gb and gb['detected']) or ta != tb:
                rows.append(dict(id=i, group=g, v02=ta, rerun=tb,
                                 v02_n=f"{len(ta)}/{ga['n_tests'] if ga else '-'}",
                                 rerun_n=f"{len(tb)}/{gb['n_tests'] if gb else '-'}"))
    return rows


def _flat(x, p=''):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from _flat(v, f'{p}{k}.')
    elif isinstance(x, list):
        for j, v in enumerate(x):
            yield from _flat(v, f'{p}{j}.')
    else:
        yield p.rstrip('.'), x


def value_diffs(v02, rr, run_id):
    """v0.3 quantities of the v0.2 run that differ (or are missing) in the re-run: only the value
    lists (not the coverage records that v0.3.0 added)."""
    out = []
    for fx in ('t1', 't2', 't3', 't4', 't5'):
        sa = v02['runs'][run_id]['fixtures'].get(fx, {}).get('summary', {}).get('v0_3')
        sb = rr['runs'][run_id]['fixtures'].get(fx, {}).get('summary', {}).get('v0_3')
        if sa is None:
            continue
        fa, fb = dict(_flat(sa)), dict(_flat(sb or {}))
        for k, x in fa.items():
            if k not in fb:
                out.append((f'{fx}.{k}', x, 'missing'))
            elif fb[k] != x:
                out.append((f'{fx}.{k}', x, fb[k]))
    return out


def new_quantity_max(rr, design):
    """Largest value of each new quantity over all runs of a design."""
    keys = {'exog_G': 't5_paths', 'exog_r': 't5_paths', 'exog_rshock': 't5_paths',
            **{k: 't3' for k in V.V_TESTS['test_v1_werning_ge_target_resid_direct'][1] +
               V.V_TESTS['test_v1_werning_ge_target_resid_reported'][1] if 'goods' in k}}
    ids = list(rr['results']['designs'][design]['table'])
    out = {}
    for k, table in keys.items():
        fx = 't5' if table == 't5_paths' else 't3'
        best = None
        for i in ids:
            vals = _values(rr['runs'][f'v0.3/{i}']['fixtures'], fx, k, table)
            xs = [_num(x) for x in (vals or []) if x is not None]
            if xs and (best is None or max(xs) > best[0]):
                best = (max(xs), i)
        out[k] = best
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--json', default=None)
    args = ap.parse_args()
    v02 = json.load(open(V02))
    rr = json.load(open(RERUN)) if os.path.exists(RERUN) else None
    report = {}
    for label, data in (('v0.2', v02), ('rerun', rr)):
        if data is None:
            continue
        for d in ('ii', 'iii'):
            mg = margins(data, d)
            pm, pt = margin_summary(mg)
            report[f'{label}/{d}'] = dict(per_mutant={i: dict(largest=list(x['largest']), smallest=list(x['smallest']), n=x['n'])
                                                      for i, x in pm.items()},
                                          smallest_tests=[list(x) for x in pt[:8]])
            print(f'\n### margins {label} ({d})')
            for i, x in pm.items():
                print(f"{i}: largest {x['largest'][0]:.4g} ({x['largest'][1]}, {x['largest'][2]}); "
                      f"smallest {x['smallest'][0]:.4g} ({x['smallest'][1]}); {x['n']} failing tests")
            print('smallest per test:', [(f'{m:.4g}', i, t) for m, i, t, q in pt[:5]])
            print('smallest largest-per-mutant:',
                  sorted((x['largest'][0], i) for i, x in pm.items())[:3])
    if rr is not None:
        for d in ('ii', 'iii'):
            rows = compare(v02, rr, d)
            report[f'changes/{d}'] = rows
            print(f'\n### changes v0.2 -> rerun ({d}): {len(rows)} cells')
            for r in rows:
                print(r)
        diffs = {}
        for rid in rr['runs']:
            if rid in v02['runs']:
                dd = value_diffs(v02, rr, rid)
                if dd:
                    diffs[rid] = dd
        report['value_diffs'] = {k: [list(x) for x in v] for k, v in diffs.items()}
        print('\n### v0.3 quantities that differ from v0.2 (value lists)')
        for rid, dd in diffs.items():
            print(rid, len(dd), dd[:6])
        for d in ('ii', 'iii'):
            nq = new_quantity_max(rr, d)
            report[f'new_quantities/{d}'] = nq
            print(f'\n### new quantities, largest over the runs of ({d})', nq)
        b = rr['results']['baselines']['v0.3']
        report['baselines'] = b
        print('\n### baselines', json.dumps(b, indent=1)[:2000])
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(report, f, indent=1)


if __name__ == '__main__':
    main()
