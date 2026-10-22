"""Compare a re-run of designs (ii) and (iii) at the fix commit of the pre-merge check with the
re-run of commit C (docs/CHECKS_v0_3_0_AMENDMENT.md, section 11.6; the rule for "same" is fixed
there before the run).

    python3 scripts/mutation_v0_3_0.py --v03-ref <fix> --design ii iii --out <scratch>/rerun_F
    python3 scripts/rerun_v0_3_0_F_compare.py --new <scratch>/rerun_F.json

Compared run by run: every test outcome (tests present in only one run are listed), every group
cell (detected, set of changed tests), every summarised number (v0.1 part and v0.3 quantities).
Differences in summarised numbers are classified as
  goods_market_path   T5 path residual max|dY - dC - dG| / s (run-to-run rounding, sections 5, 10.4)
  two_valued          t5.ssj_inconsistent_options.walras_max.rshock_r_rule (plan, section 8)
  other               anything else.
"Same" = no outcome or cell difference and no 'other' or missing number. The output JSON holds
every difference with both values, the execution record of the new run and the SHA-256 of its
raw files.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, 'out', 'mutation_v0_3_0_rerun.json')
OUT = os.path.join(ROOT, 'out', 'rerun_v0_3_0_F_compare.json')
GOODS = re.compile(r'^t5\.v0_3\.[a-z_]+\.paths\.[A-Za-z_]+\.goods_market\.\d+$')
TWO_VALUED = 't5.ssj_inconsistent_options.walras_max.rshock_r_rule'


def _flat(x, p=''):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from _flat(v, f'{p}{k}.')
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from _flat(v, f'{p}{i}.')
    else:
        yield p.rstrip('.'), x


def _kind(key):
    if GOODS.match(key):
        return 'goods_market_path'
    if key == TWO_VALUED:
        return 'two_valued'
    return 'other'


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def compare(ref, new):
    out = dict(outcomes=[], tests_only_in_new={}, tests_only_in_ref={}, cells=[], numbers=[],
               missing_numbers=[])
    for rid, rr in ref['runs'].items():
        nr = new['runs'].get(rid)
        if nr is None:
            out['outcomes'].append(dict(run=rid, problem='run missing in new'))
            continue
        oa = {k: v['outcome'] for k, v in rr['outcomes'].items()}
        ob = {k: v['outcome'] for k, v in nr['outcomes'].items()}
        for k in oa.keys() & ob.keys():
            if oa[k] != ob[k]:
                out['outcomes'].append(dict(run=rid, test=k, ref=oa[k], new=ob[k]))
        for k in ob.keys() - oa.keys():
            out['tests_only_in_new'].setdefault(k, {}).setdefault(ob[k], []).append(rid)
        for k in oa.keys() - ob.keys():
            out['tests_only_in_ref'].setdefault(k, []).append(rid)
        fa = dict(_flat({fx: r.get('summary', r.get('error')) for fx, r in rr['fixtures'].items()}))
        fb = dict(_flat({fx: r.get('summary', r.get('error')) for fx, r in nr['fixtures'].items()}))
        for k in fa.keys() | fb.keys():
            if k not in fa or k not in fb:
                out['missing_numbers'].append(dict(run=rid, key=k, ref=fa.get(k), new=fb.get(k)))
            elif fa[k] != fb[k]:
                out['numbers'].append(dict(run=rid, key=k, kind=_kind(k), ref=fa[k], new=fb[k]))
    for d, res in ref['results']['designs'].items():
        new_tab = new['results']['designs'][d]['table']
        for i, t in res['table'].items():
            for g, c in t['groups'].items():
                n = new_tab[i]['groups'][g]
                ca = (c['detected'], sorted(x['test'] for x in c['changed']))
                cb = (n['detected'], sorted(x['test'] for x in n['changed']))
                if ca != cb:
                    out['cells'].append(dict(design=d, id=i, group=g, ref=ca, new=cb))
    kinds = {}
    for x in out['numbers']:
        kinds[x['kind']] = kinds.get(x['kind'], 0) + 1
    out['number_differences_by_kind'] = kinds
    out['runs_with_number_differences'] = sorted({x['run'] for x in out['numbers']})
    out['same'] = not (out['outcomes'] or out['cells'] or out['missing_numbers'] or
                       kinds.get('other') or out['tests_only_in_ref'])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--new', required=True, help='runner JSON of the new re-run')
    ap.add_argument('--ref', default=REF)
    ap.add_argument('--out', default=OUT)
    args = ap.parse_args()
    ref, new = json.load(open(args.ref)), json.load(open(args.new))
    res = compare(ref, new)
    csv = os.path.splitext(args.new)[0] + '.csv'
    res['execution_record'] = dict(
        environment=new['environment'], suites=new['suites'], settings=new['settings'],
        reanchored=sorted(new['reanchored']),
        pytest_summaries={k: r.get('pytest_summary') for k, r in new['runs'].items()},
        wall_s={k: r.get('wall_s') for k, r in new['runs'].items()},
        baselines=new['results']['baselines'],
        suite_results={d: r['suite'] for d, r in new['results']['designs'].items()},
        controls={d: r['controls'] for d, r in new['results']['designs'].items()},
        raw_files={os.path.basename(p): dict(sha256=_sha256(p), bytes=os.path.getsize(p))
                   for p in (args.new, csv) if os.path.exists(p)},
        reference=dict(file=os.path.relpath(args.ref, ROOT), sha256=_sha256(args.ref)))
    with open(args.out, 'w') as f:
        json.dump(res, f, indent=1)
    print(f"same: {res['same']}; outcome differences {len(res['outcomes'])}; cell differences "
          f"{len(res['cells'])}; numbers by kind {res['number_differences_by_kind']} in "
          f"{len(res['runs_with_number_differences'])} runs; missing {len(res['missing_numbers'])}")
    print('tests only in new:', {k: {o: len(r) for o, r in v.items()} for k, v in res['tests_only_in_new'].items()})
    print('tests only in ref:', res['tests_only_in_ref'])
    for x in res['outcomes'][:20] + res['cells'][:20]:
        print(x)
    others = [x for x in res['numbers'] if x['kind'] != 'goods_market_path']
    for x in others[:20]:
        print(x)


if __name__ == '__main__':
    main()
