"""Model-agnostic analytic-limit checks.

Every check takes a *provider* (a function implemented by the HANK code under test; see
`hank_limits.interface`) and returns a JSON-serialisable dict with the errors along the
limit sequence and the observed convergence rates. The tests in `tests/` and the report
script call these with the SSJ providers from `hank_limits.models`.
"""
from __future__ import annotations

import copy
import time
import numpy as np

from .metrics import rel_max_err, rel_fro_err, pairwise_rates, fitted_rate
from .ra import (ra_jacobians, ra_output_fixed_real_rate, tank_output_fixed_real_rate,
                 tank_impc, zero_liquidity_output, zero_liquidity_rstar)


def _tolist(x):
    return [float(v) for v in np.asarray(x).ravel()]


def _request(x):
    """An independent deep copy of what a check asks a provider for (requested experiments, rate
    paths, chi values). Pre-merge check of v0.3.0 (docs/CHECKS_v0_3_0_AMENDMENT.md, section 11):
    each check keeps its own copy, declares the expected coverage and computes the reference
    paths and residuals from it, and gives every provider call a fresh copy. A provider that
    changes the objects it is given (drops an experiment, scales a path in place) then changes
    neither what the check compares with nor the caller's objects."""
    return copy.deepcopy(x)


# ---------------------------------------------------------------------------
# v0.3 checks (docs/CHECKS_PLAN_v0_3.md; revised for v0.3.0, docs/CHECKS_v0_3_0_AMENDMENT.md):
# GE target residual of T3 (V1), steady-state residuals (V2), income normalisation (V3) and
# specification pins (V4), computed here from primitive quantities the provider reports in the
# optional fields (`hank_limits.interface`) and from the inputs the caller passed. They are
# stored under res['v0_3'] (one list entry per setting of the check), so the v0.1 part of every
# result is unchanged.
#
# v0.3.0: coverage is declared by the check, not inferred from what the provider returned. At
# each setting the expected quantities get a placeholder None, replaced when the quantity is
# computed; res['v0_3']['at'][q] holds the setting of each value of q and
# res['v0_3']['status'][q] the coverage of q. Each quantity is computed from the fields its own
# equation needs; if one is missing, the placeholder stays and the reason is recorded in
# res['v0_3']['skipped']. `v03_verdict` turns the values into the verdict of a test: a
# non-finite, negative or too large value fails whatever else is missing, and a missing value
# without failure is "incomplete", never a pass.
# ---------------------------------------------------------------------------

def stationary_marginal(Pi):
    """Stationary distribution pi of a reported transition matrix (pi Pi = pi, sum pi = 1)."""
    Pi = np.asarray(Pi, float)
    n = Pi.shape[0]
    A = Pi.T - np.eye(n)
    A[-1] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    return np.linalg.solve(A, b)


def _pi_row_sum_err(Pi):
    return float(np.max(np.abs(Pi.sum(1) - 1)))


def _pi_negative(Pi):
    return float(max(0.0, -Pi.min()))


def _mean_e_err(pi, e):
    return float(abs(pi @ np.asarray(e, float) - 1))


def _agg_income_err(pi, y, Y_ref):
    return float(abs(pi @ np.asarray(y, float) - Y_ref) / abs(Y_ref))


def income_report(e, Pi, y=None, Y_ref=None):
    """V3: Pi row-stochastic, E_pi[e] = 1 and, if y is given, sum_e pi_e y_e = Y_ref."""
    Pi = np.asarray(Pi, float)
    pi = stationary_marginal(Pi)
    out = dict(pi_row_sum_err=_pi_row_sum_err(Pi), pi_negative=_pi_negative(Pi),
               mean_e_err=_mean_e_err(pi, e))
    if y is not None and Y_ref is not None:
        out['agg_income_err'] = _agg_income_err(pi, y, Y_ref)
    return out, pi


def household_budget_resid(C, A, r, y, pi):
    """V2: aggregate household budget at the steady state, |C - r A - sum_e pi_e y_e| / C."""
    return float(abs(C - r * A - pi @ np.asarray(y, float)) / abs(C))


def euler_errors(Pi, a_grid, a_pol, c_pol, D, r, beta, eis, bind_tol=1e-10, mass_tol=1e-12):
    """V4(b): Euler errors of a reported steady-state policy at a given (beta, eis). Same
    definitions as `hank_limits.diagnostics.euler_kkt_residuals` (numpy only): c_impl =
    (beta (1+r) E[c(e', a')^(-1/eis)])^(-eis), err = c_impl / c - 1; unit-free |err| on
    unconstrained states, max(0, -err) (KKT violation) on binding states."""
    Pi, a_grid = np.asarray(Pi, float), np.asarray(a_grid, float)
    a_pol, c_pol, D = (np.asarray(x, float) for x in (a_pol, c_pol, D))
    n_e, n_a = a_pol.shape
    up = np.empty_like(c_pol)
    for e in range(n_e):
        x = a_pol[e]
        j = np.clip(np.searchsorted(a_grid, x) - 1, 0, n_a - 2)
        w = (x - a_grid[j]) / (a_grid[j + 1] - a_grid[j])
        cn = (1 - w)[None, :] * c_pol[:, j] + w[None, :] * c_pol[:, j + 1]
        up[e] = Pi[e] @ np.maximum(cn, 1e-300) ** (-1 / eis)
    err = (beta * (1 + r) * up) ** (-eis) / c_pol - 1
    binding = a_pol <= a_grid[0] + bind_tol
    unc = ~binding & ~(a_pol > a_grid[-1])
    has_mass = D > mass_tol
    w = D[unc]
    return dict(
        euler_err_mean_weighted=float((np.abs(err[unc]) * w).sum() / max(w.sum(), 1e-300)),
        euler_err_max_with_mass=float(np.abs(err[unc & has_mass]).max()) if np.any(unc & has_mass) else 0.0,
        kkt_violation_max_with_mass=float(np.maximum(0.0, -err[binding & has_mass]).max())
        if np.any(binding & has_mass) else 0.0)


# --- records: settings, placeholders, values, skips ------------------------------------------

def _new_v03():
    return dict(skipped=[], settings=[], at={})


def _expect(v, keys, where):
    """Declare that each quantity in `keys` is expected at setting `where` (placeholder None)."""
    at = v.setdefault('at', {})
    for k in keys:
        v.setdefault(k, []).append(None)
        at.setdefault(k, []).append(where)


def _put(v, key, value, where):
    """Record a value at setting `where`: it replaces the placeholder of a quantity declared
    there, and is appended with its setting otherwise."""
    x = None if value is None else float(value)
    vals, at = v.setdefault(key, []), v.setdefault('at', {}).setdefault(key, [])
    if at and at[-1] == where and vals[-1] is None:
        vals[-1] = x
    else:
        vals.append(x)
        at.append(where)


def _guard(v, fn, *args):
    """Run one v0.3 computation; an exception is recorded (the V tests then fail) and never
    stops the v0.1 part of the check, so groups T1-T6 behave as in v0.2.1."""
    try:
        fn(v, *args)
    except Exception as exc:          # noqa: BLE001
        v.setdefault('errors', []).append(f'{fn.__name__}: {type(exc).__name__}: {exc}'[:300])


def _skip(v, what, missing):
    msg = f"{what}: provider does not report {', '.join(missing)}"
    if msg not in v['skipped']:
        v['skipped'].append(msg)


class _Fields:
    """The fields of one provider result that the V checks read: `f[name]` is the value,
    `f.missing(names)` names the fields not reported (e.g. "ss_detail['y']"). `top` are
    top-level fields; keyword arguments name a nested field and its keys
    (`ss_detail=('Y', 'C')`, `policy=('a', 'c')`)."""

    def __init__(self, out, top=(), **nested):
        self.val, self.lab = {}, {}
        for k in top:
            self.lab[k] = k
            if k in out:
                self.val[k] = out[k]
        for field, keys in nested.items():
            d = out.get(field)
            for k in keys:
                self.lab[k] = f"{field}['{k}']" if isinstance(d, dict) else field
                if isinstance(d, dict) and k in d:
                    self.val[k] = d[k]

    def __getitem__(self, k):
        return self.val[k]

    def missing(self, names):
        return list(dict.fromkeys(self.lab[n] for n in names if n not in self.val))


def _eq(v, key, where, tag, f, names, fn):
    """Compute one quantity if the provider reports every field its equation needs (R2: fields
    per equation); otherwise keep the placeholder and record the missing fields."""
    miss = f.missing(names)
    if miss:
        return _skip(v, f'{key} ({tag})', miss)
    _put(v, key, fn(), where)


def _dev(x, ref, rel=False):
    return abs(x - ref) / (abs(ref) if rel else 1.0)


def _pin(v, key, where, tag, f, name, ref, rel=False):
    _eq(v, key, where, tag, f, [name], lambda: _dev(f[name], ref, rel))


def _income(v, sfx, where, tag, f, e, Pi, y, ref_names, ref):
    """V3 for one income chain (field names e, Pi, y in f; the aggregate-income reference ref()
    needs the fields ref_names). Returns the stationary distribution of Pi, or None."""
    pi = None if f.missing([Pi]) else stationary_marginal(f[Pi])
    _eq(v, 'pi_row_sum_err' + sfx, where, tag, f, [Pi],
        lambda: _pi_row_sum_err(np.asarray(f[Pi], float)))
    _eq(v, 'pi_negative' + sfx, where, tag, f, [Pi], lambda: _pi_negative(np.asarray(f[Pi], float)))
    _eq(v, 'mean_e_err' + sfx, where, tag, f, [e, Pi], lambda: _mean_e_err(pi, f[e]))
    _eq(v, 'agg_income_err' + sfx, where, tag, f, [y, Pi, *ref_names],
        lambda: _agg_income_err(pi, f[y], ref()))
    return pi


def _coverage(at, vals):
    missing = [w for w, x in zip(at, vals) if x is None]
    nonfinite = [w for w, x in zip(at, vals) if x is not None and not np.isfinite(x)]
    n = len(vals) - len(missing)
    state = ('nonfinite' if nonfinite else 'not_reported' if n == 0 else
             'incomplete' if missing else 'complete')
    return dict(state=state, n_expected=len(vals), n_reported=n, missing=missing,
                nonfinite=nonfinite)


def _finish(v):
    """Coverage summary per quantity, res['v0_3']['status'] (T5 paths: 'paths.<label>.<q>').
    Thresholds are applied by the tests (`v03_verdict`)."""
    st = {k: _coverage(at, v[k]) for k, at in v.get('at', {}).items()}
    for lab, p in v.get('paths', {}).items():
        st.update({f'paths.{lab}.{k}': _coverage(at, p[k]) for k, at in p.get('at', {}).items()})
    v['status'] = st


# --- verdict of one V test -----------------------------------------------------------------------

def _entries(v, key, paths):
    """(setting, value) pairs of quantity `key`, or None if `v` has no record of it. Without
    coverage records (no 'at'), each list entry is one setting."""
    def pairs(d, prefix=''):
        at = d.get('at', {}).get(key) or [f'#{i}' for i in range(len(d[key]))]
        return [(prefix + str(w), x) for w, x in zip(at, d[key])]
    if not paths:
        return pairs(v) if key in v else None
    found = [(lab, p) for lab, p in (v.get('paths') or {}).items() if key in p]
    return [e for lab, p in found for e in pairs(p, f'{lab}, ')] if found else None


def _value_ok(x, threshold):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return False
    return bool(np.isfinite(x) and 0.0 <= x < threshold)


def v03_verdict(v, keys, thresholds, paths=False):
    """Verdict of one V test on the quantities `keys` of a v0.3 record `v` (res['v0_3']), with
    thresholds {quantity: threshold}; `paths=True` for T5 path quantities (all shock labels).

    'failed'     a v0.3 computation raised, or a recorded value is non-finite, negative (every V
                 quantity is an absolute value) or not below its threshold, whatever else is
                 missing;
    'skipped'    otherwise, if no quantity of the test has a value: the provider does not report
                 what the test needs (reasons in 'skipped');
    'incomplete' otherwise, if an expected value is missing (never a pass);
    'passed'     otherwise.
    A quantity the check declared at no setting does not apply ('n/a'); in a record without
    coverage records (no 'at'), every listed quantity is expected."""
    covered = 'at' in v
    per = {}
    for k in keys:
        e = _entries(v, k, paths)
        if e is None and covered:
            per[k] = dict(status='n/a')
            continue
        e = e or []
        th = thresholds[k]
        bad = [(w, x) for w, x in e if x is not None and not _value_ok(x, th)]
        missing = [w for w, x in e if x is None]
        n = len(e) - len(missing)
        status = 'failed' if bad else 'skipped' if n == 0 else 'incomplete' if missing else 'passed'
        per[k] = dict(status=status, threshold=th, failed=bad, missing=missing,
                      worst=max((float(x) for _, x in e if _value_ok(x, np.inf)), default=None))
    errors = list(v.get('errors') or [])
    st = [p['status'] for p in per.values() if p['status'] != 'n/a']
    if errors or 'failed' in st:
        status = 'failed'
    elif all(s == 'skipped' for s in st):
        status = 'skipped'
    elif any(s in ('skipped', 'incomplete') for s in st):
        status = 'incomplete'
    else:
        status = 'passed'
    return dict(status=status, quantities=per, errors=errors, skipped=list(v.get('skipped', [])))


# --- quantities per model ----------------------------------------------------------------------

_V3 = ('pi_row_sum_err', 'pi_negative', 'mean_e_err', 'agg_income_err')


def _keys_rank(spec):
    return ['hh_budget', *_V3] + (['spec_r', 'spec_eis', 'spec_A', 'spec_a_min'] if spec else [])


def _v03_rank(v, out, spec, where):
    f = _Fields(out, ss_detail=('Y', 'C', 'A', 'r', 'eis', 'a_min', 'e_grid', 'Pi', 'y'))
    pi = _income(v, '', where, 'T1', f, 'e_grid', 'Pi', 'y', ['Y'], lambda: f['Y'])
    _eq(v, 'hh_budget', where, 'T1', f, ['C', 'A', 'r', 'y', 'Pi'],
        lambda: household_budget_resid(f['C'], f['A'], f['r'], f['y'], pi))
    if spec:
        _pin(v, 'spec_r', where, 'T1', f, 'r', spec['r'])
        _pin(v, 'spec_eis', where, 'T1', f, 'eis', spec['eis'])
        _pin(v, 'spec_A', where, 'T1', f, 'A', spec['A'], rel=True)
        _eq(v, 'spec_a_min', where, 'T1', f, ['a_min', 'e_grid'], lambda: _dev(
            f['a_min'], -spec['nbl_frac'] * spec['Y'] * np.min(f['e_grid']) / spec['r'], True))


def _keys_tank(spec):
    ks = ['ss_goods', 'ss_assets', 'incidence_err'] + \
        [f'{k}_{t}' for t in 'HS' for k in (*_V3, 'hh_budget')]
    if spec:
        ks += ['spec_lam', 'spec_beta_H', 'spec_r', 'spec_eis', 'spec_a_min_H', 'spec_a_min_S']
    return ks


def _v03_tank(v, out, spec, where):
    f = _Fields(out, ss_detail=['lam', 'r', 'eis', 'Y', 'Y_H', 'Y_S'] +
                [f'{k}_{t}' for t in 'HS' for k in ('C', 'A', 'beta', 'a_min', 'e_grid', 'Pi', 'y')])
    _eq(v, 'ss_goods', where, 'T2', f, ['Y', 'lam', 'C_H', 'C_S'],
        lambda: abs(f['Y'] - f['lam'] * f['C_H'] - (1 - f['lam']) * f['C_S']) / f['Y'])
    _eq(v, 'ss_assets', where, 'T2', f, ['lam', 'A_H', 'A_S', 'Y'],
        lambda: abs(f['lam'] * f['A_H'] + (1 - f['lam']) * f['A_S']) / f['Y'])
    for t in 'HS':
        pi = _income(v, f'_{t}', where, 'T2', f, f'e_grid_{t}', f'Pi_{t}', f'y_{t}', [f'Y_{t}'],
                     lambda: f[f'Y_{t}'])
        _eq(v, f'hh_budget_{t}', where, 'T2', f, [f'C_{t}', f'A_{t}', 'r', f'y_{t}', f'Pi_{t}'],
            lambda: household_budget_resid(f[f'C_{t}'], f[f'A_{t}'], f['r'], f[f'y_{t}'], pi))
    _eq(v, 'incidence_err', where, 'T2', f, ['lam', 'Y_H', 'Y_S', 'Y'],
        lambda: abs(f['lam'] * f['Y_H'] + (1 - f['lam']) * f['Y_S'] - f['Y']) / f['Y'])
    if spec:
        _pin(v, 'spec_lam', where, 'T2', f, 'lam', spec['lam'])
        _pin(v, 'spec_beta_H', where, 'T2', f, 'beta_H', spec['beta_H'])
        _pin(v, 'spec_r', where, 'T2', f, 'r', spec['r'])
        _pin(v, 'spec_eis', where, 'T2', f, 'eis', spec['eis'])
        _pin(v, 'spec_a_min_H', where, 'T2', f, 'a_min_H', spec['a_min_H'])
        _eq(v, 'spec_a_min_S', where, 'T2', f, ['a_min_S', 'e_grid_S'], lambda: _dev(
            f['a_min_S'], -spec['nbl_frac'] * spec['Y'] * np.min(f['e_grid_S']) / spec['r'], True))


# V1 (R4): quantity names by the GE target the T3 provider reports
_V1 = {'asset_mkt': ('ge_target_resid_reported_{}', 'ge_target_resid_direct_{}'),
       'goods_mkt': ('ge_target_resid_reported_goods_{}', 'ge_target_resid_direct_goods_{}')}


def _ge_target(out):
    """(target, assumed): the GE target a T3 provider reports; 'asset_mkt' if it reports none."""
    t = out.get('ge_target')
    return ('asset_mkt', True) if t is None else (t, False)


def _known_target(t):
    # a string naming one of the two markets; anything else (another string, a number, a list)
    # is outside the convention and V1 does not apply
    return isinstance(t, str) and t in _V1


def _keys_werning(out, cyc, spec):
    tag = 'cyc' if cyc else 'acyc'
    # an unknown target keeps the asset-market placeholders missing (V1 not applicable there)
    t = _ge_target(out)[0]
    ks = [k.format(tag) for k in _V1[t if _known_target(t) else 'asset_mkt']]
    ks += [f'{k}_{tag}' for k in ('ss_assets', 'ss_goods', 'hh_budget', *_V3)]
    if cyc:
        ks.append('incidence_err_cyc')
    if spec:
        ks += [f'spec_eis_{tag}', f'spec_a_min_{tag}', f'spec_rstar_{tag}',
               'spec_zeta_low_cyc' if cyc else 'spec_zeta_acyc']
    return ks


def _v03_werning(v, out, b, dr, cyc, spec, where):
    tag = 'cyc' if cyc else 'acyc'
    T = f'T3 {tag}'
    target, assumed = _ge_target(out)
    v.setdefault('ge_target', {}).setdefault(tag, []).append(
        target + (' (assumed)' if assumed else '') if isinstance(target, str) else
        f'{target!r} (not a string)')
    if not _known_target(target):
        msg = (f"V1 ({T}, {where}): applicability condition not met: GE target {target!r} is "
               f"neither 'asset_mkt' nor 'goods_mkt'")
        if msg not in v['skipped']:
            v['skipped'].append(msg)
    else:
        # V1(a): the residual of its target the provider reports; V1(b): the same residual
        # computed here from the reported Jacobians and the rate path the check passed
        reported, direct = (k.format(tag) for k in _V1[target])
        gt = out.get('diag', {}).get('ge_target_resid')
        if gt is None:
            _skip(v, f'{reported} ({T})', ["diag['ge_target_resid']"])
        else:
            _put(v, reported, gt / b, where)
        resid = dict(
            asset_mkt=('J_A', lambda J: np.asarray(J['Y']) @ out['dY'] + np.asarray(J['r']) @ dr),
            goods_mkt=('J_C', lambda J: out['dY'] - np.asarray(J['Y']) @ out['dY']
                       - np.asarray(J['r']) @ dr))
        for market, (field, fn) in resid.items():
            J = out.get(field)
            if market == target:
                if J is None:
                    _skip(v, f'{direct} ({T})', [field])
                else:
                    _put(v, direct, np.max(np.abs(fn(J))) / b, where)
            elif J is not None:      # the other market: diagnostic, no threshold
                _put(v.setdefault('diagnostics', {}), f'{market}_resid_direct_{tag}',
                     np.max(np.abs(fn(J))) / b, where)
    f = _Fields(out, top=('e_grid', 'Pi', 'zeta', 'eis'),
                ss_detail=('Y', 'C', 'A', 'r', 'a_min', 'y', 'beta'))
    _eq(v, f'ss_assets_{tag}', where, T, f, ['A'], lambda: abs(f['A']) / b)
    _eq(v, f'ss_goods_{tag}', where, T, f, ['C', 'Y'], lambda: abs(f['C'] - f['Y']) / f['Y'])
    pi = _income(v, f'_{tag}', where, T, f, 'e_grid', 'Pi', 'y', ['Y'], lambda: f['Y'])
    _eq(v, f'hh_budget_{tag}', where, T, f, ['C', 'A', 'r', 'y', 'Pi'],
        lambda: household_budget_resid(f['C'], f['A'], f['r'], f['y'], pi))
    e = lambda: np.asarray(f['e_grid'], float)
    zeta = lambda: np.asarray(f['zeta'], float)
    if cyc:
        _eq(v, 'incidence_err_cyc', where, T, f, ['e_grid', 'Pi', 'zeta'],
            lambda: abs(pi @ (e() * zeta())))
    if spec:
        _pin(v, f'spec_eis_{tag}', where, T, f, 'eis', spec['eis'])
        _eq(v, f'spec_a_min_{tag}', where, T, f, ['a_min'], lambda: _dev(f['a_min'], -b, True))
        if cyc:
            _eq(v, 'spec_zeta_low_cyc', where, T, f, ['zeta'],
                lambda: _dev(zeta()[0], spec['zeta_low_cyc']))
        else:
            _eq(v, 'spec_zeta_acyc', where, T, f, ['zeta'], lambda: np.max(np.abs(zeta())))
        # R2: beta from ss_detail, as documented
        _eq(v, f'spec_rstar_{tag}', where, T, f, ['beta', 'eis', 'e_grid', 'Pi'], lambda: _dev(
            zero_liquidity_rstar(f['beta'], f['eis'], e(), np.asarray(f['Pi'], float)),
            spec['r_star']))


def _keys_aiyagari(spec):
    ks = ['ss_capital', 'hh_budget', *_V3]
    if spec:
        ks += ['foc_r', 'foc_w', 'ss_goods', 'spec_beta', 'spec_eis', 'spec_a_min',
               'spec_euler_err_mean_weighted', 'spec_euler_err_max_with_mass',
               'spec_kkt_violation_max_with_mass']
    return ks


def _v03_aiyagari(v, out, spec, where):
    f = _Fields(out, ss_detail=('r', 'w', 'K', 'A', 'C', 'beta', 'eis', 'a_min', 'e_grid', 'Pi', 'y'),
                policy=('a_grid', 'a', 'c', 'D'))
    L_doc = spec['L'] if spec else 1.0
    pi = _income(v, '', where, 'T4', f, 'e_grid', 'Pi', 'y', ['w'], lambda: f['w'] * L_doc)
    _eq(v, 'ss_capital', where, 'T4', f, ['A', 'K'], lambda: abs(f['A'] - f['K']) / f['K'])
    _eq(v, 'hh_budget', where, 'T4', f, ['C', 'A', 'r', 'y', 'Pi'],
        lambda: household_budget_resid(f['C'], f['A'], f['r'], f['y'], pi))
    if not spec:
        return
    al, de = spec['alpha'], spec['delta']
    L = lambda: float(pi @ np.asarray(f['e_grid'], float))      # labour = mean of the income states
    Y = lambda: f['K'] ** al * L() ** (1 - al)
    _eq(v, 'foc_r', where, 'T4', f, ['r', 'K', 'e_grid', 'Pi'],
        lambda: abs(f['r'] + de - al * (f['K'] / L()) ** (al - 1)))
    _eq(v, 'foc_w', where, 'T4', f, ['w', 'K', 'e_grid', 'Pi'],
        lambda: abs(f['w'] - (1 - al) * (f['K'] / L()) ** al) / f['w'])
    _eq(v, 'ss_goods', where, 'T4', f, ['K', 'C', 'e_grid', 'Pi'],
        lambda: abs(Y() - f['C'] - de * f['K']) / Y())
    _pin(v, 'spec_beta', where, 'T4', f, 'beta', spec['beta'])
    _pin(v, 'spec_eis', where, 'T4', f, 'eis', spec['eis'])
    _pin(v, 'spec_a_min', where, 'T4', f, 'a_min', spec['a_min'])
    miss = f.missing(['a_grid', 'a', 'c', 'D', 'Pi', 'r'])
    if miss:
        return _skip(v, 'V4(b) Euler (T4)', miss)
    for k, x in euler_errors(f['Pi'], f['a_grid'], f['a'], f['c'], f['D'], f['r'],
                             spec['beta'], spec['eis']).items():
        _put(v, 'spec_' + k, x, where)


def _lag(x):
    return np.concatenate([[0.0], np.asarray(x, float)[:-1]])


def _keys_fiscal(spec):
    ks = ['ss_government_budget', 'ss_goods', 'ss_assets', 'hh_budget', *_V3]
    if spec:
        ks += [f'spec_{k}' for k in ('r', 'eis', 'B', 'G', 'phi_B', 'a_min')]
    return ks


# T5 path quantities and the irfs fields each needs (R2); the caller's inputs give the requested
# exogenous paths (R3): G, r (fixed-rate experiments) and rshock, zero where not requested
_PATH_FIELDS = {'fiscal_rule': ('dB', 'dTax'), 'government_budget': ('dB', 'dTax', 'dr', 'dG'),
                'goods_market': ('dY', 'dC', 'dG'), 'rate_rule': ('dr', 'dY'),
                'exog_G': ('dG',), 'exog_r': ('dr',), 'exog_rshock': ('drshock',)}


def _path_keys(kind):
    """Path quantities expected in an experiment of this kind ('fixed_r' or 'r_rule')."""
    return ['fiscal_rule', 'government_budget', 'goods_market', 'exog_G', 'exog_rshock',
            'rate_rule' if kind == 'r_rule' else 'exog_r']


def _v03_fiscal(v, out, shocks, spec, where):
    f = _Fields(out, ss_detail=('Y', 'C', 'A', 'B', 'G', 'Tax', 'r', 'eis', 'phi_B', 'a_min',
                                'e_grid', 'Pi', 'y'))
    _eq(v, 'ss_government_budget', where, 'T5', f, ['Tax', 'G', 'r', 'B', 'Y'],
        lambda: abs(f['Tax'] - f['G'] - f['r'] * f['B']) / f['Y'])
    _eq(v, 'ss_goods', where, 'T5', f, ['Y', 'C', 'G'], lambda: abs(f['Y'] - f['C'] - f['G']) / f['Y'])
    _eq(v, 'ss_assets', where, 'T5', f, ['A', 'B'], lambda: abs(f['A'] - f['B']) / f['B'])
    pi = _income(v, '', where, 'T5', f, 'e_grid', 'Pi', 'y', ['Y', 'Tax'], lambda: f['Y'] - f['Tax'])
    _eq(v, 'hh_budget', where, 'T5', f, ['C', 'A', 'r', 'y', 'Pi'],
        lambda: household_budget_resid(f['C'], f['A'], f['r'], f['y'], pi))
    if not spec:
        return
    for k in ('r', 'eis', 'B', 'G', 'phi_B', 'a_min'):
        _pin(v, f'spec_{k}', where, 'T5', f, k, spec[k])
    paths = v.setdefault('paths', {})
    for label, (kind, inputs) in shocks.items():        # every experiment the caller requested
        p = paths.setdefault(label, {})
        irf = out['irfs'].get(label)
        if irf is None:
            _skip(v, f'V4(b) paths (T5, {label})', [f"irfs['{label}']"])
            continue
        n = len(next(iter(inputs.values())))
        req = {k: np.asarray(inputs[k], float) if k in inputs else np.zeros(n)
               for k in ('G', 'r', 'rshock')}
        s = max(float(np.max(np.abs(x))) for x in inputs.values())
        x = {k: np.asarray(irf[k], float) for k in ('dY', 'dC', 'dG', 'dTax', 'dr', 'dB', 'drshock')
             if k in irf}
        resid = dict(
            fiscal_rule=lambda: x['dTax'] - spec['phi_B'] * _lag(x['dB']),
            government_budget=lambda: (x['dB'] - (1 + spec['r']) * _lag(x['dB']) - spec['B'] * x['dr']
                                       - x['dG'] + x['dTax']),
            goods_market=lambda: x['dY'] - x['dC'] - x['dG'],
            rate_rule=lambda: x['dr'] - spec['phi_y'] * x['dY'] - req['rshock'],
            exog_G=lambda: x['dG'] - req['G'],
            exog_r=lambda: x['dr'] - req['r'],
            exog_rshock=lambda: x['drshock'] - req['rshock'])
        for key in _path_keys(kind):
            miss = [f"irfs['{label}']['{k}']" for k in _PATH_FIELDS[key] if k not in x]
            if miss:
                _skip(v, f'{key} (T5 path {label})', miss)
            else:
                _put(p, key, np.max(np.abs(resid[key]())) / s, where)


# ---------------------------------------------------------------------------
# T1: RANK limit
# ---------------------------------------------------------------------------

def check_rank_limit(provider, eps_list, T, T_sub=(50, 100), spec=None, **kw):
    """PIH household-block limit (partial equilibrium; not a test of the full NK GE RANK
    limit): ||J_HA(eps) - J_PIH|| for J in {dC/dY, dC/dr}, J_PIH the permanent-income closed
    form at the provider's steady state (r, eis, C, A), along the path eps -> 0 on which beta is
    recalibrated to keep mean assets at the target.

    Also reports the error on the leading T_sub x T_sub block of the SAME T x T Jacobian
    (`err_*_by_window`; `err_*_by_T` is kept as an alias). This is an observation-window
    diagnostic: it shows that the gap sits at long lags. It is not a terminal-horizon audit;
    SSJ's household Jacobian entries do not depend on the computing horizon (checked in
    scripts/diagnostics_v0_2.py, section C: T=50 and T=100 Jacobians computed separately equal
    the leading blocks of the T=200 one exactly)."""
    t0 = time.time()
    res = dict(eps=[], err_M=[], err_Jr=[], err_M_fro=[], err_Jr_fro=[], kappa=[],
               mass_constrained=[], argmax_err_M=[], err_M_by_T={str(t): [] for t in T_sub},
               err_Jr_by_T={str(t): [] for t in T_sub}, ssj_forward_resid=[])
    v03 = _new_v03()
    for eps in eps_list:
        out = provider(eps, T, **kw)
        where = f'eps={eps:g}'
        v03['settings'].append(where)
        _expect(v03, _keys_rank(spec), where)
        _guard(v03, _v03_rank, out, spec, where)
        s = out['ss']
        M_ra, Jr_ra = ra_jacobians(T, s['r'], s['eis'], s['C'], s['A'])
        M, Jr = out['J']['C']['Y'], out['J']['C']['r']
        res['eps'].append(float(eps))
        res['err_M'].append(rel_max_err(M, M_ra))
        res['err_Jr'].append(rel_max_err(Jr, Jr_ra))
        res['err_M_fro'].append(rel_fro_err(M, M_ra))
        res['err_Jr_fro'].append(rel_fro_err(Jr, Jr_ra))
        i, j = np.unravel_index(np.argmax(np.abs(M - M_ra)), M.shape)
        res['argmax_err_M'].append([int(i), int(j)])
        for t in T_sub:
            res['err_M_by_T'][str(t)].append(rel_max_err(M[:t, :t], M_ra[:t, :t]))
            res['err_Jr_by_T'][str(t)].append(rel_max_err(Jr[:t, :t], Jr_ra[:t, :t]))
        d = out.get('diag', {})
        res['kappa'].append(float(d.get('kappa', np.nan)))
        res['mass_constrained'].append(float(d.get('mass_constrained', np.nan)))
        res['ssj_forward_resid'].append(float(d.get('ssj_forward_resid', np.nan)))
    pos = [k for k, e in enumerate(res['eps']) if e > 0]
    eps_pos = [res['eps'][k] for k in pos]
    for key in ('err_M', 'err_Jr'):
        vals = [res[key][k] for k in pos]
        res[f'rates_{key}'] = pairwise_rates(eps_pos, vals)
        res[f'fitted_rate_{key}'] = fitted_rate(eps_pos, vals)
    res['err_M_by_window'] = res['err_M_by_T']
    res['err_Jr_by_window'] = res['err_Jr_by_T']
    _finish(v03)
    res['v0_3'] = v03
    res['T'] = T
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T2: TANK limit
# ---------------------------------------------------------------------------

def check_tank_limit(provider, eps_list, T, dr, chis=(1.0, 1.5), spec=None, **kw):
    """TANK benchmark with an imposed Bilbiie-style income incidence chi (reduced form; chi is
    not derived from labour-market or profit-sharing primitives). Checks, as eps -> 0 for both
    types: (i) chi = 1: aggregate iMPC = lam I + (1-lam) M_RA; (ii) any chi: the GE-relevant
    matrix lam chi M_H + (1 - lam chi) M_S -> lam chi I + (1 - lam chi) M_RA; (iii) GE output
    for a fixed real-rate path = (1-lam)/(1-lam chi) x RANK. Conditions: same steady-state
    consumption and EIS, H types exactly hand-to-mouth in the limit, 1 - lam chi > 0."""
    t0 = time.time()
    res = dict(eps=[], err_MH_vs_I=[], err_MS_vs_MRA=[], err_Magg_vs_TANK=[],
               err_Magg_chi_vs_TANK={str(c): [] for c in chis},
               err_IRF_vs_TANK={str(c): [] for c in chis}, unconstrained_H=[], A_H=[],
               kappa_S=[], mass_constrained_S=[], ge_target_resid=[],
               amplification_HA={str(c): [] for c in chis})
    I = np.eye(T)
    dr, chis = _request(dr), _request(chis)       # the check's copies (section 11)
    v03 = _new_v03()
    for eps in eps_list:
        out = provider(eps, T, _request(dr), chis=_request(chis), **kw)
        where = f'eps={eps:g}'
        v03['settings'].append(where)
        _expect(v03, _keys_tank(spec), where)
        _guard(v03, _v03_tank, out, spec, where)
        lam, s = out['lam'], out['ss']
        M_ra, _ = ra_jacobians(T, s['r'], s['eis'], s['C_S'], s['A_S'])
        M_agg = lam * out['M_H'] + (1 - lam) * out['M_S']
        res['eps'].append(float(eps))
        res['err_MH_vs_I'].append(rel_max_err(out['M_H'], I))
        res['err_MS_vs_MRA'].append(rel_max_err(out['M_S'], M_ra))
        res['err_Magg_vs_TANK'].append(rel_max_err(M_agg, tank_impc(T, s['r'], lam, 1.0)))
        for c in chis:
            # GE-relevant aggregate iMPC with income incidence chi: lam*chi*M_H + (1-lam*chi)*M_S
            M_chi = lam * c * out['M_H'] + (1 - lam * c) * out['M_S']
            res['err_Magg_chi_vs_TANK'][str(c)].append(
                rel_max_err(M_chi, tank_impc(T, s['r'], lam, c)))
            ref = tank_output_fixed_real_rate(dr, s['r'], s['eis'], lam, c, s['Y'])
            res['err_IRF_vs_TANK'][str(c)].append(rel_max_err(out['dY'][float(c)], ref))
        for c in chis:
            # impact output relative to the chi=1 (equal incidence) economy
            res['amplification_HA'][str(c)].append(
                float(out['dY'][float(c)][0] / out['dY'][float(chis[0])][0]))
        d = out.get('diag', {})
        for k in ('unconstrained_H', 'A_H', 'kappa_S', 'mass_constrained_S', 'ge_target_resid'):
            res[k].append(float(d.get(k, np.nan)))
    res['rates_err_MS_vs_MRA'] = pairwise_rates(res['eps'], res['err_MS_vs_MRA'])
    res['fitted_rate_err_MS_vs_MRA'] = fitted_rate(res['eps'], res['err_MS_vs_MRA'])
    # the H-type error vanishes exactly below a finite eps, so a single fitted slope for the
    # aggregate would be meaningless; report pairwise rates only
    res['rates_err_Magg_vs_TANK'] = pairwise_rates(res['eps'], res['err_Magg_vs_TANK'])
    for c in chis:
        e = res['err_IRF_vs_TANK'][str(c)]
        res[f'rates_err_IRF_chi{c}'] = pairwise_rates(res['eps'], e)
        res[f'fitted_rate_err_IRF_chi{c}'] = fitted_rate(res['eps'], e)
    res['amplification_factor'] = {str(c): (1 - lam) / (1 - lam * c) / ((1 - lam) / (1 - lam * chis[0]))
                                   for c in chis}
    res['lam'] = lam
    _finish(v03)
    res['v0_3'] = v03
    res['T'] = T
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T3: Werning (2015) zero / vanishing liquidity
# ---------------------------------------------------------------------------

def check_werning(provider, b_list, T, dr, asymptotic_from=1e-3, spec=None, **kw):
    """For each b (borrowing limit -b, zero net supply, b -> 0):
      acyclical (income proportional to Y): error of dY_HA vs RANK  (Werning Prop. 2);
      countercyclical incidence: error vs the generalized zero-liquidity solution, and
      the (non-vanishing) deviation from RANK."""
    t0 = time.time()
    res = dict(b=[], err_acyc_vs_RANK=[], err_cyc_vs_ZL=[], dev_cyc_vs_RANK=[],
               zl_cyc_vs_RANK=[], r_gap_over_b_acyc=[], r_gap_over_b_cyc=[],
               mass_constrained=[], max_wealth_over_b=[])
    info_cyc = None
    dr = _request(dr)                             # the check's copy (section 11)
    v03 = _new_v03()
    for b in b_list:
        a = provider(b, T, _request(dr), cyclical=False, **kw)
        c = provider(b, T, _request(dr), cyclical=True, **kw)
        where = f'b={b:g}'
        v03['settings'].append(where)
        for out, cyc in ((a, False), (c, True)):
            _expect(v03, _keys_werning(out, cyc, spec), where)
            _guard(v03, _v03_werning, out, b, dr, cyc, spec, where)
        rank_a = ra_output_fixed_real_rate(dr, a['r_ss'], a['eis'])
        rank_c = ra_output_fixed_real_rate(dr, c['r_ss'], c['eis'])
        zl_c, info_cyc = zero_liquidity_output(dr, c['r_ss'], c['eis'], c['e_grid'], c['Pi'],
                                               c['zeta'])
        res['b'].append(float(b))
        res['err_acyc_vs_RANK'].append(rel_max_err(a['dY'], rank_a))
        res['err_cyc_vs_ZL'].append(rel_max_err(c['dY'], zl_c))
        res['dev_cyc_vs_RANK'].append(rel_max_err(c['dY'], rank_c))
        res['zl_cyc_vs_RANK'].append(rel_max_err(zl_c, rank_c))
        res['r_gap_over_b_acyc'].append(float(a['diag']['r_gap_over_b']))
        res['r_gap_over_b_cyc'].append(float(c['diag']['r_gap_over_b']))
        res['mass_constrained'].append(float(a['diag']['mass_constrained']))
        res['max_wealth_over_b'].append(float(a['diag']['max_wealth_over_b']))
    asym = [k for k, b in enumerate(res['b']) if b <= asymptotic_from]
    bs = [res['b'][k] for k in asym]
    for key in ('err_acyc_vs_RANK', 'err_cyc_vs_ZL'):
        res[f'rates_{key}'] = pairwise_rates(res['b'], res[key])
        res[f'fitted_rate_{key}_asymptotic'] = fitted_rate(bs, [res[key][k] for k in asym])
    res['asymptotic_from'] = asymptotic_from
    _finish(v03)
    res['v0_3'] = v03
    res['zero_liquidity_cyc'] = {k: (float(v) if not isinstance(v, int) else v)
                                 for k, v in info_cyc.items()}
    res['T'] = T
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T4: Aiyagari steady state
# ---------------------------------------------------------------------------

def check_aiyagari(provider, eps_list, rho_limit, sd_list, rho_list, sd_for_rho, rho_for_sd,
                   spec=None, **kw):
    """No-idiosyncratic-risk limit (r -> 1/beta - 1 as risk -> 0; the code has no insurance
    markets, so this is not a complete-markets limit) and comparative statics in risk size and
    persistence for THIS calibration (regression checks, not general theorems).
    `all_below_rho` is mechanical when the provider searches r below 1/beta - 1 (the SSJ
    provider does) and is not counted as accuracy evidence. Non-mechanical diagnostics, if the
    provider reports them in diag: capital-market residual, Euler/KKT residuals, top-grid mass."""
    t0 = time.time()
    res = dict(eps=[], r=[], gap=[], mass_constrained=[], capital_market_resid_rel=[],
               euler_err_mean_weighted=[], euler_err_max_with_mass=[], kkt_violation_max_with_mass=[],
               mass_top_5pts=[], mass_policy_off_grid=[])
    v03 = _new_v03()
    for eps in eps_list:
        out = provider(eps, rho_limit, **kw)
        where = f'eps={eps:g}'
        v03['settings'].append(where)
        _expect(v03, _keys_aiyagari(spec), where)
        _guard(v03, _v03_aiyagari, out, spec, where)
        res['eps'].append(float(eps))
        res['r'].append(out['r'])
        res['gap'].append(out['rho_disc'] - out['r'])
        d = out['diag']
        res['mass_constrained'].append(d['mass_constrained'])
        ee = d.get('euler', {})
        res['capital_market_resid_rel'].append(float(d.get('capital_market_resid_rel', np.nan)))
        res['euler_err_mean_weighted'].append(float(ee.get('euler_err_mean_weighted', np.nan)))
        res['euler_err_max_with_mass'].append(float(ee.get('euler_err_max_with_mass', np.nan)))
        res['kkt_violation_max_with_mass'].append(float(ee.get('kkt_violation_max_with_mass', np.nan)))
        res['mass_top_5pts'].append(float(d.get('mass_top_5pts', np.nan)))
        res['mass_policy_off_grid'].append(float(d.get('mass_policy_off_grid', np.nan)))
    res['rho_disc'] = out['rho_disc']
    res['gap_over_eps2'] = [g / e ** 2 for g, e in zip(res['gap'], res['eps'])]
    res['rates_gap'] = pairwise_rates(res['eps'], res['gap'])
    res['fitted_rate_gap'] = fitted_rate(res['eps'], res['gap'])
    res['r_by_sd'] = dict(sd=list(map(float, sd_list)), rho=rho_for_sd,
                          r=[provider(s, rho_for_sd, **kw)['r'] for s in sd_list])
    res['r_by_rho'] = dict(rho=list(map(float, rho_list)), sd=sd_for_rho,
                           r=[provider(sd_for_rho, rh, **kw)['r'] for rh in rho_list])
    allr = res['r'] + res['r_by_sd']['r'] + res['r_by_rho']['r']
    res['all_below_rho'] = bool(all(r < res['rho_disc'] for r in allr))
    res['all_below_rho_is_mechanical'] = True     # SSJ provider: root bracket below 1/beta - 1
    res['monotone_in_sd'] = bool(np.all(np.diff(res['r_by_sd']['r']) < 0))
    res['monotone_in_rho'] = bool(np.all(np.diff(res['r_by_rho']['r']) < 0))
    _finish(v03)
    res['v0_3'] = v03
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T5: intertemporal Keynesian cross
# ---------------------------------------------------------------------------

def ikc_rhs(J, dG, dTax, dr):
    """(I - M)^{-1} [dG + J^{C,Tax} dTax + J^{C,r} dr] with M = J^{C,Y}."""
    M = J['C']['Y']
    T = M.shape[0]
    partial = dG + J['C']['Tax'] @ dTax + J['C']['r'] @ dr
    return np.linalg.solve(np.eye(T) - M, partial), partial


def check_ikc(provider, T_list, shocks_fn, spec=None, **kw):
    """shocks_fn(T) -> {label: (kind, {input: path})}. For each T and shock, compare the
    provider's GE output path with the IKC formula built from its own block Jacobians."""
    t0 = time.time()
    res = dict(T=[], err={}, walras_asset_mkt={})
    v03 = _new_v03()
    for T in T_list:
        # the check's copy of the request: expected coverage first, then the provider gets its
        # own copy, and residuals and exogenous paths are computed from this one (section 11)
        shocks = _request(shocks_fn(T))
        where = f'T={T}'
        v03['settings'].append(where)
        _expect(v03, _keys_fiscal(spec), where)
        if spec:      # path quantities: every experiment the caller requested
            for label, (kind, _) in shocks.items():
                _expect(v03.setdefault('paths', {}).setdefault(label, {}), _path_keys(kind), where)
        out = provider(T, _request(shocks), **kw)
        _guard(v03, _v03_fiscal, out, shocks, spec, where)
        res['T'].append(int(T))
        for label, irf in out['irfs'].items():
            rhs, _ = ikc_rhs(out['J'], irf['dG'], irf['dTax'], irf['dr'])
            res['err'].setdefault(label, []).append(rel_max_err(irf['dY'], rhs))
            res['walras_asset_mkt'].setdefault(label, []).append(
                float(np.max(np.abs(irf['asset_mkt']))))
    res['max_err'] = float(max(max(v) for v in res['err'].values()))
    _finish(v03)
    res['v0_3'] = v03
    res['runtime_s'] = time.time() - t0
    return res


# ---------------------------------------------------------------------------
# T6: HANK - RANK differencing identity and decomposition
# ---------------------------------------------------------------------------

def check_differencing(J, ss, irf, n_dates=12):
    """Identity (I-M_HA)^{-1} - (I-M_RA)^{-1} = (I-M_HA)^{-1}(M_HA-M_RA)(I-M_RA)^{-1}, and
    decomposition dY_HA - dY_RA = (I-M_HA)^{-1}[(M_HA-M_RA)(dY_RA - dTax) + (Jr_HA-Jr_RA) dr]
    for the same steady state (r, eis, C, A) and the same policy paths (dG, dTax, dr)."""
    M = J['C']['Y']
    T = M.shape[0]
    I = np.eye(T)
    M_ra, Jr_ra = ra_jacobians(T, ss['r'], ss['eis'], ss['C'], ss['A'])
    inv_ha = np.linalg.inv(I - M)
    inv_ra = np.linalg.inv(I - M_ra)
    lhs = inv_ha - inv_ra
    rhs = inv_ha @ (M - M_ra) @ inv_ra
    dG, dTax, dr = irf['dG'], irf['dTax'], irf['dr']
    dY_ha = irf['dY']
    dY_ra = inv_ra @ (dG - M_ra @ dTax + Jr_ra @ dr)
    first_round = (M - M_ra) @ (dY_ra - dTax) + (J['C']['r'] - Jr_ra) @ dr
    dY_gap = inv_ha @ first_round
    q = 1 / (1 + ss['r'])
    disc = q ** np.arange(T)
    out = dict(
        identity_err=rel_max_err(rhs, lhs),
        cond_I_minus_M_RA=float(np.linalg.cond(I - M_ra)),
        cond_I_minus_M_HA=float(np.linalg.cond(I - M)),
        decomposition_err=rel_max_err(dY_gap, dY_ha - dY_ra),
        dates=list(range(n_dates)),
        dY_HA=_tolist(dY_ha[:n_dates]), dY_RA=_tolist(dY_ra[:n_dates]),
        dY_gap=_tolist(dY_gap[:n_dates]),
        first_round=_tolist(first_round[:n_dates]),
        ge_propagation=_tolist((dY_gap - first_round)[:n_dates]),
        cum_multiplier_HA=float(disc @ dY_ha / (disc @ dG)) if np.any(dG) else None,
        cum_multiplier_RA=float(disc @ dY_ra / (disc @ dG)) if np.any(dG) else None,
        impact_multiplier_HA=float(dY_ha[0] / dG[0]) if dG[0] != 0 else None,
        impact_multiplier_RA=float(dY_ra[0] / dG[0]) if dG[0] != 0 else None,
    )
    out['_arrays'] = dict(inv_ha=inv_ha, M_ra=M_ra, dY_ra=dY_ra, x=dY_ra - dTax)
    return out
