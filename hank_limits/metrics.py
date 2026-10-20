"""Error norms and observed convergence rates."""
from __future__ import annotations

import numpy as np


def rel_max_err(x, ref) -> float:
    """max |x - ref| / max |ref| (entrywise, works for vectors and matrices).

    This is the largest entry difference relative to the largest reference entry; it is not
    an induced matrix norm and does not protect the relative accuracy of small (e.g. long-lag)
    entries. Use `check_close` where a reference is near zero."""
    x, ref = np.asarray(x, float), np.asarray(ref, float)
    return float(np.max(np.abs(x - ref)) / np.max(np.abs(ref)))


def rel_fro_err(x, ref) -> float:
    x, ref = np.asarray(x, float), np.asarray(ref, float)
    return float(np.linalg.norm(x - ref) / np.linalg.norm(ref))


def pairwise_rates(params, errors):
    """Observed order p_i = log(e_i / e_{i+1}) / log(h_i / h_{i+1}) between successive settings."""
    h = np.asarray(params, float)
    e = np.asarray(errors, float)
    return [float(np.log(e[i] / e[i + 1]) / np.log(h[i] / h[i + 1])) for i in range(len(h) - 1)]


def fitted_rate(params, errors):
    """Least-squares slope of log(error) on log(param)."""
    h = np.log(np.asarray(params, float))
    e = np.log(np.asarray(errors, float))
    A = np.vstack([h, np.ones_like(h)]).T
    slope, _ = np.linalg.lstsq(A, e, rcond=None)[0]
    return float(slope)


def check_close(x, ref, rtol, atol) -> dict:
    """Mixed criterion max|x - ref| <= atol + rtol * max|ref| (for references near zero)."""
    x, ref = np.asarray(x, float), np.asarray(ref, float)
    err = float(np.max(np.abs(x - ref)))
    bound = float(atol + rtol * np.max(np.abs(ref)))
    return dict(max_abs_err=err, bound=bound, ok=bool(err <= bound))
