"""Closed-form representative-agent (RA), TANK and zero-liquidity benchmarks.

Timing convention = SSJ one-asset household: r_t is the ex-post return paid at date t
on assets carried from t-1, so
    budget : c_t + a_t = (1 + r_t) a_{t-1} + y_t
    Euler  : c_t^{-1/sigma} = beta (1 + r_{t+1}) c_{t+1}^{-1/sigma}.
All Jacobians below are the infinite-horizon objects truncated to t, s < T
(the same object that SSJ's fake-news algorithm computes).
"""
from __future__ import annotations

import numpy as np


def ra_jacobians(T: int, r: float, eis: float, C: float, A: float):
    """Permanent-income RA with beta (1+r) = 1.

    Returns (M, Jr) with
      M[t, s]  = dC_t/dY_s = (1-q) q^s,                              q = 1/(1+r)
      Jr[t, s] = dC_t/dr_s = (1-q) A q^s + eis C q (1{t>=s} - q^s).
    Derivation: Euler gives dC_{t+1} - dC_t = eis C q dr_{t+1}; the intertemporal budget
    constraint gives sum_t q^t dC_t = sum_t q^t (dY_t + A dr_t).
    """
    q = 1.0 / (1.0 + r)
    s = np.arange(T)
    qs = q ** s
    M = np.tile((1 - q) * qs, (T, 1))
    step = (np.arange(T)[:, None] >= s[None, :]).astype(float)
    Jr = (1 - q) * A * qs[None, :] + eis * C * q * (step - qs[None, :])
    return M, Jr


def ra_output_fixed_real_rate(dr: np.ndarray, r: float, eis: float, Y: float = 1.0):
    """RANK output with goods-market clearing Y = C, zero net asset supply, fixed real-rate path.

    dY_t = -eis * Y * q * sum_{s>t} dr_s   (Euler iterated forward with dY -> 0).
    """
    q = 1.0 / (1.0 + r)
    tail = np.concatenate([np.cumsum(dr[::-1])[::-1][1:], [0.0]])  # sum_{s>t} dr_s
    return -eis * Y * q * tail


def tank_output_fixed_real_rate(dr, r, eis, lam, chi, Y: float = 1.0):
    """Bilbiie TANK: hand-to-mouth share lam with income elasticity chi, savers = RA,
    zero liquidity. Aggregate Euler-IS: c_t = c_{t+1} - eis (1-lam)/(1-lam chi) r_{t+1}.
    """
    return (1 - lam) / (1 - lam * chi) * ra_output_fixed_real_rate(dr, r, eis, Y)


def tank_impc(T, r, lam, chi=1.0):
    """Aggregate iMPC matrix of TANK: lam*chi*I + (1 - lam*chi) * M_RA (chi=1: lam I + (1-lam) M_RA)."""
    M, _ = ra_jacobians(T, r, 1.0, 1.0, 0.0)
    return lam * chi * np.eye(T) + (1 - lam * chi) * M


def zero_liquidity_output(dr, r, eis, e_grid, Pi, zeta=None, Y: float = 1.0):
    """Output response under zero (vanishing) liquidity, Werning (2015) Section 3.

    Income y_it = Y_t e_i (Y_t/Y)^{zeta_i}. In the zero-liquidity equilibrium everyone
    consumes income and the real rate is set by the household state most willing to
    save, h = argmax_e E[(y'/y)^{-1/eis} | e]. Linearising that household's Euler:
        (1+zeta_h) yhat_t = (1 + zbar) yhat_{t+1} - eis * dr_{t+1}/(1+r),
        zbar = sum_e' w_e' zeta_e',  w_e' = Pi[h,e'] e'^{-1/eis} / sum Pi[h,.] e^{-1/eis}.
    With zeta == 0 (income proportional to aggregate income, Werning eq. 15 /
    Proposition 2) this is exactly the RA Euler equation.
    Returns (dY, info) with dY_t = Y yhat_t, imposing yhat_T = 0.
    """
    gam = 1.0 / eis
    n = len(e_grid)
    zeta = np.zeros(n) if zeta is None else np.asarray(zeta, float)
    ratio = (Pi @ e_grid ** (-gam)) / e_grid ** (-gam)       # E[(e'/e)^{-gam} | e]
    h = int(np.argmax(ratio))
    w = Pi[h] * e_grid ** (-gam)
    w = w / w.sum()
    zbar = float(w @ zeta)
    delta = (1 + zbar) / (1 + zeta[h])
    elas = eis / (1 + zeta[h])
    T = len(dr)
    yhat = np.zeros(T + 1)
    rhat = np.asarray(dr) / (1 + r)
    for t in range(T - 1, -1, -1):
        nxt = rhat[t + 1] if t + 1 < T else 0.0
        yhat[t] = delta * yhat[t + 1] - elas * nxt
    info = dict(h=h, delta=delta, elasticity=elas, zbar=zbar)
    return Y * yhat[:T], info


def zero_liquidity_rstar(beta, eis, e_grid, Pi):
    """Highest steady-state rate consistent with no saving (Werning eq. 12 with R_t = R*_t)."""
    gam = 1.0 / eis
    ratio = (Pi @ e_grid ** (-gam)) / e_grid ** (-gam)
    return 1.0 / (beta * ratio.max()) - 1.0
