"""Minimal reproduction of finding F1 (SSJ only; no code from this repository).

Observed behaviour (sequence-jacobian 1.0.0 on PyPI, and GitHub master eb4ac4c):
HetBlock linear impulses recompute the Jacobian with the block's own differentiation
options and do not use the `Js` that are passed in.

Test: one SSJ household block (sequence_jacobian.hetblocks.hh_sim.hh).
  J2 = two-sided Jacobian (h=1e-4);  J1 = default one-sided Jacobian (h=1e-4).
  (a) impulse_linear(..., Js={'hh': J2})            vs J2 @ dr   -> differs at O(h)
  (b) same call                                      vs J1 @ dr   -> agrees to rounding
  (c) impulse_linear(..., Js={'hh': 2 * J2})         -> output unchanged (Js not used)
  (d) impulse_linear(..., Js={'hh': J2}, twosided=True) vs J2 @ dr -> agrees (workaround)
Usage: python3 scripts/repro_F1_v0_2.py      Writes out/repro_F1_v0_2.json
"""
import os, sys, json, time, platform
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
from importlib.metadata import version
import numpy as np
import sequence_jacobian as sj
from sequence_jacobian.hetblocks.hh_sim import hh

T = 300
e, pi, Pi = sj.grids.markov_rouwenhorst(rho=0.9, sigma=0.5, N=7)
calib = dict(a_grid=sj.grids.asset_grid(0.0, 100.0, 300), y=e, r=0.01, beta=0.96, eis=0.5, Pi=Pi)
ss = hh.steady_state(calib)
dr = -0.01 * 0.8 ** np.arange(T)
J2 = hh.jacobian(ss, ['r'], ['C'], T=T, twosided=True, h=1e-4)
J1 = hh.jacobian(ss, ['r'], ['C'], T=T, h=1e-4)
ref2 = np.asarray(J2['C']['r']) @ dr
ref1 = np.asarray(J1['C']['r']) @ dr


def imp(Js, **kw):
    t0 = time.time()
    out = np.asarray(hh.impulse_linear(ss, {'r': dr}, outputs=['C'], Js=Js, **kw)['C'])
    return out, time.time() - t0


a, ta = imp({'hh': J2})
c, tc = imp({'hh': sj.JacobianDict({'C': {'r': 2 * np.asarray(J2['C']['r'])}}, name='hh', T=T)})
d, td = imp({'hh': J2}, twosided=True)
rel = lambda x, y: float(np.abs(x - y).max() / np.abs(y).max())
res = dict(
    a_impulse_with_J2_vs_J2_dr=rel(a, ref2),
    b_impulse_with_J2_vs_J1_dr=rel(a, ref1),
    c_impulse_with_2xJ2_vs_impulse_with_J2=rel(c, a),
    c_impulse_with_2xJ2_vs_2xJ2_dr=rel(c, 2 * ref2),
    d_impulse_with_J2_twosided_vs_J2_dr=rel(d, ref2),
    seconds=dict(impulse_with_Js=ta, impulse_with_poisoned_Js=tc, impulse_twosided=td),
    versions=dict(python=platform.python_version(), numpy=np.__version__,
                  sequence_jacobian=version('sequence-jacobian')),
    setup=dict(T=T, block='sequence_jacobian.hetblocks.hh_sim.hh', n_e=7, n_a=300, r=0.01,
               beta=0.96, eis=0.5, rho_e=0.9, sd_e=0.5, h=1e-4))
print(json.dumps(res, indent=1))
os.makedirs(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out'), exist_ok=True)
json.dump(res, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out',
                                 'repro_F1_v0_2.json'), 'w'), indent=1)
