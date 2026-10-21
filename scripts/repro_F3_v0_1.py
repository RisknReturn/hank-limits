"""Minimal reproduction (SSJ code only): HetBlock asks for two-sided differentiation of hetinputs, but
DifferentiableCombinedExtendedFunction.diff defaults to twosided=False (not None), so the stored default is ignored and
the derivative used by backward_fakenews is one-sided. sequence-jacobian 1.0.0 (= master eb4ac4c)."""
import numpy as np
from sequence_jacobian.utilities.function import CombinedExtendedFunction


def income(Y, e_grid, zeta):
    y = e_grid * Y ** (1 + zeta)
    return y


e = np.array([0.5, 1.0, 1.5]); zeta = np.array([0.6, 0.0, -0.4]); ss = dict(Y=1.0, e_grid=e, zeta=zeta)
analytic = e * (1 + zeta)
for h in (1e-2, 1e-3, 1e-4):
    d = CombinedExtendedFunction([income]).differentiable(ss, h, True)   # as in HetBlock.jac_backward_prelim
    used = d.diff({'Y': 1.0})['y']                                       # as in HetBlock.backward_fakenews
    fwd = (income(1 + h, e, zeta) - income(1, e, zeta)) / h
    ctr = (income(1 + h, e, zeta) - income(1 - h, e, zeta)) / (2 * h)
    alt = d.diff({'Y': 1.0}, twosided=None)['y']                          # stored default honoured
    print(f"h={h:g}: |used-forward| {np.abs(used - fwd).max():.1e}  |used-central| {np.abs(used - ctr).max():.1e}  "
          f"|used-analytic| {np.abs(used - analytic).max():.1e}  |twosided=None - analytic| {np.abs(alt - analytic).max():.1e}")
