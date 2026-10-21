"""F3 end to end in HetBlock.jacobian(..., twosided=True), on SSJ's own example block hh_sim (backward function unchanged).
H: nonlinear hetinput y = w^(1+zeta_e) e_grid; Jacobian of A with respect to w.
O: linear hetinput y = w e_grid plus the hetoutput c2 = c**2; Jacobian of C2 with respect to r.
'original' = SSJ as installed; 'stored default' = DifferentiableCombinedExtendedFunction.diff called with twosided=None,
so the twosided value stored by .differentiable(...) is used. The gap between the two at the same h isolates the
one-sided hetinput/hetoutput term from every other differentiation error; it should scale with h."""
import numpy as np
from sequence_jacobian.hetblocks.hh_sim import hh, make_grids
from sequence_jacobian.utilities import function as F

cal = dict(min_a=0, max_a=1000, rho_e=0.975, sd_e=0.7, n_a=200, n_e=7, w=1.0, r=0.01/4, beta=1-0.08/4, eis=1)


def income_nl(w, e_grid):
    zeta = np.linspace(0.6, -0.4, e_grid.size)
    y = w ** (1 + zeta) * e_grid
    return y


def income_lin(w, e_grid):
    y = w * e_grid
    return y


def c2(c):
    c2 = c ** 2
    return c2


orig = F.DifferentiableCombinedExtendedFunction.diff


def stored_default(self, shock_dict, h=None, outputs=None, hide_zeros=False, twosided=None):
    return orig(self, shock_dict, h, outputs, hide_zeros, twosided)


def jac(blk, ss, i, o, h, fixed):
    F.DifferentiableCombinedExtendedFunction.diff = stored_default if fixed else orig
    try:
        return np.asarray(blk.jacobian(ss, [i], [o], T=60, twosided=True, h=h)[o][i])
    finally:
        F.DifferentiableCombinedExtendedFunction.diff = orig


cases = [("H hetinput ", hh.add_hetinputs([income_nl, make_grids]), "w", "A"),
         ("O hetoutput", hh.add_hetinputs([income_lin, make_grids]).add_hetoutputs([c2]), "r", "C2")]
for name, blk, i, o in cases:
    ss = blk.steady_state(cal)
    gaps = []
    for h in (1e-3, 1e-4, 1e-5):
        a, b = jac(blk, ss, i, o, h, False), jac(blk, ss, i, o, h, True)
        gaps.append(np.abs(a - b).max() / np.abs(b).max())
    print(name, "max|J_original - J_stored_default| / max|J_stored_default|:",
          "  ".join(f"h={h:.0e} {g:.2e}" for h, g in zip((1e-3, 1e-4, 1e-5), gaps)))
