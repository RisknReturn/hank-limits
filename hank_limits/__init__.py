"""hank_limits: analytic-limit and aggregate-identity tests for heterogeneous-agent New Keynesian (HANK) codes.

Built on top of the public `sequence-jacobian` (SSJ) package.

Thread limits are set here (before numpy/numba are imported, when possible) so that
the suite stays light on a shared machine.
"""
import os as _os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    _os.environ.setdefault(_v, "1")

__version__ = "0.3.0"
