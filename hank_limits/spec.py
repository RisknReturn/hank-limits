"""Documented specification of the suite's models (v0.3; checker side).

The values are those of the README "Exact setting" paragraphs. They are written out here and
NOT imported from `hank_limits.models` (`PIH_DEFAULTS`, `AIYAGARI_DEFAULTS`, `FISCAL_DEFAULTS`,
...): a check that read the provider's own defaults would compare a mutated default with
itself. Grid sizes, grid tops and the income-process parameters are numerical choices or enter
through the reported chain, and are not pinned (docs/CHECKS_PLAN_v0_3.md, section 3).

The checks in `hank_limits.checks` take one of these dicts as `spec=`. A code under test with a
different calibration passes its own values; the pins are regression checks for this suite.
"""

# T1: PIH household block. Borrowing limit = nbl_frac x natural limit -Y e_min / r.
T1 = dict(r=0.03, eis=0.5, A=5.0, Y=1.0, nbl_frac=0.99)

# T2: TANK. H: borrowing limit 0, beta_H 0.95. S: the T1 household with zero net supply.
T2 = dict(lam=0.3, beta_H=0.95, a_min_H=0.0, r=0.03, eis=0.5, Y=1.0, nbl_frac=0.99)

# T3: Werning vanishing liquidity. Borrowing limit -b; zeta_low = 1 in the countercyclical case.
T3 = dict(eis=0.5, r_star=0.01, zeta_low_cyc=1.0, Y=1.0)

# T4: Aiyagari. Labour 1 (mean of the income states).
T4 = dict(beta=0.96, alpha=0.36, delta=0.08, eis=0.5, a_min=0.0, L=1.0)

# T5/T6: fiscal HANK. Tax_t = Tax + phi_B (B_{t-1} - B); rule r_t = r + phi_y (Y_t - 1) + shock.
T5 = dict(r=0.02, eis=0.5, B=1.0, G=0.2, phi_B=0.1, phi_y=0.25, a_min=0.0, Y=1.0)
