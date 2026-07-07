#!/usr/bin/env python3
"""Quantify the scale-overdetermination of the Cycle Rule ratio hits.

Claim: 1/alpha_1 : 1/alpha_2 : 1/alpha_3 = 7 : 3.5 : 1 at the framework's
transparency scale. Two independent ratio conditions, one scale degree of
freedom. If the two conditions were unrelated, the scales where each holds
exactly would differ by decades (they run at very different rates). How
close are they actually?

Run: python3 consolidation/overdetermination_check.py
"""

import math

MZ = 91.1876  # GeV
# PDG-ish M_Z inputs (GUT-normalized alpha_1)
INV = {"1": 1 / 0.016946, "2": 1 / 0.033819, "3": 1 / 0.1179}
B = {"1": 41 / 10, "2": -19 / 6, "3": -7.0}   # 1-loop SM

def inv_alpha(i, Q):
    return INV[i] - B[i] / (2 * math.pi) * math.log(Q / MZ)

def solve_ratio(i, j, target):
    """Scale where inv_alpha_i / inv_alpha_j = target (Newton iteration)."""
    Q = MZ
    for _ in range(60):
        f = inv_alpha(i, Q) / inv_alpha(j, Q) - target
        h = 1e-4
        fp = (inv_alpha(i, Q * (1 + h)) / inv_alpha(j, Q * (1 + h)) - target - f) / (
            math.log(1 + h))
        Q *= math.exp(-f / fp)
    return Q

Q12 = solve_ratio("1", "2", 2.0)     # 1/a1 = 2 * 1/a2
Q23 = solve_ratio("2", "3", 3.5)     # 1/a2 = 3.5 * 1/a3
print("Scale where 1/alpha_1 : 1/alpha_2 = 2   exactly:  Q =", f"{Q12:8.2f} GeV")
print("Scale where 1/alpha_2 : 1/alpha_3 = 3.5 exactly:  Q =", f"{Q23:8.2f} GeV")
print(f"Separation: {abs(math.log(Q12/Q23)):.3f} e-folds "
      f"({max(Q12,Q23)/min(Q12,Q23):.3f}x)")

# how far COULD they have been? the ratios drift apart over the running
# range; a natural a priori window is the full desert, ~ ln(1e16/1e2) = 32
# e-folds. Probability two independent conditions land within the observed
# separation by chance ~ separation / window:
sep = abs(math.log(Q12 / Q23))
window = math.log(1e16 / 1e2)
print(f"a priori window ~ {window:.0f} e-folds -> chance coincidence "
      f"p ~ {sep/window:.4f}")
print(f"Both conditions hold at Q ~ {math.sqrt(Q12*Q23):.1f} GeV — the EW "
      f"scale (M_W = 80.4, M_Z = 91.2), the framework's pre-registered "
      f"transparency scale (alpha_s = 3/(8pi) at M_Z predates this pattern).")
