#!/usr/bin/env python3
"""Interface dilutions: how each force's 4D coupling relates to the
boundary value g^2 = 3/2 (Laws 0-1 fix it; sin^2 = 2/9 gives e^2 = 1/3).

Key reframe (2026-07-06): the alpha_s 'success' is the statement that
COLOR crosses the boundary->4D interface with dilution ~ 1. The alpha_em
problem is that the circle/fold forces do NOT. So compute the dilution
function per force, from data, no freedom:

    D_i(Q) = g_boundary^2 / g_i^2(Q) = 0.11937 / alpha_i(Q)

and ask structural questions (not factor hunts):
  1. What are D_i at M_Z?
  2. Where (if anywhere) does each force cross D = 1?
  3. Where do the fold force and channel force couplings MEET, in the
     warp coordinate u? (Tower rung check — zero freedom.)

One-loop SM running, no thresholds (adequate to ~ +-0.5 e-fold).
Run: python3 consolidation/interface_dilutions.py
"""

import math
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from boundary import SIGMA_D, M0, u as warp_u

ALPHA_BOUNDARY = 1.5 / (4 * math.pi)      # 0.11937 = boundary alpha (all forces)
MZ_GEV = 91.1876

# PDG-ish M_Z values, GUT-normalized alpha_1 = (5/3) alpha_Y
ALPHA = {"U(1)_Y (GUT norm)": 0.016946, "SU(2) fold": 0.033819, "SU(3) channels": 0.1179}
B = {"U(1)_Y (GUT norm)": 41 / 10, "SU(2) fold": -19 / 6, "SU(3) channels": -7.0}

def inv_alpha(name, Q_gev):
    return 1 / ALPHA[name] - B[name] / (2 * math.pi) * math.log(Q_gev / MZ_GEV)

print("Interface dilution D_i = alpha_boundary / alpha_i   (boundary alpha = 0.1194)")
print("-" * 76)
for name in ALPHA:
    D = ALPHA_BOUNDARY / ALPHA[name]
    print(f"  {name:<22} alpha(M_Z) = {ALPHA[name]:.5f}   D(M_Z) = {D:6.3f}")
print("""
  STRUCTURE: the force that IS the surface (channels/color) passes the
  interface at D = 1.012 — no dilution beyond its own +1.2% residual.
  The geometric forces (fold, circle) arrive diluted 3.5x / 7x. Their
  boundary values are UNREACHABLE by 4D running (asymptotically free
  couplings never get that strong above confinement) -> the dilution is
  structural, at the interface itself. THIS is the alpha problem, posed
  correctly: derive D_fold and D_circle from the anatomy. lambda-gate:
  no post-hoc factor naming; a mechanism must predict the FUNCTION D(Q).""")

# ---- where do fold and channels meet? (zero-freedom rung check) ----------
# inv_alpha_i(t) = inv_alpha_i(0) - B_i t / 2pi  ->  crossing:
# t* = 2pi * (inv2 - inv3) / (B2 - B3)
d_inv = 1 / ALPHA["SU(2) fold"] - 1 / ALPHA["SU(3) channels"]
t_cross = 2 * math.pi * d_inv / (B["SU(2) fold"] - B["SU(3) channels"])
Q_cross = MZ_GEV * math.exp(t_cross)
a_cross = 1 / inv_alpha("SU(2) fold", Q_cross)
u_cross = warp_u(Q_cross * 1000)
print(f"  alpha_2 = alpha_3 crossing: Q = {Q_cross:.3g} GeV, alpha = {a_cross:.5f}")
print(f"  in warp coordinate: u = {u_cross:.3f}   (rung 7 = u = 7, i.e. "
      f"S(7) = {M0*math.exp(7*SIGMA_D)/1e3:.3g} GeV)")
print(f"  distance from rung: {abs(u_cross-7):.3f} rung units")
print(f"  honest grade: P(random point lands this close to an integer rung)"
      f" ~ {2*abs(u_cross-7):.2f} -> suggestive, NOT survives (1-loop, no thresholds)")

# U(1) crossing with SU(3), for completeness
d1 = 1 / ALPHA["U(1)_Y (GUT norm)"] - 1 / ALPHA["SU(3) channels"]
t1 = 2 * math.pi * d1 / (B["U(1)_Y (GUT norm)"] - B["SU(3) channels"])
Q1 = MZ_GEV * math.exp(t1)
print(f"  alpha_1 = alpha_3 crossing: Q = {Q1:.3g} GeV, u = {warp_u(Q1*1000):.3f}")

# ---- the discrete-volume pattern (see 12-INTERFACE.md) --------------------
print()
print("DISCRETE-VOLUME PATTERN: D = {N_q, N_q/2, 1} (circle, fold, channels)")
print("-" * 76)
for name, vol in [("U(1)_Y (GUT norm)", 7.0), ("SU(2) fold", 3.5),
                  ("SU(3) channels", 1.0)]:
    pred = ALPHA_BOUNDARY / vol
    res = (pred / ALPHA[name] - 1) * 100
    print(f"  {name:<22} predicted alpha = {pred:.6f}  measured "
          f"{ALPHA[name]:.6f}  ({res:+.2f}%)")
aY = 3 * ALPHA_BOUNDARY / 35
a2p = ALPHA_BOUNDARY / 3.5
sin2 = aY / (a2p + aY)
aem = a2p * sin2
print(f"  sin^2(thW) = {sin2:.8f} = 3/13 exactly (measured MS-bar 0.23122, "
      f"{(sin2/0.23122-1)*100:+.2f}%)")
print(f"  1/alpha_em(M_Z) = {1/aem:.2f} (measured 128.94, "
      f"{(aem*128.94-1)*100:+.2f}%)")
print(f"  1/alpha_em(0)  ~ {1/aem+8.09:.2f} (measured 137.036, "
      f"{(137.036/(1/aem+8.09)-1)*100:+.2f}%)")
