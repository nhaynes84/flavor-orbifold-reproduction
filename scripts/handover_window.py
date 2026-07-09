#!/usr/bin/env python3
"""Scheme audit round 2, weapon 1: THE HAND-OVER WINDOW.

If a UV completion delivers the bare ratios 7 : 3.5 : 1 at scale Q, the
SM's beta functions drag the ratios off those values away from Q. The
observed lock therefore pins the hand-over scale to a computable window
— and every completion whose natural scale lies outside is PRE-EXECUTED.

One-loop SM running from measured M_Z values (two-loop moves crossings
~0.1 GeV — m3-verifier — so one-loop suffices for window edges).

Run: python3 consolidation/handover_window.py
"""
import math

B = {"1": 41/10, "2": -19/6, "3": -7.0}
INV_MZ = {"1": 59.012, "2": 29.569, "3": 8.4818}
MZ = 91.188

def inv(g, Q):
    return INV_MZ[g] - B[g]/(2*math.pi)*math.log(Q/MZ)

def ratios(Q):
    return inv("1", Q)/inv("2", Q), inv("2", Q)/inv("3", Q)

print("ratio drift rates (one-loop, at 88 GeV):")
eps = 1.001
r12a, r23a = ratios(88.0); r12b, r23b = ratios(88.0*math.e)
print(f"  r12 = 1/a1 : 1/a2 : {(r12b-r12a)/r12a*100:+.1f}% per e-fold")
print(f"  r23 = 1/a2 : 1/a3 : {(r23b-r23a)/r23a*100:+.1f}% per e-fold  <- the enforcer")

for tol in (0.01, 0.005):
    lo = hi = None
    Q = 20.0
    while Q < 1000:
        r12, r23 = ratios(Q)
        ok = abs(r12/2 - 1) < tol and abs(r23/3.5 - 1) < tol
        if ok and lo is None: lo = Q
        if ok: hi = Q
        Q *= 1.002
    if lo:
        print(f"\n  window where BOTH ratios hold within ±{tol:.1%}: "
              f"{lo:.1f} – {hi:.1f} GeV  (width x{hi/lo:.2f})")
    else:
        print(f"\n  ±{tol:.1%}: NO scale satisfies both (knife-edge)")

print("\nwhere the candidate hand-over scales sit:")
for name, Q in [("M_W", 80.38), ("crossing band", 86.3), ("M_Z", 91.19),
                ("rung 2 = v/2", 123.1), ("m_t", 172.5),
                ("1 TeV", 1000.0), ("10 TeV", 1e4), ("M_P-ish", 1e19)]:
    r12, r23 = ratios(Q)
    print(f"  {name:<14} Q = {Q:>8.4g} GeV   r12 = {r12:5.3f} ({r12/2-1:+7.1%})"
          f"   r23 = {r23:5.3f} ({r23/3.5-1:+8.1%})")

print("""
CONSEQUENCES FOR THE COMPLETION CLASS:
  1. Every completion whose natural hand-over scale is high (tower
     spacing, Planck, GUT) is PRE-EXECUTED — by 1 TeV the r23 ratio is
     off by ~25%+; at M_P the pattern is unrecognizable. 'The boundary
     delivers at its own deep scale and the ratios survive running' is
     impossible in this universe's beta functions.
  2. The window is an EW-band statement: the completion must contain
     the electroweak scale intrinsically. Rung 2's raw position (v/2 =
     123 GeV) sits OUTSIDE the +-1% window — the hand-over point is the
     gauge-boson band (~g2 x v/2), i.e., the COUPLING-DRESSED rung, not
     the bare rung. That gap (a factor g2) is real, unexplained, and
     stays inside C4's open box. Stated, not hidden.
  3. The alpha_3 knife-edge (m3-verifier) is the same fact seen from
     the input side: r23's 11%/e-fold rate is what makes the window
     tight AND makes band membership sensitive to alpha_3 at 1 sigma.""")
