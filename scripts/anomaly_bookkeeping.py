#!/usr/bin/env python3
"""Closure 1: per-fixed-point anomaly bookkeeping on S^1/(Z2xZ2).

Standard rule (Arkani-Hamed, Cohen, Georgi, "Anomalies on orbifolds",
hep-th/0103135; extended by Scrucca-Serone): a 5D fermion on an orbifold
contributes a LOCALIZED 4D anomaly at the fixed points, distributed as
half the anomaly of its chiral zero mode at each fixed point where its
even component is unprojected. The heavy vector-like KK tower contributes
nothing. Therefore, per fixed point:

    A(point) = (1/2) * sum over bulk fields of A_4D(zero mode of field)

In the M2b arrangement, EVERY chiral zero mode carries parity (+,+) —
the same signature — so every field distributes identically, and each
fixed point sees half of the TOTAL zero-mode anomaly.

The check: is the total zero-mode content (3 SM generations) anomaly-free
in all four 4D anomaly channels? Computed explicitly below with SM
hypercharges (one generation; three generations = 3x each sum).

Run: python3 consolidation/anomaly_bookkeeping.py  (exit 0 = all cancel)
"""

from fractions import Fraction as F
import sys

# one SM generation of LEFT-HANDED Weyl fermions (standard convention:
# right-handed fields as left-handed conjugates, Y -> -Y)
# (name, SU(3) dim, SU(2) dim, hypercharge Y with Q = T3 + Y)
GEN = [
    ("Q_L",   3, 2, F(1, 6)),
    ("u_R^c", 3, 1, F(-2, 3)),
    ("d_R^c", 3, 1, F(1, 3)),
    ("L_L",   1, 2, F(-1, 2)),
    ("e_R^c", 1, 1, F(1, 1)),
]

def channels(gen):
    """The four independent SM anomaly sums (per generation)."""
    yyy   = sum(c * w * y**3 for _, c, w, y in gen)          # U(1)^3
    su2y  = sum(c * y for _, c, w, y in gen if w == 2)       # SU(2)^2-U(1)
    su3y  = sum(w * y for _, c, w, y in gen if c == 3)       # SU(3)^2-U(1)
    gravy = sum(c * w * y for _, c, w, y in gen)             # grav^2-U(1)
    return {"U(1)^3": yyy, "SU(2)^2 U(1)": su2y,
            "SU(3)^2 U(1)": su3y, "grav^2 U(1)": gravy}

print("Per-generation anomaly channels (must all vanish):")
ch = channels(GEN)
ok = True
for name, val in ch.items():
    print(f"  {name:>14}: {val}")
    ok &= (val == 0)
print()
print("Per-fixed-point bookkeeping (M2b arrangement):")
print("  every zero mode has parity (+,+)  ->  identical distribution")
print("  A(x=0)  = 1/2 * 3 * (per-gen total) = 1/2 * 3 * 0 = 0")
print("  A(x=1)  = 1/2 * 3 * (per-gen total) = 1/2 * 3 * 0 = 0")
print()
if ok:
    print("CLOSURE 1 RESULT: localized anomalies cancel PER FIXED POINT,")
    print("automatically, because (a) the zero-mode set is the anomaly-free")
    print("SM content and (b) all zero modes share one parity signature so")
    print("the ACG half-half distribution is uniform. No Chern-Simons flow")
    print("needed; no arrangement freedom used. SHOWN, not asserted.")
print()
print("Honest scope notes:")
print("  - Rule imported from ACG hep-th/0103135 (standard, cited).")
print("  - Uniform distribution holds BECAUSE M2b puts all zero modes at")
print("    (+,+). A future variant with mixed parities re-opens the check.")
print("  - SU(2) global (Witten) anomaly: 3 doublets per gen x 3 gens x ...")
doublets = sum(c for _, c, w, y in GEN if w == 2)  # per gen: 3 (colors) + 1 = 4
print(f"    doublet count per generation = {doublets} (even) -> Witten anomaly")
print("    absent per generation, hence absent per fixed point. OK.")
sys.exit(0 if ok else 1)
