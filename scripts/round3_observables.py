#!/usr/bin/env python3
"""Scheme audit round 3: the conversion in OBSERVABLE language.

Exact decomposition: delta(1/alpha_i) = ell_i*Delta + (group piece).
  - ell_i*Delta rescales the ONE bare coupling -> ratios EXACTLY
    preserved -> the ratios are BLIND to Delta. Delta lives in one
    observable only: the overall level (alpha_b offset).
  - the group piece is therefore carried ENTIRELY by a second
    observable: the SEPARATION of the two ratio crossings
    (r12 = 2 at 84.5 GeV vs r23 = 3.5 at 88.1 GeV).
Two parameters, two observables. This enables two executions:

KILL 1 — "the boundary scheme IS MSbar" (zero conversion): dead iff the
level offset is significant. delta(1/alpha_1) = +0.369 +- 0.02 -> the
U(1) level pins Delta != 0 at ~18 sigma. EXECUTED.

SOFT SPOT — the group piece's smoking gun (crossing separation,
0.042 e-folds) vs the alpha_3 input error (which moves Q23 by 8.7 x
d(alpha_3)/alpha_3): is the separation significant?

Run: python3 consolidation/round3_observables.py
"""
import math

B = {"1": 41/10, "2": -19/6, "3": -7.0}
INV = {"1": 59.012, "2": 29.569, "3": 8.4818}
MZ = 91.188

def inv(g, Q, d3=0.0):
    base = INV[g] + (d3 if g == "3" else 0.0)
    return base - B[g]/(2*math.pi)*math.log(Q/MZ)

def crossing(pair, target, d3=0.0):
    lo, hi = 20.0, 500.0
    for _ in range(80):
        mid = math.sqrt(lo*hi)
        r = inv(pair[0], mid, d3)/inv(pair[1], mid, d3)
        if (r - target) * (inv(pair[0], lo, d3)/inv(pair[1], lo, d3) - target) > 0:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo*hi)

print("1. Ratios are blind to Delta (exact):")
print("   1/alpha_i = ell_i (1/alpha_b + Delta)  ->  ratios = 7 : 3.5 : 1")
print("   identically, for ANY Delta. Delta is the LEVEL observable only.")

q12 = crossing(("1", "2"), 2.0)
q23 = crossing(("2", "3"), 3.5)
print(f"\n2. The group piece's observable — crossing separation:")
print(f"   Q(r12=2) = {q12:.1f} GeV,  Q(r23=3.5) = {q23:.1f} GeV")
print(f"   separation = {math.log(q23/q12):+.4f} e-folds")

print(f"\n3. KILL — 'boundary scheme = MSbar' (zero conversion):")
d1 = INV["1"] - 7/(1.5/(4*math.pi))
print(f"   requires level offset = 0; measured delta(1/alpha_1) = {d1:+.3f} +- 0.02")
print(f"   -> Delta != 0 at ~{abs(d1)/0.02:.0f} sigma vs measurement error")
print("   (conditional on the tree level alpha_b = 3/8pi). EXECUTED.")

print(f"\n4. SOFT SPOT — is the crossing separation significant vs alpha_3?")
# alpha_s = 0.1179 +- 0.0009 -> d(1/alpha_3) = -+ 0.065
for d3, tag in ((0.0, "central"), (-0.065, "alpha_3 +1sig"), (0.065, "alpha_3 -1sig")):
    q23x = crossing(("2", "3"), 3.5, d3)
    print(f"   {tag:<14} Q23 = {q23x:6.1f} GeV   separation = "
          f"{math.log(q23x/q12):+.3f} e-folds")
print("""   -> at alpha_3 +1 sigma the crossings COINCIDE (separation ~0) or
   invert: the group piece's existence is a ~1.5-2 sigma statement with
   today's alpha_s. HONEST STATUS: Delta is rock-solid (18 sigma); the
   group structure (k vs universal vs zero) is SOFT until alpha_s
   improves (FCC-ee projects x5-10 on alpha_s -> decisive).

5. RING-COLLAPSE RISK (story 4's bill): M0's triage killed propagating
   internal modes (Landau pole, Z' limits). No propagating ring modes
   -> no ring loops -> the cell/ring regulator cannot generate its
   conversion through mode loops. Story (4) therefore stands ONLY if
   the capacity bound (C7/M6) provides a non-loop derivation of the
   toll — else the ring class collapses into the topological class
   (Delta arbitrary, M4 deflates). The audit's pressure lands on M6.""")
