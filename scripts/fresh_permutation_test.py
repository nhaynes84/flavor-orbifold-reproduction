#!/usr/bin/env python3
"""Closure 2: fresh null-test statistic for the 9-mass block.

Recomputed from scratch under the current framing (laws as in boundary.py,
current PDG-era inputs). Two nulls:

  Null A (weak, for completeness): permute measured masses among predicted
  slots. Tests ordering only; overclaims significance.

  Null B (the honest one): the DISCRETE-STRUCTURE ensemble. Freeze the
  laws' functional forms; enumerate the discrete choices the model could
  have made — sector-to-ladder-rung assignment (3! ways), gen-2 node
  indices k_d, k_u in 1..6, k_l in 1..2 — with NO continuous refitting.
  p_B = fraction of variants matching the data as well as the canonical
  assignment (log-RMS over all 9 masses).

Run: python3 consolidation/fresh_permutation_test.py
"""

import math
import itertools
import random
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from boundary import C0, M0, N_Q, N_L

E = math.e
MEASURED = {  # MeV, corpus-standard values
    "d": 4.67, "s": 93.4, "b": 4180.0,
    "u": 2.16, "c": 1270.0, "t": 172500.0,
    "e": 0.51100, "mu": 105.66, "tau": 1776.86,
}

def predict(q_map, k_d, k_u, k_l):
    """Nine masses from the frozen laws, given discrete choices.
    q_map: dict sector -> rung q in {1,2,3} (ladder g = e^{(q-1)/2})
    k_*: gen-2 interior node index per sector."""
    out = {}
    for sector, N, colored, names in (
        ("down", N_Q, True,  ("d", "s", "b")),
        ("up",   N_Q, True,  ("u", "c", "t")),
        ("lep",  N_L, False, ("e", "mu", "tau")),
    ):
        g = E ** ((q_map[sector] - 1) / 2)
        sigma = N * g - (C0 * math.sqrt(g) if colored else 0.0)
        m3 = M0 * N ** (7 / 3) * math.exp(5 * sigma / 6)
        if sector == "up":
            m3 *= (1 - C0**2)          # endpoint screening (law, up-type top)
        m1 = m3 * math.exp(-sigma)
        k = {"down": k_d, "up": k_u, "lep": k_l}[sector]
        partner = {"down": k_u, "up": k_d, "lep": None}[sector]
        if sector == "lep":
            K, sign = 6 / 5, -1        # law: toward the nu partner (x=0)
        else:
            K = partner
            sign = +1 if partner > k else -1
        dy = sign * C0 * K / (N * math.pi * k)
        y2 = k / N + dy
        m2 = m3 * math.exp(-sigma * (1 - y2))
        out[names[0]], out[names[1]], out[names[2]] = m1, m2, m3
    return out

def log_rms(pred):
    return math.sqrt(sum(
        (math.log(pred[p] / MEASURED[p])) ** 2 for p in MEASURED) / 9)

# canonical assignment
CANON = ({"down": 1, "up": 2, "lep": 3}, 3, 4, 2)
rms_true = log_rms(predict(*CANON))
print(f"Canonical assignment log-RMS = {rms_true:.4f} "
      f"(~{(math.exp(rms_true)-1)*100:.1f}% typical deviation)")

# ---- Null B: discrete-structure ensemble -----------------------------------
count, better = 0, []
for q_perm in itertools.permutations((1, 2, 3)):
    q_map = {"down": q_perm[0], "up": q_perm[1], "lep": q_perm[2]}
    for k_d in range(1, 7):
        for k_u in range(1, 7):
            for k_l in (1, 2):
                count += 1
                r = log_rms(predict(q_map, k_d, k_u, k_l))
                if r <= rms_true:
                    better.append(((q_perm, k_d, k_u, k_l), r))
pB = len(better) / count
print(f"\nNull B (discrete ensemble): {count} variants, "
      f"{len(better)} with RMS <= canonical -> p_B = {pB:.5f}")
for cfg, r in sorted(better, key=lambda x: x[1])[:5]:
    print(f"   {cfg}  RMS = {r:.4f}" + ("   <-- canonical" if r == rms_true else ""))

# ---- Null A: mass permutation (weak) ----------------------------------------
random.seed(3)
pred_c = predict(*CANON)
slots = list(MEASURED)
vals = [MEASURED[s] for s in slots]
hits, TRIALS = 0, 100000
for _ in range(TRIALS):
    perm = random.sample(vals, 9)
    r = math.sqrt(sum(
        (math.log(pred_c[s] / v)) ** 2 for s, v in zip(slots, perm)) / 9)
    if r <= rms_true:
        hits += 1
print(f"\nNull A (mass permutation, weak): p_A < {max(hits,1)/TRIALS:.0e} "
      f"({hits}/{TRIALS} hits) — ordering-only null, do not headline")

# ---- honest translation ------------------------------------------------------
def p_to_sigma(p):
    lo, hi = 0.0, 10.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2
if pB > 0:
    print(f"\nHEADLINE (Null B): the canonical assignment is a 1/{count} outlier")
    print("under the FROZEN laws — a LANDSCAPE statement. Do NOT translate to")
    print("sigma: no global significance attaches to the nine-mass fit because")
    print("the laws were developed with these masses in hand (statistician")
    print("verdict, 2026-07-07).")
print("""CAVEATS (state in the paper): the ensemble varies discrete choices
only — the laws' functional forms are frozen, so p_B answers 'given these
laws, was the discrete assignment special?' It does NOT price the laws
themselves (that is 05's grammar analysis + the ladder rigidity test).
The three quark masses carry 5-20% input fuzz; log-RMS is dominated by
structure, not precision, so the fuzz matters little here.""")
