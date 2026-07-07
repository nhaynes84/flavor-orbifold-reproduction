#!/usr/bin/env python3
"""The warp tower, correctly labeled — and the rung-2 hit.

Fixes gravity_mechanism_v2.py PART 3, which conflated Lambda_BW*e^{k sigma}
with the clean ladder. The clean ladder is anchored at the mass-scale
bridge M0 and spaced by the (pre-frozen) warp:

    S(k) = M0 * e^{k * sigma_d},   sigma_d = 34/5, M0 = 384/2500 MeV

k = 0 and k = 1 are already fixed by the mass sector (M0 anchor; k=1 is
Lambda_BW = 137.9 MeV, the QCD scale). The question with ZERO remaining
freedom: what sits at k = 2?

Run: python3 consolidation/warp_tower.py
"""

import math

C0 = 0.2
H = 25
SIGMA = 34 / 5
M0 = (4 / H) * (1 - C0**2)      # MeV

KNOWN = [  # (name, value in MeV)
    ("M0 (bridge)",        0.1536),
    ("m_pi0",              134.98),
    ("m_pi+-",             139.57),
    ("T_c/2 (lattice)",    135.0),
    ("v/2",                123_110.0),
    ("M_H",                125_250.0),
    ("M_Z",                91_188.0),
    ("v",                  246_220.0),
    ("seesaw (typ 1e10-1e15 GeV)", 1e16),   # band marker only
    ("GUT (2e16 GeV)",     2e19),
    ("M_P",                1.220890e22),
]

print("Warp tower  S(k) = M0 * e^(k*sigma_d):")
print("-" * 66)
for k in range(0, 8):
    s = M0 * math.exp(k * SIGMA)
    # nearest known landmark
    name, v = min(KNOWN, key=lambda kv: abs(math.log(s / kv[1])))
    off = (s / v - 1) * 100
    unit_s = f"{s:.4g} MeV" if s < 1e3 else (f"{s/1e3:.4g} GeV" if s < 1e9 else f"{s/1e9:.4g} PeV+")
    print(f"  k={k}:  {unit_s:>14}   nearest known: {name:<26} ({off:+.1f}%)")
print()
print("Gravity does NOT sit on an integer rung: ln(M_P/M0) / sigma ="
      f" {math.log(1.220890e22 / M0) / SIGMA:.3f}")
print("  -> consistent with v2's local-vs-global split: local scales ride the")
print("     integer tower; gravity rides the global law sigma + sigma^2")
print(f"     (exponent check: sigma + sigma^2 = {SIGMA + SIGMA**2:.2f} vs"
      f" ln(M_P/M0) = {math.log(1.220890e22 / M0):.2f}, residual = ln P)")
print()

# ---- honest surprise for the rung-2 hit -----------------------------------
# Rung 2 has zero freedom given k=0,1 fixed by the mass sector. It lands at
# 123.8 GeV. The EW landmark window (v/2 = 123.11, M_H = 125.25): distances:
s2 = M0 * math.exp(2 * SIGMA)
for name, v in [("v/2", 123_110.0), ("M_H", 125_250.0)]:
    print(f"  rung 2 = {s2/1e3:.2f} GeV vs {name} = {v/1e3:.2f} GeV: "
          f"{(s2/v-1)*100:+.2f}%")
# Surprise: the a priori window for 'the next fundamental scale above QCD'
# spans ln(S) in ~[ln 1 GeV, ln 1e19 GeV] ~ 44 e-folds. Landing within
# |ln(s2/v2)| = 0.0057 of the EW scale:
window = math.log(1e19)  # ~44
hit = 2 * abs(math.log(s2 / 123_110.0))
p_raw = hit / window
# trials: exponent k in 1..7 and sector sigma in {sigma_d, sigma_u, sigma_l}
trials = 7 * 3
print(f"  raw p (land this close in a {window:.0f} e-fold window): {p_raw:.2e}")
print(f"  trials-corrected (k x sector = {trials}): p ~ {p_raw*trials:.3f}")
print("""
  READING: ~0.6% trials-corrected — a genuine hit, sitting between the
  'marginal' and 'survives' gates. The corpus had this as 'M_H = M0 e^{2 sigma_d},
  warp tower level 2, -1.14%' and underrated it: with k=0,1 pinned by the
  mass sector, rung 2 was a zero-freedom PREDICTION of where the EW scale
  sits. Ambiguity: rung 2 is 0.6% from v/2 and 1.1% from M_H — it locates
  the EW SCALE, not a specific particle. Prospective content: rung 3 at
  ~111 TeV is the first NEW rung — nothing known lives there; any future
  fundamental scale must land on a rung or v2's tower is strained.""")
