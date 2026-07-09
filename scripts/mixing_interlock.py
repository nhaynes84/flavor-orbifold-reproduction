#!/usr/bin/env python3
"""THE GRAND INTERLOCK: every mixing parameter from masses alone.

The claim, stated so it stands independent of any model: the SM's
mixing sector (CKM + PMNS) is a FUNCTION of its mass sector. Inputs
below are MEASURED masses only — no mixing data, no model parameters.
Outputs: 7 mixing observables + 1 derived combination.

Provenance is marked per relation (literature / ours / grammar-risk),
per the trials-correction discipline of 05-LOOK-ELSEWHERE.

Run: python3 consolidation/mixing_interlock.py
"""

import math

# ---- inputs: masses only (MeV; leptons ppm-grade, quarks 5-20% scheme) ----
m_e, m_mu, m_tau = 0.51099895, 105.6584, 1776.86
m_u, m_c, m_t = 2.16, 1270.0, 172500.0
m_d, m_s, m_b = 4.67, 93.4, 4180.0

rows = []
def row(sector, name, formula, pred, obs, note):
    dev = (pred / obs - 1) * 100
    rows.append((sector, name, formula, pred, obs, dev, note))

# ---- LEPTON MIXING: from e, mu, tau alone (zero neutrino inputs) ----------
dx = 2/3 - math.log(m_mu/m_e) / math.log(m_tau/m_e)     # the muon offset
row("PMNS", "sin^2(th13)", "sqrt(3)*dx", math.sqrt(3)*dx, 0.0222,
    "ours; dx = muon node offset")
row("PMNS", "sin^2(th12)", "1/3 - 2*dx", 1/3 - 2*dx, 0.308,
    "ours; tribimaximal base + offset")
row("PMNS", "sin^2(th23)", "1/2 + (2+sqrt3)*dx", 0.5 + (2+math.sqrt(3))*dx, 0.546,
    "ours; THE LIVE BET (fits dispute: 0.470/0.561)")

# ---- QUARK MIXING: from quark masses alone --------------------------------
row("CKM", "sin(th12) [Cabibbo]", "sqrt(m_d/m_s)", math.sqrt(m_d/m_s), 0.2243,
    "LITERATURE: Gatto-Sartori-Tonin 1968")
row("CKM", "sin(th23) [V_cb]", "sqrt(m_u/m_c)", math.sqrt(m_u/m_c), 0.0415,
    "ours (alternating-sector form; Fritzsch used same-sector)")
row("CKM", "delta_CP (deg)", "arctan(m_d/m_u)",
    math.degrees(math.atan(m_d/m_u)), 65.5,
    "ours (doublet-angle mechanism candidate)")
row("CKM", "V_ub", "m_d/m_c", m_d/m_c, 0.00382,
    "GRAMMAR-RISK: from a documented ratio search; weakest, flagged")

# ---- derived combination (not independent) --------------------------------
s12, s23 = math.sqrt(m_d/m_s), math.sqrt(m_u/m_c)
s13 = m_d/m_c
c12, c23, c13 = (math.sqrt(1-s12**2), math.sqrt(1-s23**2), math.sqrt(1-s13**2))
delta = math.atan(m_d/m_u)
J = c12*c23*c13**2*s12*s23*s13*math.sin(delta)
row("CKM", "Jarlskog J (x1e-5)", "combination of the above", J*1e5, 3.08,
    "derived, not independent")

# ---- report ----------------------------------------------------------------
print(f"{'':4}{'observable':<22}{'formula':<24}{'pred':>10}{'obs':>10}{'dev':>8}")
print("-" * 86)
for sec, name, f, p, o, d, note in rows:
    print(f"{sec:<4}{name:<22}{f:<24}{p:>10.4g}{o:>10.4g}{d:>+7.1f}%")
    print(f"{'':4}  provenance: {note}")
print("-" * 86)
print("""SUMMARY: seven mixing observables + one derived combination, computed
from MASSES ALONE — no mixing inputs anywhere. Six of eight within ~1%
(quark entries limited by 5-20% scheme fuzz on their mass inputs; the
lepton entries use ppm-grade inputs). Combined with the mass model
(2 continuous dials + integers + v -> 9 masses), the ENTIRE flavor
sector -- 16 observables -- carries 2 dials.

HONEST GRADES: th12(CKM) is 1968 literature (GST) -- the pattern's
oldest confirmed member, which is corroboration, not theft. V_ub is
flagged grammar-risk (found by search). th23(PMNS) is under live
tension -- it is the package's falsifiable edge, not its weak point.
Every lepton-sector entry hangs on ONE number (the muon's node offset
dx = 0.0128), which is the same number that sets m_mu: one quantity,
four observables. That is the interlock.""")
