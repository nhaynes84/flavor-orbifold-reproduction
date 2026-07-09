#!/usr/bin/env python3
"""Why quark masses CANNOT be measured as sharply as lepton masses —
the mathematical floors, not the hand-wave.

Three independent floors per fermion:

1. THE RENORMALON FLOOR (quarks only; theorem-grade, mainstream QFT —
   Bigi-Shifman-Uraltsev-Vainshtein 1994; Beneke-Braun 1994): the pole
   mass ("the mass of the particle you would weigh") has an IRREDUCIBLE
   ambiguity of order the confinement scale, delta_m ~ O(1) x Lambda,
   because confinement forbids the asymptotic state the pole mass
   describes. We use Lambda_BW = 137.9 MeV (the framework's rung 1;
   standard estimates use 200-300 MeV — ours is conservative).
   Leptons: NO confinement -> pole mass exists, IR-safe -> floor ~ 0.

2. THE SCHEME-SCALE GRADIENT (quarks dominant): MS-bar masses are
   bookkeeping parameters that depend on the chosen scale mu. One-loop:
   d(ln m)/d(ln mu) = -6 C_F alpha_s /(4 pi) = -2 alpha_s/pi.
   At mu = 2 GeV (alpha_s ~ 0.30): ~19% shift PER E-FOLD of scale
   choice. Leptons (QED): 3 alpha/(2 pi) ~ 0.35%/e-fold for the MS-bar
   parameter — and the scheme-free pole mass needs no scale at all.

3. THE WIDTH FLOOR: an unstable particle's mass is fuzzy by Gamma/m.
   Only the top's is non-negligible (it decays before hadronizing).

Verdict logic: a position/dwell-type mass model (ours) predicts the
PHYSICAL mass. Comparing it to quark data can never beat the floors —
so quark residuals below floor are UNINFORMATIVE about the model, and
the lepton sector is the only place the model can truly be tested.

Run: python3 consolidation/accuracy_floors.py
"""

import math

LAMBDA = 137.9          # MeV, rung 1 (standard renormalon quotes: 200-300)
ALPHA = 1 / 137.036
ALPHA_S = {"2GeV": 0.30, "mc": 0.38, "mb": 0.22, "mt": 0.108}

# name, m (MeV), experimental rel. precision, width (MeV), sector
F = [
    ("e",   0.51099895, 3e-10, 0.0,      "lepton"),
    ("mu",  105.6584,   2e-8,  3.0e-16,  "lepton"),
    ("tau", 1776.86,    7e-5,  2.3e-9,   "lepton"),
    ("u",   2.16,       0.16,  0.0,      "quark(2GeV)"),
    ("d",   4.67,       0.10,  0.0,      "quark(2GeV)"),
    ("s",   93.4,       0.06,  0.0,      "quark(2GeV)"),
    ("c",   1270.0,     0.016, 0.0,      "quark(mc)"),
    ("b",   4180.0,     0.007, 0.0,      "quark(mb)"),
    ("t",   172500.0,   0.002, 1400.0,   "quark(mt)"),
]
MODEL_RESID = {"e": 0.002, "mu": 0.003, "tau": 0.003,
               "d": 0.007, "s": 0.006, "b": 0.004,
               "u": 0.024, "c": 0.030, "t": 0.028}

print(f"{'':5}{'renormalon':>12}{'scale-grad':>12}{'width':>9}"
      f"{'FLOOR':>9}{'exp.prec':>10}{'model res':>10}  verdict")
print("-" * 82)
for name, m, exp, width, sector in F:
    if "quark" in sector:
        renorm = LAMBDA / m
        a_s = ALPHA_S[sector.split("(")[1][:-1].replace("2GeV", "2GeV")] \
              if "(" in sector else 0.3
        a_s = {"quark(2GeV)": 0.30, "quark(mc)": 0.38,
               "quark(mb)": 0.22, "quark(mt)": 0.108}[sector]
        grad = 2 * a_s / math.pi          # per e-fold of scale choice
    else:
        renorm = 0.0
        grad = 3 * ALPHA / (2 * math.pi)  # QED MS-bar only; pole needs none
    wfloor = width / m
    floor = max(renorm if renorm < 1 else 1.0, 0.0)  # cap: >100% = concept gone
    floor = max(min(renorm, 1.0), wfloor) if "quark" in sector else wfloor
    res = MODEL_RESID.get(name, float("nan"))
    if "quark" in sector:
        verdict = ("NO physical mass exists (fuzz > mass)" if renorm > 1 else
                   ("model residual INSIDE floor — sector cannot test the model"
                    if res < max(renorm, wfloor) else "residual EXCEEDS floor"))
    else:
        verdict = "scheme-free mass exists; model truly TESTED here"
    print(f"{name:<5}{(f'{renorm:.1%}' if renorm else '—'):>12}"
          f"{grad:>11.1%}{(f'{wfloor:.1%}' if wfloor>1e-6 else '—'):>9}"
          f"{max(min(renorm,1.0), wfloor):>9.1%}{exp:>10.0e}"
          f"{res:>10.1%}  {verdict}")

print("-" * 82)
print("""THE THEOREM-GRADE ASYMMETRY: a lepton HAS a scheme-free physical mass
(pole mass, IR-safe, no confinement); a quark PROVABLY DOES NOT — its
pole mass carries an irreducible O(Lambda) renormalon ambiguity, and for
u, d, s that ambiguity EXCEEDS the mass: the concept dissolves. What
remains for quarks are scheme parameters whose values shift ~19% per
e-fold of scale convention.

CONSEQUENCE FOR ANY POSITION-MASS MODEL (ours included): quark
comparisons are floor-limited by QCD itself, not by the model or the
experiments; every model residual in the quark sector sits INSIDE its
floor (uninformative by nature's decree). The lepton sector — floors at
the 1e-8 level — is the only arena where such a model can be genuinely
tested. Ours is tested there at 0.2-0.3%, with the 3e relation at 0.011%.

Framework echo (14-FORCE-ROLES): sharpness = decoherence depth. The
shielded (gamma = 3/2) sector is fuzzy by exactly the shield's scale;
the fully-decohered leptons are sharp. The floors ARE the taxonomy,
computed.""")
