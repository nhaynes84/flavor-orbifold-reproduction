#!/usr/bin/env python3
"""M6 session 2: what is epsilon?

LAMBDA PREAMBLE (read first): epsilon = 0.042 +- ~0.002 is one small
number and the grammar hits anything within 1% at ~100% (doc 05). An
honest search declaration: the candidates below are EVERY expression
tried (no hidden scan); all are zero-parameter combinations of already-
frozen framework constants. The fit CANNOT establish any of them —
only the registered future bet can. Search cost: ~8 expressions.

Candidates (in delta(1/alpha) units):
    1/(4pi)      = 0.0796    one-loop unit
    1/(8pi)      = 0.0398    = alpha_b / N_c  (one color's share of bare)
    1/(2pi^2)    = 0.0507
    1/(16pi)     = 0.0199
    c0^2 * 1     = 0.0400    (fluctuation quantum squared, raw)
    c0^3 * ...   = 0.0080
    hbar_orb     = 0.0400    ( = 1/25 = c0^2, same value as above)
    g_b^2/16pi^2 = 0.0095

Note the degeneracy: 1/(8pi) = 0.0398 and c0^2 = 1/25 = 0.0400 are
numerically INDISTINGUISHABLE at current precision (0.5% apart). Two
different stories, one number — flagged, not resolved.

The zero-parameter formula under test (crease-fee + candidate quantum):
    1/alpha_i = (8pi/3) * ell_i  +  (ell_i + 2) * eps0,   eps0 = 1/(8pi)
i.e., every coupling from {7, 3.5, 1} and pi alone.

Run: python3 consolidation/m6_session2.py
"""
import math

ELL = {"U(1)": 7.0, "SU(2)": 3.5, "SU(3)": 1.0}
INV_MEAS = {"U(1)": (59.012, 0.02), "SU(2)": (29.569, 0.02),
            "SU(3)": (8.4818, 0.065)}

def test(eps0, label):
    chi2 = 0.0
    rows = []
    for g in ELL:
        pred = (8*math.pi/3)*ELL[g] + (ELL[g]+2)*eps0
        meas, sig = INV_MEAS[g]
        z = (pred - meas)/sig
        chi2 += z*z
        rows.append((g, pred, meas, z))
    return chi2, rows

print("Zero-parameter formula: 1/alpha_i = (8pi/3) ell_i + (ell_i+2) eps0")
print(f"{'eps0':<26}{'chi2 (3 dof, ZERO params)':>28}")
for eps0, label in ((1/(8*math.pi), "1/(8pi) = alpha_b/3"),
                    (1/25, "c0^2 = hbar_orb = 1/25"),
                    (1/(4*math.pi), "1/(4pi)"),
                    (1/(2*math.pi**2), "1/(2pi^2)")):
    c2, _ = test(eps0, label)
    print(f"  {label:<24}{c2:>16.2f}")

c2, rows = test(1/(8*math.pi), "")
print(f"\nDetail for eps0 = 1/(8pi) (chi2 = {c2:.2f}/3 dof):")
for g, pred, meas, z in rows:
    print(f"  {g:<6} pred {pred:8.4f}   meas {meas:8.4f}   ({z:+.2f} sigma)")

# The registered bet: what alpha_s(M_Z) does the formula demand?
inv3 = (8*math.pi/3)*1 + 3/(8*math.pi)
print(f"\nREGISTERED BET (the only honest arbiter):")
print(f"  formula demands 1/alpha_3 = {inv3:.4f}  ->  alpha_s(M_Z) = {1/inv3:.5f}")
print(f"  current: 0.1179 +- 0.0009 (formula sits {-(0.1179-1/inv3)/0.0009:+.1f} sigma)")
print(f"  FCC-ee projects +-0.0001-0.0002: DECIDES at ~5-10 sigma equivalent")
print(f"""
HONEST STATUS:
  - eps0 = 1/(8pi) and eps0 = c0^2 = 1/25 are numerically degenerate
    (0.5% apart); both give chi2 ~ 2.4-2.5 on THREE dof with ZERO
    parameters. The formula's weak spot: SU(2) at -1.4 sigma.
  - Search cost ~8 expressions: the fit alone is therefore WEAK
    evidence. What is NOT weak: a zero-parameter formula surviving
    3 dof, and a falsifiable alpha_s demand on the books.
  - Two mechanism stories for the same number, unresolved:
    (a) eps0 = alpha_b/N_c — each counting unit leaks one color-share
        of the bare noise (loop-flavored).
    (b) eps0 = c0^2 — each counting unit carries one fluctuation-quantum
        variance (Law-5-flavored: the same c0^2 that screens masses).
    If (b), the SAME constant governs mass screening AND coupling
    conversion — a cross-domain identification the framework could
    eventually test. If (a), it is bookkeeping of the bare coupling.
  - VERDICT PROTOCOL: alpha_s(M_Z) = {1/inv3:.5f} on the books. FCC-ee
    (or a 2x lattice improvement) decides. Until then: CANDIDATE.""")
