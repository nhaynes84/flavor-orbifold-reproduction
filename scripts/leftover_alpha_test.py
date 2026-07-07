#!/usr/bin/env python3
"""Test of the 'alpha is a leftover' hypothesis (Nick, 2026-07-06).

Hypothesis: alpha_em is not a coupling to be derived from channel counting
but the RESIDUE of the boundary->4D projection — 'a dirty exchange rate.'
If so, the framework's systematic residuals (bare boundary value vs dressed
4D measurement) should quantize in units of alpha.

Test 1: residual quantization. For each genuine framework residual r_i,
  compute n_i = r_i / alpha. If alpha is the dressing quantum, n_i should
  cluster near integers. Score vs a uniform null (Monte Carlo).
  Also tested: quanta alpha/pi (the corpus's O(alpha/pi) claim) and the
  framework's own screening quantum C0^3(1-C0^2).

Test 2: the leftover identity candidate. The framework's smallest colored
  screening correction is C0^3(1-C0^2) = 6/781.25 = 1/130.21 — numerically
  inside alpha's running range [1/137.04 at Q=0, 1/128.94 at M_Z]. Where on
  the running curve does 1/alpha = 130.21 sit? (One-loop leptons+quarks.)

Run: python3 consolidation/leftover_alpha_test.py
"""

import math
import random

ALPHA0 = 1 / 137.035999    # alpha(0)
C0 = 0.2

# Genuine, independent framework residuals (pred/meas - 1), from
# verify_consolidation.py output (fitted quantities excluded; one entry
# per independent structure).
RESIDUALS = {
    "V_us bare e^{-3/2}":      -0.0052,
    "sin2thW 2/9 vs on-shell": -0.0044,
    "alpha_s 3/(8pi)":         +0.0124,
    "T_c = 2 Lambda_BW":       +0.0216,
    "r_proton":                +0.0196,
    "m_proton":                -0.0040,
    "M_P route A (3/4)y_t":    +0.0145,
    "DM/baryon 26/5":          -0.0226,
    "m_b boundary formula":    -0.0043,
    "m_tau boundary formula":  +0.0031,
    "Hubble 27/25":            -0.0034,
    "lambda_H (1+C0^2)/8":     +0.0070,
    "n-p splitting":           +0.0094,
}

def quantization_score(residuals, quantum, n_max=4):
    """Mean squared distance of r/quantum to nearest integer in [-n_max,n_max].
    Lower = more quantized. Residuals beyond n_max*quantum are clipped out."""
    ds, used = [], 0
    for r in residuals:
        n = r / quantum
        if abs(n) > n_max + 0.5:
            continue
        used += 1
        ds.append((n - round(n)) ** 2)
    return (sum(ds) / len(ds)) if ds else float("nan"), used

def null_pvalue(observed_score, n_used, quantum, trials=20000, n_max=4):
    """Null: residuals uniform in the same magnitude band."""
    random.seed(11)
    lo, hi = 0.002, (n_max + 0.5) * quantum
    better = 0
    for _ in range(trials):
        rs = [random.choice((-1, 1)) * random.uniform(lo, hi)
              for _ in range(n_used)]
        s, _ = quantization_score(rs, quantum, n_max)
        if s <= observed_score:
            better += 1
    return better / trials

print("TEST 1 — residual quantization")
print("-" * 72)
for label, q in [("alpha = 1/137.036", ALPHA0),
                 ("alpha/pi (corpus O(a/pi) claim)", ALPHA0 / math.pi),
                 ("C0^3(1-C0^2) = 1/130.2 (screening quantum)",
                  C0**3 * (1 - C0**2))]:
    score, used = quantization_score(RESIDUALS.values(), q)
    p = null_pvalue(score, used, q)
    print(f"  quantum {label:<44}")
    print(f"    n_i = " + ", ".join(
        f"{r/q:+.2f}" for r in RESIDUALS.values() if abs(r/q) <= 4.5))
    print(f"    used {used}/13, mean sq dist to integer = {score:.4f}, "
          f"null p = {p:.3f}")
print()

print("TEST 2 — the leftover identity candidate")
print("-" * 72)
q_screen = C0**3 * (1 - C0**2)
print(f"  C0^3(1-C0^2) = {q_screen:.6f} = 1/{1/q_screen:.2f}")
print(f"  alpha range:   1/137.04 (Q=0)  ...  1/128.94 (M_Z)")

# one-loop QED running with lepton + quark (constituent-ish) thresholds
charged = [  # (mass GeV, Q^2 * Nc)
    (0.000511, 1), (0.10566, 1), (1.77686, 1),          # e, mu, tau
    (0.0022, 4/9 * 3), (1.27, 4/9 * 3), (172.5, 4/9 * 3),  # u, c, t
    (0.0047, 1/9 * 3), (0.095, 1/9 * 3), (4.18, 1/9 * 3),  # d, s, b
]
def inv_alpha(Q):
    s = sum(w * math.log(Q / m) for m, w in charged if Q > m)
    return 1 / ALPHA0 - (2 / (3 * math.pi)) * s
for Q in (0.001, 1.0, 2.0, 4.18, 10.0, 30.0, 91.19):
    print(f"    1/alpha({Q:>6.3f} GeV) ~ {inv_alpha(Q):7.2f}"
          + ("   <-- crosses 130.21 near here"
             if abs(inv_alpha(Q) - 1/q_screen) < 0.8 else ""))
print("""
  NOTE: current-quark masses in the one-loop sum understate hadronic
  running; treat the scale placement as +-1 decade. The point is only
  WHETHER 1/130.2 lies inside the physical running band (it does) and
  roughly where. Any claim 'alpha == screening quantum at scale X' needs
  a mechanism for X before it counts (lambda-gate).""")
