#!/usr/bin/env python3
"""Rigidity tests: Newton-mode scoring of the framework's zero-choice chains.

Complement to look_elsewhere.py. That script scores claims as if every
formula was searched (worst case). This one scores the chains that are
claimed to be FORCED by the frozen axioms, by enumerating the alternatives
the axioms could have taken and asking whether the data uniquely selects
the framework's choice.

Test 1 — the coupling ladder (Axiom A4).
  Claim: g(q) = e^{(q-1)/2} is structural. Given the ladder FAMILY
  (geometric/linear/power) and the fixed structure N_q=7, N_l=3, C0=1/5,
  the three measured sector slopes must be hit simultaneously:
    sigma_d = 7g1 - C0*sqrt(g1)   vs ln(m_b/m_d) = 6.797  +-0.7% (m_d fuzz)
    sigma_u = 7g2 - C0*sqrt(g2)   vs ln(m_t/m_u) = 11.288 +-0.5% (m_u fuzz)
    sigma_l = 3g3                 vs ln(m_tau/m_e) = 8.15399 +-0.05%
  Key fact: for ANY geometric ladder the lepton slope fixes the ratio
  rho = g(q+1)/g(q) at 1.64858 +- ~1e-4. The question becomes: which
  'natural' ratios land there, and do the other two sectors then follow
  with zero further choices?

Test 2 — the screening rule (one rule, many observables).
  Claim: corrections are C0^k (colorless) or C0^k(1-C0^{k-1}) (colored),
  sign physical. One 16-value alphabet must cover every needed correction
  across independent observables. Jointly score: how often would a random
  needed-correction set be covered this well?

Run: python3 consolidation/rigidity_tests.py
"""

import math
import random

E = math.e
C0 = 0.2
N_Q, N_L = 7, 3

# ------------------------------------------------------------------ Test 1
print("TEST 1 — coupling ladder rigidity (Axiom A4)")
print("-" * 78)

T_D, TOL_D = math.log(4180 / 4.67), 0.007        # m_d MS-bar +-5% -> 0.7%
T_U, TOL_U = math.log(172500 / 2.16), 0.005      # m_u +-5% -> 0.5%
T_L, TOL_L = math.log(1776.86 / 0.51100), 0.0005 # leptons exact; 0.05% generous

def slopes(g1, g2, g3):
    return (N_Q * g1 - C0 * math.sqrt(g1),
            N_Q * g2 - C0 * math.sqrt(g2),
            N_L * g3)

def passes(g1, g2, g3):
    sd, su, sl = slopes(g1, g2, g3)
    return (abs(sd / T_D - 1) < TOL_D and abs(su / T_U - 1) < TOL_U
            and abs(sl / T_L - 1) < TOL_L)

# The lepton slope alone pins the geometric ratio:
rho_required = math.sqrt(T_L / 3)
print(f"  geometric ladders: lepton slope forces ratio rho = {rho_required:.5f}"
      f" +- {TOL_L/2:.2%}")
print(f"  sqrt(e) = {math.sqrt(E):.5f}  -> off by "
      f"{(math.sqrt(E)/rho_required - 1):+.4%}")

# Candidate 'natural' ratios an agent could postulate for a geometric ladder
NATURAL_RATIOS = {
    "sqrt(e)": math.sqrt(E), "5/3": 5 / 3, "phi": (1 + math.sqrt(5)) / 2,
    "fine 33/20": 33 / 20, "8/5": 8 / 5, "13/8": 13 / 8, "18/11": 18 / 11,
    "sqrt(8/3)": math.sqrt(8 / 3), "sqrt(e)*(1+C0^3)": math.sqrt(E) * 1.008,
    "e/phi": E / ((1 + math.sqrt(5)) / 2), "pi/2": math.pi / 2,
    "sqrt(3)": math.sqrt(3), "5/3-1/50": 5/3 - 0.02, "1.65": 1.65,
    "cbrt(9/2)": (9 / 2) ** (1 / 3), "sqrt(pi)": math.sqrt(math.pi),
    "2^(5/7)": 2 ** (5 / 7), "e^(1/2)": E ** 0.5,
}
geo_pass = [n for n, r in NATURAL_RATIOS.items()
            if passes(1.0, r, r * r)]
print(f"  geometric ratios passing ALL THREE sectors: {geo_pass}"
      f"  ({len(geo_pass)}/{len(NATURAL_RATIOS)} tried)")

# Linear ladders g(q) = 1 + a(q-1): lepton slope forces a; check the rest
a = (T_L / 3 - 1) / 2
lin_ok = passes(1.0, 1 + a, 1 + 2 * a)
sd, su, sl = slopes(1.0, 1 + a, 1 + 2 * a)
print(f"  linear ladder (a forced by lepton = {a:.4f}): sigma_u -> {su:.3f}"
      f" vs {T_U:.3f} ({su/T_U-1:+.1%})  pass={lin_ok}")

# Power ladders g(q) = q^b
b = math.log(T_L / 3) / math.log(3)
pw_ok = passes(1.0, 2 ** b, 3 ** b)
sd, su, sl = slopes(1.0, 2 ** b, 3 ** b)
print(f"  power ladder  (b forced by lepton = {b:.4f}): sigma_u -> {su:.3f}"
      f" vs {T_U:.3f} ({su/T_U-1:+.1%})  pass={pw_ok}")
print("""  READING: the lepton sector pins the ladder ratio to 1e-4; only
  sqrt(e)-class ratios survive, and geometric form is then CONFIRMED by the
  up sector with zero further choices (linear/power fail by 7-13%). A4 is
  Newton-mode: structure fixed, two sectors predicted, both land.""")
print()

# ------------------------------------------------------------------ Test 2
print("TEST 2 — screening-rule joint rigidity (one alphabet, many targets)")
print("-" * 78)

# the rule's full alphabet: {1 +- C0^k, 1 +- C0^k(1-C0^{k-1})} for k=1..3
ALPHABET = sorted(set(
    [1 + s * C0 ** k for k in (1, 2, 3) for s in (+1, -1)] +
    [1 + s * C0 ** k * (1 - C0 ** (k - 1)) for k in (2, 3) for s in (+1, -1)]
))

# needed correction factors, from measurement / (bare framework value):
# (name, needed_factor, honest rel tolerance)
NEEDED = [
    ("y_t: m_t vs v/sqrt2 (needs 1-C0^3)", 0.9923, 0.003),
    ("M0:  eff/(4/H)      (needs 1-C0^2)", 0.154 / 0.16, 0.003),
    ("m_t: vs unscreened  (needs 1-C0^2)", 0.96, 0.01),
    ("proton: m_p/(Nq*Lambda_BW)", 938.272 / (N_Q * 137.909), 0.005),
    ("n-p: 2*1.293/(m_d-m_u)", 2 * 1.293 / (4.67 - 2.16), 0.02),
    ("lambda_H: 8*lam vs 1 (needs 1+C0^2)", 8 * 0.1291 / 1.0, 0.01),
]

def best_hit(needed, tol):
    cands = [(abs(v / needed - 1), v) for v in ALPHABET]
    d, v = min(cands)
    return (d <= tol, d, v)

n_hit = 0
for name, needed, tol in NEEDED:
    ok, d, v = best_hit(needed, tol)
    n_hit += ok
    print(f"  {name:<38} needed={needed:.4f}  nearest={v:.4f} "
          f"(off {d:+.2%}, tol {tol:.1%})  {'HIT' if ok else 'miss'}")

# Null: random needed-corrections with same magnitudes-distribution.
# Draw |correction-1| log-uniform in [0.001, 0.3], random sign; same tols.
random.seed(7)
trials, as_good = 20000, 0
for _ in range(trials):
    hits = 0
    for _, _, tol in NEEDED:
        mag = math.exp(random.uniform(math.log(0.001), math.log(0.3)))
        needed = 1 + random.choice((-1, 1)) * mag
        if best_hit(needed, tol)[0]:
            hits += 1
    if hits >= n_hit:
        as_good += 1
p_joint = as_good / trials
print(f"\n  observed: {n_hit}/{len(NEEDED)} covered by the 10-value alphabet")
print(f"  null (random corrections, same tolerances): P(>= {n_hit} hits) = {p_joint:.4f}")
sig = "genuine joint structure" if p_joint < 0.01 else (
      "suggestive" if p_joint < 0.1 else "not significant")
print(f"  -> {sig}")
