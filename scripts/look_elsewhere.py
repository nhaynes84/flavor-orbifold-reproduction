#!/usr/bin/env python3
"""Global look-elsewhere analysis for the boundary framework.

Question: given the formula grammar the corpus actually uses (small integers,
rationals, e, pi, roots, C0 powers, products of a few atoms), how surprising
is each headline match REALLY, after correcting for the searchable space?

Method:
  1. Enumerate the grammar: products of <= 3 atoms with exponents in
     {+-1, +-1/2, +-2}, complexity-costed. Dedupe numerically.
  2. Reachability: for random targets, P(some expression within eps).
  3. Per-claim: expected number of chance hits at the claimed precision
     (Poisson lambda). lambda >= 1 -> the match is grammar noise.
     lambda << 1 -> genuinely surprising, survives correction.
  4. Honest tolerances: use measurement/input uncertainty, not the
     claimed distance to a central value.

Run: python3 consolidation/look_elsewhere.py
"""

import math
import bisect
import random

E, PI = math.e, math.pi

# ----------------------------------------------------------------- grammar
# Atoms observed in actual corpus formulas. (value, cost)
ATOMS = {}
for n in range(1, 13):
    ATOMS[str(n)] = (float(n), 1 if n <= 9 else 2)
for name, v in {"25": 25.0, "26": 26.0, "27": 27.0, "13": 13.0, "15": 15.0}.items():
    ATOMS[name] = (v, 2)          # framework integers (H, beta1, ...)
ATOMS.update({
    "e": (E, 1), "pi": (PI, 1),
    "sqrt2": (math.sqrt(2), 1), "sqrt3": (math.sqrt(3), 1),
    "sqrt5": (math.sqrt(5), 2), "sqrte": (math.sqrt(E), 2),
})
EXPS = [1.0, -1.0, 0.5, -0.5, 2.0, -2.0]
EXP_COST = {1.0: 0, -1.0: 0, 0.5: 1, -0.5: 1, 2.0: 1, -2.0: 1}
MAX_COST = 5           # e.g. "3*e" costs 2; "C0*K/(N*pi*k)" ~ 5

def build_grammar():
    """All values expressible as products of <=3 (atom^exp) with cost<=MAX_COST."""
    singles = []
    for name, (v, c) in ATOMS.items():
        for x in EXPS:
            cost = c + EXP_COST[x]
            if cost <= MAX_COST:
                singles.append((v**x, cost))
    values = {}
    def add(val, cost):
        if 1e-4 < val < 1e4:
            key = round(math.log(val), 7)
            if key not in values or cost < values[key][1]:
                values[key] = (val, cost)
    for v, c in singles:
        add(v, c)
    n = len(singles)
    for i in range(n):
        v1, c1 = singles[i]
        for j in range(i, n):
            c2 = c1 + singles[j][1]
            if c2 > MAX_COST:
                continue
            v2 = v1 * singles[j][0]
            add(v2, c2)
            for k in range(j, n):
                c3 = c2 + singles[k][1]
                if c3 <= MAX_COST:
                    add(v2 * singles[k][0], c3)
    out = sorted((v, c) for v, c in values.values())
    return out

GRAMMAR = build_grammar()
G_VALS = [v for v, _ in GRAMMAR]
G_LOGS = [math.log(v) for v in G_VALS]

def hits_within(target, rel_tol, max_cost=MAX_COST):
    """Grammar values within rel_tol of target (optionally cost-limited)."""
    lo = bisect.bisect_left(G_LOGS, math.log(target * (1 - rel_tol)))
    hi = bisect.bisect_right(G_LOGS, math.log(target * (1 + rel_tol)))
    return [(v, c) for v, c in GRAMMAR[lo:hi] if c <= max_cost]

def local_density(target, max_cost=MAX_COST, window=0.30):
    """Expected chance hits per unit rel_tol near target: lambda(eps) = rate*2eps."""
    lo = bisect.bisect_left(G_LOGS, math.log(target) - window)
    hi = bisect.bisect_right(G_LOGS, math.log(target) + window)
    n = sum(1 for v, c in GRAMMAR[lo:hi] if c <= max_cost)
    return n / (2 * window)      # values per unit log(~rel) interval

def expected_hits(target, rel_tol, max_cost=MAX_COST):
    return local_density(target, max_cost) * 2 * rel_tol

# ------------------------------------------------------------ reachability
print(f"GRAMMAR SIZE: {len(GRAMMAR)} distinct values (cost <= {MAX_COST})")
random.seed(42)
for tol in (0.01, 0.001, 0.0001):
    hit = sum(
        1 for _ in range(2000)
        if hits_within(math.exp(random.uniform(math.log(0.02), math.log(50))), tol)
    )
    print(f"  P(random target in [0.02,50] hit within {tol:.2%}): {hit/2000:.1%}")
print()

# ------------------------------------------------------------- the claims
print(f"{'CLAIM':<44} {'tol used':>9} {'lam=E[chance hits]':>19}  verdict")
print("-" * 100)

results = []

def claim(name, target, tol, max_cost=MAX_COST, note=""):
    lam = expected_hits(target, tol, max_cost)
    n_actual = len(hits_within(target, tol, max_cost))
    if lam < 0.05:
        verdict = "SURVIVES"
    elif lam < 0.5:
        verdict = "marginal"
    else:
        verdict = "grammar noise"
    results.append((name, lam, verdict))
    print(f"{name:<44} {tol:>9.4%} {lam:>13.2f} ({n_actual} found)  {verdict}  {note}")
    return lam

# --- Tier candidates, honest tolerances -----------------------------------
# sigma_l = 3e vs ln(m_tau/m_e): both masses known to <1e-4 rel. The claim
# is an O(1) closed form; honest tol = actual deviation 0.011%.
claim("sigma_l = 3e vs ln(m_tau/m_e)", 8.15399, 0.00011, max_cost=2,
      note="[cost<=2 like '3e']")
claim("  same, full grammar cost<=5", 8.15399, 0.00011)

# sigma_d = 34/5 vs ln(m_b/m_d): m_d MS-bar +-5%, m_b +-0.6% -> ln-ratio
# uncertainty ~0.05 abs ~0.7% rel. Honest tol 0.7%, not 0.04%.
claim("sigma_d = 34/5 vs ln(m_b/m_d) [input-fuzz]", 6.7974, 0.007)

# delta_CP = arctan(m_d/m_u): m_u,m_d each +-5% -> ratio +-7% -> angle +-2.5%.
claim("delta_CP = arctan(m_d/m_u) [input-fuzz]", math.tan(math.radians(65.5)), 0.07)

# alpha_s = 3/(8pi) vs 0.1179 +- 0.0009 (0.8%)
claim("alpha_s = 3/(8pi) vs measured", 0.1179, 0.013, note="pred is 1.24% off; tol=dev")

# Hubble 27/25 vs 73.04/67.4 (+-~1.5% combined)
claim("Hubble ratio = 27/25", 73.04 / 67.4, 0.015)

# M_W/M_Z = sqrt(7/9): measured to 0.01%; pred lands 0.06% off. tol=dev.
claim("M_W/M_Z = sqrt(7/9)", 0.88136, 0.0007)

# y_t = 1 - C0^3 = 0.992 vs 0.9923 (m_t,v known ~0.3%). Near-1 ratios:
# grammar is DENSE near 1 (every p/q). Use dev 0.03%... measure honestly at 0.3%.
claim("y_t = 1-C0^3 [near 1: dense zone]", 0.9923, 0.003)

# M_P prefactor 11/15 at 0.001% -- the suspicious one
claim("M_P prefactor = 11/15 at observed prec", 0.733325, 0.00002,
      note="TOO exact: fitted")
claim("M_P prefactor honest window (Lam_BW 1%)", 0.733325, 0.01)

# 1/alpha_em = 132 vs 128.9 (+2.4% miss dressed as agreement)
claim("1/alpha(M_Z)=132 vs 128.9", 128.9, 0.025, note="only 'matches' at 2.5% tol")

# proton: m_p = 7*Lam*(1-C0^2(1-C0)) -- prefactor freedom test:
# needed correction factor m_p/(7*Lam) = 0.9719, claimed form hits 0.968.
claim("proton screening factor 0.9719", 0.9719, 0.004, note="near-1 dense zone")

print()

# --------------------------------------------- PMNS coefficient-triple test
# Model: s13 = c1*dx, s12 = 1/3 + c2*dx, s23 = 1/2 + c3*dx, dx = 0.0128 (rigid).
# Alphabet: plausible coefficients an agent would try.
alphabet = set()
for v in [1, 2, 3, 4, 5, 6, 1/2, 3/2, 1/3, 2/3, 4/3, 5/3, 1/4, 3/4, 5/2,
          math.sqrt(2), math.sqrt(3), 2*math.sqrt(3), math.sqrt(3)/2,
          1/math.sqrt(3), 2/math.sqrt(3), 1+math.sqrt(2), 1+math.sqrt(3),
          2+math.sqrt(3), 2-math.sqrt(3), math.sqrt(2)/2, PI/2, PI/3, PI/4,
          PI/6, E/2, E, 2*E/3]:
    alphabet.add(v); alphabet.add(-v)
alphabet = sorted(alphabet)
dx = 0.0128011
# NuFIT-era 1-sigma windows
w13 = (0.02240 - 0.00065, 0.02240 + 0.00065)
w12 = (0.307 - 0.012, 0.307 + 0.012)
w23 = (0.546 - 0.021, 0.546 + 0.021)
p1 = sum(1 for c in alphabet if w13[0] <= c * dx <= w13[1]) / len(alphabet)
p2 = sum(1 for c in alphabet if w12[0] <= 1/3 + c * dx <= w12[1]) / len(alphabet)
p3 = sum(1 for c in alphabet if w23[0] <= 0.5 + c * dx <= w23[1]) / len(alphabet)
print("PMNS coefficient-triple test (alphabet of %d plausible coefficients):" % len(alphabet))
print(f"  P(random c1 puts s13 in 1-sigma window) = {p1:.3f}")
print(f"  P(random c2 puts s12 in 1-sigma window) = {p2:.3f}")
print(f"  P(random c3 puts s23 in 1-sigma window) = {p3:.3f}")
def p_to_sigma(p):
    """two-sided p -> z, by bisection on erfc."""
    lo, hi = 0.0, 10.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

pj = p1 * p2 * p3
print(f"  JOINT (3 independent searched picks)   = {pj:.5f}  (~{p_to_sigma(pj):.1f} sigma)")
print(f"  If {{sqrt3,-2,2+sqrt3}} = {{cot,-csc,cot-half}}(pi/6) was fixed FIRST (one pick),")
print(f"  it is stronger; corpus history says c3 came last -> honest range ~{p_to_sigma(p3):.1f}-{p_to_sigma(pj):.1f} sigma.")
print()

# ------------------------------------------------- deep-rosetta demotion
# knowledge/06: 'systematic search of 10,000+ ratios' produced e.g.
# V_cb2 at 0.04%, alpha_s=sqrt(m_c/M_Z) at 0.01%. Expected best hit from
# N trials at uniform-ish density: E[hits within eps] = N * 2eps / spread.
for label, n_trials, prec in [("V_cb2 = sqrt(m_s m_b/(m_c m_t)) @0.04%", 10000, 0.0004),
                              ("alpha_s = sqrt(m_c/M_Z) @0.01%", 10000, 0.0001),
                              ("M_H = sqrt(m_t M_Z) @0.21%", 10000, 0.0021)]:
    lam = n_trials * 2 * prec / 2.0   # ratios spread over ~2 e-folds near target
    print(f"  Rosetta-probe: {label:<44} E[chance hits in 10k search] = {lam:.1f}")
print()

# ----------------------------------------------------------------- summary
print("=" * 100)
print("TIER SUMMARY (lam = expected chance hits at stated tolerance)")
for name, lam, verdict in results:
    print(f"  [{verdict:^13}] lam={lam:7.2f}  {name}")
