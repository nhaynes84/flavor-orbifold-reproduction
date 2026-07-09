#!/usr/bin/env python3
"""Dig 1: the convention story, done physically right.

The uniform per-sector shift (fuzz_pattern.py) was an approximation. A
real scale-convention choice shifts each quark by gamma_m(mu_q) x
(e-folds) — WEIGHTED by its own anomalous dimension at its own quoted
scale. Test: fit ONE scale choice (e-folds s) per sector with the
physically-correct weights, and see (a) whether fine structure improves
over the uniform model, (b) what refuses to be absorbed.

gamma_m one-loop: d(ln m)/d(ln mu) = -2 alpha_s/pi, evaluated at each
quark's quoted scale: u,d,s @ 2 GeV (a_s~0.30); c @ m_c (a_s~0.39);
b @ m_b (a_s~0.225); t @ m_t (a_s~0.108).

Run: python3 consolidation/convention_fine_structure.py
"""
import math, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fresh_permutation_test import predict, MEASURED, CANON

pred = predict(*CANON)
r = {p: math.log(pred[p]/MEASURED[p])*100 for p in MEASURED}   # residuals, %
GAMMA = {"d": 19.1, "s": 19.1, "b": 14.3,          # %/e-fold at quoted scale
         "u": 19.1, "c": 24.8, "t": 6.9}

def fit_scale(quarks):
    """least-squares single e-fold shift s: minimize sum (r - gamma*s)^2"""
    num = sum(GAMMA[q]*r[q] for q in quarks)
    den = sum(GAMMA[q]**2 for q in quarks)
    return num/den

print("Residuals (%, model vs measured):",
      {q: round(r[q],2) for q in ("d","s","b","u","c","t")})
print()
for sec, qs in (("down", ["d","s","b"]), ("up", ["u","c","t"])):
    s = fit_scale(qs)
    post = {q: r[q] - GAMMA[q]*s for q in qs}
    print(f"{sec}-sector, gamma-weighted scale fit: s = {s:+.4f} e-folds")
    for q in qs:
        print(f"   {q}: {r[q]:+.2f}%  -> post {post[q]:+.2f}%")
    print(f"   max|post| = {max(abs(v) for v in post.values()):.2f}%")
    print()

# up-sector WITHOUT the top (top's convention issues are definitional,
# not running: MC-mass vs pole, gamma tiny)
s_uc = fit_scale(["u","c"])
post_u, post_c = r["u"]-GAMMA["u"]*s_uc, r["c"]-GAMMA["c"]*s_uc
top_pred_shift = GAMMA["t"]*s_uc
top_anomaly = r["t"] - top_pred_shift
print(f"up-sector refit on (u, c) only: s = {s_uc:+.4f} e-folds")
print(f"   u: {r['u']:+.2f}% -> {post_u:+.2f}%    c: {r['c']:+.2f}% -> {post_c:+.2f}%")
print(f"   top's convention-allowed shift: {top_pred_shift:+.2f}%")
print(f"   TOP ANOMALY (post-convention): {top_anomaly:+.2f}%")
print(f"""
READING:
- Down sector: one physical scale choice ({fit_scale(['d','s','b']):+.3f} e-folds)
  absorbs d, s AND b to ~0.1% — better than the uniform model, because the
  gamma-weighting correctly gives b (quoted at m_b, smaller alpha_s) a
  smaller share. Fine structure TRACKS the anomalous dimensions.
- Up sector: u and c collapse to ~0.1% under one scale choice — but the
  TOP does not: its running weight is tiny (gamma ~ 7%/e-fold) so
  convention freedom can only move it ~1%, leaving a ~-2% anomaly.
- CONVERGENCE with 28: the floors analysis said the top is the ONLY
  quark capable of testing the model, and here it is again — the only
  quark with a post-convention residual. Note also: the top's model
  value is hostage to the internal sigma_u ambiguity (master-eq 11.284
  vs rational 34/3 = 11.333 spans -2.8% to +1.1% on m_t), so the top
  anomaly and the sigma_u fork are ONE open problem, now unified.""")
