#!/usr/bin/env python3
"""M4-partial: the matching decomposition in the physically-correct
variable, with structure discrimination.

One-loop matching is ADDITIVE in 1/alpha (not multiplicative in %).
The general form: delta(1/alpha_i) = ell_i * Delta + [matching terms],
where ell_i*Delta is a renormalization of the ONE bare coupling
(propagating through the volumes) and the matching terms carry group
structure. Three candidate structures, 3 data points:

  (i)   common bare shift only:        delta = ell_i * Delta
  (ii)  bare shift + Casimir term:     delta = ell_i * Delta + k * C_A(i)
  (iii) bare shift + universal matter: delta = ell_i * Delta + c
        (matter loops are ~universal across groups in GUT norm — that's
         why the b_i differences are pure gauge)

Run: python3 consolidation/m4_partial.py
"""
import math

alpha_b = 1.5 / (4 * math.pi)
ELL = {"U(1)": 7.0, "SU(2)": 3.5, "SU(3)": 1.0}
CA  = {"U(1)": 0.0, "SU(2)": 2.0, "SU(3)": 3.0}
INV_MEAS = {"U(1)": 59.012, "SU(2)": 29.569, "SU(3)": 8.4818}
SIG = {"U(1)": 0.02, "SU(2)": 0.02, "SU(3)": 0.065}   # 1/alpha uncertainties

delta = {g: INV_MEAS[g] - ELL[g] / alpha_b for g in ELL}
print("Additive residuals delta(1/alpha) = measured - ell/alpha_b:")
for g in ELL:
    print(f"  {g:<6} {delta[g]:+.3f}  (sigma {SIG[g]})")

def chi2(pred):
    return sum(((delta[g] - pred[g]) / SIG[g]) ** 2 for g in ELL)

# (i) common bare shift only (1 param, weighted LSQ)
num = sum(ELL[g] * delta[g] / SIG[g]**2 for g in ELL)
den = sum(ELL[g]**2 / SIG[g]**2 for g in ELL)
D1 = num / den
c2_i = chi2({g: ELL[g] * D1 for g in ELL})

# (ii) bare + Casimir (2 params, solve weighted LSQ on grid)
best_ii = None
for Dx in [x * 1e-4 for x in range(300, 700)]:
    for kx in [x * 1e-4 for x in range(0, 800)]:
        c2 = chi2({g: ELL[g] * Dx + kx * CA[g] for g in ELL})
        if best_ii is None or c2 < best_ii[0]:
            best_ii = (c2, Dx, kx)
c2_ii, D2, K2 = best_ii

# (iii) bare + universal constant (2 params)
best_iii = None
for Dx in [x * 1e-4 for x in range(200, 700)]:
    for cx in [x * 1e-4 for x in range(0, 1500)]:
        c2 = chi2({g: ELL[g] * Dx + cx for g in ELL})
        if best_iii is None or c2 < best_iii[0]:
            best_iii = (c2, Dx, cx)
c2_iii, D3, C3 = best_iii

print(f"\n(i)   common-bare only:      Delta = {D1:.4f}          chi2 = {c2_i:.1f} (2 dof)")
print(f"(ii)  bare + Casimir:        Delta = {D2:.4f}, k = {K2:.4f}/C_A  chi2 = {c2_ii:.2f} (1 dof)")
print(f"(iii) bare + universal:      Delta = {D3:.4f}, c = {C3:.4f}     chi2 = {c2_iii:.2f} (1 dof)")

print(f"""
DISCRIMINATION:
  (i) is EXCLUDED (chi2 {c2_i:.0f} on 2 dof — the residuals are not a pure
  bare-coupling renormalization).
  (ii) Casimir structure: chi2 {c2_ii:.2f} — consistent. In natural units
  the Casimir coefficient is a = 4*pi*k = {4*math.pi*K2:.2f} per C_A: O(1). The
  bare renormalization Delta = {D2:.4f} means delta(alpha_b)/alpha_b =
  {-D2*alpha_b*100:.2f}% — the ~one-loop-unit shift already noted.
  (iii) universal-matter structure: chi2 {c2_iii:.2f} — {'disfavored vs (ii)' if c2_iii > c2_ii + 1 else 'comparable to (ii)'}.

POSTED M4 TARGETS (correct variable, replaces the %-form):
  the eventual matching computation must produce
    delta(1/alpha_i) = ell_i * Delta + k * C_A(i)
  with Delta = {D2:.4f} +- ~0.005 and k = {K2:.4f} +- ~0.02 (a = {4*math.pi*K2:.1f}/4pi).
  beta_2 in the ring scheme: still blocked on the action's UV definition.""")
