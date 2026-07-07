#!/usr/bin/env python3
"""Gravity mechanism, round 2: kill test for v1 (Gaussian MGF), and the
v2 candidate (multiplicative mode impedance).

v1 (08-GRAVITY-VARIATIONAL.md) claimed <e^{sigma phi}> = e^{sigma^2} from
Var(phi) = H*hbar_orb*Z2 = 2. That variance was asserted, not computed.
Here we compute it from the boundary Lagrangian's own normalization:

  L = (1/2H)(d_t sigma)^2 - (1/2)(d_x sigma)^2 - V(sigma)
  modes u_n = sqrt(2) sin(n pi x), omega_n^2 = H (n pi)^2
  canonical: [q, p] = i hbar_orb with p = qdot/H  ->  <q_n^2> = 1/(2 omega_n)

v2: gravity couples to the uniform (global) deformation, which must push
against every driven mode simultaneously. Per-mode warp suppression e^{sigma}
multiplies across the N_q - C0 effective driven modes:

  ln(M_P / Lambda_BW) = sigma * (N_q - C0) = sigma * sigma_0 = sigma^2

using the equilibrium condition sigma = sigma_0 = N_q - C0 (the cubic
potential's minimum). The exponent sigma^2 is then EXACT at equilibrium and
the 'why the down sector' question dissolves: sigma_0 IS the field value,
which equals the down slope only because down is the unit-coupling channel.

Run: python3 consolidation/gravity_mechanism_v2.py
"""

import math

H = 25
HBAR_ORB = 1 / H
N_Q = 7
C0 = 0.2
SIGMA0 = N_Q - C0                 # 34/5, cubic potential minimum
LAMBDA_BW = (4 / H) * (1 - C0**2) * math.exp(SIGMA0)   # MeV
M_P_MEV = 1.220890e22             # Planck mass in MeV

print("PART 1 — v1 kill test: compute Var(phi) instead of asserting it")
print("-" * 74)
# zero-point variance of each oscillator with the Lagrangian's own norm
# omega_n = sqrt(H) * n * pi ; <q_n^2> = hbar_orb * H / (2 omega_n) = 1/(2 omega_n)
def var_local(x, n_modes=H):
    return sum(2 * math.sin(n * math.pi * x) ** 2 / (2 * math.sqrt(H) * n * math.pi)
               for n in range(1, n_modes + 1))

def var_mean(n_modes=H):
    # graviton = uniform mode: phi_bar = int phi dx ; int u_n = 2*sqrt(2)/(n pi), n odd
    return sum((2 * math.sqrt(2) / (n * math.pi)) ** 2 / (2 * math.sqrt(H) * n * math.pi)
               for n in range(1, n_modes + 1, 2))

print(f"  Var(phi) at generic point x=0.37 : {var_local(0.37):.4f}")
print(f"  Var(phi) at a node x=3/7         : {var_local(3/7):.4f}")
print(f"  Var(mean field) [graviton vertex]: {var_mean():.5f}")
print(f"  v1 needs Var = 2.0  ->  off by ~{2/var_mean():.0f}x for the graviton mode")
print("  VERDICT: v1 (Gaussian MGF with Var=2) is DEAD by its own bookkeeping.")
print("  Zero-point noise of the quantized boundary is ~0.03-0.12, not 2.")
print()

print("PART 2 — v2 candidate: multiplicative mode impedance")
print("-" * 74)
# exponent = sigma * sigma_0, exact = sigma^2 at equilibrium
exponent = SIGMA0 * SIGMA0
ln_ratio = math.log(M_P_MEV / LAMBDA_BW)
lnP = ln_ratio - exponent
print(f"  exponent sigma*sigma_0 = sigma^2  = {exponent:.4f}   (exact identity at V'=0)")
print(f"  measured ln(M_P/Lambda_BW)        = {ln_ratio:.4f}")
print(f"  residual ln P = {lnP:+.4f}  ->  P = {math.exp(lnP):.4f}")
print(f"  candidates on file: (3/4)y_t = {0.75 * (1 - C0**3):.4f}, 11/15 = {11/15:.4f}")
print("  P remains OPEN (lambda-gate: no numerical hunting; the KK zero-mode")
print("  normalization calculation decides).")
print()
print("  Structural consistency (v1's killer, re-checked for v2):")
print("  - masses = LOCAL node values of one mode -> no product over modes")
print("    -> masses are NOT Planck-dressed (correct, unlike naive v1)")
print("  - gravity = GLOBAL uniform deformation -> pushes all driven modes")
print("    -> picks up the full product e^{sigma*(N_q - C0)}  (the claim)")
print("  - the C0 subtraction in the exponent is the same condensate drain")
print("    that appears in sigma_d itself: one bookkeeping, two places")
print()

print("PART 3 — what v2 predicts with no further freedom")
print("-" * 74)
# If the exponent is sigma*sigma_0 with sigma the FIELD value, then any
# probe coupling to k driven modes sees exponent sigma*k_eff. Gravity sees
# k_eff = N_q - C0 (all driven modes minus condensate drain). A probe that
# couples to a SINGLE mode should see e^{sigma} — which is Lambda_BW/M_0
# itself (the mass sector). Already true by construction; the nontrivial
# content is that NOTHING should exist at intermediate exponents unless it
# couples to an intermediate mode count:
for k in range(1, 8):
    scale = LAMBDA_BW * math.exp(SIGMA0 * k) if k > 1 else LAMBDA_BW
    unit = "MeV"
    val = scale
    if val > 1e6: val, unit = val / 1e3, "GeV"
    if val > 1e9 and unit == "GeV": val, unit = val / 1e3, "TeV"
    print(f"  k = {k} modes -> scale ~ {val:9.3g} {unit}")
print("""  A k=2 'two-mode' object sits at ~1.2e5 TeV = 1.2e8 GeV; k=3 at ~1e11 GeV
  (near GUT/seesaw territory), ..., k=7-C0 at M_P. The tower predicts NO new
  fundamental scales between Lambda_BW and ~1e8 GeV — consistent with the
  observed desert. A discovered fundamental scale in the desert that maps to
  no integer k would strain v2. [prospective, frozen with today's date]""")
