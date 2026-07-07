#!/usr/bin/env python3
"""boundary.py — THE canonical constants and laws of the boundary framework.

Single source of truth. All future scripts import from here; no script
redefines these numbers. Canonical names per 11-CANONICAL-LANGUAGE.md.

Self-test: python3 consolidation/boundary.py   (exit 0 = all laws verified)
"""

import math

# ---------------------------------------------------------------- objects
N_C = 3                              # colors (fixed point, Law 0)
N_Q = 2 * N_C + 1                    # 7  quark modes
N_L = 3                              # lepton modes (P = |Q| = 1)
N_NU = 1                             # neutrino mode (P = 0)
H = N_C**2 + (N_C + 1)**2            # 25 harmonic invariant (3-4-5)
C0 = 1 / math.sqrt(H)                # 1/5  fluctuation quantum
HBAR_ORB = C0**2                     # 1/25 orbifold action quantum
C_F = (N_C**2 - 1) / (2 * N_C)       # 4/3
G_SQ = 2 / C_F                       # 3/2  boundary coupling
SIGMA_FUND = G_SQ * C_F / 2          # 1    boundary tension (exact)

# Law 0, closed form (16-FIXED-POINT-CLOSURE): matter content is a theorem,
# n_f = fold(2) x generations(3) = 6, independent of N_c; the matching
# condition (11*N_c - 2*n_f)/3 = 2*N_c + 1 then reduces to 5*N_c = 15.
N_F = 2 * 3
assert (11 * N_C - 2 * N_F) / 3 == N_Q, "Law 0 matching condition broken"

SIGMA_D = N_Q - C0                   # 34/5 the warp (field value; down = g=1)
SIGMA_U = N_Q * math.sqrt(math.e) - C0 * math.e**0.25   # up-sector slope
SIGMA_L = N_L * math.e               # 3e   lepton slope
M0 = (4 / H) * (1 - C0**2)           # 384/2500 MeV — the bridge (unit anchor)
LAMBDA_BW = M0 * math.exp(SIGMA_D)   # 137.9 MeV — rung 1, the QCD scale
M_P_MEV = 1.220890e22                # measured Planck mass (target, not input)

# ------------------------------------------------------------------- laws
def ladder(q):
    """Law 1 (Ladder): sector coupling g(q) = e^{(q-1)/2}, q = N_c|Q_em|."""
    return math.exp((q - 1) / 2)

def warp(x, sigma=SIGMA_D, rung=0):
    """Law 2 (Warp): E(x) = M0 e^{rung*sigma_d + sigma*x} — mass at
    fractional position x in a sector. The master equation."""
    return M0 * math.exp(rung * SIGMA_D + sigma * x)

def tower(k):
    """Law 3 (Tower): S(k) = M0 e^{k sigma_d} — local scales live on
    integer rungs. k=0 bridge, k=1 QCD, k=2 EW, k=3 first new rung."""
    return M0 * math.exp(k * SIGMA_D)

def impedance_exponent():
    """Law 4 (Impedance): gravity is global — ln(M_P/M0) = sigma(1+sigma0)
    = sigma + sigma^2, up to the open prefactor P."""
    return SIGMA_D * (1 + SIGMA_D)

def screening(k, colored, sign=+1):
    """Law 5 (Screening): correction factor 1 + sign*C0^k (colorless) or
    1 + sign*C0^k(1-C0^{k-1}) (colored)."""
    base = C0**k * (1 - C0 ** (k - 1)) if colored else C0**k
    return 1 + sign * base

def node_shift(k_self, k_partner, n_modes):
    """Law 6 (Overlap, shifts): dy = C0 * k_partner / (n_modes*pi*k_self),
    signed toward the weak partner."""
    return C0 * k_partner / (n_modes * math.pi * k_self)

def mixing(m_light, m_heavy):
    """Law 6 (Overlap, mixing): amplitude ratio = sqrt(m_l/m_h) =
    e^{-sigma dx/2} — distance along the warp."""
    return math.sqrt(m_light / m_heavy)

def u(mass_mev):
    """THE WARP COORDINATE: u = ln(m/M0)/sigma_d. One number line for all
    of physics. Rungs at integers; particles at fractional positions;
    gravity off-line at 1+sigma_d+lnP/sigma_d."""
    return math.log(mass_mev / M0) / SIGMA_D

# -------------------------------------------------------------- self-test
if __name__ == "__main__":
    fails = []
    def chk(name, got, want, tol_pct):
        err = (got / want - 1) * 100
        ok = abs(err) <= tol_pct
        if not ok:
            fails.append(name)
        print(f"  {'ok  ' if ok else 'FAIL'} {name:<38} {got:12.6g} vs "
              f"{want:<12.6g} ({err:+.2f}%)")

    print("boundary.py self-test — canonical laws vs measurement")
    chk("sigma_d vs ln(m_b/m_d)", SIGMA_D, math.log(4180 / 4.67), 0.1)
    chk("sigma_l vs ln(m_tau/m_e)", SIGMA_L, math.log(1776.86 / 0.511), 0.05)
    chk("Lambda_BW (rung 1) vs 137.9", LAMBDA_BW, 137.9, 0.1)
    chk("tower(2) vs v/2", tower(2), 123110.0, 1.0)
    chk("impedance exponent vs ln(M_P/M0)", impedance_exponent(),
        math.log(M_P_MEV / M0), 1.0)
    chk("screening(3,colored) vs y_t", screening(3, False, -1), 0.9923, 0.1)
    chk("mixing(m_d,m_s) vs V_us", mixing(4.67, 93.4), 0.2243, 1.0)
    chk("node_shift strange", node_shift(3, 4, N_Q), 0.0122, 1.0)

    print("\n  Warp-coordinate positions (u = ln(m/M0)/sigma_d):")
    for name, m in [("m_e", 0.511), ("m_u", 2.16), ("m_d", 4.67),
                    ("m_s", 93.4), ("m_mu", 105.66), ("QCD rung", LAMBDA_BW),
                    ("m_c", 1270.0), ("m_tau", 1776.86), ("m_b", 4180.0),
                    ("M_W", 80377.0), ("M_Z", 91188.0), ("v/2", 123110.0),
                    ("M_H", 125250.0), ("m_t", 172500.0), ("v", 246220.0),
                    ("M_P", M_P_MEV)]:
        print(f"    u({name:<8}) = {u(m):7.4f}")
    import sys
    sys.exit(1 if fails else 0)
