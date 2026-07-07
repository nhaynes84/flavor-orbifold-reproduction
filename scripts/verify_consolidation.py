#!/usr/bin/env python3
"""Numeric verification of every quantitative claim in the consolidation docs.

Run:  python3 consolidation/verify_consolidation.py
Exit code 0 = every claim in 01-03 docs reproduces at the stated accuracy.

This is the anti-backtracking gate: if a consolidation doc states a number,
it must be checked here. Companion to the full TDD suite (tdd/, 901 tests).
"""

import math
import sys

E = math.e
PI = math.pi

checks = []  # (name, predicted, measured, stated_err_pct)


def check(name, predicted, measured, max_err_pct):
    err = (predicted - measured) / measured * 100
    ok = abs(err) <= max_err_pct
    checks.append((name, predicted, measured, err, max_err_pct, ok))
    return ok


# ---------------------------------------------------------------- axioms
N_C = 3
N_Q = 2 * N_C + 1                     # 7  quark modes
N_L = 3                               # lepton modes (P=|Q|=1)
H = N_C**2 + (N_C + 1)**2             # 25 Pythagorean invariant
C0 = 1 / math.sqrt(H)                 # 1/5 condensate quantum
CF = (N_C**2 - 1) / (2 * N_C)         # 4/3
G2 = 3 / 2                            # boundary coupling g^2 = 2/CF

assert N_Q == 7 and H == 25 and abs(C0 - 0.2) < 1e-15
# fixed point: beta0 = (11*Nc - 2*nf)/3 with nf = 2*Nc equals N = 2*Nc+1 only at Nc=3
for nc in range(1, 50):
    beta0 = (11 * nc - 2 * (2 * nc)) / 3
    if beta0 == 2 * nc + 1:
        assert nc == 3, f"fixed point not unique: {nc}"

# ------------------------------------------------------- warp factors (sector slopes)
# canonical master-equation form: sigma(q) = N*g - C0*sqrt(g)*Theta_color, g = e^{(q-1)/2}
g_d, g_u, g_l = 1.0, math.sqrt(E), E
SIG_D = N_Q * g_d - C0 * math.sqrt(g_d)      # 34/5 = 6.8 exactly
SIG_U = N_Q * g_u - C0 * math.sqrt(g_u)      # 11.284 (master-equation form)
SIG_L = N_L * g_l                            # 3e = 8.1548 (no color term)
# KNOWN AMBIGUITY (see 04-GAPS): the corpus also uses the rational IR value
# sigma_u = (5/3)*sigma_d = 34/3 = 11.333 (span ratio). The boundary-mass
# formula for m_t works best with the IR value; slopes match data best with
# the master-equation value. Both are used below, labeled.
SIG_U_IR = (5 / 3) * SIG_D                   # 34/3

assert abs(SIG_D - 34 / 5) < 1e-12
check("sigma_d = 34/5 vs ln(m_b/m_d)", SIG_D, math.log(4180 / 4.67), 0.1)
check("sigma_u vs ln(m_t/m_u)", SIG_U, math.log(172500 / 2.16), 0.1)
check("sigma_l = 3e vs ln(m_tau/m_e)", SIG_L, math.log(1776.86 / 0.51100), 0.05)
# IR span ratios (approximate, NOT the exact master-eq values)
check("span sigma_u/sigma_d ~ 5/3", SIG_U / SIG_D, 5 / 3, 0.5)
check("span sigma_l/sigma_d ~ 6/5", SIG_L / SIG_D, 6 / 5, 0.1)

# ------------------------------------------------------------- mass scale + curve
M0 = (4 / H) * (1 - C0**2)                   # 384/2500 = 0.1536 MeV
check("M0 = (4/H)(1-C0^2) vs effective 0.154 MeV", M0, 0.154, 0.5)

LAMBDA_BW = M0 * math.exp(SIG_D)             # top of the down-sector curve
check("Lambda_BW = M0*e^{sigma_d} ~ 137.9 MeV", LAMBDA_BW, 137.9, 0.1)
check("T_c = 2*Lambda_BW vs lattice 270 MeV", 2 * LAMBDA_BW, 270.0, 2.5)

# gen-3 boundary masses  m3 = M0 * N^{7/3} * exp(5*sigma/6)
# NOTE: m_t needs the IR sigma_u (34/3); the master-eq value undershoots by 2.8%.
m_b = M0 * N_Q ** (7 / 3) * math.exp(5 * SIG_D / 6)
m_t = M0 * N_Q ** (7 / 3) * math.exp(5 * SIG_U_IR / 6) * (1 - C0**2)  # top screening
m_tau = M0 * N_L ** (7 / 3) * math.exp(5 * SIG_L / 6)
check("m_b", m_b, 4180.0, 1.0)
check("m_t (IR sigma_u)", m_t, 172500.0, 1.5)
check("m_tau", m_tau, 1776.86, 1.0)

# gen-1: slide down the curve (up sector anchored on measured m_t, as in
# the 4-input scorecard; master-eq sigma_u matches the measured slope)
M_T_MEAS = 172500.0
m_d = m_b * math.exp(-SIG_D)
m_u = M_T_MEAS * math.exp(-SIG_U)
m_e = m_tau * math.exp(-SIG_L)
check("m_d", m_d, 4.67, 1.0)
check("m_u (from measured m_t)", m_u, 2.16, 1.0)
check("m_e", m_e, 0.51100, 1.0)

# gen-2: interior node k/N with node shift toward the weak partner
# canonical shift: delta_y = C0 * k_partner / (N * pi * k_self), sign toward partner
dy_s = C0 * 4 / (N_Q * PI * 3)               # strange, partner charm (k=4)
dy_c = C0 * 3 / (N_Q * PI * 4)               # charm, partner strange (k=3)
dy_mu = C0 * (6 / 5) / (N_L * PI * 2)        # muon, K = span ratio (partner not on same wave)
check("gen2 shift strange +4/(105pi)", dy_s, 0.0122, 1.0)
check("gen2 shift charm -3/(140pi)", dy_c, 0.0065, 5.0)
check("gen2 shift muon -6/(150pi)", dy_mu, 0.0128, 1.0)
m_s = m_b * math.exp(-SIG_D * (1 - (3 / 7 + dy_s)))
m_c = M_T_MEAS * math.exp(-SIG_U * (1 - (4 / 7 - dy_c)))
m_mu = m_e * math.exp(SIG_L * (2 / 3 - dy_mu))
check("m_s (shifted node)", m_s, 93.4, 1.5)
check("m_c (shifted node, from measured m_t)", m_c, 1270.0, 1.0)
check("m_mu (shifted node)", m_mu, 105.66, 0.5)

# --------------------------------------------------------------- mixing = distances
V_us = math.sqrt(m_d / m_s)
V_cb = math.sqrt(m_u / m_c)
check("V_us = sqrt(m_d/m_s)", V_us, 0.2243, 1.5)
check("V_cb = sqrt(m_u/m_c)", V_cb, 0.0422, 6.5)
check("V_us bare = e^{-Nc/2}", math.exp(-1.5), 0.2243, 0.6)
delta_cp = math.degrees(math.atan(4.67 / 2.16))   # measured masses (doublet angle)
check("delta_CP = arctan(m_d/m_u) deg", delta_cp, 65.5, 1.0)

# PMNS: charged-lepton geometry (tribimaximal + muon node shift)
dx = 2 / 3 - math.log(105.66 / 0.511) / math.log(1776.86 / 0.511)
check("muon node shift dx", dx, 0.0128, 1.0)
check("sin^2(th13) = sqrt(3)*dx", math.sqrt(3) * dx, 0.0224, 1.5)
check("sin^2(th12) = 1/3 - 2*dx", 1 / 3 - 2 * dx, 0.307, 0.5)
check("sin^2(th23) = 1/2 + (2+sqrt(3))*dx", 0.5 + (2 + math.sqrt(3)) * dx, 0.546, 0.5)
# cross-sector (quarks -> PMNS, one-way bridge)
check("sin(th12) = sqrt(m_c/m_b)", math.sqrt(m_c / m_b), 0.5541, 3.5)
check("sin(th13) = sqrt(m_s/m_b)", math.sqrt(m_s / m_b), 0.1483, 3.5)

# --------------------------------------------------------------- couplings
check("alpha_s = 3/(8pi)", G2 / (4 * PI), 0.1179, 1.5)
check("sin^2(thW) tree = 2/9 vs on-shell", 2 / 9, 0.2232, 4.0)
check("M_W/M_Z = sqrt(7/9)", math.sqrt(7 / 9), 0.88136, 0.1)
check("y_t = 1 - C0^3", 1 - C0**3, 0.9923, 0.1)
check("Higgs quartic lambda = (1+C0^2)/8", (1 + C0**2) / 8, 0.1291, 1.0)
SIG_FUND = G2 * CF / 2
assert SIG_FUND == 1.0, "sigma_fund must be exactly 1"

# --------------------------------------------------------------- hadrons (junctions)
m_p = N_Q * LAMBDA_BW * (1 - C0**2 * (1 - C0))
check("m_proton = 7*Lambda_BW*screening", m_p, 938.272, 0.5)
np_split = (4.67 - 2.16) / 2 * (1 + C0**2)
check("m_n - m_p", np_split, 1.293, 1.0)
r_p = C0 * N_C * 197.327 / LAMBDA_BW         # hbar*c = 197.327 MeV fm
check("r_proton = C0*Nc*hc/Lambda_BW (fm)", r_p, 0.842, 2.0)

# --------------------------------------------------------------- gravity chain
# competing prefactors (see 04-GAPS: unresolved which is derived)
M_P_MEAS = 1.22089e19 * 1000                  # MeV
exp_warp = math.exp(SIG_D**2)
mp_a = (3 / 4) * (1 - C0**3) * LAMBDA_BW * exp_warp        # G2/CF * y_t route
mp_b = (11 / 15) * LAMBDA_BW * exp_warp                     # 11/15 route
check("M_P route A (3/4)*y_t", mp_a, M_P_MEAS, 2.0)
check("M_P route B 11/15", mp_b, M_P_MEAS, 0.2)
# G from route B
hbar_c = 197.327e-15 * 1e6                    # MeV*m ... only ratio matters:
G_ratio = (M_P_MEAS / mp_b) ** 2              # G_pred/G_meas = (Mp_meas/Mp_pred)^2
check("G accuracy via route B (ratio)", G_ratio, 1.0, 0.3)
check("hierarchy e^{2 sigma_d^2} ~ 1.5e40", math.exp(2 * SIG_D**2) / 1e40, 1.46, 2.0)

# --------------------------------------------------------------- cosmology
check("DM/baryon = (1+C0^2)/C0 = 26/5", (1 + C0**2) / C0, 5.32, 2.5)
check("Hubble ratio H0_local/H0_CMB = 27/25", 27 / 25, 73.04 / 67.4, 0.4)

# --------------------------------------------------------------- report
fails = [c for c in checks if not c[5]]
width = max(len(c[0]) for c in checks)
for name, pred, meas, err, tol, ok in checks:
    flag = "ok  " if ok else "FAIL"
    print(f"{flag} {name:<{width}}  pred={pred:<12.6g} meas={meas:<12.6g} err={err:+.3f}% (tol {tol}%)")
print(f"\n{len(checks) - len(fails)}/{len(checks)} claims verified")
sys.exit(1 if fails else 0)
