"""
Shared physical constants and measurements for the test suite.

All values from NIST CODATA 2018 and PDG 2024 (Particle Data Group)
unless otherwise noted. Uncertainties are 1σ.

This file is the single source of truth for measured values.
Tests compare framework predictions against these numbers.
"""

import pytest
import math


# ============================================================
# PYTEST MARKERS
# ============================================================

def pytest_configure(config):
    config.addinivalue_line("markers", "data: Data validation — verifies reference values")
    config.addinivalue_line("markers", "framework: Framework prediction tests")
    config.addinivalue_line("markers", "two_category: Two-category framework tests (coherence vs eigenvalue)")
    config.addinivalue_line("markers", "decoherence: Decoherence spectrum tests")


# ============================================================
# FUNDAMENTAL CONSTANTS (inputs to the framework)
# ============================================================

HBAR = 1.054571817e-34          # J·s (exact post-2019 SI)
HBAR_UNC = 0.0
C = 299792458.0                 # m/s (exact)
C_UNC = 0.0
G = 6.67430e-11                 # m³/(kg·s²)
G_UNC = 0.00015e-11
ALPHA = 7.2973525693e-3         # dimensionless
ALPHA_UNC = 0.0000000011e-3
ALPHA_INV = 1.0 / ALPHA         # 137.035999084
K_B = 1.380649e-23              # J/K (exact)
E_CHARGE = 1.602176634e-19      # C (exact)
EV = 1.602176634e-19            # J
MEV = EV * 1e6
GEV = EV * 1e9


# ============================================================
# PLANCK UNITS (derived — grid specifications)
# ============================================================

L_PLANCK = 1.616255e-35         # m
L_PLANCK_UNC = 0.000018e-35
T_PLANCK = 5.391247e-44         # s
T_PLANCK_UNC = 0.000060e-44
M_PLANCK = 2.176434e-8          # kg
M_PLANCK_UNC = 0.000024e-8
E_PLANCK = 1.220890e19          # GeV
T_PLANCK_TEMP = 1.416784e32     # K


# ============================================================
# PARTICLE MASSES (PDG 2024) — MeV/c²
# ============================================================

M_ELECTRON = 0.51099895000
M_ELECTRON_UNC = 0.00000000015
M_MUON = 105.6583755
M_MUON_UNC = 0.0000023
M_TAU = 1776.86
M_TAU_UNC = 0.12

M_NU_E_UPPER = 0.0008
M_NU_MU_UPPER = 0.19
M_NU_TAU_UPPER = 18.2
DELTA_M21_SQ = 7.53e-5          # eV²
DELTA_M32_SQ = 2.453e-3         # eV²

M_UP = 2.16
M_UP_UNC_PLUS = 0.49
M_UP_UNC_MINUS = 0.26
M_DOWN = 4.67
M_DOWN_UNC_PLUS = 0.48
M_DOWN_UNC_MINUS = 0.17
M_STRANGE = 93.4
M_STRANGE_UNC_PLUS = 8.6
M_STRANGE_UNC_MINUS = 3.4
M_CHARM = 1270.0
M_CHARM_UNC = 20.0
M_BOTTOM = 4180.0
M_BOTTOM_UNC_PLUS = 30.0
M_BOTTOM_UNC_MINUS = 20.0
M_TOP = 172760.0
M_TOP_UNC = 300.0

M_W = 80369.0
M_W_UNC = 13.0
M_Z = 91187.6
M_Z_UNC = 2.1
M_HIGGS = 125250.0
M_HIGGS_UNC = 170.0
M_PHOTON = 0.0
M_GLUON = 0.0

M_PROTON = 938.27208816
M_PROTON_UNC = 0.00000029
M_NEUTRON = 939.56542052
M_NEUTRON_UNC = 0.00000054


# ============================================================
# QCD SCALE
# ============================================================

LAMBDA_QCD = 213.0                  # MeV (PDG 2024, N_f=5, MS-bar)
LAMBDA_QCD_UNC = 8.0


# ============================================================
# CKM MIXING PARAMETERS (PDG 2024, Wolfenstein)
# ============================================================

WOLFENSTEIN_LAMBDA = 0.22453        # λ (= |V_us|)
WOLFENSTEIN_LAMBDA_UNC = 0.00044
WOLFENSTEIN_A = 0.836               # A
WOLFENSTEIN_A_UNC = 0.015

# CKM matrix elements (magnitudes)
V_UD = 0.97373
V_UD_UNC = 0.00031
V_US = 0.22453                      # = λ
V_US_UNC = 0.00044
V_UB = 0.00382
V_UB_UNC = 0.00020
V_CD = 0.22438
V_CD_UNC = 0.00044
V_CS = 0.97350
V_CS_UNC = 0.00031
V_CB = 0.04080
V_CB_UNC = 0.00080
V_TD = 0.00860
V_TD_UNC = 0.00020
V_TS = 0.04010
V_TS_UNC = 0.00090
V_TB = 0.999118
V_TB_UNC = 0.000033

# CKM mixing angles (radians, standard parametrization)
THETA_12_CKM = math.asin(V_US)     # ~12.96°
THETA_23_CKM = math.asin(V_CB)     # ~2.34°
THETA_13_CKM = math.asin(V_UB)     # ~0.219°

# CKM CP-violating phase
DELTA_CP_CKM = 1.144               # radians (~65.5°)
DELTA_CP_CKM_UNC = 0.027


# ============================================================
# PMNS MIXING PARAMETERS (PDG 2024, Normal Hierarchy)
# ============================================================

# PMNS mixing angles (radians)
THETA_12_PMNS = math.asin(math.sqrt(0.303))   # sin²θ₁₂ = 0.303, θ₁₂ ≈ 33.4°
THETA_23_PMNS = math.asin(math.sqrt(0.572))   # sin²θ₂₃ = 0.572, θ₂₃ ≈ 49.1°
THETA_13_PMNS = math.asin(math.sqrt(0.02203)) # sin²θ₁₃ = 0.02203, θ₁₃ ≈ 8.53°

SIN2_THETA_12_PMNS = 0.303
SIN2_THETA_12_PMNS_UNC = 0.012
SIN2_THETA_23_PMNS = 0.572
SIN2_THETA_23_PMNS_UNC = 0.024
SIN2_THETA_13_PMNS = 0.02203
SIN2_THETA_13_PMNS_UNC = 0.00056

# PMNS sin(theta) values (for direct comparison with predictions)
SIN_THETA_12_PMNS = math.sqrt(SIN2_THETA_12_PMNS)  # 0.5505
SIN_THETA_23_PMNS = math.sqrt(SIN2_THETA_23_PMNS)  # 0.7563
SIN_THETA_13_PMNS = math.sqrt(SIN2_THETA_13_PMNS)  # 0.1484

# PMNS CP-violating phase (NH)
DELTA_CP_PMNS = 3.86                # radians (~221°)
DELTA_CP_PMNS_UNC = 0.52

# Wolfenstein rho, eta (PDG 2024)
WOLFENSTEIN_RHO_BAR = 0.159
WOLFENSTEIN_RHO_BAR_UNC = 0.010
WOLFENSTEIN_ETA_BAR = 0.349
WOLFENSTEIN_ETA_BAR_UNC = 0.010
WOLFENSTEIN_RHO_ETA_MAG = math.sqrt(WOLFENSTEIN_RHO_BAR**2 + WOLFENSTEIN_ETA_BAR**2)  # |ρ-iη| ≈ 0.384
WOLFENSTEIN_RHO_ETA_ARG = math.atan2(WOLFENSTEIN_ETA_BAR, WOLFENSTEIN_RHO_BAR)         # arg(ρ+iη) ≈ 1.144 rad

# Jarlskog invariant (PDG 2024)
JARLSKOG_J = 3.08e-5
JARLSKOG_J_UNC = 0.14e-5

# Unitarity triangle angles (PDG 2024, degrees)
UT_ALPHA = 84.5                 # α = φ₂
UT_ALPHA_UNC = 4.5
UT_BETA = 22.2                  # β = φ₁
UT_BETA_UNC = 0.7
UT_GAMMA = 65.4                 # γ = φ₃
UT_GAMMA_UNC = 3.2


# ============================================================
# NEUTRAL MESON MIXING PARAMETERS (PDG 2024)
# ============================================================

# Mixing parameter x = Δm/Γ
X_K0 = 0.9456
X_D0 = 0.00398
X_D0_UNC = 0.00049
X_B0 = 0.769
X_B0_UNC = 0.004
X_BS = 26.89
X_BS_UNC = 0.07

# Meson masses (MeV)
M_K0 = 497.611
M_D0 = 1864.84
M_B0 = 5279.66
M_BS = 5366.92

# Meson lifetimes (seconds)
TAU_K0_L = 5.116e-8
TAU_D0 = 4.103e-13
TAU_D_PLUS = 1.033e-12
TAU_DS = 5.04e-13
TAU_B0 = 1.519e-12
TAU_B_PLUS = 1.638e-12
TAU_BS = 1.515e-12
TAU_BC = 0.510e-12

# D meson and B meson lifetime ratios (for spread calculations)
D_LIFETIME_SPREAD = TAU_D_PLUS / TAU_D0       # ~2.52
B_LIFETIME_SPREAD = TAU_B_PLUS / TAU_B0       # ~1.08


# ============================================================
# KEY MASS RATIOS
# ============================================================

RATIO_PROTON_ELECTRON = M_PROTON / M_ELECTRON
RATIO_MUON_ELECTRON = M_MUON / M_ELECTRON
RATIO_TAU_ELECTRON = M_TAU / M_ELECTRON
RATIO_TAU_ELECTRON_UNC = M_TAU_UNC / M_ELECTRON


# ============================================================
# COUPLING CONSTANTS
# ============================================================

ALPHA_S_MZ = 0.1179
ALPHA_S_MZ_UNC = 0.0009
SIN2_THETA_W = 0.23122           # MS-bar scheme (PDG 2024)
SIN2_THETA_W_UNC = 0.00004
G_FERMI = 1.1663788e-5          # GeV⁻²

# On-shell Weinberg angle (from pole masses: 1 - M_W²/M_Z²)
SIN2_THETA_W_ONSHELL = 1.0 - (M_W / M_Z)**2  # ~0.22320
# Propagated uncertainty: 2(M_W/M_Z²)√(σ_W² + (M_W σ_Z/M_Z)²)
SIN2_THETA_W_ONSHELL_UNC = 2.0 * (M_W / M_Z**2) * math.sqrt(
    M_W_UNC**2 + (M_W * M_Z_UNC / M_Z)**2
)

# Individual W mass measurements (MeV) — for axiom testing
M_W_CDF = 80433.5               # CDF II, Science 376 (2022)
M_W_CDF_UNC = 9.4
M_W_ATLAS = 80366.5             # ATLAS, Eur.Phys.J.C 84 (2024)
M_W_ATLAS_UNC = 15.9
M_W_CMS = 80360.2               # CMS, preliminary (2024)
M_W_CMS_UNC = 9.9


# ============================================================
# COSMOLOGICAL OBSERVABLES
# ============================================================

H0_CMB = 67.4                    # km/s/Mpc
H0_CMB_UNC = 0.5
H0_LOCAL = 73.04                 # km/s/Mpc
H0_LOCAL_UNC = 1.04
MPC_IN_M = 3.0857e22
H0_CMB_SI = H0_CMB * 1e3 / MPC_IN_M
H0_LOCAL_SI = H0_LOCAL * 1e3 / MPC_IN_M

T_CMB = 2.7255                   # K
T_CMB_UNC = 0.0006
A_S = 2.1e-9
A_S_UNC = 0.03e-9
SQRT_A_S = math.sqrt(A_S)

OMEGA_M = 0.315
OMEGA_M_UNC = 0.007
OMEGA_LAMBDA = 0.685
OMEGA_LAMBDA_UNC = 0.007
OMEGA_B = 0.0493
OMEGA_B_UNC = 0.0006
OMEGA_CDM = OMEGA_M - OMEGA_B

# Observable universe parameters
R_OBS_UNIVERSE = 4.4e26          # m (comoving radius)
M_OBS_UNIVERSE = 1.5e53          # kg (total mass-energy equivalent)

A0_MOND = 1.2e-10               # m/s²
A0_MOND_UNC = 0.2e-10

KBC_DELTA = 0.46
KBC_DELTA_UNC = 0.06

F_BOOST_HUBBLE = 1.72
F_BOOST_SNE = 1.742
F_AREA = 0.219
SQRT_3 = math.sqrt(3)

# DESI BAO measurements (DR2, 2025)
W0_DESI = -0.752
W0_DESI_UNC = 0.058
WA_DESI = -0.86
WA_DESI_UNC_PLUS = 0.28
WA_DESI_UNC_MINUS = 0.28

# SPARC rotation curve analysis (2026-04-05)
# Dataset: Lelli, McGaugh & Schombert 2016, 175 galaxies
SPARC_N_GALAXIES = 171          # galaxies with ≥5 data points
SPARC_N_POINTS = 3375           # total data points after quality cuts
SPARC_CHI2_BARYON = 579.05      # χ²/N for baryonic-only (no DM, no MOND)
SPARC_CHI2_FRAMEWORK = 52.46    # χ²/N for MOND with framework a₀ = cH₀_CMB/(2π)
SPARC_CHI2_EMPIRICAL = 50.20    # χ²/N for MOND with empirical a₀ = 1.2e-10
SPARC_WINS_FRAMEWORK = 81       # galaxies where framework a₀ gives best fit
SPARC_WINS_EMPIRICAL = 70       # galaxies where empirical a₀ gives best fit
SPARC_RAR_SCATTER_FW = 0.1874   # dex, RMS scatter in RAR (framework a₀)
SPARC_RAR_SCATTER_EMP = 0.1830  # dex, RMS scatter in RAR (empirical a₀)
SPARC_RAR_SCATTER_MCGAUGH = 0.13  # dex, McGaugh+2016 reported scatter
SPARC_YD = 0.5                  # disk mass-to-light ratio (fixed)
SPARC_YB = 0.7                  # bulge mass-to-light ratio (fixed)

# ISW / CMB Cold Spot observations
ISW_AISW_OBSERVED = 5.2         # Granett, Neyrinck & Szapudi 2008
ISW_AISW_UNC = 1.6
ISW_AISW_LCDM = 1.0            # ΛCDM prediction
COLD_SPOT_DEFICIT_UK = 150.0    # μK, peak deficit
COLD_SPOT_ISW_UK = 20.0         # μK, ISW accounts for this much (10-20%)
COLD_SPOT_UNEXPLAINED_FRAC = 0.85  # ~80-90% unexplained

# Schwarzschild QNM quality factors (Berti, Cardoso & Starinets 2009)
# Q(ℓ) = ω_R / (2|ω_I|) for fundamental (n=0) modes
QNM_QUALITY_FACTORS = {
    2: 2.100, 3: 3.233, 4: 4.296, 5: 5.335,
    6: 6.361, 7: 7.381, 8: 8.396, 9: 9.410, 10: 10.424,
}
# Schwarzschild QNM frequencies (Mω units)
QNM_FREQ_REAL = {
    2: 0.37367, 3: 0.59944, 4: 0.80918, 5: 1.01229,
    6: 1.21215, 7: 1.41019, 8: 1.60700, 9: 1.80299, 10: 1.99839,
}
QNM_FREQ_IMAG = {
    2: 0.08896, 3: 0.09270, 4: 0.09416, 5: 0.09487,
    6: 0.09527, 7: 0.09552, 8: 0.09568, 9: 0.09579, 10: 0.09587,
}

# Planck 2018 low-ℓ TT power spectrum (Commander, D_ℓ in μK²)
CMB_D_OBS = {
    2: 227, 3: 1017, 4: 558, 5: 779, 6: 401, 7: 1309,
    8: 873, 9: 1095, 10: 785, 11: 1038, 12: 1386, 13: 802,
    14: 1233, 15: 901, 16: 946, 17: 867, 18: 1287, 19: 1115, 20: 1127,
}
CMB_D_LCDM = {
    2: 1150, 3: 1010, 4: 960, 5: 940, 6: 950, 7: 975,
    8: 1000, 9: 1020, 10: 1040, 11: 1055, 12: 1065, 13: 1075,
    14: 1085, 15: 1095, 16: 1100, 17: 1110, 18: 1115, 19: 1120, 20: 1130,
}
CMB_L2_SUPPRESSION = CMB_D_OBS[2] / CMB_D_LCDM[2]  # ~0.197
CMB_L2_LEAKAGE = 1.0 - CMB_L2_SUPPRESSION            # ~0.803

# CMB anomaly observations
CMB_HEMI_ASYMMETRY = 0.07       # hemispheric power asymmetry A
CMB_HEMI_ASYMMETRY_UNC = 0.02
CMB_AXIS_OF_EVIL_L = 250.0      # galactic longitude, degrees
CMB_AXIS_OF_EVIL_B = 60.0       # galactic latitude, degrees

# Parent BH inferred properties (from parent_universe_inversion.mjs)
PARENT_SPIN_ESTIMATE = 0.1      # a* from hemispheric asymmetry


# ============================================================
# COSMIC RAY MUON EXCESS & STRANGENESS ENHANCEMENT
# ============================================================

# Auger 2021 (inclined showers, EPOS-LHC): R_μ = N_μ(data) / N_μ(model)
R_MU_AUGER_2021 = 1.38
R_MU_AUGER_2021_UNC = 0.04
R_MU_AUGER_2021_ENERGY_EV = 2e18

# Auger 2015 (hybrid, EPOS-LHC)
R_MU_AUGER_2015 = 1.33
R_MU_AUGER_2015_UNC = 0.16

# IceCube (GeV muons, SIBYLL 2.3)
R_MU_ICECUBE = 1.08
R_MU_ICECUBE_UNC = 0.05
R_MU_ICECUBE_ENERGY_EV = 1e15

# ALICE strangeness enhancement (Nature Physics 13, 535, 2017)
# Enhancement = yield / low-mult pp baseline
ALICE_K0S_PP_HIGH = 1.27       # K⁰_S in high-mult pp
ALICE_K0S_PBPB = 1.80          # K⁰_S in Pb-Pb central
ALICE_LAMBDA_PP_HIGH = 1.40    # Λ in high-mult pp
ALICE_XI_PP_HIGH = 1.75        # Ξ⁻ in high-mult pp (|S|=2)
ALICE_OMEGA_PP_HIGH = 2.00     # Ω⁻ in high-mult pp (|S|=3)
ALICE_OMEGA_PBPB = 8.0         # Ω⁻ in Pb-Pb central

# Pion/kaon critical energies (GeV) — where decay prob = interaction prob
E_CRIT_PION = 115              # GeV at sea level
E_CRIT_KAON = 850              # GeV at sea level


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def within_sigma(measured, predicted, uncertainty, n_sigma=3):
    """Check if predicted is within n_sigma of measured."""
    return abs(predicted - measured) <= n_sigma * uncertainty


def relative_error(measured, predicted):
    """Fractional error |predicted - measured| / measured."""
    if measured == 0:
        return float('inf') if predicted != 0 else 0.0
    return abs(predicted - measured) / abs(measured)
