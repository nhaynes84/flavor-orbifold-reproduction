"""
Fundamental constants and Planck units.

All values from NIST CODATA 2018 / PDG 2024 unless noted.
The boundary framework treats ħ, c, G as the three grid parameters:
  ħ = action per pixel-tick (resolution)
  c = Lp/Tp (rendering speed)
  G = saturation coupling (self-interaction)
"""

import math

# ============================================================
# FUNDAMENTAL CONSTANTS (NIST CODATA 2018)
# ============================================================

# Reduced Planck constant
hbar = 1.054571817e-34       # J·s (exact in 2019 SI)
hbar_eV = 6.582119569e-16    # eV·s

# Speed of light
c = 299792458.0              # m/s (exact)

# Gravitational constant
G = 6.67430e-11              # m³/(kg·s²), uncertainty ±0.00015e-11

# Fine structure constant
alpha = 7.2973525693e-3      # dimensionless, uncertainty ±0.0000000011e-3
alpha_inv = 1.0 / alpha      # ≈ 137.035999084

# Boltzmann constant
k_B = 1.380649e-23           # J/K (exact in 2019 SI)

# Elementary charge
e_charge = 1.602176634e-19   # C (exact in 2019 SI)

# Electron volt
eV = 1.602176634e-19         # J (exact)
MeV = eV * 1e6
GeV = eV * 1e9

# ============================================================
# PLANCK UNITS — Grid specifications
# ============================================================

# Planck length: pixel size
L_P = math.sqrt(hbar * G / c**3)  # ≈ 1.616255e-35 m

# Planck time: tick duration
T_P = math.sqrt(hbar * G / c**5)  # ≈ 5.391247e-44 s

# Planck mass: max energy per pixel
M_P = math.sqrt(hbar * c / G)     # ≈ 2.176434e-8 kg

# Planck energy
E_P = M_P * c**2                   # ≈ 1.956e9 J ≈ 1.22e19 GeV

# Planck temperature
T_P_temp = E_P / k_B              # ≈ 1.416784e32 K

# Planck area (one boundary pixel)
A_P = L_P**2                       # ≈ 2.612e-70 m²

# ============================================================
# DERIVED GRID RELATIONSHIPS
# ============================================================

# c = L_P / T_P (rendering speed = 1 pixel per tick)
# c² = L_P² / T_P² (boundary area element in pixel-tick space)
# L_P * M_P * c = hbar (one pixel × max momentum = one action quantum)
