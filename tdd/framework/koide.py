"""
Koide formula framework functions.

The Koide formula relates the three charged lepton masses:
  Q = (m_e + m_μ + m_τ) / (√m_e + √m_μ + √m_τ)² = 2/3

The parametrization: √m_k = M(1 + √2 cos(θ + 2πk/3))
reduces three masses to two parameters (M, θ).

Key finding: the Koide angle decomposes as θ = 2π/3 + ε,
where the excess ε = 2/9 radians for leptons (33 ppm precision).
This forms a pattern with Q = 2/3: both involve 2/3^n.

Status legend:
  IMPLEMENTED — formula exists, produces a number
  NOT_IMPLEMENTED — needs derivation, raises NotImplementedError
"""

import math


# ============================================================
# IMPLEMENTED: Koide Q parameter
# ============================================================

def koide_Q(masses: list[float]) -> float:
    """
    Koide quality factor Q = Σm / (Σ√m)².
    IMPLEMENTED.

    Q = 2/3 means the mass point sits at the midpoint of the
    barycentric simplex in √m space. This is the unique value
    where the three √m projections onto the real axis are
    equally spaced around the circle's center.

    Parameters:
        masses: list of three masses [m1, m2, m3] in any units

    Returns:
        Q value (dimensionless). Q = 2/3 for perfect Koide.
    """
    if len(masses) != 3:
        raise ValueError(f"Need exactly 3 masses, got {len(masses)}")
    s = sum(masses)
    sr = sum(math.sqrt(m) for m in masses)
    return s / (sr * sr)


def koide_params(masses: list[float]) -> tuple[float, float]:
    """
    Extract Koide parameters (M, θ) from three masses.
    IMPLEMENTED.

    General parametrization: √m_k = M + R cos(θ + 2πk/3)
    where M = (Σ√m_k)/3 and R = √(2/3 Σ(√m_k - M)²).

    When Q = 2/3 exactly, R = M√2 (Koide's special case).

    Parameters:
        masses: list of three masses [m1, m2, m3]

    Returns:
        (M, theta) where M is in √(mass units) and theta in radians.
    """
    sqrts = [math.sqrt(m) for m in masses]
    M = sum(sqrts) / 3.0

    # General radius R (works for any Q, not just Q = 2/3)
    R = math.sqrt(2.0 / 3.0 * sum((s - M) ** 2 for s in sqrts))

    if R < 1e-15:
        return (M, 0.0)  # degenerate case (all masses equal)

    # Extract theta: (√m_0 - M) / R = cos(θ)
    cos_theta = (sqrts[0] - M) / R
    # From √m_1: (√m_1 - M)/R = cos(θ + 2π/3) = -cos(θ)/2 - sin(θ)√3/2
    cos_theta_shift = (sqrts[1] - M) / R
    sin_theta = -(cos_theta_shift + cos_theta / 2.0) / (math.sqrt(3) / 2.0)
    theta = math.atan2(sin_theta, cos_theta)

    return (M, theta)


def koide_reconstruct(M: float, theta: float, R: float = None) -> list[float]:
    """
    Reconstruct three masses from Koide parameters.
    IMPLEMENTED.

    Uses √m_k = M + R cos(θ + 2πk/3).
    If R is not provided, uses R = M√2 (the Q = 2/3 Koide constraint).

    Parameters:
        M: mass scale parameter (√mass units)
        theta: phase angle (radians)
        R: circle radius. If None, uses M√2 (Koide constraint).

    Returns:
        List of three masses [m1, m2, m3].
    """
    if R is None:
        R = M * math.sqrt(2)
    masses = []
    for k in range(3):
        sqrt_mk = M + R * math.cos(theta + 2 * math.pi * k / 3)
        masses.append(sqrt_mk ** 2)
    return masses


def koide_R(masses: list[float]) -> float:
    """
    Compute the Koide circle radius R from three masses.
    IMPLEMENTED.

    R = √(2/3 Σ(√m_k - M)²). For Q = 2/3, R = M√2.

    Parameters:
        masses: list of three masses

    Returns:
        Circle radius R in √(mass units).
    """
    sqrts = [math.sqrt(m) for m in masses]
    M = sum(sqrts) / 3.0
    return math.sqrt(2.0 / 3.0 * sum((s - M) ** 2 for s in sqrts))


# ============================================================
# IMPLEMENTED: Koide excess angle
# ============================================================

def koide_excess(masses: list[float]) -> float:
    """
    Koide angle excess over 2π/3 (120°).
    IMPLEMENTED.

    The Koide angle θ decomposes as θ = 2π/3 + ε, where the
    120° part is forced by Z₃ symmetry (equal angular spacing)
    and ε encodes all the family-specific information.

    Parameters:
        masses: list of three masses

    Returns:
        Excess ε in radians.
    """
    _, theta = koide_params(masses)
    return theta - 2 * math.pi / 3


def lepton_excess_is_two_ninths(masses: list[float]) -> dict:
    """
    Test whether the lepton Koide excess equals 2/9 radians.
    IMPLEMENTED.

    Discovery: for charged leptons, ε = 2/9 radians to 33 ppm.
    Since Q = 2/3¹ and ε = 2/3² (as 2/9 = 2/3²), both the
    quality factor and the angle follow a 2/3^n pattern.

    Parameters:
        masses: [m_e, m_μ, m_τ] in MeV

    Returns:
        dict with measured excess, target 2/9, absolute and relative errors.
    """
    excess = koide_excess(masses)
    target = 2.0 / 9.0
    abs_err = abs(excess - target)
    rel_err = abs_err / target

    return {
        'excess_rad': excess,
        'target_rad': target,
        'abs_error': abs_err,
        'rel_error_ppm': rel_err * 1e6,
    }


# ============================================================
# IMPLEMENTED: D=27 quantization
# ============================================================

def quantized_excess_27(masses: list[float]) -> dict:
    """
    Check if the Koide excess is quantized in units of 1/27.
    IMPLEMENTED.

    Common denominator D=27=3³ works for all three families:
      Leptons:     ε ≈ 6/27 = 2/9  (0.003% off)
      Up quarks:   ε ≈ 2/27        (0.85% off)
      Down quarks: ε ≈ 3/27 = 1/9  (0.48% off)

    Parameters:
        masses: list of three masses

    Returns:
        dict with measured excess, nearest n/27, and error.
    """
    excess = koide_excess(masses)
    # Find nearest integer numerator n such that excess ≈ n/27
    n_float = excess * 27.0
    n_nearest = round(n_float)
    target = n_nearest / 27.0
    abs_err = abs(excess - target)
    rel_err = abs_err / target if target != 0 else float('inf')

    return {
        'excess_rad': excess,
        'n_over_27': n_nearest,
        'target_rad': target,
        'rel_error_pct': rel_err * 100,
    }


# ============================================================
# NOT IMPLEMENTED: Derivations needed
# ============================================================

def derive_koide_excess() -> float:
    """
    Derive the lepton Koide excess ε = 2/9 from first principles.
    NOT_IMPLEMENTED.

    The pattern suggests a perturbative expansion in 1/3:
      Q = 2/3¹ (zeroth order — triangle constraint)
      ε = 2/3² (first correction — ???)

    If this is a genuine expansion, the next term would be
    2/3³ = 2/27 — which IS the up quark excess.

    Possible origin:
      - Yukawa coupling structure with Z₃ selection rule
      - Radiative correction proportional to α × (color factor)
      - Boundary lattice perturbation theory
    """
    raise NotImplementedError(
        "TODO: derive ε = 2/9 from boundary geometry or gauge structure. "
        "Pattern: {Q, ε_lep, ε_up} = {2/3, 2/9, 2/27} = 2/3^{1,2,3}. "
        "Suggestive of perturbation series in 1/3."
    )


def derive_M_parameter() -> float:
    """
    Derive the Koide mass scale parameter M from first principles.
    NOT_IMPLEMENTED.

    After Q = 2/3 and ε = 2/9 fix two of three degrees of freedom,
    M is the single remaining parameter that sets the overall mass scale.

    M_leptons = 17.716 MeV^(1/2)

    This must ultimately come from the boundary's fundamental scale
    (Planck mass) via the vacuum dressing mechanism.
    """
    raise NotImplementedError(
        "TODO: derive M = 17.716 MeV^(1/2) for leptons from boundary geometry. "
        "After Q=2/3 and ε=2/9, M is the only free parameter."
    )


def predict_neutrino_masses(excess_quantized: float = None) -> list[float]:
    """
    Predict neutrino masses from quantized Koide excess.
    NOT_IMPLEMENTED.

    If neutrinos also satisfy Koide with a quantized excess n/27,
    combined with Δm² measurements, we can predict the lightest mass.
    The value of n for neutrinos is unknown.

    Returns:
        List of three neutrino masses in eV.
    """
    raise NotImplementedError(
        "TODO: determine neutrino Koide excess quantum number n. "
        "Then use Δm²₂₁ and Δm²₃₂ to predict absolute masses. "
        "JUNO and KATRIN experiments will test this."
    )
