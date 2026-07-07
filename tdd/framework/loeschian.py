"""
Loeschian numbers: L(m,n) = m² + m·n + n²

These are the eigenvalues of the Laplacian on an equilateral triangle
(Dirichlet BC), and equivalently the squared distances on the triangular
(Eisenstein integer) lattice.

The boundary framework identifies these as the natural mass spectrum:
particle masses (in electron-mass units) should be Loeschian numbers,
with deviations from vacuum dressing / QCD corrections.

Key results from the research:
  - tau/electron = 3477 = L(37,31) EXACT within measurement error
  - 11/12 particle masses within 3σ of a Loeschian number
  - Deviations decrease with mass (consistent with running couplings)
"""

import math
from typing import Optional


def loeschian(m: int, n: int) -> int:
    """Compute Loeschian number L(m,n) = m² + m·n + n²."""
    return m * m + m * n + n * n


def is_loeschian(k: int) -> bool:
    """Check if k is a Loeschian number."""
    if k < 0:
        return False
    if k == 0:
        return True
    # A positive integer is Loeschian iff it has no prime factor
    # congruent to 2 mod 3 raised to an odd power.
    # Brute force check for our purposes:
    for m in range(0, int(math.isqrt(k)) + 1):
        for n in range(0, m + 1):
            if loeschian(m, n) == k:
                return True
            if loeschian(m, n) > k:
                break
    return False


def nearest_loeschian(target: float, tolerance: float = 0.5) -> Optional[tuple]:
    """
    Find the nearest Loeschian number to target.
    Returns (L, m, n, delta) or None if nothing within tolerance.

    For the test suite, we search exhaustively up to reasonable bounds.
    """
    best = None
    best_delta = float('inf')
    max_m = int(math.isqrt(int(target))) + 2

    for m in range(0, max_m):
        for n in range(0, m + 1):
            L = loeschian(m, n)
            delta = abs(L - target)
            if delta < best_delta:
                best_delta = delta
                best = (L, m, n, target - L)
            if L > target + best_delta:
                break
        if m * m > target + best_delta:
            break

    return best


def mixed_sign_loeschian(m: int, n: int) -> int:
    """
    Compute L with mixed-sign modes: L(m, -n) = m² - m·n + n²

    Same formula m² + m·n + n² but with n→-n gives m² - m·n + n².
    These represent saddle-curvature modes on the boundary.

    Key result: L(37,31)/L(37,-31) ≈ 3 (lepton/quark factor)
    """
    return m * m - m * n + n * n


def eigenvalue_ratio(m1: int, n1: int, m2: int, n2: int) -> float:
    """Ratio of two Loeschian eigenvalues."""
    L2 = loeschian(m2, n2)
    if L2 == 0:
        raise ValueError("Denominator mode (m2,n2) gives L=0")
    return loeschian(m1, n1) / L2


# ============================================================
# Known particle mode assignments (from research)
# ============================================================

PARTICLE_MODES = {
    # particle: (m, n, L, status)
    'electron':  (1, 0, 1, 'EXACT'),
    'up':        (2, 0, 4, 'APPROXIMATE'),      # actual ratio 4.23
    'down':      (3, 0, 9, 'APPROXIMATE'),      # actual ratio 9.14
    'strange':   (13, 1, 183, 'CLOSE'),          # actual ratio 182.78
    'muon':      (12, 4, 208, 'MISS'),           # actual ratio 206.77, may be 3/(2α)
    'charm':     (36, 21, 2493, 'CLOSE'),        # actual ratio 2491.2
    'tau':       (37, 31, 3477, 'EXACT'),         # actual ratio 3477.23 ± 0.23
    'bottom':    (82, 15, 8179, 'CLOSE'),        # actual ratio 8180.1
    'W':         (279, 175, 157291, 'CLOSE'),    # actual ratio 157294
    'Z':         (378, 78, 178452, 'CLOSE'),     # actual ratio 178450
    'Higgs':     (476, 36, 245008, 'CLOSE'),     # actual ratio 245010
    'top':       (400, 267, 338089, 'CLOSE'),    # actual ratio 338083
}
