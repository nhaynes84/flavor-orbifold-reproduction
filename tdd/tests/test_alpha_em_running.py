"""
α_em test: g² = 3/2 is the SU(3) boundary coupling, NOT a universal coupling.

The naive calculation e² = g²sin²θ_W = (3/2)(2/9) = 1/3 gives
α_em = 1/(12π) = 0.027, which is 3.6× above 1/137. This is NOT a
framework failure — it's an incorrect identification: g₃ ≠ g₂.

In the SM, the three gauge couplings are distinct at M_Z:
  g₃²/(4π) ≈ 0.118 (strong)
  g₂²/(4π) ≈ 0.034 (weak SU(2))
  g₁²/(4π) ≈ 0.010 (hypercharge)

The framework derives g₃² = 3/2 from boundary string tension.
The electroweak couplings g₂ and g₁ are NOT derived — they're
separate open problems. Setting g₂ = g₃ is the error, not QED running.

QED running verification: even IF 1/(12π) were the right boundary
value, standard running cannot bridge to 1/137 at any physical scale.
This confirms that g₂ ≠ g₃ on the boundary.
"""

import math
import pytest

ALPHA_0 = 1 / 137.036
ALPHA_BOUNDARY = 1 / (12 * math.pi)
M_Z = 91.19  # GeV
M_PLANCK = 1.22e19  # GeV


class TestAlphaEmRunning:
    """Does QED running resolve the α_em discrepancy?"""

    def test_boundary_alpha_is_3p6x_measured(self):
        """The boundary value 1/(12π) is 3.6× larger than 1/137."""
        ratio = ALPHA_BOUNDARY / ALPHA_0
        assert 3.5 < ratio < 3.8

    def test_delta_alpha_needed(self):
        """Need Δα = 0.72 to bridge the gap. Known Δα(M_Z) = 0.059."""
        delta_needed = 1 - ALPHA_0 / ALPHA_BOUNDARY
        delta_at_MZ = 0.059  # measured
        assert delta_needed > 10 * delta_at_MZ, \
            "Gap requires >10× the known running — QED can't do it"

    def test_one_loop_at_planck(self):
        """One-loop QED running to Planck scale only reaches ~1/62."""
        # Sum N_c * Q_f^2 over all SM fermions = 8
        sum_NcQ2 = 8
        coeff = 2 * ALPHA_0 / (3 * math.pi) * sum_NcQ2
        ln_planck = math.log(M_PLANCK / 1.0)  # from ~1 GeV
        delta_planck = coeff * ln_planck
        alpha_planck = ALPHA_0 / (1 - delta_planck)

        # At Planck: get ~1/62, need 1/37.7
        assert 1 / alpha_planck > 55, f"alpha(Planck) too large: 1/{1/alpha_planck:.0f}"
        assert 1 / alpha_planck < 70, f"alpha(Planck) too small: 1/{1/alpha_planck:.0f}"

        # Still ~1.6-2× short of boundary value
        gap_remaining = ALPHA_BOUNDARY / alpha_planck
        assert gap_remaining > 1.5, \
            f"Gap at Planck only {gap_remaining:.1f}× — should be ~1.9×"

    def test_scale_needed_above_planck(self):
        """The scale where α = 1/(12π) is ~10^6 above Planck."""
        sum_NcQ2 = 8
        coeff = 2 * ALPHA_0 / (3 * math.pi) * sum_NcQ2
        delta_needed = 1 - ALPHA_0 / ALPHA_BOUNDARY
        ln_needed = delta_needed / coeff
        mu_needed = math.exp(ln_needed)  # in GeV (from 1 GeV)

        ratio_to_planck = mu_needed / M_PLANCK
        assert ratio_to_planck > 1e5, \
            f"Scale only {ratio_to_planck:.0e}× above Planck"

    def test_qcd_matches_but_qed_doesnt(self):
        """α_s = g²/(4π) = 3/(8π) matches at 1.2%, but α_em doesn't.
        This asymmetry is the core of the problem."""
        alpha_s_pred = 3 / (8 * math.pi)
        alpha_s_meas = 0.1179
        alpha_s_err = abs(alpha_s_pred - alpha_s_meas) / alpha_s_meas

        alpha_em_pred = ALPHA_BOUNDARY  # 1/(12π)
        alpha_em_meas = ALPHA_0
        alpha_em_err = abs(alpha_em_pred - alpha_em_meas) / alpha_em_meas

        # QCD: 1.2% match. QED: 263% miss.
        assert alpha_s_err < 0.02, "α_s should match within 2%"
        assert alpha_em_err > 2.0, "α_em should miss by >200%"
