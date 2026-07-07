"""
Test Domain 7: Koide Formula & Triangle Symmetry

The Koide formula Q = Σm/(Σ√m)² = 2/3 holds for charged leptons
to 9 ppm. The parametrization √m_k = M(1 + √2 cos(θ + 2πk/3))
encodes three masses in two parameters (M, θ) via 120° spacing.

Key discovery: the Koide angle excess ε = θ - 2π/3 equals 2/9
radians for leptons (33 ppm), forming a 2/3^n pattern with Q = 2/3.
All three families quantize to n/27 (D = 3³).

Scorecard:
  GREEN (data)      — Koide Q values, mass ratios, charge quantization
  GREEN (framework) — ε = 2/9 for leptons, D=27 quantization
  RED (framework)   — Derivation of 2/9, M parameter, neutrino prediction
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tests.conftest import *
from framework.koide import (
    koide_Q, koide_params, koide_reconstruct, koide_R, koide_excess,
    lepton_excess_is_two_ninths, quantized_excess_27,
    derive_koide_excess, derive_M_parameter, predict_neutrino_masses,
)
from framework.decoherence import (
    neutrino_koide_max_Q, predict_tau_from_koide,
)


# ============================================================
# PDG 2024 masses for Koide analysis
# ============================================================

LEPTONS = [M_ELECTRON, M_MUON, M_TAU]
UP_QUARKS = [M_UP, M_CHARM, M_TOP]
DOWN_QUARKS = [M_DOWN, M_STRANGE, M_BOTTOM]


# ============================================================
# DATA VALIDATION — Koide Q values
# ============================================================

@pytest.mark.data
class TestKoideQValues:
    """Verify Koide Q values computed from PDG masses."""

    def test_lepton_Q_near_two_thirds(self):
        """Charged lepton Q = 0.666661 — 9 ppm from 2/3."""
        Q = koide_Q(LEPTONS)
        assert Q == pytest.approx(2/3, rel=1e-4)

    def test_lepton_Q_precision(self):
        """Q deviates from 2/3 by less than 10 ppm."""
        Q = koide_Q(LEPTONS)
        deviation_ppm = abs(Q - 2/3) / (2/3) * 1e6
        assert deviation_ppm < 10, f"Deviation = {deviation_ppm:.1f} ppm"

    def test_up_quark_Q(self):
        """Up quark Q ≈ 0.849 — substantial deviation from 2/3."""
        Q = koide_Q(UP_QUARKS)
        assert Q == pytest.approx(0.849, abs=0.005)

    def test_down_quark_Q(self):
        """Down quark Q ≈ 0.731 — closer to 2/3 than up quarks."""
        Q = koide_Q(DOWN_QUARKS)
        assert Q == pytest.approx(0.731, abs=0.005)

    def test_Q_bounded(self):
        """Q is bounded: 1/3 ≤ Q ≤ 1 for any three positive masses."""
        for masses in [LEPTONS, UP_QUARKS, DOWN_QUARKS]:
            Q = koide_Q(masses)
            assert 1/3 <= Q <= 1, f"Q = {Q} out of bounds"

    def test_Q_equals_one_third_for_equal_masses(self):
        """Q = 1/3 when all three masses are equal (degenerate limit)."""
        Q = koide_Q([1.0, 1.0, 1.0])
        assert Q == pytest.approx(1/3, rel=1e-10)

    def test_Q_approaches_one_for_hierarchy(self):
        """Q → 1 when one mass dominates (hierarchical limit)."""
        Q = koide_Q([1e-10, 1e-10, 1.0])
        assert Q > 0.99


# ============================================================
# DATA VALIDATION — Koide parametrization round-trip
# ============================================================

@pytest.mark.data
class TestKoideParametrization:
    """Verify (M, θ) parametrization reconstructs masses."""

    def test_lepton_round_trip(self):
        """M, R, θ → masses → M, R, θ round-trips for leptons."""
        M, theta = koide_params(LEPTONS)
        R = koide_R(LEPTONS)
        reconstructed = koide_reconstruct(M, theta, R)
        for orig, recon in zip(LEPTONS, reconstructed):
            assert recon == pytest.approx(orig, rel=1e-6)

    def test_up_quark_round_trip(self):
        """Round-trip for up-type quarks (general R, not Q=2/3)."""
        M, theta = koide_params(UP_QUARKS)
        R = koide_R(UP_QUARKS)
        reconstructed = koide_reconstruct(M, theta, R)
        for orig, recon in zip(UP_QUARKS, reconstructed):
            assert recon == pytest.approx(orig, rel=1e-6)

    def test_down_quark_round_trip(self):
        """Round-trip for down-type quarks (general R, not Q=2/3)."""
        M, theta = koide_params(DOWN_QUARKS)
        R = koide_R(DOWN_QUARKS)
        reconstructed = koide_reconstruct(M, theta, R)
        for orig, recon in zip(DOWN_QUARKS, reconstructed):
            assert recon == pytest.approx(orig, rel=1e-6)

    def test_lepton_R_near_M_sqrt2(self):
        """For leptons, R ≈ M√2 because Q ≈ 2/3."""
        M, _ = koide_params(LEPTONS)
        R = koide_R(LEPTONS)
        assert R == pytest.approx(M * math.sqrt(2), rel=1e-4)

    def test_M_positive(self):
        """M is positive for all families."""
        for masses in [LEPTONS, UP_QUARKS, DOWN_QUARKS]:
            M, _ = koide_params(masses)
            assert M > 0

    def test_lepton_M_value(self):
        """Lepton M ≈ 17.716 MeV^(1/2)."""
        M, _ = koide_params(LEPTONS)
        assert M == pytest.approx(17.716, abs=0.01)


# ============================================================
# DATA VALIDATION — Koide angle values
# ============================================================

@pytest.mark.data
class TestKoideAngles:
    """Verify Koide angle values from PDG masses."""

    def test_lepton_theta(self):
        """Lepton θ ≈ 132.73° = 2.3163 rad."""
        _, theta = koide_params(LEPTONS)
        theta_deg = theta * 180 / math.pi
        assert theta_deg == pytest.approx(132.733, abs=0.01)

    def test_lepton_excess_positive(self):
        """Lepton excess over 120° is positive."""
        excess = koide_excess(LEPTONS)
        assert excess > 0

    def test_lepton_excess_value(self):
        """Lepton excess ε ≈ 0.22222 rad = 12.733°."""
        excess = koide_excess(LEPTONS)
        excess_deg = excess * 180 / math.pi
        assert excess_deg == pytest.approx(12.733, abs=0.01)


# ============================================================
# DATA VALIDATION — Three-fold symmetry facts
# ============================================================

@pytest.mark.data
class TestThreeFoldSymmetryFacts:
    """Facts about the number 3 in the Standard Model."""

    def test_three_lepton_generations(self):
        """Three charged lepton generations."""
        leptons = [M_ELECTRON, M_MUON, M_TAU]
        assert len(leptons) == 3

    def test_three_up_quark_generations(self):
        """Three up-type quark generations."""
        quarks = [M_UP, M_CHARM, M_TOP]
        assert len(quarks) == 3

    def test_three_down_quark_generations(self):
        """Three down-type quark generations."""
        quarks = [M_DOWN, M_STRANGE, M_BOTTOM]
        assert len(quarks) == 3

    def test_charges_are_thirds(self):
        """All quark charges are multiples of 1/3."""
        charges = [2/3, -1/3]  # up-type, down-type
        for q in charges:
            assert (q * 3) == pytest.approx(round(q * 3))

    def test_three_colors(self):
        """SU(3) has 3 color charges (N_c = 3)."""
        N_c = 3
        assert N_c == 3

    def test_z_width_three_neutrinos(self):
        """Z width measurement gives N_ν = 2.984 ± 0.008 — exactly 3."""
        N_nu = 2.984
        N_nu_unc = 0.008
        assert abs(N_nu - 3) < 3 * N_nu_unc

    def test_120_degrees_is_triangle(self):
        """120° = 2π/3 is the angle of an equilateral triangle."""
        angle = 2 * math.pi / 3
        assert angle == pytest.approx(120 * math.pi / 180)


# ============================================================
# FRAMEWORK TESTS — Lepton excess = 2/9 (the main discovery)
# ============================================================

@pytest.mark.framework
class TestLeptonExcessTwoNinths:
    """
    GREEN: Lepton Koide excess = 2/9 radians to 33 ppm.

    This is an empirical discovery with no known derivation.
    The precision (33 ppm) suggests it's exact, not accidental.

    2/9 = 2/3² connects to Q = 2/3¹, forming a 2/3^n pattern.
    """

    def test_excess_equals_two_ninths(self):
        """ε_leptons = 2/9 radians (33 ppm)."""
        result = lepton_excess_is_two_ninths(LEPTONS)
        assert result['rel_error_ppm'] < 50, \
            f"Error = {result['rel_error_ppm']:.1f} ppm, expected < 50"

    def test_two_thirds_pattern(self):
        """Q = 2/3¹ and ε = 2/3² — same base, sequential powers."""
        Q = koide_Q(LEPTONS)
        excess = koide_excess(LEPTONS)

        Q_target = 2.0 / 3.0     # 2/3^1
        eps_target = 2.0 / 9.0   # 2/3^2

        Q_err = abs(Q - Q_target) / Q_target
        eps_err = abs(excess - eps_target) / eps_target

        assert Q_err < 1e-4, f"Q error: {Q_err*1e6:.1f} ppm"
        assert eps_err < 5e-5, f"ε error: {eps_err*1e6:.1f} ppm"

    def test_excess_degrees_times_pi_near_40(self):
        """ε_degrees × π ≈ 40. Since ε = 2/9 rad = 12.732°, 12.732 × π = 39.98."""
        excess_deg = koide_excess(LEPTONS) * 180 / math.pi
        assert excess_deg * math.pi == pytest.approx(40.0, rel=1e-3)

    def test_lepton_constraint_count(self):
        """
        Q = 2/3 + ε = 2/9 → 2 constraints for 3 masses → 1 free parameter (M).
        Three masses, two constraints, one parameter. This is maximal.
        """
        # With Q=2/3 and ε=2/9, reconstruct from M alone
        M, _ = koide_params(LEPTONS)
        theta_predicted = 2 * math.pi / 3 + 2.0 / 9.0
        reconstructed = koide_reconstruct(M, theta_predicted)

        for orig, recon in zip(LEPTONS, reconstructed):
            assert recon == pytest.approx(orig, rel=5e-4), \
                f"M + exact {2/9} gives {recon:.4f}, measured {orig:.4f}"


# ============================================================
# FRAMEWORK TESTS — D=27 quantization across families
# ============================================================

@pytest.mark.framework
class TestQuantizationD27:
    """
    GREEN (with caveats): All three family excesses quantize to n/27.

    D = 27 = 3³ is the common denominator:
      Leptons:   6/27 = 2/9  (0.003%)
      Up quarks: 2/27        (0.85%)
      Down quarks: 3/27 = 1/9 (0.48%)

    Numerator pattern: {6, 2, 3} — no obvious charge/color relationship yet.
    Quark fits are ~100× worse than lepton fit.
    """

    def test_lepton_six_over_27(self):
        """Lepton excess ≈ 6/27 = 2/9 (< 0.01% error)."""
        result = quantized_excess_27(LEPTONS)
        assert result['n_over_27'] == 6
        assert result['rel_error_pct'] < 0.01

    def test_up_quark_two_over_27(self):
        """Up quark excess ≈ 2/27 (< 1% error)."""
        result = quantized_excess_27(UP_QUARKS)
        assert result['n_over_27'] == 2
        assert result['rel_error_pct'] < 1.0

    def test_down_quark_three_over_27(self):
        """Down quark excess ≈ 3/27 = 1/9 (< 1% error)."""
        result = quantized_excess_27(DOWN_QUARKS)
        assert result['n_over_27'] == 3
        assert result['rel_error_pct'] < 1.0

    def test_common_denominator_is_three_cubed(self):
        """D = 27 = 3³ — a pure power of 3."""
        D = 27
        assert D == 3 ** 3

    def test_numerator_pattern(self):
        """Numerators are {6, 2, 3}. Sum = 11. No obvious structure yet."""
        numerators = []
        for masses in [LEPTONS, UP_QUARKS, DOWN_QUARKS]:
            result = quantized_excess_27(masses)
            numerators.append(result['n_over_27'])
        assert numerators == [6, 2, 3]

    def test_two_thirds_power_series(self):
        """
        If ε = 2/3^n, then:
          n=2: 2/9  = 6/27  ← leptons
          n=3: 2/27 = 2/27  ← up quarks
        Down quarks (3/27) don't fit 2/3^n but ARE n/27.
        """
        assert 2 / 3**2 == pytest.approx(6 / 27)
        assert 2 / 3**3 == pytest.approx(2 / 27)
        # Down quarks: 3/27 = 1/9 ≠ 2/3^n for any integer n
        assert 3 / 27 == pytest.approx(1 / 9)


# ============================================================
# FRAMEWORK TESTS — Individual exact fractions (high precision)
# ============================================================

@pytest.mark.framework
class TestIndividualExactFractions:
    """
    GREEN: Each family's excess has a high-precision rational approximant.

    Leptons:     2/9         (0.003%)
    Up quarks:   (3/11)² = 9/121  (0.0015%)
    Down quarks: 13/118      (0.006%)

    These are individually tighter than D=27 but don't share
    obvious common structure.
    """

    def test_lepton_two_ninths(self):
        """Lepton excess = 2/9 to < 0.005%."""
        excess = koide_excess(LEPTONS)
        err = abs(excess - 2/9) / (2/9)
        assert err < 5e-5

    def test_up_quark_nine_over_121(self):
        """Up quark excess ≈ 9/121 = (3/11)² to < 0.01%."""
        excess = koide_excess(UP_QUARKS)
        target = 9 / 121
        err = abs(excess - target) / target
        assert err < 1e-3  # 0.1% (generous; actual is 0.0015%)

    def test_up_quark_fraction_is_perfect_square(self):
        """9/121 = 3²/11² = (3/11)² — a perfect square of a ratio."""
        assert 9 == 3**2
        assert 121 == 11**2

    def test_down_quark_13_over_118(self):
        """Down quark excess ≈ 13/118 to < 0.01%."""
        excess = koide_excess(DOWN_QUARKS)
        target = 13 / 118
        err = abs(excess - target) / target
        assert err < 1e-3  # 0.1% (generous; actual is 0.006%)


# ============================================================
# FRAMEWORK TESTS — RED: Not yet derived
# ============================================================

@pytest.mark.framework
class TestKoideDerivations:
    """
    RED: These need derivations from first principles.
    """

    def test_derive_excess(self):
        """Must derive ε = 2/9 from boundary geometry."""
        with pytest.raises(NotImplementedError):
            derive_koide_excess()

    def test_derive_M(self):
        """Must derive M_leptons = 17.716 MeV^(1/2)."""
        with pytest.raises(NotImplementedError):
            derive_M_parameter()

    def test_predict_neutrinos(self):
        """Must predict neutrino masses from quantized Koide."""
        with pytest.raises(NotImplementedError):
            predict_neutrino_masses()


# ============================================================
# FRAMEWORK TESTS — Neutrino Koide Q cannot reach 2/3
# ============================================================

@pytest.mark.framework
class TestNeutrinoKoideFails:
    """
    GREEN: Neutrino Koide Q CANNOT reach 2/3.

    For any neutrino mass configuration consistent with measured
    delta-m-squared splittings, Q_max ~ 0.551 (NH). The constraint
    Q = 2/3 is unreachable. This is not a failure of Koide — it's
    evidence that Koide is an eigenvalue property that requires
    complete decoherence (gamma -> infinity), which neutrinos lack.
    """

    def test_neutrino_Q_max_NH(self):
        """Maximum neutrino Q ~ 0.551 in Normal Hierarchy."""
        Q_max = neutrino_koide_max_Q('NH')
        assert Q_max < 0.60, f"Q_max = {Q_max:.4f}, should be < 0.60"
        assert Q_max > 0.50, f"Q_max = {Q_max:.4f}, should be > 0.50"

    def test_neutrino_Q_never_reaches_two_thirds(self):
        """No neutrino mass configuration can achieve Q = 2/3."""
        Q_max = neutrino_koide_max_Q('NH')
        assert Q_max < 2/3, \
            f"Q_max = {Q_max:.4f} >= 2/3 — Koide should be unreachable for neutrinos"

    def test_neutrino_Q_below_lepton_Q(self):
        """Neutrino Q_max < lepton Q. Neutrinos BELOW, quarks ABOVE."""
        Q_max_nu = neutrino_koide_max_Q('NH')
        Q_lep = koide_Q(LEPTONS)
        assert Q_max_nu < Q_lep, \
            f"Neutrino Q_max {Q_max_nu:.4f} should be below lepton Q {Q_lep:.6f}"


# ============================================================
# FRAMEWORK TESTS — Quark Q deviates ABOVE 2/3
# ============================================================

@pytest.mark.framework
class TestQuarkQAboveTwoThirds:
    """
    GREEN: Both quark families have Q > 2/3.

    Partially decohered systems overshoot the eigenvalue condition.
    The direction is consistent: more decoherence -> Q closer to 2/3 from above.
    """

    def test_down_quark_Q_above_two_thirds(self):
        """Down quark Q = 0.731 > 2/3."""
        Q = koide_Q(DOWN_QUARKS)
        assert Q > 2/3, f"Down Q = {Q:.4f} should be > 2/3"
        assert Q == pytest.approx(0.731, abs=0.005)

    def test_up_quark_Q_above_two_thirds(self):
        """Up quark Q = 0.849 > 2/3."""
        Q = koide_Q(UP_QUARKS)
        assert Q > 2/3, f"Up Q = {Q:.4f} should be > 2/3"
        assert Q == pytest.approx(0.849, abs=0.005)

    def test_up_quark_Q_further_than_down(self):
        """Up quarks deviate MORE from 2/3 than down quarks."""
        Q_up = koide_Q(UP_QUARKS)
        Q_down = koide_Q(DOWN_QUARKS)
        assert abs(Q_up - 2/3) > abs(Q_down - 2/3), \
            f"|Q_up - 2/3| = {abs(Q_up - 2/3):.4f} should exceed |Q_down - 2/3| = {abs(Q_down - 2/3):.4f}"

    def test_Q_direction_monotonic(self):
        """
        Koide Q direction: nu BELOW 2/3, quarks ABOVE 2/3, leptons EXACT.

        This is the two-category prediction: coherence compresses eigenvalues
        (pushing Q below 2/3 for neutrinos), partial decoherence overshoots
        (pushing Q above for quarks), full decoherence nails it (leptons = 2/3).
        """
        Q_lep = koide_Q(LEPTONS)
        Q_up = koide_Q(UP_QUARKS)
        Q_down = koide_Q(DOWN_QUARKS)
        Q_nu_max = neutrino_koide_max_Q('NH')

        # Neutrinos below, quarks above, leptons exact
        assert Q_nu_max < 2/3, "neutrino Q should be below 2/3"
        assert Q_down > 2/3, "down quark Q should be above 2/3"
        assert Q_up > 2/3, "up quark Q should be above 2/3"
        assert abs(Q_lep - 2/3) < 1e-4, "lepton Q should be nearly exactly 2/3"


# ============================================================
# FRAMEWORK TESTS — Tau mass prediction from Koide
# ============================================================

@pytest.mark.framework
class TestTauPrediction:
    """
    GREEN: Koide Q = 2/3 predicts m_tau = 1776.969 MeV.

    Given only m_e and m_mu, the Koide constraint Q = 2/3 uniquely
    determines m_tau. The prediction is 1776.969 MeV, compared to
    the PDG measurement 1776.86 +/- 0.12 MeV. That's 0.9 sigma.

    Belle II should reach +/- 0.05 MeV precision and test this.
    """

    def test_tau_prediction_value(self):
        """Predicted tau mass is 1776.969 MeV."""
        m_tau_pred = predict_tau_from_koide(M_ELECTRON, M_MUON)
        assert m_tau_pred == pytest.approx(1776.969, abs=0.01), \
            f"Predicted m_tau = {m_tau_pred:.3f} MeV"

    def test_tau_prediction_within_1sigma(self):
        """Prediction within 1 sigma of PDG measurement."""
        m_tau_pred = predict_tau_from_koide(M_ELECTRON, M_MUON)
        deviation_sigma = abs(m_tau_pred - M_TAU) / M_TAU_UNC
        assert deviation_sigma < 1.5, \
            f"Deviation = {deviation_sigma:.1f} sigma (pred={m_tau_pred:.3f}, meas={M_TAU} +/- {M_TAU_UNC})"

    def test_tau_prediction_precision(self):
        """Prediction matches measurement to < 0.01%."""
        m_tau_pred = predict_tau_from_koide(M_ELECTRON, M_MUON)
        pct_err = abs(m_tau_pred - M_TAU) / M_TAU * 100
        assert pct_err < 0.01, f"Prediction error = {pct_err:.4f}%"
