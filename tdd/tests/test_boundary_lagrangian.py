"""
Test Domain 12: Boundary Lagrangian and z-Tower (Paper 4)

The standing wave model expressed as a 1D boundary field theory.
z = √m connects boundary amplitude → mass → BH entropy.
Bekenstein bound gives S ≤ 2π at Compton radius for all particles.

Key findings:
  - z-tower: z (amplitude), z² (mass), z⁴ (BH entropy) — one variable
  - Bekenstein bound at Compton: S = 2π, mass-independent
  - Overlap matrix M_{nm} ∝ n/(n²-m²) for n+m odd, 0 otherwise
  - Scale correspondences: Ωb ≈ (sin²θ_W)², ΩΛ ≈ 2/3 (empirical)

Scorecard:
  GREEN (data)      -- Bekenstein bound, z-tower algebra
  GREEN (framework) -- overlap matrix, mixing from z-ratio
  YELLOW (empirical) -- cosmological correspondences (observations, not predictions)
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tests.conftest import *
from framework.decoherence import (
    z_amplitude, z_tower, bekenstein_at_compton, bh_entropy,
    z_tower_ratio, mixing_from_z_tower, overlap_matrix,
    scale_correspondence, chirality_overlap,
    gamma_from_color_group, cabibbo_from_qcd,
    gauge_poles_from_physics,
    pmns_cp_phase_analysis, pmns_theta23_correction,
)


# ============================================================
# DATA VALIDATION — Bekenstein bound
# ============================================================

@pytest.mark.data
class TestBekensteinBound:
    """S ≤ 2πRE/(ℏc) = 2π at Compton radius, mass-independent."""

    def test_electron_2pi(self):
        """Electron: S_bek = 2π at Compton radius."""
        result = bekenstein_at_compton(M_ELECTRON)
        assert result['S_bekenstein'] == pytest.approx(2 * math.pi, rel=1e-6)

    def test_muon_2pi(self):
        """Muon: S_bek = 2π at Compton radius."""
        result = bekenstein_at_compton(M_MUON)
        assert result['S_bekenstein'] == pytest.approx(2 * math.pi, rel=1e-6)

    def test_top_2pi(self):
        """Top quark: S_bek = 2π at Compton radius."""
        result = bekenstein_at_compton(M_TOP)
        assert result['S_bekenstein'] == pytest.approx(2 * math.pi, rel=1e-6)

    def test_mass_independent(self):
        """Bekenstein bound is the same for all particles."""
        masses = [M_ELECTRON, M_MUON, M_TAU, M_UP, M_DOWN,
                  M_STRANGE, M_CHARM, M_BOTTOM, M_TOP]
        values = [bekenstein_at_compton(m)['S_bekenstein'] for m in masses]
        for v in values:
            assert v == pytest.approx(2 * math.pi, rel=1e-6)

    def test_2pi_is_one_phase_cycle(self):
        """2π = one complete phase cycle = circumference of unit circle."""
        result = bekenstein_at_compton(M_ELECTRON)
        assert result['target'] == pytest.approx(2 * math.pi)


# ============================================================
# FRAMEWORK — z-tower
# ============================================================

@pytest.mark.framework
class TestZTower:
    """z = √m, z² = m, z⁴ ∝ S_BH — one variable, four levels."""

    def test_z_is_sqrt_m(self):
        """z = √m for all particles."""
        for m in [M_ELECTRON, M_MUON, M_TAU, M_UP, M_CHARM, M_TOP]:
            assert z_amplitude(m) == pytest.approx(math.sqrt(m))

    def test_z_squared_is_mass(self):
        """z² recovers the mass."""
        for m in [M_ELECTRON, M_CHARM, M_BOTTOM]:
            tower = z_tower(m)
            assert tower['z2'] == pytest.approx(m)
            assert tower['z'] ** 2 == pytest.approx(m)

    def test_z_fourth_is_mass_squared(self):
        """z⁴ = m² — proportional to BH entropy."""
        for m in [M_ELECTRON, M_TAU, M_TOP]:
            tower = z_tower(m)
            assert tower['z4'] == pytest.approx(m ** 2)
            assert tower['z'] ** 4 == pytest.approx(m ** 2)

    def test_z0_is_unity(self):
        """z⁰ = 1 for all particles (universal level)."""
        for m in [M_ELECTRON, M_MUON, M_TAU]:
            tower = z_tower(m)
            assert tower['z0'] == 1.0

    def test_bh_entropy_scales_as_m_squared(self):
        """S_BH(m₂)/S_BH(m₁) = (m₂/m₁)²."""
        S_e = bh_entropy(M_ELECTRON)
        S_mu = bh_entropy(M_MUON)
        S_tau = bh_entropy(M_TAU)
        assert S_mu / S_e == pytest.approx((M_MUON / M_ELECTRON) ** 2, rel=1e-4)
        assert S_tau / S_e == pytest.approx((M_TAU / M_ELECTRON) ** 2, rel=1e-4)


# ============================================================
# FRAMEWORK — z-tower ratios and mixing
# ============================================================

@pytest.mark.framework
class TestZTowerRatios:
    """Mixing angles are z-ratios; mass ratios are z²-ratios."""

    def test_gst_from_z_ratio(self):
        """V_us = z_d/z_s = √(m_d/m_s) — GST relation from z-tower."""
        v_us = mixing_from_z_tower(M_DOWN, M_STRANGE)
        assert v_us == pytest.approx(0.2236, abs=0.001)

    def test_vcb_from_z_ratio(self):
        """V_cb = z_u/z_c = √(m_u/m_c) — from z-tower."""
        v_cb = mixing_from_z_tower(M_UP, M_CHARM)
        assert v_cb == pytest.approx(0.0412, abs=0.002)

    def test_z_ratio_level0_always_one(self):
        """At level 0, all ratios = 1 (particles equal at own scale)."""
        ratio = z_tower_ratio(M_ELECTRON, M_TOP)
        assert ratio['z0_ratio'] == 1.0

    def test_z_ratio_level2_is_mass_ratio(self):
        """At level 2, ratio = m₁/m₂."""
        ratio = z_tower_ratio(M_DOWN, M_STRANGE)
        assert ratio['z2_ratio'] == pytest.approx(M_DOWN / M_STRANGE)

    def test_z_ratio_level4_is_mass_ratio_squared(self):
        """At level 4, ratio = (m₁/m₂)²."""
        ratio = z_tower_ratio(M_DOWN, M_STRANGE)
        assert ratio['z4_ratio'] == pytest.approx((M_DOWN / M_STRANGE) ** 2)


# ============================================================
# FRAMEWORK — Overlap matrix (boundary Lagrangian mass matrix)
# ============================================================

@pytest.mark.framework
class TestOverlapMatrix:
    """M_{nm} = ∫₀¹ sin(nπx) cos(mπx) dx — the mass matrix."""

    def test_diagonal_is_zero(self):
        """M_{nn} = 0 for all n (no self-coupling at same mode)."""
        M = overlap_matrix(7, 4)
        for n in range(4):
            assert M[n][n] == pytest.approx(0.0)

    def test_even_sum_is_zero(self):
        """M_{nm} = 0 when n+m is even."""
        M = overlap_matrix(7, 4)
        # n=1,m=1 (sum=2, even): 0
        assert M[0][0] == pytest.approx(0.0)
        # n=1,m=3 (sum=4, even): 0
        assert M[0][2] == pytest.approx(0.0)
        # n=2,m=2 (sum=4, even): 0
        assert M[1][1] == pytest.approx(0.0)

    def test_odd_sum_nonzero(self):
        """M_{nm} ≠ 0 when n+m is odd and n≠m."""
        M = overlap_matrix(7, 4)
        # n=1,m=2 (sum=3, odd): nonzero
        assert abs(M[0][1]) > 0.01
        # n=2,m=1 (sum=3, odd): nonzero
        assert abs(M[1][0]) > 0.01

    def test_hierarchy_from_mode_number(self):
        """Higher mode differences give smaller matrix elements."""
        M = overlap_matrix(7, 6)
        # |M_{12}| > |M_{14}| (closer modes → larger coupling)
        assert abs(M[0][1]) > abs(M[0][3])

    def test_analytic_formula(self):
        """M_{nm} = 2n/(π(n²-m²)) for n+m odd."""
        M = overlap_matrix(7, 4)
        # n=1, m=2: 2×1/(π×(1-4)) = -2/(3π) = -0.2122
        expected = 2 * 1 / (math.pi * (1 - 4))
        assert M[0][1] == pytest.approx(expected, abs=0.001)

    def test_matches_numerical_integration(self):
        """Analytic overlap matches numerical chirality_overlap."""
        M = overlap_matrix(7, 4)
        # Compare M[0][1] (n=1,m=2) with chirality_overlap(1,2,'sin','cos')
        numerical = chirality_overlap(1, 2, 'sin', 'cos')
        assert M[0][1] == pytest.approx(numerical, abs=0.002)


# ============================================================
# FRAMEWORK — Scale correspondence (empirical, YELLOW)
# ============================================================

@pytest.mark.framework
class TestScaleCorrespondence:
    """Particle↔cosmological dimensionless number matches.
    These are OBSERVATIONS, not predictions. Tolerances are generous."""

    def test_omega_b_sin2w_squared(self):
        """Ωb ≈ (sin²θ_W)² within 1%."""
        matches = scale_correspondence()
        assert matches['omega_b_sin2w']['match_pct'] < 1.0

    def test_omega_b_vus_squared(self):
        """Ωb ≈ V_us² = m_d/m_s within 2%."""
        matches = scale_correspondence()
        assert matches['omega_b_vus']['match_pct'] < 2.0

    def test_omega_lambda_near_two_thirds(self):
        """ΩΛ ≈ 2/3 within 3%."""
        matches = scale_correspondence()
        assert matches['omega_lambda_2_3']['match_pct'] < 3.0

    def test_f_boost_near_sqrt3(self):
        """f_boost ≈ √3 within 1%."""
        matches = scale_correspondence()
        assert matches['f_boost_sqrt3']['match_pct'] < 1.0

    def test_all_matches_below_10_percent(self):
        """All cataloged correspondences are within 10%."""
        matches = scale_correspondence()
        for key, m in matches.items():
            assert m['match_pct'] < 10.0, f"{key} exceeds 10%"


# ============================================================
# DATA VALIDATION — γ = N_c × T_F (the "Maxwell step")
# ============================================================

@pytest.mark.data
class TestGammaDerivation:
    """γ = N_c × T_F = 3/2: Cabibbo angle from QCD group theory."""

    def test_gamma_equals_3_over_2(self):
        """γ = 3 × 1/2 = 3/2 within 0.5% of measured."""
        result = gamma_from_color_group(N_c=3)
        assert result['gamma_predicted'] == 1.5
        assert result['gamma_error_pct'] < 0.5

    def test_lambda_from_qcd(self):
        """λ = e^{-3/2} = 0.2231 within 1% of V_us."""
        result = gamma_from_color_group()
        assert result['lambda_error_pct'] < 1.0

    def test_cabibbo_from_nc(self):
        """V_us = e^{-N_c/2} within 1%."""
        V_us_pred = cabibbo_from_qcd()
        V_us_meas = 0.22453
        assert abs(V_us_pred - V_us_meas) / V_us_meas < 0.01

    def test_ms_md_equals_e_cubed(self):
        """m_s/m_d = e³ = 20.09 within 0.5%."""
        result = gamma_from_color_group()
        assert result['ms_md_error_pct'] < 0.5

    def test_tf_is_half(self):
        """T_F = 1/2 is the SU(N) fundamental index."""
        result = gamma_from_color_group()
        assert result['T_F'] == 0.5

    def test_gamma_scales_with_nc(self):
        """γ = N_c/2 for all SU(N)."""
        for Nc in [2, 3, 4, 5]:
            result = gamma_from_color_group(Nc)
            assert result['gamma_predicted'] == Nc / 2.0


# ============================================================
# FRAMEWORK — P from gauge theory
# ============================================================

@pytest.mark.framework
class TestGaugePoles:
    """P = number of decoherence channels from gauge structure."""

    def test_quark_P_from_color(self):
        """Quarks: P = N_c = 3."""
        info = gauge_poles_from_physics('up')
        assert info['P'] == 3
        assert info['N'] == 7

    def test_lepton_P_from_em(self):
        """Charged leptons: P = |Q_EM| = 1."""
        info = gauge_poles_from_physics('lepton')
        assert info['P'] == 1
        assert info['N'] == 3

    def test_neutrino_P_zero(self):
        """Neutrinos: P = 0 (no gauge charge)."""
        info = gauge_poles_from_physics('neutrino')
        assert info['P'] == 0
        assert info['N'] == 1

    def test_three_generations_for_charged(self):
        """All charged sectors give 3 generations."""
        for sector in ['up', 'down', 'lepton']:
            info = gauge_poles_from_physics(sector)
            assert info['n_generations'] == 3

    def test_two_generations_for_neutrinos(self):
        """Neutrinos give 2 (degenerate) generations."""
        info = gauge_poles_from_physics('neutrino')
        assert info['n_generations'] == 2


# ============================================================
# FRAMEWORK — PMNS CP phase and θ₂₃ correction
# ============================================================

@pytest.mark.framework
class TestPMNSAnalysis:
    """PMNS CP phase and θ₂₃ second-order correction."""

    def test_pmns_cp_near_pi_over_2(self):
        """arctan(m_e/m_ν₁) ≈ π/2 (weak prediction)."""
        result = pmns_cp_phase_analysis()
        assert result['near_pi_over_2']

    def test_pmns_cp_mass_ratio_huge(self):
        """m_e/m_ν₁ >> 1 (ratio ~10⁷)."""
        result = pmns_cp_phase_analysis()
        assert result['mass_ratio_e_nu1'] > 1e6

    def test_theta23_deficit_positive(self):
        """Measured sin²θ₂₃ > 1/2 (needs positive correction)."""
        result = pmns_theta23_correction()
        assert result['deficit'] > 0

    def test_theta23_linear_correction_direction(self):
        """1/2 + 3δx correction goes in the right direction."""
        result = pmns_theta23_correction()
        pred = result['corrections']['linear_3dx']['predicted']
        assert pred > 0.5  # correction is positive

    def test_theta23_zeroth_order_8pct(self):
        """Tribimaximal 1/2 is ~8% off measured."""
        result = pmns_theta23_correction()
        err = abs(result['corrections']['zeroth_order']['error_pct'])
        assert 5 < err < 15
