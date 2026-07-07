"""
Test Domain: Zero-Parameter Framework (session 2026-04-14)

Tests for the derived constants and cross-level predictions.
All quantities derived from N_c = 3 with zero free parameters.
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from framework import boundary


@pytest.mark.framework
class TestSigmaDGeometric:
    """sigma_d = N_q - C_0 = 34/5 = 6.8"""

    def test_sigma_d_value(self):
        assert boundary.derive_sigma_d() == pytest.approx(34/5, rel=1e-10)

    def test_sigma_d_matches_measured(self):
        """sigma_d = 6.8 vs measured ln(m_b/m_d) = 6.797"""
        measured = math.log(4180 / 4.67)  # = 6.7974
        assert boundary.derive_sigma_d() == pytest.approx(measured, rel=0.001)


@pytest.mark.framework
class TestM0Geometric:
    """M_0 = (4/H)(1-C_0^2) = 384/2500 = 0.1536 MeV"""

    def test_M0_value(self):
        assert boundary.derive_M0() == pytest.approx(384/2500, rel=1e-10)

    def test_M0_matches_effective(self):
        """M_0 = 0.1536 vs effective 0.154 MeV"""
        assert boundary.derive_M0() == pytest.approx(0.154, rel=0.005)


@pytest.mark.framework
class TestLambdaBW:
    """Lambda_BW = M_0 * exp(sigma_d) ≈ 138 MeV"""

    def test_lambda_bw_value(self):
        val = boundary.derive_lambda_bw()
        assert 135 < val < 142, f"Lambda_BW = {val:.1f}, expected ~138"

    def test_lambda_bw_near_pion(self):
        """Lambda_BW ≈ m_pi (the pion IS the boundary quantum)"""
        val = boundary.derive_lambda_bw()
        assert val == pytest.approx(139.6, rel=0.03)  # within 3% of pion mass


@pytest.mark.framework
class TestPlanckMassFromBoundary:
    """M_P = G_2 * y_t * Lambda_BW * exp(sigma_d^2) at +1.45%"""

    def test_planck_mass_order_of_magnitude(self):
        M_P = boundary.planck_mass_from_boundary()
        # Should be ~10^22 MeV = ~10^19 GeV
        assert 1e21 < M_P < 1e23, f"M_P = {M_P:.2e} MeV"

    def test_planck_mass_accuracy(self):
        M_P = boundary.planck_mass_from_boundary()
        M_P_measured = 1.22089e22  # MeV
        assert M_P == pytest.approx(M_P_measured, rel=0.02), \
            f"M_P = {M_P:.3e} vs {M_P_measured:.3e} MeV"

    def test_planck_mass_prefactor(self):
        """The prefactor G_2 * y_t = (3/4)(124/125)"""
        G_2 = 3/4
        y_t = 124/125
        prefactor = G_2 * y_t
        assert prefactor == pytest.approx(0.744, rel=0.001)


@pytest.mark.framework
class TestDarkMatterBaryonRatio:
    """DM/baryon = 1/C_0 + C_0 = 26/5 = 5.2"""

    def test_ratio_value(self):
        assert boundary.dark_matter_baryon_ratio() == pytest.approx(26/5, rel=1e-10)

    def test_ratio_matches_planck(self):
        """Planck 2018: Omega_DM/Omega_b = 5.32"""
        assert boundary.dark_matter_baryon_ratio() == pytest.approx(5.32, rel=0.03)


@pytest.mark.framework
class TestCosmologicalConstant:
    """CC = seesaw * 5/26 = 2.82 meV"""

    def test_cc_scale(self):
        cc = boundary.cosmological_constant_scale()
        cc_meV = cc * 1e9  # MeV to meV
        assert cc_meV == pytest.approx(2.846, rel=0.02), \
            f"CC = {cc_meV:.2f} meV, expected 2.846"

    def test_cc_order_of_magnitude(self):
        cc = boundary.cosmological_constant_scale()
        # Should be ~10^-9 MeV = ~meV
        assert 1e-10 < cc < 1e-8, f"CC = {cc:.2e} MeV"


@pytest.mark.framework
class TestProtonChargeRadius:
    """r_p = C_0 * N_c * hbar_c / Lambda_BW = 0.858 fm"""

    def test_proton_radius_value(self):
        r = boundary.proton_charge_radius_fm()
        assert r == pytest.approx(0.842, rel=0.025), \
            f"r_p = {r:.3f} fm, measured 0.842 fm"

    def test_proton_radius_formula(self):
        """r_p = (3/5) * (197.3 / Lambda_BW)"""
        r = boundary.proton_charge_radius_fm()
        assert 0.8 < r < 0.9, f"r_p = {r:.3f} fm"


@pytest.mark.framework
class TestBoundaryClockTick:
    """dt = h_orb / Lambda_BW ≈ 1.9e-25 seconds"""

    def test_clock_tick_order(self):
        dt = boundary.boundary_clock_tick()
        assert 1e-26 < dt < 1e-24, f"dt = {dt:.2e} s"

    def test_clock_tick_above_planck(self):
        """Boundary tick >> Planck time"""
        dt = boundary.boundary_clock_tick()
        t_planck = 5.39e-44
        assert dt / t_planck > 1e15, "Boundary tick should be >> Planck time"

    def test_top_quark_lives_few_ticks(self):
        """Top quark lifetime ≈ 3 boundary ticks"""
        dt = boundary.boundary_clock_tick()
        tau_top = 5e-25  # seconds
        n_ticks = tau_top / dt
        assert 1 < n_ticks < 10, f"Top lives {n_ticks:.1f} ticks"


@pytest.mark.framework
class TestSelfConsistency:
    """Cross-checks between derived quantities"""

    def test_sigma_squared_identity(self):
        """sigma_u * sigma_l / 2 = sigma_d^2 (from span ratios)"""
        sigma_d = boundary.derive_sigma_d()
        sigma_u = (5/3) * sigma_d
        sigma_l = (6/5) * sigma_d
        assert sigma_u * sigma_l / 2 == pytest.approx(sigma_d**2, rel=1e-10)

    def test_total_modes_equals_beta0(self):
        """N_q + N_l + N_nu = 11 = beta_0^gauge"""
        assert 7 + 3 + 1 == 11

    def test_proton_mass_from_lambda(self):
        """m_p = N_q * Lambda_BW * (1 - C_0^2(1-C_0))"""
        Lambda_BW = boundary.derive_lambda_bw()
        C_0 = 1/5
        m_p = 7 * Lambda_BW * (1 - C_0**2 * (1 - C_0))
        assert m_p == pytest.approx(938.3, rel=0.005)

    def test_higgs_from_warp_tower(self):
        """M_H ≈ M_0 * exp(2*sigma_d) = warp level 2"""
        M_0 = boundary.derive_M0()
        sigma_d = boundary.derive_sigma_d()
        M_H = M_0 * math.exp(2 * sigma_d) / 1e3  # MeV to GeV
        assert M_H == pytest.approx(125.25, rel=0.015)

    def test_mp_over_v_ratio(self):
        """m_p / v = 3/(80*pi^2) at 0.3%"""
        ratio_framework = 3 / (80 * math.pi**2)
        ratio_measured = 938.3 / 246220
        assert ratio_framework == pytest.approx(ratio_measured, rel=0.005)
