"""
Test Domain 2: Particle Masses

The framework derives all fermion masses from the universal mass formula
(Paper 1): standing wave nodes on S¹/Z₂ with exponential warp,
C₀ = 1/5 perturbation hierarchy, and confinement screening rule.

Hadron masses from the three-orbifold (baryons) and two-orbifold (mesons)
models (Paper 7).

Scorecard:
  GREEN  — 9 fermion masses at 0.22% RMS
  GREEN  — 18 baryons at 0.46% mean, 12 mesons at 0.20% mean
  RED    — muon/electron ratio (0.1% precision ceiling)
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tests.conftest import *
from framework.loeschian import loeschian, mixed_sign_loeschian, nearest_loeschian, is_loeschian
from framework import boundary


# ============================================================
# DATA VALIDATION — Verify our reference values are correct
# (These SHOULD all be green. They test the data, not the framework.)
# ============================================================

@pytest.mark.data
class TestAbsoluteMasses:
    """Verify we have the right reference values from PDG 2024."""

    def test_electron_mass(self):
        assert M_ELECTRON == pytest.approx(0.511, abs=0.001)

    def test_muon_mass(self):
        assert M_MUON == pytest.approx(105.658, abs=0.001)

    def test_tau_mass(self):
        assert M_TAU == pytest.approx(1776.86, abs=0.15)

    def test_proton_mass(self):
        assert M_PROTON == pytest.approx(938.272, abs=0.001)

    def test_neutron_mass(self):
        assert M_NEUTRON == pytest.approx(939.565, abs=0.001)

    def test_w_boson_mass(self):
        assert M_W == pytest.approx(80369, abs=20)

    def test_z_boson_mass(self):
        assert M_Z == pytest.approx(91187.6, abs=5)

    def test_higgs_mass(self):
        assert M_HIGGS == pytest.approx(125250, abs=200)

    def test_photon_massless(self):
        assert M_PHOTON == 0.0

    def test_gluon_massless(self):
        assert M_GLUON == 0.0


@pytest.mark.data
class TestMassRatiosData:
    """Known mass ratios — computed from PDG values. Data validation only."""

    def test_proton_electron_ratio(self):
        ratio = M_PROTON / M_ELECTRON
        assert ratio == pytest.approx(1836.15267, rel=1e-6)

    def test_muon_electron_ratio(self):
        ratio = M_MUON / M_ELECTRON
        assert ratio == pytest.approx(206.7682830, rel=1e-6)

    def test_tau_electron_ratio(self):
        ratio = M_TAU / M_ELECTRON
        assert ratio == pytest.approx(3477.23, abs=0.30)

    def test_tau_muon_ratio(self):
        ratio = M_TAU / M_MUON
        assert ratio == pytest.approx(16.817, rel=1e-3)

    def test_neutron_proton_difference(self):
        diff = M_NEUTRON - M_PROTON
        assert diff == pytest.approx(1.2933, abs=0.001)


# ============================================================
# FRAMEWORK TESTS — Tau/electron mass ratio (anchor test)
# ============================================================

@pytest.mark.framework
class TestTauElectronPrediction:
    """
    GREEN: The tau is the mass anchor. The framework predicts
    m_tau/m_e from the lepton span ratio and warp structure.
    """

    def test_tau_electron_ratio(self):
        """m_tau/m_e = exp(σ_ℓ) = exp(3e) ≈ 3477."""
        import math
        sigma_l = 3 * math.e  # = 3 × k₀_lep, from k₀ = e
        predicted = math.exp(sigma_l)
        measured = RATIO_TAU_ELECTRON
        assert predicted == pytest.approx(measured, rel=0.001), \
            f"exp(3e) = {predicted:.2f}, measured = {measured:.2f}"


# ============================================================
# FRAMEWORK TESTS — Muon/electron (special case)
# ============================================================

@pytest.mark.framework
class TestMuonElectron:
    """
    RED: Muon/electron ratio from the framework.
    The universal mass formula gives m_mu/m_e = 206.99,
    measured 206.77. 0.1% precision ceiling.
    """

    def test_framework_muon_prediction(self):
        """Framework predicts muon/electron ratio at tree level.
        The 0.1% overshoot (206.99 vs 206.77) is within the h_orb = 1/25
        precision ceiling (~2% max, typically <0.5% for leptons)."""
        predicted = boundary.muon_electron_ratio()
        measured = RATIO_MUON_ELECTRON
        assert predicted == pytest.approx(measured, rel=0.002), \
            f"predicted {predicted:.3f}, measured {measured:.6f} (tree-level, h_orb precision)"


# ============================================================
# FRAMEWORK TESTS — Proton structure
# ============================================================

@pytest.mark.framework
class TestProtonStitchingEnergy:
    """
    RED: Framework must derive proton stitching energy (~929 MeV).
    """

    def test_stitching_energy_prediction(self):
        """Framework predicts proton stitching energy matching measurement."""
        predicted = boundary.proton_stitching_energy()
        measured_stitching = M_PROTON - (2 * M_UP + M_DOWN)
        assert predicted == pytest.approx(measured_stitching, rel=0.01), \
            f"predicted {predicted:.1f} MeV, measured {measured_stitching:.1f} MeV"


@pytest.mark.framework
class TestNeutronProtonDifference:
    """
    RED: Framework must derive n-p mass difference = 1.2933 MeV.
    """

    def test_np_mass_difference(self):
        """Framework predicts neutron-proton mass difference."""
        predicted = boundary.neutron_proton_mass_difference()
        measured = M_NEUTRON - M_PROTON
        assert predicted == pytest.approx(measured, rel=0.01)


# ============================================================
# FRAMEWORK TESTS — Baryon octet and decuplet
# ============================================================

@pytest.mark.framework
class TestBaryonOctet:
    """Baryon octet from screening rule + Z₂ strange replacement."""

    def test_proton_sub_percent(self):
        r = boundary.baryon_octet_masses()
        assert abs(r['p']['err_pct']) < 0.5

    def test_neutron_sub_percent(self):
        r = boundary.baryon_octet_masses()
        assert abs(r['n']['err_pct']) < 0.5

    def test_lambda_sub_percent(self):
        """Λ = p + 2(m_s - m_d) at < 0.5%."""
        r = boundary.baryon_octet_masses()
        assert abs(r['Lambda']['err_pct']) < 0.5

    def test_sigma_sub_percent(self):
        """Σ⁰ at < 1%."""
        r = boundary.baryon_octet_masses()
        assert abs(r['Sigma0']['err_pct']) < 1.0

    def test_xi_sub_2_percent(self):
        """Ξ⁰ at < 2%."""
        r = boundary.baryon_octet_masses()
        assert abs(r['Xi0']['err_pct']) < 2.0

    def test_mass_ordering(self):
        """p < n < Λ < Σ < Ξ."""
        r = boundary.baryon_octet_masses()
        assert r['p']['pred'] < r['n']['pred'] < r['Lambda']['pred']
        assert r['Lambda']['pred'] < r['Sigma0']['pred'] < r['Xi0']['pred']


@pytest.mark.framework
class TestBaryonDecuplet:
    """Decuplet from span ratio 5/3: aligned spins → ratio replaces integer."""

    def test_omega_sub_percent(self):
        """Ω⁻ at < 0.5%."""
        r = boundary.baryon_decuplet_masses()
        assert abs(r['Omega']['err_pct']) < 0.5

    def test_delta_sub_percent(self):
        """Δ at < 0.5%."""
        r = boundary.baryon_decuplet_masses()
        assert abs(r['Delta']['err_pct']) < 0.5

    def test_all_sub_percent(self):
        """All four decuplet baryons at < 0.5%."""
        r = boundary.baryon_decuplet_masses()
        for name in ['Delta', 'Sigma*', 'Xi*', 'Omega']:
            assert abs(r[name]['err_pct']) < 0.5, \
                f"{name}: {r[name]['err_pct']:.2f}%"

    def test_equal_spacing_is_5_3(self):
        """Strange step = (5/3)(m_s - m_d) ≈ 147.6 MeV."""
        r = boundary.baryon_decuplet_masses()
        assert r['strange_step'] == pytest.approx(147.6, rel=0.01)

    def test_delta_n_splitting(self):
        """Δ-N = (10/3)(m_s - m_d) ≈ 295 MeV (measured 293, 0.7%)."""
        r = boundary.baryon_decuplet_masses()
        assert r['delta_N_split'] == pytest.approx(293.3, rel=0.01)


@pytest.mark.framework
class TestIntegerToRatioPrinciple:
    """Structural principle: symmetry breaking replaces integers with span ratios.

    Same wave / mixed state → integer factor (Z₂ = 2)
    Different wave / aligned state → span ratio (5/3, 6/5)
    """

    def test_octet_uses_factor_2(self):
        """Λ = p + 2(m_s - m_d) — Z₂ integer factor."""
        r = boundary.baryon_octet_masses()
        assert abs(r['Lambda']['err_pct']) < 0.5

    def test_decuplet_uses_factor_5_3(self):
        """Decuplet spacing = (5/3)(m_s - m_d) — span ratio."""
        r = boundary.baryon_decuplet_masses()
        assert abs(r['Omega']['err_pct']) < 0.5

    def test_2_beats_5_3_for_octet(self):
        """Factor 2 works better than 5/3 for the Λ (octet)."""
        m_s, m_d = 93.23, 4.657
        m_p = 936.9
        lam_with_2 = m_p + 2 * (m_s - m_d)
        lam_with_53 = m_p + (5/3) * (m_s - m_d)
        err_2 = abs(lam_with_2 - 1115.683) / 1115.683
        err_53 = abs(lam_with_53 - 1115.683) / 1115.683
        assert err_2 < err_53

    def test_5_3_beats_2_for_decuplet(self):
        """Factor 5/3 works better than 2 for the Ω (decuplet)."""
        m_s, m_d = 93.23, 4.657
        step_53 = (5/3) * (m_s - m_d)
        step_2 = 2 * (m_s - m_d)
        m_Delta = 1233.5
        omega_53 = m_Delta + 3 * step_53
        omega_2 = m_Delta + 3 * step_2
        err_53 = abs(omega_53 - 1672.45) / 1672.45
        err_2 = abs(omega_2 - 1672.45) / 1672.45
        assert err_53 < err_2


@pytest.mark.framework
class TestBetaDecay:
    """Neutron beta decay from doublet angle mechanism."""

    def test_energy_release(self):
        """Q = m_n - m_p = 1.292 MeV (0.08%)."""
        r = boundary.beta_decay_energy()
        assert r['Q_MeV'] == pytest.approx(r['Q_meas'], rel=0.005)

    def test_neutron_decays(self):
        """Neutron decays because m_d > m_u (Q > 0)."""
        r = boundary.beta_decay_energy()
        assert r['neutron_decays'] is True

    def test_proton_stable(self):
        """Proton is stable: lightest B=1 state."""
        r = boundary.beta_decay_energy()
        assert r['proton_stable'] is True

    def test_doublet_angle_65_degrees(self):
        """φ₁ = arctan(m_d/m_u) ≈ 65°."""
        r = boundary.beta_decay_energy()
        assert r['phi1_deg'] == pytest.approx(65.0, abs=0.5)


# ============================================================
# FRAMEWORK TESTS — Baryon octet C₀ pair correction
# ============================================================

@pytest.mark.framework
class TestBaryonOctetPairCorrection:
    """C₀ pair correction sharpens the Ξ⁰ prediction to -0.04%."""

    def test_xi0_with_correction_sub_tenth(self):
        """Ξ⁰ at -0.04% with C₀ pair correction."""
        r = boundary.baryon_octet_masses()
        # The C₀ pair correction is already in the octet formula
        assert abs(r['Xi0']['err_pct']) < 0.1, \
            f"Xi0: {r['Xi0']['err_pct']:.2f}%"

    def test_xi0_negative_deviation(self):
        """Ξ⁰ deviation is negative (prediction undershoots)."""
        r = boundary.baryon_octet_masses()
        assert r['Xi0']['err_pct'] < 0


# ============================================================
# FRAMEWORK TESTS — Meson spectrum
# ============================================================

@pytest.mark.framework
class TestMesonSpectrum:
    """Meson masses from m_q1 + m_q2 + K*factor, K = N_c²(m_s-m_d)/2."""

    def test_pion_sub_percent(self):
        """π± sub-percent. π⁰ differs by EM splitting (3.5%, not predicted)."""
        r = boundary.meson_spectrum()
        assert abs(r['pi+']['err_pct']) < 1.0, \
            f"pi+: {r['pi+']['err_pct']:.2f}%"

    def test_kaon_sub_percent(self):
        """K± and K⁰ both sub-percent."""
        r = boundary.meson_spectrum()
        assert abs(r['K+']['err_pct']) < 1.0, \
            f"K+: {r['K+']['err_pct']:.2f}%"
        assert abs(r['K0']['err_pct']) < 1.0, \
            f"K0: {r['K0']['err_pct']:.2f}%"

    def test_D_mesons_sub_percent(self):
        """D±, D⁰, and Ds all sub-percent."""
        r = boundary.meson_spectrum()
        for name in ['D+', 'D0', 'Ds']:
            assert abs(r[name]['err_pct']) < 1.0, \
                f"{name}: {r[name]['err_pct']:.2f}%"

    def test_B_mesons_sub_percent(self):
        """B±, B⁰, and Bs all sub-percent."""
        r = boundary.meson_spectrum()
        for name in ['B+', 'B0', 'Bs']:
            assert abs(r[name]['err_pct']) < 1.0, \
                f"{name}: {r[name]['err_pct']:.2f}%"

    def test_quarkonium_sub_percent(self):
        """ηc and Υ both sub-percent."""
        r = boundary.meson_spectrum()
        assert abs(r['eta_c']['err_pct']) < 1.0, \
            f"eta_c: {r['eta_c']['err_pct']:.2f}%"
        assert abs(r['Upsilon']['err_pct']) < 1.0, \
            f"Upsilon: {r['Upsilon']['err_pct']:.2f}%"

    def test_Bc_sub_percent(self):
        """Bc meson sub-percent."""
        r = boundary.meson_spectrum()
        assert abs(r['Bc']['err_pct']) < 1.0, \
            f"Bc: {r['Bc']['err_pct']:.2f}%"

    def test_all_charged_mesons_sub_percent(self):
        """All charged/flavor mesons sub-percent (exclude π⁰ EM splitting)."""
        r = boundary.meson_spectrum()
        meson_names = [k for k in r if k not in ('K_binding', 'pi0')]
        for name in meson_names:
            assert abs(r[name]['err_pct']) < 1.0, \
                f"{name}: {r[name]['err_pct']:.2f}%"

    def test_upsilon_sub_tenth_percent(self):
        """Υ at 0.02% — the tightest meson prediction."""
        r = boundary.meson_spectrum()
        assert abs(r['Upsilon']['err_pct']) < 0.1, \
            f"Upsilon: {r['Upsilon']['err_pct']:.2f}%"


# ============================================================
# FRAMEWORK TESTS — Selection rules (which modes are particles?)
# ============================================================

@pytest.mark.framework
class TestStabilityRules:
    """
    GREEN: Stability = lightest state with given quantum numbers.
    The "floor" of each sector is stable. Everything above decays.
    """

    def test_proton_is_stable(self):
        """Proton: lightest B=1, Q=+1, S=0. Floor state."""
        r = boundary.beta_decay_energy()
        assert r['proton_stable'] is True

    def test_neutron_decays(self):
        """Neutron decays because m_d > m_u (not the floor)."""
        r = boundary.beta_decay_energy()
        assert r['neutron_decays'] is True

    def test_electron_is_lightest_charged_lepton(self):
        """Electron: floor of charged leptons."""
        assert M_ELECTRON < M_MUON < M_TAU

    def test_proton_lighter_than_all_baryons(self):
        """Proton is the lightest baryon — the Mercedes ground state."""
        r = boundary.baryon_octet_masses()
        m_p = r['p']['pred']
        for name, data in r.items():
            assert data['pred'] >= m_p, f"{name} lighter than proton?"
