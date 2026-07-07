"""
Test Domain: Figure-8 Spin Model — Topology, Trough/Peak Asymmetry, and Mass

Tests the figure-8 spin model formalized in figure-8-spin-model.md:
  - Spin from mode topology (figure-8, circle, dent)
  - Mass from trough/peak dwell asymmetry
  - Pauli exclusion from phase-slot counting
  - Self-trough bounded by Planck mass
  - Pair production topology conservation

Scorecard:
  GREEN  — spin classification, Pauli counting, topology conservation, Planck bounds
  RED    — anything requiring the selection rules we don't have yet
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tests.conftest import *
from framework import boundary
from framework.loeschian import loeschian


# ============================================================
# SPIN FROM TOPOLOGY
# ============================================================

class TestSpinTopology:
    """Verify spin classification from mode self-intersection count."""

    @pytest.mark.framework
    def test_electron_is_figure8(self):
        """
        Electron is a figure-8 mode (1 self-intersection) → spin 1/2.

        Derivation: A figure-8 path on the 2D boundary has one self-intersection.
        Two lobe traversals are needed per cycle → mode frequency = ω_P/2
        → angular momentum = ℏ/2 → spin 1/2.
        """
        spin = boundary.spin_from_topology(self_intersections=1)
        assert spin == 0.5, f"Figure-8 should give spin 1/2, got {spin}"

    @pytest.mark.framework
    def test_photon_is_circle(self):
        """
        Photon is a simple closed curve (0 self-intersections) → spin 1.

        Derivation: A circle on the boundary has no self-intersection.
        One traversal = one full boundary cycle → ω_mode = ω_P → spin = ℏ.
        """
        spin = boundary.spin_from_topology(self_intersections=0)
        assert spin == 1.0, f"Circle should give spin 1, got {spin}"

    @pytest.mark.framework
    def test_higgs_is_dent(self):
        """
        Higgs is a pure surface deformation (dent, no circulation) → spin 0.

        Derivation: A dent has no winding, no circulation, no angular momentum.
        It's a static boundary deformation — the vacuum expectation value
        of the Higgs field is a permanent local geometry change.
        """
        spin = boundary.spin_from_topology(self_intersections=-1)
        assert spin == 0.0, f"Dent should give spin 0, got {spin}"

    @pytest.mark.framework
    def test_all_fermions_are_figure8(self):
        """All known elementary fermions should have spin 1/2 from figure-8 topology."""
        fermion_spin = boundary.spin_from_topology(self_intersections=1)
        # Electron, muon, tau, neutrinos, quarks — all spin 1/2
        assert fermion_spin == 0.5

    @pytest.mark.framework
    def test_spin_quantized(self):
        """Spin values are discrete: only 0, 1/2, 1 for elementary particles."""
        assert boundary.spin_from_topology(-1) == 0.0   # dent
        assert boundary.spin_from_topology(0) == 1.0    # circle
        assert boundary.spin_from_topology(1) == 0.5    # figure-8

    @pytest.mark.framework
    def test_higher_spin_unstable(self):
        """
        Higher self-intersection counts give smaller spin values.
        These modes are expected to be topologically unstable on a 2D surface
        (no elementary spin-3/2+ particles observed).
        """
        # 2 self-intersections → spin 1/3 (hypothetical, not observed as elementary)
        spin_2 = boundary.spin_from_topology(self_intersections=2)
        assert spin_2 == pytest.approx(1.0/3.0)
        # Confirms the formula: spin = 1/(n+1) for n self-intersections

    @pytest.mark.framework
    def test_invalid_self_intersections(self):
        """Negative values below -1 are invalid."""
        with pytest.raises(ValueError):
            boundary.spin_from_topology(self_intersections=-2)


# ============================================================
# PAULI EXCLUSION FROM PHASE SLOTS
# ============================================================

class TestPauliExclusion:
    """Verify Pauli exclusion from phase-slot counting on oscillating boundary."""

    @pytest.mark.framework
    def test_pauli_from_phase_slots(self):
        """
        Exactly 2 fermions per quantum state.

        Derivation: The boundary oscillation has two half-cycles per period
        (trough and peak). A figure-8 mode (spin 1/2) can start on either:
          - Trough-first (spin up)
          - Peak-first (spin down)
        No third phase slot exists → maximum 2 fermions per state.

        N = 2s + 1 = 2(1/2) + 1 = 2
        """
        capacity = boundary.pauli_exclusion_capacity(spin=0.5)
        assert capacity == 2, f"Spin 1/2 should allow 2 states, got {capacity}"

    @pytest.mark.framework
    def test_boson_states(self):
        """Spin 1 bosons have 3 polarization states (2s+1 = 3)."""
        capacity = boundary.pauli_exclusion_capacity(spin=1.0)
        assert capacity == 3

    @pytest.mark.framework
    def test_scalar_single_state(self):
        """Spin 0 (Higgs-like dent) has 1 state (2×0+1 = 1)."""
        capacity = boundary.pauli_exclusion_capacity(spin=0.0)
        assert capacity == 1

    @pytest.mark.framework
    def test_pauli_rejects_negative_spin(self):
        """Negative spin is unphysical."""
        with pytest.raises(ValueError):
            boundary.pauli_exclusion_capacity(spin=-0.5)

    @pytest.mark.framework
    def test_pauli_rejects_non_half_integer(self):
        """Only integer or half-integer spins are valid."""
        with pytest.raises(ValueError):
            boundary.pauli_exclusion_capacity(spin=0.3)


# ============================================================
# FERMION MASS FROM TROUGH/PEAK ASYMMETRY
# ============================================================

class TestFigure8DwellTime:
    """Test mass from figure-8 dwell time on asymmetric boundary."""

    @pytest.mark.framework
    def test_fermion_mass_requires_two_passes(self):
        """
        Fermion mass comes from trough/peak asymmetry over TWO passes.

        A figure-8 mode traverses trough then peak (or vice versa).
        The dwell time difference between the two passes creates mass.
        With zero trough depth (flat boundary), mass = 0.

        m/M_P = L × δ × (1 + δ)
        At δ = 0: m = 0 (no asymmetry → no mass, like a massless mode)
        """
        # Zero trough depth → zero mass (no asymmetry to create dwell difference)
        mass = boundary.figure8_dwell_time(eigenvalue=1, trough_depth=0.0)
        assert mass == 0.0, "Zero trough depth should give zero mass"

    @pytest.mark.framework
    def test_boson_can_be_massless(self):
        """
        Bosons (circles on flat boundary) have zero trough depth → massless.

        Photons and gluons traverse flat boundary regions with no trough.
        Their single-pass traversal encounters no asymmetry → m = 0.
        """
        # A boson on flat boundary: any eigenvalue, but δ = 0
        mass = boundary.figure8_dwell_time(eigenvalue=100, trough_depth=0.0)
        assert mass == 0.0, "Boson on flat boundary should be massless"

    @pytest.mark.framework
    def test_mass_increases_with_eigenvalue(self):
        """Higher Loeschian eigenvalue → more mass (larger lobes, deeper trough sampling)."""
        delta = 1e-23  # typical sub-Planckian trough depth
        m1 = boundary.figure8_dwell_time(eigenvalue=1, trough_depth=delta)
        m2 = boundary.figure8_dwell_time(eigenvalue=9, trough_depth=delta)
        m3 = boundary.figure8_dwell_time(eigenvalue=3477, trough_depth=delta)
        assert m1 < m2 < m3, "Mass should increase with eigenvalue"

    @pytest.mark.framework
    def test_mass_increases_with_trough_depth(self):
        """Deeper trough → more dwell asymmetry → more mass."""
        L = loeschian(1, 0)  # electron eigenvalue = 1
        m1 = boundary.figure8_dwell_time(eigenvalue=L, trough_depth=0.01)
        m2 = boundary.figure8_dwell_time(eigenvalue=L, trough_depth=0.1)
        m3 = boundary.figure8_dwell_time(eigenvalue=L, trough_depth=0.5)
        assert m1 < m2 < m3, "Mass should increase with trough depth"

    @pytest.mark.framework
    def test_self_reinforcement_factor(self):
        """
        The (1+δ) self-reinforcement factor: mass loads boundary → deepens trough.

        At small δ: (1+δ) ≈ 1, negligible self-reinforcement.
        At δ = 0.5: (1+δ) = 1.5, 50% enhancement from self-trough.
        """
        L = 1
        d_small = 1e-20
        d_large = 0.5
        m_small = boundary.figure8_dwell_time(L, d_small)
        m_large = boundary.figure8_dwell_time(L, d_large)

        # For small δ: m ≈ L × δ (linear)
        assert m_small == pytest.approx(L * d_small, rel=1e-10)

        # For δ = 0.5: m = L × 0.5 × 1.5 = 0.75 (nonlinear enhancement)
        assert m_large == pytest.approx(L * 0.5 * 1.5)

    @pytest.mark.framework
    def test_self_trough_bounded_by_planck_mass(self):
        """
        Maximum particle mass is bounded by Planck mass (pixel saturation).

        At δ = 1 (maximum trough depth, Planck density):
          m/M_P = L × 1 × (1+1) = 2L

        For the electron (L=1), max m = 2 M_P.
        But δ > 1 is forbidden (pixel can't hold more than Planck energy).
        The trough_depth parameter is clamped to [0, 1].
        """
        # Maximum trough depth = 1.0 (Planck saturation)
        m_max = boundary.figure8_dwell_time(eigenvalue=1, trough_depth=1.0)
        assert m_max == 2.0  # 1 × 1 × (1+1) = 2 Planck masses

        # Beyond 1.0 is unphysical — should raise
        with pytest.raises(ValueError):
            boundary.figure8_dwell_time(eigenvalue=1, trough_depth=1.1)

    @pytest.mark.framework
    def test_tau_electron_mass_ratio_from_eigenvalues(self):
        """
        Tau/electron mass ratio from Loeschian eigenvalues.

        At fixed trough depth δ (same boundary, same asymmetry):
          m_tau/m_e = L(37,31)/L(1,0) = 3477/1 = 3477

        The measured ratio is 3477.23 ± 0.23, consistent within 1σ.
        """
        L_tau = loeschian(37, 31)
        L_electron = loeschian(1, 0)
        assert L_tau / L_electron == 3477


# ============================================================
# TROUGH ASYMMETRY (INVERSE: MASS → TROUGH DEPTH)
# ============================================================

class TestTroughAsymmetry:
    """Test inversion: observed mass → self-reinforcing trough depth."""

    @pytest.mark.framework
    def test_electron_trough_very_shallow(self):
        """
        Electron (0.511 MeV) creates an extremely shallow self-trough.
        δ ≈ m_e/M_P ≈ 4×10⁻²³ — negligible self-reinforcement.
        """
        delta = boundary.trough_asymmetry(M_ELECTRON)
        assert delta > 0, "Massive particle must have positive trough depth"
        assert delta < 1e-20, f"Electron trough should be ≪ 1, got {delta}"
        # Should be approximately m_e / M_P (linear regime)
        m_e_kg = M_ELECTRON * 1e6 * EV / C**2
        expected = m_e_kg / M_PLANCK
        assert delta == pytest.approx(expected, rel=1e-6)

    @pytest.mark.framework
    def test_top_quark_trough_still_small(self):
        """
        Even the top quark (173 GeV) has δ ≈ 1.4×10⁻¹⁷ — still deeply sub-Planckian.
        All known particles are in the linear regime.
        """
        delta = boundary.trough_asymmetry(M_TOP)
        assert delta < 1e-15, f"Top quark trough should be ≪ 1, got {delta}"
        assert delta > 0

    @pytest.mark.framework
    def test_trough_monotonic_with_mass(self):
        """Heavier particles → deeper self-troughs."""
        d_e = boundary.trough_asymmetry(M_ELECTRON)
        d_mu = boundary.trough_asymmetry(M_MUON)
        d_tau = boundary.trough_asymmetry(M_TAU)
        d_top = boundary.trough_asymmetry(M_TOP)
        assert d_e < d_mu < d_tau < d_top

    @pytest.mark.framework
    def test_trough_approaches_planck(self):
        """
        At Planck mass, trough depth → (√5 - 1)/2 ≈ 0.618.
        This is where self-reinforcement becomes significant.
        """
        # Planck mass in MeV
        M_P_MeV = E_PLANCK * 1e3  # E_PLANCK is in GeV
        delta = boundary.trough_asymmetry(M_P_MeV)
        golden = (math.sqrt(5) - 1) / 2  # ≈ 0.618
        assert delta == pytest.approx(golden, rel=1e-6)

    @pytest.mark.framework
    def test_massless_has_zero_trough(self):
        """Massless particles (photon) have zero trough depth."""
        delta = boundary.trough_asymmetry(0.0)
        assert delta == 0.0

    @pytest.mark.framework
    def test_roundtrip_mass_trough_mass(self):
        """
        Mass → trough depth → mass should roundtrip consistently.
        figure8_dwell_time(L=1, δ) × M_P should recover the original mass.
        """
        for mass_mev in [M_ELECTRON, M_MUON, M_TAU, M_TOP]:
            delta = boundary.trough_asymmetry(mass_mev)
            # figure8_dwell_time returns m/M_P for L=1
            m_over_Mp = boundary.figure8_dwell_time(eigenvalue=1, trough_depth=delta)
            # Convert back to MeV
            m_recovered_kg = m_over_Mp * M_PLANCK
            m_recovered_mev = m_recovered_kg * C**2 / (1e6 * EV)
            assert m_recovered_mev == pytest.approx(mass_mev, rel=1e-6), \
                f"Roundtrip failed for {mass_mev} MeV"


# ============================================================
# PAIR PRODUCTION TOPOLOGY CONSERVATION
# ============================================================

class TestPairProduction:
    """Test pair production as topological splitting."""

    @pytest.mark.framework
    def test_pair_production_conserves_topology(self):
        """
        Pair production from vacuum: net winding 0 → particle + antiparticle.
        Net winding after = 0 (conserved).
        """
        result = boundary.pair_production_topology(winding_before=0)
        assert result['topology_conserved'] is True
        assert result['winding_before'] == 0
        assert result['winding_after'] == 0
        assert result['net_charge'] == 0
        assert result['net_spin'] == 0

    @pytest.mark.framework
    def test_pair_topology_particle_properties(self):
        """Particle is CW figure-8 (trough-first), antiparticle is CCW (peak-first)."""
        result = boundary.pair_production_topology()
        p = result['particle']
        ap = result['antiparticle']

        assert p['topology'] == 'figure-8'
        assert ap['topology'] == 'figure-8'
        assert p['winding'] == +1
        assert ap['winding'] == -1
        assert p['spin'] == 0.5
        assert ap['spin'] == 0.5
        assert p['phase_start'] != ap['phase_start']  # opposite phase

    @pytest.mark.framework
    def test_pair_production_from_nonzero_winding(self):
        """Pair production conserves topology even from nonzero initial winding."""
        for w in [-2, -1, 0, 1, 2]:
            result = boundary.pair_production_topology(winding_before=w)
            assert result['topology_conserved'] is True
            assert result['winding_after'] == w


# ============================================================
# INTEGRATION: Spin + Mass + Pauli together
# ============================================================

class TestIntegration:
    """Cross-cutting tests combining spin, mass, and exclusion."""

    @pytest.mark.framework
    def test_electron_complete_picture(self):
        """
        Electron: figure-8, spin 1/2, 2 states per level,
        L(1,0) = 1 (lightest fermion), shallow self-trough.
        """
        spin = boundary.spin_from_topology(1)
        assert spin == 0.5

        states = boundary.pauli_exclusion_capacity(spin)
        assert states == 2

        L = loeschian(1, 0)
        assert L == 1  # simplest Loeschian = lightest fermion

        delta = boundary.trough_asymmetry(M_ELECTRON)
        assert delta > 0  # massive
        assert delta < 1e-20  # but very shallow

    @pytest.mark.framework
    def test_photon_complete_picture(self):
        """
        Photon: circle, spin 1, 3 polarization states (but only 2 for massless),
        no trough, massless.
        """
        spin = boundary.spin_from_topology(0)
        assert spin == 1.0

        states = boundary.pauli_exclusion_capacity(spin)
        assert states == 3  # 2s+1 = 3 (longitudinal requires mass)

        # Massless: no trough
        mass = boundary.figure8_dwell_time(eigenvalue=1, trough_depth=0.0)
        assert mass == 0.0

    @pytest.mark.framework
    def test_higgs_complete_picture(self):
        """
        Higgs: dent, spin 0, 1 state, pure surface deformation.
        """
        spin = boundary.spin_from_topology(-1)
        assert spin == 0.0

        states = boundary.pauli_exclusion_capacity(spin)
        assert states == 1

    @pytest.mark.framework
    def test_mass_hierarchy_from_eigenvalues(self):
        """
        Mass hierarchy: electron < muon < tau follows from Loeschian eigenvalues
        at fixed trough asymmetry.
        """
        delta = 1e-23  # any fixed small trough depth
        L_e = loeschian(1, 0)       # 1
        L_mu = loeschian(12, 4)     # 208 (nearest Loeschian to muon ratio)
        L_tau = loeschian(37, 31)   # 3477

        m_e = boundary.figure8_dwell_time(L_e, delta)
        m_mu = boundary.figure8_dwell_time(L_mu, delta)
        m_tau = boundary.figure8_dwell_time(L_tau, delta)

        assert m_e < m_mu < m_tau
        # Ratios should match Loeschian ratios exactly (same δ cancels)
        assert m_mu / m_e == pytest.approx(208.0 / 1.0)
        assert m_tau / m_e == pytest.approx(3477.0 / 1.0)
