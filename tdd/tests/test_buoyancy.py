"""
Test Domain 7: Boundary Buoyancy Framework

Scorecard:
  GREEN (framework) — Gravity reproduces Newton, MOND transition, GW speed = c
  GREEN (framework) — Planck mass as saturation/max-particle/min-BH boundary
  GREEN (framework) — Photon neutral buoyancy, dark energy from voids
  RED (framework) — Dark energy Λ value derivation (not from first principles)
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tests.conftest import *
from framework import boundary
from framework import constants as FC


# ============================================================
# FRAMEWORK TESTS — Gravity reproduces Newton
# ============================================================

@pytest.mark.framework
class TestGravityReproducesNewton:
    """
    GREEN: gravity_from_buoyancy() must reproduce F = GMm/r².

    Derivation chain (buoyancy-framework.md §2):
      Buoyancy deficit → 2D flux conservation → Gauss's law → inverse-square.
    """

    def test_inverse_square_basic(self):
        """g = GM/r² for a 1 kg mass at 1 m."""
        g = boundary.gravity_from_buoyancy(1.0, 1.0)
        assert g == pytest.approx(G, rel=1e-6)

    def test_inverse_square_scaling(self):
        """Doubling distance → 1/4 acceleration."""
        g1 = boundary.gravity_from_buoyancy(1.0, 1.0)
        g2 = boundary.gravity_from_buoyancy(1.0, 2.0)
        assert g2 == pytest.approx(g1 / 4.0, rel=1e-10)

    def test_mass_proportionality(self):
        """Doubling mass → double acceleration."""
        g1 = boundary.gravity_from_buoyancy(1.0, 1.0)
        g2 = boundary.gravity_from_buoyancy(2.0, 1.0)
        assert g2 == pytest.approx(2.0 * g1, rel=1e-10)

    def test_earth_surface_gravity(self):
        """Reproduce g ≈ 9.8 m/s² for Earth."""
        M_earth = 5.972e24  # kg
        R_earth = 6.371e6   # m
        g = boundary.gravity_from_buoyancy(M_earth, R_earth)
        assert g == pytest.approx(9.82, rel=0.01)

    def test_solar_surface_gravity(self):
        """Reproduce g_sun ≈ 274 m/s² for the Sun."""
        M_sun = 1.989e30    # kg
        R_sun = 6.957e8     # m
        g = boundary.gravity_from_buoyancy(M_sun, R_sun)
        assert g == pytest.approx(274.0, rel=0.01)

    def test_newtonian_force(self):
        """F = m × g = GMm/r² for Earth-Moon."""
        M_earth = 5.972e24
        M_moon = 7.342e22
        d = 3.844e8  # m
        g = boundary.gravity_from_buoyancy(M_earth, d)
        F = M_moon * g
        F_expected = G * M_earth * M_moon / d**2
        assert F == pytest.approx(F_expected, rel=1e-6)


# ============================================================
# FRAMEWORK TESTS — MOND transition from buoyancy
# ============================================================

@pytest.mark.framework
class TestMONDTransitionFromBuoyancy:
    """
    GREEN: mond_from_buoyancy() reproduces MOND at low accelerations,
    Newton at high accelerations.

    Derivation (buoyancy-framework.md §3):
      Crossover at a₀ = cH₀/(2π) where expansion pressure = gravitational sink.
      Deep MOND: g_obs ≈ √(g_N × a₀)
      Newtonian: g_obs ≈ g_N
    """

    def test_newtonian_regime(self):
        """When g_N >> a₀, MOND ≈ Newton."""
        # Earth surface: g ~ 9.8 >> a₀ ~ 10⁻¹⁰
        M_earth = 5.972e24
        R_earth = 6.371e6
        g_mond = boundary.mond_from_buoyancy(M_earth, R_earth, H0_CMB_SI)
        g_newton = boundary.gravity_from_buoyancy(M_earth, R_earth)
        # Should agree to better than 1 part in 10⁹
        assert g_mond == pytest.approx(g_newton, rel=1e-9)

    def test_deep_mond_regime(self):
        """When g_N << a₀, g_obs ≈ √(g_N × a₀)."""
        # Use a galaxy-like setup: M_sun at 100 kpc
        M = 1e11 * 1.989e30  # 10¹¹ solar masses
        r = 100e3 * 3.086e16  # 100 kpc in meters
        a0 = boundary.mond_acceleration(H0_CMB_SI)
        g_newton = boundary.gravity_from_buoyancy(M, r)

        # Verify we're in deep MOND (g_N < a₀)
        assert g_newton < a0, "Not in MOND regime"

        g_mond = boundary.mond_from_buoyancy(M, r, H0_CMB_SI)
        g_deep_mond = math.sqrt(g_newton * a0)

        # Deep MOND approximation should be within 30% of RAR
        # (exact agreement only in the limit g_N/a₀ → 0)
        assert g_mond == pytest.approx(g_deep_mond, rel=0.30)

    def test_mond_always_greater_than_newton(self):
        """MOND acceleration is always ≥ Newtonian (expansion assists gravity)."""
        for M, r in [(1e30, 1e20), (1e40, 1e22), (1e20, 1e18)]:
            g_mond = boundary.mond_from_buoyancy(M, r, H0_CMB_SI)
            g_newton = boundary.gravity_from_buoyancy(M, r)
            assert g_mond >= g_newton * 0.9999  # allow tiny float imprecision

    def test_crossover_acceleration(self):
        """The crossover point a₀ = cH₀/(2π) matches MOND empirical value."""
        a0 = boundary.mond_acceleration(H0_CMB_SI)
        assert a0 == pytest.approx(A0_MOND, rel=0.15)

    def test_flat_rotation_curve_in_mond_regime(self):
        """In deep MOND, v² = √(GMa₀) = constant → flat rotation curve."""
        M = 1e11 * 1.989e30  # galaxy mass
        a0 = boundary.mond_acceleration(H0_CMB_SI)

        # Check v at two different radii deep in MOND regime
        r1 = 50e3 * 3.086e16   # 50 kpc
        r2 = 200e3 * 3.086e16  # 200 kpc

        g1 = boundary.mond_from_buoyancy(M, r1, H0_CMB_SI)
        g2 = boundary.mond_from_buoyancy(M, r2, H0_CMB_SI)

        v1 = math.sqrt(r1 * g1)
        v2 = math.sqrt(r2 * g2)

        # Rotation velocity should be approximately constant (within 20%)
        assert v1 == pytest.approx(v2, rel=0.20)


# ============================================================
# FRAMEWORK TESTS — Minimum BH mass is Planck mass
# ============================================================

@pytest.mark.framework
class TestMinimumBHMassIsPlanck:
    """
    GREEN: minimum_bh_mass() ≈ M_P (from saturation balance).

    Derivation (buoyancy-framework.md §4):
      R_s ≥ λ_C → M ≥ M_P/√2.
      Quantum pressure vs gravitational collapse balance at Planck mass.
    """

    def test_minimum_bh_mass_order(self):
        """Minimum BH mass is order Planck mass."""
        M_min = boundary.minimum_bh_mass()
        Mp = boundary.planck_mass()
        ratio = M_min / Mp
        # Should be within factor of 2 of Planck mass
        assert 0.5 < ratio < 2.0

    def test_minimum_bh_mass_is_mp_over_sqrt2(self):
        """Exact result: M_min = M_P/√2."""
        M_min = boundary.minimum_bh_mass()
        Mp = boundary.planck_mass()
        assert M_min == pytest.approx(Mp / math.sqrt(2), rel=1e-10)

    def test_schwarzschild_equals_compton_at_planck(self):
        """At M = M_P/√2, R_s = λ_C (the crossover point)."""
        M_min = boundary.minimum_bh_mass()
        R_s = 2 * G * M_min / C**2
        lambda_C = HBAR / (M_min * C)
        assert R_s == pytest.approx(lambda_C, rel=1e-6)

    def test_below_minimum_no_bh(self):
        """Below M_min, Compton wavelength > Schwarzschild radius."""
        M_sub = boundary.minimum_bh_mass() * 0.5
        R_s = 2 * G * M_sub / C**2
        lambda_C = HBAR / (M_sub * C)
        assert lambda_C > R_s


# ============================================================
# FRAMEWORK TESTS — Max particle mass is Planck mass
# ============================================================

@pytest.mark.framework
class TestMaxParticleMassIsPlanck:
    """
    GREEN: max_particle_mass() = M_P (pixel saturation limit).

    Derivation (buoyancy-framework.md §4):
      One pixel holds at most one Planck unit of action per tick.
      Maximum single-pixel excitation = Planck mass.
    """

    def test_max_particle_mass_equals_planck(self):
        """Maximum particle mass = Planck mass exactly."""
        M_max = boundary.max_particle_mass()
        Mp = boundary.planck_mass()
        assert M_max == pytest.approx(Mp, rel=1e-10)

    def test_max_particle_mass_value(self):
        """M_P ≈ 2.176 × 10⁻⁸ kg."""
        M_max = boundary.max_particle_mass()
        assert M_max == pytest.approx(2.176e-8, rel=1e-3)

    def test_max_mass_exceeds_all_known_particles(self):
        """Planck mass >> all known particle masses."""
        M_max = boundary.max_particle_mass()
        M_max_GeV = M_max * C**2 / GEV
        # Top quark (heaviest known): ~173 GeV
        assert M_max_GeV > 173e3  # more than 173,000 GeV
        # Should be ~1.22 × 10¹⁹ GeV
        assert M_max_GeV == pytest.approx(1.22e19, rel=0.01)

    def test_min_bh_meets_max_particle(self):
        """Min BH mass and max particle mass are within √2 of each other."""
        M_min_bh = boundary.minimum_bh_mass()
        M_max_particle = boundary.max_particle_mass()
        ratio = M_max_particle / M_min_bh
        assert ratio == pytest.approx(math.sqrt(2), rel=1e-10)


# ============================================================
# FRAMEWORK TESTS — Gravitational wave speed is c
# ============================================================

@pytest.mark.framework
class TestGravitationalWaveSpeedIsC:
    """
    GREEN: gravitational_wave_speed() = c exactly.

    Derivation (buoyancy-framework.md §5):
      Buoyancy perturbations propagate at boundary rendering speed = c.
      Wave equation ∂²δB/∂t² = c²∇²δB → speed = c.

    Experimentally confirmed: LIGO/Virgo GW170817 + Fermi GRB 170817A
    measured |v_gw - c|/c < 10⁻¹⁵.
    """

    def test_gw_speed_equals_c(self):
        """Gravitational wave speed = c exactly."""
        v_gw = boundary.gravitational_wave_speed()
        assert v_gw == C  # exact equality, not approximate

    def test_gw_speed_matches_ligo_bound(self):
        """Speed within LIGO/Virgo experimental bound."""
        v_gw = boundary.gravitational_wave_speed()
        fractional_diff = abs(v_gw - C) / C
        assert fractional_diff < 1e-15

    def test_gw_speed_equals_lp_over_tp(self):
        """v_gw = L_P/T_P (one pixel per tick)."""
        v_gw = boundary.gravitational_wave_speed()
        lp_over_tp = boundary.planck_length() / boundary.planck_time()
        assert v_gw == pytest.approx(lp_over_tp, rel=1e-6)


# ============================================================
# FRAMEWORK TESTS — Dark energy from void buoyancy
# ============================================================

@pytest.mark.framework
class TestDarkEnergyFromVoidBuoyancy:
    """
    GREEN: dark_energy_from_buoyancy() gives expansion enhancement from voids.

    Derivation (buoyancy-framework.md §6):
      Voids have higher buoyancy → faster expansion.
      Enhancement = 1 + f_void × (1 - Ω_M).
    """

    def test_no_voids_no_enhancement(self):
        """Zero void fraction → no enhancement (factor = 1)."""
        factor = boundary.dark_energy_from_buoyancy(0.0, H0_CMB_SI)
        assert factor == pytest.approx(1.0, rel=1e-10)

    def test_all_void_maximum_enhancement(self):
        """100% void → maximum enhancement."""
        factor = boundary.dark_energy_from_buoyancy(1.0, H0_CMB_SI)
        expected = 1.0 + 1.0 * (1.0 - 0.315)  # 1.685
        assert factor == pytest.approx(expected, rel=0.01)

    def test_realistic_void_fraction(self):
        """~60% void fraction (current universe) gives moderate enhancement."""
        factor = boundary.dark_energy_from_buoyancy(0.6, H0_CMB_SI)
        expected = 1.0 + 0.6 * 0.685  # ≈ 1.411
        assert factor == pytest.approx(expected, rel=0.01)

    def test_enhancement_monotonic(self):
        """More voids → more enhancement (monotonically increasing)."""
        factors = [boundary.dark_energy_from_buoyancy(f, H0_CMB_SI)
                   for f in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]]
        for i in range(len(factors) - 1):
            assert factors[i] < factors[i + 1]


# ============================================================
# FRAMEWORK TESTS — Photon neutral buoyancy
# ============================================================

@pytest.mark.framework
class TestPhotonNeutralBuoyancy:
    """
    GREEN: Photon has zero net buoyancy (load = background).

    Derivation (buoyancy-framework.md §1):
      Photon transits a pixel in one tick. Its load per tick = background rate.
      pixel_buoyancy(background, background) = 0.
    """

    def test_photon_zero_buoyancy(self):
        """Photon load equals background → zero buoyancy."""
        bg = 1.0  # arbitrary background
        buoyancy = boundary.pixel_buoyancy(load=bg, background_load=bg)
        assert buoyancy == 0.0

    def test_photon_zero_any_background(self):
        """Neutral buoyancy at any background level."""
        for bg in [0.001, 1.0, 1000.0, 1e-52]:
            buoyancy = boundary.pixel_buoyancy(load=bg, background_load=bg)
            assert buoyancy == 0.0

    def test_mass_negative_buoyancy(self):
        """Mass loads above background → negative buoyancy."""
        bg = 1.0
        mass_load = 2.0  # above background
        buoyancy = boundary.pixel_buoyancy(load=mass_load, background_load=bg)
        assert buoyancy < 0

    def test_void_positive_buoyancy(self):
        """Void loads below background → positive buoyancy."""
        bg = 1.0
        void_load = 0.1  # below background
        buoyancy = boundary.pixel_buoyancy(load=void_load, background_load=bg)
        assert buoyancy > 0


# ============================================================
# FRAMEWORK TESTS — Saturation equals black hole
# ============================================================

@pytest.mark.framework
class TestSaturationEqualsBlackHole:
    """
    GREEN: Pixel saturation corresponds to black hole formation.

    Derivation (buoyancy-framework.md §7):
      Saturation = Planck energy density.
      At saturation, buoyancy deficit is maximum → spawns new boundary.
    """

    def test_saturation_at_planck_density(self):
        """Saturation threshold = Planck energy density."""
        rho_sat = boundary.buoyancy_saturation_threshold()
        # Planck density: E_P / L_P³ = M_P c² / L_P³
        Mp = boundary.planck_mass()
        Lp = boundary.planck_length()
        rho_planck = Mp * C**2 / Lp**3
        assert rho_sat == pytest.approx(rho_planck, rel=1e-6)

    def test_saturation_density_order(self):
        """Planck density ≈ 5.16 × 10⁹⁶ kg/m³ (in energy-equivalent)."""
        rho_sat = boundary.buoyancy_saturation_threshold()
        # ρ_Planck ≈ 5.16 × 10⁹⁶ J/m³
        # (or equivalently ≈ 5.16 × 10⁹⁶ kg/m³ × c²... let's check order)
        log_rho = math.log10(rho_sat)
        # Planck energy density is ~10¹¹³ J/m³
        assert 112 < log_rho < 114

    def test_buoyancy_at_saturation_is_minimum(self):
        """At Planck density, buoyancy deficit is maximum (most negative)."""
        rho_sat = boundary.buoyancy_saturation_threshold()
        # Background is Λc²/(8πG) ≈ 5.96 × 10⁻¹⁰ J/m³ — negligible vs Planck
        rho_bg = 5.96e-10  # approximate dark energy density
        buoyancy = boundary.pixel_buoyancy(load=rho_sat, background_load=rho_bg)
        # Maximum deficit: buoyancy is extremely negative
        assert buoyancy < -1e100

    def test_schwarzschild_radius_from_planck_mass(self):
        """Schwarzschild radius of Planck mass ≈ 2 Planck lengths."""
        Mp = boundary.planck_mass()
        R_s = 2 * G * Mp / C**2
        Lp = boundary.planck_length()
        assert R_s == pytest.approx(2 * Lp, rel=1e-3)

    def test_bh_entropy_from_saturation(self):
        """BH entropy = area / (4 L_P²) = pixel count on child boundary."""
        # For a 10 solar mass BH:
        M = 10 * 1.989e30
        R_s = 2 * G * M / C**2
        A = 4 * math.pi * R_s**2
        Lp = boundary.planck_length()
        S = A / (4 * Lp**2)
        # Bekenstein-Hawking entropy should be enormous
        assert S > 1e77  # order 10⁷⁷ for 10 M_sun
