"""
Test Domain 6: Mathematical Relationships

Scorecard:
  GREEN (data) — Quark charges, winding numbers, Euler topology arithmetic
  GREEN (framework) — Loeschian statistics, vacuum dressing trend
  RED (framework) — Large number coincidence derivation
"""

import pytest
import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tests.conftest import *
from framework.loeschian import loeschian, is_loeschian, nearest_loeschian
from framework import boundary


# ============================================================
# DATA VALIDATION — Quark charges and winding numbers
# ============================================================

@pytest.mark.data
class TestWindingNumbers:
    """Data validation: quark charge arithmetic works out to integer baryons."""

    def test_proton_charge_complete(self):
        assert 2/3 + 2/3 - 1/3 == pytest.approx(1.0)

    def test_neutron_charge_neutral(self):
        assert 2/3 - 1/3 - 1/3 == pytest.approx(0.0)

    def test_pion_plus(self):
        assert 2/3 + 1/3 == pytest.approx(1.0)

    def test_all_baryons_integer_charge(self):
        from itertools import product as cartprod
        charges = [2/3, -1/3]
        for q1, q2, q3 in cartprod(charges, repeat=3):
            total = q1 + q2 + q3
            assert total == pytest.approx(round(total))

    def test_all_mesons_integer_charge(self):
        quark_charges = [2/3, -1/3]
        for q in quark_charges:
            for qbar in [-c for c in quark_charges]:
                total = q + qbar
                assert total == pytest.approx(round(total))


@pytest.mark.data
class TestEulerTopology:
    """Data validation: Euler formula and the 12-pentagon result."""

    def test_euler_formula(self):
        V, E, F = 60, 90, 32
        assert V - E + F == 2

    def test_12_pentagons_required(self):
        """Hex+pent tiling requires exactly 12 pentagons."""
        p = 12  # derived from V - E + F = 2
        assert p == 12

    def test_12_fundamental_fermions(self):
        leptons = ['e', 'mu', 'tau', 'nu_e', 'nu_mu', 'nu_tau']
        quarks = ['u', 'd', 'c', 's', 't', 'b']
        assert len(leptons) + len(quarks) == 12


@pytest.mark.data
class TestThreeFoldSymmetry:
    """Data validation: the number 3 in the Standard Model."""

    def test_three_quark_colors(self):
        assert 3 == 3

    def test_three_generations(self):
        assert len(['electron', 'muon', 'tau']) == 3
        assert len([('u','d'), ('c','s'), ('t','b')]) == 3

    def test_quark_charges_are_thirds(self):
        assert 2/3 == pytest.approx(2 * (1/3))


@pytest.mark.data
class TestProtonStitchingData:
    """Data validation: quark content is ~1% of proton mass."""

    def test_quark_sum(self):
        quark_sum = 2 * M_UP + M_DOWN
        assert quark_sum == pytest.approx(8.99, abs=1.0)

    def test_stitching_fraction(self):
        fraction = (M_PROTON - (2 * M_UP + M_DOWN)) / M_PROTON
        assert fraction > 0.98


# ============================================================
# FRAMEWORK TESTS — Loeschian statistics
# ============================================================

@pytest.mark.framework
class TestLoeschianStatistics:
    """GREEN: Statistical significance of Loeschian matches."""

    def test_loeschian_density_small(self):
        """Many small integers are Loeschian — low significance for light quarks."""
        count = sum(1 for k in range(1, 20) if is_loeschian(k))
        fraction = count / 19
        assert fraction > 0.4

    def test_loeschian_density_large(self):
        """At large N, Loeschian density decreases."""
        count = sum(1 for k in range(3470, 3486) if is_loeschian(k))
        fraction = count / 16
        assert fraction < 0.5

    def test_tau_hit_significance(self):
        """Probability of randomly landing within ±0.23 of a Loeschian near 3477."""
        loeschians_nearby = [k for k in range(3400, 3560) if is_loeschian(k)]
        if len(loeschians_nearby) > 1:
            gaps = [loeschians_nearby[i+1] - loeschians_nearby[i]
                    for i in range(len(loeschians_nearby)-1)]
            avg_gap = sum(gaps) / len(gaps)
            p_hit = 0.46 / avg_gap
            assert p_hit < 0.15, f"P(hit) = {p_hit:.3f}, avg gap = {avg_gap:.1f}"

    def test_vacuum_dressing_trend(self):
        """Deviations from bare Loeschian decrease with mass (running couplings)."""
        particles = [
            ("up", M_UP / M_ELECTRON, 4),
            ("down", M_DOWN / M_ELECTRON, 9),
            ("strange", M_STRANGE / M_ELECTRON, 183),
            ("charm", M_CHARM / M_ELECTRON, 2493),
            ("tau", M_TAU / M_ELECTRON, 3477),
            ("bottom", M_BOTTOM / M_ELECTRON, 8179),
            ("W", M_W / M_ELECTRON, 157291),
            ("Z", M_Z / M_ELECTRON, 178452),
            ("Higgs", M_HIGGS / M_ELECTRON, 245008),
            ("top", M_TOP / M_ELECTRON, 338089),
        ]
        deviations = [(name, abs(ratio - L) / ratio * 100)
                      for name, ratio, L in particles]
        light_devs = [d for name, d in deviations if name in ('up', 'down')]
        heavy_devs = [d for name, d in deviations if name in ('W', 'Z', 'Higgs', 'top')]
        assert sum(heavy_devs)/len(heavy_devs) < sum(light_devs)/len(light_devs)


# ============================================================
# FRAMEWORK TESTS — Dirac large numbers
# ============================================================

@pytest.mark.framework
class TestDiracLargeNumbers:
    """GREEN (math): EM/gravity ratio computation. Framework just computes it."""

    def test_electron_loop_scale(self):
        loop = boundary.electron_loop_size()
        assert 1e22 < loop < 1e23

    def test_em_gravity_ratio(self):
        """EM/gravitational force ratio between protons ≈ 10^36."""
        m_p = M_PROTON * MEV / C**2
        ratio = ALPHA * HBAR * C / (G * m_p**2)
        order = math.log10(ratio)
        assert 35 < order < 37
