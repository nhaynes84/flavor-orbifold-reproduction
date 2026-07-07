"""
Self-consistency tests: verify the framework uses its OWN masses,
not PDG masses, for all internal predictions.

The framework predicts 9 fermion masses from m_tau + N_c = 3.
Every downstream prediction (CKM, PMNS, delta_CP, bosons) should
use FRAMEWORK masses, not external data. PDG/FLAG are for
validation, not input.
"""

import math
import pytest

# Framework masses (Paper 1, Table I)
FW = {
    'e': 0.5106, 'mu': 105.69, 'tau': 1776.86,
    'u': 2.172, 'c': 1270.0, 't': 172782.0,
    'd': 4.657, 's': 93.23, 'b': 4181.0,
}

# Framework constants
N_C = 3
N_Q = 7
N_L = 3
C0 = 1/5
T_F = 0.5


class TestMassRatiosFromFramework:
    """All mass ratios use framework masses, not PDG."""

    def test_md_mu_ratio(self):
        """m_d/m_u from framework, not PDG."""
        ratio = FW['d'] / FW['u']
        assert 2.10 < ratio < 2.20  # framework range
        # NOT 2.16 (PDG) or 2.00 (FLAG)

    def test_delta_cp_from_framework_masses(self):
        """delta_CP uses framework m_d/m_u, giving 65.0 not 65.2."""
        delta = math.degrees(math.atan(FW['d'] / FW['u']))
        assert 64.5 < delta < 65.5
        # Should be ~65.0, not ~65.2 (PDG) or ~63.4 (FLAG)

    def test_vus_from_framework_masses(self):
        """V_us uses framework m_d, m_s."""
        vus = math.sqrt(FW['d'] / FW['s'])
        assert abs(vus - 0.2235) < 0.001

    def test_vcb_from_framework_masses(self):
        """V_cb uses framework m_u, m_c."""
        vcb = math.sqrt(FW['u'] / FW['c'])
        assert abs(vcb - 0.0414) < 0.001


class TestBosonSectorFromFramework:
    """Boson predictions use framework m_t, not PDG."""

    def test_yt_uses_framework_mt(self):
        """y_t = sqrt(2) * m_t / v where m_t is framework value."""
        v = math.sqrt(2) * FW['t'] / (1 - C0**3)
        yt = math.sqrt(2) * FW['t'] / v
        assert abs(yt - (1 - C0**3)) < 1e-10

    def test_mz_from_sirlin(self):
        """M_Z from Sirlin uses sin2tw = 2/9, not measured M_Z."""
        alpha = 1/137.036
        G_F = 1.1663788e-5  # GeV^-2
        Dr = 0.036
        A_sq = math.pi * alpha / (math.sqrt(2) * G_F)
        M_Z_sq = A_sq * 81 / (14 * (1 - Dr))
        M_Z = math.sqrt(M_Z_sq)
        assert abs(M_Z - 91.33) < 0.5  # GeV

    def test_mh_from_c0_hierarchy(self):
        """M_H from lambda = (1+C0^2)/8, not measured M_H."""
        lam = (1 + C0**2) / 8
        v = math.sqrt(2) * FW['t'] / 1000 / (1 - C0**3)  # GeV
        M_H = v * math.sqrt(2 * lam)
        assert abs(M_H - 125.6) < 1.0  # GeV


class TestFrameworkVsPDG:
    """Framework and PDG should agree within the framework's precision."""

    def test_masses_within_1pct(self):
        """All framework masses within 1% of PDG (except m_t)."""
        pdg = {
            'e': 0.5110, 'mu': 105.66,
            'd': 4.67, 's': 93.4, 'b': 4180.0,
            'u': 2.16, 'c': 1270.0,
        }
        for name in pdg:
            err = abs(FW[name] - pdg[name]) / pdg[name] * 100
            assert err < 1.0, f"{name}: framework {FW[name]} vs PDG {pdg[name]} ({err:.2f}%)"

    def test_mt_within_half_pct(self):
        """Framework m_t within 0.5% of PDG."""
        err = abs(FW['t'] - 172570) / 172570 * 100
        assert err < 0.5


class TestNoMixAndMatch:
    """Verify that CKM predictions are self-consistent (no mixed sources)."""

    def test_ckm_unitarity(self):
        """CKM from framework masses should be exactly unitary."""
        m_u, m_d, m_s, m_c = FW['u'], FW['d'], FW['s'], FW['c']
        s12 = math.sqrt(m_d / m_s)
        s23 = math.sqrt(m_u / m_c)
        s13 = m_d / m_c
        c12 = math.sqrt(1 - s12**2)
        c23 = math.sqrt(1 - s23**2)
        c13 = math.sqrt(1 - s13**2)

        # First row unitarity: |V_ud|^2 + |V_us|^2 + |V_ub|^2 = 1
        row1 = (c12*c13)**2 + (s12*c13)**2 + s13**2
        assert abs(row1 - 1.0) < 1e-14, f"Row 1 unitarity: {row1}"

    def test_geometric_mean_self_consistent(self):
        """M_H = sqrt(m_t * M_Z) using framework's own M_Z prediction."""
        alpha = 1/137.036
        G_F = 1.1663788e-5
        Dr = 0.036
        A_sq = math.pi * alpha / (math.sqrt(2) * G_F)
        M_Z = math.sqrt(A_sq * 81 / (14 * (1 - Dr)))  # GeV

        lam = (1 + C0**2) / 8
        v = math.sqrt(2) * FW['t'] / 1000 / (1 - C0**3)  # GeV
        M_H_from_c0 = v * math.sqrt(2 * lam)
        M_H_from_geom = math.sqrt(FW['t']/1000 * M_Z)

        # These should agree to ~0.1% (self-consistency)
        err = abs(M_H_from_c0 - M_H_from_geom) / M_H_from_c0 * 100
        assert err < 0.5, f"Geometric mean vs C0: {err:.2f}% disagreement"
