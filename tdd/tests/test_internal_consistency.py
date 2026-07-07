"""
Internal consistency tests — the honest battleground.

These tests probe real weaknesses and internal tensions in the framework,
not just verify arithmetic. If the framework is wrong, these tests
should be the first to show cracks.
"""

import math
import pytest


# ============================================================
# Framework-predicted masses (Paper 1, Table I)
# ============================================================
M_ELECTRON = 0.5106
M_MUON = 105.69
M_TAU = 1776.86  # anchor
M_UP = 2.172
M_DOWN = 4.657
M_STRANGE = 93.23
M_CHARM = 1270.0      # m_c(m_c)
M_CHARM_2GEV = 1095.0  # m_c at 2 GeV (approximate, for scale tests)
M_BOTTOM = 4181.0
M_TOP = 172782.0

# PDG for comparison
M_UP_PDG = 2.16
M_DOWN_PDG = 4.67


class TestThetа12InternalDiscrepancy:
    """Paper 2 and Paper 3 give different sin²θ₁₂ predictions.
    This is a real internal tension, not a bug."""

    def test_paper2_theta12_from_quarks(self):
        """Paper 2: sin²θ₁₂ = m_c/m_b."""
        sin2_12_paper2 = M_CHARM / M_BOTTOM
        assert abs(sin2_12_paper2 - 0.3038) < 0.001  # ~0.3038

    def test_paper3_theta12_from_leptons(self):
        """Paper 3: sin²θ₁₂ = 1/3 - 2δx."""
        x_mu = math.log(M_MUON / M_ELECTRON) / math.log(M_TAU / M_ELECTRON)
        delta_x = 2.0 / 3 - x_mu
        sin2_12_paper3 = 1.0 / 3 - 2 * delta_x
        assert abs(sin2_12_paper3 - 0.3077) < 0.001  # ~0.3077

    def test_internal_discrepancy_quantified(self):
        """The two derivations disagree by ~1.3%. This is the framework's
        internal precision ceiling — sub-percent claims don't hold for
        all observables."""
        sin2_12_paper2 = M_CHARM / M_BOTTOM  # 0.3038
        x_mu = math.log(M_MUON / M_ELECTRON) / math.log(M_TAU / M_ELECTRON)
        delta_x = 2.0 / 3 - x_mu
        sin2_12_paper3 = 1.0 / 3 - 2 * delta_x  # 0.3077

        gap_pct = abs(sin2_12_paper2 - sin2_12_paper3) / sin2_12_paper2 * 100
        # Gap should be ~1.3% — document it honestly
        assert 0.5 < gap_pct < 3.0, f"θ₁₂ internal gap {gap_pct:.1f}% outside expected range"

    def test_both_closer_to_measurement_than_tbm(self):
        """Both derivations should be closer to measurement (0.303) than TBM (0.333)."""
        measured = 0.303
        tbm = 1.0 / 3
        sin2_12_paper2 = M_CHARM / M_BOTTOM
        x_mu = math.log(M_MUON / M_ELECTRON) / math.log(M_TAU / M_ELECTRON)
        sin2_12_paper3 = 1.0 / 3 - 2 * (2.0 / 3 - x_mu)

        err_tbm = abs(tbm - measured)
        err_p2 = abs(sin2_12_paper2 - measured)
        err_p3 = abs(sin2_12_paper3 - measured)
        assert err_p2 < err_tbm, "Paper 2 θ₁₂ worse than TBM"
        assert err_p3 < err_tbm, "Paper 3 θ₁₂ worse than TBM"


class TestCharmMassScale:
    """The charm mass at 2 GeV vs at m_c scale matters for V_ub prediction."""

    def test_charm_scale_difference(self):
        """m_c(m_c) = 1270 MeV vs m_c(2 GeV) ≈ 1095 MeV — 16% different."""
        gap_pct = abs(M_CHARM - M_CHARM_2GEV) / M_CHARM * 100
        assert gap_pct > 10, "Charm mass scale difference should be >10%"

    def test_vub_changes_with_scale(self):
        """V_ub = m_d/m_c depends on which m_c you use."""
        vub_mc = M_DOWN / M_CHARM       # 0.00368 (at m_c scale)
        vub_2gev = M_DOWN / M_CHARM_2GEV  # 0.00427 (at 2 GeV)
        measured = 0.00382

        # At m_c scale: 3.7% off (as reported)
        err_mc = abs(vub_mc - measured) / measured * 100
        # At 2 GeV: would overshoot by ~12%
        err_2gev = abs(vub_2gev - measured) / measured * 100

        assert err_mc < 5, "V_ub at m_c scale should be within 5%"
        assert err_2gev > 10, "V_ub at 2 GeV should be much worse"

    def test_mass_table_scale_mixing_documented(self):
        """The framework mixes scales: light quarks at 2 GeV, charm at m_c,
        leptons as pole masses. This is common practice but not ideal."""
        # This test just documents the fact
        scales = {
            'u': '2 GeV', 'd': '2 GeV', 's': '2 GeV',
            'c': 'm_c', 'b': 'm_b',
            'e': 'pole', 'mu': 'pole', 'tau': 'pole',
            't': 'pole',
        }
        mixed = len(set(scales.values())) > 1
        assert mixed, "Scales should be mixed (this is the known issue)"


class TestSensitivityToInputs:
    """How fragile is the framework to input mass perturbations within PDG errors?"""

    def test_span_ratio_robust_to_1sigma(self):
        """6/5 span ratio should survive ±1σ perturbation of input masses."""
        # PDG uncertainties
        m_e_1sig = 0.00000015  # negligible
        m_tau_1sig = 0.12
        m_d_1sig = 0.48
        m_b_1sig = 30.0

        target = 6.0 / 5
        worst_error = 0

        for d_sign in [-1, 1]:
            for b_sign in [-1, 1]:
                for t_sign in [-1, 1]:
                    m_d = M_DOWN + d_sign * m_d_1sig
                    m_b = M_BOTTOM + b_sign * m_b_1sig
                    m_tau = M_TAU + t_sign * m_tau_1sig

                    s_lep = math.log(m_tau / M_ELECTRON)
                    s_down = math.log(m_b / m_d)
                    ratio = s_lep / s_down
                    err = abs(ratio - target) / target * 100
                    worst_error = max(worst_error, err)

        # Worst-case 1σ perturbations give ~1.7% — the 0.03% central-value
        # precision does NOT account for input mass uncertainties.
        # The span ratio is robust (stays well within 5%) but the claimed
        # 0.03% is precision of central values, not a 0.03% prediction.
        assert worst_error < 3.0, f"6/5 ratio worst case error {worst_error:.2f}%"

    def test_ckm_angles_robust_to_1sigma(self):
        """CKM predictions from mass ratios should survive ±1σ perturbation."""
        vus_target = 0.22453
        m_d_lo, m_d_hi = M_DOWN - 0.17, M_DOWN + 0.48
        m_s_lo, m_s_hi = M_STRANGE - 3.4, M_STRANGE + 8.6

        worst_vus_err = 0
        for m_d in [m_d_lo, M_DOWN, m_d_hi]:
            for m_s in [m_s_lo, M_STRANGE, m_s_hi]:
                vus = math.sqrt(m_d / m_s)
                err = abs(vus - vus_target) / vus_target * 100
                worst_vus_err = max(worst_vus_err, err)

        # V_us from GST: worst-case ~6.5% at 1σ edges. This is dominated by
        # the huge asymmetric m_d uncertainty (+0.48/-0.17 MeV, 10%).
        # The framework's CKM predictions are NOT sharp tests until
        # light quark masses improve to ~1% precision.
        assert worst_vus_err < 10.0, f"V_us worst case error {worst_vus_err:.1f}%"


class TestAlternativeFunctionalForms:
    """Verify the framework BREAKS with different coupling forms.
    This proves N_q=7 + exponential coupling is specific, not generic."""

    def test_linear_coupling_fails(self):
        """g(q) = 1 + (q-1)*T_F instead of exp((q-1)*T_F).
        Should give wildly wrong masses."""
        T_F = 0.5
        N_c = 3
        N_q = 2 * N_c + 1  # 7
        C_0 = 1.0 / 5

        # Exponential (correct)
        g_lep_exp = math.exp(2 * T_F)
        sigma_lep_exp = 3 * g_lep_exp  # N_l * g

        # Linear (wrong)
        g_lep_lin = 1 + 2 * T_F
        sigma_lep_lin = 3 * g_lep_lin  # N_l * g

        # The span ratio S_lep/S_down depends on sigma ratios
        g_d_exp = 1.0
        sigma_d_exp = N_q * g_d_exp - C_0
        g_d_lin = 1.0
        sigma_d_lin = N_q * g_d_lin - C_0

        ratio_exp = sigma_lep_exp / sigma_d_exp
        ratio_lin = sigma_lep_lin / sigma_d_lin

        # Exponential should give ~1.2 (≈6/5), linear should be far off
        assert abs(ratio_exp - 1.2) < 0.01, f"Exponential ratio {ratio_exp} not near 6/5"
        assert abs(ratio_lin - 1.2) > 0.3, f"Linear ratio {ratio_lin} too close to 6/5"

    def test_wrong_nq_catastrophic(self):
        """Changing N_q from 7 to 6 or 8 should destroy predictions."""
        C_0 = 1.0 / 5
        T_F = 0.5

        for N_q in [5, 6, 8, 9]:
            H = (N_q // 2) ** 2 + ((N_q + 1) // 2) ** 2
            C_0_alt = 1.0 / math.sqrt(H)
            g_d = 1.0
            sigma_d = N_q * g_d - C_0_alt
            g_lep = math.exp(2 * T_F)
            sigma_lep = 3 * g_lep  # N_l always 3

            ratio = sigma_lep / sigma_d
            # Should NOT be anywhere near 6/5 for wrong N_q
            err_from_target = abs(ratio - 1.2) / 1.2 * 100
            assert err_from_target > 5, \
                f"N_q={N_q} gives ratio {ratio:.3f}, too close to 6/5"
