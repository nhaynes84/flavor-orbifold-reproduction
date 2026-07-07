#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT / "tdd"))
from framework import boundary  # noqa: E402


def cross_scale_ledger(mode_pressure: float, k_mix: float) -> dict:
    coherence_eff = k_mix / (1.0 + mode_pressure)
    farea = boundary.f_area_fraction()
    fgeom = boundary.f_boost_from_first_principles()
    H0_cmb = 67.4
    H0_local = 73.0
    delta_void = 0.30
    fobs = boundary.hubble_tension_f_boost(H0_cmb, H0_local, delta_void)
    fboost_bridged = fgeom * (1.0 + coherence_eff / farea) * (1.0 + 1.0 / 25.0)
    hubble_gap_bridged = fobs - fboost_bridged

    H0_si = 67.4 * 1000 / 3.085677581e22
    a0 = boundary.mond_acceleration(H0_si)
    a0_obs = 1.2e-10
    a0_bridged = a0 * (1.0 + coherence_eff / (1.0 - farea))
    mond_gap_bridged = (a0_bridged - a0_obs) / a0_obs

    sqrt_As = boundary.cmb_fluctuation_amplitude()
    As_pred = sqrt_As * sqrt_As
    As_obs = 2.1e-9
    cmb_frac_gap = (As_pred - As_obs) / As_obs

    score = abs(hubble_gap_bridged) + abs(mond_gap_bridged) + abs(cmb_frac_gap) + mode_pressure
    return {
        "coherence_eff": coherence_eff,
        "hubble_gap_bridged": hubble_gap_bridged,
        "mond_gap_bridged": mond_gap_bridged,
        "cmb_frac_gap": cmb_frac_gap,
        "bridged_score": score,
    }
