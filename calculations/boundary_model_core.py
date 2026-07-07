#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
from boundary_ee_network_solver import solve


def state_from_solution(st: dict) -> np.ndarray:
    common = st["F1"]
    differential = st["mismatch"] * st["F2"]
    isospin = st["F2_iso"]
    torsion = st["F2"] * st["H_disp"]
    decoh = 1.0 / (st["overlap"] + 1e-6)
    chirality = st.get("chirality", 0.0)
    anisotropy = abs(st.get("aniso_x", 0.0)) + abs(st.get("aniso_xy", 0.0))
    return np.array([1.0, common, differential, isospin, torsion, decoh, chirality, anisotropy], dtype=float)


def solve_channels(phi: float) -> dict:
    keys = {
        "uud": ["u", "u", "d"],
        "udd": ["u", "d", "d"],
        "uds_p": ["u", "d", "s"],
        "uds_m": ["u", "d", "s"],
        "uus_p": ["u", "u", "s"],
        "uus_m": ["u", "u", "s"],
        "uss_p": ["u", "s", "s"],
        "uss_m": ["u", "s", "s"],
        "sig_p": ["u", "u", "s"],
        "sig_m": ["u", "u", "s"],
        "lam_p": ["u", "d", "s"],
        "lam_m": ["u", "d", "s"],
        "x0_p": ["u", "s", "s"],
        "x0_m": ["u", "s", "s"],
    }
    out = {}
    for k, comp in keys.items():
        sign = -1 if k.endswith("_m") else +1
        out[k] = solve(comp, phi, sign)
    return out
