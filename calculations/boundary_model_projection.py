#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
from boundary_model_config import EXTRA_TARGET, F_AREA, TARGET
from boundary_model_core import state_from_solution, solve_channels


def fit_ops_weighted(phi_grid, wt_mu: float, wt_r: float, wt_w: float, ridge: float):
    X_mu, y_mu, X_r, y_r, X_w, y_w = [], [], [], [], [], []
    for phi in phi_grid:
        s = solve_channels(phi)
        X_mu += [state_from_solution(s["uud"]), state_from_solution(s["udd"])]
        y_mu += [TARGET["mu"], 1.0]
        X_r += [state_from_solution(s["uud"]), state_from_solution(s["udd"])]
        y_r += [TARGET["r"], 1.0]
        X_w += [state_from_solution(s["uds_p"]), state_from_solution(s["uus_p"])]
        y_w += [TARGET["w"], 1.0]

    def solve_w(X_rows, y_vals, wt):
        X = np.vstack(X_rows)
        y = np.array(y_vals)
        W = np.diag([wt] * len(y))
        I = np.eye(X.shape[1])
        return np.linalg.solve(X.T @ W @ X + ridge * I, X.T @ W @ y)

    return solve_w(X_mu, y_mu, wt_mu), solve_w(X_r, y_r, wt_r), solve_w(X_w, y_w, wt_w)


def mixed_split_ratio(s: dict, k_mix: float) -> float:
    d_uds = abs(s["uds_p"]["energy"] - s["uds_m"]["energy"])
    d_uss_raw = abs(s["uss_p"]["energy"] - s["uss_m"]["energy"])
    att = 1.0 / (1.0 + k_mix * (s["uss_p"]["mismatch"] + s["uss_m"]["mismatch"]) * 0.5)
    return (d_uss_raw * att) / (d_uds + 1e-12)


def eval_observables(phi: float, W_mu, W_r, W_w, k_mix: float, sat_alpha: float = 0.04):
    s = solve_channels(phi)
    mu = abs((state_from_solution(s["uud"]) @ W_mu) / ((state_from_solution(s["udd"]) @ W_mu) + 1e-12))
    r_raw = abs((state_from_solution(s["uud"]) @ W_r) / ((state_from_solution(s["udd"]) @ W_r) + 1e-12))
    sat = 1.0 + sat_alpha * r_raw * (phi / 0.56) ** 4
    r = r_raw / sat
    w = abs((state_from_solution(s["uds_p"]) @ W_w) / ((state_from_solution(s["uus_p"]) @ W_w) + 1e-12))
    ratio = mixed_split_ratio(s, k_mix)

    d_sig = abs(s["sig_p"]["energy"] - s["sig_m"]["energy"])
    d_lam = abs(s["lam_p"]["energy"] - s["lam_m"]["energy"])
    d_x0 = abs(s["x0_p"]["energy"] - s["x0_m"]["energy"])
    strange_scale = (93.23 - 4.657) * (1.0 / 5.0)
    split_sigma_lambda = abs(d_sig - d_lam) * strange_scale
    split_xi_lambda = abs(d_x0 - d_lam) * strange_scale

    def ch_scale(a, b):
        num = (a["F2_iso"] + b["F2_iso"]) * 0.5
        den = (a["F1"] + b["F1"]) * 0.5 + 1e-12
        mode_ratio = num / den
        decoh = (a["mismatch"] + b["mismatch"]) * 0.5
        return 1.0 + (k_mix / F_AREA) * mode_ratio * (1.0 + 0.5 * decoh)

    s_sl = ch_scale(s["sig_p"], s["lam_p"])
    s_xl = ch_scale(s["x0_p"], s["lam_p"])
    an_sl = abs(s["sig_p"].get("aniso_x", 0.0) - s["lam_p"].get("aniso_x", 0.0)) + abs(
        s["sig_p"].get("aniso_xy", 0.0) - s["lam_p"].get("aniso_xy", 0.0)
    )
    an_xl = abs(s["x0_p"].get("aniso_x", 0.0) - s["lam_p"].get("aniso_x", 0.0)) + abs(
        s["x0_p"].get("aniso_xy", 0.0) - s["lam_p"].get("aniso_xy", 0.0)
    )
    ch_sl = abs(s["sig_p"].get("chirality", 0.0) - s["lam_p"].get("chirality", 0.0))
    ch_xl = abs(s["x0_p"].get("chirality", 0.0) - s["lam_p"].get("chirality", 0.0))
    a_sl = 1.0 + 0.55 * an_sl + 0.25 * ch_sl
    a_xl = 1.0 + 0.55 * an_xl + 0.25 * ch_xl
    split_ratio_secondary = (split_sigma_lambda * s_sl * a_sl) / (split_xi_lambda * s_xl * a_xl + 1e-12)

    return {"ratio": ratio, "mu": mu, "r": r, "w": w, "split_ratio_secondary": split_ratio_secondary}


def loss(obs: dict) -> float:
    return (
        abs(obs["ratio"] - TARGET["ratio"])
        + 0.9 * abs(obs["mu"] - TARGET["mu"])
        + 0.8 * abs(obs["r"] - TARGET["r"])
        + 0.6 * abs(obs["w"] - TARGET["w"])
        + 0.25 * abs(obs["split_ratio_secondary"] - EXTRA_TARGET["split_ratio_secondary"])
    )
