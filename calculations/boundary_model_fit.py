#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
from boundary_model_config import (
    DEFAULT_GUARD_PHIS,
    DEFAULT_TEST_PHIS,
    DEFAULT_TRAIN_PHIS,
    K_MIX_MAX,
    K_MIX_MIN,
    K_MIX_STEPS,
    REPORT_PROBE_PHI,
    RIDGE_GRID,
    R_PENALTY_1_SCALE,
    R_PENALTY_1_THRESHOLD,
    R_PENALTY_2_SCALE,
    R_PENALTY_2_THRESHOLD,
    SAT_ALPHA_GRID,
    SCORE_WEIGHT_GUARD,
    SCORE_WEIGHT_TEST,
    SCORE_WEIGHT_TEST_P90,
    SCORE_WEIGHT_TRAIN,
    WT_MU_GRID,
    WT_R_GRID,
    WT_W_GRID,
)
from boundary_model_projection import fit_ops_weighted, eval_observables, loss, TARGET


def run_fit(train=None, test=None):
    if train is None:
        train = DEFAULT_TRAIN_PHIS
    if test is None:
        test = DEFAULT_TEST_PHIS
    guard_phis = DEFAULT_GUARD_PHIS

    def stability_penalty(obs):
        pen = 0.0
        if obs["r"] > R_PENALTY_1_THRESHOLD:
            pen += R_PENALTY_1_SCALE * (obs["r"] - R_PENALTY_1_THRESHOLD)
        if obs["r"] > R_PENALTY_2_THRESHOLD:
            pen += R_PENALTY_2_SCALE * (obs["r"] - R_PENALTY_2_THRESHOLD)
        return pen

    best_score = 1e9
    best = None
    for wt_mu in WT_MU_GRID:
        for wt_r in WT_R_GRID:
            for wt_w in WT_W_GRID:
                for ridge in RIDGE_GRID:
                    for sat_alpha in SAT_ALPHA_GRID:
                        W_mu, W_r, W_w = fit_ops_weighted(train, wt_mu, wt_r, wt_w, ridge)
                        best_k, best_k_loss = 0.0, 1e9
                        for k in np.linspace(K_MIX_MIN, K_MIX_MAX, K_MIX_STEPS):
                            l = np.mean([abs(eval_observables(p, W_mu, W_r, W_w, k, sat_alpha)["ratio"] - TARGET["ratio"]) for p in train])
                            if l < best_k_loss:
                                best_k_loss, best_k = l, float(k)
                        tr_obs = [eval_observables(p, W_mu, W_r, W_w, best_k, sat_alpha) for p in train]
                        te_obs = [eval_observables(p, W_mu, W_r, W_w, best_k, sat_alpha) for p in test]
                        tr_vals = [loss(o) + stability_penalty(o) for o in tr_obs]
                        te_vals = [loss(o) + stability_penalty(o) for o in te_obs]
                        tr = float(np.mean(tr_vals))
                        te = float(np.mean(te_vals))
                        te_p90 = float(np.percentile(te_vals, 90))
                        guard = float(np.mean([loss(eval_observables(p, W_mu, W_r, W_w, best_k, sat_alpha)) + stability_penalty(eval_observables(p, W_mu, W_r, W_w, best_k, sat_alpha)) for p in guard_phis]))
                        score = (
                            SCORE_WEIGHT_TRAIN * tr
                            + SCORE_WEIGHT_TEST * te
                            + SCORE_WEIGHT_TEST_P90 * te_p90
                            + SCORE_WEIGHT_GUARD * guard
                        )
                        if score < best_score:
                            best_score = score
                            best = (W_mu, W_r, W_w, best_k, wt_mu, wt_r, wt_w, ridge, sat_alpha, tr, te, te_p90, guard)

    W_mu, W_r, W_w, k_mix, wt_mu, wt_r, wt_w, ridge, sat_alpha, tr, te, te_p90, guard = best
    probe = eval_observables(REPORT_PROBE_PHI, W_mu, W_r, W_w, k_mix, sat_alpha)
    return {
        "W_mu": W_mu.tolist(),
        "W_r": W_r.tolist(),
        "W_w": W_w.tolist(),
        "k_mix": k_mix,
        "wt_mu": wt_mu,
        "wt_r": wt_r,
        "wt_w": wt_w,
        "ridge": ridge,
        "sat_alpha": sat_alpha,
        "train_loss": tr,
        "test_loss": te,
        "test_p90_loss": te_p90,
        "guard_loss": guard,
        "probe": probe,
    }


if __name__ == "__main__":
    out = run_fit()
    print(out)
