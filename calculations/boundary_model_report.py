#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone

from boundary_model_bridge import cross_scale_ledger
from boundary_model_config import TARGET
from boundary_model_fit import run_fit


def build_report() -> dict:
    fit = run_fit()
    probe = fit["probe"]
    mode_pressure = abs(probe["ratio"] - TARGET["ratio"]) + 0.5 * abs(probe["mu"] - TARGET["mu"])
    bridge = cross_scale_ledger(mode_pressure=mode_pressure, k_mix=fit["k_mix"])
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "fit": fit,
        "mode_pressure": mode_pressure,
        "bridge": bridge,
    }


if __name__ == "__main__":
    print(json.dumps(build_report(), indent=2))
