import sys
from pathlib import Path

import pytest


CALC_DIR = Path(__file__).resolve().parents[2] / "calculations"
if str(CALC_DIR) not in sys.path:
    sys.path.append(str(CALC_DIR))

from boundary_model_report import build_report  # noqa: E402


@pytest.mark.framework
@pytest.mark.slow
def test_unified_model_probe_and_bridge_tolerances():
    report = build_report()
    fit = report["fit"]
    probe = fit["probe"]
    bridge = report["bridge"]

    assert fit["k_mix"] == pytest.approx(0.125, abs=1e-12)
    assert probe["ratio"] == pytest.approx(0.6011, abs=0.02)
    assert probe["mu"] == pytest.approx(1.4498, abs=0.05)
    assert probe["r"] == pytest.approx(3.7209, abs=0.6)
    assert probe["w"] == pytest.approx(0.6101, abs=0.05)
    assert probe["split_ratio_secondary"] == pytest.approx(5.2895, abs=0.8)

    assert fit["test_loss"] < 0.35
    assert fit["test_p90_loss"] < 0.5
    assert fit["guard_loss"] < 0.3

    assert report["mode_pressure"] < 0.1
    assert bridge["bridged_score"] < 0.35
