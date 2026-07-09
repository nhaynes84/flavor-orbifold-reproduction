# Reproduction package: quantized-position flavor model

**Archived:** [DOI 10.5281/zenodo.21285294](https://doi.org/10.5281/zenodo.21285294) — papers A/B + prediction registry + this repo (v1.0.0 snapshot), CC-BY-4.0.

Code, tests, significance protocol, and dated prediction registry
supporting:

- **Paper A** — *A charged-lepton mass relation and a quantized-position
  flavor model on an orbifold* (N. Haynes, 2026)
- **Paper B** — *Standard Model gauge coupling ratios from orbifold
  cycle volumes* (N. Haynes, 2026)

Every number in both papers is reproducible from this repository with
Python 3 and no external dependencies beyond `pytest` for the test
suite.

## How the sausage was made

This model was developed test-first: candidate relations were committed
as runnable checks before being trusted, mechanisms that failed
adversarial re-derivation were discarded (the papers' discussion notes
this), and every quantitative claim passed a trials-corrected
significance gate before publication. The full test suite is included
deliberately — it is the development tool, not an afterthought.

## Contents

| Path | What it is |
|------|-----------|
| `scripts/boundary.py` | The model's constants and laws, single source of truth; self-testing (`python3 scripts/boundary.py`) |
| `scripts/verify_consolidation.py` | 44 quantitative claims checked against measurement in one run |
| `scripts/fresh_permutation_test.py` | Paper A's discrete-ensemble landscape statistic (432 variants, canonical = unique minimum) |
| `scripts/look_elsewhere.py` | The trials-correction protocol: expression-grammar enumeration (~29k values) and expected-chance-hit (λ) computation — the machinery behind every significance statement |
| `scripts/rigidity_tests.py` | Coupling-progression uniqueness (geometric √e vs alternative families) and screening-rule joint test |
| `scripts/anomaly_bookkeeping.py` | Per-fixed-point anomaly cancellation, exact fractions, incl. Witten anomaly |
| `scripts/interface_dilutions.py` | Paper B: per-force dilution table, volume pattern, weak angle, α_em anchors |
| `scripts/overdetermination_check.py` | Paper B: the two ratio-crossing scales (descriptive; see Paper B for the statistical status) |
| `scripts/geometry_audit.py` | Enumeration of admissible discrete orbifold substrates |
| `scripts/warp_tower.py`, `ring_running.py`, `leftover_alpha_test.py`, `gravity_mechanism_v2.py` | Supporting analyses, including negative results retained on purpose |
| `tdd/` | The full development test suite (900+ tests): `python3 -m pytest tdd` from the repository root |
| `PREDICTION_REGISTRY_v2.md` | The dated, frozen prediction registry with kill conditions and a status-annotation log (including currently-live tensions, recorded rather than reinterpreted) |

## Quick start

```bash
python3 scripts/boundary.py            # model self-test
python3 scripts/verify_consolidation.py  # all 44 claims
python3 scripts/fresh_permutation_test.py
python3 -m pytest tdd                  # full suite (~2 min)
```

## Provenance

This package is extracted from a private development repository whose
git history dates every result, including the prediction registry's
freeze dates and the discarded mechanisms. The registry's internal
dates are authoritative; the development history is available to
referees on request.
