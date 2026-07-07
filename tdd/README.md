# Boundary Framework — TDD Test Suite

Real TDD for physics research. Every test calls a framework function.
Red tests ARE the research roadmap.

## Running Tests

```bash
# Full scorecard
pytest

# Framework predictions only (the real scorecard)
pytest -m framework

# Data validation only (should be all green)
pytest -m data

# Quick summary
pytest -m framework --tb=line
```

## Current Scorecard

| Category | Passed | Failed | Total |
|----------|--------|--------|-------|
| **Data validation** (`-m data`) | 59 | 0 | 59 |
| **Framework predictions** (`-m framework`) | 144 | 26 | 170 |
| **Total** | 203 | 26 | 229 |

### What's GREEN (framework predictions that work)

- **Tau/electron = L(37,31) = 3477** — exact within 1σ (the big one)
- **MOND acceleration a₀ = cH₀/(2π)** — derived, no free parameters, within ~10%
- **SPARC rotation curves** — framework a₀ fits 171 galaxies within 4.5% of empirical (zero free params)
- **a₀ self-consistency** — local H₀ → local a₀ within 6% of empirical; same f_boost connects both
- **Hubble tension f_boost** — three independent datasets converge to ~1.73
- **QNM low quadrupole** — ℓ=2 has lowest quality factor → most suppressed CMB mode (structural)
- **Kerr spin a*≈0.1** — CMB hemispheric asymmetry A=0.07 → parent BH slowly rotating
- **Ringdown fails** — constant ω_I can't selectively suppress ℓ=2 → geometric leakage required
- **Flat Sachs-Wolfe from 1/r²** — constant damping rates across ℓ = inverse-square gravity signature
- **Planck unit grid identities** — c = Lp/Tp, Lp×Mp = ħ/c, etc.
- **Loeschian lattice math** — eigenvalues, statistics, mixed-sign modes
- **Vacuum dressing trend** — deviations decrease with mass (running couplings)
- **Electron loop size** — Compton wavelength in Planck pixels
- **Buoyancy framework** — Newton, MOND transition, GW speed, BH saturation
- **Spin topology** — figure-8 = spin ½, circle = spin 1, dent = spin 0

### What's RED (the research roadmap)

Each failed test maps to a specific formula that needs to be derived:

1. **ISW excess amplitude** — predicts ~1.7× but observe ~5× (nonlinear effects?)
2. **Vacuum dressing formula** — 9 particle mass ratios need corrections to bare Loeschian
3. **Muon/electron ratio** — suspected 3/(2α) + QED corrections
4. **Fine structure constant** — derive α = 1/137 from geometry (the holy grail)
5. **Strong coupling α_s** — from boundary stitching topology
6. **Weinberg angle** — sin²θ_W from boundary geometry
7. **Selection rules** — which Loeschian modes are stable particles?
8. **Proton stitching energy** — ~929 MeV from boundary topology
9. **Neutron-proton mass difference** — 1.2933 MeV
10. **Dark energy** — Ω_Λ, w₀, wₐ from differential boundary growth
11. **CMB fluctuation amplitude** — δT/T from boundary granularity
12. **f_boost derivation** — first-principles (not just measured)
13. **Born rule** — P(k) = n_k/N with testable predictions
14. **Cold spot boundary contribution** — quantitative boundary thinning model
15. **Hawking leak** — non-thermal boundary stress spectral distortion
16. **Quadrupole leakage fraction** — derive 80% ℓ=2 leakage from geometry
17. **Leakage suppression model** — S(ℓ) = 1 - A×exp(-βℓ) from first principles
18. **Cold Spot Y₂₀ extremum** — statistical test of Cold Spot at quadrupole maximum

## Philosophy

- **Zero `pytest.skip()` calls.** If we don't have the math, the function raises `NotImplementedError` and the test FAILS. Red means "work to do."
- **Every test calls a framework function.** No `NIST == NIST` tautologies.
- **Data validation is separate.** `@pytest.mark.data` tests confirm our reference values are right. `@pytest.mark.framework` tests are the real scorecard.
- **No faking.** If `boundary.py` doesn't have the formula, it raises `NotImplementedError`. Period.
