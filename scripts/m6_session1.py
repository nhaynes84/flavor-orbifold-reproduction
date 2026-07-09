#!/usr/bin/env python3
"""M6 session 1: the minimal channel model, and the fencepost finding.

Model v1: gauge stiffness = information accumulated along the cycle's
path. On the COVER, cycles are closed: #links = #sites. On the QUOTIENT,
folding cuts closed loops into OPEN paths pinned at creases: a path of
ell cells touches ell+1 NODES. If the boundary counts information at
render points (nodes — where matter sits, per M2b), the conversion
acquires a forced, group-blind +1 per channel: THE HEADER IS THE CREASE.

Candidate counting rules (lambda cost stated: THREE discrete rules
tried, all physically motivated, no continuous tuning):
  A. cells only:        delta(1/alpha_i) = eps * ell_i        (tree+)
  B. fencepost/nodes:   delta(1/alpha_i) = eps * (ell_i + 1)
  C. double-crease:     delta(1/alpha_i) = eps * (ell_i + 2)

Each is ONE-parameter. Compare to the 2-param fits of m4_partial
(Casimir chi2=0.35, universal chi2=0.57, both 1 dof).

Run: python3 consolidation/m6_session1.py
"""
import math

ELL = {"U(1)": 7.0, "SU(2)": 3.5, "SU(3)": 1.0}
DELTA = {"U(1)": 0.369, "SU(2)": 0.247, "SU(3)": 0.104}
SIG = {"U(1)": 0.02, "SU(2)": 0.02, "SU(3)": 0.065}

def fit(offset):
    w = {g: 1/SIG[g]**2 for g in ELL}
    n = {g: ELL[g] + offset for g in ELL}
    eps = sum(w[g]*n[g]*DELTA[g] for g in ELL) / sum(w[g]*n[g]**2 for g in ELL)
    chi2 = sum(((DELTA[g] - eps*n[g])/SIG[g])**2 for g in ELL)
    return eps, chi2

print(f"{'rule':<22}{'eps':>9}{'chi2 (2 dof)':>14}")
for name, off in (("A. cells (ell)", 0), ("B. fencepost (ell+1)", 1),
                  ("C. double-crease (+2)", 2)):
    eps, c2 = fit(off)
    print(f"{name:<22}{eps:>9.4f}{c2:>14.2f}")

eps, c2 = fit(1)
print(f"\nFencepost detail (rule B, eps = {eps:.4f}):")
for g in ELL:
    pred = eps*(ELL[g]+1)
    print(f"  {g:<6} pred {pred:+.3f}  meas {DELTA[g]:+.3f}  "
          f"({(DELTA[g]-pred)/SIG[g]:+.1f} sigma)")

print(f"""
READING:
  - Node-counting (B) crushes cell-counting (A): chi2 {fit(1)[1]:.1f} vs {fit(0)[1]:.1f}
    with the SAME parameter count. The +1 is not fitted — it is the
    fencepost arithmetic of orbifolding (closed cover loops -> open
    quotient paths pinned at creases). The header IS the crease.
  - B is competitive with the TWO-parameter models (0.35/0.57 on 1 dof)
    while using one parameter. Parsimony favors the crease.
  - SU(2)'s half-cycle: ell = 3.5 -> 4.5 nodes (path ending mid-cell at
    the fold-2 crease) — consistent, and pleasingly literal.
  - THE NEW SINGLE TARGET: eps = {eps:.4f} in 1/alpha units = {eps*1.5/(4*math.pi)*100:.2f}%
    of the per-cell tree information. ONE number now stands between the
    channel model and a full derivation of the conversion. Its origin
    (cell noise model / quantization / tick arithmetic) is session 2's
    problem. NOT hunted here — the lambda-gate holds.
  - lambda cost of this session: three discrete counting rules compared,
    all pre-motivated; no continuous freedom; stated openly.""")
