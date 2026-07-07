#!/usr/bin/env python3
"""One-loop running of the diagonal gauge coupling on the 7-site ring.

The actual math behind 18-CAPACITY-CLOSURE's lemma attempt. RING moose
(7 sites, 7 links, closed; verifier flag: earlier docstring said
"linear" — cosmetic error). PURE GAUGE — no quarks: this demonstrates
the 11 -> 70 per-cell sum rule, NOT the matter-included beta0 = 7
identity (verifier flag 2). 7 sites of SU(N), 7 complex bifundamental
link scalars with vev v; the
diagonal SU(N) survives; 6 gauge vectors become massive with the discrete
ring dispersion m_n = 2 g v sin(n pi / 7).

We do exact one-loop threshold bookkeeping of d(1/alpha_diag)/dlnQ and
verify the deconstruction sum rule state-by-state — no hand-waving:

  IR (below all masses):   b = 11N/3                      (zero mode only)
  UV (above all masses):   b = 7 * (11N/3 - N/3) = 70N/3  (7 sites, each
                             with 2 adjacent link scalars: -N/3)

State inventory that must connect them (all reps under the DIAGONAL group):
  - 1 massless vector (zero mode):             +11N/3
  - 6 massive vectors (adjoint, levels n=1,2,3 x2): each +21N/6
        (= massless vector 11N/3 + eaten real adjoint scalar -N/6)
  - link scalars: 7 complex bifund = 7 x (complex adjoint + complex singlet)
        = 14 real adjoints + 14 real singlets; 6 real adjoints are eaten
        (counted inside massive vectors) -> 8 physical real adjoints
        (each -N/6) + singlets (0)

Sum rule check: 11N/3 + 6*(21N/6) + 8*(-N/6) = (22 + 126 - 8)N/6 = 70N/3 ?

Run: python3 consolidation/ring_running.py
"""

import math
from fractions import Fraction as F

N = 3                       # fiber dimension (diagonal SU(3))
SITES = 7

# ---- exact coefficient bookkeeping (fractions, no floats) -----------------
b_massless_vector = F(11 * N, 3)
b_massive_vector = F(11 * N, 3) - F(N, 6)     # vector + eaten real adjoint
b_real_adjoint_scalar = -F(N, 6)
b_real_singlet_scalar = F(0)

ir = b_massless_vector
uv_from_states = (b_massless_vector
                  + 6 * b_massive_vector
                  + 8 * b_real_adjoint_scalar
                  + 14 * b_real_singlet_scalar)
uv_from_sites = SITES * (F(11 * N, 3) - F(N, 3))   # per-site: gauge + 2 links

print("EXACT SUM-RULE CHECK (fractions):")
print(f"  IR coefficient (zero mode):            b = {ir}")
print(f"  UV from state inventory:               b = {uv_from_states}")
print(f"  UV from deconstruction sum rule:       b = {uv_from_sites}")
assert uv_from_states == uv_from_sites, "SUM RULE FAILS"
print("  SUM RULE CLOSES EXACTLY. The bookkeeping machinery is sound.")
print()

# ---- the running staircase (numerical, thresholds at ring dispersion) -----
gv = 1.0                    # ring scale in arbitrary units
thresholds = []
for n in (1, 2, 3):         # levels, each doubly degenerate
    m = 2 * gv * math.sin(n * math.pi / SITES)
    thresholds += [(m, b_massive_vector, f"KK vector n={n}")] * 2
# physical link scalars at ~v (model-dependent; take gv, flagged)
thresholds += [(gv, b_real_adjoint_scalar, "link scalar")] * 8
thresholds.sort()

print("RUNNING STAIRCASE d(1/alpha)/dlnQ * 2pi  (b), climbing through thresholds:")
b = ir
print(f"  Q << ring scale:      b = {b}  (pure 4D zero mode = the IR theory)")
last = 0.0
for m, db, label in thresholds:
    b += db
    if m != last:
        print(f"  above Q = {m:0.3f} gv:  b = {b}   (+{label}...)")
        last = m
print(f"  Q >> ring scale:      b = {b} = 70N/3  (per-cell bookkeeping regime)")
print()

# ---- what this does and does not prove -------------------------------------
print("VERDICT ON THE LEMMA:")
print("""  PROVEN here: above the ring scale the flow is EXACTLY per-cell
  bookkeeping (sum rule closes state-by-state) — the UV half of the
  dictionary's reading is standard QFT, verified.
  NOT proven: the IR coefficient (11N/3 + matter) is set by the zero-mode
  content, and NOTHING in standard deconstruction bounds it by the cell
  count. The claim beta0 <= N_cells is NOT standard QFT — it is the
  paradigm's new physics (taxonomy: 'the boundary computes the flow'),
  with the dictionary's integer identities as its empirical support.
  So the floor's last piece is a genuine POSTULATE, now precisely located:
  kinematic fidelity = Nyquist (proven, Shannon); dynamic fidelity =
  capacity bound (postulate, evidence p ~ 1e-3..1e-2, 3 pre-existing
  integer identities). The interface action could still derive it, but
  not via one-loop deconstruction alone.""")
