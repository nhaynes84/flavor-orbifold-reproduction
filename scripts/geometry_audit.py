#!/usr/bin/env python3
"""The single-geometry audit (23-COVER-QUOTIENT's gate question).

Enumerate discrete orbifold candidates: Z_N cover with one fold (n -> -n)
or two successive folds, plus the fiber. For each, compute FORCED
properties (no fitting):

  R1: matter quotient has TWO endpoints (gen-1 at 0, gen-3 at 1 are both
      load-bearing) <=> the relevant fold has two fixed cells.
  R2: the successive-quotient volume menu contains the Cycle Rule's
      gauge volumes {7, 7/2} (fiber supplies the 1) as FORCED quotient
      volumes: Vol(Z_N/Z_2^j) = N/2^j.
  R3: the matter wave (k = N_q = 7 sin mode; corpus: sin(7 pi x)) is
      representable: on Z_N, distinct sin modes are k = 1..ceil(N/2)-1,
      so k=7 needs N >= 15 (on Z_14, sin_7 aliases to zero identically).

Run: python3 consolidation/geometry_audit.py
"""

import math

def fixed_cells(N):
    """Fixed cells of n -> -n mod N: solutions of 2n = 0 mod N."""
    return [n for n in range(N) if (2 * n) % N == 0]

def sin7_representable(N):
    """Is sin(2 pi * 7 * n / N) not identically zero on Z_N?"""
    return any(abs(math.sin(2 * math.pi * 7 * n / N)) > 1e-12 for n in range(N))

print(f"{'N':>3} {'folds':>5} {'fixed':>6} {'endpoints':>9} "
      f"{'vol menu (N/2^j)':>22} {'sin7?':>6}  R1 R2 R3")
print("-" * 78)
survivors = []
for N in range(6, 31):
    fx = fixed_cells(N)
    endpoints = len(fx)                    # first-fold fixed cells
    menu = []
    v = N
    for j in range(1, 4):
        v = v / 2
        menu.append(v)
    r1 = endpoints == 2
    r2 = (7 in menu) and (3.5 in menu)     # forced quotient volumes only
    r3 = sin7_representable(N)
    flag = ""
    if r1 and r2 and r3:
        flag = "  <-- FULL SURVIVOR"
        survivors.append((N, "discrete matter"))
    elif r1 and r2:
        flag = "  <-- survives with CONTINUUM matter only"
        survivors.append((N, "continuum matter"))
    print(f"{N:>3} {'Z2^j':>5} {len(fx):>6} {endpoints:>9} "
          f"{str(menu):>22} {str(r3):>6}  {'Y' if r1 else '.'}  "
          f"{'Y' if r2 else '.'}  {'Y' if r3 else '.'}{flag}")

print()
print("AUDIT VERDICT:")
if not survivors:
    print("  NO geometry satisfies R1+R2 (even with continuum matter).")
for N, mode in survivors:
    print(f"  N = {N} cover, successive Z2 folds, {mode}.")
print("""
NOTES (forced facts, no fitting):
- R1 (two endpoints) forces EVEN N. R2 (menu contains 7 and 3.5) forces
  N = 14 (j=1: 7, j=2: 3.5) or N = 28 (j=2: 7, j=3: 3.5). R3 (discrete
  sin_7) forces N >= 15. R1+R2+R3 discrete: N = 28 is the ONLY candidate.
- N = 14: satisfies R1+R2 but sin_7 aliases to ZERO on Z_14 (exactly at
  the Nyquist edge) — matter must live on the CONTINUUM quotient, with
  cells as gauge bookkeeping only.
- N = 28 check: fixed cells {0, 14} (two endpoints ✓); menu {14, 7, 3.5}
  — 7 and 3.5 are the j=2 and j=3 quotients; sin_7 representable ✓.
  Cost: the gauge loci become the SECOND and THIRD quotients (two extra
  folds' worth of structure to justify) and the cover volume 28 = 4N_q
  needs a meaning. More structure, all three requirements met.
- Either way the ODD 7-RING of the campaign docs (12, 13, 21) is DEAD as
  the literal substrate: it fails R1 (one endpoint). The Cycle Rule's
  VOLUMES survive (they are quotient volumes of the surviving covers);
  its 21-era locus labels (cover=7) do not.""")
