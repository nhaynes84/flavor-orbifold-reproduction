#!/usr/bin/env python3
"""
Boundary EE network solver (no particle hacks).

Maps the Y-junction to an AC network:
- arm impedance from boundary scale/mode load
- mutual admittance from bridge coupling
- solve steady-state phasors from KCL
- extract overlap/phase observables and channel splits
"""

from __future__ import annotations

import cmath
import math

C0 = 1.0 / 5.0
Q = {"u": 2.0 / 3.0, "d": -1.0 / 3.0, "s": -1.0 / 3.0}
M = {"u": 2.172, "d": 4.657, "s": 93.23}
G = {"u": 1, "d": 1, "s": 2}
I3 = {"u": 0.5, "d": -0.5, "s": 0.0}
# Relative measurement uncertainty proxy (classical volatility projection).
# Smaller = more precise / less volatile observable channel.
# Values are order-of-magnitude anchors, not a PDG ingestion pipeline.
REL_SIGMA = {"u": 0.20, "d": 0.10, "s": 0.03}


def bridge(qi: str, qj: str) -> float:
    dq = abs(Q[qi] * Q[qj])
    dg = abs(G[qi] - G[qj])
    dm = abs(math.log((M[qi] + 1.0) / (M[qj] + 1.0)))
    return (1.0 - C0) * dq * (1.0 + dg) / (1.0 + dm)


def arm_impedance(q: str, w: float) -> complex:
    # R+ jX built from screening + mode scale
    R = (1.0 + C0 * G[q]) * math.log(M[q] + 1.0)
    L = (1.0 - C0 * C0) * (1.0 + 0.4 * G[q])
    C = 1.0 / (1.0 + abs(Q[q]) + 0.2 * G[q])
    X = w * L - 1.0 / (w * C)
    return complex(R, X)


def line_phase(q: str, w: float) -> float:
    # effective arm electrical length from mode scale
    ell = (1.0 + 0.35 * G[q]) * math.log(M[q] + 1.0)
    beta = w * (1.0 - C0 * C0)
    return beta * ell


def reflection_coeff(q: str) -> complex:
    # boundary reflection from screened load mismatch
    rmag = (1.0 - C0) / (2.0 + G[q])
    rph = C0 * math.pi * (1.0 + 0.5 * abs(Q[q]))
    return cmath.rect(rmag, rph)


def solve(comp: list[str], branch_phase: float, mode_sign: int = 1) -> dict:
    n = len(comp)
    w = 1.0 + 0.15 * sum(G[q] for q in comp) / n
    Z = [arm_impedance(q, w) for q in comp]
    # composition-level renormalization from boundary screening/span:
    # higher strange occupancy and larger span increase screening of split-driving mode.
    s_count = sum(1 for q in comp if q == "s")
    g_span = max(G[q] for q in comp) - min(G[q] for q in comp)
    span_mass = max(math.log(M[q] + 1.0) for q in comp) - min(math.log(M[q] + 1.0) for q in comp)
    # Volatility transfer from measurement precision:
    # higher relative uncertainty => larger classical volatility projection.
    vol = sum(REL_SIGMA[q] for q in comp) / n
    renorm = 1.0 + C0 * (0.9 * s_count + 0.5 * g_span + 0.3 * span_mass) + 1.8 * vol
    Yarm = [1.0 / z for z in Z]
    # repeated-flavor coherence projector:
    # identical flavors share local boundary condition and suppress differential mode energy.
    rep_pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if comp[i] == comp[j]:
                rep_pairs.append((i, j))
                # enforce identical local impedance exactly for repeated flavors
                zmean = 0.5 * (Z[i] + Z[j])
                Z[i] = zmean
                Z[j] = zmean
    Yarm = [1.0 / z for z in Z]

    # source phasors at arm boundaries
    src = [cmath.exp(1j * 2 * math.pi * i / n) for i in range(n)]
    # composition-agnostic junction eigenmode perturbation (no per-particle logic)
    for i in range(n):
        sign = mode_sign if (i % 2 == 0) else -mode_sign
        src[i] *= cmath.exp(1j * sign * (branch_phase / renorm))

    # finite-line propagation + reflection dressing at each arm endpoint
    dressed = []
    for i, q in enumerate(comp):
        ph = line_phase(q, w)
        r = reflection_coeff(q)
        tf = cmath.exp(1j * ph) + r * cmath.exp(-1j * ph)
        dressed.append(src[i] * tf)

    # one-node KCL: sum_i Y_i (Vj - Vi) + sum_ij Yij(Vj - Vi)=0
    # simplify to equivalent node from arm branches + mutual bridge contribution
    num = 0 + 0j
    den = 0 + 0j
    for i in range(n):
        y_eff = Yarm[i]
        for j in range(n):
            if i == j:
                continue
            y_eff += complex(bridge(comp[i], comp[j]) / renorm, 0.0)
        num += y_eff * dressed[i]
        den += y_eff
    vj = num / den if den != 0 else 0 + 0j

    arm = [dressed[i] - vj for i in range(n)]
    # apply coherence cancellation on repeated-flavor differential mode
    for i, j in rep_pairs:
        # dynamic repeated-pair coherence:
        # stronger when pair is truly degenerate and weakly contrasted to the 3rd arm,
        # weaker when pair-vs-third mismatch is large (allows nonzero differential mode).
        kpair = bridge(comp[i], comp[j])
        others = [k for k in range(n) if k not in (i, j)]
        if others:
            k3 = 0.5 * (bridge(comp[i], comp[others[0]]) + bridge(comp[j], comp[others[0]]))
            m_pair = math.log(M[comp[i]] + 1.0)
            m_3 = math.log(M[comp[others[0]]] + 1.0)
            contrast = abs(m_pair - m_3)
        else:
            k3 = kpair
            contrast = 0.0
        # bounded 0..1, typically mid-strength
        cpair = (1.0 - C0) * (kpair / (kpair + k3 + 1e-12)) * (1.0 / (1.0 + 0.35 * contrast))
        # volatility weakens perfect locking in fuzzy channels
        vol_pair = 0.5 * (REL_SIGMA[comp[i]] + REL_SIGMA[comp[j]])
        cpair *= 1.0 / (1.0 + 2.2 * vol_pair)
        cpair = max(0.15, min(0.85, cpair))
        vcm = 0.5 * (arm[i] + arm[j])
        vdf = 0.5 * (arm[i] - arm[j])
        arm[i] = vcm + (1.0 - cpair) * vdf
        arm[j] = vcm - (1.0 - cpair) * vdf
    mags = [abs(v) for v in arm]
    ph = [cmath.phase(v) for v in arm]
    # observables
    overlap = abs(sum(arm) / n)
    mismatch = 0.0
    chirality = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = ((ph[j] - ph[i] + math.pi) % (2 * math.pi)) - math.pi
            mismatch += d * d
    for i in range(n):
        j = (i + 1) % n
        d = ((ph[j] - ph[i] + math.pi) % (2 * math.pi)) - math.pi
        chirality += math.sin(d)
    mismatch = math.sqrt(mismatch / max(1, n * (n - 1) / 2))
    # anisotropic orientation tensor components from arm phasors
    # (captures directional projection geometry at the junction)
    t_xx = 0.0
    t_yy = 0.0
    t_xy = 0.0
    for i in range(n):
        c = math.cos(ph[i])
        s = math.sin(ph[i])
        wgt = mags[i]
        t_xx += wgt * c * c
        t_yy += wgt * s * s
        t_xy += wgt * c * s
    tr = t_xx + t_yy + 1e-12
    aniso_x = (t_xx - t_yy) / tr
    aniso_xy = (2.0 * t_xy) / tr
    # Observable energy includes volatility projection weight
    energy = sum((abs(arm[i]) ** 2) * Z[i].real * (1.0 + 3.0 * REL_SIGMA[comp[i]]) for i in range(n))
    # form-factor style projections from phasor currents:
    # F1: charge-like coherent sum
    # F2: boundary-derived loop circulation (bridge-weighted graph current)
    F1 = abs(sum(Q[comp[i]] * arm[i] for i in range(n)))
    loop = 0.0
    loop_iso = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            kij = bridge(comp[i], comp[j])
            dphi = ((ph[j] - ph[i] + math.pi) % (2 * math.pi)) - math.pi
            # charge-odd circulation on pair edge
            qodd = abs(Q[comp[i]] - Q[comp[j]])
            loop += kij * qodd * (mags[i] + mags[j]) * math.sin(dphi)
            iodd = abs(I3[comp[i]] - I3[comp[j]])
            loop_iso += kij * iodd * (mags[i] + mags[j]) * math.sin(dphi)
    # screening-normalized transverse channel + dispersive transfer
    F2_raw = abs(loop) * (1.0 - C0 * C0)
    F2_iso_raw = abs(loop_iso) * (1.0 - C0 * C0)
    # effective network resonance and damping from arm impedances
    x_abs = [abs(z.imag) for z in Z]
    r_abs = [abs(z.real) for z in Z]
    xbar = sum(x_abs) / n
    rbar = sum(r_abs) / n
    w0 = 1.0 + 0.25 * sum(G[q] for q in comp) / n
    gamma = (rbar / (xbar + 1e-9)) * (0.35 + 0.65 * C0)
    dw = w - w0
    # Lorentz-like transfer magnitude
    H = 1.0 / math.sqrt((1.0 - (dw / (w0 + 1e-9)) ** 2) ** 2 + gamma**2)
    # bounded dispersive gain to avoid blow-up
    H = max(0.35, min(2.5, H))
    F2 = F2_raw * H
    F2_iso = F2_iso_raw * (0.85 + 0.15 * H)
    return {
        "energy": energy,
        "overlap": overlap,
        "mismatch": mismatch,
        "mags": mags,
        "ph": ph,
        "F1": F1,
        "F2": F2,
        "F2_iso": F2_iso,
        "H_disp": H,
        "chirality": chirality,
        "aniso_x": aniso_x,
        "aniso_xy": aniso_xy,
    }


def split_probe(name: str, comp: list[str], phi: float = 0.35) -> tuple[float, float]:
    s1 = solve(comp, phi, +1)
    s2 = solve(comp, phi, -1)
    dE = abs(s2["energy"] - s1["energy"])
    dO = abs(s2["overlap"] - s1["overlap"])
    print(f"{name}: dE={dE:.6f} dOverlap={dO:.6f} E+={s1['energy']:.6f} E-={s2['energy']:.6f}")
    return dE, dO


def main() -> None:
    print("EE boundary network split probe")
    de_uds, _ = split_probe("uds", ["u", "d", "s"])
    de_uss, _ = split_probe("uss", ["u", "s", "s"])
    if de_uds > 0:
        print(f"ratio uss/uds = {de_uss/de_uds:.3f} (target ~0.612)")

    print("\nphase sweep")
    for p in [0.15, 0.25, 0.35, 0.45, 0.55]:
        d1, _ = split_probe(f"uds(phi={p:.2f})", ["u", "d", "s"], p)
        d2, _ = split_probe(f"uss(phi={p:.2f})", ["u", "s", "s"], p)
        r = d2 / d1 if d1 > 0 else 0.0
        print(f"  phi={p:.2f} ratio={r:.3f}")


if __name__ == "__main__":
    main()
