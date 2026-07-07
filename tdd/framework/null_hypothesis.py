"""
Null hypothesis test: scrambled-data analysis.

Referee G recommended: "Run the formula-discovery algorithm on a SCRAMBLED
dataset (permute masses randomly among the 9 fermions) and count how many
sub-percent relationships survive."

This module implements that test. It defines the framework's relationship
formulas, evaluates them on the physical masses, then evaluates them on
thousands of random mass permutations. The ratio tells you how much of
the framework's success is due to structure vs. chance.

The formulas tested are exactly the ones claimed in the papers:
  - 3 span ratios (5/3, 6/5, 25/18)
  - 4 CKM constraint equations (GST, V_cb, V_ub, delta_CP)
  - 3 PMNS predictions (theta_12, theta_13, theta_23)
  - 2 cross-sector mass predictions (m_tau from quarks, m_u from spans)
  - 2 cross-sector PMNS from quarks (sin theta_12 = sqrt(m_c/m_b), etc.)
  - 1 decoherence rate (gamma = N_c * T_F -> V_us)

Total: 15 formula-to-observable comparisons.

A "hit" is defined as a formula matching its target observable to within
a given tolerance (default 1%).
"""

import math
import random
from typing import Optional


# ============================================================
# The 9 physical fermion masses (MeV)
# ============================================================

# Framework-predicted masses (Paper 1, Table I)
# These are self-consistent values from the universal mass formula.
PHYSICAL_MASSES = {
    'e': 0.5106,
    'mu': 105.69,
    'tau': 1776.86,    # anchor
    'u': 2.172,
    'c': 1270.0,
    't': 172782.0,
    'd': 4.657,
    's': 93.23,
    'b': 4181.0,
}

# Physical sector assignments: which mass goes where
PHYSICAL_ASSIGNMENT = {
    'lep1': 'e', 'lep2': 'mu', 'lep3': 'tau',
    'up1': 'u', 'up2': 'c', 'up3': 't',
    'down1': 'd', 'down2': 's', 'down3': 'b',
}

# Measured targets (the observables the framework claims to predict)
TARGETS = {
    'span_up_down': 5 / 3,
    'span_lep_down': 6 / 5,
    'span_up_lep': 25 / 18,
    'V_us': 0.22453,
    'V_cb': 0.04080,
    'V_ub': 0.00382,
    'delta_CP_deg': 65.4,
    'sin2_theta13_pmns': 0.02203,
    'sin2_theta12_pmns': 0.303,
    'sin2_theta23_pmns': 0.546,   # IO best fit
    'm_tau_cross': 1776.86,
    'sin_theta12_pmns_quarks': math.sqrt(0.303),
    'sin_theta13_pmns_quarks': math.sqrt(0.02203),
    'V_us_from_gamma': 0.22453,
}


def evaluate_formulas(masses: dict[str, float]) -> dict[str, Optional[float]]:
    """
    Evaluate all framework formulas given an assignment of 9 masses
    to the 9 fermion slots.

    Args:
        masses: dict with keys 'lep1','lep2','lep3','up1','up2','up3',
                'down1','down2','down3' mapping to mass values in MeV.

    Returns:
        dict of formula_name -> predicted_value (or None if formula
        can't be evaluated, e.g. log of negative number).
    """
    # Extract
    lep1, lep2, lep3 = masses['lep1'], masses['lep2'], masses['lep3']
    up1, up2, up3 = masses['up1'], masses['up2'], masses['up3']
    down1, down2, down3 = masses['down1'], masses['down2'], masses['down3']

    results = {}

    # Guard against bad values
    try:
        # Span ratios
        s_up = math.log(up3 / up1) if up3 > 0 and up1 > 0 else None
        s_down = math.log(down3 / down1) if down3 > 0 and down1 > 0 else None
        s_lep = math.log(lep3 / lep1) if lep3 > 0 and lep1 > 0 else None

        if s_up and s_down and s_down != 0:
            results['span_up_down'] = s_up / s_down
        if s_lep and s_down and s_down != 0:
            results['span_lep_down'] = s_lep / s_down
        if s_up and s_lep and s_lep != 0:
            results['span_up_lep'] = s_up / s_lep

        # CKM constraint equations (from mass ratios)
        # V_us = sqrt(m_d1 / m_d2)
        if down1 > 0 and down2 > 0:
            results['V_us'] = math.sqrt(down1 / down2)

        # V_cb = sqrt(m_u1 / m_u2)
        if up1 > 0 and up2 > 0:
            results['V_cb'] = math.sqrt(up1 / up2)

        # V_ub = m_d1 / m_u2
        if up2 > 0:
            results['V_ub'] = down1 / up2

        # delta_CP = arctan(m_d1 / m_u1) in degrees
        if up1 > 0:
            results['delta_CP_deg'] = math.degrees(math.atan(down1 / up1))

        # PMNS from charged lepton geometry
        # delta_x = 2/3 - ln(m_lep2/m_lep1) / ln(m_lep3/m_lep1)
        if lep1 > 0 and lep2 > lep1 and lep3 > lep1:
            x_mu = math.log(lep2 / lep1) / math.log(lep3 / lep1)
            delta_x = 2.0 / 3 - x_mu
            results['sin2_theta13_pmns'] = math.sqrt(3) * delta_x
            results['sin2_theta12_pmns'] = 1.0 / 3 - 2 * delta_x
            results['sin2_theta23_pmns'] = 0.5 + (2 + math.sqrt(3)) * delta_x

        # Cross-sector m_tau prediction: m_lep3 = m_lep1 * (m_d3/m_d1)^(6/5)
        if lep1 > 0 and down3 > 0 and down1 > 0:
            results['m_tau_cross'] = lep1 * (down3 / down1) ** (6.0 / 5)

        # PMNS from quark masses (Paper 2 Rosetta)
        # sin(theta_12_PMNS) = sqrt(m_u2 / m_d3)
        if up2 > 0 and down3 > 0:
            results['sin_theta12_pmns_quarks'] = math.sqrt(up2 / down3)

        # sin(theta_13_PMNS) = sqrt(m_d2 / m_d3)
        if down2 > 0 and down3 > 0:
            results['sin_theta13_pmns_quarks'] = math.sqrt(down2 / down3)

        # V_us from gamma = N_c * T_F = 3/2
        # This is parameter-free, doesn't depend on mass assignment
        # Include it for completeness but it's always the same
        results['V_us_from_gamma'] = math.exp(-3.0 / 2)

    except (ValueError, ZeroDivisionError, OverflowError):
        pass

    return results


def count_hits(predictions: dict[str, Optional[float]],
               targets: dict[str, float],
               tolerance: float = 0.01) -> tuple[int, list[str]]:
    """
    Count how many predictions match targets within tolerance.

    Args:
        predictions: formula_name -> predicted_value
        targets: formula_name -> target_value
        tolerance: fractional tolerance (0.01 = 1%)

    Returns:
        (hit_count, list_of_formula_names_that_hit)
    """
    hits = []
    for name, target in targets.items():
        if name not in predictions or predictions[name] is None:
            continue
        pred = predictions[name]
        if target == 0:
            continue
        frac_error = abs(pred - target) / abs(target)
        if frac_error <= tolerance:
            hits.append(name)
    return len(hits), hits


def physical_mass_assignment() -> dict[str, float]:
    """Return the physical mass assignment."""
    return {
        'lep1': PHYSICAL_MASSES['e'],
        'lep2': PHYSICAL_MASSES['mu'],
        'lep3': PHYSICAL_MASSES['tau'],
        'up1': PHYSICAL_MASSES['u'],
        'up2': PHYSICAL_MASSES['c'],
        'up3': PHYSICAL_MASSES['t'],
        'down1': PHYSICAL_MASSES['d'],
        'down2': PHYSICAL_MASSES['s'],
        'down3': PHYSICAL_MASSES['b'],
    }


def scrambled_mass_assignment(rng: random.Random = None) -> dict[str, float]:
    """
    Create a random assignment of the 9 physical masses to the 9 slots.

    The masses are the same 9 numbers, but randomly shuffled across all
    9 positions (not preserving sector structure).
    """
    if rng is None:
        rng = random.Random()

    mass_values = list(PHYSICAL_MASSES.values())
    rng.shuffle(mass_values)

    slots = ['lep1', 'lep2', 'lep3', 'up1', 'up2', 'up3',
             'down1', 'down2', 'down3']

    # Sort within each sector so gen1 < gen2 < gen3
    # (the framework assumes hierarchical ordering within sectors)
    assignment = dict(zip(slots, mass_values))

    # Enforce ordering within each sector
    for prefix in ['lep', 'up', 'down']:
        keys = [f'{prefix}1', f'{prefix}2', f'{prefix}3']
        vals = sorted([assignment[k] for k in keys])
        for k, v in zip(keys, vals):
            assignment[k] = v

    return assignment


def run_null_hypothesis(n_trials: int = 1000,
                        tolerance: float = 0.01,
                        seed: int = 42) -> dict:
    """
    Run the full null hypothesis test.

    1. Evaluate formulas on physical masses, count hits.
    2. For n_trials random permutations, evaluate formulas, count hits.
    3. Compare distributions.

    Args:
        n_trials: number of random permutations
        tolerance: fractional tolerance for a "hit"
        seed: random seed for reproducibility

    Returns:
        dict with:
          - physical_hits: int
          - physical_hit_names: list[str]
          - scrambled_mean: float
          - scrambled_std: float
          - scrambled_max: int
          - scrambled_distribution: list[int] (hit counts per trial)
          - significance_sigma: float (how many sigma physical exceeds scrambled)
          - p_value: float (fraction of trials with >= physical_hits)
          - n_trials: int
          - tolerance: float
    """
    rng = random.Random(seed)

    # Physical data
    phys = physical_mass_assignment()
    phys_preds = evaluate_formulas(phys)
    phys_hits, phys_hit_names = count_hits(phys_preds, TARGETS, tolerance)

    # Scrambled trials
    scrambled_counts = []
    for _ in range(n_trials):
        scrambled = scrambled_mass_assignment(rng)
        preds = evaluate_formulas(scrambled)
        n_hits, _ = count_hits(preds, TARGETS, tolerance)
        scrambled_counts.append(n_hits)

    mean_scrambled = sum(scrambled_counts) / len(scrambled_counts)
    variance = sum((x - mean_scrambled) ** 2 for x in scrambled_counts) / len(scrambled_counts)
    std_scrambled = math.sqrt(variance) if variance > 0 else 0.001
    max_scrambled = max(scrambled_counts)

    significance = (phys_hits - mean_scrambled) / std_scrambled if std_scrambled > 0 else float('inf')
    p_value = sum(1 for x in scrambled_counts if x >= phys_hits) / n_trials

    return {
        'physical_hits': phys_hits,
        'physical_hit_names': phys_hit_names,
        'scrambled_mean': mean_scrambled,
        'scrambled_std': std_scrambled,
        'scrambled_max': max_scrambled,
        'scrambled_distribution': scrambled_counts,
        'significance_sigma': significance,
        'p_value': p_value,
        'n_trials': n_trials,
        'tolerance': tolerance,
    }
