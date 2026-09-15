"""
Program: CHIPS_core.py
Author: Alexandre Prefontaine
Affiliation: University of Manitoba, Manitoba Centre for Proteomics and Systems Biology
Date of Last Update: 2024-06-20
Key Contributors: Oleg Krokhin, Victor Spicer, Alexandre Prefontaine

Core calculation module for CHIPs (Chromatographic Hydrophobicity Index for Peptides).

This module contains only the constants and calculation function.
It has no dependencies beyond the Python standard library and can be
imported directly into any script:

    from CHIPS_core import calculate_HI

The working equation is:

    HI = (1 - a * ln(n)) * (SUM(ci) + b)

    where:
        n    = peptide length
        a    = length correction factor (phase-dependent)
        ci   = retention coefficient for each amino acid (phase-dependent)
        b    = free term / intercept (phase-dependent)

All cysteines are assumed to be alkylated (carbamidomethylated).
PTMs are not supported — pass stripped sequences only.
"""

import math

# --- Valid amino acid codes ---
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')

# --- Common header names to exclude from peptide lists ---
KNOWN_HEADERS = {'seq', 'sequence', 'peptide', 'peptides', 'pep', 'protein'}

# --- Formic Acid (FA) constants ---
ALPHA_FA = 0.207
BETA_FA  = -0.766

RC_FA = {
    'A': 2.202,
    'C': 0.908,  # alkylated (carbamidomethyl)
    'D': 2.001,
    'E': 2.136,
    'F': 9.436,
    'G': 0.484,
    'H': -4.346,
    'I': 7.292,
    'K': -5.979,
    'L': 8.223,
    'M': 6.065,
    'N': 0.451,
    'P': 1.680,
    'Q': 0.648,
    'R': -5.614,
    'S': 0.853,
    'T': 1.486,
    'V': 4.704,
    'W': 10.394,
    'Y': 4.831,
}

# --- Trifluoroacetic Acid (TFA) constants ---
ALPHA_TFA = 0.213
BETA_TFA  = -0.766

RC_TFA = {
    'A': 2.293,
    'C': 1.111,  # alkylated (carbamidomethyl)
    'D': 1.343,
    'E': 1.799,
    'F': 9.959,
    'G': 0.739,
    'H': -1.154,
    'I': 7.555,
    'K': -1.140,
    'L': 8.591,
    'M': 6.370,
    'N': 0.516,
    'P': 1.543,
    'Q': 0.633,
    'R': -0.265,
    'S': 0.977,
    'T': 1.652,
    'V': 4.980,
    'W': 11.032,
    'Y': 5.119,
}


def calculate_HI(sequence: str, phase: str) -> float:
    """
    Calculate the Chromatographic Hydrophobicity Index (HI) of a peptide.

    Parameters
    ----------
    sequence : str
        Single-letter amino acid sequence (stripped, no PTMs).
        All cysteines assumed to be alkylated (carbamidomethyl).
    phase : str
        Ion pairing agent. Must be 'FA' or 'TFA'.

    Returns
    -------
    float
        Predicted hydrophobicity index in HI units, rounded to 5 decimal places.

    Raises
    ------
    ValueError
        If phase is not 'FA' or 'TFA'.
    ValueError
        If sequence contains characters outside the 20 standard amino acid codes.

    Examples
    --------
    >>> calculate_HI("FVTDIDELGK", "FA")
    >>> calculate_HI("VTFLGLQHWVPELAR", "TFA")
    """
    # Normalize input
    sequence = sequence.upper().strip()
    phase    = phase.upper().strip()

    # Validate phase
    if phase == 'FA':
        rc, alpha, beta = RC_FA, ALPHA_FA, BETA_FA
    elif phase == 'TFA':
        rc, alpha, beta = RC_TFA, ALPHA_TFA, BETA_TFA
    else:
        raise ValueError(f"Unknown phase '{phase}'. Must be 'FA' or 'TFA'.")

    # Validate sequence
    invalid = [aa for aa in sequence if aa not in VALID_AA]
    if invalid:
        raise ValueError(
            f"Invalid amino acid codes in sequence '{sequence}': {set(invalid)}. "
            f"Only standard single-letter codes are accepted."
        )

    # Calculate HI
    n      = len(sequence)
    sum_ci = sum(rc[aa] for aa in sequence)
    hi     = (1 - alpha * math.log(n)) * (sum_ci + beta)

    return round(hi, 5)


