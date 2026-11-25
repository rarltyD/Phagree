"""Lysogeny markers: PHROGs functions typical of temperate phages, per genome.

A genome is reported with
* ``strong`` evidence when it encodes an integrase,
* ``weak`` evidence when it encodes no integrase but a transcriptional repressor or an anti-repressor,
* ``none`` otherwise ("no markers found"; this does not mean the phage is virulent).

The marker groups were selected on one half of the relatedness blocks of the study cohort and evaluated on
the other half against the INPHARED2 lifestyle labels (PhageLeads predictions), see docs/VALIDATION_1.1.md:
integrase alone (strong) and any of the three groups (strong or weak) give different trade-offs between
sensitivity and specificity.  All candidate groups are listed per gene for transparency, but only the
selected groups enter the evidence level.
"""

from __future__ import annotations

import re

# Candidate groups (regular expressions on the PHROGs annotation of a protein with a specific function).
CANDIDATES = {
    "integrase": r"^integrase$",
    "excisionase": r"excisionase|recombination directionality factor",
    "transposase": r"^transposase$|^DNA transposition protein$",
    "CI-like repressor": r"^CI-like repressor$",
    "Cro-like repressor": r"^cro-like repressor$",
    "superinfection immunity": r"^immunity to superinfection$",
    "partition protein": r"partition protein(?! nuclease)|ParA-like|^plasmid stability$",
    "Arc-like repressor": r"^Arc-like repressor$",
    "anti-repressor": r"^anti-repressor",
    "recombinase": r"^recombinase$",
    "transcriptional repressor": r"^transcriptional repressor$",
}
STRONG = ("integrase",)
SELECTED = ("integrase", "transcriptional repressor", "anti-repressor")
_RX = {name: re.compile(pattern, re.I) for name, pattern in CANDIDATES.items()}


def marker_group(annotation):
    """Name of the candidate group of a PHROGs annotation, or None."""
    if not isinstance(annotation, str) or not annotation:
        return None
    for name, rx in _RX.items():
        if rx.search(annotation):
            return name
    return None


def evidence(groups):
    """Evidence level from the set of marker groups found in a genome."""
    groups = set(groups)
    if groups & set(STRONG):
        return "strong"
    if groups & set(SELECTED):
        return "weak"
    return "none"


def genome_markers(genes):
    """genes: iterable of dicts with keys start, end, strand, product, phrog (consensus genes of one genome).

    Returns (evidence level, list of marker genes with their group and whether the group is selected).
    """
    found = []
    for g in genes:
        if not g.get("phrog"):
            continue
        group = marker_group(g.get("product"))
        if group:
            found.append({**g, "group": group, "selected": group in SELECTED})
    return evidence(m["group"] for m in found), found
