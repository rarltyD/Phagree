"""PHROGs HMM search over the deduplicated protein catalogue (protocol section 5).

One search per exact protein sequence (``protein_sha256``); results are joined
back to models/loci via ``protein_origins``.  This module holds the search of one
chunk and the classification rules; ``phagree.phrogs_search`` runs the
chunks through a study-wide cache of immutable Parquet files, so an interrupted
search resumes without recomputing completed chunks.

Status rules (prespecified in ``configs/protocol_v1.json``; thresholds are
investigator choices, not validated biological truth):

* a hit *passes* when full-sequence E <= e_max, query coverage >= qcov_min and
  profile coverage >= pcov_min (coverages are unions of domain envelopes /
  HMM intervals of included domains);
* ``specific_function``: best passing hit belongs to a PHROG whose category is
  not ``unknown function``;
* ``homology_unknown_function``: passing hits exist, best is an unknown-function
  PHROG;
* ``ambiguous``: another passing hit with bitscore >= 0.9 x best bitscore has a
  different *known* category than the best hit;
* ``no_informative_hit``: no passing hit;
* ``not_searched_or_failed`` is assigned by callers for proteins never searched.
"""

from __future__ import annotations

import csv

import pyarrow as pa

UNKNOWN = "unknown function"
NEAR_TIE = 0.9
# Fixed search space, so that a protein gets the same hits whatever else is searched with it.  By default
# HMMER sets Z to the number of proteins searched together and domZ to the number of significant proteins
# per profile, which made borderline E-values (and functions) depend on how many genomes were processed in
# one run.  Z = 20,000 is the size of a search block of the study (E-values of full blocks are reproduced
# exactly); domZ = 13 is the hit-weighted median number of significant proteins per profile in those blocks.
SEARCH_SPACE_Z = 20000
SEARCH_SPACE_DOMZ = 13

HIT_SCHEMA = pa.schema(
    [
        ("protein_sha256", pa.string()),
        ("phrog", pa.string()),
        ("evalue", pa.float64()),
        ("bitscore", pa.float64()),
        ("best_domain_ievalue", pa.float64()),
        ("query_coverage", pa.float64()),
        ("profile_coverage", pa.float64()),
        ("n_domains_included", pa.int64()),
        ("length_aa", pa.int64()),
        ("profile_length", pa.int64()),
    ]
)


def _union_length(intervals):
    total, current = 0, None
    for start, end in sorted(intervals):
        if current is None or start > current[1]:
            if current:
                total += current[1] - current[0] + 1
            current = [start, end]
        else:
            current[1] = max(current[1], end)
    if current:
        total += current[1] - current[0] + 1
    return total


def load_annotations(path):
    table = {}
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            table["phrog_" + row["phrog"]] = {
                "annot": row["annot"] if row["annot"] not in ("", "NA") else None,
                "category": row["category"] or UNKNOWN,
            }
    return table


def search_chunk(sequences, hmm_db, cpus, e_report, z=SEARCH_SPACE_Z, domz=SEARCH_SPACE_DOMZ):
    """sequences: list of (sha, aa). Returns list of hit dicts (E <= e_report).

    ``z`` and ``domz`` fix the search space (see SEARCH_SPACE_Z); ``None`` lets HMMER set them from the
    proteins searched together (the behaviour of Phagree 1.0.0).
    """
    import pyhmmer
    from pyhmmer.easel import Alphabet, TextSequence
    from pyhmmer.plan7 import HMMFile

    alphabet = Alphabet.amino()
    block = pyhmmer.easel.DigitalSequenceBlock(
        alphabet,
        [
            TextSequence(name=sha.encode(), sequence=aa).digitize(alphabet)
            for sha, aa in sequences
        ],
    )
    lengths = {sha: len(aa) for sha, aa in sequences}
    hits = []
    with HMMFile(hmm_db) as handle:
        profiles = handle.optimized_profiles() if handle.is_pressed() else handle
        space = {k: v for k, v in (("Z", z), ("domZ", domz)) if v is not None}
        for top in pyhmmer.hmmsearch(profiles, block, cpus=cpus, E=e_report, domE=e_report, **space):
            phrog = top.query.name.decode() if isinstance(top.query.name, bytes) else top.query.name
            model_len = top.query.M
            for hit in top:
                if not hit.reported:
                    continue
                name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                domains = [d for d in hit.domains if d.included] or list(hit.domains)
                q_iv = [(d.env_from, d.env_to) for d in domains]
                p_iv = [(d.alignment.hmm_from, d.alignment.hmm_to) for d in domains]
                hits.append(
                    {
                        "protein_sha256": name,
                        "phrog": phrog,
                        "evalue": float(hit.evalue),
                        "bitscore": float(hit.score),
                        "best_domain_ievalue": float(min(d.i_evalue for d in domains)),
                        "query_coverage": _union_length(q_iv) / lengths[name],
                        "profile_coverage": _union_length(p_iv) / model_len,
                        "n_domains_included": sum(d.included for d in hit.domains),
                        "length_aa": lengths[name],
                        "profile_length": model_len,
                    }
                )
    return hits


def classify(hits, annotations, e_max=1e-5, qcov_min=0.5, pcov_min=0.5):
    """Return one status row for a single protein given all its raw hits."""
    passing = sorted(
        (
            h
            for h in hits
            if h["evalue"] <= e_max
            and h["query_coverage"] >= qcov_min
            and h["profile_coverage"] >= pcov_min
        ),
        key=lambda h: (-h["bitscore"], h["evalue"], h["phrog"]),
    )
    if not passing:
        return {"function_status": "no_informative_hit", "best_phrog": None,
                "best_category": None, "best_annot": None, "best_evalue": None,
                "best_bitscore": None, "n_passing_hits": 0}
    best = passing[0]
    info = annotations.get(best["phrog"], {"annot": None, "category": UNKNOWN})
    status = "homology_unknown_function" if info["category"] == UNKNOWN else "specific_function"
    near = {
        annotations.get(h["phrog"], {}).get("category", UNKNOWN)
        for h in passing
        if h["bitscore"] >= NEAR_TIE * best["bitscore"]
    } - {UNKNOWN}
    if status == "specific_function" and len(near) > 1:
        status = "ambiguous"
    return {
        "function_status": status,
        "best_phrog": best["phrog"],
        "best_category": info["category"],
        "best_annot": info["annot"],
        "best_evalue": best["evalue"],
        "best_bitscore": best["bitscore"],
        "n_passing_hits": len(passing),
    }
