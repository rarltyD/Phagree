from phagree.lifestyle import CANDIDATES, SELECTED, STRONG, evidence, genome_markers, marker_group


def test_marker_groups():
    assert marker_group("integrase") == "integrase"
    assert marker_group("excisionase and transcriptional regulator") == "excisionase"
    assert marker_group("anti-repressor Ant") == "anti-repressor"
    assert marker_group("ParB-like partition nuclease") is None  # common in virulent phages
    assert marker_group("ParB-like partition protein") == "partition protein"
    assert marker_group("terminase large subunit") is None
    assert marker_group(None) is None
    assert set(SELECTED) <= set(CANDIDATES) and set(STRONG) <= set(SELECTED)


def test_evidence_levels():
    assert evidence(["integrase", "excisionase"]) == "strong"
    assert evidence(["transcriptional repressor"]) == "weak"
    assert evidence(["excisionase", "transposase"]) == "none"  # listed, but not selected
    assert evidence([]) == "none"


def test_genome_markers_use_only_genes_with_a_phrogs_function():
    genes = [{"start": 1, "end": 900, "strand": 1, "product": "integrase", "phrog": None},
             {"start": 1000, "end": 1600, "strand": -1, "product": "anti-repressor", "phrog": "phrog_9"}]
    level, found = genome_markers(genes)
    assert level == "weak" and [m["group"] for m in found] == ["anti-repressor"]
