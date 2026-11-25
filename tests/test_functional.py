from phagree.functional import _union_length, classify

ANN = {
    "phrog_1": {"annot": "integrase", "category": "integration and excision"},
    "phrog_2": {"annot": "terminase", "category": "head and packaging"},
    "phrog_3": {"annot": None, "category": "unknown function"},
}


def hit(phrog, e=1e-20, bits=100.0, q=0.9, p=0.9):
    return {"phrog": phrog, "evalue": e, "bitscore": bits, "query_coverage": q,
            "profile_coverage": p}


def test_union_length_merges_overlaps_inclusive():
    assert _union_length([(1, 10), (5, 20), (30, 30)]) == 21


def test_no_hits_is_no_informative_hit():
    assert classify([], ANN)["function_status"] == "no_informative_hit"


def test_thresholds_are_inclusive_and_all_required():
    assert classify([hit("phrog_1", e=1e-5, q=0.5, p=0.5)], ANN)["function_status"] == (
        "specific_function")
    for bad in (hit("phrog_1", e=2e-5), hit("phrog_1", q=0.49), hit("phrog_1", p=0.49)):
        assert classify([bad], ANN)["function_status"] == "no_informative_hit"


def test_unknown_best_hit_is_homology_unknown():
    row = classify([hit("phrog_3", bits=200), hit("phrog_1", bits=50)], ANN)
    assert row["function_status"] == "homology_unknown_function"
    assert row["best_phrog"] == "phrog_3"


def test_near_tie_with_conflicting_known_categories_is_ambiguous():
    row = classify([hit("phrog_1", bits=100), hit("phrog_2", bits=95)], ANN)
    assert row["function_status"] == "ambiguous"
    row = classify([hit("phrog_1", bits=100), hit("phrog_2", bits=80)], ANN)
    assert row["function_status"] == "specific_function"


def test_unannotated_phrog_id_defaults_to_unknown():
    assert classify([hit("phrog_999")], ANN)["function_status"] == "homology_unknown_function"


def test_extended_keeps_study_functions():
    from phagree.functional import classify_extended

    row = classify_extended([hit("phrog_1", bits=100), hit("phrog_2", bits=300, q=1.0, p=0.2)], ANN)
    assert row["best_phrog"] == "phrog_1" and row["function_rule"] == "study"


def test_extended_rescues_protein_aligned_to_a_longer_profile():
    from phagree.functional import classify_extended

    hits = [hit("phrog_2", bits=300, q=0.99, p=0.46), hit("phrog_3", bits=200, q=0.6, p=0.6)]
    assert classify(hits, ANN)["function_status"] == "homology_unknown_function"
    row = classify_extended(hits, ANN)
    assert row["function_status"] == "specific_function"
    assert row["best_phrog"] == "phrog_2" and row["function_rule"] == "long_profile"


def test_extended_needs_most_of_query_and_part_of_profile_and_never_downgrades():
    from phagree.functional import classify_extended

    partial = [hit("phrog_2", bits=300, q=0.7, p=0.4)]
    assert classify_extended(partial, ANN)["function_status"] == "no_informative_hit"
    tiny_profile_part = [hit("phrog_2", bits=300, q=0.99, p=0.2)]
    assert classify_extended(tiny_profile_part, ANN)["function_status"] == "no_informative_hit"
    assert classify_extended(partial, ANN)["function_status"] == "no_informative_hit"
    # a better unknown-function long-profile hit does not replace the study result
    hits = [hit("phrog_3", bits=300, q=0.95, p=0.2), hit("phrog_2", bits=100, q=0.6, p=0.6)]
    row = classify_extended(hits, ANN)
    assert row["best_phrog"] == "phrog_2" and row["function_rule"] == "study"
