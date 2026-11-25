import hashlib
import json
from pathlib import Path

import pytest

from phagree.translation import translate_cds

SEQUENCES = json.loads(
    (Path(__file__).parent / "fixtures" / "synthetic" / "sequences.json").read_text()
)


@pytest.mark.parametrize(
    "key,codon_start", [("complete", 1), ("codon_start_2", 2), ("codon_start_3", 3)]
)
def test_offsets_and_terminal_stop_rule(key, codon_start):
    result = translate_cds(SEQUENCES[key], codon_start=codon_start, provided_translation="MK*")
    assert result["recomputed_translation_raw"] == "MK*"
    assert result["recomputed_translation"] == "MK"
    assert result["provided_translation"] == "MK*"
    assert result["translated_nt_length"] == 9
    assert result["translation_validity"] == "valid"
    assert result["protein_sha256"] == hashlib.sha256(b"MK").hexdigest()


def test_genbank_default_one_and_table_eleven_initiation():
    one = translate_cds(SEQUENCES["gtg"], transl_table=1)
    eleven = translate_cds(SEQUENCES["gtg"], transl_table=11)
    assert one["recomputed_translation"] == "VK"
    assert one["translation_validity"] == "invalid_start"
    assert eleven["recomputed_translation"] == "MK"
    assert eleven["translation_notes"] == ["alternative_initiation_met"]


def test_code_fifteen_tag_is_glutamine():
    alt = translate_cds(SEQUENCES["table15"], transl_table=15)
    assert alt["recomputed_translation"] == "MQ"
    assert alt["translation_validity"] == "valid"
    standard = translate_cds(SEQUENCES["table15"], transl_table=1)
    assert standard["recomputed_translation"] == "M*"
    assert standard["translation_validity"] == "internal_stop"
    assert standard["protein_sha256"] is None


def test_provided_disagreement_does_not_destroy_valid_recomputation():
    result = translate_cds(SEQUENCES["complete"], provided_translation="MM")
    assert result["provided_translation"] == "MM"
    assert result["recomputed_translation"] == "MK"
    assert result["translation_validity"] == "provided_translation_mismatch"
    assert result["protein_validation"] == "valid"
    assert result["protein_sha256"]


@pytest.mark.parametrize(
    "key,kwargs,reason",
    [
        ("invalid_symbol", {}, "invalid_nucleotide"),
        ("ambiguous", {}, "ambiguous_translation"),
        ("internal_stop", {}, "internal_stop"),
        ("complete", {"pseudo": True}, "pseudo"),
        ("complete", {"transl_except": ["synthetic exception"]}, "transl_except"),
        ("complete", {"exception": "ribosomal slippage"}, "exception"),
        ("partial", {"partial_5p": True, "partial_3p": True}, "partial"),
        ("complete", {"codon_start": 0}, "invalid_codon_start"),
        ("complete", {"transl_table": 9999}, "unsupported_table"),
    ],
)
def test_exceptional_sequences_retained_without_validated_protein(key, kwargs, reason):
    result = translate_cds(SEQUENCES[key], **kwargs)
    assert result["translation_validity"] == reason
    assert result["protein_sha256"] is None


def test_ambiguous_codons_beyond_x_are_flagged():
    # RAY (A/G, A, C/T) translates to B (D or N): the CDS must not count as valid.
    seq = "ATG" + "RAY" + "AAA" + "TAA"
    result = translate_cds(seq)
    assert "B" in result["recomputed_translation_raw"]
    assert result["translation_validity"] == "ambiguous_translation"
    assert result["protein_validation"] == "invalid"
    assert result["protein_sha256"] is None


def test_ambiguous_protein_is_masked_only_when_enabled():
    from phagree import translation

    seq = "ATG" + "RAY" + "NNN" + "AAA" + "TAA"
    translation.set_ambiguous_proteins(True)
    try:
        result = translate_cds(seq)
    finally:
        translation.set_ambiguous_proteins(False)
    assert result["protein_validation"] == "ambiguous"
    assert result["recomputed_translation"] == "MXXK"
    assert result["protein_sha256"] == hashlib.sha256(b"MXXK").hexdigest()
    assert result["translation_validity"] == "ambiguous_translation"


def test_ambiguous_protein_enters_catalogue():
    from phagree.proteins import build_catalogue

    call = {"call_id": "c1", "genome_id": "g", "source_run_id": "r", "protein_validation": "ambiguous",
            "translation_validity": "ambiguous_translation", "recomputed_translation": "MXXK",
            "protein_sha256": hashlib.sha256(b"MXXK").hexdigest()}
    cat = build_catalogue([call])
    assert cat["proteins"][0]["validation"] == "ambiguous"
    assert cat["proteins"][0]["amino_acid_sequence"] == "MXXK"
