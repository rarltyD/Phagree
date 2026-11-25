from Bio import SeqIO

from phagree.annotate import collect_inputs


def test_fasta_genomes_become_versioned_genbank_records(tmp_path):
    fasta = tmp_path / "new.fasta"
    fasta.write_text(">my phage/1 test\nATGAAATAA\n>second\nATGCCCTAA\n")
    genomes = collect_inputs([str(fasta)], tmp_path / "gb")
    assert [g for g, _ in genomes] == ["my.1", "second.1"]
    record = SeqIO.read(genomes[0][1], "genbank")
    assert record.id == "my.1" and str(record.seq) == "ATGAAATAA"
    assert record.annotations["topology"] == "linear"


def test_unknown_input_is_rejected(tmp_path):
    import pytest

    with pytest.raises(ValueError):
        collect_inputs(["not a genome"], tmp_path / "gb")


def test_input_kinds_are_reported(tmp_path):
    fasta = tmp_path / "new.fasta"
    fasta.write_text(">new_phage\nATGAAATTTTAA\n")
    kinds = {}
    genomes = collect_inputs([str(fasta)], tmp_path / "gb", kinds=kinds)
    assert kinds == {genomes[0][0]: "fasta"}
