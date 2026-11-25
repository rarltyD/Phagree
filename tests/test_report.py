import pandas as pd
from Bio.Seq import Seq

from phagree.report import coding_fraction, seq_id, translate


def test_translate_plus_minus_and_alternative_start():
    seq = Seq("GTGAAATTTTAA" + "CCC" + str(Seq("ATGGGGCCCTAG").reverse_complement()))
    assert translate(seq, 1, 12, 1) == "MKF"  # GTG start is read as M
    assert translate(seq, 16, 27, -1) == "MGP"


def test_translate_with_ambiguous_base():
    assert translate(Seq("ATGNAATTTTAA"), 1, 12, 1) == "MXF"


def test_seq_id_and_coding_fraction():
    assert seq_id("NC_001416.1") == "NC_001416_1"
    genes = pd.DataFrame({"start1": [1, 5], "end1": [10, 12]})
    assert coding_fraction(genes, 20) == 12 / 20
