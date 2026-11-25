"""PHROGs search space: hits of a protein must not depend on what else is searched with it."""

import json
import random

import pytest

pyhmmer = pytest.importorskip("pyhmmer")

from phagree.functional import search_chunk  # noqa: E402
from phagree.phrogs_search import E_REPORT, cache_init  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def _mutate(seq, rng, rate):
    return "".join(rng.choice(AA) if rng.random() < rate else c for c in seq)


@pytest.fixture(scope="module")
def tiny_db(tmp_path_factory):
    """Three profiles built from random 'families'; queries are mutated members and random decoys."""
    rng = random.Random(1)
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    path = tmp_path_factory.mktemp("hmm") / "tiny.hmm"
    families = []
    with open(path, "wb") as out:
        for i in range(3):
            seed = "".join(rng.choice(AA) for _ in range(120))
            families.append(seed)
            seq = pyhmmer.easel.TextSequence(name=f"fam{i}".encode(), sequence=seed).digitize(alphabet)
            hmm, _, _ = builder.build(seq, background)
            hmm.write(out)
    members = [(f"m{i}_{j}", _mutate(f, rng, 0.35)) for i, f in enumerate(families) for j in range(4)]
    decoys = [(f"d{k}", "".join(rng.choice(AA) for _ in range(rng.randint(60, 300)))) for k in range(300)]
    return path, members, decoys


def _key(hits):
    return sorted((h["protein_sha256"], h["phrog"], h["evalue"], h["best_domain_ievalue"], h["query_coverage"],
                   h["profile_coverage"]) for h in hits)


def test_hits_do_not_depend_on_batch(tiny_db):
    path, members, decoys = tiny_db
    alone = search_chunk(members, str(path), 1, E_REPORT)
    batch = [h for h in search_chunk(members + decoys, str(path), 1, E_REPORT)
             if h["protein_sha256"].startswith("m")]
    assert alone, "the test families must produce hits"
    assert _key(alone) == _key(batch)


def test_auto_search_space_depends_on_batch(tiny_db):
    """The 1.0.0 behaviour (z/domz from the batch) gives different E-values; kept as an explicit option."""
    path, members, decoys = tiny_db
    alone = {(h["protein_sha256"], h["phrog"]): h["evalue"]
             for h in search_chunk(members, str(path), 1, E_REPORT, z=None, domz=None)}
    batch = {(h["protein_sha256"], h["phrog"]): h["evalue"]
             for h in search_chunk(members + decoys, str(path), 1, E_REPORT, z=None, domz=None)}
    common = set(alone) & set(batch)
    assert common and any(alone[k] != batch[k] for k in common)


def test_cache_records_and_checks_search_space(tmp_path, tiny_db):
    path = tiny_db[0]
    space = {"z": 20000, "domz": 13}
    meta = cache_init(tmp_path / "c1", path, space)
    assert meta["search_space"] == space
    assert json.loads((tmp_path / "c1" / "cache_meta.json").read_text())["search_space"] == space
    assert cache_init(tmp_path / "c1", path, space)["search_space"] == space
    with pytest.raises(SystemExit):
        cache_init(tmp_path / "c1", path, {"z": None, "domz": None})


def test_legacy_cache_without_search_space_is_accepted(tmp_path, tiny_db):
    path = tiny_db[0]
    cache_init(tmp_path / "c2", path, {"z": 20000, "domz": 13})
    meta_path = tmp_path / "c2" / "cache_meta.json"
    legacy = json.loads(meta_path.read_text())
    del legacy["search_space"]
    meta_path.write_text(json.dumps(legacy))
    meta = cache_init(tmp_path / "c2", path, {"z": 20000, "domz": 13})
    assert "not recorded" in str(meta["search_space"])
