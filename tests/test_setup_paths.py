from pathlib import Path

from phagree.setup_tools import LIBRARY_NAME, db_root, make_config, run_template


def test_db_root_order(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("PHAGREE_DB", raising=False)
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    assert db_root() == tmp_path / "xdg" / "phagree"
    library = tmp_path / "data" / "db" / LIBRARY_NAME
    library.parent.mkdir(parents=True)
    library.write_bytes(b"")
    assert db_root() == Path("data/db")
    monkeypatch.setenv("PHAGREE_DB", str(tmp_path / "env"))
    assert db_root() == tmp_path / "env"
    assert db_root(tmp_path / "given") == tmp_path / "given"


def test_packaged_template_matches_study_settings():
    template = run_template()
    assert template["partial_overlap_reciprocal_primary"] == 0.8
    assert template["phanotate_command"] is None
    config = make_config(["/env/bin/python", "/env/bin/phanotate.py"])
    assert config["phanotate_command"] == ["/env/bin/python", "/env/bin/phanotate.py"]
    assert config["phanotate_version_command"][0] == "/env/bin/python"
