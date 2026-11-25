"""One-time local setup: PHROGs HMM library and a run configuration with PHANOTATE.

The PHROGs library is not redistributed.  The official PHROGs v4 alignment archive is downloaded (the
PHROGs site first, a Zenodo mirror if the site is unreachable), checked against its published SHA-256,
and rebuilt as HMMER3 profiles in parallel.  The built library is verified by a content checksum over
all profile parameters, so it is known to be identical to the library of the study even though the
file checksum differs (each profile stores its build time).
"""

from __future__ import annotations

import glob
import hashlib
import io
import json
import os
import re
import shutil
import sys
import tarfile
import time
import urllib.request
from multiprocessing import Pool
from pathlib import Path

PHROGS_FILES = {
    "MSA_phrogs.tar.gz": "1ca216aa830ebb94c925d2733c14d5870d3922bbe7a6eddffc66122cfc945923",
    "phrog_annot_v4.tsv": "502f96101597c21133bcce5711803e0b95e0c162cd4e86425c352549bd95e8c2",
}
PHROGS_SOURCES = (
    "https://phrogs.lmge.uca.fr/downloads_from_website/{name}",
    "https://zenodo.org/records/23342382/files/{name}?download=1",
)
LIBRARY_CONTENT_SHA256 = "02bff2628221973209641c29ff32f1791039121152088e2a490909af190d9981"
LIBRARY_NAME = "phrogs_v4/phrogs_v4_m50.h3m"


def db_root(db=None) -> Path:
    """Folder with phrogs_v4/, the PHROGs search cache and run_config.json.

    Order: --db, $PHAGREE_DB, data/db of a repository checkout (the study layout) when it holds
    a library, otherwise the user data folder (~/.local/share/phagree).
    """
    if db:
        return Path(db).expanduser()
    if os.environ.get("PHAGREE_DB"):
        return Path(os.environ["PHAGREE_DB"]).expanduser()
    if (Path("data/db") / LIBRARY_NAME).exists():
        return Path("data/db")
    base = os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share"
    return Path(base) / "phagree"


def run_template() -> dict:
    from importlib.resources import files

    return json.loads(files("phagree").joinpath("templates/run_template.json").read_text("utf-8"))


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download_phrogs(dest: Path, names=None, attempts=3) -> dict:
    """Download and verify PHROGs v4 inputs (alignment archive, annotation table)."""
    dest.mkdir(parents=True, exist_ok=True)
    report = {}
    for name in names or PHROGS_FILES:
        expected = PHROGS_FILES[name]
        path = dest / name
        if path.exists() and sha256_file(path) == expected:
            report[name] = "present"
            continue
        errors = []
        for attempt in range(attempts):  # transient errors (e.g. HTTP 504) are retried
            for template in PHROGS_SOURCES:
                url = template.format(name=name)
                try:
                    tmp = path.with_suffix(path.suffix + ".part")
                    with urllib.request.urlopen(url, timeout=300) as response, tmp.open("wb") as out:
                        shutil.copyfileobj(response, out)
                    if sha256_file(tmp) != expected:
                        raise ValueError("checksum mismatch")
                    tmp.replace(path)
                    report[name] = url
                    break
                except Exception as error:  # try the next source
                    errors.append(f"{url}: {type(error).__name__}: {error}")
            else:
                if attempt + 1 < attempts:
                    time.sleep(30 * (attempt + 1))
                continue
            break
        else:
            raise RuntimeError(f"could not obtain {name}:\n" + "\n".join(errors))
    return report


_BUILDER = None


def _init_builder():
    global _BUILDER
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    _BUILDER = (alphabet, pyhmmer.plan7.Builder(alphabet), pyhmmer.plan7.Background(alphabet))


def _build_one(item):
    from pyhmmer.easel import MSAFile

    name, data = item
    alphabet, builder, background = _BUILDER
    try:
        with MSAFile(io.BytesIO(data), format="afa", digital=True, alphabet=alphabet) as f:
            msa = f.read()
        msa.name = name.encode()
        hmm, _, _ = builder.build_msa(msa, background)
        buffer = io.BytesIO()
        hmm.write(buffer, binary=True)
        return name, buffer.getvalue(), None
    except Exception as error:
        return name, None, f"{type(error).__name__}: {error}"


def _alignments(msa_archive):
    with tarfile.open(msa_archive, "r:gz") as archive:
        for member in archive:
            if member.isfile() and member.name.endswith(".fma"):
                yield re.sub(r"\.fma$", "", Path(member.name).name), archive.extractfile(member).read()


def library_content_sha256(path) -> tuple[int, str]:
    """Checksum over profile names and parameters, independent of build time and command line."""
    import numpy as np
    import pyhmmer

    h, n = hashlib.sha256(), 0
    with pyhmmer.plan7.HMMFile(str(path)) as hmms:
        for hmm in hmms:
            name = hmm.name if isinstance(hmm.name, bytes) else hmm.name.encode()
            h.update(name)
            for values in (hmm.match_emissions, hmm.insert_emissions, hmm.transition_probabilities):
                h.update(np.ascontiguousarray(np.asarray(values), dtype=np.float32).tobytes())
            n += 1
    return n, h.hexdigest()


def build_library(msa_archive: Path, out: Path, workers: int) -> dict:
    """Build the HMMER3 library in parallel, in archive order, with pyhmmer Builder defaults."""
    started = time.perf_counter()
    tmp = out.with_name(out.name + ".tmp")
    failed = []
    with Pool(workers, initializer=_init_builder) as pool, tmp.open("wb") as handle:
        for name, blob, error in pool.imap(_build_one, _alignments(msa_archive), chunksize=64):
            if error:
                failed.append({"phrog": name, "error": error})
            else:
                handle.write(blob)
    tmp.replace(out)
    profiles, content = library_content_sha256(out)
    report = {"output": str(out), "profiles_built": profiles, "failed": failed, "workers": workers,
              "content_sha256": content, "identical_to_study_library": content == LIBRARY_CONTENT_SHA256,
              "wall_seconds": round(time.perf_counter() - started, 1)}
    out.with_suffix(".build.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def find_phanotate() -> list[str] | None:
    """[python, phanotate.py] from this environment, PATH or a conda environment."""
    candidates = [Path(sys.executable).parent / "phanotate.py"]
    found = shutil.which("phanotate.py")
    if found:
        candidates.append(Path(found))
    home = Path.home()
    for base in (home / "miniforge3", home / "mambaforge", home / "miniconda3", home / "anaconda3",
                 Path("/opt/conda")):
        candidates += [Path(p) for p in sorted(glob.glob(str(base / "envs" / "*" / "bin" / "phanotate.py")))]
    for script in candidates:
        if script.exists():
            python = script.parent / "python"
            return [str(python if python.exists() else Path(sys.executable)), str(script)]
    return None


def make_config(phanotate: list[str]) -> dict:
    config = run_template()
    python = phanotate[0]
    config["phanotate_command"] = phanotate
    config["phanotate_version_command"] = [python, "-c",
                                           "import importlib.metadata as m;print(m.version('phanotate'))"]
    config["phanotate_environment_command"] = [
        python, "-c", "import json,platform,importlib.metadata as m;print(json.dumps({'platform':"
        "platform.platform(),'python':platform.python_version(),'phanotate':m.version('phanotate')}))"]
    config["max_workers"] = max(1, min(os.cpu_count() or 1, config.get("max_workers", 4)))
    config["run_profile"] = "local_setup"
    return config


def setup(db=None, workers=None, rebuild=False):
    root = db_root(db)
    phrogs = root / "phrogs_v4"
    workers = workers or max(1, (os.cpu_count() or 2) - 1)
    report = {"database": str(root.resolve()),
              "phrogs_downloads": download_phrogs(phrogs, ["phrog_annot_v4.tsv"])}
    library = root / LIBRARY_NAME
    reused = None
    if library.exists() and not rebuild:
        profiles, content = library_content_sha256(library)
        reused = {"output": str(library), "profiles_built": profiles, "content_sha256": content,
                  "identical_to_study_library": content == LIBRARY_CONTENT_SHA256, "reused": True}
    if reused and reused["identical_to_study_library"]:
        report["library"] = reused  # the alignment archive is needed only to build
    else:
        report["phrogs_downloads"].update(download_phrogs(phrogs, ["MSA_phrogs.tar.gz"]))
        report["library"] = build_library(phrogs / "MSA_phrogs.tar.gz", library, workers)
    phanotate = find_phanotate()
    report["phanotate_command"] = phanotate
    if phanotate:
        config = root / "run_config.json"
        config.write_text(json.dumps(make_config(phanotate), indent=2) + "\n", encoding="utf-8")
        report["config"] = str(config)
    report["ready"] = bool(phanotate) and report["library"]["identical_to_study_library"]
    if not phanotate:
        report["next_step"] = ("PHANOTATE not found: install it in this environment "
                               "(conda install -c bioconda phanotate=1.6.7), then run setup again")
    return report
