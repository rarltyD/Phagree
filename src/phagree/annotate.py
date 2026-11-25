"""Annotate or re-annotate any set of phage genomes with the consensus rules of the study.

Inputs can be mixed: GenBank accessions (downloaded from NCBI), GenBank files (checked against their
annotation) and FASTA files (new genomes; each record is one genome).  Steps: gene calling with
PHANOTATE and Pyrodigal-gv, PHROGs search, consensus re-annotation.  CDS that contain ambiguous bases
(N or other IUPAC codes) are kept: compared by coordinates, searched as X-masked proteins and flagged.
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

from Bio import SeqIO

from .setup_tools import LIBRARY_NAME, db_root, find_phanotate, make_config

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
GENBANK_SUFFIXES = (".gb", ".gbk", ".gbff", ".genbank")
FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas")
# PHROGs hits are cached per search space; Phagree 1.0.0 used "phrogs_cache_v4" with a variable search
# space, so those hits are not reused.
CACHE_NAME = "phrogs_cache_v4_z20000_domz13"


def _fetch(accessions, dest):
    data = urllib.parse.urlencode({"db": "nuccore", "id": ",".join(accessions), "rettype": "gbwithparts",
                                   "retmode": "text", "tool": "Phagree"}).encode()
    for attempt in range(5):
        try:
            text = urllib.request.urlopen(EFETCH, data=data, timeout=300).read().decode()
            break
        except Exception:  # transient NCBI errors
            time.sleep(10 * (attempt + 1))
    else:
        raise RuntimeError("NCBI download failed")
    genomes = []
    for record in re.split(r"^//\n", text, flags=re.M):
        m = re.search(r"^VERSION\s+(\S+)", record, re.M)
        if m:
            path = dest / (m.group(1).replace(".", "_") + ".gbk")
            path.write_text(record.lstrip("\n") + "//\n")
            genomes.append((m.group(1), path))
    missing = set(accessions) - {g for g, _ in genomes}
    if missing:
        raise RuntimeError(f"not found at NCBI: {', '.join(sorted(missing))}")
    return genomes


def _genbank(path, dest):
    genomes = []
    for record in SeqIO.parse(path, "genbank"):
        out = dest / (record.id.replace(".", "_") + ".gbk")
        SeqIO.write(record, out, "genbank")
        genomes.append((record.id, out))
    return genomes


def _fasta(path, dest, topology):
    genomes = []
    for record in SeqIO.parse(path, "fasta"):
        name = re.sub(r"[^A-Za-z0-9_]", "_", record.id)[:16] or "genome"
        record.id, record.name = name + ".1", name
        record.annotations = {"molecule_type": "DNA", "topology": topology, "data_file_division": "PHG"}
        out = dest / (name + ".gbk")
        SeqIO.write(record, out, "genbank")
        genomes.append((record.id, out))
    return genomes


def collect_inputs(inputs, dest, topology="linear"):
    dest.mkdir(parents=True, exist_ok=True)
    genomes, accessions = [], []
    for item in inputs:
        path = Path(item)
        if path.suffix.lower() in GENBANK_SUFFIXES and path.exists():
            genomes += _genbank(path, dest)
        elif path.suffix.lower() in FASTA_SUFFIXES and path.exists():
            genomes += _fasta(path, dest, topology)
        elif re.fullmatch(r"[A-Z]{1,6}_?\d+(\.\d+)?", item):
            accessions.append(item)
        else:
            raise ValueError(f"cannot use input: {item}")
    if accessions:
        genomes += _fetch(accessions, dest)
    if len({g for g, _ in genomes}) != len(genomes):
        raise ValueError("duplicate genome identifiers in the input")
    return genomes


def annotate(inputs, out, db=None, config=None, workers=3, topology="linear"):
    out = Path(out)
    root = db_root(db)
    hmm = root / LIBRARY_NAME
    if not hmm.exists():
        raise FileNotFoundError(f"PHROGs library not found in {root}: run 'phagree setup' first")
    if config:
        settings = json.loads(Path(config).read_text(encoding="utf-8"))
    elif (root / "run_config.json").exists():
        settings = json.loads((root / "run_config.json").read_text(encoding="utf-8"))
    else:
        phanotate = find_phanotate()
        if not phanotate:
            raise FileNotFoundError("PHANOTATE not found: run 'phagree setup' for instructions")
        settings = make_config(phanotate)
    genomes = collect_inputs(inputs, out / "genbank", topology)
    manifest = out / "manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "accession_version", "genbank_path",
                                                    "baltimore_group", "host_domain"])
        writer.writeheader()
        for accession, path in genomes:
            writer.writerow({"accession": accession.rsplit(".", 1)[0], "accession_version": accession,
                             "genbank_path": str(path.resolve()), "baltimore_group": "Group I",
                             "host_domain": "Bacteria"})
    settings["search_ambiguous_proteins"] = True
    run_config = out / "run_config_annotate.json"
    run_config.write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
    log = out / "annotate.log"
    python = sys.executable
    steps = [
        [python, "-m", "phagree", "run", "--manifest", manifest, "--config", run_config,
         "--workers", workers, "--out", out / "run", "--raw-root", out / "callers", "--diagnostic"],
        [python, "-m", "phagree.phrogs_search", "--run", out / "run", "--hmm", hmm,
         "--cache", root / CACHE_NAME, "--annotations", root / "phrogs_v4" / "phrog_annot_v4.tsv",
         "--cpus", workers],
        [python, "-m", "phagree.reannotation", out / "run", out / "reannotation",
         "--locus-eligibility", "geometry", "--include-ambiguous"],
    ]
    for i, step in enumerate(steps):
        with log.open("a") as handle:
            rc = subprocess.call([str(x) for x in step], stdout=handle, stderr=subprocess.STDOUT)
        if i == 0 and not (out / "run" / "run_summary.json").exists():
            raise RuntimeError(f"gene calling failed; see {log}")
        if i > 0 and rc:
            raise RuntimeError(f"step {i + 1} failed; see {log}")
    summary = json.loads((out / "run" / "run_summary.json").read_text())
    result = json.loads((out / "reannotation" / "reannotation_summary.json").read_text())
    return {"genomes": len(genomes), "failed_imports": summary.get("failed_imports", []),
            "record_classes": result["record_classes"],
            "records_with_ambiguous_bases": result.get("records_with_ambiguous_bases"),
            "outputs": str(out / "reannotation")}
