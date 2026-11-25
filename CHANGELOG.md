# Changelog

## 1.0.0 (2026-11-25)

First release.

- `phagree setup`: downloads the PHROGs v4 alignments (PHROGs site, or a Zenodo mirror; transient errors are retried), checks their SHA-256, builds the HMMER3 library in parallel and verifies its content against the library of the study; finds PHANOTATE and writes `run_config.json` into the database folder (`--db`, `$PHAGREE_DB`, `data/db` of a repository checkout, or `~/.local/share/phagree`). A verified library that is already present is reused.
- `phagree annotate`: one command for any genomes (GenBank accessions, GenBank files, FASTA files of new genomes): gene calling with PHANOTATE and Pyrodigal-gv, PHROGs search and consensus re-annotation.
- Ambiguous bases (N or other IUPAC codes): CDS whose only problem is an ambiguous codon are compared by coordinates and searched as X-masked proteins; genes and records report their ambiguous bases, and `sequence_quality.csv.gz` lists them per genome.
- Conda package (recipe in `conda/meta.yaml`); PHANOTATE runs from the same environment.
- The complete study: cohort selection, gene calling of 17,882 genomes, statistics, record-age analysis and consensus re-annotation of 17,859 GenBank records (`docs/`, `results/`, `figures/`).
