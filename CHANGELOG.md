# Changelog

## 1.1.0 (2026-11-25)

- New per-genome report from `phagree annotate`: `report/genome_summary.tsv` (length, GC, coding density, genes by PHROGs category, lysogeny evidence, record class and gene counts) and `report/lysogeny_markers.tsv`.
- Lysogeny markers (`phagree.lifestyle`): `strong` evidence (integrase), `weak` (transcriptional repressor or anti-repressor only) or `none`. Marker groups were selected on one half of the relatedness blocks of the study cohort and evaluated on the other half against the INPHARED2 lifestyle labels; see `docs/VALIDATION_1.1.md`.
- NCBI submission files in `ncbi/`: GenBank file with the consensus annotation, five-column feature table and sequence for table2asn/BankIt (validated with table2asn), and protein FASTA. `--locus-tag-prefix` adds locus tags.
- `--function-rule extended`: PHROGs functions also for proteins that align with at least 80% of their length to at least 30% of a longer profile (fusions, long multi-domain families); the default `study` keeps the rule of the study. The PHROGs search writes `protein_function_extended.parquet` next to the existing tables.
- Genomes given as FASTA are reported as `draft` records, and their gene labels no longer mention GenBank.
- `reannotation/consensus_genes.csv.gz`: the consensus genes as a table, with the callers that predict each gene. Existing outputs are unchanged.
- PHROGs search: E-values are computed for a fixed search space (Z = 20,000 proteins, the size of a search block of the study; domZ = 13, the hit-weighted median number of significant proteins per profile in those blocks). Before, HMMER took both from the proteins searched together, so the hits of a genome, and borderline functions, depended on how many genomes were annotated in one run (on the 28 laboratory genomes, 6 of 3,215 proteins). Re-searching 300 proteins of a full study block reproduces the study E-values exactly; functions agree with the study for 99.8% of 2,045 proteins tested. `--z auto --domz auto` restores the 1.0.0 behaviour.
- The PHROGs cache records its search space and refuses to mix hits computed with different settings; caches without a recorded search space (the study cache, 1.0.0) are read as before. `phagree annotate` uses a new cache folder, `phrogs_cache_v4_z20000_domz13`.
- Conda recipe: test command without inner quotes (the BioContainers test wraps commands in single quotes).

## 1.0.0 (2026-11-25)

First release.

- `phagree setup`: downloads the PHROGs v4 alignments (PHROGs site, or a Zenodo mirror; transient errors are retried), checks their SHA-256, builds the HMMER3 library in parallel and verifies its content against the library of the study; finds PHANOTATE and writes `run_config.json` into the database folder (`--db`, `$PHAGREE_DB`, `data/db` of a repository checkout, or `~/.local/share/phagree`). A verified library that is already present is reused.
- `phagree annotate`: one command for any genomes (GenBank accessions, GenBank files, FASTA files of new genomes): gene calling with PHANOTATE and Pyrodigal-gv, PHROGs search and consensus re-annotation.
- Ambiguous bases (N or other IUPAC codes): CDS whose only problem is an ambiguous codon are compared by coordinates and searched as X-masked proteins; genes and records report their ambiguous bases, and `sequence_quality.csv.gz` lists them per genome.
- Conda package (recipe in `conda/meta.yaml`); PHANOTATE runs from the same environment.
- The complete study: cohort selection, gene calling of 17,882 genomes, statistics, record-age analysis and consensus re-annotation of 17,859 GenBank records (`docs/`, `results/`, `figures/`).
