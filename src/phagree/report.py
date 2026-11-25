"""Per-genome report, lysogeny markers and NCBI submission files for a ``phagree annotate`` output.

    python -m phagree.report <annotate output> [--locus-tag-prefix ABC12]

Reads ``reannotation/consensus_genes.csv.gz``, ``records_reannotation.csv.gz``, ``sequence_quality.csv.gz``
and the genome sequences in ``genbank/``, and writes

* ``report/genome_summary.tsv``: one row per genome;
* ``report/lysogeny_markers.tsv``: marker genes of temperate phages (see :mod:`phagree.lifestyle`);
* ``ncbi/<genome>.gbk``: GenBank file with the consensus annotation;
* ``ncbi/<genome>.tbl`` and ``.fsa``: feature table and sequence for NCBI table2asn / BankIt;
* ``ncbi/<genome>.faa``: proteins of the consensus genes.

Consensus genes are complete CDS (partial and origin-spanning models are not part of the consensus).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import duckdb
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .lifestyle import genome_markers

CATEGORIES = ["head and packaging", "connector", "tail", "lysis", "DNA, RNA and nucleotide metabolism",
              "transcription regulation", "integration and excision",
              "moron, auxiliary metabolic gene and host takeover", "other"]
START_CODONS_11 = {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}


def seq_id(accession):
    """Sequence identifier safe for NCBI files."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", accession)


def translate(seq, start1, end1, strand):
    """Protein of a complete CDS (translation table 11; the start codon is read as M)."""
    nt = seq[start1 - 1:end1]
    if strand == -1:
        nt = nt.reverse_complement()
    try:
        return str(nt.translate(table=11, cds=True))
    except Exception:  # ambiguous bases or an unusual start: translate codon by codon
        aa = str(nt[: len(nt) // 3 * 3].translate(table=11)).rstrip("*")
        if aa and str(nt[:3]).upper() in START_CODONS_11:
            aa = "M" + aa[1:]
        return aa


def coding_fraction(genes, length):
    covered = bytearray(length)
    for s, e in zip(genes.start1, genes.end1):
        covered[s - 1:e] = b"\x01" * (e - s + 1)
    return sum(covered) / length if length else 0.0


def tool_versions(run):
    path = (run / "parquet" / "source_runs").as_posix()
    rows = duckdb.connect().execute(f"""SELECT source, any_value(tool_version) FROM read_parquet('{path}/*.parquet',
        union_by_name=true) WHERE source IN ('phanotate', 'pyrodigal_gv') GROUP BY source""").fetchall()
    return dict(rows)


def build(out, locus_tag_prefix=None):
    out = Path(out)
    res = out / "reannotation"
    genes = pd.read_csv(res / "consensus_genes.csv.gz")
    records = pd.read_csv(res / "records_reannotation.csv.gz").set_index("accession_version")
    quality = pd.read_csv(res / "sequence_quality.csv.gz").set_index("accession_version")
    fasta = set((out / "fasta_genomes.txt").read_text().split()) if (out / "fasta_genomes.txt").exists() else set()
    manifest = pd.read_csv(out / "manifest.csv")
    versions = tool_versions(out / "run")
    inference = {"pha": f"ab initio prediction:PHANOTATE:{versions.get('phanotate', '')}".rstrip(":"),
                 "pyr": f"ab initio prediction:Pyrodigal-gv:{versions.get('pyrodigal_gv', '')}".rstrip(":")}
    (out / "report").mkdir(exist_ok=True)
    (out / "ncbi").mkdir(exist_ok=True)
    summary, markers = [], []
    tag = 0
    for acc, path in zip(manifest.accession_version, manifest.genbank_path):
        record = SeqIO.read(path, "genbank")
        seq = record.seq.upper()
        g = genes[genes.accession_version == acc].sort_values(["start1", "end1", "strand"])
        sid = seq_id(acc)
        circular = record.annotations.get("topology") == "circular"
        features = [f for f in record.features if f.type == "source"][:1]
        if not features:
            features = [SeqFeature(FeatureLocation(0, len(seq), strand=1), type="source",
                                   qualifiers={"mol_type": ["genomic DNA"]})]
        tbl = [f">Feature {sid}\n"]
        faa = []
        rows = []
        for r in g.itertuples():
            tag += 1
            lt = f"{locus_tag_prefix}_{tag:05d}" if locus_tag_prefix else None
            name = lt or f"{sid}_{r.start1}_{r.end1}_{'+' if r.strand == 1 else '-'}"
            protein = translate(seq, int(r.start1), int(r.end1), int(r.strand))
            inf = [inference[k] for k in ("pha", "pyr") if bool(getattr(r, k))]
            notes = []
            if isinstance(r.phrog, str):
                notes.append(f"PHROGs {r.phrog} ({r.category})"
                             + ("; assigned with the long-profile rule" if r.function_rule == "long_profile" else ""))
            if not inf:
                notes.append("GenBank CDS without gene caller support, kept for its PHROGs function")
            if r.ambiguous_bases:
                notes.append(f"contains {int(r.ambiguous_bases)} ambiguous bases")
            first, last = (r.start1, r.end1) if r.strand == 1 else (r.end1, r.start1)
            if lt:
                tbl.append(f"{first}\t{last}\tgene\n\t\t\tlocus_tag\t{lt}\n")
            tbl.append(f"{first}\t{last}\tCDS\n\t\t\tproduct\t{r.product}\n\t\t\ttransl_table\t11\n")
            tbl += [f"\t\t\tinference\t{x}\n" for x in inf]
            tbl += [f"\t\t\tnote\t{x}\n" for x in notes]
            q = {"product": [r.product], "codon_start": ["1"], "transl_table": ["11"], "translation": [protein]}
            if lt:
                q = {"locus_tag": [lt], **q}
                features.append(SeqFeature(FeatureLocation(int(r.start1) - 1, int(r.end1), strand=int(r.strand)),
                                           type="gene", qualifiers={"locus_tag": [lt]}))
            if inf:
                q["inference"] = inf
            if notes:
                q["note"] = notes
            features.append(SeqFeature(FeatureLocation(int(r.start1) - 1, int(r.end1), strand=int(r.strand)),
                                       type="CDS", qualifiers=q))
            faa.append(f">{name} {r.product} [{sid}]\n" + "\n".join(protein[i:i + 60] for i in range(0, len(protein), 60)) + "\n")
            rows.append({"start": int(r.start1), "end": int(r.end1), "strand": int(r.strand), "product": r.product,
                         "phrog": r.phrog if isinstance(r.phrog, str) else None, "locus": name})
        annotations = {"molecule_type": "DNA", "topology": "circular" if circular else "linear",
                       "data_file_division": "PHG"}
        for key in ("organism", "source", "taxonomy"):
            if key in record.annotations:
                annotations[key] = record.annotations[key]
        out_rec = SeqRecord(seq, id=acc, name=record.name[:16], description=record.description or sid,
                            annotations=annotations, features=features)
        SeqIO.write(out_rec, out / "ncbi" / f"{sid}.gbk", "genbank")
        (out / "ncbi" / f"{sid}.tbl").write_text("".join(tbl))
        title = record.description.strip()
        if title in ("", ".", "<unknown description>", record.id, record.name, record.id.split(".")[0]):
            title = ""  # FASTA input: no description beyond the identifier
        complete = circular and "complete genome" in title.lower()
        mods = (" [gcode=11]" + (" [topology=circular]" if circular else "")
                + (" [completeness=complete]" if complete else "") + (f" {title}" if title else ""))
        (out / "ncbi" / f"{sid}.fsa").write_text(
            f">{sid}{mods}\n" + "\n".join(str(seq[i:i + 70]) for i in range(0, len(seq), 70)) + "\n")
        (out / "ncbi" / f"{sid}.faa").write_text("".join(faa))

        level, found = genome_markers(rows)
        markers += [{"genome": acc, "group": m["group"], "selected": m["selected"], "start": m["start"],
                     "end": m["end"], "strand": "+" if m["strand"] == 1 else "-", "product": m["product"],
                     "phrog": m["phrog"], "locus": m["locus"]} for m in found]
        rec = records.loc[acc]
        qual = quality.loc[acc]
        gc = (seq.count("G") + seq.count("C")) / max(1, len(seq) - seq.count("N"))
        row = {
            "genome": acc, "input": "fasta" if acc in fasta else "genbank",
            "length_bp": len(seq), "gc_percent": round(100 * gc, 2), "topology": annotations["topology"],
            "ambiguous_bases": int(qual.n_bases + qual.other_ambiguous_bases),
            "consensus_genes": int(len(g)),
            "coding_percent": round(100 * coding_fraction(g, len(seq)), 1),
            "genes_with_function": int(g.phrog.notna().sum()),
            "genes_hypothetical": int(g.phrog.isna().sum()),
        }
        row.update({f"genes_{c}": int((g.category == c).sum()) for c in CATEGORIES})
        row.update({
            "lysogeny_evidence": level,
            "lysogeny_markers": "; ".join(f"{m['group']} {m['start']}-{m['end']}({'+' if m['strand'] == 1 else '-'})"
                                          for m in found if m["selected"]),
            "record_class": rec.reannotation_class,
            "genbank_cds_assessed": int(rec.gb_cds_assessed),
            "genbank_cds_confirmed_2_callers": int(rec.confirmed),
            "genbank_cds_supported_1_caller": int(rec.supported),
            "genbank_cds_flagged": int(rec.gb_only_unsupported),
            "proposed_additions": int(rec.add_A),
            "start_corrections": int(rec.start_corrections),
            "functional_names_for_hypothetical": int(rec.functional_relabels),
            "genes_with_ambiguous_bases": int(rec.genes_with_ambiguous_bases)
            if pd.notna(rec.get("genes_with_ambiguous_bases")) else 0,
        })
        summary.append(row)
    pd.DataFrame(summary).to_csv(out / "report" / "genome_summary.tsv", sep="\t", index=False)
    pd.DataFrame(markers, columns=["genome", "group", "selected", "start", "end", "strand", "product", "phrog",
                                   "locus"]).to_csv(out / "report" / "lysogeny_markers.tsv", sep="\t", index=False,
                                                    quoting=csv.QUOTE_MINIMAL)
    return {"genomes": len(summary), "report": str(out / "report"), "ncbi": str(out / "ncbi")}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("out", help="output folder of phagree annotate")
    p.add_argument("--locus-tag-prefix")
    a = p.parse_args(argv)
    print(json.dumps(build(a.out, a.locus_tag_prefix), indent=1))


if __name__ == "__main__":
    main()
