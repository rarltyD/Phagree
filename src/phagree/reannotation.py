"""Consensus re-annotation of GenBank phage records from two gene callers and PHROGs.

Rules were fixed before the results were inspected.  Genes are keyed by genome, strand and the
coordinates of the stop codon (the stop-locus definition of the study); only primary-eligible
calls are assessed (partial, joined, slippage and origin-spanning models are left as they are).

Gene level
  confirmed          GenBank CDS whose stop is predicted by both callers
  supported          GenBank CDS whose stop is predicted by one caller
  gb_only_function   GenBank CDS predicted by neither caller, protein has a specific PHROGs function
  gb_only_unsupported GenBank CDS predicted by neither caller, no specific function  -> review
  add_A              both callers, absent from GenBank, specific PHROGs function     -> add
  add_B              both callers, absent from GenBank, no specific function, >=100 aa (reported only)
  start_correction   shared stop, callers agree on one start, GenBank start differs  (reported only)
  functional_relabel GenBank CDS with a hypothetical-like product whose protein has a specific
                     PHROGs function (same product rule as the primary analysis; added 2026-11-14,
                     does not change gene categories or record classes)
Record level
  unannotated        no CDS in the GenBank record
  major              (gb_only_unsupported + add_A) >= 10% of GenBank CDS, or discordant with both callers
                     (stop-set Jaccard < 0.8 with each caller)
  minor              at least one gb_only_unsupported or add_A gene
  consistent         none
Outputs: per-record table, per-gene proposals, consensus GFF3 for all records, summary JSON.
Exploratory analysis added after the primary analysis (see docs/PROTOCOL.md).

    python -m phagree.reannotation results/full_cohort results/reannotation \
        --genome-table results/record_age/genome_table.csv.gz

Without --genome-table, the per-genome discordance flags are computed from the run itself.  For runs
outside the primary cohort (phagree run --diagnostic), use --locus-eligibility geometry: the
cohort flag of such loci is false by design, while their geometric eligibility is kept.

--include-ambiguous also assesses CDS whose only problem is an ambiguous codon (N or another IUPAC
code): they are compared by coordinates, and their X-masked proteins are searched when the run was made
with "search_ambiguous_proteins": true.  Every gene and record then reports the ambiguous bases it
contains, and sequence_quality.csv.gz lists them per genome.
"""
import argparse
import bisect
import gzip
import json
import re
from pathlib import Path

import duckdb
import pandas as pd

AMBIGUOUS_RUN = re.compile(r"[^ACGT]+")
# Hypothetical-like GenBank products (same rule as scripts/analysis/extra_stats.py).
HYPO_SQL = r"hypothetical|uncharacteri[sz]ed|unknown|putative protein$|^protein$|duf[0-9]+"

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("run")
    p.add_argument("out")
    p.add_argument("--genome-table", help="genome_table.csv.gz of the record-age analysis")
    p.add_argument("--locus-eligibility", choices=("primary", "geometry"), default="primary")
    p.add_argument("--include-ambiguous", action="store_true")
    a = p.parse_args(argv)
    run, out = Path(a.run), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    def pq(table):
        return (run / "parquet" / table / "*.parquet").as_posix()

    fn_path = (run / "phrogs" / "protein_function_primary.parquet").as_posix()
    con = duckdb.connect()
    con.execute("SET memory_limit='4GB'; SET threads=3; SET preserve_insertion_order=false")

    if a.genome_table:
        genomes_a = pd.read_csv(a.genome_table, usecols=["accession_version", "submitted", "period",
                                                         "unannotated", "disc_annotated", "block",
                                                         "caller", "seq_tech"])
    else:
        genomes_a = con.execute(f"""
            WITH g AS (SELECT DISTINCT genome_id, accession_version FROM read_parquet('{pq("genomes")}', union_by_name=true)),
            gb AS (SELECT genome_id, count(*) FILTER (WHERE source = 'genbank') AS n
                   FROM read_parquet('{pq("cds_calls")}', union_by_name=true) GROUP BY genome_id),
            j AS (SELECT genome_id, max(stop_jaccard) AS best FROM read_parquet('{pq("genome_metrics")}', union_by_name=true)
                  WHERE metric_scope = 'source_pair' AND 'genbank' IN (source_a, source_b) GROUP BY genome_id)
            SELECT g.accession_version, NULL AS submitted, NULL AS period,
                   (coalesce(gb.n, 0) = 0)::INT AS unannotated,
                   (coalesce(gb.n, 0) > 0 AND coalesce(j.best, 0) < 0.8)::INT AS disc_annotated,
                   NULL AS block, NULL AS caller, NULL AS seq_tech
            FROM g LEFT JOIN gb USING (genome_id) LEFT JOIN j USING (genome_id)""").df()
    con.register("cohort", genomes_a)
    locus_flag = "eligible_primary" if a.locus_eligibility == "primary" else "geometry_eligible_primary"

    # One row per primary-eligible CDS call with its stop key, start, protein function and GenBank product.
    con.execute(f"""
    CREATE TEMP TABLE calls AS
    SELECT c.genome_id, c.accession_version, c.source, c.call_id, c.strand,
           c.stop_triplet::VARCHAR AS stop_key,
           list_min(list_transform(c.segments, s -> s.start0)) AS lo,
           list_max(list_transform(c.segments, s -> s.end0)) AS hi,
           c.feature_nt_length // 3 - 1 AS length_aa, c.protein_sha256,
           CASE WHEN c.source = 'genbank' THEN json_extract_string(c.qualifiers_json, '$.product[0]') END AS product,
           f.function_status, f.best_phrog, f.best_category, f.best_annot, f.best_evalue,
           NOT c.eligible_primary AS ambiguous
    FROM read_parquet('{pq("cds_calls")}', union_by_name=true) c
    JOIN cohort USING (accession_version)
    LEFT JOIN read_parquet('{fn_path}') f USING (protein_sha256)
    WHERE c.eligible_primary
       OR ({str(a.include_ambiguous).lower()} AND c.primary_exclusion_reasons = ['ambiguous_translation'])
    """)
    con.execute("""
    CREATE TEMP TABLE calls2 AS
    SELECT *, CASE WHEN strand = 1 THEN lo ELSE hi END AS start_coord FROM calls
    """)
    # Per stop key: which sources predict it, their starts, and the best function among the models.
    con.execute("""
    CREATE TEMP TABLE keys AS
    SELECT genome_id, any_value(accession_version) AS accession_version, strand, stop_key,
           bool_or(source = 'genbank') AS gb, bool_or(source = 'phanotate') AS pha,
           bool_or(source = 'pyrodigal_gv') AS pyr,
           min(start_coord) FILTER (WHERE source = 'genbank') AS start_gb,
           min(start_coord) FILTER (WHERE source = 'phanotate') AS start_pha,
           min(start_coord) FILTER (WHERE source = 'pyrodigal_gv') AS start_pyr,
           max(length_aa) AS length_aa,
           min(lo) AS lo, max(hi) AS hi,
           any_value(product) FILTER (WHERE source = 'genbank') AS gb_product,
           arg_max(protein_sha256, length_aa) FILTER (WHERE source = 'genbank') AS gb_protein,
           arg_max(protein_sha256, length_aa) FILTER (WHERE source <> 'genbank') AS caller_protein,
           bool_or(function_status = 'specific_function') FILTER (WHERE source = 'genbank') AS gb_specific,
           bool_or(function_status = 'specific_function') FILTER (WHERE source <> 'genbank') AS caller_specific,
           arg_min(best_annot, best_evalue) FILTER (WHERE source = 'genbank' AND function_status = 'specific_function') AS gb_annot,
           arg_min(best_annot, best_evalue) FILTER (WHERE function_status = 'specific_function') AS annot,
           arg_min(best_category, best_evalue) FILTER (WHERE function_status = 'specific_function') AS category,
           arg_min(best_phrog, best_evalue) FILTER (WHERE function_status = 'specific_function') AS phrog,
           bool_or(ambiguous) AS ambiguous_codon
    FROM calls2
    GROUP BY genome_id, strand, stop_key
    """)
    # Caller support counts only at eligible stop-loci (the study's unit), so caller totals match the
    # primary analysis exactly.
    con.execute(f"""
    CREATE TEMP TABLE keys_el AS
    SELECT k.* REPLACE (k.pha AND (el.genome_id IS NOT NULL OR k.ambiguous_codon) AS pha,
                        k.pyr AND (el.genome_id IS NOT NULL OR k.ambiguous_codon) AS pyr)
    FROM keys k
    LEFT JOIN (SELECT DISTINCT genome_id, strand, stop_triplet::VARCHAR AS stop_key
               FROM read_parquet('{pq("stop_loci")}', union_by_name=true) WHERE {locus_flag}) el
      USING (genome_id, strand, stop_key)
    """)
    con.execute("DELETE FROM keys_el WHERE NOT gb AND NOT pha AND NOT pyr")
    con.execute(f"""
    CREATE TEMP TABLE genes AS
    SELECT *, (pha::INT + pyr::INT) AS n_callers,
      CASE
        WHEN gb AND pha AND pyr THEN 'confirmed'
        WHEN gb AND (pha OR pyr) THEN 'supported'
        WHEN gb AND coalesce(gb_specific, false) THEN 'gb_only_function'
        WHEN gb THEN 'gb_only_unsupported'
        WHEN pha AND pyr AND coalesce(caller_specific, false) THEN 'add_A'
        WHEN pha AND pyr AND length_aa >= 100 THEN 'add_B'
        WHEN pha AND pyr THEN 'two_caller_other'
        WHEN coalesce(caller_specific, false) THEN 'one_caller_function'
        ELSE 'one_caller_other'
      END AS category_rule,
      (gb AND pha AND pyr AND start_pha = start_pyr AND start_gb <> start_pha) AS start_correction,
      (gb AND coalesce(gb_specific, false)
          AND regexp_matches(lower(coalesce(gb_product, '')), '{HYPO_SQL}')) AS functional_relabel
    FROM keys_el
    """)

    genes = con.execute("SELECT * FROM genes").df()

    # Ambiguous bases (anything other than A, C, G, T) per genome and per gene.
    # Only genomes that contain ambiguous bases are loaded (sequences of a whole cohort do not fit in memory).
    con.execute(f"""CREATE TEMP TABLE gseq AS
        SELECT DISTINCT g.accession_version, length(g.sequence) AS length_bp,
               regexp_matches(g.sequence, '[^ACGTacgt]')::INT AS n_ambiguous
        FROM read_parquet('{pq("genomes")}', union_by_name=true) g JOIN cohort USING (accession_version)""")
    seqs = con.execute(f"""SELECT DISTINCT g.accession_version, upper(g.sequence) AS sequence
        FROM read_parquet('{pq("genomes")}', union_by_name=true) g
        JOIN (SELECT accession_version FROM gseq WHERE n_ambiguous > 0) USING (accession_version)""").df()
    amb_runs, quality = {}, []
    for acc, seq in zip(seqs.accession_version, seqs.sequence):
        runs = [(m.start(), m.end()) for m in AMBIGUOUS_RUN.finditer(seq)]
        amb_runs[acc] = ([r[0] for r in runs], [r[1] for r in runs])
        n_amb = sum(e - s_ for s_, e in runs)
        quality.append({"accession_version": acc, "length_bp": len(seq), "n_bases": seq.count("N"),
                        "other_ambiguous_bases": n_amb - seq.count("N"), "ambiguous_segments": len(runs),
                        "segments_1based": ";".join(f"{s_ + 1}-{e}" for s_, e in runs[:500])})


    def ambiguous_in(acc, lo, hi):
        """Number of ambiguous bases in [lo, hi) (0-based, half-open)."""
        starts, ends = amb_runs.get(acc, ([], []))
        i = bisect.bisect_right(ends, lo)
        n = 0
        while i < len(starts) and starts[i] < hi:
            n += min(ends[i], hi) - max(starts[i], lo)
            i += 1
        return n


    genes["ambiguous_bases"] = [ambiguous_in(a_, lo, hi) for a_, lo, hi in zip(genes.accession_version, genes.lo, genes.hi)]
    clean = con.execute("SELECT accession_version, length_bp FROM gseq WHERE n_ambiguous = 0").df()
    quality = pd.concat([pd.DataFrame(quality, columns=["accession_version", "length_bp", "n_bases",
                                                        "other_ambiguous_bases", "ambiguous_segments",
                                                        "segments_1based"]),
                         clean.assign(n_bases=0, other_ambiguous_bases=0, ambiguous_segments=0, segments_1based="")])
    quality = quality.merge(genes[genes.ambiguous_bases > 0].groupby("accession_version").size()
                            .rename("genes_with_ambiguous_bases").reset_index(), how="left").fillna({"genes_with_ambiguous_bases": 0})
    quality.sort_values("accession_version").to_csv(out / "sequence_quality.csv.gz", index=False)
    counts = genes.pivot_table(index="accession_version", columns="category_rule", values="stop_key",
                               aggfunc="size", fill_value=0)
    for c in ("confirmed", "supported", "gb_only_function", "gb_only_unsupported", "add_A", "add_B",
              "two_caller_other", "one_caller_function", "one_caller_other"):
        if c not in counts:
            counts[c] = 0
    rec = genomes_a.set_index("accession_version").join(counts, how="left").fillna(
        {c: 0 for c in counts.columns})
    for flag, col in (("start_correction", "start_corrections"), ("functional_relabel", "functional_relabels")):
        rec[col] = genes[genes[flag].eq(True)].groupby("accession_version").size()
        rec[col] = rec[col].fillna(0).astype(int)
    rec["gb_cds_assessed"] = rec[["confirmed", "supported", "gb_only_function", "gb_only_unsupported"]].sum(axis=1)
    rec["flagged"] = rec.gb_only_unsupported + rec.add_A
    rec["flagged_share"] = rec.flagged / rec.gb_cds_assessed.where(rec.gb_cds_assessed > 0)


    def record_class(r):
        if r.unannotated == 1:
            return "unannotated"
        if r.flagged_share >= 0.10 or r.disc_annotated == 1:
            return "major"
        if r.flagged > 0:
            return "minor"
        return "consistent"


    rec["reannotation_class"] = rec.apply(record_class, axis=1)
    rec = rec.join(quality.set_index("accession_version")[["n_bases", "other_ambiguous_bases",
                                                           "genes_with_ambiguous_bases"]], how="left")
    rec.sort_index().reset_index().to_csv(out / "records_reannotation.csv.gz", index=False)

    # Per-gene proposals (everything except plain confirmed/supported genes and one-caller noise).
    prop = genes[genes.category_rule.isin(["gb_only_unsupported", "gb_only_function", "add_A", "add_B"])
                 | (genes.start_correction.eq(True)) | (genes.functional_relabel.eq(True))].copy()
    prop["start1"] = prop.lo + 1
    prop["end1"] = prop.hi
    prop[["accession_version", "category_rule", "start_correction", "functional_relabel", "strand", "start1",
          "end1", "length_aa", "start_gb", "start_pha", "start_pyr", "gb_product", "phrog", "category", "annot",
          "gb_protein", "caller_protein", "ambiguous_bases"]] \
        .sort_values(["accession_version", "start1"]).to_csv(out / "gene_proposals.csv.gz", index=False)

    # Consensus GFF3: two-caller genes, one-caller genes with a specific function, GenBank-only genes
    # with a specific function.  Coordinates are those of the longest model at the stop.
    keep = genes[genes.category_rule.isin(["confirmed", "supported", "gb_only_function", "add_A", "add_B",
                                            "two_caller_other", "one_caller_function"])]
    with gzip.open(out / "consensus_annotation.gff3.gz", "wt") as g:
        g.write("##gff-version 3\n")
        for r in keep.sort_values(["accession_version", "lo"]).itertuples():
            ev = {"confirmed": "GenBank+2 callers", "supported": "GenBank+1 caller",
                  "gb_only_function": "GenBank+PHROGs", "add_A": "2 callers+PHROGs (absent from GenBank)",
                  "add_B": "2 callers (absent from GenBank)", "two_caller_other": "2 callers (absent from GenBank)",
                  "one_caller_function": "1 caller+PHROGs"}[r.category_rule]
            product = r.annot if isinstance(r.annot, str) and r.annot else "hypothetical protein"
            attrs = (f"ID={r.accession_version}_{r.lo + 1}_{'+' if r.strand == 1 else '-'};"
                     f"product={product.replace(';', ',')};evidence={ev}")
            if isinstance(r.phrog, str):
                attrs += f";phrog={r.phrog};phrog_category={str(r.category).replace(';', ',')}"
            if isinstance(r.gb_product, str):
                attrs += f";genbank_product={r.gb_product.replace(';', ',')}"
            if r.ambiguous_bases:
                attrs += f";ambiguous_bases={r.ambiguous_bases}"
            g.write(f"{r.accession_version}\tPhagree\tCDS\t{r.lo + 1}\t{r.hi}\t.\t"
                    f"{'+' if r.strand == 1 else '-'}\t0\t{attrs}\n")

    relabel = genes[genes.functional_relabel.eq(True)]
    summ = {
        "genomes": int(len(rec)),
        "gene_categories": genes.category_rule.value_counts().to_dict(),
        "start_corrections": int((genes.start_correction.eq(True)).sum()),
        "functional_relabels": int(len(relabel)),
        "functional_relabel_records": int(relabel.accession_version.nunique()),
        "functional_relabel_by_category": relabel.category.value_counts().head(12).to_dict(),
        "record_classes": rec.reannotation_class.value_counts().to_dict(),
        "records_with_add_A": int((rec.add_A > 0).sum()),
        "records_with_ambiguous_bases": int(((quality.n_bases + quality.other_ambiguous_bases) > 0).sum()),
        "genes_with_ambiguous_bases": int((genes.ambiguous_bases > 0).sum()),
        "records_with_gb_only_unsupported": int((rec.gb_only_unsupported > 0).sum()),
        "consensus_genes": int(len(keep)),
        "add_A_by_category": genes[genes.category_rule == "add_A"].category.value_counts().head(12).to_dict(),
        "gb_only_unsupported_length_aa_median": float(genes[genes.category_rule == "gb_only_unsupported"].length_aa.median()),
        "class_by_period": pd.crosstab(rec.period, rec.reannotation_class, normalize="index").round(4).to_dict(orient="index")
        if rec.period.notna().any() else None,
    }
    (out / "reannotation_summary.json").write_text(json.dumps(summ, indent=2, ensure_ascii=False, default=str))
    print(json.dumps({k: summ[k] for k in ("genomes", "gene_categories", "start_corrections", "functional_relabels",
                                            "functional_relabel_records", "record_classes", "records_with_add_A",
                                            "records_with_gb_only_unsupported", "consensus_genes")},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
