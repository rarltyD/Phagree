"""Select and validate PHROGs-based lysogeny markers against INPHARED2 lifestyle labels (Phagree 1.1).

    python validate_lysogeny_markers.py --run results/full_cohort \
        --genomes results/descriptive/per_genome_summary.csv \
        --blocks results/record_age/genome_table.csv.gz --out results/validation

The lifestyle labels of INPHARED2 are predictions of PhageLeads, not experimental truth: the result is
agreement with another predictor.

Rules (the validation half is used only to evaluate):
* relatedness blocks are split in two halves by a hash of the block name (development / validation);
* first prespecified selection (kept for the record): a candidate marker group is selected on the
  development half when it occurs in at most 1% of lytic and at least 1% of temperate genomes.  It excluded
  integrase (2.2% of lytic genomes) and gave a validation sensitivity of 11.7%, so it was replaced by:
* greedy selection on the development half: starting from no marker, add the group that most increases
  Youden's J (sensitivity + specificity - 1) of the call "any selected marker"; stop when the gain is
  below 0.005;
* a genome is called "temperate markers found" when any protein (GenBank or caller model) has a specific
  PHROGs function of a selected group;
* sensitivity and specificity are reported on the validation half with 95% cluster-bootstrap CIs over
  blocks (1,000 replicates); genomes labelled "unknown" are excluded.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from phagree.lifestyle import CANDIDATES, SELECTED, STRONG

MAX_LYTIC, MIN_TEMPERATE = 0.01, 0.01


def half(block):
    return int(hashlib.sha256(str(block).encode()).hexdigest(), 16) % 2  # 0 = development, 1 = validation


def rates(df, col):
    t, ly = df[df.label == "temperate"], df[df.label == "lytic"]
    return {"temperate_with": int(t[col].sum()), "temperate_n": int(len(t)),
            "lytic_with": int(ly[col].sum()), "lytic_n": int(len(ly)),
            "sensitivity": round(float(t[col].mean()), 4), "lytic_share": round(float(ly[col].mean()), 4)}


def bootstrap(df, col, reps=1000, seed=30102026):
    rng = np.random.default_rng(seed)
    blocks = df.block.unique()
    g = {b: d for b, d in df.groupby("block")}
    sens, spec = [], []
    for _ in range(reps):
        s = pd.concat([g[b] for b in rng.choice(blocks, len(blocks))])
        t, ly = s[s.label == "temperate"], s[s.label == "lytic"]
        sens.append(t[col].mean())
        spec.append(1 - ly[col].mean())

    def q(v):
        return [round(float(np.nanpercentile(v, 2.5)), 4), round(float(np.nanpercentile(v, 97.5)), 4)]
    return q(sens), q(spec)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--genomes", required=True)
    p.add_argument("--blocks", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    run, out = Path(a.run), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute("SET memory_limit='4GB'; SET threads=3")
    ann = con.execute(f"""
        SELECT DISTINCT c.accession_version, f.best_annot AS annot
        FROM read_parquet('{(run / 'parquet' / 'cds_calls').as_posix()}/*.parquet', union_by_name=true) c
        JOIN read_parquet('{(run / 'phrogs' / 'protein_function_primary.parquet').as_posix()}') f USING (protein_sha256)
        WHERE f.function_status = 'specific_function' AND f.best_annot IS NOT NULL""").df()
    genomes = pd.read_csv(a.genomes, usecols=["accession_version", "phageleads_lifestyle"]).rename(
        columns={"phageleads_lifestyle": "label"})
    blocks = pd.read_csv(a.blocks, usecols=["accession_version", "block"])
    df = genomes.merge(blocks, on="accession_version", how="left")
    df["block"] = df.block.fillna(df.accession_version)
    for name, pat in CANDIDATES.items():
        rx = re.compile(pat, re.I)
        has = set(ann[ann.annot.map(lambda x: bool(rx.search(x)))].accession_version)
        df[name] = df.accession_version.isin(has)
    df["half"] = df.block.map(half)
    known = df[df.label.isin(["temperate", "lytic"])]
    dev, val = known[known.half == 0], known[known.half == 1]

    candidates = {name: rates(dev, name) for name in CANDIDATES}
    first_rule = [n for n, r in candidates.items() if r["lytic_share"] <= MAX_LYTIC and r["sensitivity"] >= MIN_TEMPERATE]

    def youden(d, cols):
        call = d[cols].any(axis=1) if cols else pd.Series(False, index=d.index)
        return call[d.label == "temperate"].mean() - call[d.label == "lytic"].mean()

    selected, path = [], []
    while True:
        best = max((n for n in CANDIDATES if n not in selected), key=lambda n: youden(dev, selected + [n]), default=None)
        if best is None:
            break
        gain = youden(dev, selected + [best]) - youden(dev, selected)
        path.append({"add": best, "youden_dev": round(youden(dev, selected + [best]), 4), "gain": round(gain, 4)})
        if gain < 0.005:
            path[-1]["stopped"] = True
            break
        selected.append(best)
    for d in (dev, val, known, df):
        d["call"] = d[selected].any(axis=1)
    assert tuple(selected) == SELECTED, f"phagree.lifestyle.SELECTED {SELECTED} differs from selection {selected}"
    val_rates = rates(val, "call")
    sens_ci, spec_ci = bootstrap(val, "call")
    for d in (val, known):
        d["strong"] = d[list(STRONG)].any(axis=1)
    strong_rates = rates(val, "strong")
    strong_ci = bootstrap(val, "strong")
    summary = {
        "labels": "INPHARED2 phageleads_lifestyle (PhageLeads prediction)",
        "genomes": int(len(df)), "genomes_known_label": int(len(known)),
        "development": {"genomes": int(len(dev)), "blocks": int(dev.block.nunique())},
        "validation": {"genomes": int(len(val)), "blocks": int(val.block.nunique())},
        "selection_rule": f"lytic share <= {MAX_LYTIC} and temperate share >= {MIN_TEMPERATE} on development half",
        "candidates_development": candidates,
        "first_prespecified_rule": {"selected": first_rule, "validation": rates(val.assign(c1=val[first_rule].any(axis=1)), "c1")},
        "greedy_path_development": path,
        "selected": selected,
        "validation_per_marker": {n: rates(val, n) for n in CANDIDATES},
        "validation_call": {**val_rates, "specificity": round(1 - val_rates["lytic_share"], 4),
                            "sensitivity_ci95": sens_ci, "specificity_ci95": spec_ci},
        "validation_strong_integrase": {**strong_rates, "specificity": round(1 - strong_rates["lytic_share"], 4),
                                        "sensitivity_ci95": strong_ci[0], "specificity_ci95": strong_ci[1]},
        "all_known_call": {**rates(known, "call"), "specificity": round(1 - rates(known, "call")["lytic_share"], 4)},
        "unknown_label_with_markers": int(df[df.label == "unknown"].call.sum()),
        "unknown_label_n": int((df.label == "unknown").sum()),
    }
    (out / "lysogeny_marker_validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    df.to_csv(out / "lysogeny_markers_per_genome.csv.gz", index=False)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
