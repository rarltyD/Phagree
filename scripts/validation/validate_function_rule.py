"""Validate the long-profile rescue of PHROGs functions (Phagree 1.1) on the study cohort.

    python scripts/validation/validate_function_rule.py --run results/full_cohort \
        --cache data/db/phrogs_cache_v4 --annotations data/db/phrogs_v4/phrog_annot_v4.tsv \
        --out results/validation

1. Regression: the study rule recomputed from the raw hits must reproduce
   <run>/phrogs/protein_function_primary.parquet exactly.
2. Effect: proteins that gain a specific function with ``classify_extended`` (function_rule =
   "long_profile"); the study result of every other protein is unchanged by construction.
3. Agreement with GenBank: for GenBank CDS with a named (not hypothetical-like) product, the share whose
   product shares a function word with the PHROGs annotation, for study functions (baseline) and for
   rescued functions.  This is a coarse lexical check, not a measure of accuracy.
"""
import argparse
import json
import re
from pathlib import Path

import duckdb
import pandas as pd

from phagree.functional import classify, classify_extended, load_annotations

HYPO = re.compile(r"hypothetical|uncharacteri[sz]ed|unknown|putative protein$|^protein$|duf[0-9]+")
STOP = set("""protein proteins putative phage prophage like subunit family domain domains predicted possible
probable containing type gene product conserved and the with for from related superfamily associated
gp orf small large major minor component bacteriophage viral virus uncharacterized hypothetical""".split())


def words(text):
    return {w for w in re.findall(r"[a-z0-9]+", str(text).lower()) if len(w) >= 3 and w not in STOP}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--cache", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    run, out = Path(a.run), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute("SET memory_limit='4GB'; SET threads=3")
    hits = con.execute(f"""SELECT DISTINCT ON (protein_sha256, phrog) * FROM
        read_parquet('{Path(a.cache, 'hits').as_posix()}/*.parquet')""").df()
    ann = load_annotations(a.annotations)
    grouped = {k: g.to_dict("records") for k, g in hits.groupby("protein_sha256")}
    rows = []
    for sha, hs in grouped.items():
        s = classify(hs, ann)
        e = classify_extended(hs, ann)
        rows.append({"protein_sha256": sha, "study_status": s["function_status"], "study_phrog": s["best_phrog"],
                     "ext_status": e["function_status"], "ext_phrog": e["best_phrog"], "ext_annot": e["best_annot"],
                     "ext_category": e["best_category"], "function_rule": e["function_rule"]})
    res = pd.DataFrame(rows)

    # 1. Regression against the published table (proteins without hits are no_informative_hit there).
    # Only proteins searched in this cache (chunk member lists) can be recomputed; part of the cohort
    # was searched on another machine whose raw hits are archived elsewhere.
    members = set()
    for f in Path(a.cache, "hits").glob("*.members.json"):
        members.update(json.loads(f.read_text()))
    prim = pd.read_parquet(run / "phrogs" / "protein_function_primary.parquet",
                           columns=["protein_sha256", "function_status", "best_phrog"])
    m = prim[prim.protein_sha256.isin(members)].merge(res, on="protein_sha256", how="left")
    m["study_status"] = m.study_status.fillna("no_informative_hit")
    searched = m.function_status != "not_searched_or_failed"
    mism = m[searched & ((m.function_status != m.study_status)
                         | (m.best_phrog.fillna("") != m.study_phrog.fillna("")))]
    # Every protein with a study function must keep it.
    changed = res[res.study_status.isin(["specific_function", "ambiguous"])
                  & ((res.ext_status != res.study_status) | (res.ext_phrog != res.study_phrog))]

    # 2./3. Rescued proteins and agreement with named GenBank products.
    gb = con.execute(f"""SELECT protein_sha256, json_extract_string(qualifiers_json, '$.product[0]') AS product
        FROM read_parquet('{(run / 'parquet' / 'cds_calls').as_posix()}/*.parquet', union_by_name=true)
        WHERE source = 'genbank' AND protein_sha256 IS NOT NULL""").df()
    gb = gb[gb.protein_sha256.isin(members)].merge(res, on="protein_sha256", how="left")
    gb["ext_status"] = gb.ext_status.fillna("no_informative_hit")
    gb["function_rule"] = gb.function_rule.fillna("study")
    gb["named"] = gb["product"].notna() & ~gb["product"].fillna("").str.lower().str.contains(HYPO)

    def agreement(df):
        df = df[df.named]
        ok = [bool(words(p) & words(x)) for p, x in zip(df["product"], df.ext_annot)]
        return {"genbank_cds": int(len(df)), "share_word_agreement": round(sum(ok) / len(ok), 4) if ok else None}

    base = gb[gb.ext_status.eq("specific_function") & gb.function_rule.eq("study")]
    resc = gb[gb.function_rule.eq("long_profile")]
    rescued = res[res.function_rule.eq("long_profile")]
    hypo_named = resc[~resc.named]
    summary = {
        "proteins_searched_in_cache": int(len(members)),
        "published_proteins": int(len(prim)),
        "published_proteins_recomputed": int(len(m)),
        "proteins_with_hits": int(len(res)),
        "regression_mismatches_vs_published": int(len(mism)),
        "study_functions_changed_by_extension": int(len(changed)),
        "proteins_study_specific": int(res.study_status.eq("specific_function").sum()),
        "proteins_rescued_long_profile": int(len(rescued)),
        "rescued_by_previous_status": rescued.study_status.value_counts().to_dict(),
        "rescued_top_annotations": rescued.ext_annot.value_counts().head(15).to_dict(),
        "rescued_by_category": rescued.ext_category.value_counts().to_dict(),
        "agreement_named_genbank_study_rule": agreement(base),
        "agreement_named_genbank_rescued": agreement(resc),
        "genbank_cds_rescued_total": int(len(resc)),
        "genbank_cds_rescued_with_hypothetical_product": int(len(hypo_named)),
    }
    (out / "function_rule_validation.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    examples = resc[resc.named].assign(agree=[bool(words(p) & words(x)) for p, x in
                                              zip(resc[resc.named]["product"], resc[resc.named].ext_annot)])
    examples.drop_duplicates("protein_sha256").sample(min(60, len(examples)), random_state=1)[
        ["product", "ext_annot", "ext_phrog", "agree"]].to_csv(out / "function_rule_examples.csv", index=False)
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
