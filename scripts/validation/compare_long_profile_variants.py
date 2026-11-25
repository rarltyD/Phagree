"""Compare thresholds of the long-profile rescue on the study cohort (development step of Phagree 1.1).

    python compare_long_profile_variants.py --run results/full_cohort --cache data/db/phrogs_cache_v4 \
        --annotations data/db/phrogs_v4/phrog_annot_v4.tsv --out results/validation

Criterion fixed before running: choose the variant with the most rescued GenBank CDS among those whose
lexical agreement with named GenBank products is at least 0.78 (baseline of study functions: 0.796).
"""
import argparse
import json
import re
from pathlib import Path

import duckdb

from phagree.functional import classify_extended, load_annotations

HYPO = re.compile(r"hypothetical|uncharacteri[sz]ed|unknown|putative protein$|^protein$|duf[0-9]+")
STOP = set("""protein proteins putative phage prophage like subunit family domain domains predicted possible
probable containing type gene product conserved and the with for from related superfamily associated
gp orf small large major minor component bacteriophage viral virus uncharacterized hypothetical""".split())
VARIANTS = [(0.9, 0.0), (0.95, 0.0), (0.9, 0.2), (0.9, 0.3), (0.95, 0.3), (0.8, 0.3), (0.95, 0.4)]


# Synonyms applied identically to study and rescued functions (added after inspecting disagreements:
# the raw criterion was not met by any variant, and most disagreements were synonyms).
SYN = {"capsid": "head", "baseplate": "tail", "fibre": "fiber", "lysozyme": "endolysin", "lysin": "endolysin",
       "dutp": "dutpase", "diphosphatase": "dutpase", "pyrophosphatase": "dutpase", "spike": "tail",
       "tailspike": "tail", "portal": "portal", "neck": "connector", "connector": "connector"}


def words(text, synonyms=False):
    w = {w for w in re.findall(r"[a-z0-9]+", str(text).lower()) if len(w) >= 3 and w not in STOP}
    return {SYN.get(x, x) for x in w} if synonyms else w


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    p.add_argument("--cache", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    run, out = Path(a.run), Path(a.out)
    con = duckdb.connect()
    con.execute("SET memory_limit='3GB'; SET threads=3")
    members = set()
    for f in Path(a.cache, "hits").glob("*.members.json"):
        members.update(json.loads(f.read_text()))
    hits = con.execute(f"""SELECT DISTINCT ON (protein_sha256, phrog) * FROM
        read_parquet('{Path(a.cache, 'hits').as_posix()}/*.parquet')""").df()
    ann = load_annotations(a.annotations)
    grouped = {k: g.to_dict("records") for k, g in hits.groupby("protein_sha256")}
    del hits
    gb = con.execute(f"""SELECT protein_sha256, json_extract_string(qualifiers_json, '$.product[0]') AS product
        FROM read_parquet('{(run / 'parquet' / 'cds_calls').as_posix()}/*.parquet', union_by_name=true)
        WHERE source = 'genbank' AND protein_sha256 IS NOT NULL""").df()
    gb = gb[gb.protein_sha256.isin(members)]
    gb["named"] = gb["product"].notna() & ~gb["product"].fillna("").str.lower().str.contains(HYPO)
    def agreement(named, annot_of):
        raw = [bool(words(x) & words(annot_of[s])) for x, s in zip(named["product"], named.protein_sha256)]
        syn = [bool(words(x, True) & words(annot_of[s], True)) for x, s in zip(named["product"], named.protein_sha256)]
        return (round(sum(raw) / len(raw), 4), round(sum(syn) / len(syn), 4)) if raw else (None, None)

    from phagree.functional import classify
    study = {}
    for sha, hs in grouped.items():
        r = classify(hs, ann)
        if r["function_status"] == "specific_function":
            study[sha] = r["best_annot"]
    base = gb[gb.protein_sha256.isin(study) & gb.named]
    b_raw, b_syn = agreement(base, study)
    print("baseline (study functions):", len(base), b_raw, b_syn, flush=True)
    result = {"baseline_study_functions": {"named_genbank_cds": int(len(base)), "word_agreement": b_raw,
                                           "word_agreement_synonyms": b_syn}}
    for q, pc in VARIANTS:
        resc = {}
        for sha, hs in grouped.items():
            r = classify_extended(hs, ann, qcov_long=q, pcov_long=pc)
            if r["function_rule"] == "long_profile":
                resc[sha] = r["best_annot"]
        g = gb[gb.protein_sha256.isin(resc)]
        named = g[g.named]
        raw, syn = agreement(named, resc)
        result[f"qcov>={q},pcov>={pc}"] = {
            "proteins_rescued": len(resc), "genbank_cds_rescued": int(len(g)),
            "genbank_cds_rescued_hypothetical": int((~g.named).sum()),
            "named_genbank_cds": int(len(named)),
            "word_agreement": raw, "word_agreement_synonyms": syn}
        print(q, pc, result[f"qcov>={q},pcov>={pc}"], flush=True)
    # Prespecified criterion (raw words, >= 0.78) and, after it failed, the same criterion relative to the
    # baseline with synonyms: agreement at least baseline - 0.02.
    var = {k: v for k, v in result.items() if k != "baseline_study_functions"}
    ok = {k: v for k, v in var.items() if v["word_agreement"] is not None and v["word_agreement"] >= 0.78}
    choice = max(ok, key=lambda k: ok[k]["genbank_cds_rescued"]) if ok else None
    ok2 = {k: v for k, v in var.items() if v["word_agreement_synonyms"] is not None
           and v["word_agreement_synonyms"] >= b_syn - 0.02}
    choice2 = max(ok2, key=lambda k: ok2[k]["genbank_cds_rescued"]) if ok2 else None
    (out / "long_profile_variants.json").write_text(json.dumps(
        {"variants": result, "chosen_prespecified": choice, "chosen_synonyms": choice2}, indent=2) + "\n")
    print("chosen (prespecified):", choice, "| chosen (synonyms, baseline - 0.02):", choice2)


if __name__ == "__main__":
    main()
