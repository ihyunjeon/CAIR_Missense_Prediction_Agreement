"""
Phase 0 spike -- part 4: attach EVE scores.

EVE is the predictor that breaks the design's three-way confound: it is
MSA-based with NO structure and NO frequency supervision, so AM - EVE isolates
(structure + frequency supervision) with MSA held constant, and EVE - ESM
isolates the alignment axis.

THREE JOIN HAZARDS, all of which return rows rather than erroring:

1. EVE is keyed by UniProt ENTRY NAME (P53_HUMAN), not accession (P04637) and
   not gene symbol. TP53's entry name drops the T, so matching on the gene
   symbol silently loses the gene. We resolve accession -> entry name from the
   UniProt REST API, never by string manipulation.

2. `BRCA1_HUMAN-checkpoint` is a duplicate of `BRCA1_HUMAN` -- 37,260 rows,
   verified byte-identical. A training artifact in the release. Dropped.

3. BRCA1 additionally ships domain-only models (BRCA1_RING_HUMAN,
   BRCA1_BRCT_HUMAN) alongside the full-length fit. It is the only gene in the
   release with domain splits. We take the FULL-LENGTH model, because position
   numbering in the domain fits is not guaranteed to match the full protein and
   mixing them would put two different coordinate systems in one column.
"""
import json, urllib.request, urllib.parse
from pathlib import Path
import duckdb, pandas as pd

OUT = Path("data/spike")
EVE = "data/eve/eve_variants.parquet"
CACHE = OUT / "uniprot_entry_names.json"


def entry_names(accessions):
    """accession -> UniProt entry name, via the REST API (cached)."""
    cache = json.load(open(CACHE)) if CACHE.exists() else {}
    todo = [a for a in accessions if a not in cache]
    for acc in todo:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json?fields=id"
        with urllib.request.urlopen(url) as r:
            cache[acc] = json.load(r)["uniProtkbId"]
        print(f"    {acc} -> {cache[acc]}")
    if todo:
        json.dump(cache, open(CACHE, "w"), indent=2)
    return cache


def main():
    var = pd.read_parquet(OUT / "spike_final.parquet")
    accs = sorted(var.uniprot_id.dropna().unique())
    print(f"[1] resolving {len(accs)} accessions to UniProt entry names")
    names = entry_names(accs)
    var["eve_key"] = var.uniprot_id.map(names)

    con = duckdb.connect(); con.execute("PRAGMA threads=4")
    wanted = tuple(sorted(set(names.values())))
    print(f"[2] loading EVE rows for {wanted}")
    eve = con.execute(f"""
        SELECT gene, position, wt_aa, mt_aa,
               EVE_scores_ASM, uncertainty_ASM,
               EVE_classes_75_pct_retained_ASM AS eve_class,
               BS1, PM2
        FROM '{EVE}'
        WHERE gene IN {wanted}          -- exact match drops -checkpoint and
                                        -- the domain-only models by construction
    """).df()
    print(f"    {len(eve):,} EVE rows")

    # hazard 2/3 guard: assert we pulled exactly one model per protein
    got = set(eve.gene.unique())
    assert got == set(wanted), f"EVE model mismatch: {got ^ set(wanted)}"
    assert not any(g.endswith("-checkpoint") for g in got), "checkpoint dup leaked in"

    print("[3] joining on (entry_name, position, wt_aa, mt_aa)")
    m = var.merge(eve, how="left",
                  left_on=["eve_key", "aa_pos", "am_ref_aa", "am_alt_aa"],
                  right_on=["gene", "position", "wt_aa", "mt_aa"],
                  suffixes=("", "_eve"))
    assert len(m) == len(var), f"join fanned out: {len(m)} vs {len(var)}"

    matched = m.EVE_scores_ASM.notna()
    ok_wt = (m.loc[matched, "wt_aa"] == m.loc[matched, "am_ref_aa"]).all()
    print(f"    matched {matched.sum():,}/{len(m):,} ({100*matched.mean():.1f}%)"
          f" | wt_aa consistent: {ok_wt}")
    assert ok_wt, "EVE wt residue disagrees with the asserted reference"

    m = m.drop(columns=["gene_eve", "position", "wt_aa", "mt_aa"], errors="ignore")
    m["eve_missing"] = ~matched          # missing-not-at-random; never drop silently

    # rank-normalise within protein, same treatment as AM and ESM
    m["eve_rank"] = m.groupby("uniprot_id").EVE_scores_ASM.rank(pct=True)

    m.to_parquet(OUT / "spike_final_3pred.parquet")
    print(f"\n[out] {OUT/'spike_final_3pred.parquet'}  {len(m):,} rows x {m.shape[1]} cols")

    print("\n=== per-gene EVE coverage ===")
    cov = m.groupby("gene").agg(spike_variants=("variation_id", "size"),
                                eve_scored=("EVE_scores_ASM", "count"))
    cov["pct"] = (100 * cov.eve_scored / cov.spike_variants).round(1)
    print(cov.to_string())
    return m


if __name__ == "__main__":
    main()
