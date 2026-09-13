"""
Trial Phase 1, step 1: choose the gene set and resolve structures.

Scaled-down stand-in for the full 300-500 gene Phase 1. Selection is by ClinVar
benign count, because benign variants are the scarce class and per-gene AUROC is
meaningless without them (the Phase 0 spike's PTEN 0.997 rested on 6 negatives).

TWO FILTERS THAT MATTER:

1. prot_len <= 2700. AlphaFold DB splits longer proteins into overlapping
   fragments F1, F2, ... A glob for `AF-<acc>-F1-model_v*.pdb` therefore returns
   residues 1-1400 of TTN (34,350 aa) and silently drops everything after, with
   no error. Same silent-truncation failure class as the MANE transcript trap.
   Excluding long proteins is the honest trial-scale fix; the full run must
   handle fragments properly.

2. One UniProt accession per gene, taken as the modal accession over that gene's
   AlphaMissense rows. AM reports several isoforms per locus; picking the mode
   keeps one coordinate system per gene.

AlphaFold version is resolved per accession from the API, never hardcoded --
the _v4 URLs now 404 and the release is on v6.
"""
import json, sys, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request, urllib.error
import duckdb, pandas as pd

N_GENES   = int(sys.argv[1]) if len(sys.argv) > 1 else 60
MAX_LEN   = 2700
CANDIDATE = N_GENES * 3          # over-fetch, filters will drop some
OUT = Path("data/trial"); OUT.mkdir(parents=True, exist_ok=True)
AFDIR = Path("data/alphafold"); AFDIR.mkdir(parents=True, exist_ok=True)

TWO_STAR = ("criteria provided, multiple submitters, no conflicts",
            "reviewed by expert panel", "practice guideline")
PATH   = ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")
BENIGN = ("Benign", "Likely benign", "Benign/Likely benign")

con = duckdb.connect(); con.execute("PRAGMA threads=4")

print(f"[1] ClinVar census -> top {CANDIDATE} genes by benign count ...")
con.execute(f"""
CREATE TABLE cand AS
SELECT "GeneSymbol" AS gene,
       count(*) AS n_cv,
       sum(CASE WHEN "ClinicalSignificance" IN {BENIGN} THEN 1 ELSE 0 END) AS n_neg_cv
FROM read_csv('data/clinvar/variant_summary.txt.gz', delim='\t', header=true,
              quote='', all_varchar=true, ignore_errors=true)
WHERE "Assembly"='GRCh38' AND "Type"='single nucleotide variant'
  AND "ReviewStatus" IN {TWO_STAR}
  AND ("ClinicalSignificance" IN {PATH} OR "ClinicalSignificance" IN {BENIGN})
  AND "PositionVCF" IS NOT NULL
  AND "Name" LIKE '%(p.%' AND "GeneSymbol" NOT LIKE '%;%'
GROUP BY 1
HAVING sum(CASE WHEN "ClinicalSignificance" IN {PATH} THEN 1 ELSE 0 END) >= 10
   AND sum(CASE WHEN "ClinicalSignificance" IN {BENIGN} THEN 1 ELSE 0 END) >= 10
ORDER BY n_neg_cv DESC LIMIT {CANDIDATE}
""")
cand = con.execute("SELECT * FROM cand").df()
print(f"    {len(cand)} candidates (both classes >= 10 in ClinVar)")

print("[2] modal UniProt accession per gene, from AlphaMissense ...")
# AM has no gene symbol, so bridge through ClinVar coordinates.
con.execute(f"""
CREATE TABLE gene_acc AS
WITH cv AS (
  SELECT "GeneSymbol" AS gene, 'chr'||"Chromosome" AS chrom,
         CAST("PositionVCF" AS BIGINT) AS pos,
         "ReferenceAlleleVCF" AS ref, "AlternateAlleleVCF" AS alt
  FROM read_csv('data/clinvar/variant_summary.txt.gz', delim='\t', header=true,
                quote='', all_varchar=true, ignore_errors=true)
  WHERE "Assembly"='GRCh38' AND "Type"='single nucleotide variant'
    AND "GeneSymbol" IN (SELECT gene FROM cand) AND "PositionVCF" IS NOT NULL
), hit AS (
  SELECT cv.gene, am.uniprot_id, count(*) AS n
  FROM cv JOIN 'data/alphamissense/am_hg38.parquet' am
    ON am."CHROM"=cv.chrom AND am."POS"=cv.pos AND am."REF"=cv.ref AND am."ALT"=cv.alt
  GROUP BY 1,2
)
SELECT gene, uniprot_id, n FROM (
  SELECT *, row_number() OVER (PARTITION BY gene ORDER BY n DESC) rn FROM hit
) WHERE rn=1
""")
acc = con.execute("SELECT * FROM gene_acc").df()
acc = acc.merge(cand, on="gene").sort_values("n_neg_cv", ascending=False)
print(f"    resolved {len(acc)} gene -> accession")

print("[3] AlphaFold API: resolve version + length (skip the PDB for now) ...")
def meta(a):
    try:
        with urllib.request.urlopen(
                f"https://alphafold.ebi.ac.uk/api/prediction/{a}", timeout=30) as r:
            j = json.load(r)
        if not j: return None
        # TRAP: /prediction/<acc> may return ISOFORM entries (P21359-5, -4, -3)
        # and never the canonical one. j[0] is then a 593 aa isoform of a 2,839 aa
        # protein -- structure and ClinVar/AM positions silently disagree.
        # Accept only an entry whose accession is the bare one, no -N suffix.
        e = next((x for x in j if x.get("uniprotAccession") == a), None)
        if e is None:
            return {"uniprot_id": a, "err": "isoform-only (no canonical entry)"}
        return {"uniprot_id": a, "pdbUrl": e["pdbUrl"],
                "seq": e["uniprotSequence"], "n_frag": len(j)}
    except Exception as ex:
        return {"uniprot_id": a, "err": str(ex)}

t0 = time.time()
with ThreadPoolExecutor(16) as ex:
    metas = [m for m in ex.map(meta, acc.uniprot_id.tolist()) if m]
md = pd.DataFrame(metas)
print(f"    {len(md)} API responses in {time.time()-t0:.1f}s")
if "err" in md.columns:
    bad = md[md.err.notna()]
    print(f"    {len(bad)} unusable: "
          f"{(bad.err.str.contains('isoform')).sum()} isoform-only, "
          f"{len(bad) - (bad.err.str.contains('isoform')).sum()} API error/404 "
          f"(404 = no AlphaFold prediction; TTN and other very long proteins)")
    md = md[md.err.isna()]
md["prot_len"] = md.seq.str.len()

sel = acc.merge(md, on="uniprot_id")
too_long = (sel.prot_len > MAX_LEN).sum()
sel = sel[sel.prot_len <= MAX_LEN].head(N_GENES).reset_index(drop=True)
print(f"    dropped {too_long} proteins > {MAX_LEN} aa (AlphaFold fragments)")
print(f"    SELECTED {len(sel)} genes | median len {sel.prot_len.median():.0f} "
      f"| max {sel.prot_len.max()} | total residues {sel.prot_len.sum():,}")

print("[4] downloading PDBs ...")
def fetch(r):
    dst = AFDIR / r.pdbUrl.rsplit("/", 1)[-1]
    if dst.exists(): return 0
    try:
        urllib.request.urlretrieve(r.pdbUrl, dst); return 1
    except Exception: return -1
t0 = time.time()
with ThreadPoolExecutor(12) as ex:
    res = list(ex.map(fetch, [r for _, r in sel.iterrows()]))
print(f"    {sum(1 for x in res if x==1)} downloaded, "
      f"{sum(1 for x in res if x==0)} cached, {sum(1 for x in res if x==-1)} failed "
      f"in {time.time()-t0:.1f}s")

sel["pdb_file"] = [r.pdbUrl.rsplit("/",1)[-1] for _, r in sel.iterrows()]
sel.drop(columns=["seq"]).to_parquet(OUT/"genes.parquet")
json.dump(dict(zip(sel.uniprot_id, sel.seq)), open(OUT/"sequences.json","w"))
print(f"\n[out] {OUT/'genes.parquet'}  ({len(sel)} genes)")
print(sel[["gene","uniprot_id","prot_len","n_cv","n_neg_cv"]].head(12).to_string(index=False))
