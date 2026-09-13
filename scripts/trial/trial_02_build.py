"""
Trial Phase 1, step 2: build the joined variant table for the trial gene set.

ClinVar -> AlphaMissense (coordinate join) -> assertion -> AlphaFold structure
-> EVE. Same architecture as the Phase 0 spike, scaled and with three changes:

  * de-dup prefers the accession chosen in step 1, then MANE, then lowest ENST.
    Step 1 fixed one accession per gene; a coordinate join can still return
    other isoforms at the same locus, and mixing them puts two position systems
    in one column.
  * the EVE "one model per protein" assert is relaxed to a subset check --
    EVE covers 3,212 genes, so most of the trial set is legitimately absent.
    Absence is recorded as eve_missing, never dropped.
  * adds weighted contact number, which the design review notes is better
    behaved than RSA.

ESM-1v is deliberately NOT run here: 2.571 s/position on this laptop puts the
trial set at multiple hours. It is the Hoffman2 array's job. This table is
therefore AM + EVE; the structural Phase 2 baseline does not need ESM-1v.
"""
import json, urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, duckdb
import biotite.structure.io.pdb as pdb
import biotite.structure as struc

OUT = Path("data/trial")
AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
       "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
       "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
MAXASA = {"A":129,"R":274,"N":195,"D":193,"C":167,"Q":225,"E":223,"G":104,
          "H":224,"I":197,"L":201,"K":236,"M":224,"F":240,"P":159,"S":155,
          "T":172,"W":285,"Y":263,"V":174}   # Tien et al. 2013
TWO_STAR = ("criteria provided, multiple submitters, no conflicts",
            "reviewed by expert panel", "practice guideline")
PATH   = ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")
BENIGN = ("Benign", "Likely benign", "Benign/Likely benign")

genes = pd.read_parquet(OUT/"genes.parquet")
seqs  = json.load(open(OUT/"sequences.json"))
GENES = tuple(genes.gene.tolist())
print(f"[0] {len(genes)} genes, {genes.prot_len.sum():,} residues")

con = duckdb.connect(); con.execute("PRAGMA threads=4")
con.register("target", genes[["gene","uniprot_id"]].rename(columns={"uniprot_id":"target_acc"}))

print("[1] ClinVar >=2-star P/LP + B/LB SNVs ...")
con.execute(f"""
CREATE TABLE cv AS
SELECT "GeneSymbol" AS gene, 'chr'||"Chromosome" AS chrom,
       CAST("PositionVCF" AS BIGINT) AS pos,
       "ReferenceAlleleVCF" AS ref, "AlternateAlleleVCF" AS alt,
       "Name" AS cv_name, "ClinicalSignificance" AS cv_sig,
       CAST("VariationID" AS BIGINT) AS variation_id,
       CASE WHEN "ClinicalSignificance" IN {PATH} THEN 1 ELSE 0 END AS label_path
FROM read_csv('data/clinvar/variant_summary.txt.gz', delim='\t', header=true,
              quote='', all_varchar=true, ignore_errors=true)
WHERE "Assembly"='GRCh38' AND "Type"='single nucleotide variant'
  AND "ReviewStatus" IN {TWO_STAR}
  AND ("ClinicalSignificance" IN {PATH} OR "ClinicalSignificance" IN {BENIGN})
  AND "GeneSymbol" IN {GENES} AND "PositionVCF" IS NOT NULL
""")
n_cv = con.execute("SELECT count(*) FROM cv").fetchone()[0]
print(f"    {n_cv:,} ClinVar rows")

print("[2] MANE (transcript RANKING only -- never an equality join) ...")
con.execute("""CREATE TABLE mane AS
SELECT "symbol" AS gene, split_part("Ensembl_nuc",'.',1) AS enst_base
FROM read_csv('data/mane/mane_summary.txt.gz', delim='\t', header=true,
              all_varchar=true, ignore_errors=true)
WHERE "MANE_status"='MANE Select'""")

print("[3] coordinate join -> AlphaMissense ...")
con.execute("""
CREATE TABLE joined AS
SELECT * EXCLUDE (rn) FROM (
  SELECT cv.*, am.uniprot_id, am.transcript_id, am.protein_variant,
         am.am_pathogenicity, am.am_class,
         CASE WHEN m.enst_base IS NOT NULL THEN 1 ELSE 0 END AS is_mane,
         CASE WHEN am.uniprot_id = t.target_acc THEN 1 ELSE 0 END AS is_target,
         row_number() OVER (PARTITION BY cv.variation_id, cv.chrom, cv.pos, cv.ref, cv.alt
             ORDER BY CASE WHEN am.uniprot_id = t.target_acc THEN 0 ELSE 1 END,
                      CASE WHEN m.enst_base IS NOT NULL THEN 0 ELSE 1 END,
                      am.transcript_id) AS rn
  FROM cv
  JOIN 'data/alphamissense/am_hg38.parquet' am
    ON am."CHROM"=cv.chrom AND am."POS"=cv.pos AND am."REF"=cv.ref AND am."ALT"=cv.alt
  LEFT JOIN mane m ON m.gene=cv.gene AND m.enst_base=split_part(am.transcript_id,'.',1)
  LEFT JOIN target t ON t.gene = cv.gene
) WHERE rn=1
""")
df = con.execute("SELECT * FROM joined").df()
print(f"    {len(df):,} rows | on target accession: {df.is_target.mean()*100:.1f}%"
      f" | MANE transcript: {df.is_mane.mean()*100:.1f}%")

# keep only the accession step 1 chose, so one coordinate system per gene
df = df[df.is_target == 1].reset_index(drop=True)
print(f"    {len(df):,} rows after restricting to the per-gene accession")

print("[4] the reference-amino-acid assertion ...")
am = df.protein_variant.str.extract(r'^([A-Z])(\d+)([A-Z])$')
df["am_ref_aa"], df["aa_pos"], df["am_alt_aa"] = am[0], am[1].astype("Int64"), am[2]
cv = df.cv_name.str.extract(r'\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)')
df["cv_ref_aa"] = cv[0].map(AA3)
df["cv_pos"]    = pd.to_numeric(cv[1], errors="coerce").astype("Int64")
df["cv_alt_aa"] = cv[2].map(AA3)
df["af_ref_aa"] = [ (seqs[u][int(p)-1]
                     if (u in seqs and pd.notna(p) and 1 <= int(p) <= len(seqs[u])) else None)
                    for u, p in zip(df.uniprot_id, df.aa_pos) ]

# THE GATE -- split into a hard part and a soft part, which the 5-gene spike
# could not distinguish because it had no isoform-offset genes in it.
#
# HARD (drop the row): the residue identities must agree. If AlphaMissense's
# reference residue disagrees with the AlphaFold sequence, AM's accession and
# the structure are in different frames and every structural feature attached
# to that row is attached to the wrong residue.
#
# SOFT (keep, flag): am_pos != cv_pos while ALL THREE residues agree. This is
# not corruption -- it is ClinVar's Name field citing a RefSeq isoform whose
# protein numbering differs from the UniProt canonical that AM and AlphaFold
# use. The offset is constant within a gene (MECP2 -12, RUNX1 -27, MEN1 +5,
# MBD5 -233) and the variant identity is already pinned by the (chrom,pos,
# ref,alt) coordinate join. Dropping these would discard whole genes for a
# naming convention. Recorded as cv_frame_offset.
hard = {"AM ref == AlphaFold ref": df.am_ref_aa == df.af_ref_aa,
        "AM ref == ClinVar ref  ": df.am_ref_aa == df.cv_ref_aa,
        "AM alt == ClinVar alt  ": df.am_alt_aa == df.cv_alt_aa}
ok = hard["AM ref == AlphaFold ref"]
for _v in list(hard.values())[1:]:
    ok = ok & _v
ok = ok.fillna(False).astype(bool)
pos_ok = (df.aa_pos == df.cv_pos).fillna(False).astype(bool)
for k, v in hard.items():
    print(f"    [hard] {k}: {v.sum():>6,}/{len(df):,} = {100*v.mean():6.2f}%")
print(f"    [hard] ALL THREE (the gate)  : {ok.sum():>6,}/{len(df):,} = {100*ok.mean():6.2f}%")
print(f"    [soft] AM pos == ClinVar pos : {pos_ok.sum():>6,}/{len(df):,} = {100*pos_ok.mean():6.2f}%")

df["cv_frame_offset"] = (df.aa_pos - df.cv_pos).astype("Float64")
off = df.loc[ok & ~pos_ok]
if len(off):
    per = off.groupby("gene").cv_frame_offset.agg(["count", "median", "nunique"])
    n_const = int((per["nunique"] == 1).sum())
    print(f"    [soft] {len(off):,} rows carry a ClinVar-Name isoform offset, "
          f"{per.shape[0]} genes; offset is constant within gene for "
          f"{n_const}/{per.shape[0]} of them")
    print(per.sort_values("count", ascending=False).head(8).to_string())

if (~ok).sum():
    df[~ok].to_parquet(OUT/"assertion_failures.parquet")
    print(f"    {(~ok).sum()} HARD failures logged -> {OUT/'assertion_failures.parquet'}")
    print(df[~ok].groupby("gene").size().sort_values(ascending=False).head(8).to_string())
df["assert_ok"] = ok
df = df[ok].reset_index(drop=True)   # never model on a biased survivor set
print(f"    proceeding on {len(df):,} rows that pass the hard gate")

print("[5] structure features (pLDDT, SASA, RSA, weighted contact number) ...")
def feats_for(row):
    f = Path("data/alphafold")/row.pdb_file
    arr = pdb.PDBFile.read(str(f)).get_structure(model=1, extra_fields=["b_factor"])
    arr = arr[struc.filter_amino_acids(arr)]
    sasa = struc.sasa(arr, vdw_radii="Single", point_number=100)
    res_ids, inv = np.unique(arr.res_id, return_inverse=True)
    per_res = np.bincount(inv, weights=np.nan_to_num(sasa))
    ca = arr[arr.atom_name == "CA"]
    xyz = ca.coord
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    wcn = (1.0/np.clip(d, 1e-6, None)**2).sum(1) - 1e12*0  # self term excluded below
    np.fill_diagonal(d, np.inf)
    wcn = (1.0/d**2).sum(1)
    return row.uniprot_id, {
        "plddt": dict(zip(ca.res_id.tolist(), ca.b_factor.tolist())),
        "sasa":  dict(zip(res_ids.tolist(), per_res.tolist())),
        "wcn":   dict(zip(ca.res_id.tolist(), wcn.tolist()))}

with ThreadPoolExecutor(6) as ex:
    feats = dict(ex.map(feats_for, [r for _, r in genes.iterrows()]))
print(f"    {len(feats)} structures processed")

for key in ("plddt", "sasa", "wcn"):
    df[key] = [feats[u][key].get(int(p)) if u in feats and pd.notna(p) else None
               for u, p in zip(df.uniprot_id, df.aa_pos)]
df["rsa"] = [ (s/MAXASA[a] if (pd.notna(s) and a in MAXASA) else None)
              for s, a in zip(df.sasa, df.af_ref_aa) ]
df["buried"]   = df.rsa < 0.25
df["prot_len"] = df.uniprot_id.map({u: len(s) for u, s in seqs.items()})
df["rel_pos"]  = df.aa_pos.astype(float) / df.prot_len
df["cropped"]  = df.prot_len > 1022

print("[6] EVE ...")
CACHE = OUT/"uniprot_entry_names.json"
cache = json.load(open(CACHE)) if CACHE.exists() else {}
todo = [a for a in df.uniprot_id.unique() if a not in cache]
def ename(a):
    try:
        with urllib.request.urlopen(
                f"https://rest.uniprot.org/uniprotkb/{a}.json?fields=id", timeout=30) as r:
            return a, json.load(r)["uniProtkbId"]
    except Exception:
        return a, None
if todo:
    with ThreadPoolExecutor(12) as ex:
        for a, n in ex.map(ename, todo):
            if n: cache[a] = n
    json.dump(cache, open(CACHE, "w"), indent=2)
print(f"    {len(cache)} accessions -> entry names")
df["eve_key"] = df.uniprot_id.map(cache)

wanted = tuple(sorted(set(df.eve_key.dropna())))
eve = con.execute(f"""
    SELECT gene, position, wt_aa, mt_aa, EVE_scores_ASM, uncertainty_ASM,
           EVE_classes_75_pct_retained_ASM AS eve_class, BS1, PM2
    FROM 'data/eve/eve_variants.parquet'
    WHERE gene IN {wanted}   -- exact IN: excludes -checkpoint dup and domain-only fits
""").df()
got = set(eve.gene.unique())
assert not any(g.endswith("-checkpoint") for g in got), "checkpoint duplicate leaked in"
assert got <= set(wanted), "EVE returned an unrequested model"
print(f"    {len(eve):,} EVE rows | {len(got)}/{len(wanted)} proteins present in EVE")

n_before = len(df)
df = df.merge(eve, how="left",
              left_on=["eve_key","aa_pos","am_ref_aa","am_alt_aa"],
              right_on=["gene","position","wt_aa","mt_aa"], suffixes=("","_eve"))
assert len(df) == n_before, f"EVE join fanned out: {len(df)} vs {n_before}"
matched = df.EVE_scores_ASM.notna()
assert (df.loc[matched,"wt_aa"] == df.loc[matched,"am_ref_aa"]).all(), \
       "EVE wt residue disagrees with the asserted reference"
df = df.drop(columns=["gene_eve","position","wt_aa","mt_aa"], errors="ignore")
df["eve_missing"] = ~matched
print(f"    EVE coverage {matched.sum():,}/{len(df):,} = {100*matched.mean():.1f}%")

# rank-normalise within protein -- non-negotiable per the design review, the
# raw scales are not comparable (AM is a calibrated probability, EVE is not)
df["am_rank"]  = df.groupby("uniprot_id").am_pathogenicity.rank(pct=True)
df["eve_rank"] = df.groupby("uniprot_id").EVE_scores_ASM.rank(pct=True)

df.to_parquet(OUT/"trial_table.parquet")
print(f"\n[out] {OUT/'trial_table.parquet'}  {len(df):,} rows x {df.shape[1]} cols")
g = df.groupby("gene").agg(n=("variation_id","size"),
                           n_neg=("label_path", lambda s:(s==0).sum()))
print(f"      {df.gene.nunique()} genes | labels {int((df.label_path==1).sum()):,} P/LP "
      f"vs {int((df.label_path==0).sum()):,} B/LB")
print(f"      genes with >=10 negatives: {(g.n_neg>=10).sum()}  "
      f"(>=25: {(g.n_neg>=25).sum()})")
print(f"      unique (protein,position) [ESM-1v cost unit]: "
      f"{df[['uniprot_id','aa_pos']].drop_duplicates().shape[0]:,}")
