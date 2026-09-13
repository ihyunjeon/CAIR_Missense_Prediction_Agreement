"""
Phase 0 spike -- part 2: the reference-amino-acid assertion + structure features.

    assert af_seq[pos-1] == clinvar_ref_aa == am_ref_aa

This is the single most valuable test in the project: it is what catches
transcript-version drift, isoform mismatch and 0/1-indexing errors, all of which
fail quietly. Failures are COUNTED AND LOGGED, never warned-and-skipped.

Also attaches pLDDT and RSA (Tien et al. 2013 theoretical max-ASA) per residue.
"""
import json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
import biotite.structure.io.pdb as pdb
import biotite.structure as struc

OUT = Path("data/spike")
AA3 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
       "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
       "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"}
# Tien et al. 2013, theoretical maximum ASA (A^2)
MAXASA = {"A":129,"R":274,"N":195,"D":193,"C":167,"Q":225,"E":223,"G":104,
          "H":224,"I":197,"L":201,"K":236,"M":224,"F":240,"P":159,"S":155,
          "T":172,"W":285,"Y":263,"V":174}
THREE = {v:k for k,v in AA3.items()}

df = pd.read_parquet(OUT/"joined_pre_assert.parquet")
seqs = json.load(open("data/alphafold/sequences.json"))
print(f"[in] {len(df):,} joined rows, {df.uniprot_id.nunique()} proteins")

# --- parse AM protein_variant  e.g. R273H
am = df.protein_variant.str.extract(r'^([A-Z])(\d+)([A-Z])$')
df["am_ref_aa"], df["aa_pos"], df["am_alt_aa"] = am[0], am[1].astype("Int64"), am[2]

# --- parse ClinVar Name  e.g. NM_000546.6(TP53):c.818G>A (p.Arg273His)
cv = df.cv_name.str.extract(r'\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)')
df["cv_ref_aa"] = cv[0].map(AA3)
df["cv_pos"]    = pd.to_numeric(cv[1], errors="coerce").astype("Int64")
df["cv_alt_aa"] = cv[2].map(AA3)

# --- AlphaFold residue at that position (1-indexed -> 0-indexed)
def af_aa(r):
    s = seqs.get(r.uniprot_id)
    p = r.aa_pos
    if s is None or pd.isna(p) or p < 1 or p > len(s): return None
    return s[int(p)-1]
df["af_ref_aa"] = df.apply(af_aa, axis=1)

# ------------------------------------------------------ THE ASSERTION
ok_am_af = df.am_ref_aa == df.af_ref_aa
ok_am_cv = df.am_ref_aa == df.cv_ref_aa
ok_pos   = df.aa_pos    == df.cv_pos
ok_alt   = df.am_alt_aa == df.cv_alt_aa
ok_all   = ok_am_af & ok_am_cv & ok_pos & ok_alt
n = len(df)

print("\n=== reference-amino-acid assertion ===")
for lbl, m in [("AM ref == AlphaFold ref", ok_am_af),
               ("AM ref == ClinVar ref  ", ok_am_cv),
               ("AM pos == ClinVar pos  ", ok_pos),
               ("AM alt == ClinVar alt  ", ok_alt),
               ("ALL FOUR (the gate)    ", ok_all)]:
    print(f"  {lbl}: {m.sum():>5,}/{n:,} = {100*m.mean():6.2f}%")

fails = df[~ok_all]
if len(fails):
    print(f"\n--- {len(fails)} failing rows, by gene ---")
    print(fails.groupby("gene").size().to_string())
    fails.to_parquet(OUT/"assertion_failures.parquet")
    print(f"  logged -> {OUT/'assertion_failures.parquet'}")
    print(fails[["gene","uniprot_id","cv_name","protein_variant",
                 "am_ref_aa","cv_ref_aa","af_ref_aa","aa_pos","cv_pos"]].head(8).to_string(index=False))

# ------------------------------------------------------ structure features
print("\n=== structure features (pLDDT, RSA) ===")
feats = {}
for u in df.uniprot_id.dropna().unique():
    f = next(Path("data/alphafold").glob(f"AF-{u}-F1-model_v*.pdb"))
    arr = pdb.PDBFile.read(str(f)).get_structure(model=1, extra_fields=["b_factor"])
    arr = arr[struc.filter_amino_acids(arr)]
    sasa = struc.sasa(arr, vdw_radii="Single", point_number=200)
    res_ids = np.unique(arr.res_id)
    per_res_sasa = np.array([np.nansum(sasa[arr.res_id == r]) for r in res_ids])
    ca = arr[arr.atom_name == "CA"]
    plddt = {int(r): float(b) for r, b in zip(ca.res_id, ca.b_factor)}
    feats[u] = {"sasa": dict(zip(res_ids.tolist(), per_res_sasa.tolist())),
                "plddt": plddt}
    print(f"  {u}: {len(res_ids)} residues, mean pLDDT {np.mean(list(plddt.values())):.1f}")

def get_feat(r, key):
    d = feats.get(r.uniprot_id)
    return d[key].get(int(r.aa_pos)) if d and pd.notna(r.aa_pos) else None
df["plddt"] = df.apply(lambda r: get_feat(r, "plddt"), axis=1)
df["sasa"]  = df.apply(lambda r: get_feat(r, "sasa"),  axis=1)
df["rsa"]   = df.apply(lambda r: (r.sasa / MAXASA[r.af_ref_aa])
                       if (pd.notna(r.sasa) and r.af_ref_aa in MAXASA) else None, axis=1)
df["buried"]   = df.rsa < 0.25
df["prot_len"] = df.uniprot_id.map({u: len(s) for u, s in seqs.items()})
df["cropped"]  = df.prot_len > 1022

df["assert_ok"] = ok_all
df.to_parquet(OUT/"joined_full.parquet")
print(f"\n[out] {OUT/'joined_full.parquet'}  ({len(df):,} rows, {df.shape[1]} cols)")
print(f"      passing gate: {ok_all.sum():,}  |  cropped proteins: "
      f"{df.cropped.sum():,} rows ({df[df.cropped].gene.unique().tolist()})")
