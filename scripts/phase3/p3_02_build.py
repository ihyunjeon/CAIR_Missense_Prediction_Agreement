"""
Phase 3, step 2: the two variant tables.

  dms_table.parquet       human ProteinGym DMS singles x {AM, ESM-1v, EVE, ...}
  clinical_table.parquet  ProteinGym clinical set (ClinVar-derived P/B labels)
                          x {AM, ESM-1b, EVE, PrimateAI, ...}

Everything is put in ONE orientation before it leaves this script:
    higher = more damaging (pathogenic / lower fitness)
ProteinGym orients its model scores and DMS_score as "higher = fitter", so those
are negated; AlphaMissense and PrimateAI are already "higher = pathogenic". The
orientation is fixed by each model's published convention, not learned from the
labels, and is then CHECKED (median per-assay Spearman must be positive) -- a
sign error would otherwise invert the whole disagreement analysis silently.

The clinical set has no ESM-1v column. ESM-1b (same architecture family,
UniRef50 rather than UniRef90 training) stands in for it there; how closely the
two agree is measured on the DMS set, where both exist, so the substitution is
quantified rather than assumed.

HARD GATE (per variant, as in the trial): the wild-type residue ProteinGym
names must equal the canonical UniProt residue at the frame-shifted position,
and -- for rows with an AM score -- the AM variant string is built from that
same residue, so a successful AM join IS the second half of the assertion.
"""
import re
from pathlib import Path
import numpy as np, pandas as pd, duckdb

PG = Path("data/proteingym"); OUT = Path("data/phase3")
pm = pd.read_parquet(OUT/"protein_map.parquet")
pm = pm[pm.status == "ok"].copy()
pm["offset"] = pm.offset.astype(int)

# ProteinGym columns kept, and their sign to reach "higher = more damaging"
DMS_MODELS = {"ESM1v_ensemble": -1, "EVE_ensemble": -1, "ESM1b": -1,
              "ESM-IF1": -1, "GEMME": -1, "TranceptEVE_L": -1}
CLIN_MODELS = {"ESM1b": -1, "EVE": -1, "GEMME": -1, "TranceptEVE_L": -1,
               "PoET": +1, "PrimateAI": +1, "REVEL": +1, "CADD": +1}
# PoET's clinical column is NOT "higher = fitter": with the ProteinGym sign it
# scored macro AUROC 0.055, i.e. perfectly inverted. Its orientation in this file
# is undocumented, so +1 is inferred from that check -- PoET is a secondary
# reference here and nothing downstream depends on it. The assertion below now
# covers every column, so a flip cannot pass silently again.
# REVEL and CADD are ensembles trained on (or tuned with) clinical labels and
# allele frequency: reference only, circular against ClinVar.

MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def singles(frame):
    frame = frame[~frame.mutant.astype(str).str.contains(":")].copy()
    p = frame.mutant.str.extract(MUT)
    frame["wt"], frame["pg_pos"], frame["mt"] = p[0], p[1].astype(int), p[2]
    return frame


def gate(frame, row):
    canon = row.canon_seq
    frame["aa_pos"] = frame.pg_pos + row.offset
    inb = frame.aa_pos.between(1, len(canon))
    frame["canon_wt"] = [canon[p-1] if ok else None for p, ok in zip(frame.aa_pos, inb)]
    frame["gate_ok"] = frame.wt == frame.canon_wt
    frame["protein_variant"] = frame.wt + frame.aa_pos.astype(str) + frame.mt
    frame["uniprot_id"] = row.uniprot_id
    return frame


def attach_am(con, df):
    con.register("q", df[["uniprot_id", "protein_variant"]].drop_duplicates())
    am = con.execute("""
        SELECT a.uniprot_id, a.protein_variant, a.am_pathogenicity, a.am_class
        FROM 'data/alphamissense/am_aa_substitutions.parquet' a
        JOIN q USING (uniprot_id, protein_variant)""").df()
    n = len(df)
    df = df.merge(am, how="left", on=["uniprot_id", "protein_variant"])
    assert len(df) == n, "AM join fanned out"
    return df


con = duckdb.connect(); con.execute("PRAGMA threads=6")

# ======================================================================= DMS
print("[1] DMS table ...")
meta = pd.read_csv(PG/"DMS_substitutions.csv").set_index("DMS_id")
parts, audit = [], []
for _, r in pm[pm.set == "dms"].iterrows():
    cols = ["mutant", "DMS_score", "DMS_score_bin"] + list(DMS_MODELS)
    f = pd.read_csv(PG/"dms_scores"/r.pg_file, usecols=lambda c: c in cols)
    f = gate(singles(f), r)
    audit.append({"pg_id": r.pg_id, "n": len(f), "gate_fail": int((~f.gate_ok).sum())})
    f = f[f.gate_ok]
    m = meta.loc[r.pg_id]
    f["assay"] = r.pg_id
    f["selection"] = m.coarse_selection_type
    f["neff_l"] = m.MSA_Neff_L
    f["neff_cat"] = m.MSA_Neff_L_category
    parts.append(f)
dms = pd.concat(parts, ignore_index=True)
a = pd.DataFrame(audit)
print(f"    {len(dms):,} single-variant rows pass the hard gate; "
      f"{a.gate_fail.sum():,} dropped across {(a.gate_fail > 0).sum()} assays "
      f"(engineered construct residues)")
dms = attach_am(con, dms)
print(f"    AM coverage {dms.am_pathogenicity.notna().mean()*100:.2f}%")

dms["dms_damage"] = -dms.DMS_score
for c, s in DMS_MODELS.items():
    dms[c] = s * dms[c]
dms = dms.rename(columns={"ESM1v_ensemble": "esm1v", "EVE_ensemble": "eve",
                          "ESM1b": "esm1b", "ESM-IF1": "esmif1", "GEMME": "gemme",
                          "TranceptEVE_L": "tranceptEVE", "am_pathogenicity": "am"})

# orientation check: every predictor must correlate POSITIVELY with DMS damage
chk = {}
for c in ["am", "esm1v", "eve", "esm1b", "esmif1", "gemme", "tranceptEVE"]:
    rho = dms.groupby("assay").apply(
        lambda g: g[c].corr(g.dms_damage, method="spearman"), include_groups=False)
    chk[c] = rho.median()
print("    median per-assay Spearman vs DMS damage:",
      {k: round(v, 3) for k, v in chk.items()})
assert all(v > 0 for v in chk.values()), "orientation error"
dms.drop(columns=["mutated_sequence"], errors="ignore").to_parquet(OUT/"dms_table.parquet")
print(f"    -> {OUT/'dms_table.parquet'}: {len(dms):,} rows, "
      f"{dms.assay.nunique()} assays, {dms.uniprot_id.nunique()} proteins")

# ======================================================================= clinical
print("[2] clinical table ...")
parts, audit = [], []
for _, r in pm[pm.set == "clinical"].iterrows():
    cols = ["mutant", "DMS_bin_score"] + list(CLIN_MODELS)
    f = pd.read_csv(PG/"clin_scores"/r.pg_file, usecols=lambda c: c in cols)
    f = gate(singles(f), r)
    audit.append({"n": len(f), "gate_fail": int((~f.gate_ok).sum())})
    f = f[f.gate_ok]
    f["refseq"] = r.pg_id
    f["gene"] = r.entry_name.replace("_HUMAN", "")
    parts.append(f)
cl = pd.concat(parts, ignore_index=True)
a = pd.DataFrame(audit)
print(f"    {len(cl):,} rows pass; {a.gate_fail.sum():,} dropped by the hard gate")
cl["label_path"] = cl.DMS_bin_score.map({"Pathogenic": 1, "Benign": 0})
bad = cl.label_path.isna().sum()
assert bad == 0, f"{bad} unmapped labels: {cl.DMS_bin_score.unique()}"
n0 = len(cl)
cl = cl.drop_duplicates(["uniprot_id", "protein_variant"])
print(f"    {n0-len(cl)} duplicate (protein, variant) rows removed")
cl = attach_am(con, cl)
print(f"    AM coverage {cl.am_pathogenicity.notna().mean()*100:.2f}%")
for c, s in CLIN_MODELS.items():
    cl[c] = s * pd.to_numeric(cl[c], errors="coerce")
cl = cl.rename(columns={"ESM1b": "esm1b", "EVE": "eve", "GEMME": "gemme",
                        "TranceptEVE_L": "tranceptEVE", "PoET": "poet",
                        "PrimateAI": "primateai", "REVEL": "revel", "CADD": "cadd",
                        "am_pathogenicity": "am"})
from sklearn.metrics import roc_auc_score
for c in ["am", "esm1b", "eve", "primateai", "gemme", "tranceptEVE", "poet",
          "revel", "cadd"]:
    m = cl[c].notna()
    auc = roc_auc_score(cl.label_path[m], cl[c][m])
    print(f"    pooled AUROC {c:10s} {auc:.3f}  (coverage {m.mean()*100:.1f}%)")
    assert auc > 0.5, f"orientation error in {c}"
cl.drop(columns=["DMS_bin_score"]).to_parquet(OUT/"clinical_table.parquet")
print(f"    -> {OUT/'clinical_table.parquet'}: {len(cl):,} rows, "
      f"{cl.uniprot_id.nunique()} proteins | "
      f"{int(cl.label_path.sum()):,} P vs {int((1-cl.label_path).sum()):,} B")
