"""
Phase 0 spike -- part 3: attach ESM-1v ensemble scores and emit the final
end-to-end table. One row = one ClinVar variant with label, AlphaMissense score,
ESM-1v ensemble LLR, and AlphaFold structural context.

LLR is computed at analysis time as lp[mut] - lp[wt], per the design in
scripts/esm1v_masked_marginals.py -- per-model rows are kept so ensemble spread
survives as an uncertainty feature.
"""
from pathlib import Path
import numpy as np, pandas as pd

OUT = Path("data/spike")
var = pd.read_parquet(OUT / "joined_full.parquet")
sc  = pd.read_parquet(OUT / "esm1v_scores.parquet")
print(f"[in] {len(var):,} variants | {len(sc):,} score rows "
      f"({sc.model_idx.nunique()} models x {len(sc)//sc.model_idx.nunique():,} positions)")

# long -> per (protein,pos,model) LLR for the specific alt allele
sc = sc.rename(columns={"pos": "aa_pos"})
m = var.merge(sc, on=["uniprot_id", "aa_pos"], how="left",
              suffixes=("", "_esm"))
print(f"[join] {len(m):,} variant x model rows; "
      f"unmatched: {m.model_idx.isna().sum():,}")

# sanity: the wt residue ESM saw must equal the one we asserted on
wt_ok = (m.wt_aa == m.am_ref_aa)
print(f"[check] ESM wt_aa == AM ref_aa : {wt_ok.sum():,}/{len(m):,} "
      f"= {100*wt_ok.mean():.2f}%")

def llr(r):
    return r[f"lp_{r.am_alt_aa}"] - r[f"lp_{r.am_ref_aa}"]
m["esm_llr"] = m.apply(llr, axis=1)

ens = (m.groupby(["variation_id", "uniprot_id", "aa_pos", "am_alt_aa"])
         .agg(esm_llr_mean=("esm_llr", "mean"),
              esm_llr_sd=("esm_llr", "std"),
              n_models=("model_idx", "nunique")).reset_index())

final = var.merge(ens, on=["variation_id", "uniprot_id", "aa_pos", "am_alt_aa"],
                  how="left")

# within-protein rank normalisation -- the plan calls this non-negotiable,
# because AM is a calibrated probability and ESM is an unbounded LLR.
final["am_rank"]  = final.groupby("uniprot_id").am_pathogenicity.rank(pct=True)
final["esm_rank"] = final.groupby("uniprot_id").esm_llr_mean.rank(pct=True)
# ESM: more negative = more damaging, so flip to align with AM's direction
final["esm_rank"] = 1 - final["esm_rank"]
final["disagreement"] = final.esm_rank - final.am_rank   # signed; sign matters

final.to_parquet(OUT / "spike_final.parquet")
print(f"\n[out] {OUT/'spike_final.parquet'}  {len(final):,} rows x {final.shape[1]} cols")

cols = ["gene","cv_name","cv_sig","label_path","protein_variant","aa_pos",
        "am_pathogenicity","am_class","esm_llr_mean","esm_llr_sd",
        "am_rank","esm_rank","disagreement","plddt","rsa","buried","cropped"]
print("\n=== ONE FULLY JOINED ROW, END TO END ===")
r = final[final.gene == "TP53"].dropna(subset=["esm_llr_mean"]).iloc[0]
for c in cols:
    print(f"  {c:18s} {r[c]}")

print("\n=== coverage ===")
print(f"  variants with ESM score : {final.esm_llr_mean.notna().sum():,}/{len(final):,}")
print(f"  mean ensemble sd        : {final.esm_llr_sd.mean():.2f}")
print("\n=== does the pipeline recover known signal? ===")
for name, g in [("AM pathogenicity", "am_pathogenicity"), ("ESM LLR", "esm_llr_mean")]:
    p = final[final.label_path == 1][g].mean()
    b = final[final.label_path == 0][g].mean()
    print(f"  {name:18s} P/LP mean {p:+.3f}   B/LB mean {b:+.3f}   gap {p-b:+.3f}")
