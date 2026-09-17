"""
Phase 4 (lite): what predicts the signed disagreement between predictors?

Target: D = r_AM - r_PLM, the within-protein percentile rank of AlphaMissense
minus that of the protein language model (ESM-1b on the clinical set, ESM-1v on
DMS). Positive D = AM thinks the variant is worse than the language model does.
A second target, r_AM - r_EVE, holds the alignment axis fixed (both use MSAs),
isolating structure + frequency supervision.

Feature GROUPS, not single features -- pLDDT, RSA and WCN are collinear by
construction, so per-feature importance would be meaningless (design review):

  structure   pLDDT, RSA, WCN, secondary structure
  disorder    IDR call, UniProt disorder, windowed pLDDT / RSA
  function    ELM motif, UniProt motif, PTM, functional site, DNA/Zn, TM
  chemistry   Grantham, hydropathy/charge change, Pro/Gly/Cys/aromatic flags
  protein     length, >1022 crop flag, relative position
  alignment   EVE MSA depth per residue (NaN where EVE has no model)

Gene-held-out (GroupKFold on protein) gradient boosting. For each group g:
  R2_only(g)       model with g alone
  R2_drop(g)       model with everything except g   (LOCO, refitted)
  unique(g)      = R2_full - R2_drop(g)
  shared(g)      = R2_only(g) - unique(g)   (explained, but also available elsewhere)
That is the variance decomposition with the overlap made explicit, rather than
a Venn diagram that silently assigns shared variance to one group.

Then the disagreement is ADJUDICATED on the clinical labels by context (who is
on the right side of the label when they disagree). Unlike the DMS arm this is
label-circular for AM (frequency-tuned) and is reported as such.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import r2_score
from common import OUT, FIG, RES, CHEM, macro, boot, boot_diff, fmt_ci

GROUPS = {
    "structure": ["plddt", "rsa", "wcn", "sse_a", "sse_b"],
    "disorder": ["idr", "dis_uniprot", "plddt_w11", "rsa_w25"],
    "function": ["in_elm", "in_motif", "ptm", "func_site", "dna_or_zf", "tm"],
    "chemistry": CHEM,
    "protein": ["prot_len", "cropped", "rel_pos"],
    "alignment": ["msa_depth_per_len"],
}
BOOLS = ["idr", "dis_uniprot", "in_elm", "in_motif", "ptm", "func_site", "dna_or_zf",
         "tm", "cropped", "to_pro", "from_gly", "cys_change", "arom_loss"]

lines = []
def log(s=""):
    print(s); lines.append(s)


def prep(df, plm):
    df = df[df.has_structure.fillna(False).astype(bool)].copy()
    df = df[df[["am", plm, "eve"]].notna().all(axis=1)]
    for k in ["am", plm, "eve"]:
        df["r_" + k] = df.groupby("unit")[k].rank(pct=True)
    df["D_plm"] = df.r_am - df["r_" + plm]
    df["D_eve"] = df.r_am - df.r_eve
    df["sse_a"] = (df.sse == "a").astype(float)
    df["sse_b"] = (df.sse == "b").astype(float)
    for b in BOOLS:
        df[b] = df[b].astype(float)
    if "entry_name" not in df:
        pm = pd.read_parquet(OUT/"protein_map.parquet").dropna(subset=["uniprot_id"])
        df["entry_name"] = df.uniprot_id.map(pm.drop_duplicates("uniprot_id")
                                             .set_index("uniprot_id").entry_name)
    dep = pd.read_parquet("data/eve/eve_msa_depth.parquet").set_index("gene").msa_depth_per_len
    df["msa_depth_per_len"] = df.entry_name.map(dep)
    return df.reset_index(drop=True)


def oof_r2(df, feats, target, folds):
    pred = np.full(len(df), np.nan)
    X = df[feats].astype(float).values
    for tr, te in folds:
        h = HistGradientBoostingRegressor(max_iter=250, learning_rate=0.06,
                                          max_leaf_nodes=31, random_state=0)
        h.fit(X[tr], df[target].values[tr])
        pred[te] = h.predict(X[te])
    return r2_score(df[target], pred)


def decompose(df, target, tag):
    folds = list(GroupKFold(5).split(df, groups=df.uniprot_id))
    allf = sum(GROUPS.values(), [])
    full = oof_r2(df, allf, target, folds)
    rows = []
    for g, fs in GROUPS.items():
        only = oof_r2(df, fs, target, folds)
        drop = oof_r2(df, [f for f in allf if f not in fs], target, folds)
        uniq = full - drop
        rows.append({"dataset": tag, "target": target, "group": g, "r2_only": only,
                     "r2_drop": drop, "unique": uniq, "shared": only - uniq, "r2_full": full})
    R = pd.DataFrame(rows)
    log(f"   [{tag}] {target}: gene-held-out R2 (all groups) = {full:.3f}  (n={len(df):,})")
    for _, x in R.iterrows():
        log(f"     {x['group']:10s} alone {x['r2_only']:+.3f} | unique (LOCO) {x['unique']:+.3f} | "
            f"shared {x['shared']:+.3f}")
    return R


# ============================================================ clinical
c = pd.read_parquet(OUT/"clinical_analysis.parquet")
c["unit"] = c.uniprot_id
c = prep(c, "esm1b")
log(f"[clinical] {len(c):,} variants, {c.uniprot_id.nunique()} proteins")
log("\n1. Where does AM rank variants differently from the language model / EVE?")
for ctx, g in c.groupby("context"):
    log(f"   {ctx:16s} n={len(g):>6,}  mean D(AM-ESM1b) {g.D_plm.mean():+.3f}  "
        f"mean D(AM-EVE) {g.D_eve.mean():+.3f}  |D_plm|>=.25: {(g.D_plm.abs()>=.25).mean():.3f}")
for k in ["slim_or_ptm"]:
    g = c[(c.context == "IDR") & c[k].astype(bool)]
    log(f"   IDR & {k}   n={len(g):>6,}  mean D(AM-ESM1b) {g.D_plm.mean():+.3f}  "
        f"mean D(AM-EVE) {g.D_eve.mean():+.3f}")

log("\n2. Variance decomposition of the signed disagreement (gene-held-out)")
dec = [decompose(c, "D_plm", "clinical"), decompose(c, "D_eve", "clinical")]

log("\n3. Adjudication on clinical labels (circular for AM; see DMS arm for the clean test)")
win_rows = []
for pair, col in [("AM vs ESM-1b", "D_plm"), ("AM vs EVE", "D_eve")]:
    for ctx in ["all", "ordered-buried", "ordered-exposed", "IDR"]:
        g = c if ctx == "all" else c[c.context == ctx]
        g = g[g[col].abs() >= 0.25]
        # AM is right when it ranks the pathogenic variant higher / the benign lower
        am_right = np.where(g.label_path == 1, g[col] > 0, g[col] < 0)
        per = pd.Series(am_right, index=g.index).groupby(g.uniprot_id).mean()
        per = per[g.groupby("uniprot_id").size() >= 5]
        lo, hi = boot(per)
        win_rows.append({"pair": pair, "context": ctx, "n_var": len(g), "n_prot": len(per),
                         "am_win": per.mean(), "lo": lo, "hi": hi,
                         "am_win_P": am_right[g.label_path.values == 1].mean(),
                         "am_win_B": am_right[g.label_path.values == 0].mean()})
        log(f"   {pair:12s} {ctx:16s} n={len(g):>6,} AM right {fmt_ci(per.mean(), (lo, hi))} "
            f"(on P: {win_rows[-1]['am_win_P']:.3f}, on B: {win_rows[-1]['am_win_B']:.3f})")
pd.DataFrame(win_rows).to_csv(RES/"clin_disagreement_winrates.csv", index=False)

# ============================================================ DMS
d = pd.read_parquet(OUT/"dms_analysis.parquet")
d["unit"] = d.assay
d = prep(d, "esm1v")
# one assay per protein for the decomposition, so a protein is not counted twice
first = d.groupby("uniprot_id").assay.agg(lambda s: s.value_counts().index[0])
dd = d[d.assay.isin(first.values)].reset_index(drop=True)
log(f"\n[DMS] {len(dd):,} variants, {dd.uniprot_id.nunique()} proteins (largest assay per protein)")
dec.append(decompose(dd, "D_plm", "dms"))
pd.concat(dec).to_csv(RES/"disagreement_decomposition.csv", index=False)

# The decomposition says chemistry, not location, carries the disagreement.
# Which substitution classes, and who is right about them? Adjudicated on DMS
# (no label circularity), AM vs ESM-1v, disagreements |D| >= 0.25.
log("\n4. Disagreement by substitution class (mean signed D; DMS win rate for AM vs ESM-1v)")
def classes(df):
    w, m = df.wt, df.mt
    return {
        "to W/C/M (rare residues)": m.isin(list("WCM")),
        "to G": m == "G",
        "to P (helix/strand breaker)": m == "P",
        "to charged D/E/K/R": m.isin(list("DEKR")),
        "from F/W/Y (aromatic)": w.isin(list("FWY")),
        "from C/M (sulfur)": w.isin(list("CM")),
        "from G": w == "G",
        "everything else": ~(m.isin(list("WCMGPDEKR")) | w.isin(list("FWYCMG"))),
    }


a2p = dict(zip(d.assay, d.uniprot_id))
cls_rows = []
for (lab, mask), cmask in zip(classes(d).items(), classes(c).values()):
    g = d[mask]
    dis = g[g.D_plm.abs() >= 0.25]
    right = (dis.r_am - dis.r_dms_damage).abs() < (dis.r_esm1v - dis.r_dms_damage).abs()
    per = right.groupby(dis.assay).mean()[dis.groupby("assay").size() >= 10]
    per = per.groupby(per.index.map(a2p)).mean()
    lo, hi = boot(per)
    cd = c.loc[cmask, "D_plm"].mean()
    cls_rows.append({"class": lab, "dms_mean_D": g.D_plm.mean(), "clin_mean_D": cd,
                     "n_dis": len(dis), "n_prot": len(per), "am_win": per.mean(),
                     "lo": lo, "hi": hi})
    log(f"   {lab:28s} mean D: DMS {g.D_plm.mean():+.3f} clinical {cd:+.3f} | "
        f"DMS disagreements {len(dis):>6,}: AM closer {fmt_ci(per.mean(), (lo, hi))} n_prot={len(per)}")
pd.DataFrame(cls_rows).to_csv(RES/"disagreement_by_substitution.csv", index=False)

# ------------------------------------------------------------------ figure
D = pd.concat(dec)
fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), dpi=200, sharey=True, sharex=True)
order = list(GROUPS)
for ax, (tag, tgt, ttl) in zip(axes, [("clinical", "D_plm", "Clinical: AM - ESM-1b"),
                                     ("clinical", "D_eve", "Clinical: AM - EVE"),
                                     ("dms", "D_plm", "DMS: AM - ESM-1v")]):
    x = D[(D.dataset == tag) & (D.target == tgt)].set_index("group").loc[order]
    yy = np.arange(len(order))
    ax.barh(yy, x.unique.clip(lower=0), color="#2a78d6", height=0.6, label="unique (LOCO)")
    ax.barh(yy, x.shared.clip(lower=0), left=x.unique.clip(lower=0), color="#9dc2ee",
            height=0.6, label="shared")
    ax.set_yticks(yy); ax.set_yticklabels(order, fontsize=8); ax.invert_yaxis()
    ax.set_title(f"{ttl}\nR² = {x.r2_full.iloc[0]:.3f}", fontsize=8.5)
    ax.tick_params(axis="x", labelsize=7)
    ax.grid(axis="x", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
axes[0].legend(frameon=False, fontsize=7, loc="lower right")
axes[1].set_xlabel("gene-held-out R² of signed disagreement", fontsize=8)
fig.tight_layout(); fig.savefig(FIG/"fig_disagreement_decomposition.png"); plt.close(fig)

open(RES/"disagreement_summary.txt", "w").write("\n".join(lines) + "\n")
print(f"[out] {RES}/disagreement_*.csv")
