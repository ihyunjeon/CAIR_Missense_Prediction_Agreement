"""
Clinical arm at scale (2,344 proteins) + the disordered-region analysis.

The biochemical question: missense variants in intrinsically disordered regions
(IDRs) are mostly benign because there is no fold to destabilise -- but IDRs
carry short linear motifs and PTM switches where a single substitution can be
pathogenic. Does each predictor (a) know IDRs are mostly tolerant, and (b) still
catch the functional exceptions?

A. Macro AUROC per predictor (per protein; proteins with >=10 of each class).
B. Base rates: P fraction in IDR vs ordered, pooled AND within protein
   (Mantel-Haenszel OR stratified by protein, so gene mix cannot fake it).
   The exception class: slim_or_ptm inside IDRs.
C. Discrimination by context, and behaviour at published operating points:
   AM >= 0.564 (its 'likely pathogenic' class); ESM-1b LLR <= -7.5 (Brandes et
   al. 2023 -- background-knowledge); EVE at the cutoff that gives the same
   overall call rate as AM (label-free matching, since ProteinGym's EVE column
   is not the calibrated EVE class score).
D. Gene-held-out models: does adding biochemical context (disorder, motifs,
   substitution chemistry) to the predictors improve per-protein ranking?
   If the predictors had fully internalised "IDR = tolerant", context adds
   nothing within protein.
E. Label circularity: ACMG BS1 (benign by allele frequency) from the EVE
   release. AlphaMissense is tuned on population frequency; if its lead over
   ESM-1b is partly frequency-label leakage, the lead shrinks when the benign
   class is restricted to non-BS1 variants.
"""
import numpy as np, pandas as pd, duckdb
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import HistGradientBoostingClassifier
from common import OUT, FIG, RES, context, chemistry, CHEM, macro, boot, boot_diff, fmt_ci

PRED = ["am", "esm1b", "eve"]
EXTRA = ["primateai", "gemme", "tranceptEVE", "poet"]
CIRC = ["revel", "cadd"]
NAMES = {"am": "AlphaMissense", "esm1b": "ESM-1b", "eve": "EVE", "primateai": "PrimateAI",
         "gemme": "GEMME", "tranceptEVE": "TranceptEVE-L", "poet": "PoET",
         "revel": "REVEL (circular)", "cadd": "CADD (circular)"}
COL = {"am": "#2a78d6", "esm1b": "#eb6834", "eve": "#1baf7a"}
MIN_CLS = 10

c = pd.read_parquet(OUT/"clinical_table.parquet")
_pm = pd.read_parquet(OUT/"protein_map.parquet").dropna(subset=["uniprot_id"])
c["entry_name"] = c.uniprot_id.map(_pm.drop_duplicates("uniprot_id").set_index("uniprot_id").entry_name)
r = pd.read_parquet(OUT/"residues.parquet")
c = c.merge(r, how="left", on=["uniprot_id", "aa_pos"])
c = context(chemistry(c))
y = c.label_path
cnt = c.groupby("uniprot_id").label_path.agg(n="size", p="sum")
cnt["b"] = cnt.n - cnt.p
EVAL = set(cnt.index[(cnt.p >= MIN_CLS) & (cnt.b >= MIN_CLS)])

lines = []
def log(s=""):
    print(s); lines.append(s)

log(f"[in] {len(c):,} variants, {c.uniprot_id.nunique()} proteins; "
    f"{len(EVAL)} proteins with >={MIN_CLS} of each class ({c.uniprot_id.isin(EVAL).sum():,} variants)")
log(f"     context: {c.context.value_counts(dropna=False).to_dict()}")


def per_prot_auc(frame, col, min_cls=MIN_CLS, prots=None):
    out = {}
    f = frame if prots is None else frame[frame.uniprot_id.isin(prots)]
    for u, g in f.groupby("uniprot_id"):
        g = g[g[col].notna()]
        p = g.label_path.sum(); b = len(g) - p
        if p >= min_cls and b >= min_cls:
            out[u] = roc_auc_score(g.label_path, g[col])
    return pd.Series(out, dtype=float)


# ------------------------------------------------------------------ A
log("\nA. Macro AUROC (per protein, >=10 of each class), 95% CI over proteins")
aucs, rows = {}, []
for k in PRED + EXTRA + CIRC:
    pp = per_prot_auc(c, k)
    aucs[k] = pp
    lo, hi = boot(pp)
    pooled = roc_auc_score(y[c[k].notna()], c.loc[c[k].notna(), k])
    rows.append({"predictor": k, "n_prot": len(pp), "macro": pp.mean(), "lo": lo,
                 "hi": hi, "pooled": pooled})
    log(f"   {NAMES[k]:18s} macro {fmt_ci(pp.mean(), (lo, hi))}  pooled {pooled:.3f}  n={len(pp)}")
for x, z in [("am", "esm1b"), ("am", "eve"), ("esm1b", "eve")]:
    m, ci, p = boot_diff(aucs[x], aucs[z])
    log(f"   {NAMES[x]} - {NAMES[z]}: {fmt_ci(m, ci)} p_boot={p:.3f}")
pd.DataFrame(rows).to_csv(RES/"clin_auroc.csv", index=False)

# ------------------------------------------------------------------ B
log("\nB. Base rates (P = pathogenic fraction)")
def mh_or(frame, expo):
    num = den = 0.0
    for _, g in frame.groupby("uniprot_id"):
        e = g[expo].astype(bool); p = g.label_path == 1
        n = len(g)
        if e.all() or (~e).all(): continue
        a, b = (e & p).sum(), (e & ~p).sum()
        cc, dd = (~e & p).sum(), (~e & ~p).sum()
        num += a * dd / n; den += b * cc / n
    return num / den if den else np.nan

s = c[c.context.notna()].copy()
s["is_idr"] = s.context == "IDR"
for lab, sub in [("IDR", s[s.is_idr]), ("ordered-exposed", s[s.context == "ordered-exposed"]),
                 ("ordered-buried", s[s.context == "ordered-buried"]),
                 ("collagen-helix", s[s.context == "collagen-helix"])]:
    log(f"   {lab:16s} n={len(sub):>6,}  P fraction {sub.label_path.mean():.3f}")
s = s[s.context != "collagen-helix"]      # compare IDR with genuinely ordered residues
log(f"   IDR vs ordered, within-protein Mantel-Haenszel OR for pathogenic: {mh_or(s, 'is_idr'):.3f}")
_all = c[c.context.notna()]
_naive = _all.dis_plddt.astype(bool) & _all.dis_rsa.astype(bool)
log(f"   (naive structure-proxy IDR call, collagen included: n={_naive.sum():,}, "
    f"P fraction {_all[_naive].label_path.mean():.3f}; Gly-substitution P fraction "
    f"{_all[_naive & (_all.wt == 'G')].label_path.mean():.3f})")
s["uniprot_dis"] = s.dis_uniprot.astype(bool)
log(f"   (orthogonal definition, UniProt/MobiDB 'Disordered': OR {mh_or(s, 'uniprot_dis'):.3f})")
idr = s[s.is_idr]
for k in ["slim_or_ptm", "in_elm", "in_motif", "ptm", "func_site"]:
    sub = idr[idr[k]]
    log(f"   IDR & {k:11s} n={len(sub):>5,}  P fraction {sub.label_path.mean():.3f} "
        f"vs {idr[~idr[k]].label_path.mean():.3f} otherwise; within-protein OR "
        f"{mh_or(idr, k):.2f}")
tab = (s.assign(ctx=np.where(s.is_idr, "IDR", "ordered"))
        .groupby(["ctx", "slim_or_ptm"]).label_path.agg(["size", "mean"]))
tab.to_csv(RES/"clin_base_rates.csv")

# chemistry that matters in IDRs (within-protein logistic-free check: P fraction)
log("   IDR substitution chemistry (P fraction):")
for k in ["to_pro", "from_gly", "cys_change", "arom_loss"]:
    log(f"     {k:11s} yes {idr[idr[k]].label_path.mean():.3f} (n={idr[k].sum():,}) "
        f"no {idr[~idr[k]].label_path.mean():.3f}")
idr_ch = idr.assign(charge_change=idr.d_charge != 0)
log(f"     charge chg  yes {idr_ch[idr_ch.charge_change].label_path.mean():.3f} "
    f"no {idr_ch[~idr_ch.charge_change].label_path.mean():.3f}")

# ------------------------------------------------------------------ C
log("\nC. Discrimination by context (per-protein AUROC, >=5 of each class in stratum)")
ctx_rows = []
for ctx in ["ordered-buried", "ordered-exposed", "IDR"]:
    sub = s[s.context == ctx]
    out = []
    for k in PRED + ["primateai"]:
        pp = per_prot_auc(sub, k, min_cls=5)
        pooled = roc_auc_score(sub.label_path[sub[k].notna()], sub.loc[sub[k].notna(), k])
        lo, hi = boot(pp)
        ctx_rows.append({"context": ctx, "predictor": k, "n_prot": len(pp),
                         "macro": pp.mean(), "lo": lo, "hi": hi, "pooled": pooled})
        out.append(f"{NAMES[k]} {pp.mean():.3f} (pooled {pooled:.3f})")
    log(f"   {ctx:16s} n_prot={len(pp):3d} | " + " | ".join(out))
# the motif exceptions: can a predictor separate P from B *inside* IDR motifs?
for flag in [True, False]:
    sub = idr[idr.slim_or_ptm == flag]
    out = []
    for k in PRED:
        m = sub[k].notna()
        out.append(f"{NAMES[k]} {roc_auc_score(sub.label_path[m], sub.loc[m, k]):.3f}")
    log(f"   IDR, slim_or_ptm={flag!s:5s} n={len(sub):,} (P {sub.label_path.mean():.2f}) "
        f"pooled AUROC | " + " | ".join(out))
pd.DataFrame(ctx_rows).to_csv(RES/"clin_auroc_by_context.csv", index=False)

log("\n   Operating points (calls = predicted pathogenic):")
c["call_am"] = c.am >= 0.564
c["call_esm1b"] = c.esm1b >= 7.5
rate = c.call_am[c.am.notna()].mean()
thr_eve = c.eve.quantile(1 - rate)
c["call_eve"] = c.eve >= thr_eve
log(f"   AM call rate {rate:.3f}; EVE matched cutoff {thr_eve:.3f}; "
    f"ESM-1b call rate {c.call_esm1b.mean():.3f}")
op_rows = []
s = c[c.context.notna()]
for ctx in ["ordered-buried", "ordered-exposed", "IDR", "IDR+motif/PTM"]:
    sub = (s[(s.context == "IDR") & s.slim_or_ptm] if ctx == "IDR+motif/PTM"
           else s[s.context == ctx])
    for k in PRED:
        v = sub[sub[k].notna()]
        fpr = v.loc[v.label_path == 0, f"call_{k}"].mean()
        tpr = v.loc[v.label_path == 1, f"call_{k}"].mean()
        op_rows.append({"context": ctx, "predictor": k, "fpr": fpr, "tpr": tpr,
                        "n_b": int((v.label_path == 0).sum()), "n_p": int(v.label_path.sum())})
O = pd.DataFrame(op_rows); O.to_csv(RES/"clin_operating_points.csv", index=False)
for ctx, g in O.groupby("context", sort=False):
    log(f"   {ctx:16s} nB={g.n_b.iloc[0]:>6,} nP={g.n_p.iloc[0]:>6,} | " + " | ".join(
        f"{NAMES[k]} FPR {f:.3f} TPR {t:.3f}" for k, f, t in zip(g.predictor, g.fpr, g.tpr)))

# ------------------------------------------------------------------ D
log("\nD. Gene-held-out models (5-fold GroupKFold on protein); macro AUROC on the same proteins")
STRUCT = ["plddt", "rsa", "wcn", "plddt_w11", "rsa_w25", "rel_pos", "prot_len"]
CTXF = STRUCT + ["idr", "dis_uniprot", "in_elm", "in_motif", "ptm", "func_site",
                 "dna_or_zf", "tm"] + CHEM
m = c[c.has_structure.fillna(False).astype(bool) & c[PRED].notna().all(axis=1)].copy()
m = m.reset_index(drop=True)
for k in ["idr", "dis_uniprot", "in_elm", "in_motif", "ptm", "func_site", "dna_or_zf",
          "tm", "to_pro", "from_gly", "cys_change", "arom_loss"]:
    m[k] = m[k].astype(float)
folds = list(GroupKFold(5).split(m, m.label_path, m.uniprot_id))
SPECS = {"context_only (biochemical prior)": CTXF,
         "AM + context": ["am"] + CTXF,
         "AM + ESM-1b + EVE": PRED,
         "AM + ESM-1b + EVE + context": PRED + CTXF}
for name, feats in SPECS.items():
    oof = np.full(len(m), np.nan)
    for tr, te in folds:
        h = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06,
                                           max_leaf_nodes=31, random_state=0)
        h.fit(m.loc[tr, feats].astype(float), m.label_path[tr])
        oof[te] = h.predict_proba(m.loc[te, feats].astype(float))[:, 1]
    m[name] = oof
model_rows = []
base = {k: per_prot_auc(m, k) for k in PRED}
for name in list(SPECS) + PRED:
    pp = per_prot_auc(m, name)
    lo, hi = boot(pp)
    dm, dci, dp = boot_diff(pp, base["am"])
    row = {"model": name, "n_prot": len(pp), "macro": pp.mean(), "lo": lo, "hi": hi,
           "delta_vs_am": dm, "d_lo": dci[0], "d_hi": dci[1], "p_boot": dp}
    for ctx in ["ordered-buried", "ordered-exposed", "IDR"]:
        row[f"macro_{ctx}"] = per_prot_auc(m[m.context == ctx], name, min_cls=5).mean()
    model_rows.append(row)
    log(f"   {name:34s} {fmt_ci(pp.mean(), (lo, hi))}  vs AM {dm:+.3f} "
        f"[{dci[0]:+.3f}, {dci[1]:+.3f}] | buried {row['macro_ordered-buried']:.3f} "
        f"exposed {row['macro_ordered-exposed']:.3f} IDR {row['macro_IDR']:.3f}")
pd.DataFrame(model_rows).to_csv(RES/"clin_models.csv", index=False)
m[["uniprot_id", "protein_variant", "label_path", "context"] + list(SPECS)].to_parquet(
    OUT/"clin_model_oof.parquet")

# ------------------------------------------------------------------ E
log("\nE. Benign-label circularity: ACMG BS1 (allele-frequency benign) from the EVE release")
con = duckdb.connect()
keys = tuple(sorted(set(c.entry_name.dropna())))
eve = con.execute(f"""SELECT gene AS entry_name, position AS aa_pos, wt_aa AS wt,
    mt_aa AS mt, BS1, frequency_gv3 FROM 'data/eve/eve_variants.parquet'
    WHERE gene IN {keys}""").df()
e = c.merge(eve, how="inner", on=["entry_name", "aa_pos", "wt", "mt"])
log(f"   {len(e):,} variants over {e.uniprot_id.nunique()} proteins carry an EVE-release BS1 flag")
e["bs1"] = e.BS1.astype(str).str.lower().isin(["true", "1", "yes"])
ben = e[e.label_path == 0]
log(f"   benign labels that meet BS1: {ben.bs1.mean():.3f} overall; "
    f"IDR {ben[ben.context=='IDR'].bs1.mean():.3f} vs ordered "
    f"{ben[ben.context.isin(['ordered-buried','ordered-exposed'])].bs1.mean():.3f}")
circ_rows = []
for lab, keep in [("benign = BS1 only", True), ("benign = non-BS1 only", False)]:
    sub = e[(e.label_path == 1) | (e.bs1 == keep)]
    res = {}
    for k in PRED:
        res[k] = per_prot_auc(sub, k, min_cls=5)
    common = res["am"].index.intersection(res["esm1b"].index).intersection(res["eve"].index)
    for k in PRED:
        circ_rows.append({"benign_set": lab, "predictor": k, "n_prot": len(common),
                          "macro": res[k][common].mean()})
    g1 = boot_diff(res["am"][common], res["esm1b"][common])
    g2 = boot_diff(res["am"][common], res["eve"][common])
    log(f"   {lab:22s} n_prot={len(common):3d} | " + " | ".join(
        f"{NAMES[k]} {res[k][common].mean():.3f}" for k in PRED)
        + f" | AM-ESM1b {fmt_ci(g1[0], g1[1])} | AM-EVE {fmt_ci(g2[0], g2[1])}")
pd.DataFrame(circ_rows).to_csv(RES/"clin_bs1.csv", index=False)

# ------------------------------------------------------------------ figures
C = pd.DataFrame(ctx_rows)
ctxs = ["ordered-buried", "ordered-exposed", "IDR"]
fig, ax = plt.subplots(figsize=(7.2, 3.3), dpi=200)
for i, k in enumerate(PRED):
    sub = C[C.predictor == k].set_index("context").loc[ctxs]
    x = np.arange(3) + (i - 1) * 0.24
    ax.bar(x, sub.macro, width=0.22, color=COL[k], label=NAMES[k])
    ax.errorbar(x, sub.macro, yerr=[sub.macro - sub.lo, sub.hi - sub.macro],
                fmt="none", ecolor="#52514e", elinewidth=1)
n = C[C.predictor == "am"].set_index("context").loc[ctxs].n_prot
ax.set_xticks(range(3)); ax.set_xticklabels([f"{k}\n(n={v} proteins)" for k, v in n.items()], fontsize=8)
ax.set_ylim(0.5, 1.0); ax.set_ylabel("per-protein AUROC (macro)", fontsize=8)
ax.grid(axis="y", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper right")
fig.tight_layout(); fig.savefig(FIG/"fig_clin_auroc_context.png"); plt.close(fig)

cx = ["ordered-buried", "ordered-exposed", "IDR", "IDR+motif/PTM"]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=200, sharey=False)
for ax, met, ttl in [(axes[0], "fpr", "False-positive rate (benign called pathogenic)"),
                     (axes[1], "tpr", "Sensitivity (pathogenic called pathogenic)")]:
    for i, k in enumerate(PRED):
        sub = O[O.predictor == k].set_index("context").loc[cx]
        ax.bar(np.arange(4) + (i - 1) * 0.24, sub[met], width=0.22, color=COL[k], label=NAMES[k])
    ax.set_xticks(range(4)); ax.set_xticklabels(["buried", "exposed", "IDR", "IDR motif\n/PTM"], fontsize=8)
    ax.set_title(ttl, fontsize=8.5); ax.set_ylim(0, 1)
    ax.grid(axis="y", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
axes[0].legend(frameon=False, fontsize=7.5)
fig.tight_layout(); fig.savefig(FIG/"fig_clin_operating_points.png"); plt.close(fig)

open(RES/"clin_summary.txt", "w").write("\n".join(lines) + "\n")
c.to_parquet(OUT/"clinical_analysis.parquet")
print(f"\n[out] {RES}/clin_*.csv, figures in {FIG}")
