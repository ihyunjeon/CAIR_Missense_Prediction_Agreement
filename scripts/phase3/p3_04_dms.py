"""
Phase 3 core: adjudicate predictor disagreement against deep mutational scans.

DMS is the one label set that owes nothing to ClinVar or to population
frequency, so it is where "which predictor is right when they disagree" can be
asked without circularity. All statistics are within-assay, then averaged over
assays within a protein, then over proteins (bootstrap over proteins).

A. Overall agreement with DMS (Spearman), by assay type and alignment depth.
B. By structural context (IDR / ordered-buried / ordered-exposed).
C. Disagreement adjudication. Within an assay, every predictor and the DMS are
   converted to percentile ranks. For a pair (X, Y), a variant is a disagreement
   when |r_X - r_Y| >= 0.25 (the spike's threshold). X "wins" when its rank is
   closer to the DMS rank. Win rate 0.5 = no information about who is right.
D. IDR over-calling, calibration-free: within an assay, mean predictor rank of
   IDR variants minus mean DMS rank of the same variants. Positive = the
   predictor ranks disordered-region variants as more damaging than the
   experiment does.
E. ESM-1v vs ESM-1b agreement, which licenses ESM-1b as the language-model
   stand-in in the clinical arm.

DMS caveat, stated once: an assay measures one phenotype (abundance, binding,
growth). A variant that breaks an interaction motif can be neutral in an
abundance assay. 'Stability' assays in particular cannot see function.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from common import OUT, FIG, RES, context, chemistry, macro, boot, boot_diff, fmt_ci

PRED = ["am", "esm1v", "eve"]
REF = ["esm1b", "esmif1"]
NAMES = {"am": "AlphaMissense", "esm1v": "ESM-1v", "eve": "EVE",
         "esm1b": "ESM-1b", "esmif1": "ESM-IF1 (structure only)"}
COL = {"am": "#2a78d6", "esm1v": "#eb6834", "eve": "#1baf7a"}
MIN_N = 30

d = pd.read_parquet(OUT/"dms_table.parquet")
r = pd.read_parquet(OUT/"residues.parquet")
d = d.merge(r, how="left", on=["uniprot_id", "aa_pos"])
d = context(chemistry(d))
a2p = d.drop_duplicates("assay").set_index("assay").uniprot_id
print(f"[in] {len(d):,} variants, {d.assay.nunique()} assays, {d.uniprot_id.nunique()} proteins")
print("     context:", d.context.value_counts(dropna=False).to_dict())

for c in PRED + REF + ["dms_damage"]:
    d["r_" + c] = d.groupby("assay")[c].rank(pct=True)
d["dms_del"] = 1 - d.DMS_score_bin              # 1 = experimentally deleterious

lines = []
def log(s=""):
    print(s); lines.append(s)


def per_assay_rho(frame, col, min_n=MIN_N):
    def f(g):
        m = g[col].notna()
        return spearmanr(g[col][m], g.dms_damage[m])[0] if m.sum() >= min_n else np.nan
    return frame.groupby("assay").apply(f, include_groups=False)


# ------------------------------------------------------------------ A
log("A. Spearman with DMS damage (assay -> protein -> mean; 95% CI over proteins)")
rows, rho_p = [], {}
for c in PRED + REF:
    pp = macro(per_assay_rho(d, c), a2p)
    rho_p[c] = pp
    b = boot(pp)                     # one bootstrap per statistic: CSV and log agree
    rows.append({"stratum": "all", "predictor": c, "n_prot": len(pp),
                 "mean": pp.mean(), "lo": b[0], "hi": b[1]})
    log(f"   {NAMES[c]:26s} {fmt_ci(pp.mean(), b)}  n={len(pp)}")
for x, y in [("am", "esm1v"), ("am", "eve"), ("esm1v", "eve")]:
    m, ci, p = boot_diff(rho_p[x], rho_p[y])
    log(f"   {NAMES[x]} - {NAMES[y]}: {fmt_ci(m, ci)}  p_boot={p:.3f}")

for key in ["selection", "neff_cat"]:
    log(f"\n   by {key}:")
    for lev, sub in d.groupby(key):
        s = []
        for c in PRED:
            pp = macro(per_assay_rho(sub, c), a2p)
            b = boot(pp)
            rows.append({"stratum": f"{key}={lev}", "predictor": c, "n_prot": len(pp),
                         "mean": pp.mean(), "lo": b[0], "hi": b[1]})
            s.append(f"{NAMES[c]} {pp.mean():.3f}")
        log(f"   {lev:18s} n_prot={len(pp):2d} | " + " | ".join(s))

# ------------------------------------------------------------------ B
log("\nB. Spearman with DMS within structural context (assays with >=30 in stratum)")
for ctx in ["ordered-buried", "ordered-exposed", "IDR"]:
    sub = d[d.context == ctx]
    s = []
    for c in PRED:
        pp = macro(per_assay_rho(sub, c), a2p)
        b = boot(pp)
        rows.append({"stratum": f"context={ctx}", "predictor": c, "n_prot": len(pp),
                     "mean": pp.mean(), "lo": b[0], "hi": b[1]})
        s.append(f"{NAMES[c]} {fmt_ci(pp.mean(), b)}")
    log(f"   {ctx:16s} n_prot={len(pp):2d} | " + " | ".join(s))
pd.DataFrame(rows).to_csv(RES/"dms_spearman.csv", index=False)

# ------------------------------------------------------------------ C
log("\nC. Disagreement adjudication: who is closer to the DMS when two predictors disagree?")
PAIRS = [("am", "esm1v"), ("am", "eve"), ("esm1v", "eve")]
win_rows = []


def win_rate(frame, x, y, min_n=20):
    def f(g):
        g = g[g["r_" + x].notna() & g["r_" + y].notna()]
        dis = g[(g["r_" + x] - g["r_" + y]).abs() >= 0.25]
        if len(dis) < min_n: return np.nan
        ex = (dis["r_" + x] - dis.r_dms_damage).abs()
        ey = (dis["r_" + y] - dis.r_dms_damage).abs()
        return (ex < ey).mean()
    return frame.groupby("assay").apply(f, include_groups=False)


for x, y in PAIRS:
    frac = d.groupby("assay").apply(
        lambda g: ((g["r_" + x] - g["r_" + y]).abs() >= 0.25).mean(), include_groups=False)
    for lab, sub in [("all", d)] + [(f"context={k}", d[d.context == k])
                                    for k in ["ordered-buried", "ordered-exposed", "IDR"]] \
                   + [(f"selection={k}", g) for k, g in d.groupby("selection")]:
        pp = macro(win_rate(sub, x, y), a2p)
        m, ci, p = boot_diff(pp, pd.Series(0.5, index=pp.index))
        win_rows.append({"pair": f"{x}_vs_{y}", "stratum": lab, "n_prot": len(pp),
                         "win_x": pp.mean(), "lo": ci[0] + .5 if pp.size else np.nan,
                         "hi": ci[1] + .5 if pp.size else np.nan, "p_boot": p})
    w = [z for z in win_rows if z["pair"] == f"{x}_vs_{y}"]
    log(f"   {NAMES[x]} vs {NAMES[y]}  (disagreement rate {frac.mean()*100:.1f}% of variants)")
    for z in w:
        log(f"     {z['stratum']:30s} n_prot={z['n_prot']:2d}  {NAMES[x]} wins "
            f"{z['win_x']:.3f} [{z['lo']:.3f}, {z['hi']:.3f}]  p={z['p_boot']:.3f}")
W = pd.DataFrame(win_rows); W.to_csv(RES/"dms_winrates.csv", index=False)

# ------------------------------------------------------------------ D
log("\nD. IDR over-calling vs experiment (mean predictor rank - mean DMS rank, IDR variants)")
oc_rows = []
for c in PRED + ["esm1b"]:
    def f(g):
        g = g[(g.context == "IDR") & g["r_" + c].notna()]
        return (g["r_" + c] - g.r_dms_damage).mean() if len(g) >= MIN_N else np.nan
    pp = macro(d.groupby("assay").apply(f, include_groups=False), a2p)
    b = boot(pp)
    oc_rows.append({"predictor": c, "n_prot": len(pp), "overcall": pp.mean(),
                    "lo": b[0], "hi": b[1]})
    log(f"   {NAMES[c]:14s} {fmt_ci(pp.mean(), b)}  n_prot={len(pp)}")
# the experiment itself: are IDR variants less often deleterious?
def del_gap(g):
    i, o = g[g.context == "IDR"], g[g.context.isin(["ordered-buried", "ordered-exposed"])]
    if len(i) < MIN_N or len(o) < MIN_N: return np.nan
    return i.dms_del.mean() - o.dms_del.mean()
gap = macro(d.groupby("assay").apply(del_gap, include_groups=False), a2p)
log(f"   DMS-deleterious fraction, IDR minus ordered: {fmt_ci(gap.mean(), boot(gap))}"
    f"  n_prot={len(gap)}")
# D2. The rank-difference above is exposed to regression to the middle: most
# IDR variants are near-neutral, so their DMS ranks are largely assay noise and
# sit near the centre, while a predictor ranks them confidently low. A
# noise-robust check uses the assay's own binary call instead: of the variants
# the EXPERIMENT calls deleterious, what fraction does each predictor place in
# the top half of that assay? Compare IDR with ordered on the same assays.
log("   D2. of DMS-deleterious variants, fraction a predictor ranks in the assay's top half:")
for c in PRED:
    def hit(g, ctx):
        g = g[(g.dms_del == 1) & g["r_" + c].notna()]
        g = g[g.context == "IDR"] if ctx == "IDR" else g[g.context.isin(
            ["ordered-buried", "ordered-exposed"])]
        return (g["r_" + c] > 0.5).mean() if len(g) >= 10 else np.nan
    hi_ = macro(d.groupby("assay").apply(hit, ctx="IDR", include_groups=False), a2p)
    ho_ = macro(d.groupby("assay").apply(hit, ctx="ord", include_groups=False), a2p)
    j = hi_.dropna().index.intersection(ho_.dropna().index)
    dm, dci, _ = boot_diff(hi_[j], ho_[j])
    oc_rows[[o["predictor"] for o in oc_rows].index(c)].update(
        {"sens_idr": hi_[j].mean(), "sens_ordered": ho_[j].mean(), "sens_gap": dm,
         "gap_lo": dci[0], "gap_hi": dci[1], "n_prot_sens": len(j)})
    log(f"     {NAMES[c]:14s} IDR {hi_[j].mean():.3f} vs ordered {ho_[j].mean():.3f} "
        f"(gap {fmt_ci(dm, dci)}, n_prot={len(j)})")
pd.DataFrame(oc_rows).to_csv(RES/"dms_idr_overcall.csv", index=False)
n_idr = d[d.context == "IDR"].groupby("assay").size()
log(f"   assays with >=30 IDR variants: {(n_idr >= MIN_N).sum()} / {d.assay.nunique()}"
    f"  ({d[d.context=='IDR'].shape[0]:,} IDR variants in total)")

# ------------------------------------------------------------------ E
e = d.groupby("assay").apply(lambda g: g.esm1v.corr(g.esm1b, method="spearman"),
                             include_groups=False)
pe = macro(e, a2p)
log(f"\nE. ESM-1v vs ESM-1b within-assay Spearman: {fmt_ci(pe.mean(), boot(pe))}; "
    f"their DMS Spearman differ by {fmt_ci(*boot_diff(rho_p['esm1v'], rho_p['esm1b'])[:2])}")

# ------------------------------------------------------------------ figures
sel = ["Stability", "Expression", "Binding", "Activity", "OrganismalFitness"]
S = pd.DataFrame(rows)
fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=200)
for i, c in enumerate(PRED):
    sub = S[(S.predictor == c) & S.stratum.isin([f"selection={s}" for s in sel])]
    sub = sub.set_index("stratum").loc[[f"selection={s}" for s in sel]]
    x = np.arange(len(sel)) + (i - 1) * 0.24
    ax.bar(x, sub["mean"], width=0.22, color=COL[c], label=NAMES[c])
    ax.errorbar(x, sub["mean"], yerr=[sub["mean"] - sub.lo, sub.hi - sub["mean"]],
                fmt="none", ecolor="#52514e", elinewidth=1, capsize=0)
ax.set_xticks(np.arange(len(sel)))
ax.set_xticklabels([f"{s}\n(n={int(S[(S.stratum==f'selection={s}')].n_prot.iloc[0])})"
                    for s in sel], fontsize=8)
ax.set_ylabel("Spearman with DMS (macro over proteins)", fontsize=8)
ax.grid(axis="y", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper left")
ax.set_ylim(0, max(0.75, S["hi"].max() + 0.05))
fig.tight_layout(); fig.savefig(FIG/"fig_dms_by_selection.png"); plt.close(fig)

ctxs = ["ordered-buried", "ordered-exposed", "IDR"]
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8), dpi=200, sharey=True)
for ax, (x, y) in zip(axes, PAIRS):
    sub = W[(W.pair == f"{x}_vs_{y}")].set_index("stratum")
    labs = ["all"] + [f"context={k}" for k in ctxs]
    sub = sub.reindex(labs)
    ax.axhline(0.5, color="#8a8984", lw=0.8, ls="--")
    ax.errorbar(range(4), sub.win_x, yerr=[sub.win_x - sub.lo, sub.hi - sub.win_x],
                fmt="o", color=COL[x], ms=5, elinewidth=1.2)
    ax.set_xticks(range(4)); ax.set_xticklabels(["all", "buried", "exposed", "IDR"], fontsize=8)
    ax.set_title(f"{NAMES[x]} vs {NAMES[y]}", fontsize=9)
    ax.grid(axis="y", color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
axes[0].set_ylabel(f"P(first predictor closer to DMS)\namong disagreements", fontsize=8)
fig.tight_layout(); fig.savefig(FIG/"fig_dms_winrates.png"); plt.close(fig)

open(RES/"dms_summary.txt", "w").write("\n".join(lines) + "\n")
d.to_parquet(OUT/"dms_analysis.parquet")
print(f"\n[out] {RES}/dms_*.csv, dms_summary.txt, figures in {FIG}")
