"""
Trial Phase 2: the structure-only classifier, honestly benchmarked.

The design review is blunt about what makes this comparison dishonest, and all
four guards are implemented here:

  1. GENE-HELD-OUT splits. A structure model trained on ClinVar can reach a
     respectable AUROC by memorising "this gene's known hotspot domain".
     Random splits let it do that; GroupKFold on gene does not.
  2. MACRO per-gene AUROC, not pooled. Pooled AUROC is dominated by between-gene
     base-rate variation. Only genes with >=10 negatives are averaged -- the
     spike's PTEN 0.997 rested on 6 benign variants and was pure noise.
  3. A GENE-ID-ONLY NULL. "Does structure predict pathogenicity" is the wrong
     question (weakly, yes). The question is whether it adds anything over the
     gene-level prior.
  4. BOOTSTRAP OVER GENES for CIs, because genes, not variants, are the
     independent unit.

AlphaMissense and EVE appear as reference predictors only. They are NOT
comparable to the structure model on equal terms -- they never saw ClinVar,
and the structure model is trained on it. That asymmetry is the point of
guard 3, not something the numbers alone resolve.
"""
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import HistGradientBoostingClassifier
from pathlib import Path

rng = np.random.default_rng(0)
OUT = Path("data/trial")
MIN_NEG = 10

df = pd.read_parquet(OUT/"trial_table.parquet")
STRUCT = ["plddt", "rsa", "sasa", "wcn", "rel_pos", "prot_len"]
df = df.dropna(subset=STRUCT + ["label_path"]).reset_index(drop=True)
y = df.label_path.astype(int).values
groups = df.gene.values
X = df[STRUCT].astype(float).values
print(f"[in] {len(df):,} variants | {df.gene.nunique()} genes | "
      f"{y.sum():,} P/LP vs {(1-y).sum():,} B/LB")

# genes with enough negatives for a per-gene AUROC to mean anything
per_gene = df.groupby("gene").label_path.agg(n="size", n_pos="sum")
per_gene["n_neg"] = per_gene.n - per_gene.n_pos
EVAL = set(per_gene.index[(per_gene.n_neg >= MIN_NEG) & (per_gene.n_pos >= MIN_NEG)])
print(f"     evaluable genes (>={MIN_NEG} of BOTH classes): {len(EVAL)}")


def macro_auroc(frame, score_col, genes=EVAL):
    out = {}
    for g, sub in frame[frame.gene.isin(genes)].groupby("gene"):
        yy = sub.label_path.astype(int)
        if yy.nunique() < 2: continue
        s = sub[score_col]
        m = s.notna()
        if m.sum() < 10 or yy[m].nunique() < 2: continue
        out[g] = roc_auc_score(yy[m], s[m])
    return pd.Series(out)


def boot_ci(per_gene_auc, n=2000):
    """Bootstrap over GENES -- genes are the independent unit, not variants."""
    if len(per_gene_auc) < 2: return (np.nan, np.nan)
    v = per_gene_auc.values
    draws = [np.mean(rng.choice(v, len(v), replace=True)) for _ in range(n)]
    return tuple(np.percentile(draws, [2.5, 97.5]))


def fit_cv(X, y, groups, splitter, name):
    oof = np.full(len(y), np.nan)
    for tr, te in splitter:
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06,
                                           max_leaf_nodes=31, random_state=0)
        m.fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    return oof


print("\n[1] structure-only model, GENE-HELD-OUT (the honest protocol) ...")
gkf = list(GroupKFold(n_splits=5).split(X, y, groups))
df["struct_heldout"] = fit_cv(X, y, groups, gkf, "gene-held-out")

print("[2] same model, RANDOM splits (the circular protocol, for contrast) ...")
skf = list(StratifiedKFold(5, shuffle=True, random_state=0).split(X, y))
df["struct_random"] = fit_cv(X, y, groups, skf, "random")

print("[3] gene-ID-only null: the training-set prevalence of that gene ...")
# Under random splits the gene is seen in training, so this is available and is
# exactly the shortcut guard 1 exists to block. Under gene-held-out splits the
# held-out gene is unseen, so the best available guess is the global prior --
# a constant, which is AUROC 0.5 by construction. That contrast IS the result.
gid_rand = np.full(len(y), np.nan)
for tr, te in skf:
    prev = pd.Series(y[tr]).groupby(df.gene.values[tr]).mean()
    gid_rand[te] = pd.Series(df.gene.values[te]).map(prev).fillna(y[tr].mean()).values
df["gene_id_random"] = gid_rand
df["gene_id_heldout"] = y.mean()          # constant: unseen gene, no information
# Full-data gene prevalence: the maximally circular version, and the one whose
# POOLED AUROC is the number that matters -- it is what a model gets for free
# from the gene label alone.
df["gene_id_full"] = df.gene.map(df.groupby("gene").label_path.mean())

# NOTE on the macro AUROC of the gene-ID baselines: it is 0.5 BY CONSTRUCTION,
# because the score is constant within a gene and a constant cannot rank. The
# cross-validated variant (gene_id_random) scores slightly BELOW 0.5 for a
# mechanical reason, not a biological one: its score varies only across CV
# folds, and a gene's train-prevalence excludes the test fold, so a fold holding
# more pathogenic variants leaves a lower train prevalence -- a systematic
# anti-correlation. Verified: 5 distinct scores per gene, one per fold, with
# mean label decreasing monotonically in the score. Report it as 0.5.

print("[4] reference predictors (never trained on ClinVar) ...")
df["am_score"]  = df.am_pathogenicity
df["eve_score"] = df.EVE_scores_ASM

rows = []
for label, col in [("Structure-only, gene-held-out", "struct_heldout"),
                   ("Structure-only, RANDOM split",  "struct_random"),
                   ("Gene-ID-only (circular)",       "gene_id_full"),
                   ("Gene-ID-only, gene-held-out",   "gene_id_heldout"),
                   ("AlphaMissense (reference)",     "am_score"),
                   ("EVE (reference)",               "eve_score")]:
    pg = macro_auroc(df, col)
    lo, hi = boot_ci(pg)
    m = df[col].notna() & df.gene.isin(EVAL)
    pooled = (roc_auc_score(df.label_path[m].astype(int), df[col][m])
              if df[col][m].nunique() > 1 else 0.5)
    rows.append({"model": label, "n_genes": len(pg),
                 "macro_auroc": pg.mean() if len(pg) else np.nan,
                 "ci_lo": lo, "ci_hi": hi, "pooled_auroc": pooled})

res = pd.DataFrame(rows)
print("\n" + "="*88)
print("PHASE 2 TRIAL RESULTS   (macro = mean per-gene AUROC over evaluable genes)")
print("="*88)
print(res.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

print("\n--- what the gene-held-out penalty costs the structure model ---")
a = res.loc[res.model == "Structure-only, RANDOM split", "macro_auroc"].iat[0]
b = res.loc[res.model == "Structure-only, gene-held-out", "macro_auroc"].iat[0]
print(f"random split {a:.3f} -> gene-held-out {b:.3f}   (delta {b-a:+.3f})")
gr = res.loc[res.model == "Gene-ID-only (circular)", "pooled_auroc"].iat[0]
sp = res.loc[res.model == "Structure-only, gene-held-out", "pooled_auroc"].iat[0]
sm = res.loc[res.model == "Structure-only, gene-held-out", "macro_auroc"].iat[0]
print(f"\ngene-ID-only reaches POOLED AUROC {gr:.3f} knowing nothing but the gene name.")
print(f"The structure model's pooled AUROC is {sp:.3f} -- only {sp-gr:+.3f} over that null.")
print(f"Its MACRO AUROC is {sm:.3f} against a macro null of 0.500 ({sm-0.5:+.3f}).")
print("The two metrics tell opposite stories, and macro is the honest one:")
print("pooled credits the model for between-gene base rates it did not learn.")
print("Gene-ID macro is 0.500 by construction -- a within-gene constant cannot rank.")

print("\n--- feature importance is NOT reported ---")
print("Permutation importance on correlated structural features measures unique")
print("non-redundant contribution within one fitted model, not attribution.")
print("pLDDT, RSA and WCN are collinear by construction. That analysis is Phase 4,")
print("and needs grouped conditional permutation + LOCO refitting to mean anything.")

res.to_csv(OUT/"phase2_results.csv", index=False)
df[["gene","uniprot_id","aa_pos","label_path","struct_heldout","struct_random",
    "gene_id_random","am_score","eve_score"]].to_parquet(OUT/"phase2_predictions.parquet")
print(f"\n[out] {OUT/'phase2_results.csv'}  |  {OUT/'phase2_predictions.parquet'}")
