"""
Shared helpers for the Phase 3 analyses.

Aggregation rule used everywhere (the project's pooled-vs-macro lesson):
a statistic is computed within an assay (DMS) or gene (clinical), averaged over
assays WITHIN a protein, then averaged over proteins. Proteins are the
independent unit, so confidence intervals bootstrap over proteins.
"""
import numpy as np, pandas as pd
from pathlib import Path

OUT = Path("data/phase3")
FIG = Path("docs/figures"); FIG.mkdir(parents=True, exist_ok=True)
RES = Path("results/phase3"); RES.mkdir(parents=True, exist_ok=True)   # small, tracked
RNG = np.random.default_rng(0)

# --- substitution chemistry -------------------------------------------------
# Grantham (1974) distance, computed from his composition / polarity / volume
# formula rather than a transcribed 20x20 table.
_GC = dict(A=0, R=.65, N=1.33, D=1.38, C=2.75, Q=.89, E=.92, G=.74, H=.58, I=0,
           L=0, K=.33, M=0, F=0, P=.39, S=1.42, T=.71, W=.13, Y=.20, V=0)
_GP = dict(A=8.1, R=10.5, N=11.6, D=13.0, C=5.5, Q=10.5, E=12.3, G=9.0, H=10.4,
           I=5.2, L=4.9, K=11.3, M=5.7, F=5.2, P=8.0, S=9.2, T=8.6, W=5.4, Y=6.2, V=5.9)
_GV = dict(A=31, R=124, N=56, D=54, C=55, Q=85, E=83, G=3, H=96, I=111, L=111,
           K=119, M=105, F=132, P=32.5, S=32, T=61, W=170, Y=136, V=84)
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)                                    # Kyte-Doolittle
CHARGE = dict(D=-1, E=-1, K=1, R=1)
AROM = set("FWY")


def grantham(a, b):
    return np.sqrt(1.833 * (_GC[a] - _GC[b])**2 + 0.1018 * (_GP[a] - _GP[b])**2
                   + 0.000399 * (_GV[a] - _GV[b])**2)


def chemistry(df, wt="wt", mt="mt"):
    w, m = df[wt], df[mt]
    df["grantham"] = [grantham(a, b) for a, b in zip(w, m)]
    df["d_hydro"] = m.map(KD) - w.map(KD)
    df["d_charge"] = m.map(CHARGE).fillna(0) - w.map(CHARGE).fillna(0)
    df["to_pro"] = m == "P"
    df["from_gly"] = w == "G"
    df["cys_change"] = (w == "C") | (m == "C")
    df["arom_loss"] = w.isin(AROM) & ~m.isin(AROM)
    return df


CHEM = ["grantham", "d_hydro", "d_charge", "to_pro", "from_gly", "cys_change", "arom_loss"]


def context(df):
    """Structural context. IDR is the primary disorder call; collagen triple
    helix is its own class (AlphaFold models it as disordered, it is not)."""
    idr = df.idr.astype("boolean")
    c = np.where(idr.isna(), None,
                 np.where(df.triple_helix.fillna(False).astype(bool), "collagen-helix",
                 np.where(idr.fillna(False), "IDR",
                          np.where(df.rsa < 0.25, "ordered-buried", "ordered-exposed"))))
    df["context"] = c
    return df


def macro(per_unit, unit_to_protein):
    """per_unit: Series indexed by assay/gene -> value. Mean within protein first."""
    s = per_unit.dropna()
    return s.groupby(s.index.map(unit_to_protein)).mean()


def boot(per_protein, n=2000, stat=np.mean):
    v = np.asarray(per_protein.dropna())
    if len(v) < 3: return (np.nan, np.nan)
    d = [stat(RNG.choice(v, len(v), replace=True)) for _ in range(n)]
    return tuple(np.percentile(d, [2.5, 97.5]))


def boot_diff(a, b, n=2000):
    """Paired bootstrap over proteins of mean(a - b)."""
    j = a.dropna().index.intersection(b.dropna().index)
    v = (a[j] - b[j]).values
    if len(v) < 3: return np.nan, (np.nan, np.nan), np.nan
    d = np.array([RNG.choice(v, len(v), replace=True).mean() for _ in range(n)])
    p = 2 * min((d <= 0).mean(), (d >= 0).mean())
    return v.mean(), tuple(np.percentile(d, [2.5, 97.5])), max(p, 1 / n)


def fmt_ci(m, ci, nd=3):
    return f"{m:.{nd}f} [{ci[0]:.{nd}f}, {ci[1]:.{nd}f}]"
