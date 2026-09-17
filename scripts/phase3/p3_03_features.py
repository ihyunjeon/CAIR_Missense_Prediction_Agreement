"""
Phase 3, step 3: one row per (protein, residue) -- structure, disorder and the
functional annotations that make a disordered residue matter.

WHY THREE DISORDER DEFINITIONS. "IDR missense is almost always benign" is true
on average and false in exactly the places a clinician cares about, so the
definition has to be pinned down and its sensitivity shown:

  dis_plddt   AlphaFold pLDDT < 50. The common proxy. Misses conditionally
              folded IDRs, which AF2 often predicts confidently in their bound
              conformation.
  dis_rsa     25-residue rolling mean of AlphaFold RSA >= 0.581 (the 'AlphaFold-
              RSA' disorder predictor; window and cutoff are from the CAID
              literature -- background-knowledge, not re-derived here).
  dis_uniprot inside a UniProt 'Disordered' region (MobiDB-lite consensus;
              sequence-based, independent of AlphaFold).
  idr         the primary call: dis_plddt AND dis_rsa, i.e. both structure-
              derived proxies agree. dis_uniprot is kept as an orthogonal check.

WHY THE ANNOTATIONS. Disordered regions are not functionally inert. They carry
short linear motifs (SLiMs: degrons, docking, localisation signals), PTM sites
(phosphorylation switches), and low-complexity 'sticker' residues for phase
separation. A missense variant there can be fully pathogenic with no structural
consequence -- which is precisely the case a structure-tuned predictor may miss
and a conservation-based one may catch. So each residue is flagged:

  in_elm      inside an ELM true-positive motif instance (human)
  in_motif    UniProt 'Motif'
  ptm         UniProt Modified residue / Glycosylation / Lipidation / Cross-link
  func_site   UniProt Active site / Binding site / Site
  slim_or_ptm in_elm | in_motif | ptm   (the IDR 'exception' class)

THE COLLAGEN TRAP. A collagen alpha chain is a rigid triple helix, but only as a
trimer; AlphaFold models one chain, so the Gly-X-Y domain comes out as a long
low-pLDDT, fully exposed strand and both structure proxies call it disordered.
Glycine is obligatory at every third position (it is the only residue small
enough to sit at the helix axis), so Gly substitutions there are the classic
dominant-negative cause of osteogenesis imperfecta, Ehlers-Danlos and Alport
syndrome. Left in, collagen made the "IDR" Gly substitutions 96.5% pathogenic
(1,605 variants) and inverted the IDR biology. `triple_helix` is taken from
UniProt 'Triple-helical region' annotations OR >= 10 consecutive Gly-X-Y
repeats in sequence, and is removed from the IDR call.

AlphaFold DB splits proteins > 2700 aa into fragments (trial trap); those get no
structure features and are flagged, never silently truncated.
"""
import json, re, urllib.request, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import numpy as np, pandas as pd

OUT = Path("data/phase3"); AF = Path("data/alphafold")
MAX_LEN = 2700
MAXASA = {"A":129,"R":274,"N":195,"D":193,"C":167,"Q":225,"E":223,"G":104,
          "H":224,"I":197,"L":201,"K":236,"M":224,"F":240,"P":159,"S":155,
          "T":172,"W":285,"Y":263,"V":174}   # Tien et al. 2013
PTM_TYPES = {"Modified residue", "Glycosylation", "Lipidation", "Cross-link"}
SITE_TYPES = {"Active site", "Binding site", "Site"}

pm = pd.read_parquet(OUT/"protein_map.parquet")
pm = pm[pm.status == "ok"].drop_duplicates("uniprot_id")
seqs = dict(zip(pm.uniprot_id, pm.canon_seq))
print(f"[0] {len(seqs)} distinct accessions")


def fetch(acc):
    """Resolve the AF model from the API; never hardcode the version."""
    if len(seqs[acc]) > MAX_LEN:
        return acc, None, "too_long"
    have = sorted(AF.glob(f"AF-{acc}-F1-model_v*.pdb"))
    if have:
        return acc, have[-1], "cached"
    for i in range(4):
        try:
            with urllib.request.urlopen(
                    f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=60) as r:
                recs = json.load(r)
            rec = next((x for x in recs if x.get("uniprotAccession") == acc), None)
            if rec is None:
                return acc, None, "no_model"
            dest = AF/Path(rec["pdbUrl"]).name
            urllib.request.urlretrieve(rec["pdbUrl"], dest)
            return acc, dest, "fetched"
        except urllib.error.HTTPError as e:
            if e.code == 404: return acc, None, "no_model"
            time.sleep(2 * (i + 1))
        except Exception:
            time.sleep(2 * (i + 1))
    return acc, None, "error"


def roll(x, w):
    return pd.Series(x).rolling(w, center=True, min_periods=1).mean().values


def residue_feats(args):
    acc, path = args
    import biotite.structure.io.pdb as pdb
    import biotite.structure as struc
    from biotite.sequence import ProteinSequence
    canon = seqs[acc]
    arr = pdb.PDBFile.read(str(path)).get_structure(model=1, extra_fields=["b_factor"])
    arr = arr[struc.filter_amino_acids(arr)]
    ca = arr[arr.atom_name == "CA"]
    af_seq = "".join(ProteinSequence.convert_letter_3to1(r) for r in ca.res_name)
    if af_seq != canon:
        return acc, None                       # AF sequence != UniProt canonical: no features
    sasa = struc.sasa(arr, vdw_radii="Single", point_number=100)
    res_ids, inv = np.unique(arr.res_id, return_inverse=True)
    per_res = np.bincount(inv, weights=np.nan_to_num(sasa))
    d = np.linalg.norm(ca.coord[:, None, :] - ca.coord[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    wcn = (1.0 / d**2).sum(1)
    try:
        sse = struc.annotate_sse(arr)
    except Exception:
        sse = np.array([""] * len(ca))
    rsa = np.array([per_res[i] / MAXASA[a] for i, a in enumerate(canon)])
    plddt = ca.b_factor.astype(float)
    return acc, pd.DataFrame({
        "uniprot_id": acc, "aa_pos": ca.res_id.astype(int), "ref_aa": list(canon),
        "plddt": plddt, "rsa": np.clip(rsa, 0, 1.5), "wcn": wcn, "sse": sse,
        "plddt_w11": roll(plddt, 11), "rsa_w25": roll(np.clip(rsa, 0, 1), 25)})


def main():
    print("[1] AlphaFold models ...")
    with ThreadPoolExecutor(16) as ex:
        got = list(ex.map(fetch, list(seqs)))
    st = pd.Series([s for _, _, s in got]).value_counts()
    print("    " + st.to_string().replace("\n", "\n    "))
    jobs = [(a, p) for a, p, s in got if p is not None]

    print(f"[2] structure features for {len(jobs)} models ...")
    t0 = time.time()
    cache = OUT/"_struct_cache.parquet"          # SASA is the slow part; cache it
    if cache.exists():
        res = pd.read_parquet(cache)
        done_ = set(res.uniprot_id.unique())
        af_mismatch = [a for a, _ in jobs if a not in done_]
    else:
        frames, af_mismatch = [], []
        with ProcessPoolExecutor(8) as ex:
            for acc, f in ex.map(residue_feats, jobs, chunksize=4):
                if f is None: af_mismatch.append(acc)
                else: frames.append(f)
        res = pd.concat(frames, ignore_index=True)
        res.to_parquet(cache)
    print(f"    {res.uniprot_id.nunique()} ok, {len(af_mismatch)} AF-sequence != canonical "
          f"({time.time()-t0:.0f}s)")

    # proteins with no usable structure still get annotation rows (sequence-only)
    have_ = set(res.uniprot_id.unique())
    nostruct = [a for a in seqs if a not in have_]
    if nostruct:
        res = pd.concat([res] + [pd.DataFrame({"uniprot_id": a,
                         "aa_pos": np.arange(1, len(seqs[a]) + 1), "ref_aa": list(seqs[a])})
                         for a in nostruct], ignore_index=True)
    res["has_structure"] = res.plddt.notna()
    res["prot_len"] = res.uniprot_id.map({a: len(s) for a, s in seqs.items()})
    res["rel_pos"] = res.aa_pos / res.prot_len
    res["cropped"] = res.prot_len > 1022

    print("[3] collagen triple helix + disorder calls ...")
    feats = json.load(open(OUT/"uniprot_features.json"))
    th = set()
    for acc, fl in feats.items():
        for f in fl:
            if f["type"] == "Region" and f["desc"].startswith("Triple-helical"):
                th.update((acc, p) for p in range(int(f["start"]), int(f["end"]) + 1))
    n_ann = len(th)
    for acc, sq in seqs.items():               # sequence backup: >=10 Gly-X-Y repeats
        for m in re.finditer(r"(?:G..){10,}", sq):
            th.update((acc, p) for p in range(m.start() + 1, m.end() + 1))
    # pandas 3 strings are Arrow-backed: iterate plain numpy arrays, never the
    # Series, or each element becomes an Arrow scalar (~100x slower)
    keys = list(zip(res.uniprot_id.to_numpy(dtype=object), res.aa_pos.to_numpy()))
    res["triple_helix"] = [x in th for x in keys]
    print(f"    triple-helix residues: {res.triple_helix.sum():,} "
          f"({n_ann:,} from UniProt, rest from Gly-X-Y scan) in "
          f"{res[res.triple_helix].uniprot_id.nunique()} proteins")
    # nullable booleans: no structure -> <NA>, never a silent False
    res["dis_plddt"] = (res.plddt < 50).astype("boolean")
    res["dis_rsa"] = (res.rsa_w25 >= 0.581).astype("boolean")
    res["idr"] = res.dis_plddt & res.dis_rsa & ~res.triple_helix
    res.loc[~res.has_structure, ["dis_plddt", "dis_rsa", "idr"]] = pd.NA

    print("[4] annotations ...")
    elm = pd.read_csv("data/annot/elm_instances.tsv", sep="\t", comment="#")
    elm = elm[elm.Primary_Acc.isin(seqs)]
    flags = {k: set() for k in ["dis_uniprot", "in_motif", "ptm", "func_site",
                                "dna_or_zf", "tm", "in_elm"]}
    ptm_kind = {}
    for acc, fl in feats.items():
        for f in fl:
            span = range(int(f["start"]), int(f["end"]) + 1)
            t = f["type"]
            key = None
            if t == "Region" and f["desc"] == "Disordered": key = "dis_uniprot"
            elif t == "Motif": key = "in_motif"
            elif t in PTM_TYPES:
                key = "ptm"
                for p in span: ptm_kind[(acc, p)] = f["desc"].split(";")[0]
            elif t in SITE_TYPES: key = "func_site"
            elif t in ("DNA binding", "Zinc finger"): key = "dna_or_zf"
            elif t in ("Transmembrane", "Intramembrane"): key = "tm"
            if key:
                flags[key].update((acc, p) for p in span)
    for _, e in elm.iterrows():
        flags["in_elm"].update((e.Primary_Acc, p) for p in range(int(e.Start), int(e.End) + 1))
    for k, s in flags.items():
        res[k] = [x in s for x in keys]
    res["ptm_kind"] = [ptm_kind.get(x) for x in keys]
    res["slim_or_ptm"] = res.in_elm | res.in_motif | res.ptm

    s = res[res.has_structure]
    print(f"    residues with structure: {len(s):,} over {s.uniprot_id.nunique()} proteins")
    for k in ["dis_plddt", "dis_rsa", "idr", "dis_uniprot"]:
        print(f"    {k:12s} {s[k].astype(float).mean()*100:5.1f}% of residues")
    a, b = s.idr.astype(bool), s.dis_uniprot
    print(f"    idr vs dis_uniprot: both {np.mean(a&b)*100:.1f}%  idr-only "
          f"{np.mean(a&~b)*100:.1f}%  uniprot-only {np.mean(~a&b)*100:.1f}%  "
          f"(Cohen kappa {((np.mean(a==b)) - (a.mean()*b.mean() + (1-a.mean())*(1-b.mean()))) / (1 - (a.mean()*b.mean() + (1-a.mean())*(1-b.mean()))):.2f})")
    for k in ["in_elm", "in_motif", "ptm", "func_site"]:
        print(f"    {k:10s} {res[k].sum():>7,} residues | "
              f"{s.loc[s.idr.astype(bool), k].mean()*100:.2f}% of IDR vs "
              f"{s.loc[~s.idr.astype(bool), k].mean()*100:.2f}% of ordered")
    res.to_parquet(OUT/"residues.parquet")
    pd.Series(af_mismatch, name="uniprot_id").to_frame().to_parquet(OUT/"af_seq_mismatch.parquet")
    print(f"[out] {OUT/'residues.parquet'}  {len(res):,} residues")


if __name__ == "__main__":
    main()
