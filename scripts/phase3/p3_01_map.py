"""
Phase 3, step 1: resolve every ProteinGym protein to a UniProt canonical
accession, in the frame the AlphaMissense atlas and AlphaFold share.

Two inputs, two keys, one trap each:

  * DMS assays are keyed by UniProt ENTRY NAME (P53_HUMAN) -- the same trap as
    EVE (README trap 2). Resolved through the UniProt REST API, never by string
    edits. The assayed `target_seq` is often a domain or construct, so its
    numbering can be offset from the canonical sequence.
  * The clinical set is keyed by RefSeq protein (NP_000007.1). A RefSeq protein
    can differ from the UniProt canonical by isoform or by single residues, and
    AlphaMissense carries canonical accessions only (no '-2' isoforms).

Both are resolved by one rule: find the integer offset that maps the ProteinGym
sequence onto the canonical, and accept the frame if >= 85% of the positions
ProteinGym mutates carry the same residue at the shifted canonical position.

85%, not 99%, and the difference is deliberate. Fourteen human assays land
between 0.89 and 0.99: the offset is right, but the assayed construct carries
engineered background substitutions (most are Tsuboyama 2023 designed domains,
plus the RASK and SYUA constructs). The frame is correct and only a handful of
residues differ, so rejecting the whole assay would discard good data. This
step only chooses the frame; per-variant residue identity is a separate hard
gate in p3_02, which drops exactly the individual rows whose residue disagrees.
FRAME_MIN is therefore loose on purpose and the hard gate is not.

UniProt features are pulled in the same call: PTM sites, motifs, active and
binding sites, and UniProt's 'Disordered' regions (MobiDB-lite) -- the second
disorder definition used in the IDR analysis.
"""
import json, time, re, urllib.request, urllib.parse
from pathlib import Path
from collections import Counter
import pandas as pd

PG = Path("data/proteingym"); OUT = Path("data/phase3"); OUT.mkdir(exist_ok=True)
FRAME_MIN = 0.85
FIELDS = ("accession,id,reviewed,organism_id,length,sequence,ft_mod_res,ft_carbohyd,"
          "ft_lipid,ft_crosslnk,ft_act_site,ft_binding,ft_site,ft_motif,ft_region,"
          "ft_zn_fing,ft_dna_bind,ft_transmem,ft_signal")


def get(url, data=None, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data),
                                        timeout=300) as r:
                return r.read()
        except Exception:
            if i == tries - 1: raise
            time.sleep(4 * (i + 1))


def idmap(ids, frm):
    """UniProt ID-mapping job -> list of (from_id, entry_json)."""
    body = urllib.parse.urlencode({"from": frm, "to": "UniProtKB",
                                   "ids": ",".join(ids)}).encode()
    job = json.loads(get("https://rest.uniprot.org/idmapping/run", body))["jobId"]
    for _ in range(200):
        st = json.loads(get(f"https://rest.uniprot.org/idmapping/status/{job}"))
        if "jobStatus" not in st or st["jobStatus"] == "FINISHED": break
        if st["jobStatus"] in ("ERROR", "FAILURE"): raise RuntimeError(st)
        time.sleep(3)
    url = (f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job}"
           f"?format=json&fields={FIELDS}")
    return [(r["from"], r["to"]) for r in json.loads(get(url))["results"]]


def best_offset(pg_seq, canon, positions):
    """Offset o with canon[p+o-1] == pg_seq[p-1] over the mutated positions."""
    if pg_seq in canon:
        return canon.index(pg_seq), 1.0
    k, idx, votes = 8, {}, Counter()
    for i in range(len(canon) - k + 1):
        idx.setdefault(canon[i:i+k], []).append(i)
    for j in range(0, len(pg_seq) - k + 1, 4):
        for i in idx.get(pg_seq[j:j+k], []):
            votes[i - j] += 1
    votes[0] += 0                                   # always try the null offset
    best = (None, 0.0)
    for o, _ in votes.most_common(6):
        ok = [1 <= p+o <= len(canon) and canon[p+o-1] == pg_seq[p-1] for p in positions]
        f = sum(ok) / max(len(ok), 1)
        if f > best[1]: best = (o, f)
    return best


def mut_positions(path):
    m = pd.read_csv(path, usecols=["mutant"]).mutant.astype(str)
    m = m[~m.str.contains(":")]                     # singles only
    return sorted({int(re.match(r"[A-Z](\d+)[A-Z]", s).group(1)) for s in m})


def features(entry):
    out = []
    for f in entry.get("features", []):
        loc = f.get("location", {})
        s, e = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
        if s is None or e is None: continue
        out.append({"type": f["type"], "start": s, "end": e,
                    "desc": f.get("description", "")})
    return out


def reviewed(e):
    return e.get("entryType", "").startswith("UniProtKB reviewed")


rows, feats = [], {}

# --------------------------------------------------------------- DMS (human)
meta = pd.read_csv(PG/"DMS_substitutions.csv")
meta = meta[meta.taxon == "Human"]
names = sorted(meta.UniProt_ID.unique())
print(f"[1] DMS: {len(meta)} human assays, {len(names)} entry names -> UniProt ...")
emap = {frm: e for frm, e in idmap(names, "UniProtKB_AC-ID") if reviewed(e)}
print(f"    resolved {len(emap)}/{len(names)}")
for _, a in meta.iterrows():
    e = emap.get(a.UniProt_ID)
    if e is None:
        rows.append({"set": "dms", "pg_id": a.DMS_id, "status": "unresolved"}); continue
    canon = e["sequence"]["value"]
    pos = mut_positions(PG/"dms"/"DMS_ProteinGym_substitutions"/a.DMS_filename)
    o, f = best_offset(a.target_seq, canon, pos)
    feats[e["primaryAccession"]] = features(e)
    rows.append({"set": "dms", "pg_id": a.DMS_id, "pg_file": a.DMS_filename,
                 "uniprot_id": e["primaryAccession"], "entry_name": e["uniProtkbId"],
                 "canon_len": len(canon), "pg_len": len(a.target_seq),
                 "offset": o, "frame_match": f, "exact_seq": canon == a.target_seq,
                 "n_pos": len(pos),
                 "status": "ok" if f >= FRAME_MIN else "frame_mismatch",
                 "canon_seq": canon})

# --------------------------------------------------------------- clinical
cmeta = pd.read_csv(PG/"clinical_substitutions.csv")
ids = cmeta.DMS_id.tolist()
print(f"[2] clinical: {len(ids)} RefSeq proteins -> UniProt ...")
cands = {}
for i in range(0, len(ids), 500):
    for frm, e in idmap(ids[i:i+500], "RefSeq_Protein"):
        if reviewed(e) and e.get("organism", {}).get("taxonId") == 9606:
            cands.setdefault(frm, []).append(e)
    print(f"    batch {i//500 + 1}: {len(cands)} mapped so far")
for _, a in cmeta.iterrows():
    es = cands.get(a.DMS_id, [])
    if not es:
        rows.append({"set": "clinical", "pg_id": a.DMS_id, "status": "unresolved"}); continue
    pos = mut_positions(PG/"clin"/a.DMS_filename)
    best = None
    for e in es:                                    # several entries: best frame wins
        canon = e["sequence"]["value"]
        o, f = best_offset(a.target_seq, canon, pos)
        key = (canon == a.target_seq, f, -abs(o or 0))
        if best is None or key > best[0]: best = (key, e, o, f)
    (exact, _, _), e, o, f = best
    canon = e["sequence"]["value"]
    feats[e["primaryAccession"]] = features(e)
    rows.append({"set": "clinical", "pg_id": a.DMS_id, "pg_file": a.DMS_filename,
                 "uniprot_id": e["primaryAccession"], "entry_name": e["uniProtkbId"],
                 "canon_len": len(canon), "pg_len": len(a.target_seq),
                 "offset": o, "frame_match": f, "exact_seq": exact, "n_pos": len(pos),
                 "status": "ok" if f >= FRAME_MIN else "frame_mismatch",
                 "canon_seq": canon})

df = pd.DataFrame(rows)
df.to_parquet(OUT/"protein_map.parquet")
json.dump(feats, open(OUT/"uniprot_features.json", "w"))
print("\n[out] status by set:")
print(df.groupby(["set", "status"]).size().to_string())
ok = df[df.status == "ok"]
c = ok[ok.set == "clinical"]
print(f"clinical: exact canonical sequence {int(c.exact_seq.sum())}/{len(c)}; "
      f"accessions shared by >1 RefSeq protein: {int(c.uniprot_id.duplicated().sum())}")
d = ok[ok.set == "dms"]
print(f"dms: nonzero offsets {int((d.offset != 0).sum())}/{len(d)}; "
      f"distinct accessions {d.uniprot_id.nunique()}")
