"""
ESM-1v masked-marginal scoring for missense variants.

Emits the FULL 20-amino-acid log-probability vector at each requested position,
not a single variant score. One masked forward pass at position i yields every
substitution at i for free, so cost scales with unique (protein, position) pairs
rather than with variant count. Compute the LLR at analysis time as
lp[mut] - lp[wt].

Input : parquet/TSV with columns  uniprot_id, sequence, pos  (pos is 1-indexed)
Output: parquet with uniprot_id, pos, wt_aa, model_idx, lp_A..lp_Y,
        cropped, window_start

    python esm1v_masked_marginals.py --positions pos.parquet --out scores.parquet
"""
import argparse, time, zlib
import pandas as pd, torch, torch.nn.functional as F
from transformers import AutoTokenizer, EsmForMaskedLM

AAS = list("ACDEFGHIKLMNPQRSTVWY")
MAX_LEN = 1022                      # ESM-1v context minus BOS/EOS


def crop(seq: str, pos0: int):
    """Window of <=MAX_LEN centred on pos0. Returns (subseq, new_pos0, start, cropped)."""
    if len(seq) <= MAX_LEN:
        return seq, pos0, 0, False
    half = MAX_LEN // 2
    start = max(0, min(pos0 - half, len(seq) - MAX_LEN))
    return seq[start:start + MAX_LEN], pos0 - start, start, True


def score_model(mid, jobs, tokenizer, batch_size, device):
    """jobs: list of (uniprot, pos1, wt, subseq, local_pos0, start, cropped)."""
    model = EsmForMaskedLM.from_pretrained(mid).to(device).eval()
    aa_ids = torch.tensor([tokenizer.convert_tokens_to_ids(a) for a in AAS], device=device)
    rows = []
    for b0 in range(0, len(jobs), batch_size):
        batch = jobs[b0:b0 + batch_size]
        enc = tokenizer([j[3] for j in batch], return_tensors="pt",
                        padding=True).to(device)
        tok_idx = torch.tensor([j[4] + 1 for j in batch], device=device)  # +1 for BOS
        enc["input_ids"][torch.arange(len(batch), device=device), tok_idx] = tokenizer.mask_token_id
        with torch.no_grad():
            logits = model(**enc).logits
        sel = logits[torch.arange(len(batch), device=device), tok_idx]     # (B, vocab)
        lp = F.log_softmax(sel.float(), dim=-1)[:, aa_ids].cpu().numpy()   # (B, 20)
        for (uni, pos1, wt, _, _, start, cropped), v in zip(batch, lp):
            rows.append({"uniprot_id": uni, "pos": pos1, "wt_aa": wt,
                         "cropped": cropped, "window_start": start,
                         **{f"lp_{a}": float(x) for a, x in zip(AAS, v)}})
    del model
    if device == "mps":
        torch.mps.empty_cache()
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--positions", required=True,
                   help="parquet/tsv with uniprot_id, sequence, pos (1-indexed)")
    p.add_argument("--out", required=True)
    p.add_argument("--models", nargs="+",
                   default=[f"facebook/esm1v_t33_650M_UR90S_{i}" for i in range(1, 6)])
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--device", default="mps" if torch.backends.mps.is_available() else "cpu")
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--n-shards", type=int, default=1)
    a = p.parse_args()

    df = (pd.read_parquet(a.positions) if a.positions.endswith(".parquet")
          else pd.read_csv(a.positions, sep="\t"))
    df = df.drop_duplicates(["uniprot_id", "pos"])
    if a.n_shards > 1:                       # shard by protein, never mid-protein
        # zlib.crc32 is stable across processes; builtin hash() is randomised
        # per-process and would partition shards inconsistently between jobs.
        keep = df["uniprot_id"].map(
            lambda u: zlib.crc32(u.encode()) % a.n_shards == a.shard)
        df = df[keep]

    jobs = []
    for r in df.itertuples(index=False):
        pos0 = r.pos - 1
        if not (0 <= pos0 < len(r.sequence)):
            raise ValueError(f"{r.uniprot_id} pos {r.pos} outside sequence (len {len(r.sequence)})")
        sub, lpos, start, cropped = crop(r.sequence, pos0)
        jobs.append((r.uniprot_id, r.pos, r.sequence[pos0], sub, lpos, start, cropped))

    print(f"{len(jobs):,} (protein, position) pairs x {len(a.models)} models "
          f"on {a.device}", flush=True)

    tokenizer = AutoTokenizer.from_pretrained(a.models[0])
    out = []
    for i, mid in enumerate(a.models, start=1):
        t = time.time()
        rows = score_model(mid, jobs, tokenizer, a.batch_size, a.device)
        for r in rows:
            r["model_idx"] = i
        out.extend(rows)
        print(f"  model {i}/{len(a.models)}  {time.time()-t:.0f}s "
              f"({(time.time()-t)/max(len(jobs),1)*1000:.0f} ms/position)", flush=True)

    pd.DataFrame(out).to_parquet(a.out, index=False)
    print(f"wrote {a.out}  ({len(out):,} rows)")


if __name__ == "__main__":
    main()
