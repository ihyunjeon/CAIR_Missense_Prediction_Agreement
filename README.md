# CAIR Summer 2026 — variant-effect predictor disagreement

When AlphaMissense, ESM-1v and EVE disagree about a missense variant, which one
should you believe, and what does the disagreement itself tell you?

The three predictors differ along axes that are usually confounded. AlphaMissense
uses structure and is tuned against population frequency; ESM-1v is a protein
language model with neither; EVE is MSA-based with no structure and no frequency
supervision. Differencing them isolates one axis at a time. The project models
that disagreement as a function of alignment depth, structural context and
sequence-length cropping, then adjudicates it against deep mutational scanning
data, where there is a ground truth that owes nothing to ClinVar.

**Status: Phase 0 (feasibility spike) complete.** The join architecture is
validated end to end on five genes. Phases 1–3 are not yet built.

## What works today

`scripts/spike/` builds one fully joined variant table — ClinVar labels ×
AlphaMissense × ESM-1v × EVE × AlphaFold structure — for TP53, BRCA1, CFTR, PTEN
and MLH1, and asserts the joins are correct rather than assuming it.

The load-bearing check is a reference-amino-acid assertion: for every row,
`alphafold_seq[pos-1] == clinvar_ref_aa == alphamissense_ref_aa`. It catches
transcript-version drift, isoform mismatch and off-by-one indexing, all of which
otherwise fail silently and return plausible rows.

**It passes 1,168 / 1,168.**

### Spike results

1,149 complete-case variants (865 pathogenic / 303 benign, 2.85:1), five genes.
All figures below were recomputed from `data/spike/spike_final_3pred.parquet`.

| gene | n | n_neg | AlphaMissense | ESM-1v | EVE |
|---|---|---|---|---|---|
| BRCA1 | 342 | 186 | 0.960 | 0.913 | **0.964** |
| TP53 | 272 | 63 | 0.987 | 0.951 | 0.957 |
| MLH1 | 166 | 29 | 0.922 | 0.836 | 0.871 |
| PTEN | 194 | 6 | 0.997 | 0.997 | 0.992 |
| CFTR | 175 | 3 | 0.992 | 0.893 | 0.994 |

| AUROC | AlphaMissense | ESM-1v | EVE |
|---|---|---|---|
| pooled | 0.974 | 0.925 | 0.896 |
| macro, all 5 genes | 0.972 | 0.918 | 0.956 |
| macro, ≥10 negatives | **0.956** | **0.900** | **0.930** |

Read the bottom row only. PTEN's 0.997 rests on six benign variants and CFTR's
on three; a macro-average that weights those equally with BRCA1's 186 is
measuring noise. EVE's 0.896-pooled-to-0.956-macro swing is the largest of the
three and is a direct demonstration of between-gene base-rate inflation — which
is why every number in this project is reported macro-averaged.

None of these are generalisation estimates: five extensively studied genes, an
imbalanced label set, and ACMG BS1/BA1 contamination, where some benign labels
are benign *because* the variant is common — the same frequency signal
AlphaMissense was tuned on. The clean protocol is gene-held-out and
gene-ID-baselined (Phase 2); the honest adjudication is DMS (Phase 3).

There is real signal to model: rank-normalised disagreement has sd 0.225, and
22.3% of variants exceed |0.25|. Within-protein rank agreement is +0.700
(AM–EVE), +0.688 (AM–ESM) and +0.621 (ESM–EVE) — no pair is redundant.

## Layout

| path | what |
|---|---|
| `scripts/esm1v_masked_marginals.py` | ESM-1v scoring. Emits the full 20-AA log-prob vector per position, so every substitution at a residue costs one forward pass; compute LLR as `lp[mut] - lp[wt]` at analysis time. Cost scales with unique positions, not variants. |
| `scripts/spike/phase0_{join,assert,finalize,eve}.py` | The Phase 0 pipeline, in order. |
| `hoffman2/` | SGE submit script and cluster environment setup for the A100 run. |
| `data/README.md` | **The lab notebook.** Acquisition commands, row counts, licences, and every trap found so far. Read this before touching the data. |
| `docs/` | Progress report and design review. |

## Data

Not in this repository. ~3.6 GB of source atlases, all regenerable from the
commands in `data/README.md`. Only that notebook and the ~1.6 MB of derived
Phase 0 parquet are tracked, so the table above can be re-checked without
re-downloading anything.

`data/eve/eve_bulk.zip` is a symlink to an external drive and dangles when it is
unplugged; nothing in the normal analysis path needs it.

**AlphaMissense is CC BY-NC-SA 4.0 — non-commercial.** Record that in any
publication. The trained weights were never released, so the precomputed atlas
is the only way to obtain AM scores.

## Environment

Use the `vep` conda env, not `ml_dev` — `ml_dev` has a broken torch install and
an `esm` package that shadows `fair-esm`.

```bash
conda create -y -n vep python=3.11
~/miniconda3/envs/vep/bin/pip install torch transformers duckdb pyarrow \
    pandas biotite biopython scikit-learn lightgbm
```

Query the large parquets with DuckDB; never load them in pandas. Pin
`OMP_NUM_THREADS=1 MKL_NUM_THREADS=1` for small-matrix work.

## Three traps that have already cost days

1. **Do not join AlphaMissense on `transcript_id`.** AM was built on an older
   Ensembl release, so MANE v1.5 versions do not match and TP53 sits on a
   different transcript entirely. An equality join returned 2.9% of rows with no
   error. `(chrom, pos, ref, alt)` is the pivot; MANE only *ranks* transcripts.
2. **EVE is keyed by UniProt entry name, not accession or gene symbol.** TP53's
   entry name is `P53_HUMAN` — the T is dropped — so a gene-symbol match loses
   the gene silently. Resolve via the UniProt REST API, never by string edits.
   The release also ships a byte-identical `BRCA1_HUMAN-checkpoint` duplicate and
   domain-only BRCA1 models that a prefix join will double-count.
3. **`SGE_TASK_ID` is the literal string `"undefined"` on non-array jobs**, so
   `${SGE_TASK_ID:-1}` does not yield `1` and `argparse type=int` rejects it.
   This only bites the single-task pilot form, so it survives testing with `-t`.

## Compute

ESM-1v is the only predictor actually executed here. Measured on 734 positions ×
5 models: 2.571 s/position on an M3 Pro (MPS), 1.590 s/position on a public
A100 — **only 1.62× faster**, against the 3–4× the original plan assumed.
Budget for 30–60k positions is **13–27 A100-hours**, above the planned 5–15 h.
The likely cause is fp32 without TF32 enabled; that is a hypothesis, not a
measurement, and is worth one benchmark before the full array.

Laptop and cluster scores cross-validate to Spearman ρ = 0.9999993, so they are
interchangeable for ranking.

Phase 3 needs no GPU at all: ProteinGym v1.3 ships precomputed ESM-1v zero-shot
scores, and AlphaMissense comes from the local atlas. The DMS arm — the
minimum publishable result — is pure data engineering. Run it first.

## Next

- Phase 1: widen the gene set. Five genes cannot separate crop effects from gene
  identity (the cropped set is exactly BRCA1 + CFTR), and BRCA1 at mean pLDDT
  41.6 dominates any structural stratification.
- Phase 2: gene-held-out, macro-averaged, gene-ID-baselined evaluation.
- Phase 3: DMS adjudication via ProteinGym.
- Still to acquire: ProteinGym substitutions, PrimateAI-3D (licence registration
  — start early), ESM-IF1 scores, gnomAD constraint metrics.
