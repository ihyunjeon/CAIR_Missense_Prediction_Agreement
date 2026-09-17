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

**Status (2026-09-17): Phases 0–3 done, with a Phase 4 and disordered-region
analysis at 2,344 proteins.** Phase 3, the minimum publishable result, is run
on 94 human ProteinGym DMS assays (321,676 variants, 80 proteins). The clinical
arm is widened from 400 genes to 2,344 proteins (58,253 variants) using
ProteinGym's precomputed scores, so no GPU was needed. See
[Phase 3 results](#phase-3--dms-adjudication-and-disordered-regions) and
`data/README.md`. Still not done: fresh ESM-1v on the 400-gene ClinVar cohort
(ESM-1b stands in, see below), PrimateAI-3D, gnomAD constraint, and mmseqs
paralog clustering of the CV splits.

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

## Trial Phase 1 + 2 (400 genes)

A rehearsal of the real Phases 1 and 2, start to finish on the laptop in under
four minutes (`scripts/trial/`, outputs in `data/trial/`). 18,700 variants over
398 genes, 103 with >=10 of both classes. ESM-1v is excluded — this cohort's
15,367 unique positions cost **6.8 A100-hours**, which is the Hoffman2 array's
job, and Phase 2 is a structure-only baseline that does not need it.

Structure-only features, `HistGradientBoostingClassifier`, macro per-gene AUROC,
bootstrap over genes:

| model | macro AUROC | 95% CI | pooled |
|---|---|---|---|
| Structure-only, **gene-held-out** | **0.862** | [0.838, 0.884] | 0.839 |
| Structure-only, random split | 0.887 | [0.868, 0.904] | 0.940 |
| Gene-ID-only (circular) | 0.500 | — | **0.805** |
| AlphaMissense *(reference)* | 0.972 | [0.962, 0.981] | 0.966 |
| EVE *(reference)* | 0.951 | [0.938, 0.962] | 0.934 |

**The gene-ID null is the point.** Knowing only the gene name gives pooled AUROC
0.805; the structure model gets 0.839 pooled, just +0.034 over it. The same
model is 0.862 macro against a macro null of 0.500. Pooled and macro tell
opposite stories because pooled credits the model for between-gene base rates it
never learned. Report macro.

Holding genes out costs 0.025 macro AUROC (0.887 -> 0.862) — the memorisation
the design review warned about, now measured rather than asserted.

This was first run at 150 genes, where the structure model beat the gene-ID null
by only +0.001 pooled. **That did not replicate**: on 54 evaluable genes it was
largely noise, and the gap is 34x larger at 400. `data/README.md` records both
runs and the correction — effect sizes on a few dozen genes are not stable,
which is the confidence-interval argument for going to full scale.

Two new silent AlphaFold traps and a hard/soft split of the assertion gate came
out of the scale-up; all three are written up in `data/README.md`.

## Phase 3 — DMS adjudication and disordered regions

`scripts/phase3/`, results in `results/phase3/`, figures in `docs/figures/`.
All statistics are computed within a protein (or assay) and averaged over
proteins, with 95% CIs from a bootstrap over proteins.

**1. On experimental ground truth the three predictors are tied.** Mean
per-protein Spearman with DMS: AlphaMissense 0.453, ESM-1v 0.439, EVE 0.443;
AM − ESM-1v = +0.015 [−0.008, +0.039]. On ClinVar-derived labels AM leads
clearly (macro AUROC 0.961 vs ESM-1b 0.928 vs EVE 0.937). That a clinical
lead does not reappear against experiment is the project's central
observation. One caveat is still open: the two arms cover different protein
sets. Restricting the benign class to non-BS1 variants barely moves the lead,
so allele-frequency benign labels do not explain it.

**2. Disagreement is about the substitution, not the site.** Substitution
chemistry explains 0.15–0.16 of the gene-held-out variance in the signed
AM-minus-language-model rank difference. Structure explains ≤0.017; disorder,
motifs, protein length and MSA depth explain essentially none. When they
disagree, DMS sides with AM on introduced prolines (0.593), lost Cys/Met
(0.581) and introduced rare residues W/C/M (0.543). These are the cases where a
masked language model's amino-acid frequency prior misleads it.

**3. Who wins depends on what the assay measures.** AM beats ESM-1v in buried
cores (0.556 [0.528, 0.584]) and on activity assays (0.541), but loses to EVE
on thermodynamic-stability assays (0.440 [0.410, 0.469]).

**4. IDRs: the prior is right, the exceptions are missed.** After fixing the
collagen trap (below), 6.2% of IDR variants are pathogenic, against 42–78% in
ordered regions (within-protein OR 0.058). DMS confirms the premise with no
ClinVar circularity: 26 points fewer deleterious variants. All predictors have
learned this. Adding disorder features to AM gains +0.003 AUROC, and nothing
within IDRs. What they miss is the exceptions. At published cutoffs all three
catch only ~40% of pathogenic IDR variants. Short linear motifs (ELM) carry 11×
the pathogenic odds of other IDR residues. At IDR motif and PTM sites,
conservation-based EVE catches 65% against AM's 53%, at a higher
false-positive rate.

**The collagen trap.** AlphaFold models one collagen chain, so the Gly-X-Y
triple helix looks disordered to every structure-based disorder proxy. Left
in, it made "IDR" variants 17.5% pathogenic and IDR glycine substitutions
62.5% pathogenic, which inverted the biology. It is now its own context class.

## Layout

| path | what |
|---|---|
| `scripts/phase3/p3_0{1..6}_*.py` | Phase 3/4: UniProt mapping, ProteinGym joins, residue features (structure, disorder, motifs), DMS adjudication, clinical + IDR analysis, disagreement decomposition. |
| `results/phase3/` | The small result tables and text summaries behind every number above. |
| `scripts/esm1v_masked_marginals.py` | ESM-1v scoring. Emits the full 20-AA log-prob vector per position, so every substitution at a residue costs one forward pass; compute LLR as `lp[mut] - lp[wt]` at analysis time. Cost scales with unique positions, not variants. |
| `scripts/spike/phase0_{join,assert,finalize,eve}.py` | The Phase 0 pipeline, in order. |
| `scripts/trial/trial_0{1,2,3}_*.py` | The 400-gene trial: gene selection + structures, joined table, Phase 2 baseline. |
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

## Reproducing the trial

`data/` is untracked, so a fresh clone has the code and the derived results but
none of the source atlases. Four inputs are needed, all acquired with the
commands in `data/README.md`:

    data/clinvar/variant_summary.txt.gz       422 MB
    data/mane/mane_summary.txt.gz             1 MB
    data/alphamissense/am_hg38.parquet        2.8 GB  (converted from the TSV)
    data/eve/eve_variants.parquet             424 MB

AlphaFold structures are fetched by the pipeline itself (~285 MB for 400 genes)
and land in `data/alphafold/`. Then:

    export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    python scripts/trial/trial_01_genes.py 400   # select genes, fetch structures
    python scripts/trial/trial_02_build.py       # join, assert, structure, EVE
    python scripts/trial/trial_03_baseline.py    # Phase 2

About four minutes end to end on an M3 Pro. Step 1 takes a gene count as its
only argument.

Phase 3 additionally needs the ProteinGym v1.3 files and the ELM instance table
(commands in `data/README.md`), then:

    for s in 01_map 02_build 03_features 04_dms 05_clinical 06_disagreement; do
        python scripts/phase3/p3_$s.py; done

Roughly 30 minutes the first time (UniProt API, 2,055 AlphaFold downloads,
SASA), about 10 minutes after that.

## Licence

Code is MIT (see `LICENSE`). **This does not extend to the data.** AlphaMissense
is CC BY-NC-SA 4.0 (non-commercial) and EVE, ClinVar, MANE and AlphaFold each
carry their own terms — none of that data is redistributed here, and the
derived parquet in `data/` inherits the source terms. Record the AlphaMissense
licence in any publication.

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

## Phases

The design review (`docs/plan_review.html`) stages the work 0-5. Phase 0 is
done; 4 and 5 are explicitly upside rather than load-bearing, and there is a
defensible result at the end of Phase 3.

| # | phase | owner | deliverable |
|---|---|---|---|
| 0 | Spike | Ihyun | One joined row end to end. **Done** — assertion 1,168/1,168 |
| 1 | Data | Labels owner · Ihyun | The joined table at 300-500 genes: 4 predictors, DSSP/RSA/pLDDT/contact number/Pfam, gnomAD AF + LOEUF + Neff + crop flag, plus ProteinGym DMS as a second label set. **Done at 2,344 proteins** via ProteinGym clinical (ESM-1b for ESM-1v; no gnomAD LOEUF, no Pfam) |
| 2 | Baseline | Baseline owner | Structure-only classifier, gene-held-out CV, macro AUROC, bootstrap-over-genes CIs, against a gene-ID-only null. **Done** (400-gene trial; context-only prior at 2,344 proteins: macro 0.906) |
| 3 | Core | DMS owner | **Minimum publishable result.** On DMS, restrict to disagreement cases and ask which predictor tracks experimental fitness, stratified by structure. **Done** — 94 assays, 80 proteins |
| 4 | Model | Modelling owner | Predict the signed disagreement residual. Grouped conditional permutation, LOCO with refitting, variance decomposition as a Venn with the overlap made explicit. **Lite done** — grouped LOCO + single-group R², unique/shared split |
| 5 | Stress | Ihyun | Stratify by gnomAD-observed vs unobserved; replicate on DMS. **Partial** — BS1 benign-label split only; no gnomAD-observed split |

**Phase 3 does not depend on Phase 1** — its inputs are ProteinGym (which ships
its own ESM-1v scores) and the local AlphaMissense atlas. It needs no GPU and is
unblocked now. Run it first and let Phase 1's ClinVar arm proceed behind it.

The plan assigns no durations beyond "week 1" for Phase 0, and judges the
unstaged scope at roughly 2x what five students finish in a summer.

## Immediate next actions

- **Explain the clinical-vs-DMS gap.** AM leads on ClinVar-derived labels but is
  tied on DMS, and the BS1 split does not account for it. Candidates to test:
  label curation that used earlier predictors, DMS assay-phenotype mismatch,
  and protein-set differences (80 DMS proteins vs 287 clinical ones — rerun the
  clinical comparison on the overlap).
- **IDR motifs are the open clinical problem.** Build a motif-aware test set
  (ELM, phospho-switches, degrons) and ask whether a conservation-weighted
  score recovers the ~60% of pathogenic IDR variants every predictor misses.
- Run fresh ESM-1v on the 400-gene ClinVar cohort (6.8 A100-h) to replace the
  ESM-1b stand-in; test the TF32 hypothesis first.
- Handle obligate oligomers beyond collagen (coiled coils, keratins) in the
  disorder call, and add gnomAD LOEUF + Pfam for the full Phase 1 table.
- **Start PrimateAI-3D licence registration.** The ProteinGym column is the
  original PrimateAI (macro 0.869), not -3D.
- Record which GPU models each teammate can reach. A non-member's `highp` job
  waits forever with no error.
