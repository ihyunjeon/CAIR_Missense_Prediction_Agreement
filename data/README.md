# Data acquisition log

Acquired 2026-08-27. Everything here is regenerable from the commands below —
do not back it up, do not commit it.

## Storage layout

The internal disk is 460 GB and chronically tight. Split as of 2026-09-01:

| tier | lives on | what goes here |
|---|---|---|
| hot | internal `data/` | anything joined or queried per-run — the parquets, AlphaMissense |
| cold | LaCie `/Volumes/LaCie/data/` | write-once archives read rarely, reached by symlink |

`data/eve/eve_bulk.zip` is a **symlink** to `/Volumes/LaCie/data/eve/eve_bulk.zip`.
Code that opens it by path keeps working, but only while the LaCie is mounted;
with the drive unplugged the symlink dangles and `zipfile` raises `FileNotFoundError`.
The derived parquets are local, so nothing in the normal analysis path needs the drive.

Do not move AlphaMissense to the LaCie. `am_hg38` is the join pivot for the whole
project and gets hit constantly; it belongs on internal disk.

## AlphaMissense — `alphamissense/`

Precomputed scores. **AlphaMissense cannot be run**: DeepMind released the model
and loss code (Apache 2.0) as an implementation reference but withheld the
trained weights. This atlas is the only way to obtain AM scores.

    curl -L -O https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz
    curl -L -O https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
    curl -L -O https://storage.googleapis.com/dm_alphamissense/AlphaMissense_gene_hg38.tsv.gz

| file | rows | distinct UniProt | keyed by |
|---|---|---|---|
| `am_aa_substitutions.parquet` | 216,175,351 | 20,516 | `uniprot_id` + `protein_variant` (all isoforms) |
| `am_hg38.parquet` | 71,697,556 | 19,117 | genomic coords **+ `uniprot_id` + `transcript_id`** |

`am_hg38` carries chrom/pos/ref/alt *and* uniprot_id *and* transcript_id in one
row. That is the pivot for the whole project: it lets ClinVar (genomic coords)
join to AlphaFold (UniProt accession) without a separate mapping step.

Row count 71,697,556 matches the ~71M figure in Cheng et al. 2023.

`am_class` over all possible substitutions: 43.4% pathogenic, 42.2% benign,
14.4% ambiguous. Note this is the space of *possible* variants, not observed
ones — it is not a prior over anything real.

**License: CC BY-NC-SA 4.0 (non-commercial).** Fine for academic work; record it
in any publication.

### Regenerating the parquet

TSVs have 3 comment lines before the header, hence `skip=3`:

    duckdb -c "COPY (SELECT * FROM read_csv('AlphaMissense_hg38.tsv.gz',
                 delim='\t', header=true, skip=3))
               TO 'am_hg38.parquet' (FORMAT parquet, COMPRESSION zstd)"

~10s per file. Query the parquet with DuckDB, never load it in pandas.

## EVE — `eve/`

`eve_bulk.zip` (8.9 GB, integrity verified) from <https://evemodel.org/download/bulk>.
**Stored on the LaCie**, symlinked in as `eve/eve_bulk.zip`; md5
`e012cde632cfcb4cda119802dd6b27d5`, re-verified on the 2026-09-01 move:

    curl -L -o eve_bulk.zip https://evemodel.org/api/proteins/bulk/download/

**Do not unzip it.** Contents are 63.4 GB uncompressed against ~30 GB free:

| group | files | uncompressed | use |
|---|---|---|---|
| `vcf_files_missense_mutations/` | 2,951 | 28.70 GB | skip — duplicates the CSVs in genomic coords |
| `MSAs/*.a2m` | 3,220 | 24.15 GB | **the alignments EVE was trained on** |
| `variant_files/*.csv` | 3,212 | 9.93 GB | the scores |
| ROC/PRC `*.png` | 5,550 | 0.58 GB | skip — plots |

Everything is streamed out of the zip directly into parquet; nothing is
extracted. See `scripts/`-adjacent conversion code in the job log, or rerun the
snippet in this repo's history.

### `eve_variants.parquet` — 424 MB, 44,514,580 rows, 3,212 genes

Curated subset of the 43 source columns (all 43 remain in the zip). Coverage:

| column | non-null | |
|---|---|---|
| `EVE_scores_ASM` | 36,013,930 | **80.9%** — 8.5M positions have no EVE score |
| `uncertainty_ASM` | 35,997,704 | 80.9% — EVE's own uncertainty |
| `frequency_gv2` / `gv3` | 1,024,834 | 2.3% — gnomAD AF; only observed variants |
| `ClinVar_ClinicalSignificance` | 283,048 | 0.6% |
| `BS1`, `PM2` | 44,514,580 | 100% — precomputed ACMG criteria |

`EVE_scores_ASM`: min 0.006, median 0.521, mean 0.502, max 1.000.

**Three things to know before using this.**

1. **19% of positions have no EVE score.** Any join must treat missing EVE as
   missing-not-at-random — EVE fails where its MSA is thin, which correlates
   with exactly the low-Neff proteins where ESM-1v also degrades. Do not drop
   those rows silently.

2. **The bundled ClinVar is a 2021 snapshot and is mostly VUS.** Of 283,048
   labelled rows, 203,021 (72%) are "Uncertain significance"; only ~48k carry a
   confident P/LP/B/LB call, and those run ~1.7:1 pathogenic-to-benign. Pull
   current ClinVar separately; use this only as a cross-check.

3. **`BS1` is precomputed for every row.** BS1 is the ACMG *allele-frequency*
   benign criterion — the mechanism by which ClinVar benign labels are
   frequency-derived. Having it per-variant makes the label-contamination
   analysis direct rather than inferred.

**Join key:** files are named by UniProt *entry name* (`RPC22_HUMAN`), not
accession (`P0CG48`). A UniProt idmapping step is required to join to
AlphaMissense or AlphaFold.

### `eve_msa_depth.parquet` — per-gene alignment depth (3,220 genes, 81 KB)

Sequence count and query length per `.a2m`, computed streaming from the zip.
MSA depth is plausibly the largest single driver of AM/ESM-1v disagreement and
belongs to neither the structural nor the frequency block. EVE ships its own
alignments, so we get this without running jackhmmer.

**Depth spans 5.9 orders of magnitude** — 1 to 793,992 sequences, median 3,499,
p10/p90 = 749 / 39,685 (a 53x interdecile ratio). 16.3% of genes have fewer than
1,000 homologs. This is a large enough dynamic range that leaving it out of the
disagreement model would let its variance leak into the structural block.

**A confound chain worth modelling explicitly.** Longer proteins have shallower
alignments *and* get cropped by ESM-1v's 1022-token context — two degradations
that co-occur:

| query length | genes | median MSA depth |
|---|---|---|
| <250 | 699 | 4,895 |
| 250-500 | 1,250 | 5,632 |
| 500-1000 | 863 | 2,528 |
| 1000-1022 | 15 | 2,712 |
| **>1022 (cropped)** | **393 (12.2%)** | **2,032** |

The depth/length correlation is statistically overwhelming but **modest in
magnitude**: Spearman rho = -0.219 (p = 2.2e-36), Pearson on log10-log10 =
-0.210. That is ~5% of variance — the extreme cases (deepest MSAs are ~240 aa,
shallowest are 2,000-4,000 aa) suggest a stronger relationship than actually
holds. Model depth and crop status as separate features; do not treat one as a
proxy for the other.

12.2% of genes exceed the ESM-1v context and will be cropped. One alignment is
degenerate (`TENX_HUMAN`, depth 1 — no homologs found, 4,244 aa); EVE cannot
have fit a meaningful model there. Filter it.

## ESM-1v — HuggingFace cache (not in this directory)

Five-model ensemble at `~/.cache/huggingface/hub/models--facebook--esm1v_t33_650M_UR90S_{1..5}`,
~2.4 GB each. The only predictor in this project that is actually executed.

Model 1 was already present from prior work.

**Duplicate revisions removed 2026-09-01 (12 GB reclaimed; cache 27 GB -> 15 GB).**
All five models had fetched two revisions, not just model 1 as previously recorded
here: `pytorch_model.bin` on `main`, plus a PR-branch snapshot holding *only*
`model.safetensors` — no config, no tokenizer, so `from_pretrained` could never
resolve to it. Each model is now 2.4 GB.

Verified after deletion: model 1 reproduces the documented ubiquitin `I44V`
LLR of **-6.27** exactly. The `esm.contact_head.regression.*` MISSING warning on
load is expected — the contact head is absent from the checkpoint and unused here.

If a future `from_pretrained` re-fetches the safetensors branch, it is safe to
drop again; prefer `huggingface-cli delete-cache` over deleting blobs by hand.

Scoring: `scripts/esm1v_masked_marginals.py`. Emits the full 20-AA log-prob
vector per position, so every substitution at a residue costs one forward pass
total. Compute LLR at analysis time as `lp[mut] - lp[wt]`.

## Environment

Use the `vep` conda env, **not `ml_dev`**. `ml_dev` has a broken torch (conda
`pytorch 2.5.1` and pip `torch 2.6.0` both installed; the pip `_C.so` links
against the conda `libtorch_cpu.dylib` and fails with a missing
`__ZN2at3cpu20is_arm_sve_supportedEv` symbol). `ml_dev` also carries
`esm 3.2.1`, which is EvolutionaryScale's ESM-3 and shadows `fair-esm` — another
reason to stay out of it. We load ESM-1v through `transformers`, so neither
package is needed.

    conda create -y -n vep python=3.11
    ~/miniconda3/envs/vep/bin/pip install torch transformers duckdb pyarrow \
        pandas biotite biopython scikit-learn lightgbm

Verified: torch 2.13.0 (MPS available), transformers 5.16.1, duckdb 1.5.5.

## Validation performed

All five ensemble models scored human ubiquitin (76 aa) at F4, I44, G76.
Per-model mean LLR over the 19 non-WT substitutions:

| site | models 1-5 | ensemble | sd |
|---|---|---|---|
| F4  | -8.93, -9.76, -8.96, -11.26, -9.28  | -9.64  | 0.87 |
| I44 | -12.50, -10.82, -8.40, -12.04, -10.32 | -10.81 | 1.44 |
| G76 | -8.43, -8.85, -7.84, -7.79, -7.87   | -8.16  | 0.42 |

Biologically correct: I44 is the most constrained site, as expected for
ubiquitin's buried hydrophobic-patch residue, central to receptor recognition.
`I44V` (conservative hydrophobic) is the most tolerated substitution there at
-6.27 in model 1; `I44W` the worst at -15.53. All sites strongly negative, as
expected for one of the most conserved proteins known.

**Use all five models.** They are independent seeds and disagree materially: at
I44, model 1 alone gives -12.50 against the ensemble's -10.81. Ensemble spread
also varies 3.4x across just these three positions (sd 0.42 to 1.44), which is
why `esm1v_masked_marginals.py` writes per-model rows rather than collapsing to
a mean -- that spread is a candidate uncertainty feature for the disagreement
model.

### Throughput is NOT yet measured

The observed ~1.1-1.6 s/position is dominated by 2.4 GB model loads amortised
over only 3 positions. It is not a usable throughput figure and the working
estimate of 5-15 GPU-hours for a full run is **unvalidated**. Benchmark on a few
thousand positions before sizing any Hoffman2 array job.

## Phase 0 spike — results (2026-09-01)

Scripts: `scripts/spike/phase0_join.py`, `scripts/spike/phase0_assert.py`.
Genes: TP53, BRCA1, CFTR, PTEN, MLH1. Output: `data/spike/joined_full.parquet`.

**The assertion passes 100%.** On 1,168 joined ClinVar x AlphaMissense rows,
`am_ref == af_ref == clinvar_ref`, and position and alt allele agree too — 4/4
checks at 1,168/1,168. The join architecture is sound; Phase 1 can proceed on it.

### The trap that cost the most, and will recur

**Do not join AlphaMissense to MANE on `transcript_id`.** AM was built on an
older Ensembl release. Against MANE v1.5:

| gene | MANE v1.5 | AlphaMissense | |
|---|---|---|---|
| PTEN | ENST00000371953.8 | ENST00000371953.8 | match |
| BRCA1 | ENST00000357654.9 | ENST00000357654.8 | version drift |
| CFTR | ENST00000003084.11 | ENST00000003084.10 | version drift |
| MLH1 | ENST00000231790.8 | ENST00000231790.6 | version drift |
| TP53 | ENST00000269305.9 | ENST00000445888.6 | **different transcript** |

An equality join on `transcript_id` returned 194 rows from one gene (PTEN) and
silently dropped the rest — 2.9% survival, no error. `(chrom,pos,ref,alt)` is the
pivot; MANE is used only to *rank* transcripts when AM reports several at a
locus. Coordinate join: 1,168 rows, all five genes, 6x the rows.

TP53 is not a coverage gap — AM has all 2,569 of its missense SNVs, just under
the old Ensembl canonical rather than the MANE transcript.

### Join yield is expected, not lossy

1,168 / 6,677 ClinVar rows = 17.5%, which looks alarming and is not. Of those
6,677 two-star P/LP+B/LB SNVs only 1,185 are missense (985 nonsense, 4,507
splice/UTR/other). Against the AM-eligible denominator the join captures
**1,168 / 1,185 = 98.6%**.

### Corrections to the plan's acquisition table

- **AlphaFold DB is on v6, not v4.** The `..._v4.pdb` URLs now 404. Resolve the
  version per accession from `https://alphafold.ebi.ac.uk/api/prediction/<acc>`
  (`pdbUrl` / `latestVersion`) rather than hardcoding it.
- **ClinVar `variant_summary.txt.gz` is 422 MB**, not the estimated ~100 MB.
- MANE release in use: v1.5.

### Measured ESM-1v throughput — replaces the unvalidated estimate

Five-model ensemble over 734 unique (protein, position) pairs, M3 Pro MPS,
batch size 16, `scripts/esm1v_masked_marginals.py`:

| model | wall | per position |
|---|---|---|
| 1 | 365 s | 498 ms |
| 2 | 368 s | 501 ms |
| 3 | 367 s | 500 ms |
| 4 | 373 s | 508 ms |
| 5 | 414 s | 565 ms |
| **ensemble** | **1,887 s** | **2.571 s** |

Extrapolating to the plan's 30-60k unique positions for 300-500 genes:
**21-43 hours on this laptop's MPS.** The plan's 5-15 h figure was for an A100.
**MEASURED 2026-09-02 (Hoffman2 job 14625415): the A100 is 1.62x faster than MPS**,
giving 1.590 s/position and a revised budget of **13-27 A100-hours** -- above the
plan's 5-15 h. See `hoffman2/README.md`. The pilot also cross-validated the two
stacks: per-variant ensemble LLRs agree to Spearman rho 0.9999993 between laptop
and cluster, so scores from either are interchangeable for ranking.

Two things make this cohort *pessimistic* as a benchmark: 377 of 734 positions
come from proteins over the 1022-token cap (BRCA1, CFTR), and the mean scored
window is 771 tokens. A gene set with typical-length proteins will be faster.

### Spike results — the pipeline recovers known signal

| | pooled AUROC | macro per-gene AUROC |
|---|---|---|
| AlphaMissense | 0.974 | **0.972** |
| ESM-1v (5-model LLR) | 0.926 | **0.920** |

Per gene: PTEN is saturated (AM 0.997 / ESM 0.997), MLH1 hardest (0.921 / 0.837).

**Do not read these as generalisation estimates.** Three reasons, all structural:
five extensively-studied genes; a 2.9:1 P/LP-to-B/LB imbalance; and the ACMG
BS1/BA1 problem — a share of these benign labels are benign *because* they are
common, which is the frequency signal AlphaMissense was tuned on. Clean numbers
require the gene-held-out, macro-averaged, gene-ID-baselined protocol in Phase 2,
and the DMS replication in Phase 3.

Signed rank-normalised disagreement: mean -0.004, sd 0.225; 22.3% of variants
exceed |0.25|. So there is real disagreement to model at this cohort size.

**A confound to note now:** mean |disagreement| is 0.183 on cropped vs 0.158 on
uncropped rows, which looks like the crop artifact the plan predicts — but in
this cohort *cropped is perfectly confounded with gene identity* (the cropped set
is exactly BRCA1 + CFTR). It is not evidence of a crop effect. Testing that needs
a gene set containing both long and short proteins with the crop flag varying
within gene-length strata.

### EVE added to the spike (2026-09-03) — `scripts/spike/phase0_eve.py`

EVE is the predictor that breaks the design's three-way confound (MSA-based, no
structure, no frequency supervision). Output: `data/spike/spike_final_3pred.parquet`.

**Coverage: 1,149/1,168 = 98.4%**, wt residue consistent on every matched row.
Per gene: PTEN 100%, BRCA1 99.1%, CFTR 98.9%, MLH1 98.8%, TP53 95.8%.

#### Three join hazards, all of which return rows rather than erroring

1. **EVE is keyed by UniProt entry name, not accession and not gene symbol.**
   TP53's entry name is `P53_HUMAN` — the T is dropped — so matching on the gene
   symbol silently loses the gene entirely. Resolve accession -> entry name from
   the UniProt REST API (`/uniprotkb/<acc>.json?fields=id`), never by string
   manipulation. Same failure class as the MANE transcript trap.
2. **`BRCA1_HUMAN-checkpoint` is a duplicate** of `BRCA1_HUMAN`: 37,260 rows,
   verified byte-identical via `EXCEPT`. A training artifact in the release. A
   prefix/LIKE join double-counts BRCA1.
3. **BRCA1 ships domain-only models** (`BRCA1_RING_HUMAN`, `BRCA1_BRCT_HUMAN`)
   beside the full-length fit — the only gene in the release with domain splits.
   We take full-length; domain-fit position numbering is not guaranteed to match
   the full protein, and mixing them puts two coordinate systems in one column.

An exact `gene IN (...)` match on resolved entry names avoids all three, and the
script asserts one model per protein.

### CORRECTION to the earlier per-gene AUROC table

The earlier table reported per-gene AUROCs without checking how many negatives
each rested on. **PTEN's 0.997 rests on 6 benign variants and CFTR's 0.992 on 3.**
Those two are noise, and the macro-average weighted them equally with BRCA1's 345.

| gene | n | **n_neg** | AM | ESM-1v | EVE | |
|---|---|---|---|---|---|---|
| BRCA1 | 342 | 186 | 0.960 | 0.913 | **0.964** | |
| TP53 | 272 | 63 | 0.987 | 0.951 | 0.957 | |
| MLH1 | 166 | 29 | 0.922 | 0.836 | 0.871 | |
| PTEN | 194 | **6** | 0.997 | 0.997 | 0.992 | unreliable |
| CFTR | 175 | **3** | 0.992 | 0.893 | 0.994 | unreliable |

| macro | AM | ESM-1v | EVE |
|---|---|---|---|
| all 5 genes | 0.972 | 0.918 | 0.956 |
| **3 genes with >=10 negatives** | **0.956** | **0.900** | **0.930** |

Pooled AUROC on identical complete-case rows (n=1,149), 95% bootstrap CI over
2,000 resamples:

| predictor | pooled AUROC | 95% CI |
|---|---|---|
| AlphaMissense | 0.974 | [0.964, 0.982] |
| ESM-1v | 0.925 | [0.909, 0.939] |
| EVE | 0.896 | [0.874, 0.917] |

**EVE moves from 0.896 pooled to 0.956 macro** — the largest pooled-vs-macro gap
of the three, and a direct demonstration of the between-gene base-rate inflation
the design review warns about. On BRCA1, the one gene with a substantial number
of negatives, **EVE beats AlphaMissense** (0.964 vs 0.960). Report macro, never
pooled.

Rank agreement (within-protein Spearman, **variant-count weighted** across the
five proteins): AM-EVE +0.700, AM-ESM +0.688, ESM-EVE +0.621. No pair is
redundant; all three carry distinct signal. The weighting matters and was not
recorded originally -- an unweighted mean over proteins gives +0.664 / +0.668 /
+0.590, about 0.03 lower throughout. State the aggregation whenever a per-gene
statistic is averaged; it is the same pooled-vs-macro distinction that decides
the Phase 2 result below.

### EVE missingness is strongly non-random — do not drop those rows silently

Only 19 spike variants lack an EVE score, but they are not a random 19:

| | missing EVE | has EVE |
|---|---|---|
| mean pLDDT | **51.0** | 82.3 |
| fraction P/LP | **0.16** | 0.75 |

Missing rows are disordered *and* mostly benign. At project scale (EVE is 80.9%
complete overall) this is a large missing-not-at-random block that will bias any
complete-case analysis. `spike_final_3pred.parquet` carries an explicit
`eve_missing` flag; keep it in the model rather than dropping the rows.

### Precomputed-score audit (the plan's Phase 0 question)

ProteinGym **v1.3** ships precomputed zero-shot scores, and **ESM-1v is among the
scored models**:

| resource | size | file |
|---|---|---|
| DMS benchmark, substitutions | 1.0 GB | `DMS_ProteinGym_substitutions.zip` |
| zero-shot model scores, substitutions | 4.4 GB | `zero_shot_substitutions_scores.zip` |

Base URL: `https://marks.hms.harvard.edu/proteingym/ProteinGym_v1.3/<file>`

**Consequence for the compute plan: Phase 3 needs no GPU at all.** The DMS
adjudication experiment — the minimum publishable result — can be assembled from
ProteinGym's own ESM-1v scores plus AlphaMissense from the local atlas, joined on
UniProt + `protein_variant` via `am_aa_substitutions.parquet`. ProteinGym does
*not* ship AlphaMissense, which is exactly the gap our local copy fills.

GPU work remains only for the ClinVar arm (ESM-1v over 300-500 genes), which does
not overlap ProteinGym's assay set. Run the DMS phase first: it is the
load-bearing result and it is now pure data engineering.

Provenance: **metadata-only** — read from the ProteinGym GitHub README, not from
downloading the archives. Verify sizes and the v1.3 URL before committing quota.

### Spike cohort characteristics

| | |
|---|---|
| variants / unique (protein,position) | 1,168 / **734** — 1.59 variants per forward pass |
| label balance | 865 P/LP vs 303 B/LB (2.9:1) |
| cropped (>1022 aa) | 522 rows — BRCA1 (1,863 aa), CFTR (1,480 aa) |
| mean pLDDT | PTEN 83.0, MLH1 77.3, CFTR 75.6, TP53 75.1, **BRCA1 41.6** |

BRCA1 at mean pLDDT 41.6 is largely disordered. With 246 of the 734 positions,
it alone is a third of the cohort — structural stratification on this five-gene
set will be dominated by one poorly-ordered protein. Widen the gene set before
reading anything into buried/exposed contrasts.

## Acquisition status

Acquired and in use:

- ClinVar `variant_summary.txt.gz` (labels, ≥2-star) — 422 MB
- MANE Select summary — v1.5, used only to *rank* transcripts
- AlphaMissense hg38 + aa_substitutions parquet
- EVE `eve_variants.parquet` + `eve_msa_depth.parquet`
- AlphaFold structures — fetched per gene by `scripts/trial/trial_01_genes.py`,
  **v6, resolved per accession from the API**. Not the bulk proteome tarball;
  the plan's "human proteome v4" is superseded (the _v4 URLs 404).

Still to acquire:

- ProteinGym substitutions (DMS ground truth, ships ESM-1v + AM scores).
  Phase 3 depends on this and it needs no GPU — highest-value next download.
- PrimateAI-3D (**license registration required — start early**, it is the
  long pole on a full Phase 1)
- ESM-IF1 scores (via ProteinGym) — fills the structure-only cell
- gnomAD constraint metrics (LOEUF, gene-level; ~10 MB, do NOT pull gnomAD VCFs)

## Trial Phase 1 + 2 (2026-09-12) — 150-gene scale-up

A deliberately undersized rehearsal of Phases 1 and 2, run end to end on the
laptop in about three minutes. Scripts: `scripts/trial/trial_0{1,2,3}_*.py`.
Outputs: `data/trial/`. **ESM-1v is not in this table** — at the measured
2.571 s/position it would need ~5.5 h for this cohort's 7,673 unique positions,
which is the Hoffman2 array's job, not the laptop's. Phase 2 is a structure-only
baseline and does not need it.

### Gene selection

Selection is by ClinVar benign count, because benign is the scarce class and a
per-gene AUROC without negatives is noise (the Phase 0 spike's PTEN 0.997 rested
on 6). ClinVar is far richer than the spike implied: **2,906 genes have >=10
two-star benign variants**, so ClinVar is not the binding constraint on gene-set
size. AlphaFold availability is.

| stage | genes |
|---|---|
| candidates (>=10 of both classes in ClinVar) | 360 |
| resolved to a modal UniProt accession via AM | 356 |
| lost: AlphaFold returned **isoform-only** entries | 27 |
| lost: AlphaFold **404**, no prediction at all | 33 |
| selected for the trial | 150 |

### Two new AlphaFold traps, both silent

1. **`/api/prediction/<acc>` can return only ISOFORM entries.** For NF1 the
   endpoint returns `P21359-5` (593 aa), `P21359-4` (1,598 aa) and `P21359-3`
   (551 aa) — and never canonical `P21359`, which is 2,839 aa. Taking `j[0]`
   hands you a 593-residue structure while ClinVar and AlphaMissense positions
   are on the 2,839-residue canonical, so every structural feature attaches to
   the wrong residue with no error. **Accept only an entry whose
   `uniprotAccession` equals the bare accession, with no `-N` suffix.**
2. **Very long proteins 404 rather than returning fragments.** TTN (Q8WZ42,
   34,350 aa) has no prediction at this endpoint at all. A `prot_len <= 2700`
   filter intended to avoid AlphaFold's F1/F2 fragment splitting was therefore a
   no-op — the long proteins had already been dropped as "API errors". Do not
   mistake that for the filter working. The full run must handle fragments.

### The assertion gate has to be split into hard and soft

At five genes the reference-AA assertion passed 1,168/1,168 and looked binary.
At 150 genes it drops to 93.6% — and most of that is not corruption.

| check | rows | |
|---|---|---|
| AM ref == AlphaFold ref | 9,600 / 9,768 | 98.28% |
| AM ref == ClinVar ref | 9,759 / 9,768 | 99.91% |
| AM alt == ClinVar alt | 9,759 / 9,768 | 99.91% |
| **all three (the hard gate)** | **9,595 / 9,768** | **98.23%** |
| AM pos == ClinVar pos | 9,147 / 9,768 | 93.64% |

The position check fails on 449 rows **where all three amino acids agree**. The
offset is constant within a gene — MECP2 −12, RUNX1 −27, MEN1 +5, MBD5 −233,
SPTAN1 −5 (constant in 9 of 12 affected genes) — because ClinVar's `Name` field
cites a RefSeq isoform whose protein numbering differs from the UniProt
canonical that AM and AlphaFold use. MECP2's e1/e2 isoforms differ by exactly 12
N-terminal residues. The variant's identity is already pinned by the
`(chrom,pos,ref,alt)` coordinate join, so these rows are sound.

**Treatment: hard-gate on residue identity, soft-flag the offset.** Dropping on
position would have discarded MECP2 and RUNX1 almost entirely for a naming
convention. The 173 genuine hard failures (ARID1B 92, MEN1 55, ALPK3 23) are
rows where AM and AlphaFold disagree on the residue itself — real frame
mismatch — and are dropped and logged. `cv_frame_offset` is kept on every row.

### Trial cohort

| | |
|---|---|
| variants | 9,595 |
| genes | 148 |
| label balance | 5,757 P/LP vs 3,838 B/LB (**1.50:1**, vs the spike's 2.85:1) |
| genes with >=10 of both classes | 54 |
| EVE coverage | 8,306 / 9,595 = 86.6% (135 of 148 proteins in EVE) |
| unique (protein, position) | 7,673 — the ESM-1v cost unit |

### Phase 2 result — the gene-ID null is the whole story

Structure-only features (pLDDT, RSA, SASA, weighted contact number, relative
position, protein length), `HistGradientBoostingClassifier`, 5-fold. Macro =
mean per-gene AUROC over the 54 evaluable genes; CI = 2,000-sample bootstrap
**over genes**, not variants.

| model | macro AUROC | 95% CI | pooled |
|---|---|---|---|
| Structure-only, **gene-held-out** | **0.844** | [0.816, 0.870] | 0.822 |
| Structure-only, random split | 0.908 | [0.888, 0.928] | 0.956 |
| Gene-ID-only (circular) | 0.500 | — | **0.821** |
| Gene-ID-only, gene-held-out | 0.500 | — | 0.500 |
| AlphaMissense *(reference)* | 0.966 | [0.948, 0.980] | 0.970 |
| EVE *(reference)* | 0.948 | [0.926, 0.965] | 0.931 |

**Pooled and macro tell opposite stories, and that is the finding.** Knowing
nothing but the gene name gives pooled AUROC 0.821; the structure model's pooled
AUROC is 0.822 — **+0.001 over the null**. Macro-averaged, the same model is
0.844 against a macro null of 0.500. The design review predicted this shape
("if gene-ID-only gets 0.70 and your structural model gets 0.74, you have
learned almost nothing"); the real gap is far tighter than that hypothetical.
Pooled AUROC credits the model for between-gene base rates it never learned.
**Report macro. Always.**

Gene-ID macro is 0.500 *by construction* — a within-gene constant cannot rank.
A cross-validated gene-prevalence baseline scores slightly *below* 0.5 (0.367)
for a mechanical reason worth knowing: its score varies only across CV folds,
and a gene's train-prevalence excludes the test fold, so a fold holding more
pathogenic variants leaves a lower train prevalence — systematic
anti-correlation. Verified as 5 distinct scores per gene, one per fold, with
mean label decreasing monotonically in the score. It is an artifact; the null
is 0.5.

Circularity costs the structure model **0.064 macro AUROC** (0.908 random ->
0.844 gene-held-out). That is the memorisation the review warned about, measured.

AM and EVE are shown as references only and are **not** on equal terms — neither
trained on ClinVar, while the structure model did. That asymmetry is why the
gene-ID null exists; the numbers alone do not resolve it.

**Feature importance is deliberately not reported.** pLDDT, RSA and WCN are
collinear by construction, and permutation importance on correlated features
measures unique non-redundant contribution within one fitted model, not
attribution. That is Phase 4's job and needs grouped conditional permutation
plus LOCO refitting.

### What this trial does not establish

Still five-gene problems at 150-gene scale: no DMS replication, no frequency
stratification, no paralog clustering (gene-held-out splits still leak across
paralogs — mmseqs at 30% identity is the fix and is not applied here), and no
ESM-1v, so no disagreement analysis at all.

### Scaled to 400 genes (same day) — and one 150-gene result does not replicate

The 150-gene run above was rerun at 400 genes with no code changes, to check
which findings were real and which were small-sample. Runtime: gene selection
35 s, table build 2 min 57 s, Phase 2 baseline 5 s.

The candidate pool is the ceiling: **827 genes have >=10 of both classes** in
two-star ClinVar, and after the ~17% AlphaFold loss (isoform-only + 404) that
allows roughly 690. 400 is comfortable; much beyond ~690 is not reachable on
this label definition without relaxing the two-star or >=10-of-both criteria.

| | 150 genes | 400 genes |
|---|---|---|
| variants | 9,595 | 18,700 |
| genes in table | 148 | 398 |
| label balance | 1.50:1 | 1.62:1 |
| evaluable genes (>=10 of both) | 54 | **103** |
| EVE coverage | 86.6% | 84.9% |
| unique (protein, position) | 7,673 | **15,367** |
| hard-gate failures | 173 | 193 |
| isoform-offset rows kept | 449 | 951 |

**ESM-1v budget for this cohort:** 15,367 unique positions x 5 models is
**6.8 A100-hours** at the measured 1.590 s/position, or 11.0 h on the laptop.
That is the number to size the Hoffman2 array against, and it is why ESM-1v is
absent from this table.

| model | macro AUROC | 95% CI | pooled |
|---|---|---|---|
| Structure-only, **gene-held-out** | **0.862** | [0.838, 0.884] | 0.839 |
| Structure-only, random split | 0.887 | [0.868, 0.904] | 0.940 |
| Gene-ID-only (circular) | 0.500 | — | 0.805 |
| AlphaMissense *(reference)* | 0.972 | [0.962, 0.981] | 0.966 |
| EVE *(reference)* | 0.951 | [0.938, 0.962] | 0.934 |

**CORRECTION to the 150-gene write-up.** That run reported the structure model
beating the gene-ID null by **+0.001** pooled (0.822 vs 0.821) and presented it
as a stark confirmation of the review's warning. It does not replicate. At 400
genes the gap is **+0.034** (0.839 vs 0.805) — still small, but 34x larger, and
the +0.001 was substantially small-sample noise on 54 evaluable genes. Treat the
150-gene figure as superseded.

Two other quantities moved in the same direction, and for the same reason:

- The circularity penalty **halved**, from −0.064 to −0.025 (random 0.887 ->
  gene-held-out 0.862). With 398 genes there is far less to gain from memorising
  any one gene's hotspot domain, so holding genes out costs less.
- The gene-held-out macro AUROC **rose**, 0.844 -> 0.862, with a tighter CI
  ([0.816, 0.870] -> [0.838, 0.884]) on 103 evaluable genes instead of 54.

What survives the scale-up unchanged: pooled and macro still tell opposite
stories (gene-ID-only reaches pooled 0.805 while its macro is 0.500 by
construction), so **report macro**. AM 0.972 and EVE 0.951 macro remain far
above the structure-only model, and remain non-comparable — neither trained on
ClinVar.

The lesson generalises beyond this table: **effect sizes estimated on a few
dozen evaluable genes are not stable**, which is exactly the confidence-interval
argument the design review made for 300-500 genes. Anything measured on the
five-gene spike or on 150 genes should be re-estimated at full scale before it
is written down as a finding.
