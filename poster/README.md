# Poster figures

Four figures for a poster, built the same way as `Lab_Kappel/symposium_2026/poster/`:
every plotted number is re-derived from `data/trial/` or `data/spike/` and
asserted against the value recorded in `data/README.md` before the figure is
saved (`make_figs.py`). Categorical comparisons use a fixed blue/orange/gray
set; the one heatmap uses viridis. No red-against-green anywhere.

## Rebuild

    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ~/miniconda3/envs/vep/bin/python poster/make_figs.py

Needs the `vep` conda env (duckdb, pandas, scikit-learn, matplotlib) and
`data/trial/`, `data/spike/` populated (see the project's `data/README.md`).

## Figures and captions

**Figure 1 — `fig1_gene_id_null.png`.** Knowing only the gene name pools to
AUROC 0.805; the structure-only model, gene-held-out, adds only +0.034 pooled
(0.839), while the same model is 0.862 macro against a macro null of 0.500.
Pooled and macro tell opposite stories because pooled credits the model for
between-gene base rates it never learned. AlphaMissense and EVE are shown as
references, not on equal terms — neither trained on ClinVar labels.
*(400-gene trial, 103 evaluable genes, `data/trial/phase2_results.csv`.)*

**Figure 2 — `fig2_replication.png`.** The 150-gene run's gene-ID-null gap
(+0.001 pooled) does not replicate at 400 genes (+0.034) — effect sizes on a
few dozen evaluable genes are not stable. The gene-held-out macro AUROC itself
rises with scale (0.844 → 0.862) and its bootstrap CI tightens (54 → 103
evaluable genes). *(Left: pooled AUROC minus gene-ID-null pooled AUROC. Right:
structure-only gene-held-out macro AUROC, 95% CI bootstrapped over genes.)*

**Figure 3 — `fig3_rank_agreement.png`.** No pair of predictors is redundant:
variant-count-weighted within-protein Spearman rank agreement is 0.700
(AlphaMissense–EVE), 0.696 (AlphaMissense–ESM-1v) and 0.622 (ESM-1v–EVE).
*(Phase 0 spike, 5 genes, 1,149 complete-case variants.)*

**Figure 4 — `fig4_pergene_auroc.png`.** PTEN and CFTR look saturated
(AUROC ≈ 0.99–1.00 for all three predictors) but rest on only 6 and 3 benign
variants respectively — open markers flag them as unreliable. On BRCA1, the
gene with the most negatives (186), EVE edges out AlphaMissense (0.964 vs
0.960). *(Phase 0 spike, per-gene AUROC, complete-case rows only; ESM-1v LLR
is sign-flipped since more negative = more damaging.)*

## Notes

- Figure 4's AUROC values are on the 1,149-row complete-case subset (all three
  predictors scored), matching `data/README.md`'s per-gene table — not the
  full 1,168-row spike, which pads a few genes' AUROC by a few thousandths.
- Figure 2's "34x larger" language in `data/README.md` is computed on the
  rounded +0.001/+0.034 values; the full-precision ratio is ~27x. The figure
  states the rounded deltas rather than repeating the ratio.
- No ESM-1v in Figures 1–2: the 400-gene Phase 2 baseline is structure-only by
  design (ESM-1v needs the Hoffman2 array, see the top-level README).
