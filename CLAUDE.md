# CAIR_Summer

Summer 2026 CAIR project: variant-effect predictor disagreement
(AlphaMissense vs ESM-1v vs EVE). See `README.md` for scope and current status.

Predecessor: the Winter 2026 ECG work (shared notebooks, `utils.py`) was kept at
`_archive/CAIR_W2026/`. **That directory is not present in this checkout** — if
you need it, find it before citing it. Either way it is a separate project; do
not treat those notebooks as the current codebase.

## Layout

| path | what |
|---|---|
| `scripts/` | `esm1v_masked_marginals.py` (the only predictor actually executed) |
| `scripts/spike/` | Phase 0 pipeline: `phase0_join` -> `_assert` -> `_finalize` -> `_eve` |
| `hoffman2/` | SGE submit script + cluster env setup for the A100 run |
| `data/` | Untracked (~3.6 GB). `data/README.md` IS tracked and is the lab notebook |
| `docs/` | Progress report, design review |

## Conventions

- Launch Claude from this folder, not from `Coding_Project`.
- Shared resources live one level up: `../Reference_Repos/`, `../tools/`.
- Working preferences that apply everywhere are in `~/.claude/CLAUDE.md`.
- Use the `vep` conda env, **not `ml_dev`** (broken torch; `esm` shadows
  `fair-esm`). See `data/README.md`.
- Query the large parquets with DuckDB, never pandas.
- Findings, row counts and traps go in `data/README.md` as they are discovered.
  That file is the project's memory — keep it current in the same commit as the
  code that produced the finding.
