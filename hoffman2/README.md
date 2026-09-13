# Hoffman2 setup — CAIR Summer VEP project

Account `ijeon`, `$SCRATCH=/u/scratch/i/ijeon`. Adapted from
`Reference_Repos/condenseq-esm-handoff/hoffman2/` — read `HOWTO_GPU_DISPATCH.md`
there before changing any submit line.

## Queue policy for this project

**Public queue only. Do not add `highp`.** `ijeon` *is* in the `kappel` group, so
`highp` would dispatch to the group-owned node g20 — but this is CAIR work, not
Kappel lab work, so it does not belong on the lab's owned allocation. Teammates
on their own accounts are not in that ACL anyway, and a non-member's `highp` job
waits forever with no error.

## The two dispatch rules that cost days last time

1. **Never pass `-q` on a GPU job.** A queue list combined with a generic
   `-l gpu` is a conflicting request; the job sits in `qw` indefinitely with no
   error message. (IDRE ticket #100868426.) Request the GPU through `-l`.
2. **`cuda=N` is the RSMAP consumable `qw` actually waits on** — not CPU slots.
   `qstat -g c` shows free CPU slots and will mislead you. Check GPUs with:

   ```bash
   ssh hoffman2 "qhost -F gpu_model,cuda | awk '/^g/{h=\$1} /gpu_model=/{m=\$0} /hc:cuda=[1-9]/{print h,m,\$0}'"
   ```

## Environment

`$HOME/venvs/vep` — a conda env, built by `setup_vep_env.sh` as a **clone of the
condenseq env** rather than a fresh install, because that env already carries a
working `torch 2.6.0+cu118` and re-resolving the CUDA build is the risky part.

Verified: torch 2.6.0+cu118 (cuda 11.8), transformers 4.57.6, protobuf 3.20.3.

**Two traps found building it (2026-09-02):**

- `conda create --clone` leaves a **broken pip** — both `bin/pip` and
  `python -m pip` fail with `ImportError: cannot import name 'get_runnable_pip'`.
  Fix with `conda install -p $HOME/venvs/vep pip protobuf` (see
  `fix_vep_env.sh`), which bypasses pip entirely.
- `protobuf` is required by `EsmTokenizer.from_pretrained` and **the error
  message does not say so**. Install it. Do NOT install `sentencepiece` — no
  wheel builds on this cluster.

## `SGE_TASK_ID` is the string "undefined" on non-array jobs

Submitting without `-t` does **not** leave `SGE_TASK_ID` unset — SGE sets it to
the literal string `undefined`. So `${SGE_TASK_ID:-1}` yields `undefined`, not
`1`, and anything typed (`argparse type=int`) rejects it:

    esm1v_masked_marginals.py: error: argument --shard: invalid int value: 'undefined'

Job 14625404 died on this within seconds of dispatch. `esm1v_score.qsub`
normalises it explicitly:

```bash
SHARD="${SGE_TASK_ID:-1}"
[ "$SHARD" = "undefined" ] && SHARD=1
```

Worth knowing because it only bites the single-task pilot form, never the array
form — so it survives testing with `-t` and fails the first time someone runs one
task.

## Model cache

`$SCRATCH/hf_cache/hub` — the five ESM-1v checkpoints are **pushed from the
laptop by rsync**, never downloaded on the login node (it OOMs on the
downloader). Jobs run with `HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1`, so a
missing checkpoint is a hard failure rather than a silent download.

Note: macOS ships rsync 2.6.9, which does **not** support `--info=progress2`.
Use plain `rsync -a`.

## Submitting

```bash
# pilot: one task, all positions, 5-model ensemble (measures the A100 rate)
ssh hoffman2 'cd $SCRATCH/vep && qsub -l gpu,A100,cuda=1 esm1v_score.qsub'

# full run: shard by protein; the lockfile de-dups if you broaden to a
# second array on another GPU model
ssh hoffman2 'cd $SCRATCH/vep && qsub -t 1-20 -v N_SHARDS=20 -l gpu,A100,cuda=1 esm1v_score.qsub'
```

Diagnostics: `qstat -u ijeon`; `qalter -w v <jobid>` answers "will this job ever
run?" — `no suitable queues` means a disabled queue or an ACL you are not in.

## Pilot result — job 14625415 (2026-09-02)

734 positions x 5 models, one task, public A100-SXM4-40GB on g15. Queued ~2 min.

| model | wall | per position |
|---|---|---|
| 1 | 265 s | 361 ms (includes warmup) |
| 2 | 208 s | 284 ms |
| 3 | 266 s | 362 ms |
| 4 | 199 s | 270 ms |
| 5 | 229 s | 312 ms |
| **ensemble** | **1,167 s** | **1.590 s** |

**The A100 is only 1.62x faster than the M3 Pro laptop** (1,167 s vs 1,887 s).
That is far less than the 3-4x the plan's 5-15 GPU-hour estimate implicitly
assumed.

Revised budget for 30-60k unique positions: **13-27 A100-hours**, against the
plan's 5-15 h. As a `-t 1-20` array that is 40-80 min wall, so it is not a
schedule risk — but it is the real GPU-hour cost, and it is up to ~1.8x the
planned figure at the high end.

### Cross-validation: the cluster stack produces the laptop's numbers

Same 734 positions, different hardware *and* different transformers (4.57.6 on
the cluster vs 5.16.1 locally):

| | |
|---|---|
| per-variant ensemble LLR, max abs diff | 0.0092 |
| mean abs diff | 0.00168 |
| Spearman rho | 0.99999926 |
| max percentile-rank shift | 0.0034 |

Float32 kernel nondeterminism only. Scores from either machine are
interchangeable for ranking, which is what the disagreement analysis uses.

### Why only 1.62x — worth testing before the full run

**Not padding.** Batches are length-homogeneous by construction: positions are
grouped by protein, and within a protein every crop is the same length (the full
sequence, or exactly 1022).

**Most likely fp32 without TF32.** PyTorch defaults `allow_tf32=False` for matmul,
so an A100 runs this at ~19.5 TFLOPS instead of ~156. Worth testing as a one-line
change before committing to the full array:

```python
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True
```

Also untested: raising `--batch-size` above 16 (40 GB is not close to full at
650M / 1022 tokens), and bf16 autocast. **This is a hypothesis, not a
measurement** — benchmark it on the same 734 positions before believing it, and
re-run the cross-validation above if you change precision, since TF32 will move
the log-probs more than the 1e-3 seen here.
