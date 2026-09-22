# Batch-confound audit for public transcriptomic datasets

Reference implementation, derived results and figures for the manuscript
**"A diagnostic audit for batch and confounding artifacts in reuse of public
bulk transcriptomic data: a two-domain study in perioperative neurocognitive
disorder and traumatic brain injury"** (under review).

[![Figure contract](https://img.shields.io/badge/figures-%E2%89%A47.2in%20%7C%20%E2%89%A57pt%20%7C%20300dpi-blue)](#figure-specification)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.13+](https://img.shields.io/badge/python-3.13%2B-blue.svg)](requirements.txt)

**Source:** <https://github.com/eagle-wk/pocd-batch-audit> ·
**Archived release (cite this):** [10.5281/zenodo.22899702](https://doi.org/10.5281/zenodo.22899702)

---

## What this is

A public dataset's group labels are usually taken at face value. If the
treatment arm was also the batch that ran on a different day, every downstream
number — DEG lists, module–trait correlations, immune deconvolution — is
confounded, and no amount of downstream statistics detects it.

This repository implements the audit we propose instead: **four defences applied
in order, ending in one of four verdicts**, and it ships the result of running
that audit over 14 contrast-units drawn from 12 GEO series in two disease
domains.

The headline finding is negative and stays negative:

> Of 14 contrast-units, **2 are unconditionally analysable, 6 only under
> documented restrictions, and 6 are categorically unusable**. The initial
> biological hypothesis (an m6A × programmed-cell-death signature) did **not**
> replicate across datasets — direction agreement was 11/21 = 52.4 %
> (p = 1.000), i.e. chance. The one entity that survived every test is a
> down-regulation of interferon *signalling* components, and even there only the
> sensing/transduction layer moves.

---

## Quick start

No R, no Bioconductor, no network access required to verify the shipped results.

```bash
git clone https://github.com/eagle-wk/pocd-batch-audit.git
cd pocd-batch-audit
pip install -r requirements.txt

# 1. Check every figure against the printed-size contract (width, font, dpi, PDF)
python checks/verify_figures.py

# 2. Recompute the manuscript's headline numbers from the shipped derived results
python checks/verify_numbers.py

# 3. Regenerate all six figures from the shipped intermediate results
cd figures
python make_fig1_fig2_cases.py     # Fig 1, 2, 5
python make_fig3.py                # Fig 3
python make_fig4_fig6.py           # Fig 4, 6
```

Step 3 re-creates byte-identical PNGs to the ones in `figures/output/` — that is
the reproducibility claim, and it is checkable in one command. If any figure
differs, something in the environment has drifted.

---

## The four defences

| Stage | Question | Kills | Implemented in |
|---|---|---|---|
| **0** | Are the group labels independently verifiable? | Datasets with no usable labels | `audit/audit_one.py` |
| **1** | Is the batch real, and is it entangled with the group? | Structural confounding; technically indistinguishable blocks; nesting | `audit/batch_audit_generic.py` |
| **2** | Does the signal survive batch correction without moving globally? | Global drift presented as a finding | `analysis/geneset_shift_diagnostic.py` |
| **3** | Is the module's identity asserted or tested externally? | Circular identities read off hub genes | `analysis/module_geneset_overlap.py` |
| **4** | Does the direction replicate across independent datasets? | Single-dataset artefacts | `analysis/cross_dataset_meta.py` |

Stage 1 is the operative one. It applies **three gates that must all fire**:

1. **Structural** — how many groups sit 100 % inside one batch level.
2. **Technical separability** — can the batch be told apart from technical
   covariates alone? Leave-one-out accuracy must beat the random baseline by
   `gain ≥ 0.15` **and** a label-permutation test must give `p < 0.05`.
   Accuracy alone is not enough; a block label is not evidence of a real batch.
3. **Nesting** — a block difference that is fully explained by a stratified
   variable (time point, age) is not independent batch evidence.

The four verdicts are `stop` / `correct` / `proceed` / `normalize`, assigned per
**contrast-unit** — not per GEO series. One series can contain both a usable and
an unusable comparison.

---

## Layout

```
audit/        the audit engine            → Stage 0-1
analysis/     post-audit analysis         → Stage 2-4
benchmark/    simulation benchmark + tool comparison (ComBat / sva / RUV)
figures/      figure generators + figstyle.py (the size contract)
checks/       pre-submission verifiers
data/genesets/  pre-defined gene sets (literature-derived, not data-derived)
results/      derived result tables the manuscript cites
  derived/    small intermediates (verdict JSONs, benchmark JSONs)
docs/         data sources, pipeline notes, limitations, Zenodo guide
```

`figures/output/` holds the six published figures as PNG + PDF.

---

## Script → figure / result map

| Output | Script | Reads |
|---|---|---|
| Fig 1 — audit framework | `figures/make_fig1_fig2_cases.py` | schematic, no data |
| Fig 2 — simulation benchmark | `figures/make_fig1_fig2_cases.py` | `results/derived/benchmark/benchmark_sim.json` |
| Fig 3 — ComBat / sva / RUV | `figures/make_fig3.py` | `results/derived/benchmark/benchmark_tool_compare.json` |
| Fig 4 — 14-unit audit map | `figures/make_fig4_fig6.py` | `results/master_audit_table.tsv`, `results/derived/verdicts/*.json` |
| Fig 5 — three real cases | `figures/make_fig1_fig2_cases.py` | `master_audit_table.tsv`, `GSE304956_cross_batch_replication.csv` |
| Fig 6 — separability gate | `figures/make_fig4_fig6.py` | `results/derived/verdicts/*.json` |
| Master audit table | `audit/audit_one.py --verdict-dir ...` | per-unit verdict JSONs |
| Simulation benchmark | `benchmark/benchmark_audit_sim.py` | synthetic, seeded |
| Tool comparison | `benchmark/benchmark_tool_compare.py` | synthetic, seeded |

---

## Reproducing the benchmarks

Both benchmarks are self-contained and seeded (`--seed 20260918`), so they run
in a few minutes and reproduce exactly:

```bash
cd benchmark
python benchmark_audit_sim.py   --reps 100 --perm 200 --outdir ../results/derived/benchmark
python benchmark_tool_compare.py --reps 100 --perm 200 --outdir ../results/derived/benchmark
```

`benchmark_audit_sim.py` contains an assertion pinning the overall accuracy to
0.994; the source file it expects is the one cited in the manuscript, and the
assertion makes a wrong pick fail loudly rather than silently producing a figure
that disagrees with the text.

Re-running the audit itself against the **real** data needs the GEO downloads
(~527 MB, not shipped). See [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md):

```bash
python audit/audit_one.py --input <series_matrix.txt.gz> --outroot <workdir>
```

---

## Results files

Everything the manuscript cites lives in `results/`:

- `master_audit_table.tsv` — one row per contrast-unit: sample and gene counts,
  batch provenance, confounded groups, both separability statistics, verdict.
  **This is the authoritative table.**
- `audit_summary.json` — the 2 / 6 / 6 distortion split.
- `final_numbers.json` — every headline number in the manuscript, with its
  provenance, so text and data can be diffed mechanically.
- `derived/verdicts/*.json` — the raw per-unit verdict records.
- `derived/` — benchmark JSONs, cross-dataset meta-analysis, and the GSE304956
  evidence tables (cross-batch replication, PCA attribution, variance
  decomposition).

**On two statistics that intentionally differ.** `results/derived/cross_dataset_meta.tsv`
reports `p_bh` computed from the **one-sided** Fisher input (`0.0046` for
`IFN_signaling`), whereas the manuscript quotes the **two-sided** figure
(`BH = 0.0079`), because it is the more conservative of the two. Both are
correct; see [`results/README.md`](results/README.md). `final_numbers.json` is
the authority for anything quoted in the text.

---

## Figure specification

Every figure is generated through `figures/figstyle.py`, which **asserts** the
printing contract at save time and exits non-zero on violation:

- physical width ≤ **7.2 in** (double-column, full width)
- smallest text ≥ **7 pt**
- raster ≥ **300 dpi**
- PDF alongside each PNG, with TrueType fonts embedded (`pdf.fonttype = 42`)

Measured, not assumed: `save()` re-reads the written PNG's resolution and walks
the figure's text artists for the smallest visible point size. A layout that
only looks fine on screen fails at generation time rather than at proof stage.

---

## Data availability

All datasets are public on GEO; accessions, platforms and sample counts are
tabulated in [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md). No new sequencing
was performed and no restricted-access data was used.

## Citation

See [`CITATION.cff`](CITATION.cff). If you use the audit framework, please cite
the manuscript; if you use the code directly, the archived release is at
`https://doi.org/10.5281/zenodo.22899702`. The development repository is
<https://github.com/eagle-wk/pocd-batch-audit>.

## License

MIT — see [`LICENSE`](LICENSE).

---

## Releasing a new version

The DOI quoted above is the **concept DOI**: it always resolves to the newest
version of this archive. To cut a new version, follow
[`docs/ZENODO_RELEASE_GUIDE.md`](docs/ZENODO_RELEASE_GUIDE.md).
