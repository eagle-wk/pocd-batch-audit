# results/ — what each file is, and which number is authoritative

Short version: **`final_numbers.json` is the authority for anything quoted in
the manuscript.** `master_audit_table.tsv` is the authority for per-unit
verdicts. Everything else is supporting evidence or an intermediate.

---

## Files

| File | Content | Authority |
|---|---|---|
| `master_audit_table.tsv` | One row per contrast-unit (14 rows): sample/gene counts, batch provenance, confounded groups, separability statistics, verdict, free-text recommendation | **per-unit verdicts** |
| `audit_summary.json` | The three-way distortion split and the headline sentence | derived from the table |
| `final_numbers.json` | Every headline number cited in the manuscript, with provenance | **all quoted numbers** |
| `master_audit_report.md` | Human-readable pass over the same table | — |
| `derived/verdicts/*.json` | Raw per-unit audit records (13 files; `GSE276942` is a design-derived case documented in the table instead) | source for Fig 4 / Fig 6 |
| `derived/benchmark/benchmark_sim.json` | Simulation benchmark, 5 scenarios × 100 reps | source for Fig 2 |
| `derived/benchmark/benchmark_tool_compare.json` | ComBat / sva / RUV head-to-head | source for Fig 3 |
| `derived/GSE304956_cross_batch_replication.csv` | Per-contrast cross-batch replication | source for Fig 5B |
| `derived/GSE304956_*.csv` | Confound crosstab, PCA attribution, variance decomposition | supporting evidence |
| `derived/cross_dataset_meta.tsv` | 21 sets × 2 datasets direction consistency | source for the concordance finding |

---

## The distortion split: 2 / 6 / 6

Counted from the `cat` column of `master_audit_table.tsv`:

| Verdict | *n* | Units |
|---|---:|---|
| `usable` | **2** | GSE80174_Cortex, GSE80174_Thalamus |
| `caveat` | **6** | GSE31357, GSE115614, GSE330865, GSE297195, GSE283401, GSE163943 |
| `unusable` | **6** | GSE59645, GSE64978, GSE68207, GSE80174_Hippocampus, GSE304956, GSE276942 |

Note that `usable` (the boolean column) is `True` for 8 units while `cat` is
`usable` for only 2. That is intentional and not a contradiction: the boolean
means "the audit does not forbid analysis", whereas `cat` separates
unconditional use from use **under a documented restriction**. Only the latter
is counted as clean.

> An earlier revision of `audit_summary.json` carried a stale `2 / 5 / 7` split
> from before a verdict was reclassified. It was found during repository
> assembly and recomputed. `checks/verify_numbers.py` now asserts that the
> summary, the table and the manuscript agree, so the two cannot drift apart
> again silently.

---

## Two BH values for the same finding (not an inconsistency)

The interferon result is reported with two different Benjamini–Hochberg values:

| Source | Fisher input | `p_BH` | Used for |
|---|---|---|---|
| `derived/cross_dataset_meta.tsv` | one-sided | `0.0046` | intermediate analysis |
| manuscript; `final_numbers.json` | two-sided | **`0.0079`** | **what is quoted** |

Both are computed correctly; they differ because BH depends on the input
p-values across the whole family, and the two-sided Fisher input is the larger
of the two. The manuscript deliberately quotes the **more conservative**
two-sided value. If you are diffing the text against `cross_dataset_meta.tsv`
and see `0.0046`, this is why.

The underlying meta-analysis statistics agree across both inputs:

- Stouffer `Z = 3.328`, one-sided `p = 4.38 × 10⁻⁴`
- per-dataset Δz: `−0.225` (p = 0.037, GSE330865), `−0.611` (p = 0.0076, GSE297195)

---

## The one finding that survived

`IFN_signaling`, reported in `derived/cross_dataset_meta.tsv` and
`derived/verdicts/*`:

- down-regulated in **both** datasets (same sign)
- meta `Z = 3.33`, `BH = 0.0079`
- **but only the sensing / transduction layer moves.** `ISG_core` (p = 0.286),
  `MHC_I_AP` (p = 0.74) and `ISG_effector` (opposite direction) do not follow.
  The defensible statement is *"interferon signalling components are
  down-regulated"*, **not** *"the interferon response is suppressed"*.
  Counter-examples `Jak1`, `Irf3`, `Irf5` behave against the aggregate and are
  reported rather than smoothed over.

The set has **13 defined members; 12 are analysed.** `Ifnb1` (interferon-β
itself) is absent from the expression matrices in both datasets — it is an
inducible ligand with near-zero baseline expression, so it is filtered out
upstream. The two datasets use exactly the same 12 genes, so the comparison is
like-for-like. See `data/genesets/interferon_mhc_genesets_mouse.tsv`.

---

## Negative results are kept negative

- The original m6A × PCD hypothesis did **not** replicate: direction agreement
  **11/21 = 52.4 %**, binomial `p = 1.000`. This is the primary result.
- The two domains do not differ in unusable fraction (2/6 vs 4/8, Fisher
  `p = 0.63`) — with power this low that is *"no difference detected"*, not
  *"no difference exists"*.
- `mod_tan` has the smallest BH (4.9 × 10⁻⁹) but is **excluded as untrustworthy**:
  one side p = 0.154, single-dataset-driven, no external identity.
- The WGCNA reconstruction of GSE297195 is **permanently discarded** (β = 18,
  96.3 % grey, largest module 79 genes at n = 23).
