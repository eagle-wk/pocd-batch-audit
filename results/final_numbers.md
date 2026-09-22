# Final recomputed numbers for manuscript v0.4

## 1. Batch provenance — corrected grading

| Provenance class | n | datasets |
|---|---:|---|
| inferred | 7 | GSE283401, GSE31357, GSE64978, GSE68207, GSE80174_Cortex, GSE80174_Hippocampus, GSE80174_Thalamus |
| measured | 3 | GSE115614, GSE304956, GSE59645 |
| GSM-numbering block | 2 | GSE297195, GSE330865 |
| design-derived | 1 | GSE276942 |
| no batch label | 1 | GSE163943 |

> **Correction A7b.** The manuscript claims "measured batch metadata for 6/14". In fact GSE330865 and GSE297195 have `batch = UNKNOWN` for every sample; their "batch" is the GSM numbering block (B1-B6), which is a submitter numbering artefact of exactly the kind the manuscript warns about. Only **3/14 carry truly measured batch metadata** (scan date or declared processing block).

## 2. Corrected verdict table (n = 14 contrast-units, 12 GEO series)

| # | Unit | Domain | n | Provenance | Verdict | Basis |
|---:|---|---|---:|---|---|---|
| 1 | GSE330865 | POCD | 30 | GSM-numbering block | **caveat** | 6/6 groups single GSM block; separable (0.43 vs 0.17, p<0.001) but fully explained by nesting in `age` → not independent batch evidence |
| 2 | GSE297195 | POCD | 33 | GSM-numbering block | **caveat** | 6/6 groups single GSM block; separable (0.52 vs 0.21, p<0.001) but fully explained by nesting in `timepoint` → not independent batch evidence |
| 3 | GSE283401 | POCD | 96 | inferred | **caveat** | inferred batch inseparable from group (0.92 vs baseline 0.97) |
| 4 | GSE304956 | POCD | 29 | measured | **unusable** | E: cross-batch replication failed in 3/3 contrasts |
| 5 | GSE276942 | POCD | 17 | design-derived | **unusable** | batch ≡ timepoint by design; PC1 36.1% captured by batch |
| 6 | GSE163943 | POCD | 8 | no batch label | **caveat** | no batch label; n=8 ceiling; not an independent validation cohort |
| 7 | GSE59645 | TBI | 16 | measured | **unusable** | B: 6/6 groups single scan date; separable (0.94 vs 0.50, p=0.005) and not nested → genuine red line; PC1 R2_batch=0.961, R2_group=0.996 |
| 8 | GSE115614 | TBI | 14 | measured | **caveat** | 5/5 confounded but not separable (0.86 vs 0.86, p=0.030) → downgraded |
| 9 | GSE31357 | TBI | 32 | inferred | **caveat** | 7/8 confounded but not separable (0.72 vs 0.91, p=0.050) → downgraded |
| 10 | GSE68207 | TBI | 8 | inferred | **unusable** | C: 4703/5000 high-variance genes batch-dominated; + E failed |
| 11 | GSE64978 | TBI | 10 | inferred | **unusable** | E: cross-batch replication failed |
| 12 | GSE80174_Cortex | TBI | 10 | inferred | **usable** | batch separable from group; E passed (r=0.77, sign 86%) |
| 13 | GSE80174_Thalamus | TBI | 10 | inferred | **usable** | batch separable from group; **E not computable** (controls in a single batch level) |
| 14 | GSE80174_Hippocampus | TBI | 9 | inferred | **unusable** | B: injured group 100% single batch; C: 4820/5000 |


### Distortion rate (corrected)

| Tier | k/n | % | 95% CI (Wilson) |
|---|---:|---:|---|
| usable | 2/14 | 14.3% | [4.0%, 39.9%] |
| caveat | 6/14 | 42.9% | [21.4%, 67.4%] |
| unusable | 6/14 | 42.9% | [21.4%, 67.4%] |

| Domain | usable | caveat | unusable | n | unusable % |
|---|---:|---:|---:|---:|---:|
| POCD | 0 | 4 | 2 | 6 | 33.3% |
| TBI | 2 | 2 | 4 | 8 | 50.0% |

- Fisher exact (domain × unusable): **OR = 0.500, p = 0.627**

> **Correction B3.** With GSE297195 reclassified, the two domains are no longer identical (33% vs 50%). The Fisher test has essentially no power (p = 0.63). State: "no domain difference detected", never "the rate reproduces".

> **New positive finding.** The nesting threshold (third gate) is what reclassified GSE297195: its blocks are technically separable (0.52 vs 0.21, permutation p < 0.001) but the block difference is entirely explained by `timepoint`, so it is not independent batch evidence. Same for GSE330865 (`age`). The manuscript previously demonstrated only the separability gate; it can now demonstrate all three gates acting on real data.

## 3. Simulation & tool benchmark

| Run | genes | reps | correct | 95% CI |
|---|---:|---:|---:|---|
| Simulation benchmark | 1500 | 500 | 497/500 (99.4%) | [98.3%, 99.8%] |
| Tool comparison | 2000 | 500 | 496/500 (99.2%) | [98.0%, 99.7%] |

### Scenario C (batch ≡ group), true action = stop

| Tool | retention | note |
|---|---:|---|
| ComBat | 1.74e-16 | machine precision — true signal erased |
| RUV | 0.0185 (1.9%) | signal erased |
| sva (supervised) | 1.0000 (100.0%) | **corrected nothing** — see note |
| our audit | stop (100/100) | warning issued |

> **Reframing B1 (mandatory).** sva retaining ~100% is *not* success. Run supervised with the group protected, sva cannot estimate a surrogate variable orthogonal to the group when batch ≡ group, so it degenerates to correcting nothing: 100% of the signal is retained because **100% of the batch contamination is also retained**, with no warning. ComBat and RUV by contrast actively remove the group difference as if it were batch. Three tools, three different silent failures, one correct action: stop.

- Scenario B (legitimately correctable) ComBat retention = 112.6% — noise re-expansion, as expected. The tools are sound when their premise holds; that is why a gate is needed, not a replacement.

## 4. Cross-dataset meta-analysis (§3.5) — recomputed

- Entities tested: 21; directionally concordant **11/21 = 52.4%**, binomial two-sided **p = 1.000**

- **IFN_signaling** (n_genes = 12): GSE330865 Δz = -0.225 (p = 0.037); GSE297195 Δz = -0.611 (p = 0.0076)
- Stouffer (concordance-signed): Z = +3.328, **one-sided p = 4.38e-04**, two-sided p = 8.76e-04
- Fisher, two-sided p as input (standard): χ²₄ = 16.35, **p = 2.58e-03**
- Fisher, one-sided p as input (as in the original script): χ²₄ = 19.13, p = 7.43e-04
- BH across 21 entities (two-sided p input): **q = 0.0079**

> **Correction A2/A3.** Report Stouffer and Fisher side by side, and label the one-/two-sided convention. The conclusion is robust: q < 0.05 under every variant (Stouffer one-sided 4.4e-4; Fisher two-sided-input 2.6e-3; BH 0.0079).

> **Sign convention (must state).** Z is positive while both Δz are negative because Z is signed by *directional concordance*, not by up/down. Without this sentence a reader reads "+Z" as up-regulation.

- **mod_tan** is the *most* significant entity (BH = 4.9e-09), yet is absent from the manuscript. Excluded because: GSE297195 side p = 0.1542 (not significant → single-dataset driven) and no external gene-set identity. Both reasons must be stated in §3.5.

- **p-value underflow.** mod_tan reports p_865 = 0 exactly. A p of exactly 0 is impossible from a permutation test; it means "0 of N permutations exceeded the statistic" and must be reported as `< 1/N`.
