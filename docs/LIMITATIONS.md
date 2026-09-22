# Limitations and known issues

Kept deliberately blunt. Everything here is either a limit we could not remove
or a place where an earlier revision was wrong and the correction is recorded so
it is not rediscovered as a surprise.

---

## 1. What the audit cannot decide

**The audit judges design admissibility, not biology.** A `usable` verdict means
"this contrast can be analysed without the batch explaining the group effect". It
says nothing about whether the contrast is interesting, well-powered, or
correctly normalised for anything else.

**Sample sizes were not a selection criterion we could afford.** Roughly half the
corpus is under 20 samples per group. Effects at that size are estimable but not
robust, which is one reason the primary result is negative.

**The two domains cannot be compared for distortion rate.** 2/6 unusable in the
POCD domain versus 4/8 in the TBI domain gives Fisher `p = 0.63`. With these
counts that is *"no difference detected"*. It is **not** evidence that the
domains are equally affected, and it is written that way in the manuscript,
which says explicitly that the comparison has near-zero power.

**Human coverage is one dataset.** GSE163943, n = 8, 4 POD versus 4 control, has
been used in at least four prior publications and is at its ceiling. There is no
second human bulk resource in this space. This is a property of the field, not of
the search strategy — a wider search will not fix it.

**Disease entity is not uniform.** The human-facing resources measure
*postoperative delirium* (POD), not POCD. They are labelled as such throughout
and are not treated as a POCD validation cohort. GSE304956 is a burn/excision
injury model; it is retained as a systemic-insult comparison and explicitly
labelled rather than silently dropped.

---

## 2. Batch provenance is graded, and two labels are weaker than they look

The `batch_source` field is one of `measured`, `GSM block`, `inferred`,
`design-derived`, `none`. Two entries deserve flagging:

- **GSE330865 and GSE297195 are `GSM block`, not `measured`.** The submitter's
  batch column reads `UNKNOWN` for every sample; the block structure was
  recovered from GSM numbering. Anything describing these as experimentally
  measured batch assignments is wrong.
- **All eight `inferred` labels come from a median split of library size in a
  unimodal distribution.** `libsize_bimodal` is `False` for all of them. An
  earlier revision described the provenance as "inferred from bimodal library
  size", which overstates the evidence; the correct wording is *"inferred by
  median split of library size (no bimodality detected)"*. This is the kind of
  overstatement the audit exists to catch, and it had to be caught in our own
  text first.

---

## 3. A criterion applied inconsistently across the corpus

GSE297195, GSE330865 and GSE283401 were judged by an **older version of the audit
script that had no permutation test**. They carry `loocv_accuracy_on_technical`
but no `loocv_perm_p`, so the technical-separability gate cannot be evaluated for
them on the same footing as the rest.

Consequences, both handled explicitly:

- **Figure 6** plots only the units with a complete dual gate. This is why the
  selection is declared by domain and field presence in
  `figures/make_fig4_fig6.py` rather than by directory location.
- The manuscript states the version difference where the three-way comparison
  is discussed, rather than presenting all units as judged by one uniform
  criterion.

Re-running those three with the current script would remove the inconsistency and
is the right thing to do before any extension of this work.

---

## 4. `GSE59645` ↔ PMID 31442236 — verification record

GEO links both **GSE59645** and **GSE115614** to PMID 31442236, whose title is
about *microRNA profiling* of an antidepressant compound. Table 2 describes
GSE59645 as a six-arm design with two scan dates, which does not obviously match
that title. Because the scan-date confound is the strongest single claim in the
paper, the link was verified by hand before submission.

**The association holds.** Three independent lines of evidence:

1. `!Series_pubmed_id = 31442236` in the GEO records for **both** series.
2. The paper's data-availability statement lists GSE59645, GSE59646 and GSE115614.
3. Overlapping authors and affiliation (GenUs BioSystems).

The title mismatch is expected rather than suspicious: the three series belong to
superseries **GSE59647**, which contains an mRNA arm (GSE59645) and an miRNA arm
(GSE59646) over the same 16 samples. The paper leads with the miRNA result;
GSE59645 is its mRNA counterpart.

**Correction to an earlier draft.** It claimed that *both* comparisons in
GSE59645 are zero-overlap between arms and batches. That is wrong. The scan-date
layer separates the treated arms (all 2012-04-01, slides C/D) from the control
side (all 2012-05-04, slides A/B), which is what makes the **treated-versus-
control contrast** unusable at both resolutions. The **sham-versus-TBI-only**
contrast shares 2012-05-04 and separates only at the slide layer, so it is *not*
date-confounded. The manuscript carries the corrected statement in both
languages.

**Reproduced arm sizes differ from the series description.** GEO's
`Series_overall_design` says three replicates per arm; the GSM records give
**3/2/3/3/3/2** (sham and E33 have two each). The GSM records were used.

---

## 5. Statistics that look stronger than they are

**`sva` retaining 100 % of the true signal in scenario C is a failure, not a
win.** In a fully collinear design `sva` is blind to the confound, issues no
warning, and effectively corrects nothing — so the signal passes through
untouched. Reading its 100 % as "the best tool" inverts the result.

**ComBat's and RUV's 0.0 % retentions are floating-point zero** (order 1e-16),
not a measured near-zero. "0.0 %" is defensible; it should be qualified as "at
machine precision".

**None of the three tools produces a design-admissibility warning.** Only the
audit returns a `stop` verdict. This is a statement about the tools' interfaces,
and it belongs in text and legend — it is **not** a measured quantity. An earlier
version of Figure 3 plotted it as a bar chart, which asserted an interface
property as if it were data; that panel was replaced by the audit's own recorded
verdict distribution.

**`mod_tan` is excluded despite having the smallest BH (4.9 × 10⁻⁹).** It is
driven by a single dataset (the other side gives p = 0.154) and has no external
identity. The exclusion rationale is stated in the manuscript rather than the
module being quietly omitted.

**The WGCNA reconstruction of GSE297195 is permanently discarded** (β = 18,
96.3 % of genes in grey, largest module 79 genes at n = 23). It cannot support
any claim.

---

## 6. Gene-set bookkeeping

`IFN_signaling` has **13 defined members but 12 are analysed**. `Ifnb1` —
interferon-β itself — is absent from the expression matrices of both datasets,
being an inducible ligand with near-zero baseline expression that is filtered
upstream. Both datasets use the identical 12 genes, so the cross-dataset
comparison is like-for-like, but the difference is documented here because
counting the definition file gives 13.

The defensible claim from this set is narrow, and the manuscript keeps it narrow:

> the **sensing and transduction** components of interferon signalling are
> down-regulated

not "the interferon response is suppressed". `ISG_core` (p = 0.286),
`MHC_I_AP` (p = 0.74) and `ISG_effector` (opposite direction) do not follow, and
`Jak1`, `Irf3` and `Irf5` run against the aggregate. Those counter-examples are
reported.

---

## 7. Where an earlier revision was wrong

Recorded so these do not get reintroduced:

| Claim in an earlier draft | Status |
|---|---|
| Module `royalblue` = interferon / MHC-I | **Wrong.** Overlap with the pre-defined interferon family is 5/121 = 4.1 %, and none of the 13 genes are in the WGCNA network at all — structurally it cannot be that module (it is `Olfr*` olfactory receptors plus epithelial/choroid plexus). The actual immune module is `darkturquoise`. |
| Hypothesis direction: m6A ↓ with interferon/MHC-I ↑ | **Wrong direction.** Interferon *signalling* is down-regulated. |
| Both domain comparisons in GSE59645 are zero-overlap | **Wrong.** Only the treated-versus-control contrast is; see §4. |
| Batch provenance for GSE330865 / GSE297195 is "measured" | **Wrong.** GSM-numbering block; submitter column is `UNKNOWN`. |
| Inferred batches come from bimodal library size | **Wrong.** Unimodal; median split. See §2. |
| The N = 599 cohort is a POCD cohort | **Wrong.** It is a POD cohort (PMID 41544911), and its raw data are not public. |
| `sva` is the best-correting tool in scenario C | **Wrong reading.** Its 100 % retention means it corrected nothing. See §5. |
