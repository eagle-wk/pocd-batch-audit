# Data sources

Every dataset behind the audit is public on GEO. Nothing proprietary, nothing
restricted — a reviewer can re-download all of it.

**Note on sample counts.** The `n` column below is the **analysed subset**, not
the size of the GEO series. Several series were built for a broader question and
we deliberately used a slice of them (a single tissue, a matched time point, one
assay arm). The GEO total is given alongside so the difference is not mistaken
for an error. Platform IDs were pulled from the NCBI `gds` database on
2026-09-22; see `fetch_geo_metadata.py` for the exact query.

---

## Analysed units (14 contrast-units across 12 series)

### TBI / neurotrauma domain (8 units)

| Unit | GEO | Platform | Species | Tissue | GEO *n* | Analysed *n* | Contrast | Verdict |
|---|---|---|---|---|---|---|---|---|
| GSE59645 | [GSE59645](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59645) | GPL14746 | *Rattus norvegicus* | hippocampus | 16 | 16 | sham vs TBI vs TBI+PMI/JM6/E33 | **unusable** |
| GSE64978 | [GSE64978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64978) | GPL18694 | *Rattus norvegicus* | — | 10 | 10 | TBI vs control | **unusable** |
| GSE68207 | [GSE68207](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68207) | GPL18694 | *Rattus norvegicus* | — | 16 | 8 | TBI vs control | **unusable** |
| GSE80174_Cortex | [GSE80174](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80174) | GPL14844 | *Rattus norvegicus* | cortex | 30 | 10 | injured vs sham | **usable** |
| GSE80174_Hippocampus | GSE80174 | GPL14844 | *Rattus norvegicus* | hippocampus | 30 | 9 | injured vs sham | **unusable** |
| GSE80174_Thalamus | GSE80174 | GPL14844 | *Rattus norvegicus* | thalamus | 30 | 10 | injured vs sham | **usable** |
| GSE31357 | [GSE31357](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31357) | GPL4478 | *Rattus norvegicus* | — | 32 | 32 | 4 h / 24 h, control vs TBI vs drug | **caveat** |
| GSE115614 | [GSE115614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115614) | GPL15084 | *Rattus norvegicus* | — | 28 | 14 | sham vs TBI vs antidepressant arms | **caveat** |

`GSE80174` is one series read as three tissue-wise units; GPL14844 is the
platform actually used in the analysed slice.

### POCD / neurodegeneration domain (6 units)

| Unit | GEO | Platform | Species | Tissue | GEO *n* | Analysed *n* | Contrast | Verdict |
|---|---|---|---|---|---|---|---|---|
| GSE330865 | [GSE330865](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330865) | GPL34290 | *Mus musculus* | hippocampus | 120 | 30 | 3 m / 17 m / 27 m, sham vs surgery | **caveat** |
| GSE297195 | [GSE297195](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297195) | GPL34290 | *Mus musculus* | hippocampus | 66 | 33 | young / aged, sham vs surgery | **caveat** |
| GSE283401 | [GSE283401](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283401) | GPL19057 | *Mus musculus* | sorted microglia | 96 | 96 | young / old × surgery × 6 h / 48 h | **caveat** |
| GSE304956 | [GSE304956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304956) | GPL24247 | *Mus musculus* | — | 71 | 29 | burn vs excision vs sham | **unusable** |
| GSE276942 | [GSE276942](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276942) | GPL24247 | *Mus musculus* | hippocampus | 17 | 17 | post-operative time course | **unusable** |
| GSE163943 | [GSE163943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163943) | GPL26963 | *Homo sapiens* | peripheral blood | 8 | 8 | 4 POD vs 4 control | **caveat** |

> **On disease entity.** The two human-facing resources measure *postoperative
> delirium* (POD), not postoperative cognitive dysfunction (POCD). They are
> labelled as such throughout and are not used as a POCD validation cohort.

---

## Candidate series examined and dropped (5)

Listed so the selection is auditable — these were retrieved and inspected, then
excluded for a stated reason rather than silently.

| GEO | Platform | Species | *n* | Reason dropped |
|---|---|---|---|---|
| [GSE181804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181804) | GPL16791 | *Homo sapiens* | 24 | no group labels in GEO metadata |
| [GSE184942](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184942) | GPL11154 | *Homo sapiens* | 10 | no group labels in GEO metadata |
| [GSE2392](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2392) | GPL85 / GPL81 | rat, mouse | 61 | RAW archive is an HTML error page |
| [GSE45997](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45997) | GPL1355 | *Rattus norvegicus* | 9 | Affymetrix CEL; offline annotation unavailable |
| [GSE86579](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86579) | GPL11534 | *Rattus norvegicus* | 11 | Affymetrix CEL; offline annotation unavailable |

---

## Primary publications for each series

| GEO | PMID | Note |
|---|---|---|
| GSE163943 | 34093168 | Song Y et al., *Front Aging Neurosci* 2021 |
| GSE276942 | 39955470 | — |
| GSE283401 | 39696348 | — |
| GSE297195 | 40540720 | — |
| GSE330865 | 42126809 | — |
| GSE31357 | 23326402 | — |
| GSE64978, GSE68207 | 28174132 | same study |
| GSE80174 | 27530814 | — |
| GSE59645, GSE115614 | 31442236 | shared study; GSE59645 = mRNA, GSE59646 = miRNA |
| GSE304956 | — | not yet published |

`GSE59645` and `GSE115614` both resolve to PMID 31442236. GEO records that link
independently (`!Series_pubmed_id` in both series), the paper's data-availability
section lists them, and the submitting lab matches. The miRNA-sounding title of
that paper is expected: the superseries **GSE59647** contains both an mRNA arm
(GSE59645) and an miRNA arm (GSE59646). This was verified by hand before
submission — see `docs/LIMITATIONS.md`.

---

## How the raw data is fetched

Raw GEO archives are **not** shipped in this repository (~527 MB). Re-download:

```bash
python audit/geo_download_and_audit.py --candidates <candidate-list.tsv> \
                                       --outroot <workdir>
```

or fetch a single series by accession and audit it offline:

```bash
python audit/audit_one.py --input <series_matrix.txt.gz> --outroot <workdir>
```

`audit/parse_agilent_raw.py` additionally reads scan metadata (`Scan_Date`,
`Scan_OriginalGUID`) out of Agilent RAW files, which GEO's sample annotation
layer does **not** expose. That field is what identifies the GSE59645
confound — see the manuscript.
