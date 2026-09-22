# Getting a DOI: step by step

Target: a citable, permanent DOI for the audit code and results, with GitHub as
the working copy and Zenodo as the archive. Roughly 20 minutes of actual work,
most of it waiting.

> **⚠️ This project already has a reserved DOI.**
>
> **`10.5281/zenodo.22899702`** was reserved on Zenodo before packaging, and it is
> already written into the manuscript's availability statement, `README.md`,
> `CITATION.cff` and `.zenodo.json`.
>
> **So the route is fixed: manual upload.** Send
> `pocd_data/zenodo_upload/pocd-batch-audit-1.0.0.zip` (+ `README.md`) to that
> reserved record and publish it. `build_zenodo_upload.py` produces exactly those
> files, and `pocd_data/zenodo_upload/_先读我_如何上传.md` walks through the upload.
>
> **Do NOT switch on Zenodo's GitHub integration.**
>
> Zenodo's integration archives **each new GitHub release as a brand-new record
> with a brand-new DOI**. There is no option to attach a release to a DOI that was
> already reserved. Turning it on would give you:
>
> ```
> 10.5281/zenodo.22899702   <- reserved, stays empty forever
> 10.5281/zenodo.<other>    <- the integration's record, the one with the files
> ```
>
> Two DOIs for one paper -- and the manuscript cites the empty one. Untangling
> that afterwards means editing the manuscript, the cover letter and the record
> relations, and a published DOI cannot be recycled.
>
> The rest of this file documents the manual route. Where it describes the git
> side (steps 1-3, 5), that is about getting the repository published so reviewers
> can browse it -- **not** about how the archive gets its DOI.

If you have *not* reserved a DOI yet, the integration is a reasonable choice: see
Zenodo's own documentation. But do not mix the two routes on one version -- that
is how projects end up with two DOIs for the same thing.

**What is and is not reversible** (Zenodo policy, checked 2026-09):

- **Metadata** — title, creators, description, keywords — can be edited at any
  time, before *or after* publication. Editing and re-publishing does not change
  the DOI.
- **Files** — after publication, minor corrections are allowed for **30 days**
  (and the corrective draft must itself be published within 45 days of the
  original). Beyond that, or for anything substantial, you need a **new version**
  — which mints a second DOI linked to the first — or a support request.
- **Deleting the record** — by the owner only, within 30 days of publication;
  afterwards only in justified cases.

So a wrong author list is recoverable; a wrong ZIP is not (cheaply). Treat the
package — not the metadata — as the thing that has to be right the first time.

---

## 0. Pre-flight

```bash
cd pocd_data/release_repo

# No placeholders may survive into the archive -- .zenodo.json especially.
python ../setup_repo_identity.py --dry-run

# Figures must satisfy the printing contract
python checks/verify_figures.py

# Text and data must agree
python checks/verify_numbers.py

# And the figures must actually regenerate from the shipped intermediates
cd figures && python make_fig1_fig2_cases.py && python make_fig3.py \
            && python make_fig4_fig6.py && cd ..
```

All four must pass. If the last one produces figures that differ byte-wise from
`figures/output/`, stop and find out why before archiving — a DOI pointing at
figures that cannot be reproduced is worse than no DOI.

> `.zenodo.json` is shipped as the machine-readable statement of the intended
> metadata, and `build_zenodo_upload.py` renders the paste-ready form checklist
> from it. But on the manual route **Zenodo never reads the file itself** -- the
> web form is authoritative. Keep it for reference and for the packager; nothing
> breaks if the form and the file disagree, though they should not.

---

## 1. Fill in your identity

```bash
python ../setup_repo_identity.py
```

It asks only for the fields the repository is actually still missing, then
rewrites every affected file and tells you which ones it touched. If it reports
that no placeholders are left, the identity fields are already final — skip to
step 2. The DOI is deliberately left for step 6.

---

## 2. Make it a git repository

From inside `release_repo/`:

```bash
git init -b main
git add .
git status          # check nothing large or unintended is staged
git commit -m "Batch-confound audit: pipeline, derived results and figures"
```

**Check the staged size.** Raw GEO data is intentionally excluded; the whole
repository should be around 2.5 MB. If `git status` shows hundreds of MB, the
`.gitignore` is not doing its job — find out why before committing:

```bash
git count-objects -vH     # after committing
du -sh .git
```

The large files are a common git mistake and are painful to undo after a push.

---

## 3. Create the GitHub repository

Either through the web UI, or:

```bash
gh repo create pocd-batch-audit --public --source=. --remote=origin --push
```

Without `gh`, create an empty repository on github.com (**do not** initialise it
with a README or licence — you already have both) and then:

```bash
git remote add origin https://github.com/eagle-wk/pocd-batch-audit.git
git push -u origin main
```

Public is the point: a reviewer must be able to reach it without an account. If
the manuscript is still under double-blind review, see "Blind review" below.

---

## 4. Publish the reserved record — manual upload, **not** the integration

Open the reserved record on Zenodo (**Upload → My uploads**), or create it if the
reservation is not there yet:

1. **Upload → New upload.** If the DOI is already reserved, the form shows
   `DOI: 10.5281/zenodo.22899702` in a *reserve* state. Leave it alone;
   publishing activates it.
2. Drag in **two** files: `pocd-batch-audit-1.0.0.zip` and `README.md`.
3. Fill the metadata from `pocd_data/zenodo_upload/Zenodo表单填写清单.md` — every
   field has a paste-ready value. **Resource type is `Software`**, not Dataset.
   Related works: `Is supplement to` / `Software` / `URL` /
   `https://github.com/eagle-wk/pocd-batch-audit`.
4. **Publish.** The DOI goes live.

> **Do not** also upload the ZIP's contents as loose files — Zenodo expands the ZIP
> into a browsable directory tree on the record page, so uploading both duplicates
> the whole tree.

> **Do not** connect this repository under *Settings → GitHub*. See the warning at
> the top of this file: the integration mints a separate DOI and orphans the
> reserved one. That is the single most expensive mistake available here.

> The checklist in the archive was rendered from `.zenodo.json`, but on manual
> upload Zenodo **does not read** that file — the web form is authoritative.

---

## 5. Tag a release

The version string should match `CITATION.cff` (`version: 1.0.0`) and
`.zenodo.json`. Mismatched version strings are the most common metadata wart.

```bash
git tag -a v1.0.0 -m "v1.0.0 -- manuscript under review"
git push origin v1.0.0
```

Then on GitHub: **Releases → Draft a new release**, pick tag `v1.0.0`, title it
`v1.0.0`, and publish. (Alternatively `gh release create v1.0.0 --notes "..."`.)

> Tagging here is purely provenance: it records which commit the archived snapshot
> corresponds to. It is safe **only because the Zenodo integration stays off** --
> if it were ever switched on, every release would mint another DOI.


---

## 6. Collect the DOI

You already know it — it was reserved before packaging:

```
concept DOI   10.5281/zenodo.22899702      <- cite this; always the newest version
version DOI   10.5281/zenodo.<n>           <- appears alongside, pins v1.0.0 exactly
```

**Cite the concept DOI in the manuscript.** It keeps resolving after you publish
v1.1. Use the version DOI only when you mean this exact snapshot. Both are listed
in the record's **Versions** section — picking the wrong one is the most common
mistake in the whole process.

Because the DOI was reserved *before* packaging, the ZIP already carries it in
`README.md`, `CITATION.cff` and `.zenodo.json`, and the manuscript was regenerated
with it. There is nothing to back-fill.

> If the reserved DOI ever has to change, rewrite it in **one** place and re-run
> the packagers:
>
> ```bash
> python ../setup_repo_identity.py --doi 10.5281/zenodo.<new>
> python ../make_v05.py              # regenerates both manuscript drafts
> python ../build_zenodo_upload.py   # rebuilds the ZIP
> ```
>
> Editing `.zenodo.json` alone does nothing to an already-published record: on
> manual upload Zenodo never read that file. To correct **metadata** on a live
> record, use the **Edit** button on the record — the DOI does not change.

---

## 7. Back-fill the manuscript

The availability statement already names both the repository and the DOI
(`10.5281/zenodo.22899702`), because the DOI was reserved *before* packaging. If
either changes, update it through the generator, not by hand:

```
pocd_data/make_v05.py        # the single source of truth for the drafts
```

Re-run it so both the English and Chinese drafts are regenerated, then confirm:

```bash
MANUSCRIPT_DIR=../.. python checks/verify_figures.py
```

The DOI in the statement must be the **concept DOI**, and the repository URL must
match the repository exactly as it exists on GitHub. A renamed repository, or one
that was never pushed, leaves a 404 inside a file that freezes 30 days after
publication.

> The manuscript is generated, so the availability text lives in `make_v05.py`.
> Editing the HTML directly works until the next re-run, then silently reverts.

---

## 8. Verify it resolves

```bash
curl -sI https://doi.org/10.5281/zenodo.22899702 | head -3   # expect 302
```

Also open the Zenodo record once and confirm:

- licence reads **MIT**; creators are the **six** real authors in manuscript order,
  with Song Keqin marked as contact; the description is the real one, not a template;
- the file list shows the ZIP with its directory tree expanded;
- Related works carries `Is supplement to` / `Software` pointing at
  `https://github.com/eagle-wk/pocd-batch-audit`.

Then confirm that repository link actually opens. The packager probes it for you:

```bash
cd ../.. && python build_zenodo_upload.py --check-links-only    # want: ✔ 可达
```

> If the probe reports `? 网络不可达，未能判定`, that is **not** a failure: it means
> the machine could not reach GitHub. This is common behind a proxy the browser
> uses but the shell does not -- `curl` returns `000` (a timeout code) rather than
> a 404. Open the URL by hand in that case.

---

## Blind review

If the journal requires the code to be anonymised during review, do **not**
publish the DOI in the manuscript yet. Two options:

- **Recommended:** put the repository on a private-but-shareable service that
  issues an anonymous link (e.g. an anonymised Zenodo deposit), cite that in the
  manuscript, then swap in the real DOI at acceptance.
- Keep the repository private and state *"code will be released upon
  acceptance"* — weaker, and some journals now reject it.

Either way, prepare the repository now so that acceptance only requires publishing
the record that is already reserved.

---

## Common ways this goes wrong

| Symptom | Cause |
|---|---|
| **The record carries a DOI you never reserved** | The Zenodo GitHub integration was switched on. It creates its own record and leaves the reserved DOI permanently empty — the most expensive mistake available here |
| DOI metadata shows `AUTHOR_FAMILY_NAME` | The identity step was skipped, or the packet was built before it was filled |
| A code link in the record or the manuscript 404s | The repository was never pushed, or was renamed / made private / deleted after the archive froze |
| `git push` rejected for large files | Raw GEO data got staged; fix `.gitignore` and rewrite history before pushing |
| DOI resolves to an old version | The version DOI was cited instead of the concept DOI |
| Figures differ from a fresh checkout | A dependency not pinned, or a path still pointing at a local directory |
| Only part of the tree appears on the record | The ZIP was built from the wrong root; the archive must contain exactly one top-level directory |
| The ZIP differs byte-wise from `MANIFEST_sha256.txt` | Line endings drifted — a Windows `open(p, "w")` writes CRLF. Rebuild; do not hand-upload |
| The manuscript points at a dead repo | The repository was renamed. Four files plus `make_v05.py` must change together, then re-package |

---

## After acceptance

1. Fill in `journal` and the article `doi` in the `preferred-citation` block of
   `CITATION.cff`.
2. Add the article DOI to the Zenodo record's **Related works** as
   `isSupplementTo`. Metadata can be edited after publication and the DOI does not
   change, so this needs no new version. It is the only related identifier the
   record should carry beyond the code repository.
3. Bump to `v1.0.1` and re-run `build_zenodo_upload.py` if the code changed during
   revision, so the archived version matches what the paper describes. On the
   published record use **New version**; the concept DOI stays the same.
4. Put the concept DOI in the manuscript's availability statement **and** in the
   journal's submission form — they are two separate places.
