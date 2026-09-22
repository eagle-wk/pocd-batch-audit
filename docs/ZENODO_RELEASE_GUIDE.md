# Getting a DOI: step by step

Target: a citable, permanent DOI for the audit code and results, with GitHub as
the working copy and Zenodo as the archive. Roughly 30 minutes of actual work,
most of it waiting.

> **Two routes to the same kind of DOI.**
>
> *This guide* uses Zenodo's GitHub integration: connect the repository once,
> then every tagged release is archived automatically. Zenodo also mirrors the
> code to Software Heritage, so the snapshot stays readable even if GitHub
> disappears.
>
> If you would rather not deal with git at all, run
> `python build_zenodo_upload.py` (from `pocd_data/`) and it produces a folder
> with a ready-to-upload ZIP. Follow
> `pocd_data/zenodo_upload/_先读我_如何上传.md`. The two routes are independent
> — you can upload by hand now and turn the integration on later, but do **not**
> archive the same version twice or you will end up with two DOIs for one thing.

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

> Keep a file named `.zenodo.json` only if you want Zenodo to read metadata from
> it. It is read **only** for the release that triggers archiving. If you would
> rather type the metadata into the Zenodo web form, delete the file and fill the
> form by hand; nothing else depends on it.

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

## 4. Connect Zenodo to GitHub

1. Sign in at [zenodo.org](https://zenodo.org) with GitHub (or link GitHub under
   *Settings → Applications*).
2. Go to **Settings → GitHub** and find the repository in the list.
3. Flip its switch **on**.

That is the whole integration. From now on, **every GitHub release** in that
repository is archived automatically and gets its own DOI.

> Zenodo will not archive a repository that is already archived elsewhere, and it
> will not archive one it cannot see. If the repository does not appear in the
> list, it is either private or the GitHub app lacks access to it.

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

---

## 6. Collect the DOI

Give Zenodo a minute. Under **Zenodo → Upload → My uploads** a new record
appears (it may take a few minutes to show; the badge on the GitHub repo page
appears faster).

The record shows two DOIs:

- a **version DOI** — `10.5281/zenodo.XXXXXXX`, specific to v1.0.0
- a **concept DOI** — `10.5281/zenodo.YYYYYYY`, always the newest version

**Cite the concept DOI in the manuscript.** It keeps resolving after you publish
v1.1. Use the version DOI only when you mean this exact snapshot.

Back-fill it:

```bash
python ../setup_repo_identity.py --doi 10.5281/zenodo.YYYYYYY
git add -A
git commit -m "Add Zenodo DOI"
git push
```

This rewrites `README.md`, `CITATION.cff` and `.zenodo.json`. Note that editing
`.zenodo.json` after archiving does **not** retroactively change the record — on
the GitHub-integration route the file is read at release time, so it only affects
the next release. To correct a record's metadata, use the **Edit** button on the
record itself; the DOI does not change.

---

## 7. Back-fill the manuscript

The manuscript's Data and code availability statement must point at the
repository and the DOI. Update it through the generator, not by hand:

```
pocd_data/make_v05.py        # the single source of truth for the drafts
```

Re-run it so both the English and Chinese drafts are regenerated, then confirm:

```bash
MANUSCRIPT_DIR=../.. python checks/verify_figures.py
```

The DOI in the statement should be the **concept DOI**.

---

## 8. Verify it resolves

```bash
curl -sI https://doi.org/10.5281/zenodo.YYYYYYY | head -3   # expect 302
```

Also open the Zenodo record once and confirm: licence reads MIT, creators are the
real authors, and the description is not the placeholder text.

---

## Blind review

If the journal requires the code to be anonymised during review, do **not**
publish the DOI in the manuscript yet. Two options:

- **Recommended:** put the repository on a private-but-shareable service that
  issues an anonymous link (e.g. an anonymised Zenodo deposit), cite that in the
  manuscript, then swap in the real DOI at acceptance.
- Keep the repository private and state *"code will be released upon
  acceptance"* — weaker, and some journals now reject it.

Either way, prepare the repository now so that acceptance only requires flipping
a switch.

---

## Common ways this goes wrong

| Symptom | Cause |
|---|---|
| DOI metadata shows `AUTHOR_FAMILY_NAME` | Step 1 was skipped, or `.zenodo.json` was edited after tagging |
| No Zenodo record appears | The repo switch in *Settings → GitHub* is off, or the release is a draft |
| `git push` rejected for large files | Raw GEO data got staged; fix `.gitignore` and rewrite history before pushing |
| DOI resolves to an old version | The version DOI was cited instead of the concept DOI |
| Figures differ on a fresh clone | A dependency not pinned, or a path still pointing at a local directory |
| Zenodo created a new record instead of a version | `CITATION.cff` / `.zenodo.json` version strings disagree with the tag |

---

## After acceptance

1. Fill in `journal` and `doi` in the `preferred-citation` block of
   `CITATION.cff`.
2. Tag `v1.0.1` if anything in the code changed during revision, so the archived
   version matches what the paper describes.
3. Put the concept DOI in the manuscript's availability statement and in the
   journal's submission form — they are two separate places.
