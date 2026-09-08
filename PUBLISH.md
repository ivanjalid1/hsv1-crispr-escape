# PUBLISH.md — taking this repository public and minting its DOI

This is the author's checklist. It goes from a local repository with no remote to a
public GitHub repository, an archived Zenodo record with a citable DOI, and a
manuscript whose last unresolved citation (`manuscript/manuscript.md:376`) is filled.

Work through it in order. Steps 0–2 happen before anything is public and are the ones
that are expensive to undo.

Line numbers below are as of the commit that added this file. If a file has been
edited since, find the placeholder instead of trusting the number:

```bash
grep -rn "\[AUTHOR NAME\]\|\[AUTHOR GIVEN NAME\]\|\[AUTHOR FAMILY NAME\]\|\[AFFILIATION\|\[ORCID\]\|\[GITHUB-USER\]\|\[REPO-NAME\]\|\[BIORXIV-DOI\|\[ZENODO-\|\[YYYY-MM-DD\]" \
  README.md LICENSE CITATION.cff .zenodo.json manuscript/manuscript.md
```

---

## Step 0 — DECIDE THIS FIRST: what name goes on the git history

**This is the one irreversible privacy decision, and it must be made before the first
push.**

Every one of the nine commits in this repository was authored and committed under a
single identity — your display name and your **institutional email address**. See it
for yourself before deciding:

```bash
git log --all --format='%an <%ae>' | sort -u
```

That comes from your local `git config --global user.name` / `user.email`. It is baked
into the commit objects, it is **not** covered by any `.gitignore`, and once the
repository is public it is visible to anyone on every commit page and through the API.
GitHub does not let you edit it after the fact without rewriting history and
force-pushing. Search engines and scrapers index it.

The manuscript deliberately carries `[AUTHOR NAME]` and `[AFFILIATION]` placeholders,
so the repository would currently disclose an identity the manuscript does not.
Decide which of these you want:

**(a) Keep it.** Nothing to do. Your name and institutional address appear on every
commit. This is completely normal for academic code and is what most researchers do.

**(b) Replace it with a GitHub no-reply address**, keeping your display name but not
your institutional email. GitHub issues you one at
`https://github.com/settings/emails` in the form `ID+username@users.noreply.github.com`.
Set it for future commits and rewrite the existing nine:

```bash
git config user.name  "Your Name"
git config user.email "ID+username@users.noreply.github.com"

git filter-branch --env-filter '
  export GIT_AUTHOR_NAME="Your Name"
  export GIT_AUTHOR_EMAIL="ID+username@users.noreply.github.com"
  export GIT_COMMITTER_NAME="Your Name"
  export GIT_COMMITTER_EMAIL="ID+username@users.noreply.github.com"
' --tag-name-filter cat -- --branches --tags
```

**(c) Replace it with the identity you intend to publish under**, if the name on the
paper differs from the name on your git config. Same commands as (b) with different
values.

Verify whichever you chose, and confirm nothing else is left:

```bash
git log --all --format='%an <%ae> | %cn <%ce>' | sort -u    # expect exactly one line
git log -p --all | grep -nEI '[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}' \
  | grep -v 'you@example.org'
```

The second command will still print the corresponding-author addresses of Amrani et
al. and the affiliations inside `refs/amrani2024_full.xml` and
`refs/amrani2024_full.txt`. Those are the published article's own contact details in a
verbatim CC BY-NC-ND copy of it; they are not yours and they are meant to be there.

- [ ] Step 0 decided and, if (b) or (c), executed and verified

---

## Step 1 — Fill the placeholders

Nothing here is guessable, and nothing should be invented. Fill what you know now;
the DOI placeholders are filled later, in Steps 5 and 7, because those DOIs do not
exist yet.

### 1a. Fill now, before the first push

| # | file:line | placeholder | what to put |
|---|---|---|---|
| 1 | `LICENSE:3` | `Copyright (c) 2026 [AUTHOR NAME]` | your name as the copyright holder |
| 2 | `CITATION.cff:25` | `given-names: "[AUTHOR GIVEN NAME]"` | given name(s) |
| 3 | `CITATION.cff:26` | `family-names: "[AUTHOR FAMILY NAME]"` | family name |
| 4 | `CITATION.cff:27` | `affiliation: "[AFFILIATION - independent researcher]"` | affiliation, or `Independent researcher` |
| 5 | `CITATION.cff:31` | `# orcid: "https://orcid.org/0000-0000-0000-0000"` | **uncomment** and put your real ORCID URL, e.g. `https://orcid.org/0000-0002-1825-0097`. It ships commented out on purpose: the CFF schema validates the ORCID format strictly, so a `[ORCID]` placeholder there would make the whole file invalid and GitHub would silently drop the "Cite this repository" button. **No ORCID yet? Get one free at <https://orcid.org/register> — Zenodo and bioRxiv both use it, and it is the only durable way to keep this work attached to you.** If you truly want none, leave the line commented and delete line 9 of `.zenodo.json`. |
| 6 | `CITATION.cff:33-34` | `repository-code:` / `url:` | `https://github.com/<user>/<repo>` — the URL you create in Step 3 |
| 7 | `CITATION.cff:39` | `# date-released: 2026-01-01` | **uncomment** and set to the date you cut the release in Step 4, `YYYY-MM-DD`. Commented out for the same schema reason as #5. |
| 8 | `CITATION.cff:69-72` | the same four author fields again | duplicate of #2–#5, inside `preferred-citation` |
| 9 | `.zenodo.json:7` | `"name": "[AUTHOR FAMILY NAME], [AUTHOR GIVEN NAME]"` | **`Family, Given` order** — Zenodo requires it |
| 10 | `.zenodo.json:8` | `"affiliation"` | same as #4 |
| 11 | `.zenodo.json:9` | `"orcid": "[ORCID]"` | **bare digits with hyphens, no URL**: `0000-0002-1825-0097`. Zenodo rejects the URL form and rejects an invalid ORCID outright, so a leftover placeholder here will make the deposition fail — which is the intended behaviour, not a bug. |
| 12 | `.zenodo.json:40` | `"identifier": "https://github.com/[GITHUB-USER]/[REPO-NAME]"` | the repository URL from #6 |
| 13 | `README.md:83-84` | `git clone https://github.com/[GITHUB-USER]/[REPO-NAME].git` | the repository URL from #6 |
| 14 | `manuscript/manuscript.md:3` | `**[AUTHOR NAME]**` | your name, as it will appear on the preprint |
| 15 | `manuscript/manuscript.md:5` | `[AFFILIATION — independent researcher]` | affiliation |
| 16 | `manuscript/manuscript.md:7` | `ORCID: [ORCID]` | ORCID |
| 17 | `manuscript/manuscript.md:9` | `Correspondence: [AUTHOR NAME]` | name and the contact address you want printed on the preprint. **This one is published deliberately — do not use an address you would not want scraped.** |
| 18 | `manuscript/manuscript.md:384` | author-contributions line | name, affiliation, ORCID |

### 1b. Fill in Step 5, once the Zenodo DOI exists

| # | file:line | placeholder |
|---|---|---|
| 19 | `CITATION.cff:44` | `doi: "10.5281/zenodo.[ZENODO-CONCEPT-RECORD-ID]"` |
| 20 | `README.md:389` | `doi:`10.5281/zenodo.[ZENODO-CONCEPT-RECORD-ID]`` |
| 21 | `manuscript/manuscript.md:376` | the whole `[CITATION NEEDED: …]` bracket |

### 1c. Fill in Step 7, once the preprint is posted

| # | file:line | placeholder |
|---|---|---|
| 22 | `CITATION.cff:75-76` | `doi:` / `url:` under `preferred-citation` |
| 23 | `.zenodo.json:34` | `"identifier": "10.1101/[BIORXIV-DOI-SUFFIX]"` |
| 24 | `README.md:385-387` | the preprint line under **How to cite** |

- [ ] 1a done
- [ ] 1b done (Step 5)
- [ ] 1c done (Step 7)

---

## Step 2 — Pre-flight, locally, before anything is public

```bash
# 1. The full test suite: 123 tests across eight files, all offline.
for t in tests/test_*.py; do .venv/Scripts/python.exe "$t" || echo "FAILED: $t"; done

# 2. The headline numbers, from the version-controlled results. < 1 s, offline.
.venv/Scripts/python.exe verify.py

# 3. What is actually about to be published. Read this list.
git status --short
git ls-files

# 4. Nothing large or unintended is staged.
git ls-files -z | xargs -0 ls -l | sort -k5 -n -r | head -10

# 5. Credentials, one last time.
git ls-files -z | xargs -0 grep -nIiE 'api[_-]?key *=|NCBI_EMAIL *=|secret|password|token' \
  | grep -v 'you@example.org' | grep -v '^refs/amrani2024_full'
```

Expected results: 123 passed; `verify.py` reports `ALL CHECKS PASSED`; the largest
tracked file is `figures/fig4.png` at ~184 KB; and the credential grep returns only
the documentation of the two environment variables in `README.md`, `src/common.py`
and `.gitignore` — never a value.

**Do not add `.env`, `data/raw/`, `data/genome/` or the bulk `results/` tables.** They
are gitignored on purpose and the reasoning is in `.gitignore` itself.

- [ ] 123 tests pass
- [ ] `verify.py` passes
- [ ] `git ls-files` reviewed line by line
- [ ] credential grep clean

---

## Step 3 — Create the GitHub repository and push

Pick a name. Something descriptive and stable, because it becomes part of the citation
and of the Zenodo record: `hsv1-crispr-conservation`, `hsv1-guide-escape-audit`, or
similar. Avoid renaming it later.

### With the `gh` CLI

```bash
gh auth status || gh auth login

gh repo create <REPO-NAME> \
  --public \
  --description "Alignment-free conserved CRISPR-Cas9 guide discovery, multiplex escape modelling and human off-target screening for HSV-1" \
  --source . \
  --remote origin \
  --push
```

### Or by hand

Create an **empty** public repository at <https://github.com/new> — no README, no
`.gitignore`, no licence, or the first push will conflict. Then:

```bash
git remote add origin https://github.com/<user>/<REPO-NAME>.git
git branch -M main
git push -u origin main
```

Then, on the repository page:

- [ ] the "Cite this repository" button appears on the right. This confirms
      `CITATION.cff` both parsed as YAML *and* validated against the CFF schema. If it
      does not appear, the usual cause is a placeholder left in a field the schema
      type-checks — `orcid` and `date-released` are the two that bite, which is why
      both ship commented out. Paste the file into <https://citation-file-format.github.io/cff-initializer-javascript/>
      to see the actual error.
- [ ] the licence is detected as **MIT** in the sidebar
- [ ] set the About description and topics: `crispr`, `cas9`, `hsv-1`, `guide-rna`,
      `bioinformatics`, `reproducible-research`
- [ ] Settings → General → Features: turn **Issues** on. A reviewer with a question
      needs somewhere to put it.

---

## Step 4 — Connect Zenodo and mint the DOI

**Order matters. Zenodo only archives releases created *after* the repository is
switched on, so this step must come before the release.**

1. Sign in at <https://zenodo.org> **with GitHub** (Log in → GitHub), and authorise
   the `admin:repo_hook` and `read:org` scopes it asks for.
2. Go to <https://zenodo.org/account/settings/github/>. Click **Sync now** if the new
   repository is not listed.
3. Flip the toggle **ON** for `<user>/<REPO-NAME>`.
4. Confirm `.zenodo.json` is committed and its placeholders from Step 1a are filled.
   Zenodo reads this file at release time; if the ORCID or an identifier is malformed
   the deposition fails and you will have to delete the release and redo it.
5. Create the release:

```bash
git tag -a v1.0.0 -m "First public release: analysis code and archived results for the HSV-1 CRISPR guide audit"
git push origin v1.0.0
gh release create v1.0.0 \
  --title "v1.0.0 — first public release" \
  --notes "Analysis code and archived results for the preprint. Pinned corpus: 183 complete HSV-1 genomes retrieved 2026-09-07T21:56:25Z. See README.md for a reproduce-from-clone walkthrough and verify.py for a fast check that the headline numbers reproduce."
```

6. Within a few minutes the record appears at
   <https://zenodo.org/account/settings/github/> with a DOI badge. Open it and check
   the title, description, creator, ORCID and licence came through from
   `.zenodo.json`. Fields can still be edited on Zenodo afterwards; the DOI cannot.

7. **Write down two DOIs, and use the right one.** Zenodo mints both:
   - a **concept DOI**, which always resolves to the newest version — this is the one
     that goes in `CITATION.cff`, the README and the manuscript;
   - a **version DOI**, specific to `v1.0.0`.

   The record page shows the concept DOI as "Cite all versions". Take that one.

- [ ] Zenodo GitHub integration enabled **before** the release
- [ ] `v1.0.0` tagged, pushed and released
- [ ] Zenodo record exists and its metadata is correct
- [ ] concept DOI recorded: `10.5281/zenodo.__________`

---

## Step 5 — Feed the DOI back into the repository and the manuscript

Fill placeholders 19–21 from Step 1b, then commit and push. This is the step that
closes the manuscript's last open citation.

`manuscript/manuscript.md:376` currently reads:

> **Code.** [CITATION NEEDED: public repository URL and archived DOI (e.g. Zenodo) for
> the analysis code, to be inserted on deposition.] The pipeline is dependency-light …

Replace the bracket with something of this form:

> **Code.** All analysis code is available at
> `https://github.com/<user>/<REPO-NAME>` and archived at Zenodo,
> doi:`10.5281/zenodo.<CONCEPT-ID>`. The pipeline is dependency-light …

Then update the checklist entry in `manuscript/references.md` (item 14, around line
331) from **UNRESOLVED** to **RESOLVED**, and the count in the paragraph at line 226
from "Thirteen are now RESOLVED; one remains" to "All fourteen are now RESOLVED". The
sentence at line 10 saying one placeholder remains needs the same treatment.

```bash
git add -A
git commit -m "Record the public repository URL and the archived Zenodo DOI"
git push
```

Optionally cut `v1.0.1` so the archived copy on Zenodo also contains its own DOI. This
is cosmetic — the concept DOI already resolves to it — but it makes the deposited
snapshot self-describing.

- [ ] placeholders 19–21 filled
- [ ] `manuscript/references.md` checklist item 14 marked RESOLVED
- [ ] pushed

---

## Step 6 — Post the preprint

Submit `manuscript/manuscript.md` to bioRxiv with the code URL and DOI now present in
Section 7. Category: Bioinformatics, or Genomics. Declare no competing interests, as
Section 9 already states.

- [ ] preprint submitted
- [ ] preprint DOI received: `10.1101/__________`

---

## Step 7 — Close the loop

Fill placeholders 22–24 from Step 1c with the bioRxiv DOI, commit, push, and cut
`v1.1.0`. Zenodo will archive it automatically and link the preprint as a related
identifier, so the record, the repository and the paper all point at each other.

```bash
git add -A
git commit -m "Link the posted preprint from the citation metadata"
git push
git tag -a v1.1.0 -m "Preprint posted; citation metadata now complete"
git push origin v1.1.0
gh release create v1.1.0 --title "v1.1.0 — preprint linked" --notes "Citation metadata now carries the posted preprint DOI."
```

- [ ] placeholders 22–24 filled
- [ ] `v1.1.0` released
- [ ] Zenodo record shows the preprint under Related identifiers

---

## Quick reference: every placeholder, one grep

```bash
grep -rn "\[AUTHOR NAME\]\|\[AUTHOR GIVEN NAME\]\|\[AUTHOR FAMILY NAME\]\|\[AFFILIATION\|\[ORCID\]\|\[GITHUB-USER\]\|\[REPO-NAME\]\|\[BIORXIV-DOI\|\[ZENODO-\|\[YYYY-MM-DD\]\|CITATION NEEDED" \
  README.md LICENSE CITATION.cff .zenodo.json PUBLISH.md manuscript/
```

When that returns nothing outside `PUBLISH.md` itself and the two `references.md`
lines that *describe* the placeholder convention, the repository is fully filled in.

**Never invent a DOI or a URL to make the grep quiet.** A visible placeholder is
correct; a plausible-looking fabricated identifier is a catastrophic failure, for
exactly the reason `manuscript/references.md:8-11` gives about fabricated references.
