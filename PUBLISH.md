# PUBLISH.md — taking this repository public and minting its DOI

This is the author's checklist. It goes from a local repository with no remote to a
public GitHub repository, an archived Zenodo record with a citable DOI, and a
manuscript whose last unresolved citation (`manuscript/manuscript.md:376`) is filled.

Work through it in order. Steps 0–2 happen before anything is public and are the ones
that are expensive to undo.

Line numbers below were refreshed by the commit that recorded the author identity.
If a file has been edited since, find the placeholder instead of trusting the number:

```bash
grep -rn "\[AUTHOR NAME\]\|\[AUTHOR GIVEN NAME\]\|\[AUTHOR FAMILY NAME\]\|\[AFFILIATION\|\[ORCID\]\|\[GITHUB-USER\]\|\[REPO-NAME\]\|\[BIORXIV-DOI\|\[ZENODO-\|\[YYYY-MM-DD\]" \
  README.md LICENSE CITATION.cff .zenodo.json manuscript/manuscript.md
```

---

## Step 0 — DECIDE THIS FIRST: what name goes on the git history — **DONE**

**This was the one irreversible privacy decision, and it had to be made before the
first push. It has been made, executed and verified.**

### The decision

Option **(c)** below: the whole history was rewritten to carry the identity the paper
will be published under.

| field | value |
|---|---|
| Name | `Ivan Heredia Jalid` |
| Email | `ivanjalid@gmail.com` |
| Affiliation | `Independent Researcher, Córdoba, Argentina` |
| ORCID | **not registered yet** — see Step 1b. Not invented, so the placeholder stays. |

Two surnames, no hyphen. `Heredia Jalid` is one family name and is recorded as a
single field everywhere the metadata is structured, so that indexers cannot split it
and attribute half the work to someone else.

### What was done

Every commit in this repository had originally been authored **and** committed under a
display name and an **institutional email address**. All of them were rewritten:

```bash
git config user.name  "Ivan Heredia Jalid"
git config user.email "ivanjalid@gmail.com"

git filter-branch -f --env-filter '
  export GIT_AUTHOR_NAME="Ivan Heredia Jalid"
  export GIT_AUTHOR_EMAIL="ivanjalid@gmail.com"
  export GIT_COMMITTER_NAME="Ivan Heredia Jalid"
  export GIT_COMMITTER_EMAIL="ivanjalid@gmail.com"
' --tag-name-filter cat -- master
```

Both the author and the committer fields were set. Setting only `GIT_AUTHOR_*` leaves
the old address in the committer field, where GitHub still displays it.

The rewrite was verified to be **identity-only** before anything was discarded: the
commit count was unchanged, every commit message and every tree hash was identical to
the pre-rewrite fingerprint, the tree hash of `HEAD` was unchanged, author and
committer dates were preserved, and `git diff` between a `backup-pre-rewrite` ref and
the rewritten branch was empty. Only then were the backup ref and `refs/original/`
deleted and the old objects expired:

```bash
git branch -D backup-pre-rewrite
git for-each-ref --format='%(refname)' refs/original/ | while read r; do git update-ref -d "$r"; done
git reflog expire --expire=now --all
git gc --prune=now --aggressive
```

The pre-rewrite commit objects are now unreachable and `git cat-file -e` on them
fails, so the institutional address is not recoverable from the object store. All of
this happened while the repository still had **no remote**, so nothing under the old
identity was ever pushed anywhere.

### Confirm it at any time

```bash
git log --all --format='%an <%ae> | %cn <%ce>' | sort -u    # expect exactly one line
git log -p --all | grep -nEI '[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}' \
  | grep -v 'you@example.org' | grep -v 'ivanjalid@gmail.com'
```

`ivanjalid@gmail.com` is excluded above because it is now *meant* to be in the tree:
it is the correspondence address on the preprint and the `email` field in
`CITATION.cff`. It is a personal address chosen for publication, not an institutional
one disclosed by accident.

The second command will still print the corresponding-author addresses of Amrani et
al. and the affiliations inside `refs/amrani2024_full.xml` and
`refs/amrani2024_full.txt`. Those are the published article's own contact details in a
verbatim CC BY-NC-ND copy of it; they are not yours and they are meant to be there.

- [x] Step 0 decided, executed and verified — option (c), identity as above

---

## Step 1 — Fill the placeholders

Nothing here is guessable, and nothing should be invented. The author identity is now
filled in everywhere (1a). Everything still open — 1b to 1e — waits on an identifier
that **does not exist yet**: an ORCID, a repository URL, a Zenodo DOI, a bioRxiv DOI.
Each is filled at the step that creates it.

### 1a. Author identity — **DONE**

Filled from the Step 0 decision. Recorded here so the values can be checked, not
re-decided.

| file:line | field | value now in the file |
|---|---|---|
| `LICENSE:3` | copyright holder | `Copyright (c) 2026 Ivan Heredia Jalid` |
| `CITATION.cff:27` | `given-names` | `Ivan` |
| `CITATION.cff:28` | `family-names` | `Heredia Jalid` — **one** field, both surnames, no hyphen |
| `CITATION.cff:29` | `email` | `ivanjalid@gmail.com` |
| `CITATION.cff:30` | `affiliation` | `Independent Researcher, Córdoba, Argentina` |
| `CITATION.cff:74-77` | the same four fields inside `preferred-citation` | as above |
| `.zenodo.json:7` | `creators[0].name` | `Heredia Jalid, Ivan` — Zenodo's **`Family, Given`** order |
| `.zenodo.json:8` | `creators[0].affiliation` | as above |
| `manuscript/manuscript.md:3` | author line | `**Ivan Heredia Jalid**` |
| `manuscript/manuscript.md:5` | affiliation | `Independent Researcher, Córdoba, Argentina` |
| `manuscript/manuscript.md:9` | correspondence | `Correspondence: Ivan Heredia Jalid <ivanjalid@gmail.com>` — published deliberately |
| `manuscript/manuscript.md:384` | Author Contributions | name and affiliation; the ORCID there is still a placeholder |
| `README.md:385,388` | the two **How to cite** entries | `Ivan Heredia Jalid` |

Section 9 of the manuscript, Competing Interests, needed no edit: it already declares
no competing interests for "the author" without naming them.

- [x] 1a done

### 1b. ORCID — **open**. Needed before Step 4, not before Step 3.

There is no ORCID yet. **Do not invent one.** Get one free at
<https://orcid.org/register> — Zenodo and bioRxiv both use it, and it is the only
durable way to keep this work attached to you. That matters more than usual for a
two-surname name, which indexers routinely split.

| # | file:line | placeholder | what to put |
|---|---|---|---|
| 1 | `.zenodo.json:9` | `"orcid": "[ORCID]"` | **bare digits with hyphens, no URL**: `0000-0002-1825-0097`. Left live on purpose: Zenodo rejects a malformed ORCID outright, so this placeholder is a hard gate that fails the deposition rather than quietly publishing a record with no ORCID. That is the intended behaviour, not a bug. |
| 2 | `CITATION.cff:36` | `# orcid: "https://orcid.org/0000-0000-0000-0000"` | **uncomment** and put the real ORCID **URL**. It ships commented out, rather than carrying an `[ORCID]` placeholder, because the CFF schema type-checks this field: a placeholder would invalidate the whole file and GitHub would silently drop the "Cite this repository" button. |
| 3 | `CITATION.cff:79` | the same line inside `preferred-citation` | same value, same reason |
| 4 | `manuscript/manuscript.md:7` | `ORCID: [ORCID]` | the ORCID as printed on the preprint |
| 5 | `manuscript/manuscript.md:384` | `ORCID [ORCID]` in Author Contributions | same |

If you truly want no ORCID at all: leave both `CITATION.cff` lines commented and
**delete** line 9 of `.zenodo.json` rather than leaving `[ORCID]` in it, or the
deposition in Step 4 will fail.

- [ ] 1b done — ORCID registered and filled in all five places

### 1c. Repository URL — fill in Step 3, once the GitHub repository exists

| # | file:line | placeholder |
|---|---|---|
| 6 | `CITATION.cff:38` | `repository-code: "https://github.com/[GITHUB-USER]/[REPO-NAME]"` |
| 7 | `CITATION.cff:39` | `url: "https://github.com/[GITHUB-USER]/[REPO-NAME]"` |
| 8 | `.zenodo.json:40` | `"identifier": "https://github.com/[GITHUB-USER]/[REPO-NAME]"` |
| 9 | `README.md:83-84` | `git clone https://github.com/[GITHUB-USER]/[REPO-NAME].git` and `cd [REPO-NAME]` |

- [ ] 1c done (Step 3)

### 1d. Release date and Zenodo DOI — fill in Steps 4 and 5

| # | file:line | placeholder |
|---|---|---|
| 10 | `CITATION.cff:44` | `# date-released: 2026-01-01` — **uncomment** and set to the release date, `YYYY-MM-DD`. Commented out for the same schema reason as `orcid`. |
| 11 | `CITATION.cff:49` | `doi: "10.5281/zenodo.[ZENODO-CONCEPT-RECORD-ID]"` |
| 12 | `README.md:389` | `doi:` `10.5281/zenodo.[ZENODO-CONCEPT-RECORD-ID]` |
| 13 | `manuscript/manuscript.md:376` | the whole `[CITATION NEEDED: …]` bracket |

- [ ] 1d done (Steps 4-5)

### 1e. bioRxiv DOI — fill in Step 7, once the preprint is posted

| # | file:line | placeholder |
|---|---|---|
| 14 | `CITATION.cff:82` | `doi: "10.1101/[BIORXIV-DOI-SUFFIX]"` under `preferred-citation` |
| 15 | `CITATION.cff:83` | `url: "https://doi.org/10.1101/[BIORXIV-DOI-SUFFIX]"` |
| 16 | `.zenodo.json:34` | `"identifier": "10.1101/[BIORXIV-DOI-SUFFIX]"` |
| 17 | `README.md:387` | `doi:` `[BIORXIV-DOI]` under **How to cite** |

- [ ] 1e done (Step 7)

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
4. Confirm `.zenodo.json` is committed and that its ORCID (Step 1b) and repository
   URL (Step 1c) placeholders are filled.
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

Fill the Step 1d placeholders (10–13), then commit and push. This is the step that
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

- [ ] Step 1d placeholders (10–13) filled
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

Fill the Step 1e placeholders (14–17) with the bioRxiv DOI, commit, push, and cut
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

- [ ] Step 1e placeholders (14–17) filled
- [ ] `v1.1.0` released
- [ ] Zenodo record shows the preprint under Related identifiers

---

## Quick reference: every placeholder, one grep

```bash
grep -rn "\[AUTHOR NAME\]\|\[AUTHOR GIVEN NAME\]\|\[AUTHOR FAMILY NAME\]\|\[AFFILIATION\|\[ORCID\]\|\[GITHUB-USER\]\|\[REPO-NAME\]\|\[BIORXIV-DOI\|\[ZENODO-\|\[YYYY-MM-DD\]\|CITATION NEEDED" \
  README.md LICENSE CITATION.cff .zenodo.json PUBLISH.md manuscript/
```

The author-identity placeholders — `[AUTHOR NAME]`, `[AUTHOR GIVEN NAME]`,
`[AUTHOR FAMILY NAME]`, `[AFFILIATION …]` — are gone, and must never come back; they
are kept in the pattern above purely as a regression check. What that grep still
returns, correctly, is `[ORCID]` (Step 1b), `[GITHUB-USER]` / `[REPO-NAME]`
(Step 1c), `[ZENODO-CONCEPT-RECORD-ID]` and the `CITATION NEEDED` bracket (Step 1d),
and `[BIORXIV-DOI…]` (Step 1e) — plus `PUBLISH.md` itself and the two
`references.md` lines that *describe* the placeholder convention. When it returns
nothing but those last two, the repository is fully filled in.

**Never invent a DOI or a URL to make the grep quiet.** A visible placeholder is
correct; a plausible-looking fabricated identifier is a catastrophic failure, for
exactly the reason `manuscript/references.md:8-11` gives about fabricated references.
