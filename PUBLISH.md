# PUBLISH.md — taking this repository public and minting its DOI

This is the author's checklist. It goes from a local repository with no remote to a
public GitHub repository, an archived Zenodo record with a citable DOI, and a
manuscript whose last unresolved citation (`manuscript/manuscript.md:376`) is filled.
Steps 0-5 are done. **The only identifier still missing anywhere in this repository is
the bioRxiv DOI**, and only bioRxiv can create it (Steps 6-7).

Work through it in order. Steps 0–2 happen before anything is public and are the ones
that are expensive to undo.

Line numbers below were refreshed by the commit that recorded the Zenodo DOI (Steps 4
and 5). If a file has been edited since, find the placeholder instead of trusting the
number:

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
| ORCID | `0009-0003-6702-1295` |

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

Nothing here is guessable, and nothing should be invented. The author identity (1a),
the ORCID (1b) and the repository URL (1c) are now filled in everywhere. What is still
open — **1e only** — waits on an identifier that **does not exist yet**: the bioRxiv
DOI. It is filled at the step that creates it.

| still open | identifier | created by |
|---|---|---|
| ~~1d~~ | ~~Zenodo concept DOI and release date~~ — **DONE**, Steps 4-5 | Step 4 |
| 1e | bioRxiv DOI (`10.1101/…`) | Step 6 |

### 1a. Author identity — **DONE**

Filled from the Step 0 decision. Recorded here so the values can be checked, not
re-decided.

| file:line | field | value now in the file |
|---|---|---|
| `LICENSE:3` | copyright holder | `Copyright (c) 2026 Ivan Heredia Jalid` |
| `CITATION.cff:32` | `given-names` | `Ivan` |
| `CITATION.cff:33` | `family-names` | `Heredia Jalid` — **one** field, both surnames, no hyphen |
| `CITATION.cff:34` | `email` | `ivanjalid@gmail.com` |
| `CITATION.cff:35` | `affiliation` | `Independent Researcher, Córdoba, Argentina` |
| `CITATION.cff:79-82` | the same four fields inside `preferred-citation` | as above |
| `.zenodo.json:7` | `creators[0].name` | `Heredia Jalid, Ivan` — Zenodo's **`Family, Given`** order |
| `.zenodo.json:8` | `creators[0].affiliation` | as above |
| `manuscript/manuscript.md:3` | author line | `**Ivan Heredia Jalid**` |
| `manuscript/manuscript.md:5` | affiliation | `Independent Researcher, Córdoba, Argentina` |
| `manuscript/manuscript.md:9` | correspondence | `Correspondence: Ivan Heredia Jalid <ivanjalid@gmail.com>` — published deliberately |
| `manuscript/manuscript.md:384` | Author Contributions | name, affiliation and ORCID |
| `README.md:387,390` | the two **How to cite** entries | `Ivan Heredia Jalid` |

Section 9 of the manuscript, Competing Interests, needed no edit: it already declares
no competing interests for "the author" without naming them.

- [x] 1a done

### 1b. ORCID — **DONE**

`0009-0003-6702-1295`. Registered at ORCID; the ISO 7064 MOD 11-2 check digit was
verified before it was written anywhere (computed 5, declared 5).

**The two files want two different forms. This is not a style choice — each schema
rejects the other's form.**

| file:line | form required | value now in the file |
|---|---|---|
| `CITATION.cff:40` | **full URI** | `orcid: "https://orcid.org/0009-0003-6702-1295"` |
| `CITATION.cff:83` | **full URI**, inside `preferred-citation` | same |
| `.zenodo.json:9` | **bare identifier, no URL** | `"orcid": "0009-0003-6702-1295"` |
| `manuscript/manuscript.md:7` | bare, after the `ORCID:` label | `ORCID: 0009-0003-6702-1295` |
| `manuscript/manuscript.md:384` | bare, in Author Contributions | `ORCID 0009-0003-6702-1295` |

Both `CITATION.cff` lines are now **uncommented and live**. They shipped commented out
only because a placeholder would have failed the schema's type check; with a real value
there is nothing left to protect against, and leaving them commented would have thrown
away the ORCID's whole purpose. The file was re-validated against the CFF 1.2.0 schema
with the field live — `cffconvert --validate -i CITATION.cff` reports
*"Citation metadata are valid according to schema version 1.2.0"* — so GitHub's
"Cite this repository" button still works.

That check also confirms the two-surname handling survives the round trip: the
generated BibTeX is `author = {Heredia Jalid, Ivan}` and the APA form is
`Heredia Jalid I.`, with the surname intact rather than split.

`date-released` was commented out for the same reason and is now live at
`CITATION.cff:49` — see 1d below.

- [x] 1b done — ORCID registered, verified and filled in all five places

### 1c. Repository URL — **DONE**

The repository is live at <https://github.com/ivanjalid1/hsv1-crispr-escape>, public,
with `master` tracking `origin/master`. **The four places do not all want the same
form** — two schemas want the bare browse URL, the README wants the clone URL:

| # | file:line | form required | value now in the file |
|---|---|---|---|
| 6 | `CITATION.cff:42` | bare URL, no `.git` (CFF type-checks it as a URI) | `repository-code: "https://github.com/ivanjalid1/hsv1-crispr-escape"` |
| 7 | `CITATION.cff:43` | bare URL, no `.git` | `url: "https://github.com/ivanjalid1/hsv1-crispr-escape"` |
| 8 | `.zenodo.json:40` | bare URL — the entry declares `"scheme": "url"`, so a `.git` suffix would be archived as the identifier | `"identifier": "https://github.com/ivanjalid1/hsv1-crispr-escape"` |
| 9 | `README.md:85-86` | **git clone URL, with `.git`**, and the bare directory name | `git clone https://github.com/ivanjalid1/hsv1-crispr-escape.git` / `cd hsv1-crispr-escape` |

`cffconvert --validate -i CITATION.cff` still reports *"Citation metadata are valid
according to schema version 1.2.0"* with the URLs live, and the generated BibTeX now
carries `url = {https://github.com/ivanjalid1/hsv1-crispr-escape}`.

- [x] 1c done — repository created, pushed, and its URL filled in all four places

### 1d. Release date and Zenodo DOI — **DONE**

The record is live. **Zenodo minted two DOIs and they are not interchangeable:**

| DOI | what it is | where it belongs |
|---|---|---|
| **`10.5281/zenodo.22664837`** | **CONCEPT DOI** — the parent record; always resolves to the newest version | **everywhere**: `CITATION.cff`, the README badge and citation, the manuscript, `references.md` |
| `10.5281/zenodo.22664838` | version DOI — pinned to `v1.0.0` alone | nowhere in this repository. Use it only if a field ever asks specifically for a version-pinned identifier |

Verified against the API rather than the UI: `https://zenodo.org/api/records/22664838`
reports `conceptdoi: 10.5281/zenodo.22664837` and `conceptrecid: 22664837`.
**Zenodo's own GitHub settings page displays the *version* DOI in the repository row.**
That is the trap; the row is not the citation. Re-check at any time:

```bash
curl -s https://zenodo.org/api/records/22664838 | python -c "import json,sys;d=json.load(sys.stdin);print(d['conceptdoi'],d['doi'])"

# The version DOI must appear in no file but this one, which records it on purpose.
git ls-files -z | xargs -0 grep -ln "zenodo.22664838" | grep -v '^PUBLISH.md$'
```

Where the concept DOI now is, and in which form — **the forms differ and are not
interchangeable either**:

| # | file:line | form | value now in the file |
|---|---|---|---|
| 10 | `CITATION.cff:49` | bare date, uncommented | `date-released: 2026-09-08` |
| 11 | `CITATION.cff:54` | **bare DOI string**, no `https://doi.org/` prefix (CFF's `doi` field is type-checked as a DOI, not a URL) | `doi: "10.5281/zenodo.22664837"` |
| 12 | `README.md:3` | **badge**, Markdown image link, target = concept DOI | `[![DOI](https://zenodo.org/badge/1361724608.svg)](https://doi.org/10.5281/zenodo.22664837)` |
| 13 | `README.md:391-392` | bare DOI in the `doi:` entry, plus the resolvable URL | `doi:` `10.5281/zenodo.22664837` — <https://doi.org/10.5281/zenodo.22664837> |
| 14 | `manuscript/manuscript.md:376` | prose: repository URL, bare DOI **and** the resolvable `https://doi.org/…` URL; the `[CITATION NEEDED]` bracket is gone | see Step 5 |
| 15 | `manuscript/references.md:61-65` | full software-citation entry (author, title, version, `[software]`, publisher Zenodo, year, DOI, repository) | reference 10, `Heredia Jalid I.` |
| 16 | `manuscript/references.md:344-353` | checklist item 14, flipped `UNRESOLVED` → `RESOLVED` | — |

`cffconvert --validate -i CITATION.cff` reports *"Citation metadata are valid according
to schema version 1.2.0"* with both fields live, and the generated BibTeX carries
`doi = {10.5281/zenodo.22664837}`.

- [x] 1d done (Steps 4-5)

### 1e. bioRxiv DOI — fill in Step 7, once the preprint is posted

| # | file:line | placeholder |
|---|---|---|
| 17 | `CITATION.cff:86` | `doi: "10.1101/[BIORXIV-DOI-SUFFIX]"` under `preferred-citation` |
| 18 | `CITATION.cff:87` | `url: "https://doi.org/10.1101/[BIORXIV-DOI-SUFFIX]"` |
| 19 | `.zenodo.json:34` | `"identifier": "10.1101/[BIORXIV-DOI-SUFFIX]"` |
| 20 | `README.md:389` | `doi:` `[BIORXIV-DOI]` under **How to cite** |

- [ ] 1e done (Step 7)

---

## Step 1f — Why `LICENSE` is pristine and `NOTICE.md` exists — **DONE**

GitHub reported the licence as `NOASSERTION`, not `MIT`:

```bash
gh api repos/ivanjalid1/hsv1-crispr-escape --jq .license.spdx_id   # -> NOASSERTION
```

The cause was not a wrong licence — it was the shape of the file. GitHub detects
licences with [`licensee`](https://github.com/licensee/licensee), which only matches a
file against a known licence when the normalised text is a near-exact match. `LICENSE`
had a `THIRD-PARTY MATERIAL` section appended to it, and that appended text pushed it
below the similarity threshold, so licensee refused to name the licence at all.

The fix was to **move, not delete**. `LICENSE` is now the unmodified MIT text and
nothing else, and the third-party notice — a real licensing constraint on `refs/`,
covering the CC BY-NC-ND 4.0 Amrani et al. 2024 full text and the CC BY 4.0
Ramadoss et al. 2025 derived table — now lives verbatim in **[`NOTICE.md`](NOTICE.md)**.

**No pointer line was added to `LICENSE`.** Even one extra sentence risks putting the
file back below licensee's threshold, and re-breaking detection to gain a cross-
reference is a bad trade. `NOTICE.md` is instead referenced from three places that are
read by humans and by machines:

| file | how it points at `NOTICE.md` |
|---|---|
| `README.md` (Licence section, and both repository-layout listings) | prose link, with the reason the notice is not inside `LICENSE` |
| `CITATION.cff` (`message`) | one sentence in the field GitHub shows under "Cite this repository" |
| `.zenodo.json` (`notes`) | so the Zenodo record carries it too |

If `LICENSE` is ever edited again, re-check detection rather than assuming:

```bash
gh api repos/ivanjalid1/hsv1-crispr-escape --jq .license.spdx_id   # must print MIT
```

- [x] `LICENSE` restored to pristine MIT; notice moved verbatim to `NOTICE.md`
- [x] `spdx_id` confirmed `MIT` against the live API after the push

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

## Step 3 — Create the GitHub repository and push — **DONE**

The repository exists, is **public**, and the whole history is pushed:

| field | value |
|---|---|
| URL | <https://github.com/ivanjalid1/hsv1-crispr-escape> |
| owner | `ivanjalid1` |
| visibility | public |
| default branch | `master`, tracking `origin/master` |

The name is now part of the citation and will be part of the Zenodo record.
**Do not rename it.** A rename leaves a redirect on GitHub but silently invalidates the
URL already written into `CITATION.cff`, `.zenodo.json` and the README.

Re-confirm the state at any time:

```bash
gh repo view ivanjalid1/hsv1-crispr-escape --json name,visibility,url,defaultBranchRef
git status -sb            # expect: ## master...origin/master  (no ahead/behind)
```

Then, on the repository page:

- [x] the "Cite this repository" button appears on the right. This confirms
      `CITATION.cff` both parsed as YAML *and* validated against the CFF schema. If it
      ever stops appearing, the usual cause is a placeholder left in a field the schema
      type-checks — `orcid` and `date-released` are the two that bite, which is why
      `date-released` still ships commented out. Paste the file into
      <https://citation-file-format.github.io/cff-initializer-javascript/>
      to see the actual error.
- [x] the licence is detected as **MIT** in the sidebar. This did **not** work on the
      first push — it read "Other" / `NOASSERTION` — and was fixed by making `LICENSE`
      a pristine MIT text and moving the third-party notice to `NOTICE.md`. The full
      reasoning is in **Step 1f** above. Verify it from the API, not by eye:
      `gh api repos/ivanjalid1/hsv1-crispr-escape --jq .license.spdx_id` → `MIT`.
- [ ] set the About description and topics: `crispr`, `cas9`, `hsv-1`, `guide-rna`,
      `bioinformatics`, `reproducible-research`
- [ ] Settings → General → Features: turn **Issues** on. A reviewer with a question
      needs somewhere to put it.

---

## Step 4 — Connect Zenodo and mint the DOI — **DONE**

Done. The Zenodo GitHub integration was enabled first, `v1.0.0` was tagged, pushed and
released, and Zenodo archived it and minted the DOIs. **The only identifier still
missing from this project is the bioRxiv DOI, and only bioRxiv can create it.**

The record and its two DOIs:

| field | value |
|---|---|
| concept record | `22664837` — **concept DOI `10.5281/zenodo.22664837`** |
| v1.0.0 record | `22664838` — version DOI `10.5281/zenodo.22664838` |
| release date | 2026-09-08 |
| badge | `https://zenodo.org/badge/1361724608.svg` |

**Use the concept DOI.** It is the one in the paper, the citation metadata, the README
badge and `references.md`. The version DOI appears nowhere in this repository and
should stay that way unless some future field asks specifically for a version-pinned
identifier. See 1d above for how each was verified, and for the trap where Zenodo's own
GitHub settings page shows the *version* DOI in the repository row.

The steps below are kept as the record of what was done, in case a later release has to
repeat them.

**Order matters. Zenodo only archives releases created *after* the repository toggle
is switched on, so the OAuth link and the toggle must both come before the tag.** A
release cut first is not archived, and the fix is to delete the release and redo it.

1. Sign in at <https://zenodo.org> **with GitHub** (Log in → GitHub), and authorise
   the `admin:repo_hook` and `read:org` scopes it asks for.
2. Go to <https://zenodo.org/account/settings/github/>. The repository is already
   public, so it should be listed; click **Sync now** if it is not.
3. Flip the toggle **ON** for **`ivanjalid1/hsv1-crispr-escape`**.
4. `.zenodo.json` is committed and its repository-URL placeholder is filled — nothing
   to do here beyond confirming it, because Zenodo reads this file at release time and
   a malformed ORCID or identifier makes the deposition fail, after which you have to
   delete the release and redo it. To re-check before tagging:

   ```bash
   python -c "import json;json.load(open('.zenodo.json',encoding='utf-8'));print('ok')"
   git status --short          # expect empty: .zenodo.json must be committed, not just saved
   ```

   Note that at this point in the sequence `CITATION.cff` still carried its
   `10.5281/zenodo.[…]` placeholder. **That was correct and expected** — it is a
   chicken-and-egg field, filled in Step 5 from the DOI this step mints. It must never
   be given an invented value to make it look finished. It is filled now.
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

- [x] Zenodo GitHub integration enabled **before** the release
- [x] `v1.0.0` tagged, pushed and released
- [x] Zenodo record exists and its metadata is correct
- [x] concept DOI recorded: `10.5281/zenodo.22664837`
- [x] version DOI recorded and deliberately unused: `10.5281/zenodo.22664838`

---

## Step 5 — Feed the DOI back into the repository and the manuscript — **DONE**

The Step 1d placeholders (10–16) are filled and pushed. This is the step that closed the
manuscript's last open citation: **`manuscript/manuscript.md` now contains zero
`[CITATION NEEDED]` brackets.**

`manuscript/manuscript.md:376` now reads:

> **Code.** All analysis code is public at
> `https://github.com/ivanjalid1/hsv1-crispr-escape` and permanently archived at Zenodo
> (Heredia Jalid, 2026), doi:`10.5281/zenodo.22664837`, which resolves at
> `https://doi.org/10.5281/zenodo.22664837`. That is the concept DOI and always resolves
> to the most recent archived version; the release reported here is v1.0.0. The pipeline
> is dependency-light …

The author–year form matches the manuscript's citation style throughout, and the
corresponding software entry was added to `manuscript/references.md` as **reference 10,
`Heredia Jalid I.`**, inserted alphabetically between Hanley and Hsu (the entries that
followed were renumbered 10–23 → 11–24). A note on its usage was added to the
*Notes on citation usage* list, checklist item 14 was flipped from **UNRESOLVED** to
**RESOLVED**, and the two summary sentences — the one in the header at
`references.md:10` and the count at `references.md:239` — now read that all fourteen are
RESOLVED.

Optionally cut `v1.0.1` so the archived copy on Zenodo also contains its own DOI. This
is cosmetic — the concept DOI already resolves to it — but it makes the deposited
snapshot self-describing.

- [x] Step 1d placeholders (10–16) filled
- [x] `manuscript/references.md` checklist item 14 marked RESOLVED
- [x] software reference added to `manuscript/references.md` (reference 10)
- [x] Zenodo DOI badge added at `README.md:3`, targeting the **concept** DOI
- [x] pushed

---

## Step 6 — Post the preprint ← **YOU ARE HERE**

Submit `manuscript/manuscript.md` to bioRxiv with the code URL and DOI now present in
Section 7. Category: Bioinformatics, or Genomics. Declare no competing interests, as
Section 9 already states.

- [ ] preprint submitted
- [ ] preprint DOI received: `10.1101/__________`

---

## Step 7 — Close the loop

Fill the Step 1e placeholders (17–20) with the bioRxiv DOI, commit, push, and cut
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

- [ ] Step 1e placeholders (17–20) filled
- [ ] `v1.1.0` released
- [ ] Zenodo record shows the preprint under Related identifiers

---

## Quick reference: every placeholder, one grep

```bash
grep -rn "\[AUTHOR NAME\]\|\[AUTHOR GIVEN NAME\]\|\[AUTHOR FAMILY NAME\]\|\[AFFILIATION\|\[ORCID\]\|\[GITHUB-USER\]\|\[REPO-NAME\]\|\[BIORXIV-DOI\|\[ZENODO-\|\[YYYY-MM-DD\]\|CITATION NEEDED" \
  README.md LICENSE NOTICE.md CITATION.cff .zenodo.json PUBLISH.md manuscript/
```

The author-identity placeholders — `[AUTHOR NAME]`, `[AUTHOR GIVEN NAME]`,
`[AUTHOR FAMILY NAME]`, `[AFFILIATION …]` — and `[ORCID]` are gone, and so now are
`[GITHUB-USER]`, `[REPO-NAME]` and `[ZENODO-CONCEPT-RECORD-ID]`. None of them must ever
come back; they are kept in the pattern above purely as a regression check.

**Exactly one identifier remains unfilled, in four places, and it is waiting on an
external service:**

| what | where it still appears | filled by |
|---|---|---|
| bioRxiv DOI | `CITATION.cff:86`, `CITATION.cff:87`, `.zenodo.json:34`, `README.md:389` | Step 7, from the DOI Step 6 receives |

The Zenodo concept DOI and the release date are **done** and are no longer placeholders:
`CITATION.cff:49`, `CITATION.cff:54`, `README.md:3` (badge), `README.md:391-392`,
`manuscript/manuscript.md:376`, `manuscript/references.md:61-65` and
`manuscript/references.md:344-353` all carry `10.5281/zenodo.22664837`. The version DOI
`10.5281/zenodo.22664838` is intentionally absent from every file except this one, which
records it so it is not re-confused later; the check is in Step 1d.

Everything else that grep returns is `PUBLISH.md` itself and the two `references.md`
lines that *describe* the placeholder convention. When it returns nothing but those,
the repository is fully filled in.

**Never invent a DOI or a URL to make the grep quiet.** A visible placeholder is
correct; a plausible-looking fabricated identifier is a catastrophic failure, for
exactly the reason `manuscript/references.md:8-11` gives about fabricated references.
