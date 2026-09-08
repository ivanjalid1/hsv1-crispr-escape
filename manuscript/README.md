# `manuscript/`

Source and build for the preprint *"What a genome corpus can and cannot certify about
CRISPR antiviral escape: resolution floors, joint coverage, and a worked audit of an
HSV-1 guide pair"*.

## What each file is

| file | what it is |
|---|---|
| `manuscript.md` | The manuscript itself, and the only place its prose is edited. Sections 1-10 plus the abstract and the figure captions in Section 6. Every number in it is produced by the pipeline; nothing here is hand-entered. |
| `references.md` | The reference list, plus the citation audit that produced it: notes on how each entry was verified against a primary record, and the closed checklist of the fourteen `[CITATION NEEDED]` placeholders the draft once carried. Only the **Verified references** section goes into the typeset PDF; the audit apparatus is repository-internal. |
| `build_pdf.py` | Builds the submission PDF from the two files above and `../figures/*.png`. |
| `verify_pdf.py` | Checks the built PDF against those sources. Run it after every build. |
| `hsv1-crispr-escape-preprint.pdf` | The build product: the submission-ready A4 PDF, tracked in the repository as a release artefact. |

## Regenerating the PDF

```bash
# once
.venv/Scripts/python.exe -m pip install pymupdf websockets

# build, then check
.venv/Scripts/python.exe manuscript/build_pdf.py
.venv/Scripts/python.exe manuscript/verify_pdf.py
```

`build_pdf.py --keep-html` also leaves the intermediate
`hsv1-crispr-escape-preprint.html` in place, which is the fastest way to inspect a
layout problem in a browser. That file is a build artefact and is not tracked.

## How the build works, and why

Markdown to styled HTML to headless Chrome. There is no LaTeX distribution, no pandoc
and no wkhtmltopdf on the build machine, and installing one was out of scope; Chrome is
present, its print pipeline embeds font subsets properly, and it is driven here over the
DevTools Protocol (`Page.printToPDF`) rather than the `--print-to-pdf` command-line flag
because only the protocol call accepts a custom footer template — that is where the page
numbers come from.

Specifics worth knowing before editing the CSS:

- **Page geometry.** A4, 25 mm text margins left and right, 22/20 mm top and bottom.
  The body is a fixed 170 mm wide with a 10 mm left gutter, so the text block is exactly
  160 mm and screen layout is identical to print layout at scale 1. That identity is
  what makes the line numbering correct; do not make the body width relative.
- **Line numbers.** Measured in the browser, then each *visual* line is promoted to its
  own block element with its number absolutely positioned in the gutter. The obvious
  implementation — one absolutely positioned number per line box inside the paragraph —
  silently mis-numbers every paragraph that crosses a page break, because a number whose
  offset lands in the leftover strip at the foot of a page stays there while its line
  moves to the next page. `verify_pdf.py` asserts the numbers run 1..N unbroken in
  reading order, which is the check that catches a regression here.
- **Figures.** Placed at the paragraph that cites them (Figure 1 in Section 2.1,
  Figure 2 in 2.2, Figures 3 and 5 in 2.3, Figure 4 in 2.5), each with its full caption
  from Section 6, and `break-inside: avoid` so neither figure nor caption is ever split.
  The 300 dpi PNGs are used and embedded as base64, not the PDF versions, because Chrome
  will not render a PDF inside an `<img>`. Section 6 keeps its heading and numbering and
  states where the figures are, so section numbers in the PDF match `manuscript.md`.
- **Underscores are not emphasis.** The Markdown converter deliberately ignores `_`,
  because the manuscript uses it for subscripts (`q_i`, `tol_invariant`) and guide
  identifiers (`RL2_5335+`).

## What `verify_pdf.py` checks

Every paragraph, list item, heading, code line, table cell and reference in the sources
is present in the PDF text verbatim; the character multiset of the PDF covers that of
the sources, so nothing can have been silently dropped; no Markdown syntax,
`[CITATION NEEDED]` placeholder or build marker leaked through; all fonts are embedded
subsets; all five figures are present; nothing is drawn outside the text block on any
page; and the line numbers are contiguous.

It exits non-zero on any failure, so it can be wired into a pre-release check.
