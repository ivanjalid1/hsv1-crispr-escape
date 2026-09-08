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
| `hsv1-crispr-escape-preprint.pdf` | The build product: the submission-ready US Letter PDF, tracked in the repository as a release artefact. |

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

- **Page size: US Letter, not A4.** 8.5 x 11 in = 215.9 x 279.4 mm, portrait.
  bioRxiv's submission form asks for it in as many words — *"Please use the standard
  paper size of 8.5 x 11 inches (21.6 x 27.9 cm) in portrait orientation to ensure
  successful PDF conversion"* — and a conversion failure is a rejected submission, so
  the build targets Letter even though the authors are not in a Letter country. **Do not
  change it back to A4.** The paper size lives in exactly one place, `PAGE_W_MM` /
  `PAGE_H_MM` in `build_pdf.py`; the CSS `@page` rule, the body width, the Chrome
  `printToPDF` call and `verify_pdf.py`'s margin checks are all derived from it, and
  `verify_pdf.py` asserts the MediaBox of every page really is 612 x 792 pt (Chrome is
  told `preferCSSPageSize: false`, so the CSS rule alone would not settle it).
- **Page geometry.** 25 mm text margins left and right, 22/20 mm top and bottom. The
  body is a fixed 175.9 mm wide with a 10 mm left gutter, so the text block is exactly
  165.9 mm and screen layout is identical to print layout at scale 1. That identity is
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
  `FIG_MAX_H_MM` caps the image at five eighths of the *text height*, a fraction rather
  than a constant, because that is what keeps 45-70 mm of the page free for text
  underneath a figure. Fixing it in millimetres is the trap: A4's 160 mm carried over to
  Letter's 18 mm shorter page made every figure block too tall to share a page, so each
  one started a fresh page, the document grew to 28 pages and its whitespace rose by a
  third. Page-relative, it is 27 pages at the same 9-10% whitespace A4 had.
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
subsets; all five figures are present; every page's MediaBox is the US Letter 612 x 792
pt the build asked for; nothing is drawn outside the text block on any page; and the
line numbers are contiguous.

It exits non-zero on any failure, so it can be wired into a pre-release check.
