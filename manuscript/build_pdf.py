#!/usr/bin/env python3
"""Build a submission-ready US Letter PDF of the manuscript.

Pipeline:  manuscript.md + references.md  ->  self-contained HTML  ->  headless
Chrome (DevTools Protocol Page.printToPDF)  ->  manuscript/hsv1-crispr-escape-preprint.pdf

No LaTeX, no pandoc, no wkhtmltopdf.  Requires only Google Chrome plus the
`websockets` package (CDP transport).  Figures are embedded as base64 PNGs, so
the intermediate HTML is fully self-contained.

Usage:  .venv/Scripts/python.exe manuscript/build_pdf.py [--keep-html]
"""

from __future__ import annotations

import argparse
import base64
import html as htmllib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MS = ROOT / "manuscript" / "manuscript.md"
REFS = ROOT / "manuscript" / "references.md"
FIGDIR = ROOT / "figures"
OUT_PDF = ROOT / "manuscript" / "hsv1-crispr-escape-preprint.pdf"
OUT_HTML = ROOT / "manuscript" / "hsv1-crispr-escape-preprint.html"

CHROME_CANDIDATES = [
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
]

# ---------------------------------------------------------------- page geometry
# US Letter, 8.5 x 11 in = 215.9 x 279.4 mm, portrait.  bioRxiv's submission
# form asks for this size explicitly ("Please use the standard paper size of
# 8.5 x 11 inches ... to ensure successful PDF conversion"), so the build
# targets Letter rather than A4; do not switch it back.
#
# These are the single source of truth for the page: the CSS @page rule, the
# body width, the Chrome printToPDF call and verify_pdf.py's margin checks are
# all derived from them, so changing the paper size means changing it here only.
#
# Margins below are the *paper* margins handed to Chrome; the body then carries
# a 10 mm left gutter that holds the line numbers, so the text block itself
# starts 25 mm from the left edge and ends 25 mm from the right.
PAGE_W_MM, PAGE_H_MM = 215.9, 279.4
MARGIN_TOP_MM, MARGIN_BOTTOM_MM = 22.0, 20.0
MARGIN_LEFT_MM, MARGIN_RIGHT_MM = 15.0, 25.0
GUTTER_MM = 10.0
CONTENT_W_MM = PAGE_W_MM - MARGIN_LEFT_MM - MARGIN_RIGHT_MM   # 175.9 mm
TEXT_W_MM = CONTENT_W_MM - GUTTER_MM                          # 165.9 mm
TEXT_H_MM = PAGE_H_MM - MARGIN_TOP_MM - MARGIN_BOTTOM_MM      # 237.4 mm

# A figure and its caption carry `break-inside: avoid`, so the whole block moves
# to the next page if it does not fit.  Cap the image at five eighths of the text
# height: captions run 20-45 mm, which leaves 45-70 mm of the page for text, so a
# figure can still follow the paragraph that cites it instead of starting a fresh
# page and stranding the bottom of the previous one.  (This is a *fraction of the
# page*, not a constant -- Letter's text block is 18 mm shorter than A4's, and at
# A4's old fixed 160 mm every figure began its own page and the whitespace in the
# document rose by a third.)
FIG_MAX_H_MM = TEXT_H_MM * 5 / 8                              # 148.4 mm

MM_PER_IN = 25.4

# Where each figure goes.  Two anchor kinds:
#   ("before", heading-prefix)      -- immediately before that heading
#   ("after", section, marker)      -- after the first paragraph inside
#                                      `section` that contains `marker`
# Figure 1 -> 2.1, Figure 2 -> 2.2, Figures 3 and 5 -> 2.3, Figure 4 -> 2.5.
FIG_PLACEMENT = [
    ("fig1", ("before", "### 2.2 ")),
    ("fig2", ("after", "### 2.2 ", "(Figure 2)")),
    ("fig3", ("after", "### 2.3 ", "(Figure 3)")),
    ("fig5", ("after", "### 2.3 ", "(Figure 5)")),
    ("fig4", ("after", "### 2.5 ", "Figure 4)")),
]

FIGURES_SECTION_NOTE = (
    "In this typeset version the five figures are placed inline at the sections "
    "they support, each printed with its full caption: Figure 1 in Section 2.1, "
    "Figure 2 in Section 2.2, Figures 3 and 5 in Section 2.3, and Figure 4 in "
    "Section 2.5."
)


# =============================================================== markdown -> html
def esc(text: str) -> str:
    return htmllib.escape(text, quote=False)


def inline(text: str) -> str:
    """Inline markdown: code spans, then bold, then italic.

    Underscores are deliberately NOT treated as emphasis: the manuscript uses
    them as mathematical subscripts (q_i, tol_invariant, Sigma_genomes).
    Link syntax is not parsed because the sources contain none.
    """
    holes: list[str] = []

    def stash(rendered: str) -> str:
        holes.append(rendered)
        return f"\x00{len(holes) - 1}\x00"

    # code spans first, so their contents are never touched by emphasis rules
    text = re.sub(r"`([^`]+)`", lambda m: stash("<code>" + esc(m.group(1)) + "</code>"), text)
    text = esc(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text, flags=re.S)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", text, flags=re.S)
    text = re.sub(r"\x00(\d+)\x00", lambda m: holes[int(m.group(1))], text)
    return text


def render_table(rows: list[str]) -> str:
    def cells(line: str) -> list[str]:
        line = line.strip()
        if line.startswith("|"):
            line = line[1:]
        if line.endswith("|"):
            line = line[:-1]
        return [c.strip() for c in line.split("|")]

    header = cells(rows[0])
    body = [cells(r) for r in rows[2:]]
    ncol = len(header)
    out = ['<div class="tablewrap"><table>', "<thead><tr>"]
    for h in header:
        out.append(f"<th>{inline(h)}</th>")
    out.append("</tr></thead><tbody>")
    for row in body:
        row = (row + [""] * ncol)[:ncol]
        out.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in row) + "</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


TABLE_CAPTION = re.compile(
    r'(<p><strong>Table \d+\.</strong>.*?</p>)\s*(<div class="tablewrap">.*?</div>)',
    re.S,
)


def bind_table_captions(html: str) -> str:
    """Keep a "Table N." caption on the same page as the table it introduces."""
    return TABLE_CAPTION.sub(
        lambda m: f'<div class="tableblock">{m.group(1)}{m.group(2)}</div>', html
    )


LIST_UL = re.compile(r"^([-*+])\s+(.*)$")
LIST_OL = re.compile(r"^(\d+)\.\s+(.*)$")


def md_to_html(md: str, *, drop_hr: bool = True) -> str:
    lines = md.replace("\r\n", "\n").split("\n")
    out: list[str] = []
    i, n = 0, len(lines)
    para: list[str] = []

    def flush_para() -> None:
        if para:
            out.append("<p>" + inline(" ".join(para)) + "</p>")
            para.clear()

    def peek_list_continues(j: int, ordered: bool) -> bool:
        """After a blank line, does the same list carry on?  (loose lists)"""
        k = j
        while k < n and not lines[k].strip():
            k += 1
        if k >= n:
            return False
        pat = LIST_OL if ordered else LIST_UL
        return bool(pat.match(lines[k])) and not lines[k].startswith("---")

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # fenced code
        if stripped.startswith("```"):
            flush_para()
            i += 1
            buf = []
            while i < n and not lines[i].strip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            out.append("<pre><code>" + esc("\n".join(buf)) + "</code></pre>")
            continue

        # headings
        m = re.match(r"^(#{1,6})\s+(.*)$", line)
        if m:
            flush_para()
            lvl = len(m.group(1))
            out.append(f"<h{lvl}>{inline(m.group(2).strip())}</h{lvl}>")
            i += 1
            continue

        # horizontal rule
        if re.fullmatch(r"-{3,}|\*{3,}|_{3,}", stripped):
            flush_para()
            if not drop_hr:
                out.append("<hr/>")
            i += 1
            continue

        # table
        if stripped.startswith("|") and i + 1 < n and re.match(r"^\|[\s:|-]+\|?\s*$", lines[i + 1].strip()):
            flush_para()
            rows = []
            while i < n and lines[i].strip().startswith("|"):
                rows.append(lines[i])
                i += 1
            out.append(render_table(rows))
            continue

        # lists (loose or tight, ordered or unordered)
        m_ul, m_ol = LIST_UL.match(line), LIST_OL.match(line)
        if (m_ul or m_ol) and not stripped.startswith("---"):
            flush_para()
            ordered = m_ol is not None
            tag = "ol" if ordered else "ul"
            pat = LIST_OL if ordered else LIST_UL
            start = m_ol.group(1) if ordered else None
            items: list[list[str]] = []
            while i < n:
                cur = lines[i]
                if not cur.strip():
                    if peek_list_continues(i + 1, ordered):
                        while i < n and not lines[i].strip():
                            i += 1
                        continue
                    break
                mm = pat.match(cur)
                if mm:
                    items.append([mm.group(2).strip()])
                elif cur.startswith((" ", "\t")) and items:
                    items[-1].append(cur.strip())
                else:
                    break
                i += 1
            attr = f' start="{start}"' if ordered and start and start != "1" else ""
            out.append(f"<{tag}{attr}>")
            for item in items:
                out.append("<li>" + inline(" ".join(item)) + "</li>")
            out.append(f"</{tag}>")
            continue

        if not stripped:
            flush_para()
            i += 1
            continue

        para.append(stripped)
        i += 1

    flush_para()
    return "\n".join(out)


# ======================================================================= sources
def split_sections(md: str) -> tuple[str, str, str, str]:
    """front matter, body up to section 6, section 6 (figures), rest."""
    lines = md.replace("\r\n", "\n").split("\n")
    i_abs = next(i for i, l in enumerate(lines) if l.startswith("## Abstract"))
    i_fig = next(i for i, l in enumerate(lines) if l.startswith("## 6. Figures"))
    i_after = next(i for i, l in enumerate(lines) if l.startswith("## 7. "))
    return (
        "\n".join(lines[:i_abs]),
        "\n".join(lines[i_abs:i_fig]),
        "\n".join(lines[i_fig + 1 : i_after]),
        "\n".join(lines[i_after:]),
    )


def parse_captions(fig_md: str) -> dict[str, str]:
    """Section 6 -> {'fig1': '<markdown caption paragraph>', ...}"""
    caps: dict[str, str] = {}
    for block in re.split(r"\n\s*\n", fig_md.strip()):
        block = block.strip()
        if not block or block.startswith("---"):
            continue
        m = re.match(r"\*\*Figure (\d+)\.", block)
        if m:
            caps[f"fig{m.group(1)}"] = " ".join(l.strip() for l in block.split("\n"))
    return caps


def figure_html(key: str, caption_md: str) -> str:
    png = FIGDIR / f"{key}.png"
    data = base64.b64encode(png.read_bytes()).decode("ascii")
    from PIL import Image  # noqa: PLC0415

    with Image.open(png) as im:
        aspect = im.size[1] / im.size[0]
    width_mm = min(TEXT_W_MM, FIG_MAX_H_MM / aspect)
    return (
        f'<figure class="fig" id="{key}">'
        f'<img style="width:{width_mm:.1f}mm" alt="{key}" '
        f'src="data:image/png;base64,{data}"/>'
        f'<figcaption>{inline(caption_md)}</figcaption>'
        f"</figure>"
    )


def insert_figures(body_md: str) -> str:
    blocks = re.split(r"\n\s*\n", body_md.replace("\r\n", "\n"))

    for key, anchor in FIG_PLACEMENT:
        marker = f"@@FIGURE:{key}@@"
        if anchor[0] == "before":
            idx = next(i for i, b in enumerate(blocks) if b.startswith(anchor[1]))
            blocks.insert(idx, marker)
        else:
            _, section, needle = anchor
            start = next(i for i, b in enumerate(blocks) if b.startswith(section))
            idx = next(
                i
                for i in range(start, len(blocks))
                if needle in blocks[i] and not blocks[i].startswith("#")
            )
            blocks.insert(idx + 1, marker)
    return "\n\n".join(blocks)


def extract_reference_list(refs_md: str) -> str:
    """Only the verified reference list; the repo-internal citation audit
    apparatus (which quotes the string [CITATION NEEDED]) is not part of the
    manuscript and is excluded."""
    lines = refs_md.replace("\r\n", "\n").split("\n")
    i0 = next(i for i, l in enumerate(lines) if l.startswith("## Verified references"))
    i1 = next(i for i, l in enumerate(lines) if i > i0 and l.startswith("### Notes on citation"))
    block = "\n".join(lines[i0 + 1 : i1])
    block = re.sub(r"^\s*---\s*$", "", block, flags=re.M)
    return block.strip()


# ========================================================================== CSS
CSS = f"""
@page {{ size: {PAGE_W_MM}mm {PAGE_H_MM}mm; }}

html {{ -webkit-print-color-adjust: exact; print-color-adjust: exact; }}

body {{
  width: {CONTENT_W_MM}mm;
  padding-left: {GUTTER_MM}mm;
  margin: 0;
  font-family: Cambria, Georgia, "Times New Roman", serif;
  font-size: 11pt;
  line-height: 1.5;
  color: #111;
  background: #fff;
  text-align: left;
  hyphens: none;
}}

pre, .vline {{ position: relative; }}
.vline {{ margin: 0; padding: 0; }}

p {{ margin: 0 0 .55em 0; }}

h1, h2, h3, h4 {{
  font-family: "Segoe UI", Calibri, Arial, sans-serif;
  color: #000;
  break-after: avoid; page-break-after: avoid;
  line-height: 1.25;
}}
h2 {{ font-size: 14pt; margin: 1.5em 0 .5em; padding-bottom: .15em;
      border-bottom: .6pt solid #bbb; }}
h3 {{ font-size: 11.5pt; margin: 1.2em 0 .4em; }}

ul, ol {{ margin: 0 0 .6em 0; padding-left: 1.5em; }}
li {{ margin-bottom: .4em; }}

strong {{ font-weight: 700; }}
code {{ font-family: Consolas, "Courier New", monospace; font-size: .87em;
        background: #f2f2f2; padding: 0 .18em; border-radius: 2px; }}
pre {{ background: #f5f5f5; border-left: 2pt solid #bbb; padding: .5em .7em;
       margin: .6em 0 .8em; white-space: pre-wrap; word-break: break-word;
       break-inside: avoid; page-break-inside: avoid; }}
pre code {{ background: none; padding: 0; font-size: 8.5pt; line-height: 1.35; }}

/* ------------------------------------------------------------------ title */
.titleblock {{ margin: 0 0 1.4em 0; padding-bottom: .8em;
               border-bottom: 1.2pt solid #333; }}
.titleblock h1 {{ font-size: 17pt; line-height: 1.28; margin: 0 0 .55em 0;
                  font-weight: 600; letter-spacing: -.1pt; }}
.author {{ font-size: 12pt; font-weight: 700; margin: 0 0 .15em; }}
.meta {{ font-size: 9.5pt; line-height: 1.45; margin: 0; color: #222;
         font-family: "Segoe UI", Calibri, Arial, sans-serif; }}
.meta a {{ color: #222; text-decoration: none; }}
.archive {{ margin-top: .6em; font-size: 9pt; color: #333;
            font-family: "Segoe UI", Calibri, Arial, sans-serif; }}
.notpeer {{ margin-top: .5em; font-size: 9pt; font-style: italic; color: #555;
            font-family: "Segoe UI", Calibri, Arial, sans-serif; }}

/* ----------------------------------------------------------------- tables */
.tablewrap {{ break-inside: avoid; page-break-inside: avoid; margin: .5em 0 1em; }}
.tableblock {{ break-inside: avoid; page-break-inside: avoid; }}
table {{ border-collapse: collapse; width: 100%; table-layout: auto;
         font-family: "Segoe UI", Calibri, Arial, sans-serif;
         font-size: 8.3pt; line-height: 1.28; }}
th, td {{ border: .5pt solid #999; padding: 2.6pt 3.4pt; text-align: left;
          vertical-align: top; overflow-wrap: break-word; word-break: normal; }}
th {{ background: #ececec; font-weight: 700; }}
tbody tr:nth-child(even) td {{ background: #fafafa; }}
td code, th code {{ font-size: .95em; background: none; padding: 0; }}

/* ---------------------------------------------------------------- figures */
figure.fig {{ margin: 1.1em 0 1.2em; break-inside: avoid; page-break-inside: avoid;
              text-align: center; }}
figure.fig img {{ display: block; margin: 0 auto .5em; height: auto; }}
figcaption {{ font-size: 9pt; line-height: 1.42; text-align: left;
              font-family: "Segoe UI", Calibri, Arial, sans-serif; color: #1a1a1a; }}

/* ------------------------------------------------------------- references */
#references ol {{ padding-left: 1.9em; }}
#references li {{ font-size: 9.5pt; line-height: 1.4; margin-bottom: .5em;
                  overflow-wrap: break-word; }}

.figures-note {{ font-size: 10pt; }}

/* ----------------------------------------------------------- line numbers */
.lineno {{
  position: absolute;
  left: -{GUTTER_MM}mm;
  width: {GUTTER_MM - 2.2}mm;
  text-align: right;
  font-family: "Segoe UI", Calibri, Arial, sans-serif;
  font-size: 7pt;
  line-height: 1;
  color: #9a9a9a;
  user-select: none;
}}
"""

# JS: number every visual line of body text.
#
# Naive approach -- one absolutely positioned number per line box inside the
# paragraph -- breaks down whenever a paragraph is split across a page: the
# number whose offset falls in the leftover strip at the foot of the page stays
# there while its line moves to the next page, shifting every number in the
# continuation by one line.  Instead each *visual* line is promoted to its own
# block element, measured first at exactly the printed text width, so page
# breaks can only ever fall between numbered lines and a number can never come
# adrift from its text.  Line breaks are unchanged: the body is a fixed
# millimetre width and the print scale is 1, so screen and print layout agree.
LINENO_JS = r"""
(function () {
  var n = 0;
  // Every number must sit at the same absolute x, whatever the indent of the
  // block it belongs to (list items sit inside the <ol> padding).
  var bodyCS = getComputedStyle(document.body);
  var bodyLeft = document.body.getBoundingClientRect().left +
                 parseFloat(bodyCS.paddingLeft);
  var GUTTER_PX = parseFloat(bodyCS.paddingLeft);

  function lineStarts(el, lh) {
    var walker = document.createTreeWalker(el, NodeFilter.SHOW_TEXT, null);
    var starts = [], lineTop = null, r = document.createRange(), node;
    while ((node = walker.nextNode())) {
      var t = node.nodeValue;
      for (var i = 0; i < t.length; i++) {
        if (/\s/.test(t[i])) { continue; }
        r.setStart(node, i); r.setEnd(node, i + 1);
        var rect = r.getBoundingClientRect();
        if (rect.height <= 0) { continue; }
        if (lineTop === null || rect.top - lineTop > lh * 0.6) {
          starts.push([node, i]);
          lineTop = rect.top;
        }
      }
    }
    return starts;
  }

  function splitIntoLines(el) {
    var cs = getComputedStyle(el);
    var lh = parseFloat(cs.lineHeight);
    if (!isFinite(lh)) { lh = parseFloat(cs.fontSize) * 1.5; }
    var starts = lineStarts(el, lh);
    if (starts.length === 0) { return; }

    // Ranges are built against the untouched DOM, then cloned, so partially
    // selected <strong>/<em>/<code> spans keep their markup on both sides.
    var frags = [];
    for (var i = 0; i < starts.length; i++) {
      var r = document.createRange();
      r.setStart(starts[i][0], starts[i][1]);
      if (i + 1 < starts.length) { r.setEnd(starts[i + 1][0], starts[i + 1][1]); }
      else { r.setEnd(el, el.childNodes.length); }
      // cloneContents() only preserves markup from commonAncestorContainer
      // downwards, so a line that falls entirely inside a <strong>/<em>/<code>
      // would come back as bare text.  Re-wrap it in shallow clones of the
      // ancestors between that container and the block itself.
      var frag = r.cloneContents();
      var anc = r.commonAncestorContainer;
      if (anc.nodeType === 3) { anc = anc.parentNode; }
      while (anc && anc !== el && el.contains(anc)) {
        var w = anc.cloneNode(false);
        w.appendChild(frag);
        frag = document.createDocumentFragment();
        frag.appendChild(w);
        anc = anc.parentNode;
      }
      frags.push(frag);
    }

    var indent = el.getBoundingClientRect().left - bodyLeft;
    while (el.firstChild) { el.removeChild(el.firstChild); }
    for (var j = 0; j < frags.length; j++) {
      var d = document.createElement('div');
      d.className = 'vline';
      var s = document.createElement('span');
      s.className = 'lineno';
      s.textContent = String(++n);
      s.style.left = (-GUTTER_PX - indent) + 'px';
      d.appendChild(s);
      d.appendChild(frags[j]);
      el.appendChild(d);
    }
    // emulate orphans:2 / widows:2, which block-level splitting would lose
    if (frags.length > 2) {
      el.firstChild.style.breakAfter = 'avoid';
      el.firstChild.style.pageBreakAfter = 'avoid';
      el.lastChild.style.breakBefore = 'avoid';
      el.lastChild.style.pageBreakBefore = 'avoid';
    }
  }

  var blocks = Array.prototype.slice.call(
    document.querySelectorAll('p, li, figcaption')
  ).filter(function (el) {
    return !el.closest('.titleblock') && !el.closest('table') &&
           el.textContent.trim().length > 0;
  });
  blocks.forEach(splitIntoLines);
  return n;
})();
"""


# ======================================================================== build
def build_html(with_line_numbers: bool) -> str:
    md = MS.read_text(encoding="utf-8")
    front, body_md, fig_md, tail_md = split_sections(md)

    title = re.match(r"^#\s+(.*)$", md.split("\n", 1)[0]).group(1).strip()
    caps = parse_captions(fig_md)
    missing = {f"fig{i}" for i in range(1, 6)} - set(caps)
    if missing:
        raise SystemExit(f"missing captions for {sorted(missing)}")

    body_html = bind_table_captions(md_to_html(insert_figures(body_md)))
    for key in sorted(caps):
        body_html = body_html.replace(
            f"<p>@@FIGURE:{key}@@</p>", figure_html(key, caps[key])
        )
    if "@@FIGURE" in body_html:
        raise SystemExit("figure placeholder left unsubstituted")

    figures_section = (
        "<h2>6. Figures</h2>"
        f'<p class="figures-note">{inline(FIGURES_SECTION_NOTE)}</p>'
    )

    tail_html = md_to_html(tail_md)
    refs_html = md_to_html(extract_reference_list(REFS.read_text(encoding="utf-8")))

    titleblock = f"""
<header class="titleblock">
  <h1>{inline(title)}</h1>
  <p class="author">Ivan Heredia Jalid</p>
  <p class="meta">Independent Researcher, C&oacute;rdoba, Argentina<br/>
     ORCID 0009-0003-6702-1295<br/>
     Correspondence: ivanjalid@gmail.com</p>
  <p class="archive">Archived code and data-processing pipeline:
     doi:10.5281/zenodo.22664837
     (https://doi.org/10.5281/zenodo.22664837);
     repository https://github.com/ivanjalid1/hsv1-crispr-escape</p>
  <p class="notpeer">Preprint. Not peer reviewed.</p>
</header>
"""

    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"/>
<title>{esc(title)}</title>
<style>{CSS}</style>
</head><body>
{titleblock}
{body_html}
{figures_section}
{tail_html}
<h2 id="references-h">11. References</h2>
<div id="references">{refs_html}</div>
</body></html>
"""
    return doc


# ===================================================================== chrome/cdp
def find_chrome() -> str:
    for c in CHROME_CANDIDATES:
        if os.path.exists(c):
            return c
    raise SystemExit("no Chrome/Edge binary found")


def free_port() -> int:
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def print_pdf(html_path: Path, pdf_path: Path, run_lineno: bool) -> None:
    import asyncio  # noqa: PLC0415

    import websockets  # noqa: PLC0415

    chrome = find_chrome()
    port = free_port()
    profile = Path(tempfile.mkdtemp(prefix="chrome-pdf-"))
    proc = subprocess.Popen(
        [
            chrome,
            "--headless=new",
            "--disable-gpu",
            "--no-sandbox",
            "--no-first-run",
            "--no-default-browser-check",
            "--disable-extensions",
            "--hide-scrollbars",
            "--force-device-scale-factor=1",
            "--window-size=1400,2000",
            f"--user-data-dir={profile}",
            f"--remote-debugging-port={port}",
            "about:blank",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    try:
        ws_url = None
        for _ in range(120):
            try:
                with urllib.request.urlopen(f"http://127.0.0.1:{port}/json/version", timeout=1) as r:
                    ws_url = json.load(r)["webSocketDebuggerUrl"]
                break
            except Exception:
                time.sleep(0.25)
        if not ws_url:
            raise SystemExit("chrome devtools endpoint never came up")

        async def run() -> None:
            async with websockets.connect(ws_url, max_size=200 * 1024 * 1024) as ws:
                nid = 0

                async def cmd(method, params=None, session=None):
                    nonlocal nid
                    nid += 1
                    msg = {"id": nid, "method": method, "params": params or {}}
                    if session:
                        msg["sessionId"] = session
                    await ws.send(json.dumps(msg))
                    while True:
                        resp = json.loads(await ws.recv())
                        if resp.get("id") == nid:
                            if "error" in resp:
                                raise SystemExit(f"{method}: {resp['error']}")
                            return resp.get("result", {})

                t = await cmd("Target.createTarget", {"url": "about:blank"})
                sess = (await cmd("Target.attachToTarget",
                                  {"targetId": t["targetId"], "flatten": True}))["sessionId"]
                await cmd("Page.enable", session=sess)
                await cmd("Runtime.enable", session=sess)
                url = html_path.resolve().as_uri()
                await cmd("Page.navigate", {"url": url}, session=sess)

                # wait for load + fonts + images
                deadline = time.time() + 60
                while time.time() < deadline:
                    r = await cmd(
                        "Runtime.evaluate",
                        {
                            "expression": (
                                "(document.readyState==='complete') && "
                                "Array.from(document.images).every(i=>i.complete) ? "
                                "document.fonts.ready.then(()=>true) : false"
                            ),
                            "awaitPromise": True,
                            "returnByValue": True,
                        },
                        session=sess,
                    )
                    if r.get("result", {}).get("value") is True:
                        break
                    time.sleep(0.25)
                time.sleep(0.6)

                if run_lineno:
                    r = await cmd(
                        "Runtime.evaluate",
                        {"expression": LINENO_JS, "returnByValue": True},
                        session=sess,
                    )
                    print(f"  line numbers placed: {r['result'].get('value')}")

                footer = (
                    '<div style="width:100%;font-family:Georgia,serif;font-size:8px;'
                    'color:#555;text-align:center;margin:0 0 0 0;">'
                    '<span class="pageNumber"></span> of '
                    '<span class="totalPages"></span></div>'
                )
                res = await cmd(
                    "Page.printToPDF",
                    {
                        "landscape": False,
                        "printBackground": True,
                        "scale": 1,
                        "paperWidth": PAGE_W_MM / MM_PER_IN,
                        "paperHeight": PAGE_H_MM / MM_PER_IN,
                        "marginTop": MARGIN_TOP_MM / MM_PER_IN,
                        "marginBottom": MARGIN_BOTTOM_MM / MM_PER_IN,
                        "marginLeft": MARGIN_LEFT_MM / MM_PER_IN,
                        "marginRight": MARGIN_RIGHT_MM / MM_PER_IN,
                        "displayHeaderFooter": True,
                        "headerTemplate": "<div></div>",
                        "footerTemplate": footer,
                        "preferCSSPageSize": False,
                        "transferMode": "ReturnAsBase64",
                    },
                    session=sess,
                )
                pdf_path.write_bytes(base64.b64decode(res["data"]))

        asyncio.run(run())
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=15)
        except Exception:
            proc.kill()
        shutil.rmtree(profile, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-html", action="store_true", help="keep the intermediate HTML")
    ap.add_argument("--no-line-numbers", action="store_true")
    args = ap.parse_args()

    print("building HTML ...")
    doc = build_html(not args.no_line_numbers)
    OUT_HTML.write_text(doc, encoding="utf-8")
    print(f"  {OUT_HTML.name}: {len(doc):,} bytes")

    print("printing with headless Chrome ...")
    print_pdf(OUT_HTML, OUT_PDF, run_lineno=not args.no_line_numbers)
    print(f"  {OUT_PDF.name}: {OUT_PDF.stat().st_size:,} bytes")

    if not args.keep_html:
        OUT_HTML.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
