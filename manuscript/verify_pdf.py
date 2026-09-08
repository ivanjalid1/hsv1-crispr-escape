#!/usr/bin/env python3
"""Verify the built PDF against the Markdown sources.

Checks:
  1. every paragraph, list item, heading, code line and table cell of
     manuscript.md is present in the PDF text, verbatim after normalisation;
  2. every verified reference in references.md is present;
  3. no Markdown syntax, placeholder or build marker leaked through;
  4. all fonts are embedded subsets, every figure is present exactly once,
     nothing is drawn outside the intended text block, and the line numbers
     form one unbroken 1..N sequence in reading order.

Usage:  .venv/Scripts/python.exe manuscript/verify_pdf.py
"""

from __future__ import annotations

import re
import sys
import unicodedata
from pathlib import Path

import pymupdf

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parent.parent
PDF = ROOT / "manuscript" / "hsv1-crispr-escape-preprint.pdf"
MS = ROOT / "manuscript" / "manuscript.md"
REFS = ROOT / "manuscript" / "references.md"

MM = 72 / 25.4
TEXT_LEFT_MM, TEXT_RIGHT_MM = 24.0, 186.0
TEXT_TOP_MM, TEXT_BOTTOM_MM = 20.0, 278.0
GUTTER_MAX_MM = 23.0


def norm(s: str) -> str:
    """Normalise for comparison: strip markdown, unify dashes/quotes/space."""
    s = unicodedata.normalize("NFKC", s)
    s = re.sub(r"[`*]", "", s)
    s = s.replace("’", "'").replace("‘", "'")
    s = s.replace("“", '"').replace("”", '"')
    s = re.sub(r"[‐-―−]", "-", s)
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def nows(s: str) -> str:
    """Whitespace-free form.  PDF text extraction inserts and drops spaces at
    span and line boundaries, so comparison is done on the character sequence
    with all whitespace removed; line breaks only ever fall at an existing
    space or hyphen (hyphenation is off), so this is lossless for the check."""
    return re.sub(r"\s+", "", norm(s))


def pdf_body_text(doc: pymupdf.Document) -> str:
    """All page text except the line-number gutter and the page footer."""
    chunks = []
    for page in doc:
        for block in page.get_text("dict")["blocks"]:
            if block["type"] != 0:
                continue
            for line in block["lines"]:
                for span in line["spans"]:
                    x0, y0, x1, y1 = span["bbox"]
                    if x1 / MM <= GUTTER_MAX_MM:      # line number
                        continue
                    if y1 / MM > TEXT_BOTTOM_MM + 4:  # footer
                        continue
                    chunks.append(span["text"])
                chunks.append(" ")
    return norm(" ".join(chunks))


def source_fragments() -> list[tuple[str, str]]:
    """(label, text) fragments that must all appear in the PDF."""
    md = MS.read_text(encoding="utf-8").replace("\r\n", "\n")
    lines = md.split("\n")
    i_abs = next(i for i, l in enumerate(lines) if l.startswith("## Abstract"))
    i_fig = next(i for i, l in enumerate(lines) if l.startswith("## 6. Figures"))
    i_after = next(i for i, l in enumerate(lines) if l.startswith("## 7. "))
    # Section 6 holds the figure captions; they are printed with the figures.
    body = lines[i_abs:i_fig] + lines[i_fig + 1 : i_after] + lines[i_after:]

    frags: list[tuple[str, str]] = []
    para: list[str] = []
    in_code = False
    for ln in body:
        s = ln.strip()
        if s.startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            if s:
                frags.append(("code", s))
            continue
        if s.startswith("|"):
            if re.match(r"^\|[\s:|-]+\|?$", s):
                continue
            for cell in s.strip("|").split("|"):
                if cell.strip():
                    frags.append(("cell", cell.strip()))
            continue
        if not s or s.startswith("#") or re.fullmatch(r"-{3,}", s):
            if para:
                frags.append(("para", " ".join(para)))
                para = []
            if s.startswith("#"):
                frags.append(("heading", s.lstrip("# ")))
            continue
        m = re.match(r"^(?:[-*+]|\d+\.)\s+(.*)$", s)
        if m:                       # a new list item starts its own fragment
            if para:
                frags.append(("para", " ".join(para)))
            para = [m.group(1)]
        else:
            para.append(s)
    if para:
        frags.append(("para", " ".join(para)))

    # references: the verified list only
    rl = REFS.read_text(encoding="utf-8").replace("\r\n", "\n").split("\n")
    r0 = next(i for i, l in enumerate(rl) if l.startswith("## Verified references"))
    r1 = next(i for i, l in enumerate(rl) if i > r0 and l.startswith("### Notes on citation"))
    item: list[str] = []
    for ln in rl[r0 + 1 : r1]:
        if re.match(r"^\d+\.\s", ln.strip()) or not ln.strip():
            if item:
                frags.append(("ref", " ".join(item)))
                item = []
        if ln.strip() and not re.fullmatch(r"-{3,}", ln.strip()):
            item.append(re.sub(r"^\d+\.\s*", "", ln.strip()))
    if item:
        frags.append(("ref", " ".join(item)))

    return [(k, norm(v)) for k, v in frags if len(norm(v)) >= 12]


def main() -> int:
    if not PDF.exists():
        print("FAIL: PDF missing")
        return 1

    doc = pymupdf.open(PDF)
    body = pdf_body_text(doc)
    problems: list[str] = []

    # 1/2 -- content completeness
    frags = source_fragments()
    body_ns = nows(body)

    def present(text: str) -> bool:
        """Containment, tolerant of page-break interleaving.

        PyMuPDF emits blocks page by page, so a paragraph continuing over a
        page break is interrupted in the flat stream by that page's heading or
        by the list markers PyMuPDF groups into one block.  The fragment is
        therefore consumed greedily: take the longest prefix that occurs in the
        PDF, drop it, repeat.  Text that is really present survives in a few
        long chunks; text that is missing leaves a short unmatched piece.
        """
        rest, chunks = nows(text), 0
        while rest:
            lo, hi = 0, len(rest)
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if rest[:mid] in body_ns:
                    lo = mid
                else:
                    hi = mid - 1
            if lo < 20 and len(rest) >= 20:
                return False
            rest = rest[max(lo, 1):]
            chunks += 1
            if chunks > 4:
                return False
        return True

    missing = [(k, v) for k, v in frags if not present(v)]
    for k, v in missing:
        problems.append(f"missing {k}: {v[:110]!r}")
    print(f"content fragments checked : {len(frags)}  missing: {len(missing)}")

    # 1b -- nothing dropped: order-independent character multiset.  Any text
    # lost in conversion would lower some character's count below the source.
    from collections import Counter
    src_chars = Counter(nows(" ".join(v for _, v in frags)))
    pdf_chars = Counter(body_ns)
    short = {c: n - pdf_chars[c] for c, n in src_chars.items() if pdf_chars[c] < n}
    print(f"characters (source/pdf)   : {sum(src_chars.values()):,} / "
          f"{sum(pdf_chars.values()):,}  deficient classes: {len(short)}")
    for c, d in list(short.items())[:10]:
        problems.append(f"character {c!r} appears {d} time(s) fewer in the PDF")

    # 3 -- no leaked syntax or placeholders
    raw = "".join(p.get_text() for p in doc)
    for bad in ("[CITATION NEEDED]", "**", "@@FIGURE", "|---", "\x00"):
        if bad in raw:
            problems.append(f"leaked marker in PDF: {bad!r}")
    if re.search(r"^\s*#{1,6}\s", raw, re.M):
        problems.append("leaked ATX heading syntax")

    # 4a -- fonts
    fonts = {(f[3], f[1]) for p in doc for f in p.get_fonts(full=True)}
    unembedded = [n for n, ext in fonts if not re.match(r"^[A-Z]{6}\+", n) or ext == "n/a"]
    print(f"fonts                     : {len(fonts)} all subset-embedded: {not unembedded}")
    for n in unembedded:
        problems.append(f"font not embedded: {n}")

    # 4b -- figures
    imgs = sum(len(p.get_images()) for p in doc)
    print(f"embedded figure images    : {imgs}")
    if imgs != 5:
        problems.append(f"expected 5 figure images, found {imgs}")
    for n in range(1, 6):
        if f"Figure {n}." not in raw:
            problems.append(f"caption for Figure {n} missing")

    # 4c -- geometry
    for i, page in enumerate(doc, 1):
        boxes = [b[:4] for b in page.get_text("blocks") if b[3] / MM < TEXT_BOTTOM_MM + 4]
        boxes += [im["bbox"] for im in page.get_image_info()]
        if not boxes:
            problems.append(f"page {i} is empty")
            continue
        if max(b[2] for b in boxes) / MM > TEXT_RIGHT_MM:
            problems.append(f"page {i} overflows the right margin")
        if min(b[0] for b in boxes) / MM < 15.0:
            problems.append(f"page {i} overflows the left margin")
        if max(b[3] for b in boxes) / MM > TEXT_BOTTOM_MM:
            problems.append(f"page {i} overflows the bottom margin")

    # 4d -- line numbers form 1..N in reading order
    nums: list[int] = []
    for page in doc:
        page_nums = []
        for block in page.get_text("dict")["blocks"]:
            if block["type"] != 0:
                continue
            for line in block["lines"]:
                for span in line["spans"]:
                    if span["bbox"][2] / MM <= GUTTER_MAX_MM and span["text"].strip().isdigit():
                        page_nums.append((span["bbox"][1], int(span["text"])))
        nums += [n for _, n in sorted(page_nums)]
    expected = list(range(1, len(nums) + 1))
    print(f"line numbers              : {len(nums)} contiguous 1..N: {nums == expected}")
    if nums != expected:
        bad = next((i for i, (a, b) in enumerate(zip(nums, expected)) if a != b), None)
        problems.append(f"line numbering breaks at index {bad}: {nums[max(0,(bad or 0)-2):(bad or 0)+3]}")

    print(f"pages                     : {doc.page_count}")
    print(f"file size                 : {PDF.stat().st_size:,} bytes "
          f"({PDF.stat().st_size / 1048576:.2f} MB)")

    if problems:
        print(f"\n{len(problems)} PROBLEM(S):")
        for p in problems[:40]:
            print("  -", p)
        return 1
    print("\nOK: all checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
