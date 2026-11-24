#!/usr/bin/env python3
"""Build 6 Korean-language review PDFs from versioned Markdown sources."""
import argparse
import json
import re
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, SimpleDocTemplate

FONT_FILES = (
    "/usr/share/fonts/truetype/nanum/NanumGothic.ttf",
    "/usr/share/fonts/truetype/nanum/NanumBarunGothic.ttf",
)


def register_font():
    for path in FONT_FILES:
        if Path(path).is_file():
            pdfmetrics.registerFont(TTFont("PacketKorean", path))
            return
    raise FileNotFoundError("Korean font missing; on Ubuntu install fonts-nanum")


def render(text: str, path: Path):
    body = ParagraphStyle("body", fontName="PacketKorean", fontSize=9.5, leading=15, spaceAfter=7, wordWrap="CJK")
    title = ParagraphStyle("title", parent=body, fontSize=16, leading=23, spaceAfter=12)
    heading = ParagraphStyle("heading", parent=body, fontSize=12, leading=19, spaceBefore=13, spaceAfter=6, textColor=colors.HexColor("#183548"))
    bullet = ParagraphStyle("bullet", parent=body, leftIndent=14, firstLineIndent=-8)
    story = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("<!--"):
            continue
        if line.startswith("# "):
            style, line = title, line[2:]
        elif line.startswith("## "):
            style, line = heading, line[3:]
        elif line.startswith("- "):
            style, line = bullet, "• " + line[2:]
        else:
            style = body
        line = re.sub(r"\*\*(.*?)\*\*", r"\1", line)
        story.append(Paragraph(escape(line), style))
    if not story:
        raise ValueError(f"empty rendered proposal: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(str(path), pagesize=A4, topMargin=42, bottomMargin=42,
                            leftMargin=41, rightMargin=41, title=path.stem)
    doc.build(story)


def main():
    args = argparse.ArgumentParser()
    args.add_argument("--root", type=Path, default=Path("."))
    args.add_argument("--output", type=Path, default=Path("output/proposal-packet"))
    a = args.parse_args()
    root = a.root.resolve()
    entries = json.loads((root / "proposals/2026-11/manifest.json").read_text(encoding="utf-8"))["proposals"]
    register_font()
    for e in entries:
        output = root / a.output / f"{e['id']}_{e['slug']}.pdf"
        render((root / e["proposal"]).read_text(encoding="utf-8"), output)
        print(f"CREATED {output.name} : {output.stat().st_size} bytes")


if __name__ == "__main__":
    main()
