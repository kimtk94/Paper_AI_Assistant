#!/usr/bin/env python3
"""Organize downloaded PDFs by disease and publication year.

Features
- Robust text extraction with multiple backends (PyPDF, pdfplumber, PyMuPDF, pdftotext).
- Disease classification from extracted text using keyword rules.
- Year inference from filename + extracted text.
- File organization into <output>/<disease>/<year>/.
- Machine-friendly JSON report + human-friendly Markdown summary.

Example
    python src/pdf_disease_year_organizer.py \
      --pdf-dir data/pdfs \
      --output-dir data/organized_pdfs \
      --report-json data/pdf_organize_report.json \
      --summary-md data/pdf_organize_summary.md
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable
import importlib
import importlib.util


DEFAULT_DISEASE_KEYWORDS: dict[str, list[str]] = {
    "type_2_diabetes": [
        "type 2 diabetes",
        "t2d",
        "insulin resistance",
        "hyperglycemia",
        "metformin",
        "hba1c",
    ],
    "alzheimers_disease": [
        "alzheimer",
        "amyloid",
        "tauopathy",
        "mci",
        "dementia",
        "adni",
    ],
    "parkinsons_disease": [
        "parkinson",
        "alpha-synuclein",
        "dopaminergic",
        "pd patients",
    ],
    "cardiovascular_disease": [
        "cardiovascular",
        "atherosclerosis",
        "myocardial",
        "heart failure",
        "coronary artery",
    ],
    "cancer": [
        "cancer",
        "tumor",
        "oncology",
        "carcinoma",
        "metastasis",
        "chemotherapy",
    ],
}


@dataclass
class ExtractionResult:
    backend: str
    text: str


@dataclass
class PdfClassification:
    filename: str
    source_path: str
    backend: str
    text_length: int
    year: int | None
    diseases: list[str]
    destination_paths: list[str]
    extraction_error: str | None


def _load_backend_function() -> tuple[str, Callable[[Path], str]]:
    if importlib.util.find_spec("pypdf") is not None:
        module = importlib.import_module("pypdf")

        def extract_with_pypdf(pdf_path: Path) -> str:
            reader = module.PdfReader(str(pdf_path))
            return "\n".join((page.extract_text() or "") for page in reader.pages)

        return "pypdf", extract_with_pypdf

    if importlib.util.find_spec("pdfplumber") is not None:
        module = importlib.import_module("pdfplumber")

        def extract_with_pdfplumber(pdf_path: Path) -> str:
            texts: list[str] = []
            with module.open(str(pdf_path)) as pdf:
                for page in pdf.pages:
                    texts.append(page.extract_text() or "")
            return "\n".join(texts)

        return "pdfplumber", extract_with_pdfplumber

    if importlib.util.find_spec("fitz") is not None:
        module = importlib.import_module("fitz")

        def extract_with_pymupdf(pdf_path: Path) -> str:
            texts: list[str] = []
            doc = module.open(str(pdf_path))
            try:
                for page in doc:
                    texts.append(page.get_text("text") or "")
            finally:
                doc.close()
            return "\n".join(texts)

        return "pymupdf", extract_with_pymupdf

    pdftotext_bin = shutil.which("pdftotext")
    if pdftotext_bin:

        def extract_with_pdftotext(pdf_path: Path) -> str:
            result = subprocess.run(
                [pdftotext_bin, "-layout", str(pdf_path), "-"],
                check=False,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                stderr = (result.stderr or "").strip()
                raise RuntimeError(f"pdftotext failed (code={result.returncode}): {stderr}")
            return result.stdout or ""

        return "pdftotext", extract_with_pdftotext

    raise RuntimeError(
        "No PDF extractor available. Install one: pypdf (recommended), pdfplumber, pymupdf, or pdftotext CLI."
    )


def extract_text(pdf_path: Path) -> ExtractionResult:
    backend, extractor = _load_backend_function()
    text = extractor(pdf_path)
    return ExtractionResult(backend=backend, text=_normalize_whitespace(text))


def _normalize_whitespace(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def infer_year(pdf_path: Path, text: str) -> int | None:
    current_year = datetime.now().year
    valid_years = set(range(1900, current_year + 1))

    filename_matches = re.findall(r"(?:19|20)\d{2}", pdf_path.name)
    for cand in filename_matches:
        value = int(cand)
        if value in valid_years:
            return value

    text_slice = text[:20000]
    text_matches = re.findall(r"(?:19|20)\d{2}", text_slice)
    for cand in sorted(text_matches, reverse=True):
        value = int(cand)
        if value in valid_years:
            return value
    return None


def infer_diseases(text: str, disease_keywords: dict[str, list[str]]) -> list[str]:
    lowered = text.lower()
    hits: list[str] = []
    for disease, keywords in disease_keywords.items():
        if any(keyword.lower() in lowered for keyword in keywords):
            hits.append(disease)
    return sorted(set(hits))


def load_disease_keywords(custom_path: Path | None) -> dict[str, list[str]]:
    if not custom_path:
        return DEFAULT_DISEASE_KEYWORDS
    payload = json.loads(custom_path.read_text(encoding="utf-8"))
    normalized: dict[str, list[str]] = {}
    for disease, keywords in payload.items():
        if not isinstance(keywords, list):
            raise ValueError(f"Keyword list must be array: disease={disease}")
        normalized[str(disease)] = [str(v) for v in keywords]
    return normalized


def organize_pdfs(
    *,
    pdf_dir: Path,
    output_dir: Path,
    text_dir: Path,
    disease_keywords: dict[str, list[str]],
) -> list[PdfClassification]:
    output_dir.mkdir(parents=True, exist_ok=True)
    text_dir.mkdir(parents=True, exist_ok=True)

    records: list[PdfClassification] = []
    for pdf_path in sorted(pdf_dir.glob("*.pdf")):
        try:
            extraction = extract_text(pdf_path)
            text = extraction.text
            backend = extraction.backend
            error_message = None
        except Exception as exc:
            records.append(
                PdfClassification(
                    filename=pdf_path.name,
                    source_path=str(pdf_path),
                    backend="n/a",
                    text_length=0,
                    year=None,
                    diseases=[],
                    destination_paths=[],
                    extraction_error=str(exc),
                )
            )
            continue

        year = infer_year(pdf_path, text)
        diseases = infer_diseases(text, disease_keywords)
        if not diseases:
            diseases = ["unknown_disease"]

        year_label = str(year) if year is not None else "unknown_year"

        txt_out = text_dir / f"{pdf_path.stem}.txt"
        txt_out.write_text(text, encoding="utf-8")

        destinations: list[str] = []
        for disease in diseases:
            disease_dir = output_dir / disease / year_label
            disease_dir.mkdir(parents=True, exist_ok=True)
            dst = disease_dir / pdf_path.name
            shutil.copy2(pdf_path, dst)
            destinations.append(str(dst))

        records.append(
            PdfClassification(
                filename=pdf_path.name,
                source_path=str(pdf_path),
                backend=backend,
                text_length=len(text),
                year=year,
                diseases=diseases,
                destination_paths=destinations,
                extraction_error=error_message,
            )
        )

    return records


def write_reports(records: list[PdfClassification], report_json: Path, summary_md: Path) -> None:
    report_json.parent.mkdir(parents=True, exist_ok=True)
    summary_md.parent.mkdir(parents=True, exist_ok=True)

    report_json.write_text(
        json.dumps([asdict(record) for record in records], ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    total = len(records)
    success = sum(1 for r in records if r.extraction_error is None)
    failed = total - success

    disease_counter: dict[str, int] = {}
    year_counter: dict[str, int] = {}
    for record in records:
        for disease in record.diseases:
            disease_counter[disease] = disease_counter.get(disease, 0) + 1
        year_key = str(record.year) if record.year is not None else "unknown_year"
        year_counter[year_key] = year_counter.get(year_key, 0) + 1

    lines = [
        "# PDF 정리 요약",
        "",
        f"- 전체 PDF: **{total}**",
        f"- 텍스트 추출 성공: **{success}**",
        f"- 텍스트 추출 실패: **{failed}**",
        "",
        "## 질병별 분포",
        "",
        "| disease | count |",
        "|---|---:|",
    ]
    for disease, count in sorted(disease_counter.items(), key=lambda x: x[0]):
        lines.append(f"| {disease} | {count} |")

    lines.extend(["", "## 연도별 분포", "", "| year | count |", "|---|---:|"])
    for year, count in sorted(year_counter.items(), key=lambda x: x[0]):
        lines.append(f"| {year} | {count} |")

    lines.extend(["", "## 실패 파일", ""])
    failed_rows = [r for r in records if r.extraction_error]
    if not failed_rows:
        lines.append("- 없음")
    else:
        for item in failed_rows:
            lines.append(f"- {item.filename}: {item.extraction_error}")

    summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Organize PDFs by disease and publication year.")
    parser.add_argument("--pdf-dir", type=Path, required=True, help="Input directory containing downloaded PDF files.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/organized_pdfs"),
        help="Output directory. PDFs are copied into <output>/<disease>/<year>/.",
    )
    parser.add_argument(
        "--text-dir",
        type=Path,
        default=Path("data/extracted_text"),
        help="Directory where extracted text (*.txt) will be stored.",
    )
    parser.add_argument(
        "--report-json",
        type=Path,
        default=Path("data/pdf_organize_report.json"),
        help="JSON report path.",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=Path("data/pdf_organize_summary.md"),
        help="Markdown summary path.",
    )
    parser.add_argument(
        "--disease-keywords-json",
        type=Path,
        default=None,
        help="Optional disease keyword mapping JSON. Format: {\"disease\": [\"keyword1\", ...]}",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    disease_keywords = load_disease_keywords(args.disease_keywords_json)

    records = organize_pdfs(
        pdf_dir=args.pdf_dir,
        output_dir=args.output_dir,
        text_dir=args.text_dir,
        disease_keywords=disease_keywords,
    )
    write_reports(records, args.report_json, args.summary_md)

    print(f"[done] processed={len(records)}")
    print(f"[done] json report={args.report_json}")
    print(f"[done] markdown summary={args.summary_md}")


if __name__ == "__main__":
    main()
