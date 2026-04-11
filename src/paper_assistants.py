#!/usr/bin/env python3
"""Three-assistant pipeline for disease-focused paper retrieval and review.

Assistants:
1) RetrievalAssistant: fetches papers from OpenAlex for target diseases.
2) SummarizerAssistant: creates structured summaries from metadata + abstract.
3) RelevanceReviewerAssistant: evaluates topic fit against requested diseases.

The retrieval result includes a boolean flag indicating whether any author
has an affiliation that appears to be in Korea.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from urllib.error import URLError, HTTPError
from typing import Any

OPENALEX_WORKS_URL = "https://api.openalex.org/works"
DEFAULT_MAILTO = "paper-assistant@example.com"
TARGET_DISEASES = {
    "type_2_diabetes": [
        '"type 2 diabetes"',
        '"T2D"',
        '"insulin resistance"',
    ],
    "alzheimers_disease": [
        '"Alzheimer disease"',
        '"Alzheimer\'s disease"',
        '"AD dementia"',
    ],
}


@dataclass
class PaperResult:
    paper_id: str
    title: str
    year: int | None
    doi: str | None
    openalex_id: str
    source: str
    authors: list[str]
    abstract: str | None
    concepts: list[str]
    primary_location_url: str | None
    korea_affiliation_present: bool
    korea_affiliation_evidence: list[str]
    disease_hits: list[str]


@dataclass
class PaperSummary:
    paper_id: str
    summary: str
    methods_hint: str
    possible_limitations: str


@dataclass
class RelevanceReview:
    paper_id: str
    is_relevant: bool
    matched_diseases: list[str]
    reason: str


class RetrievalAssistant:
    def __init__(self, *, mailto: str = DEFAULT_MAILTO, per_page: int = 25, timeout: int = 30) -> None:
        self.mailto = mailto
        self.per_page = per_page
        self.timeout = timeout

    def fetch(self, *, max_results: int = 20) -> list[PaperResult]:
        collected: list[PaperResult] = []
        seen_ids: set[str] = set()

        for disease_name, query_terms in TARGET_DISEASES.items():
            for term in query_terms:
                if len(collected) >= max_results:
                    break
                works = self._query_openalex(term)
                for work in works:
                    parsed = self._parse_work(work)
                    if parsed.paper_id in seen_ids:
                        continue
                    parsed.disease_hits = self._disease_hits(parsed)
                    if not parsed.disease_hits:
                        continue
                    seen_ids.add(parsed.paper_id)
                    collected.append(parsed)
                    if len(collected) >= max_results:
                        break
            if len(collected) >= max_results:
                break

        return collected

    def _query_openalex(self, search_term: str) -> list[dict[str, Any]]:
        params = {
            "search": search_term,
            "per-page": self.per_page,
            "sort": "cited_by_count:desc",
            "mailto": self.mailto,
        }
        query = urlencode(params)
        url = f"{OPENALEX_WORKS_URL}?{query}"
        try:
            with urlopen(url, timeout=self.timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
            return payload.get("results", [])
        except (URLError, HTTPError):
            return []

    def _parse_work(self, work: dict[str, Any]) -> PaperResult:
        paper_id = str(work.get("id", ""))
        title = str(work.get("title") or "(untitled)")
        year = work.get("publication_year")
        doi = work.get("doi")
        source = "openalex"
        authors = []
        korea_evidence: list[str] = []

        authorships = work.get("authorships") or []
        for authorship in authorships:
            author_obj = authorship.get("author") or {}
            name = author_obj.get("display_name")
            if name:
                authors.append(str(name))

            raw_affiliations = authorship.get("raw_affiliation_strings") or []
            for aff in raw_affiliations:
                if self._is_korea_affiliation(str(aff)):
                    korea_evidence.append(str(aff))

            institutions = authorship.get("institutions") or []
            for inst in institutions:
                country_code = str(inst.get("country_code") or "").upper()
                display_name = str(inst.get("display_name") or "")
                if country_code == "KR":
                    korea_evidence.append(f"{display_name} (country_code=KR)")

        abstract = self._reconstruct_abstract(work.get("abstract_inverted_index"))
        concepts = [str(c.get("display_name")) for c in (work.get("concepts") or []) if c.get("display_name")]
        primary_location = work.get("primary_location") or {}
        primary_location_url = primary_location.get("landing_page_url") or primary_location.get("pdf_url")

        return PaperResult(
            paper_id=paper_id,
            title=title,
            year=year,
            doi=doi,
            openalex_id=paper_id,
            source=source,
            authors=authors,
            abstract=abstract,
            concepts=concepts,
            primary_location_url=primary_location_url,
            korea_affiliation_present=bool(korea_evidence),
            korea_affiliation_evidence=sorted(set(korea_evidence)),
            disease_hits=[],
        )

    @staticmethod
    def _reconstruct_abstract(inverted_index: dict[str, list[int]] | None) -> str | None:
        if not inverted_index:
            return None
        positions: dict[int, str] = {}
        for token, pos_list in inverted_index.items():
            for pos in pos_list:
                positions[pos] = token
        if not positions:
            return None
        return " ".join(positions[i] for i in sorted(positions))

    @staticmethod
    def _is_korea_affiliation(affiliation_text: str) -> bool:
        text = affiliation_text.lower()
        keywords = ["korea", "republic of korea", "south korea", "korean"]
        return any(key in text for key in keywords)

    @staticmethod
    def _disease_hits(paper: PaperResult) -> list[str]:
        haystack = " ".join(
            [
                paper.title,
                paper.abstract or "",
                " ".join(paper.concepts),
            ]
        ).lower()

        hits: list[str] = []
        t2d_terms = ["type 2 diabetes", "t2d", "insulin resistance"]
        ad_terms = ["alzheimer", "ad dementia", "amyloid"]
        if any(term in haystack for term in t2d_terms):
            hits.append("type_2_diabetes")
        if any(term in haystack for term in ad_terms):
            hits.append("alzheimers_disease")
        return hits


class SummarizerAssistant:
    def summarize(self, papers: list[PaperResult]) -> list[PaperSummary]:
        summaries: list[PaperSummary] = []
        for paper in papers:
            abstract = paper.abstract or "Abstract unavailable."
            snippet = abstract[:380].strip()
            summary = f"{paper.title} ({paper.year})는 {', '.join(paper.disease_hits)} 관련 연구로 보이며, 핵심 초록은 다음과 같습니다: {snippet}"
            methods_hint = self._infer_methods(paper)
            limitations = "초록 기반 자동 요약이므로 실험 설정/통계 검정의 세부 한계는 원문 확인이 필요합니다."
            summaries.append(
                PaperSummary(
                    paper_id=paper.paper_id,
                    summary=summary,
                    methods_hint=methods_hint,
                    possible_limitations=limitations,
                )
            )
        return summaries

    @staticmethod
    def _infer_methods(paper: PaperResult) -> str:
        joined = " ".join([paper.title, paper.abstract or "", " ".join(paper.concepts)]).lower()
        if "mendelian" in joined:
            return "Mendelian randomization 가능성이 높음"
        if "genome" in joined or "gwas" in joined:
            return "GWAS/유전연관 분석 계열 가능성이 높음"
        if "single-cell" in joined:
            return "single-cell 전사체 분석 가능성이 높음"
        return "방법론은 원문 Methods 섹션 확인 필요"


class RelevanceReviewerAssistant:
    def review(self, papers: list[PaperResult]) -> list[RelevanceReview]:
        reviews: list[RelevanceReview] = []
        for paper in papers:
            matched = paper.disease_hits
            relevant = bool(matched)
            reason = (
                f"질병 키워드 매칭: {', '.join(matched)}"
                if matched
                else "type 2 diabetes / Alzheimer 관련 키워드 매칭이 약함"
            )
            reviews.append(
                RelevanceReview(
                    paper_id=paper.paper_id,
                    is_relevant=relevant,
                    matched_diseases=matched,
                    reason=reason,
                )
            )
        return reviews


def build_report(
    papers: list[PaperResult],
    summaries: list[PaperSummary],
    reviews: list[RelevanceReview],
) -> dict[str, Any]:
    summary_map = {s.paper_id: s for s in summaries}
    review_map = {r.paper_id: r for r in reviews}

    rows: list[dict[str, Any]] = []
    for paper in papers:
        row = asdict(paper)
        row["summary"] = asdict(summary_map[paper.paper_id])
        row["review"] = asdict(review_map[paper.paper_id])
        rows.append(row)

    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "query_targets": ["type_2_diabetes", "alzheimers_disease"],
        "count": len(rows),
        "results": rows,
    }


def run_pipeline(*, max_results: int, output_path: Path, mailto: str) -> dict[str, Any]:
    retriever = RetrievalAssistant(mailto=mailto)
    summarizer = SummarizerAssistant()
    reviewer = RelevanceReviewerAssistant()

    papers = retriever.fetch(max_results=max_results)
    summaries = summarizer.summarize(papers)
    reviews = reviewer.review(papers)
    report = build_report(papers, summaries, reviews)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Disease paper retrieval + summary + relevance review assistants")
    parser.add_argument("--max-results", type=int, default=12, help="Maximum number of papers to store")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/paper_assistant_report.json"),
        help="Path to output JSON report",
    )
    parser.add_argument(
        "--mailto",
        default=DEFAULT_MAILTO,
        help="Contact email sent to OpenAlex API as polite pool identifier",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = run_pipeline(max_results=args.max_results, output_path=args.output, mailto=args.mailto)
    print(json.dumps({"saved": str(args.output), "count": report["count"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
