#!/usr/bin/env python3
"""Three-assistant pipeline for disease-focused paper retrieval and review.

Assistants:
1) RetrievalAssistant: fetches papers from OpenAlex for target diseases.
2) SummarizerAssistant: creates structured summaries from metadata + abstract.
3) RelevanceReviewerAssistant: evaluates topic fit against requested diseases.

The retrieval result includes a boolean flag indicating whether any author
has an affiliation that appears to be in Korea.

Quick start:
    python src/paper_assistants.py --max-results 20 --mailto kimtk7830@naver.com
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
DEFAULT_MAILTO = "kimtk7830@naver.com"
TARGET_DISEASES = {
    "type_2_diabetes": [
        "type 2 diabetes",
        "T2D",
        "insulin resistance",
    ],
    "alzheimers_disease": [
        "Alzheimer disease",
        "Alzheimer's disease",
        "AD dementia",
    ],
}


@dataclass
class PaperResult:
    paper_id: str
    title: str
    year: int | None
    doi: str | None
    pmid: str | None
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
    def __init__(
        self,
        *,
        mailto: str = DEFAULT_MAILTO,
        per_page: int = 25,
        timeout: int = 30,
        max_pages_per_term: int = 3,
        show_progress: bool = True,
    ) -> None:
        self.mailto = mailto
        self.per_page = per_page
        self.timeout = timeout
        self.max_pages_per_term = max_pages_per_term
        self.show_progress = show_progress

    def fetch(self, *, max_results: int = 20, skip_index: "ProcessedIndex | None" = None) -> list[PaperResult]:
        collected: list[PaperResult] = []
        seen_ids: set[str] = set()

        for disease_name, query_terms in TARGET_DISEASES.items():
            self._progress(f"[retrieve] disease={disease_name} terms={len(query_terms)}")
            for term in query_terms:
                if len(collected) >= max_results:
                    break
                self._progress(f"[retrieve] term={term} start")
                for page in range(1, self.max_pages_per_term + 1):
                    self._progress(f"[retrieve] term={term} page={page}")
                    works = self._query_openalex(term, page=page)
                    if not works:
                        self._progress(f"[retrieve] term={term} page={page} no_results")
                        break
                    for work in works:
                        parsed = self._parse_work(work)
                        if parsed.paper_id in seen_ids:
                            continue
                        if skip_index and skip_index.contains(parsed):
                            self._progress(
                                f"[skip] already_processed openalex={parsed.openalex_id} doi={parsed.doi} pmid={parsed.pmid}"
                            )
                            continue
                        parsed.disease_hits = self._disease_hits(parsed)
                        if not parsed.disease_hits:
                            continue
                        seen_ids.add(parsed.paper_id)
                        collected.append(parsed)
                        self._progress(
                            f"[collect] count={len(collected)}/{max_results} id={parsed.paper_id} diseases={parsed.disease_hits}"
                        )
                        if len(collected) >= max_results:
                            break
                    if len(collected) >= max_results:
                        break
            if len(collected) >= max_results:
                break

        return collected

    def _query_openalex(self, search_term: str, *, page: int) -> list[dict[str, Any]]:
        params = {
            "search": search_term,
            "per-page": self.per_page,
            "page": page,
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
        pmid = self._extract_pmid(work)
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
            pmid=pmid,
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
    def _extract_pmid(work: dict[str, Any]) -> str | None:
        ids = work.get("ids") or {}
        pmid = ids.get("pmid")
        if not pmid:
            return None
        pmid_text = str(pmid).strip().rstrip("/")
        if "/" in pmid_text:
            pmid_text = pmid_text.split("/")[-1]
        return pmid_text or None

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

    def _progress(self, message: str) -> None:
        if self.show_progress:
            print(message)


@dataclass
class ProcessedIndex:
    openalex_ids: set[str]
    dois: set[str]
    pmids: set[str]

    @classmethod
    def empty(cls) -> "ProcessedIndex":
        return cls(openalex_ids=set(), dois=set(), pmids=set())

    def contains(self, paper: PaperResult) -> bool:
        if paper.openalex_id and paper.openalex_id in self.openalex_ids:
            return True
        if paper.doi and paper.doi in self.dois:
            return True
        if paper.pmid and paper.pmid in self.pmids:
            return True
        return False

    def add(self, paper: PaperResult) -> None:
        if paper.openalex_id:
            self.openalex_ids.add(paper.openalex_id)
        if paper.doi:
            self.dois.add(paper.doi)
        if paper.pmid:
            self.pmids.add(paper.pmid)

    def extend(self, papers: list[PaperResult]) -> None:
        for paper in papers:
            self.add(paper)

    def to_json(self) -> dict[str, list[str]]:
        return {
            "openalex_ids": sorted(self.openalex_ids),
            "dois": sorted(self.dois),
            "pmids": sorted(self.pmids),
        }

    @classmethod
    def from_json(cls, payload: dict[str, Any]) -> "ProcessedIndex":
        return cls(
            openalex_ids=set(str(v) for v in payload.get("openalex_ids", [])),
            dois=set(str(v) for v in payload.get("dois", [])),
            pmids=set(str(v) for v in payload.get("pmids", [])),
        )


def _load_processed_index(index_path: Path, *, show_progress: bool) -> ProcessedIndex:
    if not index_path.exists():
        if show_progress:
            print(f"[state] no existing index file: {index_path}")
        return ProcessedIndex.empty()
    try:
        payload = json.loads(index_path.read_text(encoding="utf-8"))
        index = ProcessedIndex.from_json(payload)
        if show_progress:
            print(
                f"[state] loaded index: openalex={len(index.openalex_ids)} doi={len(index.dois)} pmid={len(index.pmids)} from {index_path}"
            )
        return index
    except json.JSONDecodeError:
        if show_progress:
            print(f"[state] invalid JSON index, starting fresh: {index_path}")
        return ProcessedIndex.empty()


def _save_processed_index(index_path: Path, index: ProcessedIndex, *, show_progress: bool) -> None:
    index_path.parent.mkdir(parents=True, exist_ok=True)
    index_path.write_text(json.dumps(index.to_json(), ensure_ascii=False, indent=2), encoding="utf-8")
    if show_progress:
        print(
            f"[state] saved index: openalex={len(index.openalex_ids)} doi={len(index.dois)} pmid={len(index.pmids)} to {index_path}"
        )


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


def run_pipeline(
    *,
    max_results: int,
    output_path: Path,
    index_path: Path,
    mailto: str,
    max_pages_per_term: int,
    show_progress: bool,
) -> dict[str, Any]:
    retriever = RetrievalAssistant(mailto=mailto, max_pages_per_term=max_pages_per_term, show_progress=show_progress)
    summarizer = SummarizerAssistant()
    reviewer = RelevanceReviewerAssistant()

    processed_index = _load_processed_index(index_path, show_progress=show_progress)
    papers = retriever.fetch(max_results=max_results, skip_index=processed_index)
    summaries = summarizer.summarize(papers)
    reviews = reviewer.review(papers)
    report = build_report(papers, summaries, reviews)
    processed_index.extend(papers)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    _save_processed_index(index_path, processed_index, show_progress=show_progress)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Disease paper retrieval + summary + relevance review assistants")
    parser.add_argument("--max-results", type=int, default=12, help="Maximum number of papers to store")
    parser.add_argument(
        "--max-pages-per-term",
        type=int,
        default=3,
        help="Maximum OpenAlex pages to fetch for each query term",
    )
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
    parser.add_argument(
        "--index-path",
        type=Path,
        default=Path("outputs/paper_assistant_seen_ids.json"),
        help="Path to persistent processed-id index used for skip behavior",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Disable pipeline progress logs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = run_pipeline(
        max_results=args.max_results,
        output_path=args.output,
        index_path=args.index_path,
        mailto=args.mailto,
        max_pages_per_term=args.max_pages_per_term,
        show_progress=not args.quiet,
    )
    print(json.dumps({"saved": str(args.output), "count": report["count"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
