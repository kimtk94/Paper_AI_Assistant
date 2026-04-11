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
import re
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
import xml.etree.ElementTree as ET
from typing import Any

OPENALEX_WORKS_URL = "https://api.openalex.org/works"
IDCONV_URL = "https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/"
OA_FCGI_URL = "https://pmc.ncbi.nlm.nih.gov/utils/oa/oa.fcgi"
UNPAYWALL_URL = "https://api.unpaywall.org/v2"
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


@dataclass
class MethodDataReview:
    paper_id: str
    method: str
    data_db: str
    dataset: str
    checklist: list[str]


@dataclass
class DownloadResult:
    input_id: str
    normalized_id: str
    doi: str | None
    pmid: str | None
    pmcid: str | None
    downloaded: bool
    filename: str | None
    source_url: str | None
    reason: str


@dataclass
class MissingPdfCandidate:
    paper_id: str
    title: str
    doi: str | None
    pmid: str | None
    identifier: str
    miss_reason: str


@dataclass
class ScihubDownloadAttempt:
    command: list[str]
    returncode: int
    input_file: str
    candidate_count: int
    stdout_tail: str
    stderr_tail: str


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


class MethodDataReviewerAssistant:
    def review(self, papers: list[PaperResult]) -> list[MethodDataReview]:
        reviews: list[MethodDataReview] = []
        for paper in papers:
            joined = " ".join([paper.title, paper.abstract or "", " ".join(paper.concepts)]).lower()
            method = self._infer_method(joined)
            data_db = self._infer_data_db(paper)
            dataset = self._infer_dataset(joined)
            checklist = [
                "train/val/test split 명시 여부 확인 필요",
                "external validation 여부 확인 필요",
                "data leakage 위험(환자/문서 단위 split) 확인 필요",
            ]
            reviews.append(
                MethodDataReview(
                    paper_id=paper.paper_id,
                    method=method,
                    data_db=data_db,
                    dataset=dataset,
                    checklist=checklist,
                )
            )
        return reviews

    @staticmethod
    def _infer_method(joined_text: str) -> str:
        if "mendelian" in joined_text:
            return "Mendelian randomization"
        if "gwas" in joined_text or "genome" in joined_text:
            return "GWAS / 유전연관 분석"
        if "single-cell" in joined_text:
            return "Single-cell omics 분석"
        if "meta-analysis" in joined_text or "systematic review" in joined_text:
            return "Systematic review / meta-analysis"
        if "deep learning" in joined_text or "neural network" in joined_text:
            return "딥러닝 기반 예측 모델"
        return "원문 Methods 섹션 확인 필요"

    @staticmethod
    def _infer_data_db(paper: PaperResult) -> str:
        if paper.source:
            return f"{paper.source} (retrieval source)"
        return "데이터 출처 확인 필요"

    @staticmethod
    def _infer_dataset(joined_text: str) -> str:
        if "uk biobank" in joined_text:
            return "UK Biobank"
        if "adni" in joined_text:
            return "ADNI"
        if "tcga" in joined_text:
            return "TCGA"
        if "mimic" in joined_text:
            return "MIMIC"
        return "논문 원문에서 dataset 이름/버전 확인 필요"


class PdfDownloader:
    def __init__(self, *, outdir: Path, timeout: int = 30, unpaywall_email: str | None = None) -> None:
        self.outdir = outdir
        self.timeout = timeout
        self.unpaywall_email = unpaywall_email
        self.outdir.mkdir(parents=True, exist_ok=True)

    def download(self, raw_identifier: str) -> DownloadResult:
        normalized = normalize_identifier(raw_identifier)
        mapping = self._id_convert(normalized)
        doi = mapping.get("doi")
        pmid = mapping.get("pmid")
        pmcid = mapping.get("pmcid")

        if pmcid:
            pmc_pdf_url = self._resolve_pmc_pdf_url(pmcid)
            if pmc_pdf_url:
                filename = self._safe_filename(normalized, pmcid=pmcid, doi=doi) + ".pdf"
                dst = self.outdir / filename
                if self._download_file(pmc_pdf_url, dst):
                    return DownloadResult(
                        input_id=raw_identifier,
                        normalized_id=normalized,
                        doi=doi,
                        pmid=pmid,
                        pmcid=pmcid,
                        downloaded=True,
                        filename=filename,
                        source_url=pmc_pdf_url,
                        reason="Downloaded from PubMed Central OA endpoint.",
                    )

        if doi and self.unpaywall_email:
            pdf_url = self._resolve_unpaywall_pdf_url(doi)
            if pdf_url:
                filename = self._safe_filename(normalized, pmcid=pmcid, doi=doi) + ".pdf"
                dst = self.outdir / filename
                if self._download_file(pdf_url, dst):
                    return DownloadResult(
                        input_id=raw_identifier,
                        normalized_id=normalized,
                        doi=doi,
                        pmid=pmid,
                        pmcid=pmcid,
                        downloaded=True,
                        filename=filename,
                        source_url=pdf_url,
                        reason="Downloaded from Unpaywall open-access link.",
                    )

        if doi and not self.unpaywall_email:
            reason = "No PMC OA PDF found. Provide --unpaywall-email to try legal OA mirrors by DOI."
        else:
            reason = "No open-access PDF found from configured sources."

        return DownloadResult(
            input_id=raw_identifier,
            normalized_id=normalized,
            doi=doi,
            pmid=pmid,
            pmcid=pmcid,
            downloaded=False,
            filename=None,
            source_url=None,
            reason=reason,
        )

    def _id_convert(self, identifier: str) -> dict[str, str | None]:
        qid = quote(identifier, safe="")
        url = f"{IDCONV_URL}?ids={qid}&format=json"
        try:
            payload = http_get_json(url, timeout=self.timeout)
        except (URLError, HTTPError, json.JSONDecodeError):
            return {"doi": extract_doi(identifier), "pmid": extract_pmid(identifier), "pmcid": extract_pmcid(identifier)}

        records = payload.get("records") or []
        if not records:
            return {"doi": extract_doi(identifier), "pmid": extract_pmid(identifier), "pmcid": extract_pmcid(identifier)}

        first = records[0]
        return {
            "doi": text_or_none(first.get("doi")) or extract_doi(identifier),
            "pmid": digits_only(text_or_none(first.get("pmid"))) or extract_pmid(identifier),
            "pmcid": normalize_pmcid(text_or_none(first.get("pmcid")) or extract_pmcid(identifier)),
        }

    def _resolve_pmc_pdf_url(self, pmcid: str) -> str | None:
        url = f"{OA_FCGI_URL}?id={quote(pmcid, safe='')}"
        try:
            xml_text = http_get_text(url, timeout=self.timeout)
        except (URLError, HTTPError):
            return None

        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return None

        for link in root.findall(".//link"):
            fmt = (link.attrib.get("format") or "").lower()
            href = link.attrib.get("href")
            if fmt == "pdf" and href:
                return href
        return None

    def _resolve_unpaywall_pdf_url(self, doi: str) -> str | None:
        email = quote(self.unpaywall_email or "", safe="")
        doi_q = quote(doi, safe="")
        url = f"{UNPAYWALL_URL}/{doi_q}?email={email}"
        try:
            payload = http_get_json(url, timeout=self.timeout)
        except (URLError, HTTPError, json.JSONDecodeError):
            return None
        best = payload.get("best_oa_location") or {}
        pdf = text_or_none(best.get("url_for_pdf"))
        if pdf:
            return pdf
        return text_or_none(best.get("url"))

    def _download_file(self, url: str, dst: Path) -> bool:
        req = Request(url, headers={"User-Agent": "PaperAI-Assistant/1.0"})
        try:
            with urlopen(req, timeout=self.timeout) as response:
                data = response.read()
        except (URLError, HTTPError):
            return False
        if not data:
            return False
        dst.write_bytes(data)
        return True

    @staticmethod
    def _safe_filename(identifier: str, *, pmcid: str | None, doi: str | None) -> str:
        base = pmcid or doi or identifier
        base = base.replace("/", "_")
        base = re.sub(r"[^A-Za-z0-9._-]", "_", base)
        base = re.sub(r"_+", "_", base).strip("_")
        return base[:160] or "paper"


def build_report(
    papers: list[PaperResult],
    summaries: list[PaperSummary],
    reviews: list[RelevanceReview],
    method_data_reviews: list[MethodDataReview],
) -> dict[str, Any]:
    summary_map = {s.paper_id: s for s in summaries}
    review_map = {r.paper_id: r for r in reviews}
    method_data_map = {m.paper_id: m for m in method_data_reviews}

    rows: list[dict[str, Any]] = []
    for paper in papers:
        row = asdict(paper)
        row["summary"] = asdict(summary_map[paper.paper_id])
        row["review"] = asdict(review_map[paper.paper_id])
        row["method_data_review"] = asdict(method_data_map[paper.paper_id])
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
    method_data_reviewer = MethodDataReviewerAssistant()

    processed_index = _load_processed_index(index_path, show_progress=show_progress)
    papers = retriever.fetch(max_results=max_results, skip_index=processed_index)
    summaries = summarizer.summarize(papers)
    reviews = reviewer.review(papers)
    method_data_reviews = method_data_reviewer.review(papers)
    report = build_report(papers, summaries, reviews, method_data_reviews)
    processed_index.extend(papers)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    _save_processed_index(index_path, processed_index, show_progress=show_progress)
    return report


def _load_papers_from_report(report_path: Path) -> list[PaperResult]:
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    rows = payload.get("results", [])
    papers: list[PaperResult] = []
    for row in rows:
        papers.append(
            PaperResult(
                paper_id=str(row.get("paper_id", "")),
                title=str(row.get("title", "(untitled)")),
                year=row.get("year"),
                doi=row.get("doi"),
                pmid=row.get("pmid"),
                openalex_id=str(row.get("openalex_id", row.get("paper_id", ""))),
                source=str(row.get("source", "unknown")),
                authors=[str(a) for a in row.get("authors", [])],
                abstract=row.get("abstract"),
                concepts=[str(c) for c in row.get("concepts", [])],
                primary_location_url=row.get("primary_location_url"),
                korea_affiliation_present=bool(row.get("korea_affiliation_present", False)),
                korea_affiliation_evidence=[str(v) for v in row.get("korea_affiliation_evidence", [])],
                disease_hits=[str(v) for v in row.get("disease_hits", [])],
            )
        )
    return papers


def _write_review_markdown(report: dict[str, Any], md_path: Path) -> None:
    lines = [
        "# Paper Review (Method / Data DB / Dataset)",
        "",
        f"- Generated at (UTC): {report.get('generated_at_utc', '')}",
        f"- Paper count: {report.get('count', 0)}",
        "",
        "| Title | Method | Data DB | Dataset | Relevance |",
        "|---|---|---|---|---|",
    ]
    for row in report.get("results", []):
        review = row.get("method_data_review", {})
        rel = row.get("review", {})
        title = str(row.get("title", "")).replace("|", "/")
        method = str(review.get("method", "")).replace("|", "/")
        data_db = str(review.get("data_db", "")).replace("|", "/")
        dataset = str(review.get("dataset", "")).replace("|", "/")
        relevance = "relevant" if rel.get("is_relevant") else "check-needed"
        lines.append(f"| {title} | {method} | {data_db} | {dataset} | {relevance} |")
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def http_get_json(url: str, *, timeout: int) -> dict[str, Any]:
    req = Request(url, headers={"User-Agent": "PaperAI-Assistant/1.0"})
    with urlopen(req, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def http_get_text(url: str, *, timeout: int) -> str:
    req = Request(url, headers={"User-Agent": "PaperAI-Assistant/1.0"})
    with urlopen(req, timeout=timeout) as response:
        return response.read().decode("utf-8")


def text_or_none(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def digits_only(text: str | None) -> str | None:
    if not text:
        return None
    digits = "".join(ch for ch in text if ch.isdigit())
    return digits or None


def normalize_pmcid(text: str | None) -> str | None:
    if not text:
        return None
    cleaned = text.upper().strip()
    if not cleaned.startswith("PMC"):
        if cleaned.isdigit():
            cleaned = f"PMC{cleaned}"
        else:
            return None
    return cleaned


def extract_doi(text: str) -> str | None:
    match = re.search(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+", text)
    return match.group(0) if match else None


def extract_pmid(text: str) -> str | None:
    if "pmid" in text.lower():
        nums = re.findall(r"\d+", text)
        return nums[0] if nums else None
    if text.isdigit() and not text.upper().startswith("PMC"):
        return text
    return None


def extract_pmcid(text: str) -> str | None:
    match = re.search(r"PMC\d+", text.upper())
    if match:
        return match.group(0)
    return None


def normalize_identifier(text: str) -> str:
    raw = text.strip()
    if not raw:
        return raw
    doi = extract_doi(raw)
    if doi:
        return doi
    pmcid = extract_pmcid(raw)
    if pmcid:
        return pmcid
    pmid = extract_pmid(raw)
    if pmid:
        return pmid
    return raw


def _download_pdfs_for_papers(
    papers: list[PaperResult],
    *,
    outdir: Path,
    report_path: Path,
    timeout: int,
    unpaywall_email: str | None,
) -> dict[str, Any]:
    downloader = PdfDownloader(outdir=outdir, timeout=timeout, unpaywall_email=unpaywall_email)
    raw_identifiers: list[str] = []
    for paper in papers:
        if paper.pmid:
            raw_identifiers.append(paper.pmid)
            continue
        if paper.doi:
            raw_identifiers.append(paper.doi)

    seen: set[str] = set()
    identifiers: list[str] = []
    for identifier in raw_identifiers:
        key = identifier.strip()
        if key and key not in seen:
            seen.add(key)
            identifiers.append(key)

    results: list[DownloadResult] = []
    for identifier in identifiers:
        result = downloader.download(identifier)
        results.append(result)
        status = "OK" if result.downloaded else "MISS"
        print(f"[pdf:{status}] {identifier} -> {result.reason}")

    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "total": len(results),
        "downloaded": sum(1 for r in results if r.downloaded),
        "missed": sum(1 for r in results if not r.downloaded),
        "results": [asdict(r) for r in results],
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"pdf_report_saved": str(report_path), "downloaded": payload["downloaded"]}, ensure_ascii=False))
    return payload


def _collect_missing_pdf_candidates(
    papers: list[PaperResult],
    *,
    pdf_payload: dict[str, Any],
) -> list[MissingPdfCandidate]:
    doi_map = {p.doi: p for p in papers if p.doi}
    pmid_map = {p.pmid: p for p in papers if p.pmid}
    candidates: list[MissingPdfCandidate] = []

    for row in pdf_payload.get("results", []):
        if row.get("downloaded"):
            continue
        doi = text_or_none(row.get("doi"))
        pmid = text_or_none(row.get("pmid"))
        miss_reason = str(row.get("reason") or "")

        paper = None
        if doi and doi in doi_map:
            paper = doi_map[doi]
        elif pmid and pmid in pmid_map:
            paper = pmid_map[pmid]
        if paper is None:
            continue

        identifier = doi or pmid or str(row.get("normalized_id") or "")
        if not identifier:
            continue
        candidates.append(
            MissingPdfCandidate(
                paper_id=paper.paper_id,
                title=paper.title,
                doi=paper.doi,
                pmid=paper.pmid,
                identifier=identifier,
                miss_reason=miss_reason,
            )
        )

    return candidates


def _run_scihub_cli(
    candidates: list[MissingPdfCandidate],
    *,
    scihub_bin: str,
    outdir: Path,
    input_file: Path,
    email: str | None,
    extra_args: list[str],
) -> dict[str, Any]:
    seen: set[str] = set()
    ids: list[str] = []
    for candidate in candidates:
        key = candidate.identifier.strip()
        if key and key not in seen:
            seen.add(key)
            ids.append(key)

    input_file.parent.mkdir(parents=True, exist_ok=True)
    input_file.write_text("\n".join(ids) + "\n", encoding="utf-8")

    command = [scihub_bin, str(input_file), "-o", str(outdir)]
    if email:
        command.extend(["--email", email])
    command.extend(extra_args)

    proc = subprocess.run(command, capture_output=True, text=True)
    stdout_tail = proc.stdout[-1000:]
    stderr_tail = proc.stderr[-1000:]
    status = "OK" if proc.returncode == 0 else "FAIL"
    print(f"[scihub:{status}] input={input_file} ids={len(ids)} outdir={outdir}")

    attempt = ScihubDownloadAttempt(
        command=command,
        returncode=proc.returncode,
        input_file=str(input_file),
        candidate_count=len(ids),
        stdout_tail=stdout_tail,
        stderr_tail=stderr_tail,
    )

    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "total_runs": 1,
        "success_runs": 1 if proc.returncode == 0 else 0,
        "failed_runs": 0 if proc.returncode == 0 else 1,
        "results": [asdict(attempt)],
    }


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
    parser.add_argument(
        "--from-report",
        type=Path,
        default=None,
        help="Re-process an existing JSON report's results list instead of fetching new papers",
    )
    parser.add_argument(
        "--review-md-output",
        type=Path,
        default=None,
        help="Optional markdown output path for method/data/dataset review table",
    )
    parser.add_argument(
        "--download-pdfs",
        action="store_true",
        help="Download legal OA PDFs for collected papers using PMID/DOI.",
    )
    parser.add_argument(
        "--pdf-outdir",
        type=Path,
        default=Path("outputs/pdfs"),
        help="Directory where downloaded PDFs are stored.",
    )
    parser.add_argument(
        "--pdf-report",
        type=Path,
        default=Path("outputs/pdf_download_report.json"),
        help="JSON report path for PDF download results.",
    )
    parser.add_argument(
        "--unpaywall-email",
        default=None,
        help="Email for Unpaywall API (optional DOI OA fallback).",
    )
    parser.add_argument(
        "--pdf-timeout",
        type=int,
        default=30,
        help="HTTP timeout seconds for PDF download requests.",
    )
    parser.add_argument(
        "--missing-pdf-report",
        type=Path,
        default=Path("outputs/pdf_missing_candidates.json"),
        help="JSON report path for papers that still miss PDF after legal OA download attempts.",
    )
    parser.add_argument(
        "--retry-missing-with-scihub-cli",
        action="store_true",
        help="Retry missing-PDF candidates via scihub-cli batch mode.",
    )
    parser.add_argument(
        "--scihub-cli-bin",
        default="scihub-cli",
        help="scihub-cli executable name or full path.",
    )
    parser.add_argument(
        "--scihub-input-file",
        type=Path,
        default=Path("outputs/scihub_missing_input.txt"),
        help="Input text file path generated for scihub-cli (one identifier per line).",
    )
    parser.add_argument(
        "--scihub-email",
        default=None,
        help="Optional email passed to scihub-cli --email (for Unpaywall integration).",
    )
    parser.add_argument(
        "--scihub-extra-args",
        nargs="*",
        default=[],
        help="Additional arguments passed through to scihub-cli (e.g. --verbose --no-fast-fail).",
    )
    parser.add_argument(
        "--scihub-report",
        type=Path,
        default=Path("outputs/scihub_download_report.json"),
        help="JSON report path for scihub command execution results.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.from_report:
        papers = _load_papers_from_report(args.from_report)
        summarizer = SummarizerAssistant()
        reviewer = RelevanceReviewerAssistant()
        method_data_reviewer = MethodDataReviewerAssistant()
        report = build_report(
            papers,
            summarizer.summarize(papers),
            reviewer.review(papers),
            method_data_reviewer.review(papers),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    else:
        report = run_pipeline(
            max_results=args.max_results,
            output_path=args.output,
            index_path=args.index_path,
            mailto=args.mailto,
            max_pages_per_term=args.max_pages_per_term,
            show_progress=not args.quiet,
        )
        papers = _load_papers_from_report(args.output)
    if args.review_md_output:
        _write_review_markdown(report, args.review_md_output)
    if args.download_pdfs:
        pdf_payload = _download_pdfs_for_papers(
            papers,
            outdir=args.pdf_outdir,
            report_path=args.pdf_report,
            timeout=args.pdf_timeout,
            unpaywall_email=args.unpaywall_email,
        )
        missing_candidates = _collect_missing_pdf_candidates(papers, pdf_payload=pdf_payload)
        args.missing_pdf_report.parent.mkdir(parents=True, exist_ok=True)
        args.missing_pdf_report.write_text(
            json.dumps(
                {
                    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                    "total": len(missing_candidates),
                    "results": [asdict(c) for c in missing_candidates],
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )
        print(json.dumps({"missing_pdf_report_saved": str(args.missing_pdf_report), "total": len(missing_candidates)}))

        if args.retry_missing_with_scihub_cli and missing_candidates:
            scihub_payload = _run_scihub_cli(
                missing_candidates,
                scihub_bin=args.scihub_cli_bin,
                outdir=args.pdf_outdir,
                input_file=args.scihub_input_file,
                email=args.scihub_email,
                extra_args=args.scihub_extra_args,
            )
            args.scihub_report.parent.mkdir(parents=True, exist_ok=True)
            args.scihub_report.write_text(
                json.dumps(scihub_payload, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            print(
                json.dumps(
                    {"scihub_report_saved": str(args.scihub_report), "success_runs": scihub_payload["success_runs"]}
                )
            )
    print(json.dumps({"saved": str(args.output), "count": report["count"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
