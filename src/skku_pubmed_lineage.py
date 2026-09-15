#!/usr/bin/env python3
"""Build an SKKU-centered PubMed publication lineage graph.

The pipeline has two layers:
1) Seed papers: PubMed papers with Sungkyunkwan/SKKU in author affiliation.
2) Connections: citation links, shared SKKU authors, and topic similarity.

Important:
- PubMed affiliation proves publication-time affiliation, not alumni status.
- Author-name matching can be ambiguous. ORCID is retained when PubMed provides it.
"""

from __future__ import annotations

import argparse
import calendar
import csv
import hashlib
import http.client
import json
import math
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "skku_pubmed_lineage"
SKKU_RE = re.compile(r"\bsungkyunkwan\b|\bskku\b", re.I)


@dataclass
class Author:
    name: str
    last_name: str
    fore_name: str
    initials: str
    orcid: str
    affiliations: list[str]
    is_skku: bool


@dataclass
class Paper:
    pmid: str
    title: str
    journal: str
    publication_date: str
    year: int
    doi: str
    pmcid: str
    pubmed_url: str
    authors: list[Author]
    skku_authors: list[str]
    skku_orcids: list[str]
    skku_affiliation_evidence: list[str]
    abstract: str
    mesh_terms: list[str]
    keywords: list[str]
    publication_types: list[str]


@dataclass
class Edge:
    source: str
    target: str
    relation: str
    weight: float
    evidence: str


def clean_text(node: Optional[ET.Element]) -> str:
    if node is None:
        return ""
    return " ".join("".join(node.itertext()).split())


def child_text(root: ET.Element, path: str) -> str:
    return clean_text(root.find(path))


def parse_year(article: ET.Element) -> tuple[str, int]:
    pub = article.find("./MedlineCitation/Article/Journal/JournalIssue/PubDate")
    if pub is None:
        return "", 0
    year_s = child_text(pub, "Year")
    month = child_text(pub, "Month")
    day = child_text(pub, "Day")
    medline = child_text(pub, "MedlineDate")
    if not year_s and medline:
        m = re.search(r"(19|20)\d{2}", medline)
        year_s = m.group(0) if m else ""
    try:
        year = int(year_s)
    except ValueError:
        year = 0
    date = " ".join(x for x in [year_s, month, day] if x) or medline
    return date, year


def parse_ids(article: ET.Element) -> tuple[str, str]:
    doi = ""
    pmcid = ""
    for node in article.findall("./PubmedData/ArticleIdList/ArticleId"):
        value = clean_text(node)
        id_type = (node.attrib.get("IdType") or "").lower()
        if id_type == "doi":
            doi = value
        elif id_type == "pmc":
            pmcid = value
    return doi, pmcid


def parse_author(node: ET.Element) -> Author:
    last = child_text(node, "LastName")
    fore = child_text(node, "ForeName")
    initials = child_text(node, "Initials")
    collective = child_text(node, "CollectiveName")
    name = " ".join(x for x in [fore, last] if x).strip() or collective or last or initials or "Unknown"

    affs = []
    for aff in node.findall("./AffiliationInfo/Affiliation"):
        value = clean_text(aff)
        if value:
            affs.append(value)

    orcid = ""
    for ident in node.findall("./Identifier"):
        if (ident.attrib.get("Source") or "").upper() == "ORCID":
            orcid = clean_text(ident).replace("https://orcid.org/", "")
            break

    return Author(
        name=name,
        last_name=last,
        fore_name=fore,
        initials=initials,
        orcid=orcid,
        affiliations=affs,
        is_skku=any(SKKU_RE.search(a) for a in affs),
    )


def parse_abstract(article: ET.Element) -> str:
    out = []
    for node in article.findall("./MedlineCitation/Article/Abstract/AbstractText"):
        value = clean_text(node)
        if not value:
            continue
        label = (node.attrib.get("Label") or "").strip()
        out.append(f"{label}: {value}" if label else value)
    return "\n".join(out)


def parse_paper(article: ET.Element) -> Paper:
    citation = article.find("./MedlineCitation")
    art = article.find("./MedlineCitation/Article")
    if citation is None or art is None:
        raise ValueError("Invalid PubMed XML record")

    pmid = child_text(citation, "PMID")
    date, year = parse_year(article)
    doi, pmcid = parse_ids(article)
    authors = [parse_author(x) for x in art.findall("./AuthorList/Author")]
    skku_authors = [a.name for a in authors if a.is_skku]
    skku_orcids = sorted({a.orcid for a in authors if a.is_skku and a.orcid})

    evidence = []
    for author in authors:
        if not author.is_skku:
            continue
        for affiliation in author.affiliations:
            if SKKU_RE.search(affiliation):
                evidence.append(f"{author.name}: {affiliation}")

    mesh = [
        clean_text(x)
        for x in citation.findall("./MeshHeadingList/MeshHeading/DescriptorName")
        if clean_text(x)
    ]
    keywords = [
        clean_text(x)
        for x in citation.findall("./KeywordList/Keyword")
        if clean_text(x)
    ]
    pub_types = [
        clean_text(x)
        for x in art.findall("./PublicationTypeList/PublicationType")
        if clean_text(x)
    ]

    return Paper(
        pmid=pmid,
        title=child_text(art, "ArticleTitle"),
        journal=child_text(art, "Journal/Title"),
        publication_date=date,
        year=year,
        doi=doi,
        pmcid=pmcid,
        pubmed_url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
        authors=authors,
        skku_authors=skku_authors,
        skku_orcids=skku_orcids,
        skku_affiliation_evidence=evidence,
        abstract=parse_abstract(article),
        mesh_terms=mesh,
        keywords=keywords,
        publication_types=pub_types,
    )


class PubMedClient:
    def __init__(self, email: str, api_key: str = "", timeout: int = 60):
        self.email = email.strip()
        self.api_key = api_key.strip()
        self.timeout = timeout
        self.last_request = 0.0
        self.min_interval = 0.12 if self.api_key else 0.36

    def _wait(self) -> None:
        elapsed = time.monotonic() - self.last_request
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)

    def _request(
        self,
        endpoint: str,
        params: dict[str, str],
        post: bool = False,
        max_attempts: int = 6,
    ) -> bytes:
        """Call NCBI E-utilities with throttling and transient-error retries.

        Full-corpus crawls make many E-utility calls, so a single 429/5xx/timeout
        must not abort the entire run. Invalid API keys are detected and the
        client falls back to the unauthenticated NCBI rate limit.
        """
        url = f"{EUTILS}/{endpoint}"
        retryable_http = {429, 500, 502, 503, 504}
        dropped_bad_api_key = False
        last_error: Exception | None = None

        for attempt in range(1, max_attempts + 1):
            common = {"tool": TOOL, "email": self.email}
            if self.api_key:
                common["api_key"] = self.api_key
            body = urllib.parse.urlencode({**params, **common}).encode("utf-8")
            self._wait()

            if post:
                req = urllib.request.Request(
                    url,
                    data=body,
                    headers={"User-Agent": f"{TOOL}/1.0 ({self.email})"},
                    method="POST",
                )
            else:
                req = urllib.request.Request(
                    f"{url}?{body.decode('utf-8')}",
                    headers={"User-Agent": f"{TOOL}/1.0 ({self.email})"},
                )

            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as response:
                    payload = response.read()
                self.last_request = time.monotonic()
                return payload
            except urllib.error.HTTPError as exc:
                self.last_request = time.monotonic()
                last_error = exc
                try:
                    detail = exc.read().decode("utf-8", errors="replace")[:500]
                except Exception:
                    detail = ""

                invalid_key = (
                    bool(self.api_key)
                    and exc.code in {400, 401, 403}
                    and "api" in detail.lower()
                    and "key" in detail.lower()
                )
                if invalid_key and not dropped_bad_api_key:
                    print(
                        "[ncbi] API key was rejected; retrying without NCBI_API_KEY",
                        file=sys.stderr,
                    )
                    self.api_key = ""
                    self.min_interval = 0.36
                    dropped_bad_api_key = True
                    continue

                if exc.code not in retryable_http or attempt >= max_attempts:
                    raise RuntimeError(
                        f"NCBI {endpoint} failed with HTTP {exc.code}: "
                        f"{detail or exc.reason}"
                    ) from exc

                retry_after = exc.headers.get("Retry-After") if exc.headers else None
                try:
                    delay = float(retry_after) if retry_after else min(30.0, 2 ** (attempt - 1))
                except ValueError:
                    delay = min(30.0, 2 ** (attempt - 1))
                print(
                    f"[ncbi] HTTP {exc.code} on {endpoint}; "
                    f"retry {attempt}/{max_attempts} after {delay:g}s",
                    file=sys.stderr,
                )
                time.sleep(delay)
            except (
                urllib.error.URLError,
                TimeoutError,
                http.client.IncompleteRead,
                http.client.RemoteDisconnected,
                ConnectionResetError,
                BrokenPipeError,
            ) as exc:
                self.last_request = time.monotonic()
                last_error = exc
                if attempt >= max_attempts:
                    raise RuntimeError(
                        f"NCBI {endpoint} failed after {max_attempts} attempts: {exc}"
                    ) from exc
                delay = min(30.0, 2 ** (attempt - 1))
                print(
                    f"[ncbi] network/timeout on {endpoint}: {exc}; "
                    f"retry {attempt}/{max_attempts} after {delay:g}s",
                    file=sys.stderr,
                )
                time.sleep(delay)

        raise RuntimeError(f"NCBI {endpoint} failed: {last_error}")

    def search_page(
        self,
        query: str,
        retstart: int = 0,
        retmax: int = 500,
    ) -> tuple[int, list[str]]:
        if retstart < 0:
            raise ValueError("retstart must be >= 0")
        if retmax < 0 or retmax > 10000:
            raise ValueError("retmax must be between 0 and 10000")
        if retstart + retmax > 10000:
            raise ValueError(
                "PubMed ESearch only exposes the first 10,000 records per query; "
                "partition the query by date before paging further."
            )
        payload = self._request(
            "esearch.fcgi",
            {
                "db": "pubmed",
                "term": query,
                "retmode": "json",
                "retstart": str(retstart),
                "retmax": str(retmax),
                "sort": "pub_date",
            },
        )
        result = json.loads(payload.decode("utf-8")).get("esearchresult", {})
        return int(result.get("count", 0)), list(result.get("idlist", []))

    def search(self, query: str, max_results: int) -> tuple[int, list[str]]:
        # Backward-compatible one-page search used by author follow-up.
        return self.search_page(query, retstart=0, retmax=max_results)

    def fetch(
        self,
        pmids: Iterable[str],
        batch_size: int = 150,
        cache_dir: Path | None = None,
        resume: bool = False,
    ) -> list[Paper]:
        """Fetch PubMed records in resumable batches.

        When cache_dir is provided, each successful efetch XML batch is written
        atomically to disk. With resume=True, matching cached batches are parsed
        locally instead of being downloaded again. The batch filename includes
        a hash of its PMID list so stale caches cannot be silently reused.
        """
        ids = list(dict.fromkeys(str(x) for x in pmids if str(x).strip()))
        papers: list[Paper] = []
        total_batches = math.ceil(len(ids) / batch_size) if ids else 0

        cache_path = Path(cache_dir) if cache_dir is not None else None
        if cache_path is not None:
            cache_path.mkdir(parents=True, exist_ok=True)

        cache_hits = 0
        downloaded = 0

        for batch_index, start in enumerate(range(0, len(ids), batch_size), 1):
            batch = ids[start : start + batch_size]
            batch_token = hashlib.sha1(",".join(batch).encode("utf-8")).hexdigest()[:16]
            cache_file = (
                cache_path / f"batch_{batch_index:05d}_{batch_token}.xml"
                if cache_path is not None
                else None
            )

            payload: bytes | None = None
            root: ET.Element | None = None
            source = "download"

            if resume and cache_file is not None and cache_file.exists():
                try:
                    payload = cache_file.read_bytes()
                    root = ET.fromstring(payload)
                    cache_hits += 1
                    source = "cache"
                except (OSError, ET.ParseError):
                    # A partial/corrupt cache is never trusted.
                    try:
                        cache_file.unlink()
                    except OSError:
                        pass
                    payload = None
                    root = None

            if payload is None:
                payload = self._request(
                    "efetch.fcgi",
                    {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"},
                    post=len(batch) > 20,
                )
                try:
                    root = ET.fromstring(payload)
                except ET.ParseError as exc:
                    preview = payload[:300].decode("utf-8", errors="replace")
                    raise RuntimeError(
                        f"NCBI efetch returned malformed XML for batch "
                        f"{batch_index}/{total_batches}: {preview}"
                    ) from exc

                downloaded += 1
                if cache_file is not None:
                    tmp_file = cache_file.with_suffix(cache_file.suffix + ".tmp")
                    tmp_file.write_bytes(payload)
                    tmp_file.replace(cache_file)

            assert root is not None
            papers.extend(parse_paper(x) for x in root.findall("./PubmedArticle"))

            if total_batches > 1 and (
                batch_index == 1
                or batch_index == total_batches
                or batch_index % 10 == 0
            ):
                print(
                    f"[fetch] batch {batch_index}/{total_batches}; "
                    f"parsed papers={len(papers):,}; source={source}; "
                    f"cached={cache_hits:,}; downloaded={downloaded:,}",
                    file=sys.stderr,
                )

        if cache_path is not None:
            print(
                f"[fetch] complete; cached batches={cache_hits:,}; "
                f"downloaded batches={downloaded:,}; cache={cache_path}",
                file=sys.stderr,
            )
        return papers

    def links(self, pmid: str, linkname: str) -> set[str]:
        payload = self._request(
            "elink.fcgi",
            {
                "dbfrom": "pubmed",
                "db": "pubmed",
                "id": pmid,
                "linkname": linkname,
                "retmode": "xml",
            },
        )
        root = ET.fromstring(payload)
        result = set()
        for db in root.findall(".//LinkSetDb"):
            if child_text(db, "LinkName") != linkname:
                continue
            for node in db.findall("./Link/Id"):
                value = clean_text(node)
                if value:
                    result.add(value)
        return result


def build_query(start_year: int, end_year: int, topic: str, extra: str) -> str:
    parts = ["(Sungkyunkwan[ad] OR SKKU[ad])"]
    if start_year and end_year:
        parts.append(f"{start_year}:{end_year}[dp]")
    elif start_year:
        parts.append(f"{start_year}:3000[dp]")
    elif end_year:
        parts.append(f"1800:{end_year}[dp]")
    if topic.strip():
        parts.append(f"({topic.strip()})")
    if extra.strip():
        parts.append(f"({extra.strip()})")
    return " AND ".join(parts)


def build_query_date_window(
    start_date: str,
    end_date: str,
    topic: str,
    extra: str,
) -> str:
    parts = [
        "(Sungkyunkwan[ad] OR SKKU[ad])",
        f"{start_date}:{end_date}[dp]",
    ]
    if topic.strip():
        parts.append(f"({topic.strip()})")
    if extra.strip():
        parts.append(f"({extra.strip()})")
    return " AND ".join(parts)


def _checkpoint_payload(
    *,
    spec: dict,
    pmids: list[str],
    completed_windows: list[dict],
) -> dict:
    return {
        "version": 1,
        "spec": spec,
        "retrieved_pmids": pmids,
        "completed_windows": completed_windows,
        "updated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def search_partitioned_pubmed(
    client: PubMedClient,
    *,
    start_year: int,
    end_year: int,
    topic: str = "",
    extra: str = "",
    page_size: int = 500,
    max_results: int = 0,
    checkpoint_path: Path | None = None,
    resume: bool = False,
) -> tuple[int, list[str], list[dict]]:
    """Retrieve an arbitrary-size PubMed result set by date partitioning.

    PubMed ESearch only exposes the first 10,000 records for a single query.
    We therefore partition by publication year; any year above 10,000 hits is
    automatically partitioned by calendar month. Each partition is then paged
    with retstart/retmax. max_results=0 means no user cap.
    """
    if start_year <= 0 or end_year <= 0:
        raise ValueError("--all-results requires explicit positive start/end years")
    if start_year > end_year:
        raise ValueError("start_year must be <= end_year")
    if page_size <= 0 or page_size > 10000:
        raise ValueError("page_size must be between 1 and 10000")
    if max_results < 0:
        raise ValueError("max_results must be >= 0")

    spec = {
        "start_year": start_year,
        "end_year": end_year,
        "topic": topic,
        "extra": extra,
        "page_size": page_size,
        "max_results": max_results,
    }
    pmids: list[str] = []
    seen: set[str] = set()
    completed_windows: list[dict] = []
    completed_labels: set[str] = set()

    if resume and checkpoint_path and checkpoint_path.exists():
        saved = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if saved.get("spec") == spec:
            for pmid in saved.get("retrieved_pmids", []):
                value = str(pmid)
                if value and value not in seen:
                    seen.add(value)
                    pmids.append(value)
            completed_windows = list(saved.get("completed_windows", []))
            completed_labels = {
                str(item.get("label", ""))
                for item in completed_windows
                if item.get("label")
            }
            print(
                f"[resume] checkpoint PMIDs={len(pmids):,}; "
                f"windows={len(completed_windows):,}",
                file=sys.stderr,
            )
        else:
            print("[resume] checkpoint spec differs; starting a new crawl", file=sys.stderr)

    overall_query = build_query(start_year, end_year, topic, extra)
    total_hits, _ = client.search_page(overall_query, retstart=0, retmax=0)

    def persist() -> None:
        if checkpoint_path is None:
            return
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
        checkpoint_path.write_text(
            json.dumps(
                _checkpoint_payload(
                    spec=spec,
                    pmids=pmids,
                    completed_windows=completed_windows,
                ),
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )

    def remaining_cap() -> int | None:
        if max_results == 0:
            return None
        return max(0, max_results - len(pmids))

    def collect_window(query: str, label: str, expected_count: int | None = None) -> int:
        if label in completed_labels:
            prior = next(
                (x for x in completed_windows if x.get("label") == label),
                {},
            )
            return int(prior.get("pubmed_hits", expected_count or 0) or 0)

        cap = remaining_cap()
        if cap == 0:
            return 0

        first_retmax = page_size if cap is None else min(page_size, cap)
        count, first_ids = client.search_page(query, retstart=0, retmax=first_retmax)
        if count > 10000:
            raise RuntimeError(
                f"Partition {label} still has {count:,} hits (>10,000); "
                "split it into smaller date windows."
            )
        if expected_count is not None and count != expected_count:
            print(
                f"[crawl] {label}: count changed {expected_count:,} -> {count:,}",
                file=sys.stderr,
            )

        target = count if cap is None else min(count, cap)
        page_ids = first_ids[:target]
        for pmid in page_ids:
            value = str(pmid)
            if value and value not in seen:
                seen.add(value)
                pmids.append(value)

        retstart = len(first_ids)
        while retstart < target:
            cap_now = remaining_cap()
            if cap_now == 0:
                break
            retmax = min(page_size, target - retstart)
            if cap_now is not None:
                retmax = min(retmax, cap_now)
            _, ids = client.search_page(
                query,
                retstart=retstart,
                retmax=retmax,
            )
            if not ids:
                break
            for pmid in ids:
                value = str(pmid)
                if value and value not in seen:
                    seen.add(value)
                    pmids.append(value)
            retstart += len(ids)
            print(
                f"[crawl] {label}: {min(retstart, target):,}/{target:,}",
                file=sys.stderr,
            )

        completed_windows.append(
            {
                "label": label,
                "query": query,
                "pubmed_hits": count,
                "retrieved_target": target,
                "retrieved_total_after_window": len(pmids),
            }
        )
        completed_labels.add(label)
        persist()
        return count

    # Newest years first, matching PubMed pub_date ordering at the partition level.
    for year in range(end_year, start_year - 1, -1):
        if remaining_cap() == 0:
            break

        year_query = build_query(year, year, topic, extra)
        year_label = f"{year}"
        if year_label in completed_labels:
            prior = next(x for x in completed_windows if x.get("label") == year_label)
            continue

        year_count, _ = client.search_page(year_query, retstart=0, retmax=0)
        print(f"[crawl] year={year} hits={year_count:,}", file=sys.stderr)

        if year_count <= 10000:
            collect_window(year_query, year_label, expected_count=year_count)
            continue

        # Annual partition is still above the PubMed 10k ceiling: split monthly.
        for month in range(12, 0, -1):
            if remaining_cap() == 0:
                break
            last_day = calendar.monthrange(year, month)[1]
            start_date = f"{year}/{month:02d}/01"
            end_date = f"{year}/{month:02d}/{last_day:02d}"
            label = f"{year}-{month:02d}"
            month_query = build_query_date_window(start_date, end_date, topic, extra)
            month_count, _ = client.search_page(month_query, retstart=0, retmax=0)
            if month_count > 10000:
                raise RuntimeError(
                    f"Monthly partition {label} has {month_count:,} hits. "
                    "A daily partition is required for this query."
                )
            collect_window(month_query, label, expected_count=month_count)

    persist()
    return total_hits, pmids, completed_windows


def normalize_term(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def topic_set(paper: Paper) -> set[str]:
    values = paper.mesh_terms + paper.keywords
    return {normalize_term(x) for x in values if normalize_term(x)}


def jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def direction(a: Paper, b: Paper) -> tuple[Paper, Paper]:
    if a.year and b.year and a.year != b.year:
        return (a, b) if a.year < b.year else (b, a)
    try:
        return (a, b) if int(a.pmid) < int(b.pmid) else (b, a)
    except ValueError:
        return (a, b)


def build_edges(
    papers: list[Paper],
    client: PubMedClient,
    include_citations: bool,
    topic_threshold: float,
    max_pairs: int,
) -> list[Edge]:
    by_id = {p.pmid: p for p in papers}
    paper_ids = set(by_id)
    edges: dict[tuple[str, str, str], Edge] = {}

    if include_citations:
        for idx, paper in enumerate(papers, 1):
            refs = client.links(paper.pmid, "pubmed_pubmed_refs")

            for older_id in refs & paper_ids:
                if older_id == paper.pmid:
                    continue
                key = (older_id, paper.pmid, "citation")
                edges[key] = Edge(
                    source=older_id,
                    target=paper.pmid,
                    relation="citation",
                    weight=1.0,
                    evidence=f"{paper.pmid} cites {older_id}",
                )

            if idx % 25 == 0:
                print(f"[citation] {idx}/{len(papers)} papers checked", file=sys.stderr)

    pair_count = 0
    topic_cache = {p.pmid: topic_set(p) for p in papers}
    for i, a in enumerate(papers):
        for b in papers[i + 1 :]:
            pair_count += 1
            if max_pairs and pair_count > max_pairs:
                break

            older, newer = direction(a, b)
            shared_names = sorted(set(a.skku_authors) & set(b.skku_authors))
            shared_orcids = sorted(set(a.skku_orcids) & set(b.skku_orcids))

            if shared_orcids or shared_names:
                weight = 0.95 if shared_orcids else 0.80
                evidence_parts = []
                if shared_orcids:
                    evidence_parts.append("ORCID=" + ", ".join(shared_orcids))
                if shared_names:
                    evidence_parts.append("authors=" + ", ".join(shared_names))
                key = (older.pmid, newer.pmid, "shared_skku_author")
                edges[key] = Edge(
                    source=older.pmid,
                    target=newer.pmid,
                    relation="shared_skku_author",
                    weight=weight,
                    evidence="; ".join(evidence_parts),
                )

            sim = jaccard(topic_cache[a.pmid], topic_cache[b.pmid])
            if sim >= topic_threshold:
                key = (older.pmid, newer.pmid, "topic_similarity")
                edges[key] = Edge(
                    source=older.pmid,
                    target=newer.pmid,
                    relation="topic_similarity",
                    weight=round(sim, 4),
                    evidence=f"MeSH/keyword Jaccard={sim:.3f}",
                )

        if max_pairs and pair_count > max_pairs:
            break

    return sorted(
        edges.values(),
        key=lambda e: (e.source, e.target, e.relation),
    )


def flat_paper(p: Paper) -> dict:
    return {
        "pmid": p.pmid,
        "year": p.year,
        "title": p.title,
        "journal": p.journal,
        "doi": p.doi,
        "pmcid": p.pmcid,
        "pubmed_url": p.pubmed_url,
        "all_authors": "; ".join(a.name for a in p.authors),
        "skku_authors": "; ".join(p.skku_authors),
        "skku_orcids": "; ".join(p.skku_orcids),
        "skku_affiliation_evidence": " | ".join(p.skku_affiliation_evidence),
        "mesh_terms": "; ".join(p.mesh_terms),
        "keywords": "; ".join(p.keywords),
        "abstract": p.abstract,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def connected_components(nodes: set[str], edges: list[Edge]) -> list[list[str]]:
    adj = {n: set() for n in nodes}
    for edge in edges:
        if edge.source in adj and edge.target in adj:
            adj[edge.source].add(edge.target)
            adj[edge.target].add(edge.source)

    seen = set()
    comps = []
    for node in nodes:
        if node in seen:
            continue
        stack = [node]
        seen.add(node)
        comp = []
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nxt in adj[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        comps.append(comp)
    return sorted(comps, key=len, reverse=True)


def write_report(
    path: Path,
    query: str,
    total_hits: int,
    papers: list[Paper],
    edges: list[Edge],
) -> None:
    by_id = {p.pmid: p for p in papers}
    rel_counts = Counter(e.relation for e in edges)
    author_counts = Counter(a for p in papers for a in p.skku_authors)
    comps = connected_components(set(by_id), edges)

    lines = [
        "# SKKU PubMed Research Lineage",
        "",
        f"- Generated UTC: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"- PubMed query: {query}",
        f"- PubMed hits: {total_hits:,}",
        f"- Verified SKKU-affiliated papers: {len(papers):,}",
        f"- Graph edges: {len(edges):,}",
        f"- Connected components: {len(comps):,}",
        "",
        "## Edge types",
        "",
        "| Relation | Count | Meaning |",
        "|---|---:|---|",
        f"| citation | {rel_counts['citation']} | Direct PubMed reference/cited-by link |",
        f"| shared_skku_author | {rel_counts['shared_skku_author']} | Same SKKU-affiliated author/ORCID across papers |",
        f"| topic_similarity | {rel_counts['topic_similarity']} | Similar MeSH/keyword profile |",
        "",
        "## Top SKKU authors",
        "",
        "| Author | Papers |",
        "|---|---:|",
    ]
    for name, count in author_counts.most_common(30):
        lines.append(f"| {name} | {count} |")

    lines += ["", "## Largest research chains", ""]
    for rank, comp in enumerate(comps[:20], 1):
        comp_papers = sorted(
            (by_id[x] for x in comp),
            key=lambda p: (p.year or 9999, p.pmid),
        )
        if len(comp_papers) < 2:
            continue
        lines.append(f"### Chain {rank} ({len(comp_papers)} papers)")
        for p in comp_papers[:40]:
            authors = ", ".join(p.skku_authors[:4]) or "SKKU affiliation"
            lines.append(
                f"- {p.year or 'Unknown'} · PMID {p.pmid} · {p.title} · [{authors}]({p.pubmed_url})"
            )
        lines.append("")

    lines += [
        "## Interpretation",
        "",
        "- Citation edges are the strongest evidence of direct research continuation.",
        "- Shared-author edges can reveal a laboratory/researcher trajectory even without direct citation.",
        "- Topic-similarity edges are exploratory and should not be interpreted as proof of intellectual lineage.",
        "- PubMed affiliation is publication-time affiliation; it does not establish alumni/graduation status.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_html(path: Path, papers: list[Paper], edges: list[Edge]) -> None:
    nodes = []
    for p in papers:
        authors = ", ".join(p.skku_authors[:5])
        label = f"{p.year or '?'} | {p.title[:72]}"
        title = (
            f"PMID {p.pmid}<br><b>{p.title}</b><br>{p.journal}<br>"
            f"SKKU: {authors}<br><a href='{p.pubmed_url}' target='_blank'>PubMed</a>"
        )
        nodes.append(
            {
                "id": p.pmid,
                "label": label,
                "title": title,
                "year": p.year,
                "url": p.pubmed_url,
            }
        )

    vis_edges = []
    for e in edges:
        color = {
            "citation": "#d62728",
            "shared_skku_author": "#1f77b4",
            "topic_similarity": "#7f7f7f",
        }.get(e.relation, "#7f7f7f")
        dashes = e.relation == "topic_similarity"
        vis_edges.append(
            {
                "from": e.source,
                "to": e.target,
                "arrows": "to",
                "label": e.relation,
                "title": e.evidence,
                "value": max(1, int(math.ceil(e.weight * 5))),
                "color": color,
                "dashes": dashes,
            }
        )

    node_json = json.dumps(nodes, ensure_ascii=False).replace("</", "<\\/")
    edge_json = json.dumps(vis_edges, ensure_ascii=False).replace("</", "<\\/")
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SKKU PubMed Research Lineage</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
<style>
body {{ font-family: Arial, sans-serif; margin: 0; }}
header {{ padding: 12px 16px; border-bottom: 1px solid #ddd; }}
#network {{ width: 100vw; height: calc(100vh - 82px); }}
.legend span {{ margin-right: 18px; }}
</style>
</head>
<body>
<header>
  <b>SKKU PubMed Research Lineage</b>
  <div class="legend">
    <span>red: citation</span>
    <span>blue: shared SKKU author</span>
    <span>gray dashed: topic similarity</span>
  </div>
</header>
<div id="network"></div>
<script>
const nodes = new vis.DataSet({node_json});
const edges = new vis.DataSet({edge_json});
const network = new vis.Network(
  document.getElementById("network"),
  {{nodes, edges}},
  {{
    physics: {{stabilization: true}},
    interaction: {{hover: true, navigationButtons: true}},
    edges: {{smooth: {{type: "dynamic"}}}},
    nodes: {{shape: "dot", size: 12, font: {{size: 11}}}}
  }}
);
network.on("doubleClick", function(params) {{
  if (params.nodes.length) {{
    const n = nodes.get(params.nodes[0]);
    if (n && n.url) window.open(n.url, "_blank");
  }}
}});
</script>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build an SKKU PubMed research lineage graph.")
    parser.add_argument("--email", default=os.getenv("NCBI_EMAIL", ""))
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY", ""))
    parser.add_argument("--start-year", type=int, default=2018)
    parser.add_argument("--end-year", type=int, default=datetime.now().year)
    parser.add_argument("--topic", default="")
    parser.add_argument("--extra-query", default="")
    parser.add_argument(
        "--max-results",
        type=int,
        default=500,
        help="Maximum PMIDs to retrieve. With --all-results, 0 means no cap.",
    )
    parser.add_argument(
        "--all-results",
        action="store_true",
        help="Crawl the full date range using year/month partitions and pagination.",
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=500,
        help="ESearch page size for --all-results (1-10000).",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume --all-results from crawl_checkpoint.json when the crawl spec matches.",
    )
    parser.add_argument("--topic-threshold", type=float, default=0.30)
    parser.add_argument("--max-pairs", type=int, default=250000)
    parser.add_argument("--skip-citations", action="store_true")
    parser.add_argument("--output-dir", default="outputs/skku_pubmed_lineage")
    args = parser.parse_args()

    if not args.email:
        parser.error("Provide --email or set NCBI_EMAIL. NCBI asks E-utility clients to identify themselves.")

    query = build_query(args.start_year, args.end_year, args.topic, args.extra_query)
    client = PubMedClient(args.email, args.api_key)

    print(f"[search] {query}", file=sys.stderr)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    crawl_windows: list[dict] = []
    if args.all_results:
        total_hits, pmids, crawl_windows = search_partitioned_pubmed(
            client,
            start_year=args.start_year,
            end_year=args.end_year,
            topic=args.topic,
            extra=args.extra_query,
            page_size=args.page_size,
            max_results=args.max_results,
            checkpoint_path=out / "crawl_checkpoint.json",
            resume=args.resume,
        )
        retrieval_mode = "partitioned_all_results"
    else:
        total_hits, pmids = client.search(query, args.max_results)
        retrieval_mode = "single_page"

    print(
        f"[search] PubMed hits={total_hits:,}; retrieving={len(pmids):,}; "
        f"mode={retrieval_mode}",
        file=sys.stderr,
    )

    fetched = client.fetch(
        pmids,
        cache_dir=out / "efetch_cache",
        resume=args.resume,
    )
    papers = [p for p in fetched if p.skku_affiliation_evidence]
    papers.sort(key=lambda p: (p.year, p.pmid), reverse=True)
    print(f"[verify] SKKU affiliation verified={len(papers):,}", file=sys.stderr)

    edges = build_edges(
        papers=papers,
        client=client,
        include_citations=not args.skip_citations,
        topic_threshold=args.topic_threshold,
        max_pairs=args.max_pairs,
    )
    print(f"[graph] edges={len(edges):,}", file=sys.stderr)

    (out / "query.txt").write_text(query + "\n", encoding="utf-8")
    (out / "papers.json").write_text(
        json.dumps([asdict(p) for p in papers], ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (out / "edges.json").write_text(
        json.dumps([asdict(e) for e in edges], ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (out / "crawl_windows.json").write_text(
        json.dumps(crawl_windows, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_csv(out / "papers.csv", [flat_paper(p) for p in papers])
    write_csv(out / "edges.csv", [asdict(e) for e in edges])
    write_report(out / "lineage.md", query, total_hits, papers, edges)
    write_html(out / "lineage.html", papers, edges)

    summary = {
        "query": query,
        "pubmed_total_hits": total_hits,
        "retrieval_mode": retrieval_mode,
        "retrieved_pmids": len(pmids),
        "crawl_windows": len(crawl_windows),
        "verified_papers": len(papers),
        "edges": len(edges),
        "edge_types": dict(Counter(e.relation for e in edges)),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    (out / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
