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
import csv
import json
import math
import os
import re
import sys
import time
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

    def _request(self, endpoint: str, params: dict[str, str], post: bool = False) -> bytes:
        common = {"tool": TOOL, "email": self.email}
        if self.api_key:
            common["api_key"] = self.api_key
        body = urllib.parse.urlencode({**params, **common}).encode("utf-8")
        url = f"{EUTILS}/{endpoint}"
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
        with urllib.request.urlopen(req, timeout=self.timeout) as response:
            payload = response.read()
        self.last_request = time.monotonic()
        return payload

    def search(self, query: str, max_results: int) -> tuple[int, list[str]]:
        payload = self._request(
            "esearch.fcgi",
            {
                "db": "pubmed",
                "term": query,
                "retmode": "json",
                "retmax": str(max_results),
                "sort": "pub_date",
            },
        )
        result = json.loads(payload.decode("utf-8")).get("esearchresult", {})
        return int(result.get("count", 0)), list(result.get("idlist", []))

    def fetch(self, pmids: Iterable[str], batch_size: int = 150) -> list[Paper]:
        ids = list(dict.fromkeys(str(x) for x in pmids if str(x).strip()))
        papers = []
        for start in range(0, len(ids), batch_size):
            batch = ids[start : start + batch_size]
            payload = self._request(
                "efetch.fcgi",
                {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"},
                post=len(batch) > 20,
            )
            root = ET.fromstring(payload)
            papers.extend(parse_paper(x) for x in root.findall("./PubmedArticle"))
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
    parser.add_argument("--max-results", type=int, default=500)
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
    total_hits, pmids = client.search(query, args.max_results)
    print(f"[search] PubMed hits={total_hits:,}; retrieving={len(pmids):,}", file=sys.stderr)

    fetched = client.fetch(pmids)
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

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    (out / "query.txt").write_text(query + "\n", encoding="utf-8")
    (out / "papers.json").write_text(
        json.dumps([asdict(p) for p in papers], ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (out / "edges.json").write_text(
        json.dumps([asdict(e) for e in edges], ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_csv(out / "papers.csv", [flat_paper(p) for p in papers])
    write_csv(out / "edges.csv", [asdict(e) for e in edges])
    write_report(out / "lineage.md", query, total_hits, papers, edges)
    write_html(out / "lineage.html", papers, edges)

    summary = {
        "query": query,
        "pubmed_total_hits": total_hits,
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
