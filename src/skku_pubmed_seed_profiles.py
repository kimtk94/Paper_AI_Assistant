#!/usr/bin/env python3
"""Build scalable researcher profiles directly from the full SKKU PubMed seed corpus.

This path is intended for institution-wide mapping. It avoids issuing a new PubMed
search for every researcher. Instead, it derives researcher identity, publication
profiles, and co-publication membership from SKKU-affiliated papers already fetched
by skku_pubmed_lineage.py.

Identity policy:
- High confidence: ORCID.
- Medium confidence: a name-only record is mapped to an ORCID only when that exact
  normalized name maps uniquely to one ORCID anywhere in the seed corpus.
- Low confidence: remaining exact normalized full-name identity.

Outputs are shaped to be compatible with skku_pubmed_researcher_network.py.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

from skku_pubmed_annotation import annotate_paper
from skku_pubmed_lineage import Author, Paper, is_skku_author, skku_author_match_status


@dataclass
class SeedResearcherProfile:
    key: str
    name: str
    orcid: str
    confidence: str
    paper_count: int
    first_year: int
    last_year: int
    top_diseases: list[str]
    top_methods: list[str]
    top_data_types: list[str]
    stage_path: list[str]
    strong_lineage_count: int
    trajectory_summary: str


def norm_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def load_seed_papers(path: Path) -> list[Paper]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, list):
        raise ValueError("papers.json must contain a list")

    papers: list[Paper] = []
    for row in raw:
        authors = []
        for a in row.get("authors", []):
            name = str(a.get("name", ""))
            initials = str(a.get("initials", ""))
            affiliations = list(a.get("affiliations", []) or [])
            authors.append(
                Author(
                    name=name,
                    last_name=str(a.get("last_name", "")),
                    fore_name=str(a.get("fore_name", "")),
                    initials=initials,
                    orcid=str(a.get("orcid", "")).replace("https://orcid.org/", ""),
                    affiliations=affiliations,
                    # Recompute from raw affiliation text. Do not trust the legacy
                    # is_skku boolean embedded in papers.json.
                    is_skku=is_skku_author(name, initials, affiliations),
                )
            )
        papers.append(
            Paper(
                pmid=str(row.get("pmid", "")),
                title=str(row.get("title", "")),
                journal=str(row.get("journal", "")),
                publication_date=str(row.get("publication_date", "")),
                year=int(row.get("year", 0) or 0),
                doi=str(row.get("doi", "")),
                pmcid=str(row.get("pmcid", "")),
                pubmed_url=str(row.get("pubmed_url", "")),
                authors=authors,
                skku_authors=[a.name for a in authors if a.is_skku],
                skku_orcids=sorted({a.orcid for a in authors if a.is_skku and a.orcid}),
                skku_affiliation_evidence=[
                    f"{a.name}: {aff}"
                    for a in authors if a.is_skku
                    for aff in a.affiliations
                    if "sungkyunkwan" in aff.lower() or re.search(r"\\bskku\\b", aff, re.I)
                ],
                abstract=str(row.get("abstract", "")),
                mesh_terms=list(row.get("mesh_terms", []) or []),
                keywords=list(row.get("keywords", []) or []),
                publication_types=list(row.get("publication_types", []) or []),
            )
        )
    return papers


def resolve_identity_map(papers: list[Paper]) -> dict[tuple[str, str], tuple[str, str]]:
    """Return (pmid, normalized-name) -> (identity-key, confidence)."""
    name_orcids: dict[str, set[str]] = defaultdict(set)
    preferred_name: dict[str, str] = {}

    for paper in papers:
        for author in paper.authors:
            if not author.is_skku:
                continue
            name_key = norm_name(author.name)
            if not name_key:
                continue
            preferred_name.setdefault(name_key, author.name)
            if author.orcid:
                name_orcids[name_key].add(author.orcid)

    unique_name_orcid = {
        name: next(iter(orcids))
        for name, orcids in name_orcids.items()
        if len(orcids) == 1
    }

    identity: dict[tuple[str, str], tuple[str, str]] = {}
    for paper in papers:
        for author in paper.authors:
            if not author.is_skku:
                continue
            name_key = norm_name(author.name)
            if not name_key:
                continue
            if author.orcid:
                identity[(paper.pmid, name_key)] = (f"orcid:{author.orcid}", "high")
            elif name_key in unique_name_orcid:
                identity[(paper.pmid, name_key)] = (
                    f"orcid:{unique_name_orcid[name_key]}",
                    "medium",
                )
            else:
                identity[(paper.pmid, name_key)] = (f"name:{name_key}", "low")
    return identity


def build_seed_profiles(
    papers: list[Paper],
) -> tuple[list[SeedResearcherProfile], list[dict], dict]:
    identity_map = resolve_identity_map(papers)
    annotations = {paper.pmid: annotate_paper(paper) for paper in papers}

    per_key_papers: dict[str, set[str]] = defaultdict(set)
    per_key_names: dict[str, Counter] = defaultdict(Counter)
    per_key_orcid: dict[str, str] = {}
    per_key_confidence: dict[str, set[str]] = defaultdict(set)

    lineage_papers: list[dict] = []

    for paper in papers:
        keys = []
        names = []
        seen = set()
        for author in paper.authors:
            if not author.is_skku:
                continue
            name_key = norm_name(author.name)
            if not name_key:
                continue
            key, confidence = identity_map[(paper.pmid, name_key)]
            if key in seen:
                continue
            seen.add(key)
            keys.append(key)
            names.append(author.name)
            per_key_papers[key].add(paper.pmid)
            per_key_names[key][author.name] += 1
            per_key_confidence[key].add(confidence)
            if key.startswith("orcid:"):
                per_key_orcid[key] = key.split(":", 1)[1]

        ann = annotations[paper.pmid]
        lineage_papers.append(
            {
                "pmid": paper.pmid,
                "year": paper.year,
                "title": paper.title,
                "journal": paper.journal,
                "doi": paper.doi,
                "pubmed_url": paper.pubmed_url,
                "tracked_authors": names,
                "tracked_author_keys": keys,
                "current_affiliations": paper.skku_affiliation_evidence,
                "skku_current": True,
                "disease_terms": ann.disease_terms,
                "methods": ann.methods,
                "data_types": ann.data_types,
                "research_stage": ann.research_stage,
                "research_question": ann.research_question,
            }
        )

    paper_by_id = {paper.pmid: paper for paper in papers}
    profiles: list[SeedResearcherProfile] = []

    for key, pmids in per_key_papers.items():
        ordered = sorted(
            (paper_by_id[pmid] for pmid in pmids),
            key=lambda p: (p.year or 9999, int(p.pmid) if p.pmid.isdigit() else 0),
        )
        disease_counts = Counter(
            value
            for paper in ordered
            for value in annotations[paper.pmid].disease_terms
        )
        method_counts = Counter(
            value
            for paper in ordered
            for value in annotations[paper.pmid].methods
        )
        data_counts = Counter(
            value
            for paper in ordered
            for value in annotations[paper.pmid].data_types
        )
        stage_path = []
        for paper in ordered:
            stage = annotations[paper.pmid].research_stage
            if not stage_path or stage_path[-1] != stage:
                stage_path.append(stage)

        confidence_set = per_key_confidence[key]
        if "high" in confidence_set:
            confidence = "high"
        elif "medium" in confidence_set:
            confidence = "medium"
        else:
            confidence = "low"

        name = per_key_names[key].most_common(1)[0][0]
        years = [paper.year for paper in ordered if paper.year]
        top_diseases = [x for x, _ in disease_counts.most_common(6)]
        top_methods = [x for x, _ in method_counts.most_common(6)]
        top_data = [x for x, _ in data_counts.most_common(6)]

        parts = []
        if years:
            parts.append(f"{min(years)}–{max(years)}")
        if stage_path:
            parts.append("stage: " + " → ".join(stage_path[:8]))
        if top_methods:
            parts.append("methods: " + ", ".join(top_methods[:3]))
        if top_diseases:
            parts.append("disease: " + ", ".join(top_diseases[:2]))

        profiles.append(
            SeedResearcherProfile(
                key=key,
                name=name,
                orcid=per_key_orcid.get(key, ""),
                confidence=confidence,
                paper_count=len(pmids),
                first_year=min(years) if years else 0,
                last_year=max(years) if years else 0,
                top_diseases=top_diseases,
                top_methods=top_methods,
                top_data_types=top_data,
                stage_path=stage_path,
                strong_lineage_count=0,
                trajectory_summary=" · ".join(parts),
            )
        )

    profiles.sort(key=lambda x: (-x.paper_count, x.name.lower()))

    affiliation_status_counts = Counter(
        skku_author_match_status(author.name, author.initials, author.affiliations)
        for paper in papers
        for author in paper.authors
    )
    ambiguous_papers = sum(
        any(
            skku_author_match_status(author.name, author.initials, author.affiliations)
            == "ambiguous"
            for author in paper.authors
        )
        for paper in papers
    )
    tracked_counts = [len(x["tracked_author_keys"]) for x in lineage_papers]

    summary = {
        "profile_version": 2,
        "author_affiliation_policy": "conservative_shared_block_initial_attribution",
        "seed_papers": len(papers),
        "researchers": len(profiles),
        "orcid_researchers": sum(bool(x.orcid) for x in profiles),
        "name_only_researchers": sum(not x.orcid for x in profiles),
        "high_confidence": sum(x.confidence == "high" for x in profiles),
        "medium_confidence": sum(x.confidence == "medium" for x in profiles),
        "low_confidence": sum(x.confidence == "low" for x in profiles),
        "papers_with_multiple_skku_researchers": sum(
            len(x["tracked_author_keys"]) >= 2 for x in lineage_papers
        ),
        "papers_without_attributed_skku_researcher": sum(x == 0 for x in tracked_counts),
        "max_tracked_researchers_per_paper": max(tracked_counts, default=0),
        "ambiguous_shared_affiliation_papers": ambiguous_papers,
        "affiliation_status_counts": dict(affiliation_status_counts),
    }
    return profiles, lineage_papers, summary


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build institution-wide researcher profiles directly from SKKU seed papers."
    )
    parser.add_argument("--papers-json", default="outputs/skku_pubmed_lineage/papers.json")
    parser.add_argument("--output-dir", default="outputs/skku_pubmed_seed_map")
    args = parser.parse_args()

    papers = load_seed_papers(Path(args.papers_json))
    profiles, lineage_papers, summary = build_seed_profiles(papers)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    profile_json = [asdict(x) for x in profiles]
    (out / "researcher_summary.json").write_text(
        json.dumps(profile_json, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (out / "lineage_papers.json").write_text(
        json.dumps(lineage_papers, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (out / "lineage_edges.json").write_text("[]\n", encoding="utf-8")

    profile_rows = []
    for item in profiles:
        row = asdict(item)
        for key in ["top_diseases", "top_methods", "top_data_types", "stage_path"]:
            row[key] = "; ".join(row[key])
        profile_rows.append(row)
    write_csv(out / "researcher_summary.csv", profile_rows)

    paper_rows = []
    for item in lineage_papers:
        row = dict(item)
        for key in [
            "tracked_authors", "tracked_author_keys", "current_affiliations",
            "disease_terms", "methods", "data_types",
        ]:
            row[key] = "; ".join(row[key])
        paper_rows.append(row)
    write_csv(out / "lineage_papers.csv", paper_rows)
    (out / "lineage_edges.csv").write_text("", encoding="utf-8")

    (out / "seed_profile_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
