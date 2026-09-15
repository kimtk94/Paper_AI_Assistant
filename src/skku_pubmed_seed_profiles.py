#!/usr/bin/env python3
"""Build scalable researcher profiles directly from the full SKKU PubMed seed corpus.

This path is intended for institution-wide mapping. It avoids issuing a new PubMed
search for every researcher. Instead, it derives researcher identity, publication
profiles, and co-publication membership from SKKU-affiliated papers already fetched
by skku_pubmed_lineage.py.

Identity policy:
- High confidence: ORCID.
- Medium confidence: a name-only record maps to a unique ORCID only when its
  affiliation fingerprint is compatible with that ORCID's observed affiliations.
- Low-affiliation confidence: remaining name-only records are split by a
  conservative SKKU department/institute/hospital fingerprint instead of being
  merged institution-wide by name alone.

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
    top_domains: list[str]
    top_topics: list[str]
    top_diseases: list[str]
    top_methods: list[str]
    top_data_types: list[str]
    stage_path: list[str]
    strong_lineage_count: int
    trajectory_summary: str


def norm_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def affiliation_fingerprint(author: Author) -> str:
    """Return a conservative SKKU sub-affiliation fingerprint for name-only IDs."""
    candidates = [
        aff for aff in author.affiliations
        if "sungkyunkwan" in aff.lower() or re.search(r"\bskku\b", aff, re.I)
    ]
    if not candidates:
        return ""

    hospital_anchors = (
        "samsung medical center",
        "kangbuk samsung hospital",
        "samsung changwon hospital",
    )
    for aff in candidates:
        low = re.sub(r"\s+", " ", aff.lower())
        parts = []
        for hospital in hospital_anchors:
            if hospital in low:
                parts.append(hospital.replace(" ", "_"))
                break

        match = re.search(
            r"\b(department|division|school|college|institute|center|centre)\s+of\s+([^,;]{2,80})",
            low,
        )
        if match:
            unit = re.sub(r"[^a-z0-9]+", "_", match.group(2)).strip("_")
            if unit:
                parts.append(f"{match.group(1)}_{unit}")

        if parts:
            return "|".join(parts[:2])

    return "skku_general"


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
                    if "sungkyunkwan" in aff.lower() or re.search(r"\bskku\b", aff, re.I)
                ],
                abstract=str(row.get("abstract", "")),
                mesh_terms=list(row.get("mesh_terms", []) or []),
                keywords=list(row.get("keywords", []) or []),
                publication_types=list(row.get("publication_types", []) or []),
            )
        )
    return papers


def resolve_identity_map(papers: list[Paper]) -> dict[tuple[str, str], tuple[str, str]]:
    """Return (pmid, normalized-name) -> (identity-key, confidence).

    ORCID remains authoritative. Name-only records are only attached to a unique
    ORCID when affiliation evidence is compatible; otherwise they are split by a
    conservative SKKU affiliation fingerprint.
    """
    name_orcids: dict[str, set[str]] = defaultdict(set)
    name_orcid_fingerprints: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )

    for paper in papers:
        for author in paper.authors:
            if not author.is_skku:
                continue
            name_key = norm_name(author.name)
            if not name_key:
                continue
            if author.orcid:
                name_orcids[name_key].add(author.orcid)
                fp = affiliation_fingerprint(author)
                if fp:
                    name_orcid_fingerprints[name_key][author.orcid].add(fp)

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
                identity[(paper.pmid, name_key)] = (
                    f"orcid:{author.orcid}",
                    "high",
                )
                continue

            fp = affiliation_fingerprint(author)
            if name_key in unique_name_orcid:
                orcid = unique_name_orcid[name_key]
                known_fps = name_orcid_fingerprints[name_key].get(orcid, set())
                if not fp or not known_fps or fp in known_fps:
                    identity[(paper.pmid, name_key)] = (
                        f"orcid:{orcid}",
                        "medium",
                    )
                    continue

            fp_key = fp or "unknown"
            identity[(paper.pmid, name_key)] = (
                f"name:{name_key}|aff:{fp_key}",
                "low_affiliation" if fp else "low",
            )

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
                "primary_domain": ann.primary_domain,
                "research_domains": ann.research_domains,
                "topic_terms": ann.topic_terms,
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
        domain_counts = Counter(
            value
            for paper in ordered
            for value in annotations[paper.pmid].research_domains
        )
        topic_counts = Counter(
            value
            for paper in ordered
            for value in annotations[paper.pmid].topic_terms
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
        elif "low_affiliation" in confidence_set:
            confidence = "low_affiliation"
        else:
            confidence = "low"

        name = per_key_names[key].most_common(1)[0][0]
        years = [paper.year for paper in ordered if paper.year]
        top_domains = [x for x, _ in domain_counts.most_common(6)]
        top_topics = [x for x, _ in topic_counts.most_common(8)]
        top_diseases = [x for x, _ in disease_counts.most_common(6)]
        top_methods = [x for x, _ in method_counts.most_common(6)]
        top_data = [x for x, _ in data_counts.most_common(6)]

        parts = []
        if years:
            parts.append(f"{min(years)}–{max(years)}")
        if top_domains:
            parts.append("domain: " + ", ".join(top_domains[:2]))
        if top_topics:
            parts.append("topics: " + ", ".join(top_topics[:2]))
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
                top_domains=top_domains,
                top_topics=top_topics,
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
    summary_name_orcids: dict[str, set[str]] = defaultdict(set)
    for paper in papers:
        for author in paper.authors:
            if author.is_skku and author.orcid:
                name_key = norm_name(author.name)
                if name_key:
                    summary_name_orcids[name_key].add(author.orcid)

    summary = {
        "profile_version": 3,
        "author_affiliation_policy": "conservative_shared_block_plus_affiliation_fingerprint_identity",
        "seed_papers": len(papers),
        "researchers": len(profiles),
        "orcid_researchers": sum(bool(x.orcid) for x in profiles),
        "name_only_researchers": sum(not x.orcid for x in profiles),
        "high_confidence": sum(x.confidence == "high" for x in profiles),
        "medium_confidence": sum(x.confidence == "medium" for x in profiles),
        "low_confidence": sum(x.confidence.startswith("low") for x in profiles),
        "affiliation_disambiguated_name_profiles": sum(
            x.key.startswith("name:") and "|aff:" in x.key and not x.key.endswith("|aff:unknown")
            for x in profiles
        ),
        "ambiguous_names_with_multiple_orcids": sum(
            len(orcids) >= 2 for orcids in summary_name_orcids.values()
        ),
        "papers_with_multiple_skku_researchers": sum(
            len(x["tracked_author_keys"]) >= 2 for x in lineage_papers
        ),
        "papers_without_attributed_skku_researcher": sum(x == 0 for x in tracked_counts),
        "max_tracked_researchers_per_paper": max(tracked_counts, default=0),
        "ambiguous_shared_affiliation_papers": ambiguous_papers,
        "affiliation_status_counts": dict(affiliation_status_counts),
        "primary_domain_counts": dict(
            Counter(annotations[p.pmid].primary_domain for p in papers)
        ),
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
        for key in ["top_domains", "top_topics", "top_diseases", "top_methods", "top_data_types", "stage_path"]:
            row[key] = "; ".join(row[key])
        profile_rows.append(row)
    write_csv(out / "researcher_summary.csv", profile_rows)

    paper_rows = []
    for item in lineage_papers:
        row = dict(item)
        for key in [
            "tracked_authors", "tracked_author_keys", "current_affiliations",
            "research_domains", "topic_terms", "disease_terms", "methods", "data_types",
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
