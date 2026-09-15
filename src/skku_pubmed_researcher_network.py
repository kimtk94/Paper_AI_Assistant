#!/usr/bin/env python3
"""Build a researcher-level network from SKKU PubMed lineage outputs.

This layer converts paper-level lineage into a researcher graph using:
- co-publication among tracked researchers,
- direct PubMed citation relationships between tracked researchers,
- overlap in disease, method, and data-type profiles.

It is deterministic and does not call external APIs.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path


@dataclass
class ResearcherNode:
    key: str
    name: str
    orcid: str
    paper_count: int
    first_year: int
    last_year: int
    top_diseases: list[str]
    top_methods: list[str]
    top_data_types: list[str]
    stage_path: list[str]
    strong_lineage_count: int


@dataclass
class ResearcherEdge:
    source: str
    target: str
    source_name: str
    target_name: str
    relation: str
    score: float
    shared_papers: int
    direct_citations: int
    topic_similarity: float
    shared_diseases: list[str]
    shared_methods: list[str]
    shared_data_types: list[str]
    evidence: str


def load_json(path: Path) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"Expected a JSON list: {path}")
    return data


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def profile_similarity(a: ResearcherNode, b: ResearcherNode) -> tuple[float, list[str], list[str], list[str]]:
    disease_a, disease_b = set(a.top_diseases), set(b.top_diseases)
    method_a, method_b = set(a.top_methods), set(b.top_methods)
    data_a, data_b = set(a.top_data_types), set(b.top_data_types)

    shared_disease = sorted(disease_a & disease_b)
    shared_methods = sorted(method_a & method_b)
    shared_data = sorted(data_a & data_b)

    score = (
        0.40 * jaccard(disease_a, disease_b)
        + 0.35 * jaccard(method_a, method_b)
        + 0.25 * jaccard(data_a, data_b)
    )
    return round(score, 4), shared_disease, shared_methods, shared_data


def build_nodes(researchers: list[dict]) -> list[ResearcherNode]:
    nodes = []
    for row in researchers:
        nodes.append(
            ResearcherNode(
                key=str(row.get("key", "")),
                name=str(row.get("name", "")),
                orcid=str(row.get("orcid", "")),
                paper_count=int(row.get("paper_count", 0) or 0),
                first_year=int(row.get("first_year", 0) or 0),
                last_year=int(row.get("last_year", 0) or 0),
                top_diseases=list(row.get("top_diseases", []) or []),
                top_methods=list(row.get("top_methods", []) or []),
                top_data_types=list(row.get("top_data_types", []) or []),
                stage_path=list(row.get("stage_path", []) or []),
                strong_lineage_count=int(row.get("strong_lineage_count", 0) or 0),
            )
        )
    return [x for x in nodes if x.key and x.name]


def build_researcher_network(
    researchers: list[dict],
    papers: list[dict],
    lineage_edges: list[dict],
    topic_threshold: float = 0.30,
    include_thematic: bool = True,
) -> tuple[list[ResearcherNode], list[ResearcherEdge]]:
    nodes = build_nodes(researchers)
    by_key = {x.key: x for x in nodes}
    name_to_key = {x.name: x.key for x in nodes}

    collaboration: dict[tuple[str, str], set[str]] = {}
    for paper in papers:
        keys = sorted(
            {
                str(x)
                for x in (paper.get("tracked_author_keys", []) or [])
                if str(x) in by_key
            }
        )
        pmid = str(paper.get("pmid", ""))
        for left, right in combinations(keys, 2):
            collaboration.setdefault((left, right), set()).add(pmid)

    citation_counts: dict[tuple[str, str], int] = {}
    citation_examples: dict[tuple[str, str], list[str]] = {}
    for edge in lineage_edges:
        if edge.get("relation") != "cross_researcher_citation":
            continue
        names = sorted({str(x) for x in (edge.get("tracked_authors", []) or []) if str(x) in name_to_key})
        keys = sorted({name_to_key[name] for name in names})
        for left, right in combinations(keys, 2):
            pair = (left, right)
            citation_counts[pair] = citation_counts.get(pair, 0) + 1
            citation_examples.setdefault(pair, []).append(
                f'{edge.get("source", "")}->{edge.get("target", "")}'
            )

    edges: list[ResearcherEdge] = []
    keys = sorted(by_key)

    if include_thematic:
        candidate_pairs = [
            (left, right)
            for i, left in enumerate(keys)
            for right in keys[i + 1 :]
        ]
    else:
        candidate_pairs = sorted(set(collaboration) | set(citation_counts))

    for left, right in candidate_pairs:
        a, b = by_key[left], by_key[right]
        pair = (left, right)
        shared_pmids = sorted(x for x in collaboration.get(pair, set()) if x)
        citations = citation_counts.get(pair, 0)
        similarity, shared_disease, shared_methods, shared_data = profile_similarity(a, b)

        if not shared_pmids and citations == 0 and (
            not include_thematic or similarity < topic_threshold
        ):
            continue

        if shared_pmids and citations:
            relation = "collaboration+citation"
            score = min(1.0, 0.95 + 0.05 * similarity)
        elif citations:
            relation = "citation"
            score = min(1.0, 0.85 + 0.10 * similarity)
        elif shared_pmids:
            relation = "collaboration"
            score = min(0.94, 0.75 + 0.15 * similarity + 0.02 * min(len(shared_pmids), 5))
        else:
            relation = "thematic_overlap"
            score = min(0.84, 0.45 + 0.45 * similarity)

        evidence_parts = []
        if shared_pmids:
            evidence_parts.append(
                f"shared papers={len(shared_pmids)} ({', '.join(shared_pmids[:5])})"
            )
        if citations:
            examples = citation_examples.get(pair, [])
            evidence_parts.append(
                f"cross-researcher citations={citations}"
                + (f" ({', '.join(examples[:5])})" if examples else "")
            )
        if shared_disease:
            evidence_parts.append("shared disease=" + ", ".join(shared_disease[:4]))
        if shared_methods:
            evidence_parts.append("shared method=" + ", ".join(shared_methods[:4]))
        if shared_data:
            evidence_parts.append("shared data=" + ", ".join(shared_data[:4]))
        evidence_parts.append(f"profile similarity={similarity:.3f}")

        edges.append(
            ResearcherEdge(
                source=left,
                target=right,
                source_name=a.name,
                target_name=b.name,
                relation=relation,
                score=round(score, 4),
                shared_papers=len(shared_pmids),
                direct_citations=citations,
                topic_similarity=similarity,
                shared_diseases=shared_disease,
                shared_methods=shared_methods,
                shared_data_types=shared_data,
                evidence="; ".join(evidence_parts),
            )
        )

    edges.sort(key=lambda e: (-e.score, -e.direct_citations, -e.shared_papers, e.source_name, e.target_name))
    return nodes, edges


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def render_html(path: Path, nodes: list[ResearcherNode], edges: list[ResearcherEdge]) -> None:
    vis_nodes = []
    for item in nodes:
        tooltip = (
            f"<b>{item.name}</b><br>"
            f"ORCID: {item.orcid or '-'}<br>"
            f"Papers: {item.paper_count}<br>"
            f"Years: {item.first_year or '?'}–{item.last_year or '?'}<br>"
            f"Strong lineage: {item.strong_lineage_count}<br>"
            f"Disease: {', '.join(item.top_diseases[:4]) or '-'}<br>"
            f"Methods: {', '.join(item.top_methods[:4]) or '-'}<br>"
            f"Data: {', '.join(item.top_data_types[:4]) or '-'}<br>"
            f"Stage path: {' → '.join(item.stage_path) or '-'}"
        )
        vis_nodes.append(
            {
                "id": item.key,
                "label": item.name,
                "title": tooltip,
                "value": max(1, item.paper_count + 2 * item.strong_lineage_count),
            }
        )

    vis_edges = []
    for idx, item in enumerate(edges):
        vis_edges.append(
            {
                "id": f"e{idx}",
                "from": item.source,
                "to": item.target,
                "value": max(1, round(item.score * 6, 2)),
                "title": f"{item.relation} | score={item.score:.2f}<br>{item.evidence}",
                "relation": item.relation,
                "score": item.score,
            }
        )

    nodes_json = json.dumps(vis_nodes, ensure_ascii=False)
    edges_json = json.dumps(vis_edges, ensure_ascii=False)
    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>SKKU PubMed Researcher Network</title>
<script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
<style>
body {{ margin:0; font-family:Arial,sans-serif; background:#fafafa; }}
header {{ padding:14px 18px; background:white; border-bottom:1px solid #ddd; }}
.controls {{ display:flex; gap:16px; flex-wrap:wrap; align-items:center; }}
#network {{ height:82vh; background:white; }}
select,input {{ padding:6px; }}
</style>
</head>
<body>
<header>
<h2 style="margin:0 0 10px 0">SKKU-seeded PubMed Researcher Network</h2>
<div class="controls">
<label>Relation:
<select id="relation">
<option value="">All</option>
<option>collaboration+citation</option>
<option>citation</option>
<option>collaboration</option>
<option>thematic_overlap</option>
</select>
</label>
<label>Minimum score: <input id="minScore" type="number" min="0" max="1" step="0.05" value="0.60"></label>
<span id="summary"></span>
</div>
</header>
<div id="network"></div>
<script>
const allNodes = {nodes_json};
const allEdges = {edges_json};
const network = new vis.Network(
  document.getElementById("network"),
  {{nodes:new vis.DataSet([]), edges:new vis.DataSet([])}},
  {{
    physics: {{stabilization:true, barnesHut:{{gravitationalConstant:-10000}}}},
    interaction: {{hover:true, navigationButtons:true}},
    nodes: {{shape:"dot", scaling:{{min:10,max:38}}}},
    edges: {{smooth:{{type:"dynamic"}}, scaling:{{min:1,max:8}}}}
  }}
);

function redraw() {{
  const relation = document.getElementById("relation").value;
  const minScore = parseFloat(document.getElementById("minScore").value || "0");
  const edges = allEdges.filter(e => (!relation || e.relation === relation) && e.score >= minScore);
  const ids = new Set();
  for (const e of edges) {{ ids.add(e.from); ids.add(e.to); }}
  const nodes = allNodes.filter(n => ids.has(n.id));
  network.setData({{nodes:new vis.DataSet(nodes), edges:new vis.DataSet(edges)}});
  document.getElementById("summary").textContent = "researchers=" + nodes.length + " · edges=" + edges.length;
  if (nodes.length) network.fit();
}}
document.getElementById("relation").addEventListener("change", redraw);
document.getElementById("minScore").addEventListener("change", redraw);
redraw();
</script>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build researcher-level network from SKKU PubMed lineage outputs.")
    parser.add_argument("--input-dir", default="outputs/skku_pubmed_followup")
    parser.add_argument("--topic-threshold", type=float, default=0.30)
    parser.add_argument(
        "--skip-thematic",
        action="store_true",
        help="Only evaluate observed collaboration/citation pairs; avoids O(R^2) thematic comparisons.",
    )
    parser.add_argument("--output-dir", default="")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    researchers = load_json(input_dir / "researcher_summary.json")
    papers = load_json(input_dir / "lineage_papers.json")
    lineage_edges = load_json(input_dir / "lineage_edges.json")

    nodes, edges = build_researcher_network(
        researchers,
        papers,
        lineage_edges,
        topic_threshold=args.topic_threshold,
        include_thematic=not args.skip_thematic,
    )

    node_json = [asdict(x) for x in nodes]
    edge_json = [asdict(x) for x in edges]
    (output_dir / "researcher_network_nodes.json").write_text(
        json.dumps(node_json, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (output_dir / "researcher_network_edges.json").write_text(
        json.dumps(edge_json, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    node_rows = []
    for item in nodes:
        row = asdict(item)
        for key in ["top_diseases", "top_methods", "top_data_types", "stage_path"]:
            row[key] = "; ".join(row[key])
        node_rows.append(row)
    edge_rows = []
    for item in edges:
        row = asdict(item)
        for key in ["shared_diseases", "shared_methods", "shared_data_types"]:
            row[key] = "; ".join(row[key])
        edge_rows.append(row)

    write_csv(output_dir / "researcher_network_nodes.csv", node_rows)
    write_csv(output_dir / "researcher_network_edges.csv", edge_rows)
    render_html(output_dir / "researcher_network.html", nodes, edges)

    summary = {
        "researchers": len(nodes),
        "network_edges": len(edges),
        "collaboration_citation": sum(e.relation == "collaboration+citation" for e in edges),
        "citation": sum(e.relation == "citation" for e in edges),
        "collaboration": sum(e.relation == "collaboration" for e in edges),
        "thematic_overlap": sum(e.relation == "thematic_overlap" for e in edges),
        "high_confidence_edges": sum(e.score >= 0.85 for e in edges),
        "thematic_enabled": not args.skip_thematic,
    }
    (output_dir / "researcher_network_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
