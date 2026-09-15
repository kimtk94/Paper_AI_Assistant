#!/usr/bin/env python3
"""Cluster the SKKU researcher network into interpretable research communities.

Communities are inferred from the researcher-level weighted graph. This module does
not claim that a graph community is an administrative laboratory or that the hub
researcher is a PI. It reports a reproducible network community and its central hub.

Inputs (from skku_pubmed_researcher_network.py):
  researcher_network_nodes.json
  researcher_network_edges.json

Outputs:
  research_communities.csv / research_communities.json
  researcher_community_membership.csv / researcher_community_membership.json
  community_edges.csv / community_edges.json
  research_community_map.html
  research_community_summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict

try:
    import networkx as nx
except ImportError:  # pragma: no cover - surfaced with a clear runtime message.
    nx = None
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass
class CommunityMember:
    community_id: str
    researcher_key: str
    researcher_name: str
    orcid: str
    paper_count: int
    weighted_degree: float
    within_community_degree: float
    role: str


@dataclass
class ResearchCommunity:
    community_id: str
    label: str
    researcher_count: int
    researchers: list[str]
    hub_researcher: str
    hub_key: str
    total_papers: int
    first_year: int
    last_year: int
    top_diseases: list[str]
    top_methods: list[str]
    top_data_types: list[str]
    stage_path: list[str]
    internal_edge_count: int
    internal_weight: float
    strong_edge_count: int
    evidence_summary: str


@dataclass
class CommunityEdge:
    source: str
    target: str
    source_label: str
    target_label: str
    researcher_edge_count: int
    max_score: float
    mean_score: float
    collaboration_edges: int
    citation_edges: int
    thematic_edges: int
    evidence: str


def load_json(path: Path) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError(f"Expected a JSON list: {path}")
    return data


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _edge_weight(edge: dict) -> float:
    return float(edge.get("score", 0.0) or 0.0)


def _build_adjacency(
    node_keys: set[str],
    edges: list[dict],
    min_edge_score: float,
) -> dict[str, dict[str, float]]:
    adjacency: dict[str, dict[str, float]] = {key: {} for key in node_keys}
    for edge in edges:
        source = str(edge.get("source", ""))
        target = str(edge.get("target", ""))
        score = _edge_weight(edge)
        if source not in node_keys or target not in node_keys or source == target:
            continue
        if score < min_edge_score:
            continue
        adjacency[source][target] = max(score, adjacency[source].get(target, 0.0))
        adjacency[target][source] = max(score, adjacency[target].get(source, 0.0))
    return adjacency


def _connected_component_labels(
    node_keys: set[str],
    adjacency: dict[str, dict[str, float]],
) -> dict[str, str]:
    visited: set[str] = set()
    labels: dict[str, str] = {}
    for start in sorted(node_keys):
        if start in visited:
            continue
        stack = [start]
        component: list[str] = []
        visited.add(start)
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in sorted(adjacency.get(current, {})):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        canonical = min(component)
        for member in component:
            labels[member] = canonical
    return labels


def _build_networkx_graph(
    node_keys: set[str],
    edges: list[dict],
    min_edge_score: float,
):
    if nx is None:
        raise RuntimeError(
            "networkx is required for Louvain community detection. "
            "Install networkx>=3.2 or use --algorithm connected."
        )
    graph = nx.Graph()
    graph.add_nodes_from(sorted(node_keys))
    for edge in edges:
        source = str(edge.get("source", ""))
        target = str(edge.get("target", ""))
        score = _edge_weight(edge)
        if (
            source not in node_keys
            or target not in node_keys
            or source == target
            or score < min_edge_score
        ):
            continue
        if graph.has_edge(source, target):
            if score > graph[source][target].get("weight", 0.0):
                graph[source][target]["weight"] = score
        else:
            graph.add_edge(source, target, weight=score)
    return graph


def detect_communities(
    nodes: list[dict],
    edges: list[dict],
    min_edge_score: float = 0.65,
    algorithm: str = "louvain",
    resolution: float = 1.0,
    seed: int = 42,
) -> dict[str, str]:
    """Return researcher -> deterministic raw community label.

    Louvain is the default because connected components collapse institution-scale
    collaboration graphs into giant components through a handful of bridge edges.
    The connected mode is retained only as a diagnostic/backward-compatible option.
    """
    node_keys = {str(node.get("key", "")) for node in nodes if node.get("key")}
    adjacency = _build_adjacency(node_keys, edges, min_edge_score)

    if algorithm == "connected":
        return _connected_component_labels(node_keys, adjacency)
    if algorithm != "louvain":
        raise ValueError(f"Unsupported community algorithm: {algorithm}")

    graph = _build_networkx_graph(node_keys, edges, min_edge_score)
    communities = nx.algorithms.community.louvain_communities(
        graph,
        weight="weight",
        resolution=resolution,
        threshold=1e-7,
        seed=seed,
    )

    ordered = sorted(
        (sorted(group) for group in communities),
        key=lambda group: (-len(group), group[0] if group else ""),
    )
    labels: dict[str, str] = {}
    for group in ordered:
        if not group:
            continue
        canonical = group[0]
        for key in group:
            labels[key] = canonical

    # networkx includes isolated nodes, but keep this defensive fallback explicit.
    for key in sorted(node_keys):
        labels.setdefault(key, key)
    return labels


def calculate_modularity(
    nodes: list[dict],
    edges: list[dict],
    membership_rows: list["CommunityMember"],
    min_edge_score: float,
) -> float:
    if nx is None:
        return 0.0
    node_keys = {str(node.get("key", "")) for node in nodes if node.get("key")}
    graph = _build_networkx_graph(node_keys, edges, min_edge_score)
    if graph.number_of_edges() == 0:
        return 0.0
    grouped: dict[str, set[str]] = defaultdict(set)
    for member in membership_rows:
        grouped[member.community_id].add(member.researcher_key)
    partition = [group for group in grouped.values() if group]
    return round(
        nx.algorithms.community.modularity(graph, partition, weight="weight"),
        6,
    )


def _counter_top(values: list[list[str]], n: int = 6) -> list[str]:
    counter = Counter(item for group in values for item in (group or []) if item)
    return [item for item, _ in counter.most_common(n)]


def _community_label(
    diseases: list[str],
    methods: list[str],
    data_types: list[str],
    community_id: str,
) -> str:
    pieces = []
    if diseases:
        pieces.append(diseases[0])
    if methods:
        pieces.append(methods[0])
    elif data_types:
        pieces.append(data_types[0])
    return " / ".join(pieces) if pieces else f"Research community {community_id}"


def build_research_communities(
    nodes: list[dict],
    edges: list[dict],
    min_edge_score: float = 0.65,
    strong_edge_score: float = 0.85,
    algorithm: str = "louvain",
    resolution: float = 1.0,
    seed: int = 42,
) -> tuple[list[ResearchCommunity], list[CommunityMember], list[CommunityEdge]]:
    raw_labels = detect_communities(
        nodes,
        edges,
        min_edge_score=min_edge_score,
        algorithm=algorithm,
        resolution=resolution,
        seed=seed,
    )
    node_map = {str(node.get("key", "")): node for node in nodes if node.get("key")}

    grouped: dict[str, list[str]] = defaultdict(list)
    for key, raw in raw_labels.items():
        grouped[raw].append(key)

    # Assign compact, stable community IDs by largest group then lexical key.
    ordered_groups = sorted(
        grouped.items(),
        key=lambda item: (-len(item[1]), min(item[1])),
    )
    compact_id: dict[str, str] = {
        raw: f"C{idx:03d}" for idx, (raw, _) in enumerate(ordered_groups, 1)
    }
    membership = {
        key: compact_id[raw]
        for key, raw in raw_labels.items()
    }

    weighted_degree: dict[str, float] = defaultdict(float)
    within_degree: dict[str, float] = defaultdict(float)
    internal_edges_by_community: dict[str, list[dict]] = defaultdict(list)
    cross_edges: dict[tuple[str, str], list[dict]] = defaultdict(list)

    for edge in edges:
        source = str(edge.get("source", ""))
        target = str(edge.get("target", ""))
        if source not in membership or target not in membership:
            continue
        score = _edge_weight(edge)
        weighted_degree[source] += score
        weighted_degree[target] += score

        source_c = membership[source]
        target_c = membership[target]
        if source_c == target_c:
            within_degree[source] += score
            within_degree[target] += score
            internal_edges_by_community[source_c].append(edge)
        else:
            pair = tuple(sorted((source_c, target_c)))
            cross_edges[pair].append(edge)

    communities: list[ResearchCommunity] = []
    members: list[CommunityMember] = []

    for raw, keys in ordered_groups:
        cid = compact_id[raw]
        group_nodes = [node_map[key] for key in sorted(keys)]
        diseases = _counter_top([list(x.get("top_diseases", []) or []) for x in group_nodes])
        methods = _counter_top([list(x.get("top_methods", []) or []) for x in group_nodes])
        data_types = _counter_top([list(x.get("top_data_types", []) or []) for x in group_nodes])

        stage_counter = Counter(
            stage
            for node in group_nodes
            for stage in (node.get("stage_path", []) or [])
            if stage
        )
        stages = [stage for stage, _ in stage_counter.most_common(8)]

        years_first = [int(x.get("first_year", 0) or 0) for x in group_nodes if int(x.get("first_year", 0) or 0)]
        years_last = [int(x.get("last_year", 0) or 0) for x in group_nodes if int(x.get("last_year", 0) or 0)]

        hub_key = max(
            keys,
            key=lambda key: (
                within_degree[key],
                weighted_degree[key],
                int(node_map[key].get("strong_lineage_count", 0) or 0),
                int(node_map[key].get("paper_count", 0) or 0),
                str(node_map[key].get("name", "")),
            ),
        )
        hub_name = str(node_map[hub_key].get("name", ""))
        internal = internal_edges_by_community.get(cid, [])
        internal_weight = round(sum(_edge_weight(edge) for edge in internal), 4)
        strong_count = sum(_edge_weight(edge) >= strong_edge_score for edge in internal)
        label = _community_label(diseases, methods, data_types, cid)

        evidence_parts = [
            f"researchers={len(group_nodes)}",
            f"internal edges={len(internal)}",
            f"strong edges={strong_count}",
        ]
        if diseases:
            evidence_parts.append("disease=" + ", ".join(diseases[:3]))
        if methods:
            evidence_parts.append("method=" + ", ".join(methods[:3]))
        if data_types:
            evidence_parts.append("data=" + ", ".join(data_types[:3]))

        communities.append(
            ResearchCommunity(
                community_id=cid,
                label=label,
                researcher_count=len(group_nodes),
                researchers=sorted(str(x.get("name", "")) for x in group_nodes),
                hub_researcher=hub_name,
                hub_key=hub_key,
                total_papers=sum(int(x.get("paper_count", 0) or 0) for x in group_nodes),
                first_year=min(years_first) if years_first else 0,
                last_year=max(years_last) if years_last else 0,
                top_diseases=diseases,
                top_methods=methods,
                top_data_types=data_types,
                stage_path=stages,
                internal_edge_count=len(internal),
                internal_weight=internal_weight,
                strong_edge_count=strong_count,
                evidence_summary="; ".join(evidence_parts),
            )
        )

        for key in sorted(keys):
            node = node_map[key]
            members.append(
                CommunityMember(
                    community_id=cid,
                    researcher_key=key,
                    researcher_name=str(node.get("name", "")),
                    orcid=str(node.get("orcid", "")),
                    paper_count=int(node.get("paper_count", 0) or 0),
                    weighted_degree=round(weighted_degree[key], 4),
                    within_community_degree=round(within_degree[key], 4),
                    role="hub" if key == hub_key else "member",
                )
            )

    community_map = {x.community_id: x for x in communities}
    community_edges: list[CommunityEdge] = []
    for (source, target), group in sorted(cross_edges.items()):
        scores = [_edge_weight(edge) for edge in group]
        relation_counts = Counter(str(edge.get("relation", "")) for edge in group)
        community_edges.append(
            CommunityEdge(
                source=source,
                target=target,
                source_label=community_map[source].label,
                target_label=community_map[target].label,
                researcher_edge_count=len(group),
                max_score=round(max(scores), 4),
                mean_score=round(sum(scores) / len(scores), 4),
                collaboration_edges=relation_counts["collaboration"] + relation_counts["collaboration+citation"],
                citation_edges=relation_counts["citation"] + relation_counts["collaboration+citation"],
                thematic_edges=relation_counts["thematic_overlap"],
                evidence=(
                    f"researcher edges={len(group)}; max score={max(scores):.3f}; "
                    f"mean score={sum(scores) / len(scores):.3f}"
                ),
            )
        )

    communities.sort(key=lambda x: (-x.researcher_count, -x.internal_weight, x.community_id))
    members.sort(key=lambda x: (x.community_id, x.role != "hub", -x.within_community_degree, x.researcher_name))
    community_edges.sort(key=lambda x: (-x.max_score, -x.researcher_edge_count, x.source, x.target))
    return communities, members, community_edges


def render_html(
    path: Path,
    communities: list[ResearchCommunity],
    community_edges: list[CommunityEdge],
) -> None:
    vis_nodes = []
    for item in communities:
        tooltip = (
            f"<b>{item.label}</b><br>"
            f"Community: {item.community_id}<br>"
            f"Researchers: {item.researcher_count}<br>"
            f"Hub researcher: {item.hub_researcher}<br>"
            f"Years: {item.first_year or '?'}–{item.last_year or '?'}<br>"
            f"Disease: {', '.join(item.top_diseases[:5]) or '-'}<br>"
            f"Methods: {', '.join(item.top_methods[:5]) or '-'}<br>"
            f"Data: {', '.join(item.top_data_types[:5]) or '-'}<br>"
            f"Researchers: {', '.join(item.researchers[:10])}"
        )
        vis_nodes.append(
            {
                "id": item.community_id,
                "label": f"{item.community_id} | {item.label}\n{item.researcher_count} researchers",
                "title": tooltip,
                "value": max(1, item.researcher_count + item.strong_edge_count),
            }
        )

    vis_edges = [
        {
            "from": edge.source,
            "to": edge.target,
            "value": max(1, edge.max_score * 6),
            "title": (
                f"researcher edges={edge.researcher_edge_count}<br>"
                f"max score={edge.max_score:.2f}<br>"
                f"mean score={edge.mean_score:.2f}"
            ),
        }
        for edge in community_edges
    ]

    nodes_json = json.dumps(vis_nodes, ensure_ascii=False)
    edges_json = json.dumps(vis_edges, ensure_ascii=False)
    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>SKKU PubMed Research Communities</title>
<script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
<style>
body {{ margin:0; font-family:Arial,sans-serif; background:#fafafa; }}
header {{ padding:14px 18px; background:white; border-bottom:1px solid #ddd; }}
#network {{ height:82vh; background:white; }}
.note {{ color:#555; font-size:12px; margin-top:6px; }}
</style>
</head>
<body>
<header>
<h2 style="margin:0">SKKU-seeded PubMed Research Communities</h2>
<div class="note">Communities are graph-derived research groups. Hub researcher is a network-central researcher, not a confirmed PI/lab head.</div>
</header>
<div id="network"></div>
<script>
const nodes = new vis.DataSet({nodes_json});
const edges = new vis.DataSet({edges_json});
const network = new vis.Network(
  document.getElementById("network"),
  {{nodes, edges}},
  {{
    physics: {{stabilization:true, barnesHut:{{gravitationalConstant:-10000}}}},
    interaction: {{hover:true, navigationButtons:true}},
    nodes: {{shape:"box", margin:10, scaling:{{min:12,max:40}}}},
    edges: {{smooth:{{type:"dynamic"}}, scaling:{{min:1,max:8}}}}
  }}
);
</script>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Cluster SKKU researcher network into research communities.")
    parser.add_argument("--input-dir", default="outputs/skku_pubmed_followup")
    parser.add_argument("--output-dir", default="")
    parser.add_argument("--min-edge-score", type=float, default=0.65)
    parser.add_argument("--strong-edge-score", type=float, default=0.85)
    parser.add_argument(
        "--algorithm",
        choices=["louvain", "connected"],
        default="louvain",
        help="Community detection algorithm. Louvain is recommended at institution scale.",
    )
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_json(input_dir / "researcher_network_nodes.json")
    edges = load_json(input_dir / "researcher_network_edges.json")
    communities, members, community_edges = build_research_communities(
        nodes,
        edges,
        min_edge_score=args.min_edge_score,
        strong_edge_score=args.strong_edge_score,
        algorithm=args.algorithm,
        resolution=args.resolution,
        seed=args.seed,
    )

    community_json = [asdict(x) for x in communities]
    member_json = [asdict(x) for x in members]
    edge_json = [asdict(x) for x in community_edges]

    (output_dir / "research_communities.json").write_text(
        json.dumps(community_json, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (output_dir / "researcher_community_membership.json").write_text(
        json.dumps(member_json, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (output_dir / "community_edges.json").write_text(
        json.dumps(edge_json, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    community_rows = []
    for item in communities:
        row = asdict(item)
        for key in ["researchers", "top_diseases", "top_methods", "top_data_types", "stage_path"]:
            row[key] = "; ".join(row[key])
        community_rows.append(row)
    write_csv(output_dir / "research_communities.csv", community_rows)
    write_csv(output_dir / "researcher_community_membership.csv", member_json)
    write_csv(output_dir / "community_edges.csv", edge_json)
    render_html(output_dir / "research_community_map.html", communities, community_edges)

    largest_size = max((x.researcher_count for x in communities), default=0)
    summary = {
        "community_version": 2,
        "algorithm": args.algorithm,
        "resolution": args.resolution,
        "seed": args.seed,
        "communities": len(communities),
        "researchers": len(nodes),
        "multi_researcher_communities": sum(x.researcher_count >= 2 for x in communities),
        "singleton_communities": sum(x.researcher_count == 1 for x in communities),
        "inter_community_edges": len(community_edges),
        "largest_community_size": largest_size,
        "largest_community_fraction": round(largest_size / len(nodes), 6) if nodes else 0.0,
        "modularity": calculate_modularity(
            nodes,
            edges,
            members,
            min_edge_score=args.min_edge_score,
        ),
        "min_edge_score": args.min_edge_score,
        "strong_edge_score": args.strong_edge_score,
    }
    (output_dir / "research_community_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
