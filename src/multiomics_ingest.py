#!/usr/bin/env python3
"""Batch ingest utility for multi-omics paper/idea bundles.

- Accepts a bundle JSON with top-level keys: `papers`, `ideas`
- Optionally validates with Pydantic models (if installed)
- Appends normalized records into JSONL stores for incremental accumulation
- Saves immutable timestamped snapshot for audit/replay
"""

from __future__ import annotations

import argparse
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from uuid import uuid4


def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _validate_bundle(bundle: Any, strict: bool) -> dict[str, Any]:
    """Validate via pydantic models when available.

    Envelope validation always runs. Only a missing optional Pydantic dependency
    permits a non-strict fallback; invalid records must never bypass validation.
    """
    if not isinstance(bundle, dict):
        raise ValueError("Bundle must be a JSON object with papers and/or ideas arrays.")
    unknown_keys = set(bundle) - {"papers", "ideas"}
    if not bundle or unknown_keys:
        raise ValueError(
            "Bundle must contain only papers and/or ideas arrays; "
            "paper_assistants reports with results are not ingest bundles."
        )
    for key, rows in bundle.items():
        if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
            raise ValueError(f"Bundle {key} must be an array of record objects.")

    try:
        from multiomics_models import parse_bundle  # local module
    except ModuleNotFoundError as exc:
        if exc.name != "pydantic":
            raise
        if strict:
            raise RuntimeError(
                "Strict validation requires Pydantic. Install dependencies first: "
                "pip install -r requirements-multiomics.txt"
            ) from exc
        warnings.warn(
            "Pydantic is not installed; only the bundle envelope was validated. "
            "Install requirements-multiomics.txt for record validation.",
            RuntimeWarning,
            stacklevel=2,
        )
        return bundle

    model = parse_bundle(bundle)
    return model.model_dump()


def _write_json_exclusive(path: Path, payload: dict[str, Any]) -> None:
    """Create a new audit file without replacing any existing run's evidence."""
    with path.open("x", encoding="utf-8") as f:
        f.write(json.dumps(payload, ensure_ascii=False, indent=2))


def _append_jsonl(path: Path, rows: list[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        for row in rows:
            f.write(json.dumps(row, ensure_ascii=False) + "\n")
    return len(rows)


def ingest_bundle(
    input_json: Path,
    data_dir: Path,
    *,
    strict_validation: bool = False,
) -> dict[str, Any]:
    bundle_raw = _load_json(input_json)
    bundle = _validate_bundle(bundle_raw, strict=strict_validation)

    ts = _utc_now()
    run_id = f"{ts}_{uuid4().hex}"
    snapshots_dir = data_dir / "snapshots"
    snapshots_dir.mkdir(parents=True, exist_ok=True)
    snapshot_path = snapshots_dir / f"bundle_{run_id}.json"
    _write_json_exclusive(snapshot_path, bundle)

    papers = bundle.get("papers", [])
    ideas = bundle.get("ideas", [])

    papers_rows = [{**p, "ingested_at": ts} for p in papers]
    ideas_rows = [{**i, "ingested_at": ts} for i in ideas]

    paper_n = _append_jsonl(data_dir / "paper_records.jsonl", papers_rows)
    idea_n = _append_jsonl(data_dir / "idea_records.jsonl", ideas_rows)

    manifest = {
        "run_id": run_id,
        "ingested_at": ts,
        "input": str(input_json),
        "snapshot": str(snapshot_path),
        "papers_appended": paper_n,
        "ideas_appended": idea_n,
    }

    manifests = data_dir / "manifests"
    manifests.mkdir(parents=True, exist_ok=True)
    _write_json_exclusive(manifests / f"manifest_{run_id}.json", manifest)

    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Ingest multi-omics bundle JSON into append-only JSONL stores")
    parser.add_argument("--input", required=True, help="Path to bundle JSON file")
    parser.add_argument("--data-dir", default="data/multiomics", help="Output data directory")
    parser.add_argument(
        "--strict-validation",
        action="store_true",
        help="Require pydantic validation success (fails if dependencies are missing)",
    )
    args = parser.parse_args()

    result = ingest_bundle(
        input_json=Path(args.input),
        data_dir=Path(args.data_dir),
        strict_validation=args.strict_validation,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
