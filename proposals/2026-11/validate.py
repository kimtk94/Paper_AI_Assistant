#!/usr/bin/env python3
"""Validate the six research proposals and public-safe unsent email drafts."""
import argparse
import json
import re
from pathlib import Path

DIRECTORY = Path("proposals/2026-11")
PRIVATE_TERMS = re.compile(r"(?i)\b(?:SFDC|Salesforce|CRM)\b|api[_ -]?key|Bearer\s+\S+|password")
EMAIL_ADDRESS = re.compile(r"(?i)\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b")


def inspect_text(text: str, kind: str) -> list[str]:
    issues = []
    if len(text) < (650 if kind == "proposal" else 200):
        issues.append("text too short")
    if not re.search(r"^# .+", text, re.M):
        issues.append("top-level heading missing")
    if PRIVATE_TERMS.search(text) or EMAIL_ADDRESS.search(text):
        issues.append("potentially private text or email address")
    if kind == "proposal":
        for section in ("연구 배경", "연구 질문", "파일럿", "활용 가능한", "참고문헌"):
            if section not in text:
                issues.append(f"missing section: {section}")
        if not re.search(r"https://doi\.org/10\.\d{4,9}/\S+", text):
            issues.append("missing scholarly DOI")
        if "확정 연구계획·견적이 아닙니다" not in text:
            issues.append("missing scope disclaimer")
    elif "[담당자]" not in text or "첨부:" not in text:
        issues.append("missing public sender placeholder or attachment")
    return issues


def validate(root: Path) -> list[str]:
    problems = []
    base = root / DIRECTORY
    f = base / "manifest.json"
    try:
        entries = json.loads(f.read_text(encoding="utf-8"))["proposals"]
    except (OSError, ValueError, KeyError) as exc:
        return [f"manifest unavailable: {exc}"]
    if len(entries) != 6 or [e.get("id") for e in entries] != [f"{i:02d}" for i in range(1, 7)]:
        problems.append("manifest must have six ordered unique entries")
    for entry in entries:
        ident = entry.get("id", "?")
        for kind in ("proposal", "email"):
            path = entry.get(kind)
            if not isinstance(path, str) or not path.startswith(str(DIRECTORY) + "/") or ".." in path or not path.endswith(".md"):
                problems.append(f"{ident} invalid {kind} path")
                continue
            p = root / path
            if not p.is_file():
                problems.append(f"{ident} missing {path}")
                continue
            problems.extend(f"{ident}/{kind}: {issue}" for issue in inspect_text(p.read_text(encoding="utf-8"), kind))
        if not str(entry.get("source_doc", "")).startswith("https://docs.google.com/document/d/"):
            problems.append(f"{ident} invalid source_doc")
        if not str(entry.get("delivery_pdf", "")).startswith("https://drive.google.com/file/d/"):
            problems.append(f"{ident} invalid delivery_pdf")
        if not entry.get("doi") or any(not str(d).startswith("10.") for d in entry["doi"]):
            problems.append(f"{ident} DOI list missing or malformed")
    expected = {e["proposal"] for e in entries if "proposal" in e}
    actual = {str(p.relative_to(root)) for p in base.glob("*.md") if p.name not in {"README.md", "REFERENCES.md"}}
    if actual != expected:
        problems.append(f"unregistered or missing Markdown proposal files: {sorted(actual ^ expected)}")
    return problems


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path("."))
    args = p.parse_args()
    errors = validate(args.root.resolve())
    for error in errors:
        print("FAIL:", error)
    print("PASS: 6 public proposals, 6 unsent emails, DOI links and Drive provenance" if not errors else f"FAILED: {len(errors)} issues")
    raise SystemExit(0 if not errors else 2)
