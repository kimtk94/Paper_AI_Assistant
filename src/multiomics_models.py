"""Pydantic models and OpenClaw-friendly JSON parser helpers for multi-omics paper/idea records."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, conint, confloat

DiseaseDomain = Literal["metabolic", "brain"]
DiseaseName = Literal[
    "type_2_diabetes",
    "obesity",
    "nafld",
    "metabolic_syndrome",
    "alzheimers_disease",
    "parkinsons_disease",
    "schizophrenia",
    "major_depressive_disorder",
]
Modality = Literal["GWAS", "TWAS", "PWAS", "scRNA_seq", "eQTL", "pQTL", "ATAC_seq"]


class Population(BaseModel):
    model_config = ConfigDict(extra="forbid")
    ancestry: list[str] = Field(default_factory=list)
    sample_size: conint(ge=0) | None = None
    sex_info: str | None = None
    age_range: str | None = None


class DataResource(BaseModel):
    model_config = ConfigDict(extra="forbid")
    name: str
    type: Literal["summary_stats", "single_cell", "qtl", "proteomics", "other"]
    access: Literal["open", "controlled", "unknown"]
    url: str | None = None


class Reproducibility(BaseModel):
    model_config = ConfigDict(extra="forbid")
    code_available: bool
    code_url: str | None = None
    data_available: bool
    workflow_clarity_score: confloat(ge=0, le=10) | None = None


class EvidenceStrength(BaseModel):
    model_config = ConfigDict(extra="forbid")
    overall_score: confloat(ge=0, le=10)
    notes: str | None = None


class PaperRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    paper_id: str
    title: str
    authors: list[str] = Field(default_factory=list)
    journal_or_venue: str | None = None
    year: conint(ge=1990, le=2100)
    disease_domain: DiseaseDomain
    disease_name: DiseaseName
    population: Population | None = None
    modalities: list[Modality] = Field(min_length=1)
    study_type: Literal["method", "association", "integrative", "benchmark", "review"]
    main_question: str
    methods: list[str] = Field(default_factory=list)
    key_findings: list[str] = Field(min_length=1)
    limitations: list[str] = Field(min_length=1)
    open_problems: list[str] = Field(min_length=1)
    data_resources: list[DataResource]
    reproducibility: Reproducibility
    evidence_strength: EvidenceStrength
    tags: list[str] = Field(default_factory=list)


class Rationale(BaseModel):
    model_config = ConfigDict(extra="forbid")
    gap_statement: str
    biological_plausibility: str
    novelty_statement: str


class RequiredDataset(BaseModel):
    model_config = ConfigDict(extra="forbid")
    name: str
    modality: str
    access: Literal["open", "controlled", "unknown"]
    url: str | None = None


class AnalysisStep(BaseModel):
    model_config = ConfigDict(extra="forbid")
    step: conint(ge=1)
    description: str
    tool_hint: str | None = None


class RiskAndFallback(BaseModel):
    model_config = ConfigDict(extra="forbid")
    main_risks: list[str] = Field(default_factory=list)
    fallback_plan: list[str] = Field(default_factory=list)


class Scores(BaseModel):
    model_config = ConfigDict(extra="forbid")
    novelty: confloat(ge=0, le=10)
    feasibility: confloat(ge=0, le=10)
    impact: confloat(ge=0, le=10)
    overall: confloat(ge=0, le=10)


class IdeaRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    idea_id: str
    title: str
    disease_domain: DiseaseDomain
    disease_name: DiseaseName
    hypothesis: str
    target_modalities: list[Modality] = Field(min_length=2)
    rationale: Rationale
    supporting_papers: list[str] = Field(min_length=1)
    required_datasets: list[RequiredDataset]
    analysis_plan: list[AnalysisStep] = Field(min_length=3)
    evaluation_metrics: list[str] = Field(min_length=2)
    risk_and_fallback: RiskAndFallback
    scores: Scores
    priority_rank: conint(ge=1)


class MultiomicsBundle(BaseModel):
    model_config = ConfigDict(extra="forbid")
    papers: list[PaperRecord] = Field(default_factory=list)
    ideas: list[IdeaRecord] = Field(default_factory=list)


def parse_bundle(payload: str | bytes | dict) -> MultiomicsBundle:
    """Parse and validate an OpenClaw input payload into typed records."""
    if isinstance(payload, (str, bytes)):
        raw = json.loads(payload)
    else:
        raw = payload
    return MultiomicsBundle.model_validate(raw)


def load_bundle(path: str | Path) -> MultiomicsBundle:
    """Load a JSON file and parse into a validated model bundle."""
    data = Path(path).read_text(encoding="utf-8")
    return parse_bundle(data)


def dump_bundle(bundle: MultiomicsBundle, *, indent: int = 2) -> str:
    """Serialize validated models for downstream OpenClaw tools."""
    return bundle.model_dump_json(indent=indent)


if __name__ == "__main__":
    sample_path = Path(__file__).resolve().parents[1] / "examples" / "multiomics_records.example.json"
    bundle = load_bundle(sample_path)
    print(f"Validated bundle: papers={len(bundle.papers)}, ideas={len(bundle.ideas)}")
