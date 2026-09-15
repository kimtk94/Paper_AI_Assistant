#!/usr/bin/env python3
"""Deterministic PubMed research-profile annotation for SKKU lineage papers.

This module intentionally avoids LLM/API dependencies so the same annotations can
be reproduced in CI, Colab, and local runs. It uses title, abstract, MeSH terms,
keywords, and publication types already returned by PubMed.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from typing import Iterable

from skku_pubmed_lineage import Paper


@dataclass
class PaperAnnotation:
    pmid: str
    disease_terms: list[str]
    methods: list[str]
    data_types: list[str]
    research_stage: str
    research_question: str
    evidence_terms: list[str]


DISEASE_PATTERNS: dict[str, tuple[str, ...]] = {
    "cancer": ("cancer", "carcinoma", "neoplasm", "tumor", "tumour", "leukemia", "lymphoma", "melanoma"),
    "stroke": ("stroke", "cerebrovascular", "brain ischem", "cerebral infarct"),
    "cardiovascular disease": ("cardiovascular", "coronary", "myocardial infar", "heart failure", "atherosclero"),
    "diabetes": ("diabetes", "diabetic", "hyperglyc"),
    "obesity/metabolic disease": ("obesity", "metabolic syndrome", "dyslipid", "fatty liver", "steat"),
    "Alzheimer disease/dementia": ("alzheimer", "dementia", "cognitive impairment"),
    "Parkinson disease": ("parkinson",),
    "psychiatric disease": ("depression", "schizophrenia", "bipolar", "anxiety disorder", "psychiatric"),
    "kidney disease": ("kidney disease", "renal disease", "nephro", "chronic kidney"),
    "liver disease": ("liver disease", "hepatic disease", "hepatitis", "cirrhosis"),
    "respiratory disease": ("asthma", "copd", "pulmonary disease", "respiratory disease"),
    "infectious disease": ("infection", "infectious", "sepsis", "covid", "sars-cov", "influenza", "tuberculosis"),
    "autoimmune/inflammatory disease": ("autoimmune", "arthritis", "lupus", "inflammatory bowel", "crohn", "colitis"),
    "osteoporosis/bone disease": ("osteoporosis", "bone disease", "osteoarthritis"),
    "rare/genetic disease": ("rare disease", "genetic disorder", "congenital disorder"),
}

METHOD_PATTERNS: dict[str, tuple[str, ...]] = {
    "GWAS": ("genome-wide association", "gwas"),
    "Mendelian randomization": ("mendelian random",),
    "genetic association": ("genetic association", "polygenic risk", "prs", "variant association"),
    "WGS/WES": ("whole genome sequencing", "whole-genome sequencing", "whole exome", "whole-exome", "wgs", "wes"),
    "bulk RNA-seq": ("rna-seq", "rna sequencing", "transcriptome sequencing"),
    "single-cell RNA-seq": ("single-cell rna", "single cell rna", "scrna", "snrna", "single-nucleus rna"),
    "single-cell multiome": ("multiome", "single-cell atac", "single cell atac", "snatac"),
    "spatial transcriptomics": ("spatial transcript", "visium", "xenium", "stereo-seq"),
    "proteomics": ("proteomic", "protein profiling", "mass spectrometry proteom"),
    "metabolomics": ("metabolomic", "metabolite profiling"),
    "epigenomics": ("methylation", "epigenom", "atac-seq", "chip-seq"),
    "microbiome": ("microbiome", "microbiota", "16s rrna", "metagenom"),
    "imaging": ("mri", "magnetic resonance imaging", "ct imaging", "pet imaging", "radiomics"),
    "prospective cohort": ("prospective cohort",),
    "retrospective cohort": ("retrospective cohort",),
    "case-control": ("case-control", "case control"),
    "randomized trial": ("randomized controlled", "randomised controlled", "clinical trial"),
    "meta-analysis": ("meta-analysis", "meta analysis", "systematic review"),
    "machine learning": ("machine learning", "deep learning", "artificial intelligence", "neural network"),
    "survival analysis": ("survival analysis", "cox regression", "cox proportional"),
    "causal inference": ("causal inference", "instrumental variable", "target trial"),
}

DATA_PATTERNS: dict[str, tuple[str, ...]] = {
    "genomics": ("genome", "genomic", "gwas", "variant", "snp", "whole exome", "whole genome"),
    "transcriptomics": ("transcriptom", "rna-seq", "gene expression", "single-cell rna", "single cell rna"),
    "epigenomics": ("epigen", "methylation", "chromatin", "atac-seq", "chip-seq"),
    "proteomics": ("proteom", "protein profiling"),
    "metabolomics": ("metabolom", "metabolite"),
    "microbiome": ("microbiome", "microbiota", "metagenom", "16s rrna"),
    "single-cell": ("single-cell", "single cell", "single-nucleus", "single nucleus"),
    "spatial": ("spatial transcript", "visium", "xenium", "stereo-seq"),
    "clinical/EHR": ("electronic health", "ehr", "clinical data", "medical record", "patient record"),
    "imaging": ("mri", "magnetic resonance", "computed tomography", "radiomics", "imaging"),
    "questionnaire/phenotype": ("questionnaire", "survey", "phenotype", "phenotypic"),
}

STAGE_PATTERNS: list[tuple[str, tuple[str, ...]]] = [
    ("intervention", ("randomized controlled", "randomised controlled", "clinical trial", "intervention")),
    ("prediction", ("predict", "prediction", "prognostic", "risk score", "machine learning", "deep learning")),
    ("causal inference", ("mendelian random", "causal", "instrumental variable")),
    ("validation", ("validation", "replication", "validate", "external cohort")),
    ("mechanism", ("mechanism", "pathway", "functional", "cellular", "molecular mechanism")),
    ("discovery/association", ("association", "identify", "discover", "characterize", "profiling", "landscape")),
]

GENERIC_MESH = {
    "humans", "male", "female", "adult", "aged", "middle aged", "young adult",
    "animals", "mice", "rats", "retrospective studies", "prospective studies",
    "case-control studies", "cohort studies", "risk factors", "biomarkers",
}


def _norm(value: str) -> str:
    return re.sub(r"\s+", " ", value.lower().replace("‑", "-")).strip()


def _paper_text(paper: Paper) -> str:
    chunks = [
        paper.title,
        paper.abstract,
        " ".join(paper.mesh_terms),
        " ".join(paper.keywords),
        " ".join(paper.publication_types),
    ]
    return _norm(" ".join(x for x in chunks if x))


def _match_labels(text: str, patterns: dict[str, tuple[str, ...]]) -> tuple[list[str], list[str]]:
    labels: list[str] = []
    evidence: list[str] = []
    for label, needles in patterns.items():
        hits = [needle for needle in needles if needle in text]
        if hits:
            labels.append(label)
            evidence.extend(hits[:2])
    return labels, evidence


def infer_disease_terms(paper: Paper, text: str) -> tuple[list[str], list[str]]:
    labels, evidence = _match_labels(text, DISEASE_PATTERNS)

    # Preserve informative MeSH disease/phenotype descriptors even when they are
    # outside the curated broad disease dictionary.
    mesh_candidates = []
    for term in paper.mesh_terms:
        norm = _norm(term)
        if not norm or norm in GENERIC_MESH:
            continue
        if any(method in norm for method in (
            "sequencing", "analysis", "study", "model", "algorithm", "machine learning",
            "gene expression", "genomics", "proteomics", "metabolomics",
        )):
            continue
        if any(token in norm for token in (
            "disease", "syndrome", "cancer", "carcinoma", "neoplasm", "stroke",
            "diabetes", "obesity", "infection", "infarction", "failure", "disorder",
            "arthritis", "dementia", "alzheimer", "parkinson", "asthma",
        )):
            mesh_candidates.append(term)

    combined = []
    for item in labels + mesh_candidates:
        if item not in combined:
            combined.append(item)
    return combined[:8], evidence


def infer_stage(text: str, methods: Iterable[str]) -> str:
    for stage, needles in STAGE_PATTERNS:
        if any(needle in text for needle in needles):
            return stage
    method_set = set(methods)
    if method_set & {"GWAS", "genetic association", "bulk RNA-seq", "single-cell RNA-seq",
                     "single-cell multiome", "spatial transcriptomics", "proteomics",
                     "metabolomics", "microbiome", "epigenomics"}:
        return "discovery/association"
    if "meta-analysis" in method_set:
        return "evidence synthesis"
    return "descriptive/observational"


def make_research_question(
    disease_terms: list[str],
    methods: list[str],
    data_types: list[str],
    stage: str,
) -> str:
    disease = ", ".join(disease_terms[:2]) if disease_terms else "the target phenotype"
    method = ", ".join(methods[:2]) if methods else "observational analysis"
    data = ", ".join(data_types[:2]) if data_types else "biomedical data"

    templates = {
        "intervention": f"Does an intervention improve outcomes related to {disease}?",
        "prediction": f"Can {data} predict risk or outcomes for {disease}?",
        "causal inference": f"Is there evidence for a causal relationship involving {disease} using {method}?",
        "validation": f"Can previously reported findings for {disease} be validated using {method}?",
        "mechanism": f"What biological mechanisms underlying {disease} are revealed by {method}?",
        "evidence synthesis": f"What is the combined evidence regarding {disease} across prior studies?",
        "discovery/association": f"Which features are associated with {disease} using {method} and {data}?",
        "descriptive/observational": f"How is {disease} characterized in {data}?",
    }
    return templates.get(stage, f"How is {disease} studied using {method}?")


def annotate_paper(paper: Paper) -> PaperAnnotation:
    text = _paper_text(paper)
    disease_terms, disease_evidence = infer_disease_terms(paper, text)
    methods, method_evidence = _match_labels(text, METHOD_PATTERNS)
    data_types, data_evidence = _match_labels(text, DATA_PATTERNS)
    stage = infer_stage(text, methods)
    question = make_research_question(disease_terms, methods, data_types, stage)

    evidence = []
    for item in disease_evidence + method_evidence + data_evidence:
        if item not in evidence:
            evidence.append(item)

    return PaperAnnotation(
        pmid=paper.pmid,
        disease_terms=disease_terms,
        methods=methods,
        data_types=data_types,
        research_stage=stage,
        research_question=question,
        evidence_terms=evidence[:20],
    )


def progression_summary(source: PaperAnnotation, target: PaperAnnotation) -> str:
    changes = []

    if source.research_stage != target.research_stage:
        changes.append(f"stage: {source.research_stage} → {target.research_stage}")

    source_methods = set(source.methods)
    target_methods = set(target.methods)
    new_methods = [x for x in target.methods if x not in source_methods]
    if new_methods:
        changes.append("new method: " + ", ".join(new_methods[:3]))

    source_data = set(source.data_types)
    new_data = [x for x in target.data_types if x not in source_data]
    if new_data:
        changes.append("new data: " + ", ".join(new_data[:3]))

    source_disease = set(source.disease_terms)
    target_disease = set(target.disease_terms)
    if source_disease and target_disease and source_disease != target_disease:
        new_disease = [x for x in target.disease_terms if x not in source_disease]
        if new_disease:
            changes.append("expanded phenotype: " + ", ".join(new_disease[:3]))

    if not changes:
        shared_methods = sorted(source_methods & target_methods)
        shared_data = sorted(set(source.data_types) & set(target.data_types))
        stable = []
        if shared_methods:
            stable.append("method=" + ", ".join(shared_methods[:2]))
        if shared_data:
            stable.append("data=" + ", ".join(shared_data[:2]))
        if stable:
            return "continuation with stable " + "; ".join(stable)
        return "research continuation; no major deterministic profile shift detected"

    return "; ".join(changes)


def annotation_dict(annotation: PaperAnnotation) -> dict:
    return asdict(annotation)
