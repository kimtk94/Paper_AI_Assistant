#!/usr/bin/env python3
"""Deterministic cross-domain classification for SKKU PubMed papers.

The SKKU corpus spans medicine, life science, engineering, chemistry, materials,
energy, and digital-health research. This module provides a coarse domain layer
before downstream disease/method annotation so non-biomedical papers are not forced
into clinical labels.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable

from skku_pubmed_lineage import Paper


@dataclass
class DomainProfile:
    primary_domain: str
    research_domains: list[str]
    topic_terms: list[str]
    evidence_terms: list[str]


DOMAIN_PATTERNS: dict[str, tuple[tuple[str, int], ...]] = {
    "Energy / Photovoltaics": (
        ("perovskite solar cell", 10),
        ("perovskite photovoltaic", 9),
        ("photovoltaic", 7),
        ("solar cell", 7),
        ("power conversion efficiency", 6),
        ("photoelectric conversion", 5),
        ("dye-sensitized", 6),
        ("tandem solar", 6),
    ),
    "Materials Science": (
        ("materials science", 8),
        ("perovskite", 5),
        ("thin film", 5),
        ("semiconductor", 5),
        ("nanomaterial", 5),
        ("nanostructure", 4),
        ("quantum dot", 6),
        ("memristor", 7),
        ("ferroelectric", 6),
        ("crystal structure", 4),
        ("heterostructure", 5),
        ("electrode material", 4),
        ("surface modification", 3),
    ),
    "Chemical Engineering": (
        ("chemical engineering", 9),
        ("electrocatal", 7),
        ("photocatal", 7),
        ("catalysis", 6),
        ("catalyst", 5),
        ("reaction engineering", 6),
        ("separation membrane", 6),
        ("membrane separation", 6),
        ("adsorption", 4),
        ("reactor", 4),
    ),
    "Chemistry": (
        ("organic synthesis", 7),
        ("chemical synthesis", 6),
        ("synthetic chemistry", 6),
        ("coordination complex", 5),
        ("organometallic", 6),
        ("ligand", 3),
        ("nmr spectroscopy", 5),
        ("raman spectroscopy", 4),
        ("x-ray diffraction", 4),
        ("spectroscopy", 3),
    ),
    "Digital Health / AI": (
        ("digital health", 9),
        ("machine learning", 7),
        ("deep learning", 7),
        ("artificial intelligence", 7),
        ("neural network", 6),
        ("large language model", 7),
        ("electronic health record", 6),
        ("ehr", 5),
        ("radiomics", 5),
        ("computer-aided diagnosis", 5),
    ),
    "Genomics / Omics": (
        ("genome-wide association", 8),
        ("gwas", 7),
        ("rna-seq", 7),
        ("single-cell rna", 8),
        ("single cell rna", 8),
        ("single-nucleus rna", 8),
        ("spatial transcript", 8),
        ("proteomic", 6),
        ("metabolomic", 6),
        ("epigenom", 6),
        ("methylation", 5),
        ("microbiome", 5),
        ("metagenom", 5),
        ("whole genome sequencing", 6),
        ("whole exome", 6),
    ),
    "Neuroscience": (
        ("neuroscience", 8),
        ("neuron", 5),
        ("neuronal", 5),
        ("brain", 4),
        ("synaptic plasticity", 6),
        ("neurodegener", 6),
        ("alzheimer", 6),
        ("parkinson", 6),
        ("stroke", 5),
        ("cerebrovascular", 5),
    ),
    "Public Health / Epidemiology": (
        ("public health", 8),
        ("epidemiolog", 7),
        ("population-based", 6),
        ("population based", 6),
        ("prevalence", 5),
        ("incidence", 5),
        ("cohort", 4),
        ("case-control", 4),
        ("risk factor", 4),
        ("registry", 3),
        ("mortality", 3),
    ),
    "Pharmaceutical Science": (
        ("pharmaceutical", 7),
        ("pharmacokinetic", 7),
        ("drug delivery", 7),
        ("drug formulation", 6),
        ("medicinal chemistry", 7),
        ("therapeutic drug monitoring", 6),
        ("drug concentration", 4),
        ("nanomedicine", 5),
    ),
    "Clinical Medicine": (
        ("clinical trial", 6),
        ("randomized controlled", 6),
        ("randomised controlled", 6),
        ("patient", 3),
        ("patients", 3),
        ("diagnosis", 4),
        ("treatment", 4),
        ("surgery", 4),
        ("therapy", 3),
        ("prognosis", 4),
        ("hospital", 2),
    ),
    "Biomedical / Life Science": (
        ("molecular biology", 6),
        ("cell signaling", 6),
        ("cell signalling", 6),
        ("protein expression", 5),
        ("gene expression", 5),
        ("animal model", 5),
        ("mouse model", 5),
        ("cell line", 4),
        ("pathway", 3),
        ("biomarker", 3),
        ("enzyme", 3),
        ("receptor", 3),
    ),
}

TOPIC_PATTERNS: dict[str, tuple[str, ...]] = {
    "perovskite solar cells": ("perovskite solar cell", "perovskite photovoltaic"),
    "photovoltaics": ("photovoltaic", "solar cell", "power conversion efficiency"),
    "memristors / neuromorphic devices": ("memristor", "neuromorphic", "artificial synapse"),
    "batteries / energy storage": ("battery", "lithium-ion", "lithium ion", "energy storage", "supercapacitor"),
    "catalysis": ("catalyst", "catalysis", "electrocatal", "photocatal"),
    "semiconductors / electronics": ("semiconductor", "transistor", "heterostructure", "electronic device"),
    "nanomaterials": ("nanomaterial", "nanostructure", "nanoparticle", "quantum dot"),
    "polymers": ("polymer", "copolymer"),
    "sensors": ("sensor", "biosensor"),
    "drug delivery": ("drug delivery", "nanomedicine", "formulation"),
    "cancer": ("cancer", "carcinoma", "tumor", "tumour", "neoplasm"),
    "stroke": ("stroke", "cerebrovascular", "brain ischem"),
    "single-cell omics": ("single-cell", "single cell", "single-nucleus", "single nucleus"),
    "genomics": ("genome-wide", "gwas", "whole genome", "whole exome", "genomic"),
    "AI / machine learning": ("machine learning", "deep learning", "artificial intelligence", "neural network"),
    "epidemiology": ("epidemiolog", "population-based", "prevalence", "incidence"),
    "neuroscience": ("neuron", "neuronal", "brain", "synaptic plasticity", "neurodegener"),
}

NON_BIOMEDICAL_DOMAINS = {
    "Energy / Photovoltaics",
    "Materials Science",
    "Chemical Engineering",
    "Chemistry",
}


def _norm(value: str) -> str:
    return re.sub(r"\s+", " ", (value or "").lower().replace("‑", "-")).strip()


def paper_text(paper: Paper) -> str:
    return _norm(
        " ".join(
            x
            for x in [
                paper.title,
                paper.abstract,
                " ".join(paper.mesh_terms),
                " ".join(paper.keywords),
                " ".join(paper.publication_types),
            ]
            if x
        )
    )


def _contains(text: str, needle: str) -> bool:
    if len(needle) <= 4 and needle.isalnum():
        return bool(re.search(rf"\b{re.escape(needle)}\b", text))
    return needle in text


def classify_domain(paper: Paper) -> DomainProfile:
    text = paper_text(paper)
    scores: dict[str, int] = {}
    evidence: dict[str, list[str]] = {}

    for domain, patterns in DOMAIN_PATTERNS.items():
        score = 0
        hits: list[str] = []
        for needle, weight in patterns:
            if _contains(text, needle):
                score += weight
                hits.append(needle)
        if score:
            scores[domain] = score
            evidence[domain] = hits

    # Avoid clinical/psychiatric semantics dominating engineering papers that use
    # terms such as artificial synapse or long-term depression.
    engineering_strength = max(
        (
            scores.get("Energy / Photovoltaics", 0),
            scores.get("Materials Science", 0),
            scores.get("Chemical Engineering", 0),
            scores.get("Chemistry", 0),
        ),
        default=0,
    )
    if engineering_strength >= 7:
        for biomedical in (
            "Clinical Medicine",
            "Public Health / Epidemiology",
            "Neuroscience",
            "Biomedical / Life Science",
        ):
            if scores.get(biomedical, 0) < engineering_strength:
                scores.pop(biomedical, None)
                evidence.pop(biomedical, None)

    ranked = sorted(scores, key=lambda d: (-scores[d], d))
    if not ranked:
        ranked = ["Other / Multidisciplinary"]
        evidence["Other / Multidisciplinary"] = []

    top_score = scores.get(ranked[0], 0)
    research_domains = [
        domain
        for domain in ranked
        if domain == ranked[0] or scores.get(domain, 0) >= max(4, int(top_score * 0.55))
    ][:4]

    topic_terms = [
        label
        for label, needles in TOPIC_PATTERNS.items()
        if any(_contains(text, needle) for needle in needles)
    ][:8]

    evidence_terms: list[str] = []
    for domain in research_domains:
        for term in evidence.get(domain, []):
            if term not in evidence_terms:
                evidence_terms.append(term)

    return DomainProfile(
        primary_domain=ranked[0],
        research_domains=research_domains,
        topic_terms=topic_terms,
        evidence_terms=evidence_terms[:20],
    )


def is_non_biomedical_domain(domain: str) -> bool:
    return domain in NON_BIOMEDICAL_DOMAINS
