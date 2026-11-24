#!/usr/bin/env python3
"""Fail-closed, offline check for a proposed 2023 HS GeoMx 24-ROI count matrix.

A file with correct columns and numeric values is NOT proof of its scientific
origin or compatibility with the 2026 study. Do not publish candidate matrices.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

SOURCE_DOI = '10.1093/bjd/ljad274'
EXPECTED_BY_DONOR = {'HS41': 6, 'HS42': 6, 'HS43': 12}
EXPECTED_ROIS = 24
MIN_GENES = 1000
MAX_GENES = 100000


def read_expected_manifest(path: Path) -> set[str]:
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        required = {'ROI_composite_key', 'patient_id', 'author_ROI_label',
                    'source_paper_doi', 'expression_value_available'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('STOP: official S2-derived ROI manifest schema incomplete')
        rows = list(reader)
    keys = [r['ROI_composite_key'].strip() for r in rows]
    donors = Counter(r['patient_id'].strip() for r in rows)
    if len(rows) != EXPECTED_ROIS or len(set(keys)) != EXPECTED_ROIS or donors != EXPECTED_BY_DONOR:
        raise ValueError('STOP: historical S2 manifest ROI/donor count or uniqueness changed')
    for r in rows:
        if (r['source_paper_doi'].strip().lower() != SOURCE_DOI or
            r['ROI_composite_key'].strip() != r['patient_id'].strip()+':'+r['author_ROI_label'].strip() or
            not r['author_ROI_label'].strip() or
            r['expression_value_available'].strip().upper() != 'NO'):
            raise ValueError('STOP: historical ROI source identity or unverified status changed')
    return set(keys)


def scan_candidate(path: Path, expected: set[str], *, minimum_genes: int = MIN_GENES) -> dict:
    """Validate raw integer count format only; never certify source authenticity."""
    extension = path.name.lower()
    if extension.endswith(('.csv.gz', '.csv')):
        separator = ','
    elif extension.endswith(('.tsv.gz', '.tsv')):
        separator = '\t'
    else:
        raise ValueError('STOP: only plain/gzip CSV or TSV raw integer count candidates are accepted')
    opener = gzip.open if extension.endswith('.gz') else open
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(data)
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f, delimiter=separator, strict=True)
        header = next(reader, [])
        if (len(header) != EXPECTED_ROIS+1 or
            header[0] not in {'gene_id', 'gene_symbol', 'target_name'} or
            len(set(header)) != len(header) or
            set(header[1:]) != expected):
            raise ValueError('STOP: candidate must contain exactly the 24 official HS ROI keys')
        library_totals = {roi: 0 for roi in header[1:]}
        detected = {roi: 0 for roi in header[1:]}
        seen_genes = set()
        genes = 0
        for row in reader:
            if len(row) != len(header):
                raise ValueError('STOP: ragged raw count matrix')
            gene = row[0].strip()
            if not gene or gene in seen_genes:
                raise ValueError('STOP: blank or duplicate gene identifier')
            seen_genes.add(gene)
            genes += 1
            if genes > MAX_GENES:
                raise ValueError('STOP: too many candidate features')
            for roi, value in zip(header[1:], row[1:]):
                v = value.strip()
                # No coercion from normalized floats, missing values or negatives.
                if not v.isascii() or not v.isdecimal():
                    raise ValueError('STOP: matrix includes missing, noninteger, or negative values')
                n = int(v)
                if n > 10**12:
                    raise ValueError('STOP: implausibly large raw count')
                library_totals[roi] += n
                detected[roi] += bool(n)
    if genes < minimum_genes:
        raise ValueError('STOP: candidate gene count too small')
    if any(v == 0 for v in library_totals.values()):
        raise ValueError('STOP: one or more ROI has zero counts')
    return {
        'candidate_input_sha256': sha.hexdigest(),
        'n_gene_features': genes,
        'n_roi_columns': EXPECTED_ROIS,
        'donor_roi_counts': dict(EXPECTED_BY_DONOR),
        'min_roi_library_total': min(library_totals.values()),
        'max_roi_library_total': max(library_totals.values()),
        'min_detected_gene_features': min(detected.values()),
        'max_detected_gene_features': max(detected.values()),
        'candidate_format': 'NONNEGATIVE_INTEGER_COUNT_MATRIX',
    }


def evaluate(manifest: Path, matrix: Path | None = None) -> dict:
    expected = read_expected_manifest(manifest)
    result = {
        'source_paper_doi': SOURCE_DOI,
        'expected_roi': len(expected),
        'expected_donors': len(EXPECTED_BY_DONOR),
        'scientific_replication': False,
        'source_authenticity_verified': False,
        'combined_156_roi_analysis_allowed': False,
    }
    if matrix is None:
        result['decision'] = 'HOLD_HISTORICAL_24_COUNT_MATRIX_NOT_SUPPLIED'
    else:
        result.update(scan_candidate(matrix, expected))
        result['decision'] = 'CANDIDATE_FORMAT_PASS__HOLD_INDEPENDENT_SOURCE_AND_METHOD_REVIEW'
    result['next_gate'] = (
        'Independently verify author/official repository provenance, accession, '
        '24 explicit ROI identities, dataset SHA-256, probe platform, count units, '
        'processing steps and cross-study normalization before any pooling.'
    )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--matrix', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    try:
        result = evaluate(args.manifest, args.matrix)
    except (ValueError, OSError, csv.Error, UnicodeError, EOFError) as exc:
        result = {'decision': 'ERROR_CANDIDATE_SOURCE_OR_FORMAT',
                  'error_type': type(exc).__name__, 'error': str(exc)[:350],
                  'scientific_replication': False,
                  'combined_156_roi_analysis_allowed': False}
        (args.out/'HISTORIC_24_QUANT_GATE.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result))
        raise SystemExit(1)
    (args.out/'HISTORIC_24_QUANT_GATE.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
