"""Strict processed GEO count-table identity and numeric guards."""
import re
import numpy as np


def resolve_sample(col, groups, lookup):
    c = col.lower()
    hits = set()
    for gsm in groups:
        if re.search(r"(?<![A-Za-z0-9])"+re.escape(gsm.lower())+r"(?![0-9])", c):
            hits.add(gsm)
    for gsm, identifier in lookup.items():
        if re.search(r"(?<![0-9])"+re.escape(identifier)+r"(?![0-9])", c):
            hits.add(gsm)
    for word, idx in re.findall(r"(?:ems[\W_]*)?(none|untreated|visanne|dienogest)[\W_]*([1-4])(?!\d)",c):
        group = "dienogest" if word in ("visanne","dienogest") else "untreated"
        hits.update(gsm for gsm, pair in groups.items() if pair==(group,int(idx)))
    if len(hits)>1:
        raise ValueError("STOP: a column contains conflicting sample identities")
    if not hits:
        return None
    gsm=next(iter(hits))
    return gsm,*groups[gsm]


def classify_metric(col, groups, lookup):
    c=col.lower()
    # pandas may suffix duplicate headers with .1/.2. Classify them in the
    # same bucket so the audit rejects duplicate sample/count mappings.
    m=re.fullmatch(r"(gsm\d+|\d+)_read_count(?:\.\d+)?",c)
    if m and resolve_sample(col,groups,lookup):
        return "count"
    for metric in ("fpkm","rpkm","tpm"):
        if metric in c:
            return metric
    return "unspecified"


def validate_raw_counts(values):
    arr=np.asarray(values,dtype=float)
    if arr.ndim!=2 or not arr.size:
        raise ValueError("STOP: empty count matrix")
    if not np.isfinite(arr).all() or (arr<0).any():
        raise ValueError("STOP: nonfinite, missing or negative raw counts")
    if (arr!=np.floor(arr)).any():
        raise ValueError("STOP: fractional values are not raw integer counts")
    if (arr>1e12).any():
        raise ValueError("STOP: implausibly large raw count")
    if (arr.sum(axis=0)<=0).any():
        raise ValueError("STOP: zero-count sample")
