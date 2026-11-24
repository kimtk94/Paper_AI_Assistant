#!/usr/bin/env python3
"""GSE283755 public RNA-seq exploratory audit. No treatment-causality inference.
Downloads public GEO Excel, stores SHA256, identifies known GSM group mapping, and
performs only exploratory metrics if group/sample columns can be unambiguously found.
"""
from __future__ import annotations
import argparse, datetime as dt, hashlib, io, json, math, re, sys, time, urllib.request
from pathlib import Path

from pathlib import Path as _Path
sys.path.insert(0, str(_Path(__file__).resolve().parent))
from count_evidence import resolve_sample, classify_metric, validate_raw_counts

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

ACCESSION = "GSE283755"
BASENAME = "GSE283755_Expression_Profile.GRCh38.gene.xlsx"
URLS = [
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE283755&format=file&file=" + BASENAME,
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283755/suppl/" + BASENAME,
]
# Confirmed from the NCBI GEO GSM entries, no inferred matching.
SAMPLE_GROUPS = {
    "GSM8671652": ("untreated", 1), "GSM8671653": ("dienogest", 1),
    "GSM8671654": ("untreated", 2), "GSM8671655": ("dienogest", 2),
    "GSM8671656": ("dienogest", 3), "GSM8671657": ("untreated", 3),
    "GSM8671658": ("dienogest", 4), "GSM8671659": ("untreated", 4),
}
# Numeric identifiers appear in GSE283755 processed XLSX; cross-check against each GEO GSM Description.
SAMPLE_ID_LOOKUP = {
    "GSM8671652": "10250742", "GSM8671653": "10532202",
    "GSM8671654": "10638608", "GSM8671655": "10925077",
    "GSM8671656": "11019271", "GSM8671657": "2449152",
    "GSM8671658": "7139439", "GSM8671659": "8127689",
}



def download(dest: Path, manifest: dict) -> None:
    if dest.exists() and dest.stat().st_size > 2_000_000:
        manifest["download"] = {"mode": "existing", "bytes": dest.stat().st_size}
        return
    errors = []
    for url in URLS:
        for attempt in range(2):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Research-QC/1.0 (public GEO study reanalysis)"})
                with urllib.request.urlopen(req, timeout=90) as resp:
                    length = int(resp.headers.get("Content-Length", "0") or 0)
                    if length > 100_000_000:
                        raise ValueError("Refusing oversized GEO supplementary file")
                    n = 0
                    with open(dest, "wb") as out:
                        while True:
                            chunk = resp.read(1024 * 1024)
                            if not chunk:
                                break
                            n += len(chunk)
                            if n > 100_000_000:
                                raise ValueError("Refusing oversized GEO file")
                            out.write(chunk)
                if n < 2_000_000 or not dest.read_bytes()[:2] == b"PK":
                    raise ValueError("Expected xlsx ZIP header and >2MB, got incompatible payload")
                manifest["download"] = {"mode": "remote", "source_url": url, "bytes": n}
                return
            except Exception as e:
                errors.append(f"{url}: attempt={attempt+1}: {type(e).__name__}: {e}")
                if dest.exists(): dest.unlink()
                time.sleep(2)
    raise RuntimeError("Unable to download GEO file: " + " | ".join(errors))


def canonical_sample(col: str):
    return resolve_sample(str(col), SAMPLE_GROUPS, SAMPLE_ID_LOOKUP)


def kind_of(col):
    return classify_metric(str(col), SAMPLE_GROUPS, SAMPLE_ID_LOOKUP)


def find_gene_col(df):
    candidates = []
    for col in df.columns:
        c = str(col).lower().replace(" ", "_")
        if re.search(r"gene.?symbol|symbol|gene.?name|hgnc", c):
            candidates.append((0, col))
        elif c in ("gene", "gene_id", "ensembl_gene_id", "geneid", "id"):
            candidates.append((1, col))
    return sorted(candidates, key=lambda x: x[0])[0][1] if candidates else None


def audit(input_path: Path, dest: Path, manifest: dict):
    excel = pd.ExcelFile(input_path, engine="openpyxl")
    if len(SAMPLE_ID_LOOKUP) != 8 or len(set(SAMPLE_ID_LOOKUP.values())) != 8:
        raise ValueError("Sample ID mapping is not one-to-one")
    manifest["sample_map"] = [{"gsm": gsm, "internal_id": SAMPLE_ID_LOOKUP[gsm], "group": group, "replicate": rep} for gsm,(group,rep) in SAMPLE_GROUPS.items()]
    candidates = []
    inventory = []
    for sheet in excel.sheet_names:
        df = pd.read_excel(excel, sheet_name=sheet, nrows=8)
        cols = [str(c) for c in df.columns]
        rec = {"sheet": sheet, "columns": cols, "first_8_rows": len(df), "genes_column": str(find_gene_col(df)), "matched": []}
        sets = {}
        for c in cols:
            mat = canonical_sample(c)
            if mat:
                gsm, group, idx = mat
                k = kind_of(c)
                group_cols = sets.setdefault(k, {})
                # Never let a second column silently replace the first GSM match.
                # Multiple measurements for one sample/metric require manual review.
                if gsm in group_cols:
                    rec.setdefault("ambiguous_matches", []).append(
                        {"gsm": gsm, "kind": k, "columns": [group_cols[gsm], c]}
                    )
                else:
                    group_cols[gsm] = c
                rec["matched"].append({"column": c, "gsm": gsm, "group": group, "kind": k})
        inventory.append(rec)
        for kind, colmap in sets.items():
            if len(colmap) == 8 and not any(
                a["kind"] == kind for a in rec.get("ambiguous_matches", [])
            ):
                candidates.append((sheet, kind, colmap))
    (dest / "sheet_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2), encoding="utf8")
    # Audit explicitly records ambiguity rather than silently picking arbitrary sample columns.
    if not candidates:
        ambiguous = [
            {"sheet": row["sheet"], "matches": row["ambiguous_matches"]}
            for row in inventory if row.get("ambiguous_matches")
        ]
        manifest["analysis"] = {"status": "BLOCKED_AMBIGUOUS_SAMPLE_MAPPING" if ambiguous else "BLOCKED_MATRIX_MAPPING", "reason": "Duplicate columns match the same GSM/metric; inspect inventory before analysis" if ambiguous else "No sheet contained 8 matched GSM/title columns for a common numeric metric", "candidate_count": 0, "ambiguous_sheets": ambiguous}
        return
    priority = {"count": 0, "tpm": 1, "fpkm": 2, "rpkm": 3, "unspecified": 4}
    candidates.sort(key=lambda x: priority.get(x[1], 10))
    best = priority.get(candidates[0][1], 10)
    if sum(priority.get(c[1],10)==best for c in candidates)>1:
        manifest["analysis"]={"status":"BLOCKED_AMBIGUOUS_SAMPLE_MAPPING",
                              "reason":"Multiple equally preferred complete matrices; no arbitrary sheet choice"}
        return
    sheet, kind, colmap = candidates[0]
    manifest["analysis"] = {"status": "PROFILE_COMPLETE", "all_candidates": [{"sheet": s, "metric": k} for s,k,_ in candidates], "chosen_sheet": sheet, "chosen_metric": kind}
    df = pd.read_excel(excel, sheet_name=sheet)
    gene_col = find_gene_col(df)
    if gene_col is None:
        manifest["analysis"].update(status="BLOCKED_GENE_IDENTIFIER", reason="Gene identifier column not recognized")
        return
    order = list(SAMPLE_GROUPS)
    selected = [colmap[gsm] for gsm in order]
    numbers = df[selected].apply(pd.to_numeric, errors="coerce")
    if kind == "count":
        if any(not re.fullmatch(r"(?:GSM\d+|\d+)_Read_Count", c, re.I) for c in selected):
            raise ValueError("STOP: raw count columns require exact published Read_Count headers")
        validate_raw_counts(numbers.to_numpy())
    valid = numbers.notna().all(axis=1) & (numbers >= 0).all(axis=1) & np.isfinite(numbers).all(axis=1)
    numbers = numbers.loc[valid]
    identifiers = df.loc[valid, gene_col].astype(str).str.strip()
    nonempty = identifiers.ne("") & identifiers.ne("nan")
    numbers = numbers.loc[nonempty]
    identifiers = identifiers.loc[nonempty]
    numbers.columns = order
    # Gene symbol duplicates: sum counts, instead of averaging, for count data; leave non-count metrics averaged.
    manifest["analysis"].update(source_rows=int(len(df)), invalid_numeric_rows=int((~valid).sum()),
                                missing_gene_identifier_rows=int((~nonempty).sum()),
                                gene_identifier_column=str(gene_col),
                                duplicate_identifier_rows=int(identifiers.duplicated().sum()))
    if kind == "count":
        numbers = numbers.groupby(identifiers.values, sort=False).sum()
        validate_raw_counts(numbers.to_numpy())
        numbers = numbers.astype("int64")
    else:
        numbers = numbers.groupby(identifiers.values, sort=False).mean()
    if len(numbers) < 1000:
        manifest["analysis"].update(status="BLOCKED_TOO_FEW_GENES", n_genes=int(len(numbers)))
        return
    props = []
    for gsm in order:
        group, idx = SAMPLE_GROUPS[gsm]
        s = numbers[gsm]
        props.append({"gsm": gsm, "group": group, "replicate": idx, "metric": kind, "sum": float(s.sum()), "median": float(s.median()), "zero_fraction": float((s == 0).mean()), "n_features": int(s.size)})
    pd.DataFrame(props).to_csv(dest / "sample_qc.csv", index=False)
    if kind == "count":
        numbers.to_csv(dest / "raw_counts_by_gene_symbol.csv", index_label="gene")
    # For count-like values use library-size normalized logCPM; otherwise log2(value+1).
    # This is NOT DESeq2 or confirmatory differential expression.
    if kind == "count":
        log = np.log2((numbers + .5) / (numbers.sum(axis=0) + 1) * 1_000_000 + 1)
    else:
        log = np.log2(numbers + 1)
    groups = [SAMPLE_GROUPS[x][0] for x in order]
    untreated = [x for x in order if SAMPLE_GROUPS[x][0] == "untreated"]
    treated = [x for x in order if SAMPLE_GROUPS[x][0] == "dienogest"]
    x, y = log[untreated].to_numpy(), log[treated].to_numpy()
    diff = y.mean(axis=1) - x.mean(axis=1)
    _, pval = ttest_ind(y, x, axis=1, equal_var=False, nan_policy="omit")
    validp = np.isfinite(pval)
    qvals = np.full(len(pval), np.nan)
    if validp.any(): qvals[validp] = multipletests(pval[validp], method="fdr_bh")[1]
    effects = pd.DataFrame({"gene": numbers.index, "delta_mean_log2": diff, "exploratory_welch_p": pval, "exploratory_bh_q": qvals})
    effects = effects.sort_values(["exploratory_bh_q", "exploratory_welch_p"], na_position="last")
    effects.to_csv(dest / "exploratory_gene_effects.csv", index=False)
    # Standard SVD PCA from sample-level expression. This cannot account for case mix or confounding.
    v = log.to_numpy().T
    vv = v - v.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(vv, full_matrices=False)
    variance = s**2
    variance = variance / max(float(variance.sum()), 1e-12)
    pd.DataFrame({"gsm": order, "group": groups, "pc1": (u*s)[:,0], "pc2": (u*s)[:,1], "pc1_variance":variance[0], "pc2_variance":variance[1]}).to_csv(dest / "pca_samples.csv", index=False)
    checks=[]
    for gene in ("LRRC15", "MMP3", "CXCL8", "IL6", "CCL2", "COL1A1"):
        hit=effects[effects.gene.str.upper().eq(gene)]
        checks.append({"gene":gene,"found":len(hit)>0,"delta_mean_log2":float(hit.iloc[0].delta_mean_log2) if len(hit) else None,"exploratory_q":float(hit.iloc[0].exploratory_bh_q) if len(hit) and np.isfinite(hit.iloc[0].exploratory_bh_q) else None})
    pd.DataFrame(checks).to_csv(dest / "candidate_direction_checks.csv", index=False)
    manifest["analysis"].update(status="EXPLORATORY_COMPLETE", n_genes=int(len(numbers)), n_samples=8, n_untreated=4, n_dienogest=4, n_exploratory_q_lt_0_05=int(np.nansum(qvals < .05)), caveat="Post-surgery sample cross-sectional association, not pre-treatment predictor. Welch test NOT DESeq2.")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output", default="output/endometriosis_gse283755")
    ap.add_argument("--input-xlsx", default="")
    args=ap.parse_args()
    out=Path(args.output); out.mkdir(parents=True, exist_ok=True)
    manifest={"accession":ACCESSION,"doi":"10.3349/ymj.2024.0518", "fetched_at_utc":dt.datetime.now(dt.timezone.utc).isoformat(), "study_design":"4 untreated vs 4 pre-operative dienogest-exposed ovarian endometriosis tissue; observational and unpaired", "download":None, "analysis":None, "constraints":["DO NOT treat treated and untreated as randomized", "DO NOT infer pretreatment response biomarker from post-treatment RNA", "GSE266265 is eutopic endometrium of untreated participants, not the ovarian endometrioma dataset", "No patient-level treatment response outcome in GSE283755"]}
    try:
        path=Path(args.input_xlsx) if args.input_xlsx else (out/BASENAME)
        if args.input_xlsx:
            manifest["download"]={"mode":"user_local_input","bytes":path.stat().st_size}
        else:
            download(path,manifest)
        h=hashlib.sha256()
        with open(path,"rb") as fh:
            for chunk in iter(lambda:fh.read(1024*1024),b""):h.update(chunk)
        manifest["sha256"]=h.hexdigest()
        audit(path,out,manifest)
    except Exception as e:
        manifest["analysis"]={"status":"FAILED", "error_type":type(e).__name__, "error":str(e)}
        print(f"[ERROR] {e}",file=sys.stderr)
    (out/"audit.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf8")
    with open(out/"REPORT.md","w",encoding="utf8") as fh:
        a=manifest["analysis"] or {}
        fh.write("# GSE283755 first-pass public data audit\n\n")
        fh.write("Status: **"+a.get("status","NOT_STARTED")+"**\n\n")
        fh.write("Study: Seok Kyo Seo et al. *Yonsei Med J* (2025), DOI 10.3349/ymj.2024.0518.\n\n")
        fh.write("Groups: untreated n=4; pre-operative dienogest-exposed n=4, different patients, ovarian endometrioma tissues.\n\n")
        fh.write("Published result: 406 reported DEGs; this pipeline does not presume reproduction.\n\n")
        fh.write("Analysis: "+json.dumps(a,ensure_ascii=False)+"\n\n")
        fh.write("Limits: observational confounding, low n, no pretreatment-response phenotype, FDR robustness uncertain. "
                 "Exploratory Welch q-values are NOT a DESeq2 reproduction and should not be used for definitive claims.\n\n")
        fh.write("Next gate: verify original raw counts and gene IDs, DESeq2 reproduction with design ~ treatment only, "
                 "independent cell-source contextualization and planned wet-lab validation.\n")
    print(json.dumps(manifest,ensure_ascii=False,indent=2))
    if manifest["analysis"]["status"] != "EXPLORATORY_COMPLETE" or manifest["analysis"].get("chosen_metric") != "count":
        raise SystemExit(1)

if __name__=="__main__":main()