#!/usr/bin/env python3
"""Audit original public ISA/MAF text, not the API's annotation-only projection.

Bounded metadata reads only. Never download instrument files or emit sample-level
measurements. Numeric columns alone do not establish phenotype, pairing or units.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import html.parser
import io
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

BASE="https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS13630/"
ARTICLE_ARCHIVE="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13098509/supplementaryFiles"
MODES=("GC-MS___metabolite_profiling","LC-MS_positive_reverse-phase_metabolite_profiling",
       "LC-MS_negative_reverse-phase_metabolite_profiling")
MISSING={"","na","n/a","nan","null","nd","not detected"}
TRANSIENT={429,500,502,503,504}


def fetch(url, cap=8_000_000):
    attempts=[]
    for attempt in range(2):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"Public-MetaboLights-Source-Audit/1.0"})
            with urllib.request.urlopen(req,timeout=25) as res:
                data=res.read(cap+1)
                if len(data)>cap:
                    return {"status":"TOO_LARGE","url":url,"limit":cap},None
                return {"status":"OK","url":url,"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest(),
                        "http":res.status,"attempts":attempts},data
        except urllib.error.HTTPError as e:
            attempts.append({"http":e.code})
            if e.code not in TRANSIENT:break
        except (OSError,TimeoutError,urllib.error.URLError) as e:
            attempts.append({"error_type":type(e).__name__})
        if attempt==0:time.sleep(2)
    return {"status":"UNAVAILABLE","url":url,"attempts":attempts},None


def tsv(data):
    text=data.decode("utf-8-sig",errors="strict")
    if text.lstrip().startswith(("<","{","[")):
        raise ValueError("STOP: not an original TSV table")
    rows=list(csv.reader(io.StringIO(text),delimiter="\t",strict=True))
    if len(rows)<2 or len(rows[0])<2:
        raise ValueError("STOP: empty or malformed TSV")
    if any(len(r)!=len(rows[0]) for r in rows[1:]):
        raise ValueError("STOP: ragged TSV")
    return [x.strip() for x in rows[0]],rows[1:]


def field(header,name):
    hits=[i for i,x in enumerate(header) if x==name]
    if len(hits)!=1:raise ValueError("STOP: missing or ambiguous identity field: "+name)
    return hits[0]


def sample_scope(data):
    header,rows=tsv(data)
    sample=field(header,"Sample Name")
    source=field(header,"Source Name")
    tissue=field(header,"Characteristics[Sample type]")
    variant=field(header,"Characteristics[Variant]")
    ids=[r[sample].strip() for r in rows]
    if not all(ids) or len(ids)!=len(set(ids)) or any(not r[source].strip() for r in rows):
        raise ValueError("STOP: blank or duplicate sample identity")
    # A Source Name denotes a specimen source, not a verified independent subject.
    groups=collections.Counter(r[tissue].strip().lower() for r in rows)
    phenotype_fields=[(i,h) for i,h in enumerate(header) if re.search(r"(?:variant|diagnosis|disease|clinical group|phenotype)",h,re.I)]
    summary={"sample_rows":len(rows),"unique_sample_names":len(ids),
             "unique_source_names":len({r[source] for r in rows}),
             "sample_type_counts":dict(groups),"blank_variant_rows":sum(not r[variant].strip() for r in rows),
             "explicit_phenotype_field_nonempty_counts":{h:sum(bool(r[i].strip()) for r in rows) for i,h in phenotype_fields},
             "independent_participant_mapping_verified":False,"case_control_mapping_verified":False}
    return set(ids),summary


def assay_scope(data,samples):
    h,rows=tsv(data); pos=field(h,"Sample Name"); maf=field(h,"Metabolite Assignment File")
    ids=[r[pos].strip() for r in rows]
    if any(x not in samples for x in ids) or len(set(ids))!=len(ids) or not all(ids):
        raise ValueError("STOP: assay contains missing, duplicate or unlinked specimen identity")
    refs={r[maf].strip().rsplit("/",1)[-1] for r in rows}
    if len(refs)!=1 or "" in refs:
        raise ValueError("STOP: ambiguous MAF reference")
    return set(ids),next(iter(refs))


def maf_scope(data, samples, assay_samples):
    h,rows=tsv(data)
    # Direct MAF quantitative columns are exact ISA sample names. Do not guess
    # from position, sample prefixes, raw filenames or aggregate abundance fields.
    positions=[(i,x) for i,x in enumerate(h) if x in samples]
    labels=[x for i,x in positions]
    if len(labels)!=len(set(labels)) or not set(labels)<=assay_samples:
        raise ValueError("STOP: duplicated or unlinked exact MAF sample column")
    numeric=0; missing=0; invalid=0; feature_rows=0; coverage=collections.Counter()
    for row in rows:
        found=False
        for i,name in positions:
            raw=row[i].strip()
            if raw.lower() in MISSING:missing+=1;continue
            try:v=float(raw)
            except ValueError:invalid+=1;continue
            if not math.isfinite(v):invalid+=1;continue
            numeric+=1;coverage[name]+=1;found=True
        feature_rows+=found
    unexpected=[x for x in h if x not in samples and (x.lower().startswith("s_") or re.search(r"sample_abundance",x,re.I))]
    return {"annotation_rows":len(rows),"columns":len(h),"exact_sample_columns":len(labels),
            "assay_samples":len(assay_samples),"sample_columns_with_finite_numeric_values":len(coverage),
            "finite_numeric_sample_cells":numeric,"missing_sample_cells":missing,"invalid_sample_cells":invalid,
            "features_with_sample_values":feature_rows,"unmatched_sample_like_headers":len(unexpected),
            "status":"CANDIDATE_VALUES_NEED_UNITS_QC_AND_PHENOTYPE" if numeric and not invalid else
                     "HOLD_INVALID_OR_UNLINKED_VALUES" if invalid or unexpected else "NO_EXACT_LINKED_NUMERIC_VALUES"}


class Links(html.parser.HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag=="a":
            href=dict(attrs).get("href")
            if href:self.links.append(href)


def directory_scope(data):
    p=Links();p.feed(data.decode("utf-8",errors="strict"))
    # Count only immediate relative file/directory names. Do not follow external
    # URLs or parent directories, and never download any listed instrument files.
    names=set()
    for h in p.links:
        s=urllib.parse.unquote(h)
        if not re.fullmatch(r"[^/?#\\]+/?",s) or s.strip("/") in (".",".."):
            continue
        names.add(s)
    if not names:raise ValueError("STOP: no directory entries parsed")
    dirs=[n for n in names if n.endswith("/")]
    files=[n for n in names if not n.endswith("/")]
    processed=[n for n in files if re.search(r"\.(?:csv|tsv|txt|xlsx|xls|parquet|mzTab)(?:\.gz)?$",n,re.I)
               and not re.match(r"^[aims]_MTBLS13630_|^i_Investigation",n,re.I)]
    return {"immediate_files":len(files),"immediate_subdirectories":len(dirs),
            "archive_or_instrument_files":sum(bool(re.search(r"\.(?:zip|raw|mzml|mzxml|d|cdf|tar|gz)$",n,re.I)) for n in files),
            "potential_processed_file_count":len(processed),
            "potential_processed_file_extensions":dict(collections.Counter(n.rsplit(".",1)[-1].lower() for n in processed)),
            "listing_scope":"immediate entries only; contents and nested directories not verified"}


def supplement_scope(data):
    if not zipfile.is_zipfile(io.BytesIO(data)):
        raise ValueError("STOP: official supplemental response is not a ZIP")
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        members=[i for i in archive.infolist() if not i.is_dir()]
        if len(members)>200:raise ValueError("STOP: excessive supplement member count")
        return {"member_count":len(members),
                "members":[{"basename":i.filename.rsplit("/",1)[-1],"bytes":i.file_size} for i in members],
                "tabular_member_count":sum(bool(re.search(r"\.(xlsx|xls|csv|tsv|txt)$",i.filename,re.I)) for i in members),
                "contents_verified":False,"note":"Inventory only; no spreadsheet, PDF or raw data extracted."}


def evaluate():
    result={"study":"MTBLS13630","source":"original public TSV and directory metadata",
            "case_control_analysis_allowed":False,"clinical_prediction_validated":False,
            "participant_identity_verified":False,"assays":{},"access":{}}
    def read(key,url,fn,cap=8_000_000):
        prov,data=fetch(url,cap);result["access"][key]=prov
        if data is None:return None
        try:return fn(data)
        except (ValueError,UnicodeError,csv.Error,zipfile.BadZipFile) as e:
            result["access"][key].update(status="PARSE_ERROR",error=str(e)[:180])
            return None
    sample=read("sample_tsv",BASE+"s_MTBLS13630.txt",sample_scope)
    if sample:
        ids,result["sample_scope"]=sample
        for mode in MODES:
            assay="a_MTBLS13630_"+mode+".txt"
            maf="m_MTBLS13630_"+mode+"_v2_maf.tsv"
            scope=read(assay,BASE+assay,lambda b:assay_scope(b,ids))
            if not scope:continue
            assay_ids,ref=scope
            if ref!=maf:
                result["assays"][mode]={"status":"HOLD_MAF_REFERENCE_MISMATCH"}
                continue
            report=read(maf,BASE+maf,lambda b:maf_scope(b,ids,assay_ids))
            if report:result["assays"][mode]=report
    result["root_listing"]=read("root_listing",BASE,directory_scope,2_000_000)
    result["files_listing"]=read("files_listing",BASE+"FILES/",directory_scope,2_000_000)
    result["article_supplement"]=read("article_supplement",ARTICLE_ARCHIVE,supplement_scope,25_000_000)
    a=result["assays"]
    if any(x.get("finite_numeric_sample_cells",0)>0 for x in a.values()):
        decision="HOLD_CANDIDATE_NUMERIC_VALUES_NEED_SOURCE_UNITS_QC_AND_PHENOTYPE"
    elif len(a)==3 and all(x.get("status")=="NO_EXACT_LINKED_NUMERIC_VALUES" for x in a.values()):
        decision="HOLD_NO_EXACT_LINKED_NUMERIC_VALUES_IN_THREE_DIRECT_MAFS"
    else:
        decision="HOLD_DIRECT_SOURCE_ACCESS_OR_SCHEMA_UNVERIFIED"
    result["decision"]=decision
    result["limitations"]=["No instrument raw archive fetched or processed; no quantitative matrix outside the three exact MAFs inspected.",
        "Directory/supplement listing alone does not establish file contents or completeness.",
        "PFO-associated plasma is not paired longitudinal Yonsei serum SCFA data.",
        "Independent participant, clinical phenotype, specimen pairing and measurement units require explicit authoritative mapping."]
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    r=evaluate();a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(r,ensure_ascii=False,allow_nan=False))
