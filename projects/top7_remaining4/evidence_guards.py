"""Fail-closed source and biological-unit guards for public HS/MG/NMOSD audits.
These checks establish table identity and cohort scope, not treatment prediction.
"""
import collections
import gzip
import io
import re
import tarfile

GSM=re.compile(r"GSM[0-9]+")
TRACK=re.compile(r"L[0-9]+_Track-[0-9]+")
MEMBER=re.compile(r"(GSM[0-9]+)_(.+)\.txt(?:\.gz)?",re.I)


def unique_keys(values, what):
    values=list(values)
    if any(not str(x).strip() for x in values) or len(set(values))!=len(values):
        raise ValueError("STOP: blank or duplicate "+what)
    return values


def library_key(name):
    # R make.names can prefix a numeric column with X and replace punctuation.
    value=str(name).strip().strip('"').lower()
    if re.match(r"^x[0-9]",value):
        value=value[1:]
    return re.sub(r"[.\-_]","",value)


def parse_raw_count_member(data):
    if len(data)>5_000_000:
        raise ValueError("STOP: individual count member too large")
    lines=data.decode("utf-8-sig").splitlines()
    if not lines:
        raise ValueError("STOP: empty raw count file")
    header=None
    if len(lines[0].strip().split())==1:
        header=lines.pop(0).strip().strip('"')
    counts={}
    for line in lines:
        fields=line.strip().split()
        if len(fields)!=2:
            raise ValueError("STOP: invalid raw count row")
        gene,n=(s.strip('"') for s in fields)
        if not gene.isascii() or not gene.isdecimal() or gene in counts or not n.isascii() or not n.isdecimal():
            raise ValueError("STOP: invalid or duplicated gene/count row")
        counts[gene]=int(n)
    if not counts:
        raise ValueError("STOP: no gene counts")
    return header,counts


def verify_archive_mapping(blob, columns, genes, vector_for_column, expected_gsms):
    """Verify all individual source vectors against the exact combined matrix.

    No substring match, positional match or partial signature is sufficient.
    vector_for_column returns the combined matrix values in the given gene order.
    """
    columns=unique_keys(columns,"matrix columns")
    genes=unique_keys(map(str,genes),"matrix gene identifiers")
    expected=set(unique_keys(expected_gsms,"GSM identifiers"))
    if len(columns)!=len(expected) or any(not GSM.fullmatch(s) for s in expected):
        raise ValueError("STOP: matrix/metadata sample identities incomplete")
    aliases={}
    for col in columns:
        key=library_key(col)
        if not key or key in aliases:
            raise ValueError("STOP: ambiguous normalized library identifier")
        aliases[key]=col
    mapping={}
    with tarfile.open(fileobj=io.BytesIO(blob),mode="r:*") as archive:
        files=[m for m in archive.getmembers() if m.isfile()]
        if len(files)!=len(columns):
            raise ValueError("STOP: raw archive member count mismatch")
        for m in files:
            match=MEMBER.fullmatch(m.name.rsplit("/",1)[-1])
            if not match or m.size>5_000_000:
                raise ValueError("STOP: unexpected individual source filename or size")
            gsm,stem=match.groups()
            col=aliases.get(library_key(stem))
            if gsm not in expected or col is None or col in mapping or gsm in mapping.values():
                raise ValueError("STOP: ambiguous or unknown raw library/GSM mapping")
            with archive.extractfile(m) as stream:
                if m.name.lower().endswith(".gz"):
                    with gzip.GzipFile(fileobj=stream) as gz:
                        data=gz.read(5_000_001)
                else:
                    data=stream.read(5_000_001)
            header,raw=parse_raw_count_member(data)
            if header is not None and library_key(header)!=library_key(col):
                raise ValueError("STOP: source library header conflicts with filename")
            if set(raw)!=set(genes):
                raise ValueError("STOP: individual/combined gene universe differs")
            combined=list(vector_for_column(col))
            if len(combined)!=len(genes) or any(raw[g]!=v for g,v in zip(genes,combined)):
                raise ValueError("STOP: full source count vector differs from combined matrix")
            mapping[col]=gsm
    if set(mapping)!=set(columns) or set(mapping.values())!=expected:
        raise ValueError("STOP: incomplete source mapping")
    return mapping


def paired_lesions(rows, expected_patients=14, expected_each=7):
    """Validate the biological pairing and response strata before emitting counts."""
    unique_keys([r["gsm"] for r in rows],"mapped GSMs")
    unique_keys([r["column"] for r in rows],"mapped columns")
    by={}
    for r in rows:
        if r["tissue"] not in {"lesional","nonlesional"} or r["visit"] not in {"V1","V12"} or not r["patient"]:
            raise ValueError("STOP: unexpected patient/tissue/visit")
        key=(r["patient"],r["tissue"],r["visit"])
        if key in by:
            raise ValueError("STOP: duplicate patient/tissue/visit")
        by[key]=r
    pairs=[(r,by[(p,t,"V12")]) for (p,t,v),r in by.items() if v=="V1" and (p,t,"V12") in by]
    lesion=[(a,b) for a,b in pairs if a["tissue"]=="lesional"]
    counts=collections.Counter()
    for a,b in lesion:
        if a["responder"] not in {"YES","NO"} or a["responder"]!=b["responder"]:
            raise ValueError("STOP: missing or inconsistent paired response label")
        counts[a["responder"]]+=1
    if len(lesion)!=expected_patients or counts!={"YES":expected_each,"NO":expected_each}:
        raise ValueError("STOP: unexpected paired patient/response group count")
    return pairs,lesion


def nmosd_track_mapping(samples, header):
    unique_keys([s["accession"] for s in samples],"NMOSD GSMs")
    tracks={}
    for s in samples:
        hits=re.findall(r"\[(L[0-9]+_Track-[0-9]+)\]",s["title"])
        if len(hits)!=1 or hits[0] in tracks:
            raise ValueError("STOP: missing or duplicate NMOSD Track identity")
        tracks[hits[0]]=s
    if header[:4]!=["","Gene_id","gene_biotype","description"]:
        raise ValueError("STOP: unexpected NMOSD matrix annotation header")
    cols=unique_keys(header[4:],"NMOSD Track columns")
    if any(not TRACK.fullmatch(col) or col not in tracks for col in cols):
        raise ValueError("STOP: unknown NMOSD matrix Track")
    return {col:tracks[col] for col in cols},sorted(set(tracks)-set(cols))
