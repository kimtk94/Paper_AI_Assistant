#!/usr/bin/env python3
"""Rebuild all public GSE294009 WTA endogenous counts for S2-linked 132 HS ROIs.

This program DOES NOT reproduce original 156-ROI published expression or normalized DE.
It pins the original raw GeoMx DCC/PKC source to an earlier verified action,
keeps 132 unique GSMs and 15 donors, and exports gene x ROI integer counts.
"""
import argparse,collections,csv,gzip,hashlib,io,json,re,tarfile
from array import array
from pathlib import Path
from geo_file_probe import ACC,BASE,get,parse_series_matrix

CROSS=Path("projects/hs_spatial_gse294009/author_s2_geo_crosswalk.csv")
TARGETS={"JCHAIN","CXCL13","MS4A1","CD79A","CCL20","IL17F","TNF"}
def rows(p):
    with p.open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))
def write(p,table):
    if not table:return
    fields=list(table[0])
    with p.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(table)
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",default="output/hs_spatial_full_matrix")
    ap.add_argument("--source-audit",default="output/hs_spatial_gse294009/RECONSTRUCTION_AUDIT.json")
    args=ap.parse_args()
    dest=Path(args.out);dest.mkdir(parents=True,exist_ok=True)
    old=json.loads(Path(args.source_audit).read_text())
    if old.get("accession")!="GSE294009" or old.get("GSMs_matched")!=204:
        raise RuntimeError("STOP verified 204 ROI raw source audit missing")
    cross=rows(CROSS)
    if len(cross)!=146 or len({r["gsm"] for r in cross})!=146:
        raise RuntimeError("STOP author S2 crosswalk uniqueness")
    hs=sorted((x for x in cross if x["patient"].startswith("HS")),key=lambda v:int(v["gsm"][3:]))
    ac=[x for x in cross if x["patient"].startswith("AC")]
    if len(hs)!=132 or len(ac)!=14 or len({r["patient"] for r in hs})!=15:
        raise RuntimeError("STOP S2 study sample definition")
    target_gsms={r["gsm"] for r in hs}
    gsm_order=[r["gsm"] for r in hs]
    index={gsm:i for i,gsm in enumerate(gsm_order)}
    prior_stage={r["patient"] for r in hs if r["clinical_stage_from_s1"]!=r["author_stage_from_roi"]}
    if prior_stage!={"HS12"}:raise RuntimeError("STOP discordance differs from original")
    source,matrix_meta=get(BASE+"/matrix/"+ACC+"_series_matrix.txt.gz",timeout=105)
    pkc,pkc_meta=get(BASE+"/suppl/GSE294009_Hs_R_NGS_WTA_v1.0.pkc.gz",limit=35_000_000,timeout=125)
    tar_blob,tar_meta=get(BASE+"/suppl/GSE294009_RAW.tar",limit=20_000_000,timeout=155)
    if source is None or pkc is None or tar_blob is None:
        raise RuntimeError("STOP raw GEO source network unavailable: "+str((matrix_meta,pkc_meta,tar_meta)))
    for key,meta in (("series_sha256",matrix_meta),("pkc_sha256",pkc_meta),("raw_tar_sha256",tar_meta)):
        if old.get(key)!=meta.get("sha256"):
            raise RuntimeError("STOP raw source checksum changed for "+key)
    series={s["gsm"]:s for s in parse_series_matrix(source)}
    if len(series)!=204 or not target_gsms<=set(series):
        raise RuntimeError("STOP GEO series 204 sample integrity")
    pkc_obj=json.loads(gzip.decompress(pkc))
    probes={};gene_probe=collections.Counter();classes=collections.Counter()
    for target in pkc_obj.get("Targets",[]):
        gene=str(target.get("DisplayName","")).strip()
        codeclass=str(target.get("CodeClass","")).lower().strip()
        classes[codeclass]+=1
        if not gene:continue
        for pp in target.get("Probes",[]):
            rid=str(pp.get("RTS_ID","")).strip()
            if not rid:continue
            if rid in probes and probes[rid]!=(gene,codeclass):
                raise RuntimeError("STOP ambiguous RTS ID "+rid)
            probes[rid]=(gene,codeclass)
            if codeclass.startswith("endo"):gene_probe[gene]+=1
    names=sorted(gene_probe)
    if not 15000<=len(names)<=25000:
        raise RuntimeError("STOP unexpected WTA gene count "+str(len(names)))
    if not TARGETS<=set(names):raise RuntimeError("STOP required target gene missing")
    count_by_gene={g:array("I",[0])*len(gsm_order) for g in names}
    seen=set();raw_qc=[];unknown=0;member_total=0
    with tarfile.open(fileobj=io.BytesIO(tar_blob),mode="r:*") as archive:
        members=[x for x in archive.getmembers() if x.isfile() and ".dcc" in x.name.lower()]
        if len(members)!=204:raise RuntimeError("STOP DCC member count "+str(len(members)))
        for member in members:
            match=re.search(r"GSM\d+",member.name)
            if not match or match.group() not in series:
                raise RuntimeError("STOP DCC without author GEO GSM "+member.name)
            gsm=match.group()
            if gsm in seen:raise RuntimeError("STOP duplicate DCC GSM "+gsm)
            seen.add(gsm);member_total+=1
            if gsm not in target_gsms:continue
            b=archive.extractfile(member).read()
            if member.name.lower().endswith(".gz"):b=gzip.decompress(b)
            txt=b.decode("utf8",errors="replace")
            match=re.search(r"<Code_Summary>\s*(.*?)\s*</Code_Summary>",txt,re.S|re.I)
            if not match:raise RuntimeError("STOP missing DCC codes "+gsm)
            i=index[gsm];processed=0
            seen_codes=set();neg=0
            for line in match.group(1).splitlines():
                fields=line.strip().split(",")
                if len(fields)!=2 or not re.fullmatch(r"RTS\d+",fields[0].strip()) or not re.fullmatch(r"\d+",fields[1].strip()):
                    continue
                rid=fields[0].strip();val=int(fields[1].strip())
                if rid in seen_codes:raise RuntimeError("STOP duplicate RTS in "+gsm)
                seen_codes.add(rid);processed+=1
                rec=probes.get(rid)
                if rec is None:unknown+=1;continue
                gene,codeclass=rec
                if codeclass.startswith("endo"):
                    if count_by_gene[gene][i]+val>2**32-1:raise RuntimeError("STOP integer overflow")
                    count_by_gene[gene][i]+=val
                elif "negative" in codeclass:
                    neg+=1
            if processed<12000:raise RuntimeError("STOP too few RTS counts "+gsm)
            raw_qc.append({"gsm":gsm,"rts_observed":processed,"n_negative_probes":neg})
    if seen!=set(series) or len(raw_qc)!=132:
        raise RuntimeError("STOP expected 204 total GEO DCC and 132 author S2 HS DCC")
    filename=dest/"HS_132_HS_GENE_WTA_COUNTS.csv.gz"
    with gzip.open(filename,"wt",encoding="utf8",newline="",compresslevel=6) as f:
        w=csv.writer(f);w.writerow(["gene"]+gsm_order)
        for gene in names:w.writerow([gene]+list(count_by_gene[gene]))
    slide={}
    old_qc=Path(args.source_audit).parent/"roi_qc.csv"
    if old_qc.is_file():
        prev=rows(old_qc)
        slide={x["gsm"]:x["slide_label_not_donor"] for x in prev}
        if len(slide)!=204:raise RuntimeError("STOP source slide QC mismatch")
    manifest=[]
    for x in hs:
        label=x["author_roi_classification"].split("_HS_",1)
        if len(label)!=2:raise RuntimeError("STOP invalid source ROI compartment")
        manifest.append({
          "gsm":x["gsm"],"donor":x["patient"],
          "stage_s2":x["author_stage_from_roi"],
          "stage_s1":x["clinical_stage_from_s1"],
          "discordant_stage":x["discordant_stage"],
          "site":x["institution"],"preparation":x["slide_preparation"],
          "compartment":label[1],"author_roi_classification":x["author_roi_classification"],
          "slide_label_not_donor":slide.get(x["gsm"],"UNKNOWN")
        })
    meta_file=dest/"HS_132_AUTHOR_S2_MANIFEST.csv";write(meta_file,manifest)
    write(dest/"HS_132_RAW_ROI_PROBE_QC.csv",raw_qc)
    write(dest/"HS_WTA_GENE_PROBE_COUNTS.csv",
          [{"gene":g,"n_RTS_probes":gene_probe[g]} for g in names])
    comp=collections.defaultdict(lambda:collections.defaultdict(set))
    for m in manifest:
        comp[m["compartment"]][m["stage_s2"]].add(m["donor"])
    support=[{"compartment":k,"early_donors":len(c["Early"]),"late_donors":len(c["Late"]),
              "early_ROI":sum(x["compartment"]==k and x["stage_s2"]=="Early" for x in manifest),
              "late_ROI":sum(x["compartment"]==k and x["stage_s2"]=="Late" for x in manifest)}
             for k,c in sorted(comp.items())]
    write(dest/"HS_S2_GENE_COMPARTMENT_SUPPORT.csv",support)
    stats={
      "status":"FULL_WTA_PUBLIC_132_HS_SOURCE_COUNTS_VERIFIED__NO_DE",
      "accession":ACC,
      "original_GEO_DCC":204,"included_S2_GEO_HS_ROI":132,
      "excluded_GEO_AC_ROI":14,
      "historical_2023_roi_unavailable":24,
      "source_series_sha256":matrix_meta["sha256"],
      "source_pkc_sha256":pkc_meta["sha256"],
      "source_RAW_tar_sha256":tar_meta["sha256"],
      "matched_prior_reconstruction_sha256":True,
      "all_source_DCC_GSMs_matched":len(seen),
      "endogenous_gene_rows":len(names),
      "n_independent_HS_donors":15,
      "n_S1_S2_discordant_donors":len(prior_stage),
      "n_gene_probes_non_one":sum(1 for k in names if gene_probe[k]!=1),
      "n_unknown_RTS_instances":unknown,
      "matrix_csv_gz_sha256":sha(filename),
      "manifest_sha256":sha(meta_file),
      "source_gene_target_zero_warning":"All endogenous genes present in PKC; missing per-ROI RTS ID treated as zero, review original NanoString probe QC and imputation.",
      "count_model":"Raw probe counts, NOT standR-normalized; no author-equivalent 156 ROI DE.",
      "compartment_support":support
    }
    (dest/"HS_132_FULL_WTA_SOURCE_AUDIT.json").write_text(json.dumps(stats,indent=2)+"\n",encoding="utf8")
    print(json.dumps({"status":stats["status"],"gene_rows":len(names),
                      "donors":stats["n_independent_HS_donors"],
                      "ROI":stats["included_S2_GEO_HS_ROI"],
                      "support":support,
                      "counts_sha":stats["matrix_csv_gz_sha256"],
                      "count_probe_non_one":stats["n_gene_probes_non_one"]}))
if __name__=="__main__":main()
