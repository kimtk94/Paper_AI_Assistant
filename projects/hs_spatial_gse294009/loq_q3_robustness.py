#!/usr/bin/env python3
"""GeoMx GSE294009 132-HS-ROI negative-probe/Q3/Q90 sensitivity.

This is a donor-aware, exploratory detectability audit, not author LOQ replication.
No inferential DEG p-values, pseudo-replicated ROIs, or treatment predictions.
"""
import argparse,collections,csv,datetime,json,math,statistics
from pathlib import Path

CROSS=Path("projects/hs_spatial_gse294009/author_s2_geo_crosswalk.csv")
SITES={"Yonsei University"}
GENES=("CXCL13","JCHAIN","MS4A1","CD79A","CXCL8","CCL20","IL17F","TNF")
COMPARTMENTS=("Inflammatory Aggregate_DeepDermis","Inflammatory Aggregate_SuperficialDermis","Epithelium_Epidermis")
VARIANTS=("Q3_all","Q90_all","Q3_LOQx1_zero","Q3_LOQx2_zero","Q3_SNR3_subset","Q3_SNR5_subset")
PSEUDOCOUNT=0.1
def load(path):
    with path.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))
def save(path,rows):
    if not rows:return
    names=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=names)
        w.writeheader();w.writerows(rows)
def number(value):
    try:
        x=float(value)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None
def sign(value):
    if value is None:return "NO_COMPARISON"
    return "LATE_HIGHER" if value>1e-9 else "EARLY_HIGHER" if value<-1e-9 else "TIE"
def ratio(e,l):
    return math.log2((l+PSEUDOCOUNT)/(e+PSEUDOCOUNT))
def make_donor_comparison(expr,meta,qc,gene,comp,variant):
    patient_values=collections.defaultdict(list)
    patient_detection=collections.defaultdict(list)
    eligible_roi=collections.defaultdict(set)
    missing_neg=0
    skipped_multi=0
    for r in expr:
        if r["gene"]!=gene:continue
        gsm=r["gsm"]
        m=meta.get(gsm)
        if not m or not m["patient"].startswith("HS"):continue
        if m["institution"] not in SITES or m["patient"]=="HS12":continue
        if m["slide_preparation"]!="Formalin-Fixed Paraffin-Embedded slide":continue
        pieces=m["author_roi_classification"].split("_HS_",1)
        if len(pieces)!=2 or pieces[1]!=comp:continue
        if m["author_stage_from_roi"] not in ("Early","Late"):continue
        q=qc.get(gsm)
        if q is None:raise RuntimeError("STOP: ROI QC missing GSM "+gsm)
        stage=m["author_stage_from_roi"]
        p=(m["patient"],stage)
        q3=number(r.get("q3_scaled_ratio"))
        q90=number(r.get("q90_scaled_ratio"))
        raw=number(r.get("raw_count"))
        loq=number(q.get("exploratory_LOQ_clamp1_geomean_geosd2"))
        snr=number(q.get("Q3_div_negative_geomean"))
        nprobe=number(r.get("gene_probe_count"))
        if raw is None or nprobe is None:raise RuntimeError("STOP: no gene raw count or probe count")
        if q3 is None:continue
        if variant=="Q90_all":
            v=q90
        elif variant=="Q3_all":
            v=q3
        elif variant in ("Q3_LOQx1_zero","Q3_LOQx2_zero"):
            # An aggregated multi-probe gene count is NOT directly comparable to
            # the background distribution of individual negative probes.
            if nprobe!=1:
                skipped_multi+=1
                continue
            if loq is None or loq<=0:
                missing_neg+=1
                continue
            threshold=loq*(2 if variant=="Q3_LOQx2_zero" else 1)
            v=q3 if raw>=threshold else 0.0
        else:
            if snr is None:
                missing_neg+=1
                continue
            threshold=3 if variant=="Q3_SNR3_subset" else 5
            if snr<threshold:continue
            v=q3
        if v is None:continue
        patient_values[p].append(v)
        eligible_roi[p].add(gsm)
        if nprobe==1 and loq is not None and loq>0:
            patient_detection[p].append(int(raw>=loq))
    e=[statistics.median(values) for (patient,stage),values in patient_values.items() if stage=="Early"]
    l=[statistics.median(values) for (patient,stage),values in patient_values.items() if stage=="Late"]
    e_det=[statistics.mean(vals) for (p,stage),vals in patient_detection.items() if stage=="Early" and vals]
    l_det=[statistics.mean(vals) for (p,stage),vals in patient_detection.items() if stage=="Late" and vals]
    ee=statistics.median(e) if e else None
    ll=statistics.median(l) if l else None
    effect=ratio(ee,ll) if ee is not None and ll is not None else None
    return {
        "gene":gene,"compartment":comp,"variant":variant,
        "early_donors":len(e),"late_donors":len(l),
        "early_ROIs":sum(len(v) for (p,stage),v in patient_values.items() if stage=="Early"),
        "late_ROIs":sum(len(v) for (p,stage),v in patient_values.items() if stage=="Late"),
        "early_donor_median":round(ee,6) if ee is not None else "",
        "late_donor_median":round(ll,6) if ll is not None else "",
        "log2_Late_over_Early_pseudo0p1":round(effect,6) if effect is not None else "",
        "direction":sign(effect),
        "early_donor_median_LOQ_detect_frac":round(statistics.median(e_det),4) if e_det else "",
        "late_donor_median_LOQ_detect_frac":round(statistics.median(l_det),4) if l_det else "",
        "neg_background_missing_rows":missing_neg,
        "non_single_probe_skipped_rows":skipped_multi,
        "inference":"DESCRIPTIVE_ONLY_NO_PVALUES"
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",default="output/hs_spatial_gse294009")
    ap.add_argument("--out",default="output/hs_loq_sensitivity")
    a=ap.parse_args()
    src=Path(a.input);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    source=json.loads((src/"RECONSTRUCTION_AUDIT.json").read_text())
    if source["accession"]!="GSE294009" or source["GSMs_matched"]!=204:raise RuntimeError("STOP source accession")
    if "exploratory_loq_method" not in source:raise RuntimeError("STOP missing fresh raw negative-probe reconstruction")
    m=load(CROSS);meta={x["gsm"]:x for x in m}
    if len(m)!=146 or len(meta)!=146:raise RuntimeError("STOP S2 donor crosswalk")
    hs=[x for x in m if x["patient"].startswith("HS")]
    ac=[x for x in m if x["patient"].startswith("AC")]
    if len(hs)!=132 or len(ac)!=14 or len({x["patient"] for x in hs})!=15:
        raise RuntimeError("STOP S2 cohort mismatch")
    discordant=[x for x in hs if x["clinical_stage_from_s1"]!=x["author_stage_from_roi"]]
    if len(discordant)!=4 or {x["patient"] for x in discordant}!={"HS12"}:
        raise RuntimeError("STOP unexpected stage disagreement")
    qc_rows=load(src/"roi_qc.csv")
    qc={x["gsm"]:x for x in qc_rows}
    if len(qc)!=204:raise RuntimeError("STOP QC source missing ROI")
    ex=load(src/"roi_gene_targeted_q3.csv")
    if len({x["gsm"] for x in ex})!=204:raise RuntimeError("STOP ROI matrix missing")
    available={x["gene"] for x in ex}
    missing=[x for x in GENES if x not in available]
    if missing:raise RuntimeError("STOP target missing in PKC "+",".join(missing))
    uniq=set()
    for x in ex:
        pair=(x["gsm"],x["gene"])
        if pair in uniq:raise RuntimeError("STOP duplicated target per ROI "+str(pair))
        uniq.add(pair)
    site_roi=[x for x in hs if x["institution"]=="Yonsei University" and x["patient"]!="HS12"]
    if len({x["patient"] for x in site_roi if x["author_stage_from_roi"]=="Early"})!=7:
        raise RuntimeError("STOP early donors")
    if len({x["patient"] for x in site_roi if x["author_stage_from_roi"]=="Late"})!=3:
        raise RuntimeError("STOP late donors")
    rows=[]
    for comp in COMPARTMENTS:
        for gene in GENES:
            for variant in VARIANTS:
                rows.append(make_donor_comparison(ex,meta,qc,gene,comp,variant))
    save(out/"HS_LOQ_Q3_Q90_DONOR_SENSITIVITY.csv",rows)
    index={(r["gene"],r["compartment"],r["variant"]):r for r in rows}
    ranks=[]
    for comp in COMPARTMENTS:
        for gene in GENES:
            z={v:index[(gene,comp,v)] for v in VARIANTS}
            base=z["Q3_all"]
            loq=z["Q3_LOQx1_zero"]
            loq2=z["Q3_LOQx2_zero"]
            q90=z["Q90_all"]
            snr3=z["Q3_SNR3_subset"]
            snr5=z["Q3_SNR5_subset"]
            d=base["direction"]
            observed=all(z[v]["direction"]==d for v in ("Q90_all","Q3_LOQx1_zero","Q3_LOQx2_zero"))
            valid_base=int(base["early_donors"])>=3 and int(base["late_donors"])>=2 and d not in ("TIE","NO_COMPARISON")
            has_loq=int(loq["early_donors"])>=3 and int(loq["late_donors"])>=2
            dete=number(loq["early_donor_median_LOQ_detect_frac"])
            detl=number(loq["late_donor_median_LOQ_detect_frac"])
            detectable=(dete is not None and detl is not None and dete>=.5 and detl>=.5)
            if not valid_base:decision="HOLD_INSUFFICIENT_DONORS"
            elif not has_loq:decision="HOLD_LOQ_NOT_EVALUABLE"
            elif not detectable:decision="HOLD_LOW_DETECTABILITY"
            elif observed:decision="RETAIN_EXPLORATORY_FOR_ORTHOGONAL_VALIDATION"
            else:decision="HOLD_QC_SENSITIVE_DIRECTION"
            ranks.append({
              "gene":gene,"compartment":comp,
              "decision":decision,"early_donors":base["early_donors"],"late_donors":base["late_donors"],
              "Q3_effect":base["log2_Late_over_Early_pseudo0p1"],
              "Q90_effect":q90["log2_Late_over_Early_pseudo0p1"],
              "LOQ1x_effect":loq["log2_Late_over_Early_pseudo0p1"],
              "LOQ2x_effect":loq2["log2_Late_over_Early_pseudo0p1"],
              "SNR3_effect":snr3["log2_Late_over_Early_pseudo0p1"],
              "SNR5_effect":snr5["log2_Late_over_Early_pseudo0p1"],
              "baseline_direction":d,
              "Q3_Q90_LOQ1_LOQ2_same_direction":observed,
              "median_early_LOQ_detect_fraction":dete if dete is not None else "",
              "median_late_LOQ_detect_fraction":detl if detl is not None else "",
              "SNR3_donors":str(snr3["early_donors"])+"/"+str(snr3["late_donors"]),
              "SNR5_donors":str(snr5["early_donors"])+"/"+str(snr5["late_donors"]),
              "LOQ_method":"EXPLORATORY_GEOMEAN_x_GEOSD2_CLAMP1"
            })
    save(out/"HS_QC_CANDIDATE_DECISION.csv",ranks)
    qcsite=collections.defaultdict(list)
    for x in site_roi:
        q=qc[x["gsm"]]
        qcsite[(x["author_stage_from_roi"],x["institution"])].append(q)
    def med(rows,k):
        vals=[number(x.get(k)) for x in rows]
        vals=[x for x in vals if x is not None]
        return round(statistics.median(vals),6) if vals else ""
    site_summary=[]
    for (stage,site),xs in sorted(qcsite.items()):
        site_summary.append({
          "stage":stage,"site":site,"n_ROI":len(xs),
          "n_donors":len({x["patient"] for x in site_roi if x["author_stage_from_roi"]==stage and x["institution"]==site}),
          "median_negative_probe_count":med(xs,"n_negative_probes"),
          "median_Q3":med(xs,"Q3_endogenous_raw"),
          "median_negative_geomean":med(xs,"negative_probe_geomean_clamp1"),
          "median_exploratory_LOQ":med(xs,"exploratory_LOQ_clamp1_geomean_geosd2"),
          "median_Q3_to_negative_geomean":med(xs,"Q3_div_negative_geomean"),
          "n_SNR_lt_3":sum((number(x["Q3_div_negative_geomean"]) or -1)<3 for x in xs),
          "n_SNR_lt_5":sum((number(x["Q3_div_negative_geomean"]) or -1)<5 for x in xs)
        })
    save(out/"HS_ROI_BACKGROUND_BY_STAGE.csv",site_summary)
    all_roi=[qc[x["gsm"]] for x in hs]
    nwith=sum(number(x["exploratory_LOQ_clamp1_geomean_geosd2"]) is not None for x in all_roi)
    key=[x for x in ranks if x["compartment"]=="Inflammatory Aggregate_DeepDermis" and x["gene"] in ("CXCL13","JCHAIN","MS4A1","CD79A")]
    audit={
      "status":"SOURCE_LOQ_PROXY_AUDITED__PATIENT_DESCRIPTIVE_ONLY",
      "timestamp":datetime.datetime.now(datetime.timezone.utc).isoformat(),
      "source_accession":"GSE294009","source_raw_tar_sha256":source.get("raw_tar_sha256"),
      "source_pkc_sha256":source.get("pkc_sha256"),
      "n_geo_roi":204,"n_hs_s2_roi":132,"n_ac_s2_roi":14,"n_hs_donors":15,
      "n_hs_roi_with_negative_probe_loq":nwith,
      "n_focus_gene_comparisons":len(ranks),
      "n_exploratory_retained":sum(x["decision"]=="RETAIN_EXPLORATORY_FOR_ORTHOGONAL_VALIDATION" for x in ranks),
      "n_qc_sensitive":sum(x["decision"].startswith("HOLD") for x in ranks),
      "candidate_evidence":key,
      "analysis":{"primary":"same-site Yonsei FFPE, HS12 excluded; median of ROI gene/Q3 per donor then median across donors",
                  "secondary":["Q90 endogenous raw scale","set to 0 if gene raw less than exploratory LOQ x1 or x2","keep only Q3/negative geomean >=3 or >=5"],
                  "loq":"max(2, negGeoMean * (negGeoSD^2)) from ROI DCC raw negative controls; clamp counts<1 before log; sample SD",
                  "probe_count_guard":"only 1-RTS gene compared with negative-probe LOQ",
                  "effect":"log2((late donor median + .1)/(early donor median + .1))",
                  "decision":"sign consistency Q3/Q90/LOQ masks and >=50% donor median LOQ detectability as exploratory heuristic",
                  "all_effects_exploratory":True,"p_values_computed":False,"FDR_computed":False},
      "limitations":[
        "Per-target Q3 fraction is not publisher's batch-adjusted LOQ normalized transcriptome.",
        "Negative control geomean*GeoSD^2 uses clamped zeros and sample SD, unlike potentially different original pipeline; not exact official LOQ.",
        "Some genes may have multiple probes; LOQ comparison blocked for those genes.",
        "Q3 and Q90 calculated from all endogenous PKC targets, no original gene filtering/LOQ batch model.",
        "Same-site n early<=7, late<=3 and fewer donors by compartment; retained candidates are NOT DEG hits.",
        "Q3/negative signal cutoff 3/5 is exploratory, not a manufacturer pass/fail threshold.",
        "No 24 historical source ROI raw data; full 156 ROI quantitative recreation HOLD.",
        "No treatment-response patient pairing; do not infer causal drug response."]
    }
    (out/"HS_LOQ_Q3_Q90_AUDIT.json").write_text(json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    lines=["# HS WTA DCC negative probe / Q3 / Q90 exploratory audit","",
           "**Status: patient-aware hypothesis generation only. No formal DEG/therapy response.**","",
           "- Original GEO ROI: 204; S2 HS ROI: 132 / 15 donors; AC control: 14 / 1 donor.",
           "- Negative-probe LOQ available for "+str(nwith)+"/132 HS ROIs. The LOQ is **an exploratory proxy**, not the author pipeline's LOQ.",
           "- Candidate genes are assessed in identical site (Yonsei) and preparation (FFPE), HS12 excluded.",
           "- Pseudocount=0.1 and all per-donor Q3/RNA normalisation proxies must accompany numerical results.",
           "- Retained: "+str(audit["n_exploratory_retained"])+"/"+str(len(ranks))+" gene-compartment comparisons under predeclared heuristic.",
           "",
           "## Targeted deep inflammatory dermis QC table","",
           "|Gene|Decision|Q3 effect|Q90 effect|LOQ1x|LOQ2x|Early/Late LOQ detectable median|",
           "|---|---|---:|---:|---:|---:|---|"]
    for x in key:lines.append("|"+str(x["gene"])+"|"+str(x["decision"])+"|"+str(x["Q3_effect"])+"|"+str(x["Q90_effect"])+"|"+str(x["LOQ1x_effect"])+"|"+str(x["LOQ2x_effect"])+"|"+str(x["median_early_LOQ_detect_fraction"])+"/"+str(x["median_late_LOQ_detect_fraction"])+"|")
    lines+=["","## Critical caveats","",
            "- LOQ threshold per gene is only validly compared with counts if exactly one RTS probe maps the gene.",
            "- Zero-censor at LOQ is a robustness stress test, not a recommended imputation method.",
            "- Q90 changes denominator compared with Q3 and effect magnitude is not directly comparable as official fold-change.",
            "- If the QC result is retained, seek orthogonal patient-level RNA/protein validation before making biomarker or treatment claims.",
            "- Detailed CSV and JSON are in this same artifact."]
    (out/"HS_LOQ_Q3_Q90_REPORT.md").write_text("\n".join(lines)+"\n",encoding="utf8")
    print(json.dumps({"status":audit["status"],"n_roi_with_loq":nwith,
                     "retained":audit["n_exploratory_retained"],"comparisons":len(ranks),
                     "primary_candidates":key,"site_summary":site_summary},ensure_ascii=False))
if __name__=="__main__":main()
