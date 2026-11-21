#!/usr/bin/env python3
"""HS GeoMx 132 ROI targeted patient-aware descriptive sensitivity.
Source: public GSE294009 DCC/PKC Q3 reconstruction and published S2 ROI donor mapping.
Never treat ROIs as independent patients; no p-values, FDR, or treatment-response inference.
"""
import argparse,csv,collections,json,statistics,math
from pathlib import Path

REPO=Path("projects/hs_spatial_gse294009/author_s2_geo_crosswalk.csv")
PANEL={
 "B_cell":("CXCL13","MS4A1","CD79A","JCHAIN"),
 "Type17":("IL17A","IL17F","CCL20","TNF"),
 "Neutrophil":("S100A8","S100A9","CXCL8"),
 "Epithelial":("KRT17",),
}
FAMILIES={g:p for p,gs in PANEL.items() for g in gs}
PSEUDO=0.1
EXPECTED={"GSM":146,"HS_ROI":132,"AC_ROI":14,"HS_DONORS":15}
STAGE_SCENARIOS={
 "Yonsei_FFPE_exclude_HS12":(lambda x:x["institution"]=="Yonsei University" and x["patient"]!="HS12","author_stage_from_roi"),
 "Yonsei_FFPE_S2_include_HS12":(lambda x:x["institution"]=="Yonsei University","author_stage_from_roi"),
 "Yonsei_FFPE_S1_include_HS12":(lambda x:x["institution"]=="Yonsei University","clinical_stage_from_s1"),
 "All_sites_exclude_HS12":(lambda x:x["patient"]!="HS12","author_stage_from_roi"),
}
FOCUS={
 ("Inflammatory Aggregate_SuperficialDermis","CCL20"):"LATE_HIGHer_in_source",
 ("Inflammatory Aggregate_SuperficialDermis","IL17F"):"LATE_HIGHER_in_source",
 ("Inflammatory Aggregate_SuperficialDermis","TNF"):"LATE_HIGHER_in_source",
}
def load(path):
    with open(path,newline="",encoding="utf8") as handle:return list(csv.DictReader(handle))

def write(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",encoding="utf8",newline="") as handle:
        w=csv.DictWriter(handle,fieldnames=fields);w.writeheader();w.writerows(rows)

def parse(v):
    try:
        n=float(v)
        return n if math.isfinite(n) and n>=0 else None
    except (ValueError,TypeError):return None

def ratio(a,b):
    return math.log2((b+PSEUDO)/(a+PSEUDO))

def direction(v):
    return "LATE_HIGHER" if v>1e-9 else "EARLY_HIGHER" if v< -1e-9 else "TIE"

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",default="output/hs_spatial_gse294009")
    p.add_argument("--out",default="output/hs_same_site_sensitivity")
    a=p.parse_args()
    src=Path(a.input);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    m=load(REPO);meta={x["gsm"]:x for x in m}
    if len(meta)!=146 or len(m)!=146:raise RuntimeError("STOP S2 GSM mapping not unique 146")
    hs=[x for x in m if x["patient"].startswith("HS")]
    ac=[x for x in m if x["patient"].startswith("AC")]
    if len(hs)!=132 or len(ac)!=14 or len({x["patient"] for x in hs})!=15:raise RuntimeError("STOP unexpected HS source cohort")
    conflict=[x for x in hs if x["author_stage_from_roi"]!=x["clinical_stage_from_s1"]]
    if len(conflict)!=4 or {x["patient"] for x in conflict}!={"HS12"}:raise RuntimeError("STOP unexpected S1/S2 disagreement")
    expected_rois={x["gsm"] for x in hs}
    expr=load(src/"roi_gene_targeted_q3.csv")
    if not expr:raise RuntimeError("STOP missing expression")
    observed={x["gsm"] for x in expr}
    if len(observed)!=204 or not set(meta).issubset(observed):
        raise RuntimeError("STOP reconstructed DCC and S2 GSM mismatch")
    probe_genes={x["gene"] for x in expr}
    available=set(FAMILIES)&probe_genes
    if not available:raise RuntimeError("STOP no target genes mapped")
    # Each GSM/gene should only have one value; q3 is per ROI.
    uniq=set();pairs=collections.defaultdict(list);gene_nonmissing=collections.Counter()
    for row in expr:
        if row["gsm"] not in expected_rois or row["gene"] not in available:continue
        key=(row["gsm"],row["gene"])
        if key in uniq:raise RuntimeError("STOP duplicate ROI/gene")
        uniq.add(key)
        v=parse(row.get("q3_scaled_ratio"))
        if v is None:continue
        sub=meta[row["gsm"]]
        family_comp=sub["author_roi_classification"].split("_HS_",1)
        if len(family_comp)!=2:raise RuntimeError("STOP unexpected author ROI designation")
        comp=family_comp[1]
        if sub["slide_preparation"] not in ("Formalin-Fixed Paraffin-Embedded slide","Fresh Frozen slide"):
            raise RuntimeError("STOP unexpected sample preparation")
        pairs[(sub["patient"],row["gene"],comp)].append(v)
        gene_nonmissing[row["gene"]]+=1
    donor_rows=[]
    for (patient,gene,comp),values in pairs.items():
        m0=next(x for x in hs if x["patient"]==patient)
        donor_rows.append({**m0,"gene":gene,"compartment":comp,
                           "n_roi":len(values),"donor_median":statistics.median(values),
                           "roi_detection_fraction":sum(x>0 for x in values)/len(values)})
    outrows=[];support=[];primary={}
    for scenario,(predicate,stagecol) in STAGE_SCENARIOS.items():
        select=[x for x in donor_rows if predicate(x) and x[stagecol] in ("Early","Late")]
        keys=sorted({(x["gene"],x["compartment"]) for x in select})
        for gene,comp in keys:
            selected=[x for x in select if x["gene"]==gene and x["compartment"]==comp]
            e=[x for x in selected if x[stagecol]=="Early"]
            l=[x for x in selected if x[stagecol]=="Late"]
            if len({x["patient"] for x in e})!=len(e) or len({x["patient"] for x in l})!=len(l):
                raise RuntimeError("STOP duplicate donor within comparison")
            if scenario=="Yonsei_FFPE_exclude_HS12":
                support.append({"compartment":comp,"gene":gene,"panel":FAMILIES[gene],
                                "early_donors":len(e),"late_donors":len(l),
                                "n_early_ROI":sum(x["n_roi"] for x in e),"n_late_ROI":sum(x["n_roi"] for x in l)})
            if not e or not l:continue
            ev=[x["donor_median"] for x in e]; lv=[x["donor_median"] for x in l]
            v=ratio(statistics.median(ev),statistics.median(lv))
            deletion=[]
            if len(e)>1 and len(l)>1:
                for i in range(len(ev)):deletion.append(ratio(statistics.median(ev[:i]+ev[i+1:]),statistics.median(lv)))
                for i in range(len(lv)):deletion.append(ratio(statistics.median(ev),statistics.median(lv[:i]+lv[i+1:])))
            d=direction(v)
            stable=bool(deletion) and all(direction(z)==d for z in deletion) and d!="TIE"
            raw_directions={direction(z) for z in deletion}
            row={"scenario":scenario,"compartment":comp,"gene":gene,"panel":FAMILIES[gene],
                 "early_donors":len(e),"late_donors":len(l),"early_ROI":sum(x["n_roi"] for x in e),
                 "late_ROI":sum(x["n_roi"] for x in l),
                 "early_median_donor_Q3_ratio":round(statistics.median(ev),6),
                 "late_median_donor_Q3_ratio":round(statistics.median(lv),6),
                 "late_vs_early_log2_ratio_pseudo0p1":round(v,4),
                 "direction":d,"LOO_n":len(deletion),"LOO_direction_stable":stable,
                 "LOO_min":round(min(deletion),4) if deletion else "",
                 "LOO_max":round(max(deletion),4) if deletion else "",
                 "early_median_donor_detection_fraction":round(statistics.median(x["roi_detection_fraction"] for x in e),4),
                 "late_median_donor_detection_fraction":round(statistics.median(x["roi_detection_fraction"] for x in l),4),
                 "low_detection_flag":statistics.median(x["roi_detection_fraction"] for x in e)<.2 or statistics.median(x["roi_detection_fraction"] for x in l)<.2,
                 "inference":"DESCRIPTIVE_ONLY_NO_PVALUES"}
            outrows.append(row)
            if scenario=="Yonsei_FFPE_exclude_HS12":primary[(comp,gene)]=row
    if not primary:raise RuntimeError("STOP no within-site contrasts")
    # The original manuscript's claim is early immune activation vs comparator,
    # NOT a universal early > late contrast; check late directions only for explicit matches.
    directional=[]
    for (comp,gene),expect in FOCUS.items():
        r=primary.get((comp,gene))
        directional.append({"compartment":comp,"gene":gene,"paper_late_v_early_direction":"LATE_HIGHER",
                            "observed_direction":r["direction"] if r else "NO_MATCHED_DONORS",
                            "matching":("CONSISTENT" if r and r["direction"]=="LATE_HIGHER" else
                                        "OPPOSITE" if r and r["direction"]=="EARLY_HIGHER" else "UNRESOLVED"),
                            "source_context":"PMC13552446 late-vs-early expression changes in superficial dermis",
                            "status":"DIRECTIONAL_DESCRIPTIVE_NOT_QUANT_REPLICATION"})
    write(out/"HS_SAME_SITE_DONOR_PANEL.csv",outrows)
    write(out/"HS_PRIMARY_GENE_SUPPORT.csv",support)
    write(out/"HS_PAPER_DIRECTIONAL_AUDIT.csv",directional)
    num_primary=len(primary)
    high=[x for x in primary.values() if x["LOO_direction_stable"] and not x["low_detection_flag"] and x["early_donors"]>=3 and x["late_donors"]>=2]
    audit={"accession":"GSE294009","source":"raw DCC + PKC, source-linked official S2",
           "n_total_GEO_DCC":len(observed),"n_S2_GEO_HS_ROI":len(hs),"n_S2_GEO_AC_ROI":len(ac),
           "n_HS_donors":len({x["patient"] for x in hs}),"S1_S2_HS12_conflicted_ROIs":len(conflict),
           "primary":"Yonsei FFPE exclude HS12; median of ROI-level Q3 ratio by donor and compartment",
           "panels":{k:list(v) for k,v in PANEL.items()},
           "measured_target_genes":sorted(available),"n_primary_comparisons":num_primary,
           "n_high_support_stable_non_sparse":len(high),
           "status":"TARGETED_DESCRIPTIVE_GO_CONFIRMED_DE_HOLD",
           "limitations":["Source ratio not official author LOQ/Q3 normalization",
                          "Preselected targeted panel is not a transcriptome-wide test",
                          "Concordance to publication is descriptive direction only",
                          "No ROI-level p/FDR or donor-inference with these sample counts",
                          "Donor and site/fixation confounded in all-site comparison",
                          "HS12 discrepancy; no treatment response outcomes",
                          "24 previous GeoMx ROI do not have public matched quant data"]}
    (out/"HS_SENSITIVITY_AUDIT.json").write_text(json.dumps(audit,indent=2)+"\n",encoding="utf8")
    lines=["# HS spatial GeoMx — same-site targeted donor-level sensitivity","",
           "**Verdict: descriptive only. No formal DE/FDR.**","",
           "Primary = Yonsei FFPE; exclude discrepant HS12; median Q3 per patient and compartment;",
           "effect is log2[(late median + 0.1)/(early median + 0.1)]. No optimized pseudocount.",
           "Leave-one-donor-out direction checks omit each early or late donor in turn.",
           "A stable sign is NOT statistical significance or replication of official normalisation.","",
           "## High-support, LOO direction-stable candidate contrasts","",
           "|Compartment|Gene|Panel|Donors early/late|log2 late/early|LOO range|",
           "|---|---|---|---:|---:|---:|"]
    for x in sorted(high,key=lambda r:(r["compartment"],r["panel"],r["gene"])):
        lines.append(f'|{x["compartment"]}|{x["gene"]}|{x["panel"]}|{x["early_donors"]}/{x["late_donors"]}|{x["late_vs_early_log2_ratio_pseudo0p1"]}|{x["LOO_min"]} to {x["LOO_max"]}|')
    lines+=["","## Original manuscript targeted direction check","",
            "|Compartment|Gene|Paper late-vs-early|Our result|Consistency|","|---|---|---|---|---|"]
    for x in directional:lines.append(f'|{x["compartment"]}|{x["gene"]}|{x["paper_late_v_early_direction"]}|{x["observed_direction"]}|{x["matching"]}|')
    lines+=["","**Limits:** Early B-cell/Type17 immune activation is documented in the paper",
            "but not every gene is expected to be higher in early than late.",
            "Patient-level descriptive robustness alone cannot validate 156-ROI differential expression.",""]
    (out/"HS_SENSITIVITY_REPORT.md").write_text("\n".join(lines),encoding="utf8")
    print(json.dumps({"status":audit["status"],"primary_n":num_primary,
                      "high_support_stable":len(high),"high_support_genes":[{"gene":r["gene"],"compartment":r["compartment"],"n":str(r["early_donors"])+"/"+str(r["late_donors"]),"direction":r["direction"],"log2_ratio":r["late_vs_early_log2_ratio_pseudo0p1"]} for r in high],
                      "paper_directional":directional},ensure_ascii=False))

if __name__=="__main__":main()
