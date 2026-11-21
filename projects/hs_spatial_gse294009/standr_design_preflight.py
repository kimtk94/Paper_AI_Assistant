#!/usr/bin/env python3
"""Independent-donor and design-rank feasibility before any standR/RUV4 model.

Reject ROI-as-replicate inference, inspect stage vs slide/site collinearity,
freeze the 132 S2 HS source population and early/late overlap per compartment.
"""
import argparse,collections,csv,gzip,json,hashlib
from pathlib import Path

CORE="Inflammatory Aggregate_DeepDermis"
def read(p):
    with p.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))
def rank(rows,tol=1e-9):
    a=[[float(z) for z in row] for row in rows]
    if not a:return 0
    n=len(a);m=len(a[0]);r=0
    for j in range(m):
        pivot=max(range(r,n),key=lambda i:abs(a[i][j]),default=None)
        if pivot is None or abs(a[pivot][j])<tol:continue
        a[r],a[pivot]=a[pivot],a[r]
        f=a[r][j]
        for k in range(j,m):a[r][k]/=f
        for i in range(n):
            if i==r:continue
            fac=a[i][j]
            if abs(fac)<tol:continue
            for k in range(j,m):a[i][k]-=fac*a[r][k]
        r+=1
        if r>=n:break
    return r
def cols_stage_slider(ms):
    slides=sorted(set(x["slide_label_not_donor"] for x in ms))
    donors=sorted(set(x["donor"] for x in ms))
    late=[int(x["stage_s2"]=="Late") for x in ms]
    base=[[1]+[int(x["slide_label_not_donor"]==z) for z in slides[1:]] for x in ms]
    with_stage=[a+[b] for a,b in zip(base,late)]
    dono=[[1]+[int(x["donor"]==p) for p in donors[1:]] for x in ms]
    return {
      "n_slides":len(slides),"slide_labels":slides,
      "rank_intercept_slide":rank(base),
      "rank_plus_stage":rank(with_stage),
      "stage_estimable_after_slide_fixed_effect":rank(with_stage)>rank(base),
      "n_donors":len(donors),
      "rank_intercept_donor":rank(dono),
      "rank_donor_plus_stage":rank([a+[b] for a,b in zip(dono,late)]),
      "stage_estimable_with_donor_fixed_effect":rank([a+[b] for a,b in zip(dono,late)])>rank(dono),
      "n_late_independent_donors":len({x["donor"] for x in ms if x["stage_s2"]=="Late"}),
      "n_early_independent_donors":len({x["donor"] for x in ms if x["stage_s2"]=="Early"}),
    }
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",default="output/hs_spatial_full_matrix")
    ap.add_argument("--out",default="output/hs_standr_design_preflight")
    a=ap.parse_args()
    source=Path(a.input);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    audit=json.loads((source/"HS_132_FULL_WTA_SOURCE_AUDIT.json").read_text())
    if audit["status"]!="FULL_WTA_PUBLIC_132_HS_SOURCE_COUNTS_VERIFIED__NO_DE":
        raise RuntimeError("STOP source matrix absent")
    ms=read(source/"HS_132_AUTHOR_S2_MANIFEST.csv")
    if len(ms)!=132 or len({x["gsm"] for x in ms})!=132 or len({x["donor"] for x in ms})!=15:
        raise RuntimeError("STOP 132/15 source S2 cohort mismatch")
    full_matrix=source/"HS_132_HS_GENE_WTA_COUNTS.csv.gz"
    if hashlib.sha256(full_matrix.read_bytes()).hexdigest()!=audit["matrix_csv_gz_sha256"]:
        raise RuntimeError("STOP source matrix checksum changed")
    with gzip.open(full_matrix,"rt",encoding="utf8",newline="") as f:
        h=next(csv.reader(f))
    if h[0]!="gene" or set(h[1:])!={x["gsm"] for x in ms}:raise RuntimeError("STOP matrix and sample manifest mismatch")
    if any(x["slide_label_not_donor"] in ("","UNKNOWN") for x in ms):
        raise RuntimeError("STOP missing slide metadata for model")
    full=cols_stage_slider(ms)
    primary=[x for x in ms if x["site"]=="Yonsei University" and
             x["preparation"]=="Formalin-Fixed Paraffin-Embedded slide" and
             x["donor"]!="HS12"]
    if len(primary)!=82 or len({x["donor"] for x in primary})!=10:
        raise RuntimeError("STOP primary 82 ROI / 10 patients not recovered")
    if any(x["discordant_stage"]!="0" for x in primary):
        raise RuntimeError("STOP S1/S2 conflict not excluded")
    if (len({x["donor"] for x in primary if x["stage_s2"]=="Early"}),
        len({x["donor"] for x in primary if x["stage_s2"]=="Late"}))!=(7,3):
        raise RuntimeError("STOP primary 7 Early / 3 Late independent donors")
    supports=[]
    for comp in sorted({x["compartment"] for x in primary}):
        sub=[x for x in primary if x["compartment"]==comp]
        desc=cols_stage_slider(sub)
        supports.append({
           "compartment":comp,"n_roi":len(sub),
           "n_early_roi":sum(x["stage_s2"]=="Early" for x in sub),
           "n_late_roi":sum(x["stage_s2"]=="Late" for x in sub),
           "n_early_donors":desc["n_early_independent_donors"],
           "n_late_donors":desc["n_late_independent_donors"],
           "n_slide_labels":desc["n_slides"],
           "rank_intercept_slide":desc["rank_intercept_slide"],
           "rank_plus_stage":desc["rank_plus_stage"],
           "stage_estimable_after_slide_fixed_effect":desc["stage_estimable_after_slide_fixed_effect"],
           "note":"Slide fixed-effect alias check does not substitute for donor-blocked limma and RUV4"}
        )
    deep=[x for x in supports if x["compartment"]==CORE]
    if len(deep)!=1 or (deep[0]["n_early_donors"],deep[0]["n_late_donors"])!=(7,3):
        raise RuntimeError("STOP deep dermis independent donor coverage")
    within=cols_stage_slider(primary)
    # Stage is constant for each patient. Never fit stage simultaneously
    # with donor fixed effects: the stage coefficient is aliased.
    if within["stage_estimable_with_donor_fixed_effect"]:
        raise RuntimeError("STOP expected nesting of stage within donor failed")
    decision={
      "status":"PRIMARY_FFPE_SOURCE_AND_DESIGN_READY__CONFIRMATORY_STAGE_DE_HOLD",
      "source_accession":"GSE294009",
      "n_raw_GEO_DCC":audit["original_GEO_DCC"],
      "n_gene_rows":audit["endogenous_gene_rows"],
      "n_geo_hs_ROI":132,"n_HS_donors":15,
      "primary_site":"Yonsei University",
      "primary_preparation":"FFPE",
      "excluded_conflicted_donor":"HS12",
      "primary_ROIs":82,"primary_early_ROIs":57,"primary_late_ROIs":25,
      "primary_early_independent_donors":7,"primary_late_independent_donors":3,
      "full_source_model_rank":full,
      "same_site_primary_model_rank":within,
      "compartments":supports,
      "ruv_k_required":2,
      "original_normalization_method_parameter":"UNSPECIFIED__TEST_TMM_AND_UPPERQUARTILE_NOT_PAPER_EQUIVALENCE",
      "NCG_selection":"REQUIRES_DEFINED_200_LEAST_VARIABLE_GENES_BY_SLIDE_NOT_AUTHORS_EXACT",
      "limma_stage_design_with_donor_fixed_effect":"FORBIDDEN_ALIASED",
      "limma_stage_design_with_donor_block_correlation":"PLANNED_EXPLORATORY_NOT_CONFIRMED",
      "full_156_expression":"HOLD_MISSING_2023_24_ROI",
      "DE_inference":"HOLD_UNTIL_VALID_RUV_AND_DONOR_BLOCK_PROVEN",
      "geographic_batch_caution":"Late UC Davis Frozen and Yonsei FFPE confounded with preparation; do not estimate joint site/fixation/stage causal effect"
    }
    (out/"HS_132_STANDR_MODEL_PREFLIGHT.json").write_text(json.dumps(decision,indent=2)+"\n",encoding="utf8")
    with (out/"HS_132_COMPARTMENT_DONOR_SUPPORT.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(supports[0]));w.writeheader();w.writerows(supports)
    print(json.dumps({"status":decision["status"],"primary_n_roi":len(primary),
      "primary_n_donors":10,"primary_stage_after_slide_estimable":within["stage_estimable_after_slide_fixed_effect"],
      "donor_fixed_alias":not within["stage_estimable_with_donor_fixed_effect"],
      "deep":deep[0],"compartments":supports},ensure_ascii=False))
if __name__=="__main__":main()
