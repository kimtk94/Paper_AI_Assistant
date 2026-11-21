#!/usr/bin/env python3
"""Independently audit completed HS GSE294009 full-nuisance donor LOO results.

This validates the actual Actions artifact, NOT only workflow completion.
44 successful model fits establish directional robustness, not DE significance.
"""
import argparse,json,math,hashlib
from pathlib import Path
NORMS=("TMM","upperquartile")
MODELS=("Stage_RUV","Stage_Slide_RUV")
EXPECTED_STATUS="FULL_NUISANCE_REESTIMATED_DONOR_LOO_AND_SLIDE_DESIGN_AUDIT__NOT_PAPER_REPLICATION"
def finite(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",default="output/hs_standr_full_refit/HS_FULL_REFIT_SLIDE_DONOR_LOO_AUDIT.json")
    p.add_argument("--out",default="output/hs_full_refit_final_gate")
    a=p.parse_args()
    path=Path(a.input)
    if not path.is_file():raise RuntimeError("STOP: full-refit result artifact missing")
    raw=path.read_bytes()
    j=json.loads(raw)
    if j["status"]!=EXPECTED_STATUS:raise RuntimeError("STOP: wrong experimental model status")
    if j["source_accession"]!="GSE294009":raise RuntimeError("STOP: wrong GEO source")
    expected=(j["n_source_HS_roi"],j["n_primary_HS_roi"],j["n_primary_donors"],
              j["n_primary_early_donors"],j["n_primary_late_donors"],j["n_primary_deep_roi"])
    if expected!=(132,82,10,7,3,24):
        raise RuntimeError("STOP: unexpected independently sampled cohort")
    if j["source_matrix_sha256"]!="943dff9ffcd6a9b75b0356a1c388863b48ec95804cf1d004ca17571a4ce5c169":
        raise RuntimeError("STOP: source matrix checksum drift")
    if any(j[x] is not False for x in ("original_156_ROI_reproduced",
                                       "confirmed_DE_biomarker",
                                       "paper_equivalent_normalization")):
        raise RuntimeError("STOP: unverified published effect claim")
    if j["full_reestimation_per_heldout_donor"] is not True or j["slide_fixed_effect_explicitly_tested"] is not True:
        raise RuntimeError("STOP: full nuisance refit or slide not tested")
    if set(j["results"])!=set(NORMS) or len(j["summaries"])!=4:
        raise RuntimeError("STOP: incomplete two-by-two model grid")
    rows=[];success_count=0;loo_count=0;negative_count=0;all_withheld=set()
    for norm in NORMS:
        block=j["results"][norm]
        if len(block)!=22:raise RuntimeError(f"STOP: expected 22 model fits for {norm}")
        for model in MODELS:
            label=f"{norm}__{model}"
            summary=j["summaries"][label]
            if summary["normalization"]!=norm or summary["model"]!=model:
                raise RuntimeError("STOP: summary key disagreement")
            full_key="FULL__"+model
            base=block[full_key]
            if base["donor_excluded"]!="NONE" or base["normalization"]!=norm or base["model"]!=model:
                raise RuntimeError("STOP: full-model bookkeeping mismatch")
            if base["status"]!="PASS_MODEL_FIT" or not finite(base.get("coef")):
                raise RuntimeError("STOP: primary slide/non-slide model did not fit")
            if base["n_upstream_roi"]!=82 or base["n_deep_roi"]!=24:
                raise RuntimeError("STOP: original model used wrong ROI")
            if (base["n_early_donors"],base["n_late_donors"])!=(7,3):
                raise RuntimeError("STOP: original model wrong biological replicates")
            if not 100<=base["n_NCG"]<=200:
                raise RuntimeError("STOP: invalid NCG support")
            if base["rank"]!=base["n_cols"] or base["full_nuisance_reestimation"] is not True:
                raise RuntimeError("STOP: model collinear or nuisance reused")
            success_count+=1
            withheld_rows=[block[x] for x in block if x.startswith("WITHOUT_") and x.endswith("__"+model)]
            if len(withheld_rows)!=10:raise RuntimeError("STOP: incomplete held-out patient coverage")
            effects=[]
            withheld=set()
            for z in withheld_rows:
                if z["status"]!="PASS_MODEL_FIT" or z["normalization"]!=norm or z["model"]!=model:
                    raise RuntimeError("STOP: unreported failed donor fold")
                if z["donor_excluded"]=="NONE" or z["donor_excluded"] in withheld:
                    raise RuntimeError("STOP: duplicated or missing donor holdout")
                withheld.add(z["donor_excluded"])
                all_withheld.add(z["donor_excluded"])
                if not 0<z["n_upstream_roi"]<82 or not 100<=z["n_NCG"]<=200:
                    raise RuntimeError("STOP: donor exclusion did not change upstream cohort")
                if z["rank"]!=z["n_cols"] or z["full_nuisance_reestimation"] is not True:
                    raise RuntimeError("STOP: held-out design quality")
                if z["n_late_donors"]<2 or z["n_early_donors"]<3:
                    raise RuntimeError("STOP: insufficient patient replication")
                if not finite(z.get("coef")):
                    raise RuntimeError("STOP: no stage effect from successful fold")
                effects.append(float(z["coef"]))
                loo_count+=1
            if not all(v>0 for v in effects) or float(base["coef"])<=0:
                negative_count+=sum(v<=0 for v in effects)
                raise RuntimeError(f"STOP: observed direction reversal in {label}")
            if (summary["base_status"]!="PASS_MODEL_FIT" or
                summary["n_donors_loo"]!=10 or
                summary["n_successful_refits"]!=10 or
                summary["n_HOLD"]!=0 or
                summary["all_success_and_same_direction"] is not True):
                raise RuntimeError("STOP: summary disagrees with donor-wise rows")
            rng=summary["loo_effect_range"]
            if len(rng)!=2 or max(abs(min(effects)-rng[0]),abs(max(effects)-rng[1]))>0.015:
                raise RuntimeError("STOP: reported effect range disagrees with original values")
            rows.append({"normalization":norm,"model":model,
                         "n_model_success":11,"n_donor_refits":10,
                         "baseline_effect":round(base["coef"],4),
                         "loo_min":round(min(effects),4),"loo_max":round(max(effects),4),
                         "n_neg_or_zero":0})
    if success_count!=4 or loo_count!=40 or len(all_withheld)!=10:
        raise RuntimeError("STOP: 4+40 models/10 independent patients not proven")
    out=Path(a.out);out.mkdir(exist_ok=True,parents=True)
    result={
      "status":"EXPLORATORY_FULL_REFIT_AND_SLIDE_DIRECTION_PASS__PAPER_AND_BIOMARKER_HOLD",
      "verified_input_sha256":hashlib.sha256(raw).hexdigest(),
      "source_matrix_sha256":j["source_matrix_sha256"],
      "n_full_model_success":success_count,
      "n_fold_model_success":loo_count,
      "n_unique_heldout_donors":len(all_withheld),
      "n_fold_direction_reversal":negative_count,
      "strictly_exploratory":True,
      "study_replication":False,
      "differential_expression_significance":False,
      "note":"Folds reestimate standR/NCG/RUV4 but estimate coefficients only from retained patients; not test-set prediction or independent biological validation. n Late donors=3.",
      "model_summary":rows
    }
    (out/"HS_FINAL_QC_GATE.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(result,ensure_ascii=False))
    print("PASS_FULL_NUISANCE_REESTIMATED_SLIDE_44_FITS_40_POSITIVE_NO_FDR_CLAIM")
if __name__=="__main__":main()
