#!/usr/bin/env python3
"""Compare verified author GeoMx methods with exploratory Q3/LOQ audit.

Fail closed: robust source retrieval and technical exploratory JCHAIN direction
do not establish exact normalization, within-donor stage effect, or novel biomarker.
"""
import argparse,datetime as dt,json,csv,re
from pathlib import Path

PDF_SHA="02546a7a94e71e6c1a7551d32f9a97a1e9adbc13ef2d9a38825810e29b8e113d"
PDF_NAME="NIHMS2203138-supplement-1.pdf"
METHODS_DOI="10.1016/j.jid.2026.01.041"
EXTERNAL_DOI="10.1038/s41467-026-70883-5"
def load(p):
    if not p.is_file():raise RuntimeError("STOP absent source artifact "+str(p))
    return json.loads(p.read_text(encoding="utf8"))
def excerpts(x):
    # All strings are extracted from the DOI-matched original PDF with page labels.
    rows=[]
    for collection in (x.get("method_contexts",{}),x.get("author_pipeline_anchor_context",{})):
        for hits in collection.values():
            for item in hits:
                rows.append(str(item.get("context") or item.get("text") or ""))
    return " ".join(rows)
def write_csv(path,rows):
    if not rows:return
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--methods",default="output/hs_geomx_published_methods/ORIGINAL_METHODS_AUDIT.json")
    ap.add_argument("--qc",default="output/hs_loq_sensitivity/HS_LOQ_Q3_Q90_AUDIT.json")
    ap.add_argument("--out",default="output/hs_geomx_methods_equivalence")
    a=ap.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    methods=load(Path(a.methods));qc=load(Path(a.qc))
    if methods.get("source_file")!=PDF_NAME or methods.get("pdf_sha256")!=PDF_SHA:
        raise RuntimeError("STOP published Supplementary Methods digest / provenance drift")
    if methods.get("source_publication_doi")!=METHODS_DOI:
        raise RuntimeError("STOP wrong supplementary paper")
    if methods.get("pdf_pages")!=8:raise RuntimeError("STOP PDF page count changed")
    if qc.get("status")!="SOURCE_LOQ_PROXY_AUDITED__PATIENT_DESCRIPTIVE_ONLY":
        raise RuntimeError("STOP previous exploratory QC evidence unavailable")
    if (qc.get("n_geo_roi"),qc.get("n_hs_s2_roi"),qc.get("n_hs_donors"))!=(204,132,15):
        raise RuntimeError("STOP cohort provenance mismatch")
    if qc.get("n_hs_roi_with_negative_probe_loq")!=132:
        raise RuntimeError("STOP negative-control Q3 QC missing")
    source=excerpts(methods).lower()
    required={
      "standR_geomxNorm":bool(re.search(r"geo-\s*mxnorm|geomxnorm",source)),
      "standR_RUV4":("ruv4" in source),
      "standR_geomxBatchCorrection":("geomxbatchcorrection" in source),
      "standR_RUV_k2":bool(re.search(r"k,\s*which\s+represents\s+unwanted\s+variation\s+factors.{0,80}set\s+to\s+2",source)),
      "limma_voom":("limma-voom" in source),
      "limma_RUV_design_covariates":bool(re.search(r"design\s+matrix.{0,200}ruv",source)),
      "robust_eBayes":("ebayes" in source and "robust" in source),
      "BH_adjusted_p":("benjami" in source and "hochberg" in source)
    }
    absent=[k for k,v in required.items() if not v]
    if absent:raise RuntimeError("STOP source pipeline incomplete: "+", ".join(absent))
    if not all(x in qc.get("analysis",{}).get("loq","").lower() for x in ("neggeomean","neggeosd")):
        raise RuntimeError("STOP exploratory LOQ source method altered")
    bykey={(x["gene"],x["compartment"]):x for x in qc.get("candidate_evidence",[])}
    comp="Inflammatory Aggregate_DeepDermis"
    target=bykey.get(("JCHAIN",comp))
    if not target or target.get("decision")!="RETAIN_EXPLORATORY_FOR_ORTHOGONAL_VALIDATION":
        raise RuntimeError("STOP JCHAIN baseline screen status changed")
    comparisons=[
      {"step":"GeoMx raw DCC/PKC mapping","paper":"GeoMx NGS Pipeline DCC","current":"204 DCC RTS-to-WTA PKC mapping","equivalence":"PARTIAL_GEO_SOURCE"},
      {"step":"Expression normalization","paper":"standR::geomxNorm; referenced quantile-normalized expression","current":"raw gene counts / ROI Q3 and Q90 exploratory denominators","equivalence":"NOT_REPRODUCED"},
      {"step":"Batch correction","paper":"standR::geomxBatchCorrection, RUV4 factors k=2; additional limma batch correction mentioned for technical variance","current":"same-site FFPE restriction; no RUV4 fitted","equivalence":"NOT_REPRODUCED"},
      {"step":"Statistical design","paper":"limma-voom DGEList, disease status + RUV factors; voom weights, robust eBayes and BH","current":"independent-donor ROI-median comparison; no adjusted p/FDR","equivalence":"NOT_REPRODUCED"},
      {"step":"LOQ threshold","paper":"No explicit LOQ equation located in extracted original supplement PDF; this does not prove absent QC","current":"exploratory negGeoMean x (negGeoSD)^2 from clamped negative-probe counts","equivalence":"UNKNOWN_ORIGINAL_CONFIGURATION"},
      {"step":"Study population","paper":"132 source GEO HS ROI plus 24 historic early HS ROI in integrated 156","current":"132 source GEO HS ROI; missing original 24 expression","equivalence":"PARTIAL_132_ONLY"},
      {"step":"Patient correlation","paper":"authors describe latent donor/slide RUV4 variation mitigation; no independent donor-level proof from this text alone","current":"compare donor medians in same site; n<=7 early/3 late","equivalence":"DIFFERENT_METHOD_WITH_DONOR_GUARD"}
    ]
    # Literature assay validation is for HS vs normal tissue, not early-vs-late
    # deep inflammatory dermis. Do not upgrade cross-study evidence to stage proof.
    verdict={
      "run_date_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
      "source": {"original_doi":METHODS_DOI,"original_pdf":PDF_NAME,"sha256":PDF_SHA,
                 "original_page_methods":[1,2]},
      "source_pipeline_presence_verified":required,
      "matched_original_normalization":False,
      "matched_original_batch_correction":False,
      "matched_original_DE_inference":False,
      "verified_public_full_156_roi_quantification":False,
      "technical_exploratory_JCHAIN":"RETAIN_BY_PROXY_QC_ONLY",
      "external_JCHAIN_HS_lesion_association":{
          "status":"SUPPORTED_BY_PUBLISHED_ORTHOGONAL_ASSAYS_NOT_NEW_RESULT",
          "doi":EXTERNAL_DOI,
          "reported_assays":["lesion-versus-control JCHAIN qPCR","J chain protein Western blot",
                             "immunofluorescence IgA localization","lesion versus non-lesion RNA-seq GSE155176",
                             "scRNA expression GSE249793 B/plasma cells"],
          "clinical_stage_specific_deep_dermis_validation":False,
          "novel_gene_biomarker_claim_allowed":False
      },
      "stage_specific_JCHAIN":"HOLD_UNVERIFIED",
      "patient_level_independence_or_therapy_prediction":"HOLD",
      "scientific_decision":"GO_FOR_RESEARCH_HYPOTHESIS_ONLY__HOLD_CLAIM_OF_PAPER_REPRODUCTION",
      "comparison":comparisons,
      "next_steps":["Request original normalized, RUV-adjusted matrices or replicate standR on 132 with locked sample design",
                    "Verify original 2023 24-ROI source before 156-ROI reconstruction",
                    "Compare JCHAIN within compartment by individual donors under full probe/LOQ and stage model",
                    "Do not present JCHAIN HS expression or protein elevation as novel discovery"],
    }
    (out/"HS_METHOD_EQUIVALENCE_DECISION.json").write_text(json.dumps(verdict,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    write_csv(out/"HS_METHOD_PIPELINE_GAPS.csv",comparisons)
    lines=["# Original GSE294009 GeoMx Methods vs our source-data exploratory analysis","",
      "**Decision: HYPOTHESIS GO / exact original-normalization and inferential replication HOLD.**","",
      "- Original Methods DOI: https://doi.org/"+METHODS_DOI,
      "- Verified original official supplement: "+PDF_NAME+", 8 pages, SHA256 "+PDF_SHA+".",
      "- Source Methods page 1 describes standR geomxNorm and RUV4 via geomxBatchCorrection with **k=2**.",
      "- Page 2 describes limma-voom with disease status and RUV factors, robust eBayes and BH correction.",
      "- The retrieved original supplement does not give an explicit LOQ equation in searchable text; a typical GeoMxTools formula is not automatically the authors' method.",
      "- Our raw/Q3/Q90 and negative-control LOQ computations are exploratory, not equivalent to original integrated or batch-adjusted values.",
      "","## Method differences","","|Step|Original source|Our analysis|Decision|","|---|---|---|---|"]
    for x in comparisons:lines.append("|"+x["step"]+"|"+x["paper"]+"|"+x["current"]+"|"+x["equivalence"]+"|")
    lines+=["","## JCHAIN independent validation and novelty","",
       "- An independent 2026 Nature Communications study (https://doi.org/"+EXTERNAL_DOI+") already reports HS skin JCHAIN elevated by qPCR, J chain Western blot, and IGHA/JCHAIN expression in B/plasma cells using public single-cell data.",
       "- Reported qPCR cohort: HS lesion n=15, healthy controls n=3. Bulk paired lesions/nonlesional RNA cohorts GSE155176: n=15 per group; nonlesional and normal expression differ.",
       "- **This is independent support for JCHAIN in HS lesions, not an independent stage-specific deep-dermis test.**",
       "- JCHAIN's identity as a HS lesion-associated gene is not novel; potential research question is within-same-compartment early/late trajectory, but this dataset lacks confirmatory statistical design.","",
       "## Final gates","",
       "1. Keep JCHAIN as technical follow-up/hypothesis only; no validated early/late biomarker claim.",
       "2. Keep CXCL13, MS4A1 and CD79A stage effects HOLD after exploratory LOQ robustness failures.",
       "3. Do not report original paper's DEG FDR or compare ratios as if author-equivalent.",
       "4. Preserve paper SHA, patient/ROI provenance and analysis mode in downstream artifacts.",""]
    (out/"HS_METHOD_EQUIVALENCE_DECISION.md").write_text("\n".join(lines),encoding="utf8")
    print(json.dumps({"decision":verdict["scientific_decision"],"original_methods_present":required,
                      "n_method_gaps":len(comparisons),
                      "jchain":verdict["technical_exploratory_JCHAIN"],
                      "stage_specific":verdict["stage_specific_JCHAIN"],
                      "external_DOI":EXTERNAL_DOI},ensure_ascii=False))
if __name__=="__main__":main()
