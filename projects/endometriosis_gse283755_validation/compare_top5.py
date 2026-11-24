#!/usr/bin/env python3
"""Compare five article Table 3 effects; never certify the full 406-gene set."""
import argparse
import csv
import json
import math
from pathlib import Path

SOURCE="https://doi.org/10.3349/ymj.2024.0518"
# Published numeric facts, Table 3, to the printed 2-decimal precision.
TOP5={"HLA-DOB":-1.81,"ELOVL2":-1.43,"NSFP1":-1.55,"GSTM1":3.22,"LRRC15":-2.10}

def compare_results(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if not {"gene","log2FoldChange","padj"}.issubset(reader.fieldnames or []):
            raise ValueError("required DESeq2 result columns absent")
        found={}
        for row in reader:
            gene=row["gene"]
            if not gene or gene in found:
                raise ValueError("blank or duplicate gene identifier")
            found[gene]=row
    rows=[]
    for gene, published in TOP5.items():
        row=found.get(gene)
        try:
            lfc=float(row["log2FoldChange"]); q=float(row["padj"])
            if not math.isfinite(lfc) or not math.isfinite(q) or not 0<=q<=1:
                raise ValueError("invalid effect or FDR")
        except (TypeError,ValueError,KeyError):
            lfc=q=None
        direction=lfc is not None and lfc*published>0
        match=direction and abs(lfc-published)<0.005
        significant=lfc is not None and q<0.05 and abs(lfc)>1
        rows.append({"gene":gene,"published_log2fc":published,"current_log2fc":lfc,
                     "current_padj":q,"direction_matches":direction,
                     "printed_effect_matches":match,"current_deg_threshold":significant})
    n=sum(x["printed_effect_matches"] and x["current_deg_threshold"] for x in rows)
    return {"source":SOURCE,"scope":"Table 3 top five only, same cohort",
            "status":"PASS_TOP5_PRINTED_EFFECTS__HOLD_FULL_DEG_SET" if n==5 else "HOLD_TOP5_DIFFERENCE",
            "matched":n,"expected":5,"complete_406_gene_reproduction":False,
            "independent_clinical_validation":False,"genes":rows}

def main():
    p=argparse.ArgumentParser();p.add_argument("--results",required=True);p.add_argument("--out",required=True)
    a=p.parse_args()
    report=compare_results(a.results)
    Path(a.out).write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report))
if __name__=="__main__":main()
