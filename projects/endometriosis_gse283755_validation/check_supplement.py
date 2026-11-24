#!/usr/bin/env python3
"""Fetch original article Supplementary Table 2, extract evidence without fabricating 406-gene overlap."""
import argparse, hashlib, json, os, re, shutil, subprocess, time, urllib.request
from pathlib import Path
from collections import Counter

URLS=[
 'https://pmc.ncbi.nlm.nih.gov/articles/instance/12568438/bin/ymj-66-780-s002.pdf',
 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12568438/bin/ymj-66-780-s002.pdf',
 'https://www.eymj.org/Synapse/Data/PDFData/0069YMJ/ymj-66-780-s002.pdf',
 'https://www.ymj.or.kr/Synapse/Data/PDFData/0069YMJ/ymj-66-780-s002.pdf',
 'https://ymj.or.kr/Synapse/Data/PDFData/0069YMJ/ymj-66-780-s002.pdf',
 'https://synapse.koreamed.org/upload/SynapseData/PDFData/0069YMJ/ymj-66-780-s002.pdf',
]
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args()
 dest=Path(a.out);dest.mkdir(parents=True,exist_ok=True)
 report={'accession':'GSE283755','supplement':'Supplementary Table 2, original DEGs','fetched':False,'attempts':[],'status':'HOLD_SUPPLEMENT_ACCESS'}
 for url in URLS:
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 (scientific supplementary data audit)'})
   with urllib.request.urlopen(req,timeout=35) as resp:
    content=resp.read(2_000_000)
   if not content.startswith(b'%PDF') or len(content)<30000:raise ValueError('non-PDF or unexpectedly short file')
   path=dest/'original_supplement_2.pdf';path.write_bytes(content)
   report.update({'fetched':True,'url':url,'pdf_sha256':hashlib.sha256(content).hexdigest(),'pdf_bytes':len(content),'status':'PDF_FETCHED_NOT_PARSED'})
   break
  except Exception as e:
   report['attempts'].append({'url':url,'error':str(e)[:180]})
   time.sleep(.5)
 if report['fetched']:
  pdf=str(dest/'original_supplement_2.pdf');txt=str(dest/'supplement_table2_text.txt')
  try:
   if not shutil.which('pdftotext'):raise RuntimeError('pdftotext binary unavailable')
   subprocess.run(['pdftotext','-layout',pdf,txt],check=True,timeout=40)
   lines=Path(txt).read_text(encoding='utf-8',errors='replace').splitlines()
   report.update({'text_lines':len(lines),'text_chars':sum(map(len,lines)),
       'text_preview':[x[:220] for x in lines if x.strip()][:65],
       'status':'TEXT_EXTRACTED_REQUIRE_TABLE_SCHEMA_CHECK'})
   if len(lines)<30:report['status']='HOLD_TEXT_OCR_OR_PARSE'
  except Exception as e:report['status']='HOLD_PDF_TEXT_EXTRACTION';report['error']=str(e)[:300]
 (dest/'supplement_audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
 print(json.dumps({'status':report['status'],'fetched':report['fetched'],'preview':report.get('text_preview',[])[:10]}))
if __name__=='__main__':main()