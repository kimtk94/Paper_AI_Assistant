#!/usr/bin/env python3
"""Fetch original article Supplementary Table 2, extract evidence without fabricating 406-gene overlap."""
import argparse, hashlib, io, json, os, re, shutil, subprocess, time, urllib.request, zipfile
from pathlib import Path
from collections import Counter

ARCHIVE_URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12568438/supplementaryFiles"
MEMBER="ymj-66-780-s002.pdf"
MAX_ARCHIVE=50_000_000
MAX_PDF=2_000_000

def pdf_from_archive(content):
    if len(content)>MAX_ARCHIVE or not zipfile.is_zipfile(io.BytesIO(content)):
        raise ValueError("invalid or oversized official supplementary ZIP")
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        matches=[x for x in archive.infolist() if x.filename.rsplit("/",1)[-1]==MEMBER]
        if len(matches)!=1:
            raise ValueError("official Table 2 member missing or ambiguous")
        if matches[0].file_size>MAX_PDF:
            raise ValueError("official Table 2 PDF too large")
        data=archive.read(matches[0])
    if not data.startswith(b"%PDF-") or len(data)<30000:
        raise ValueError("official Table 2 member is not a valid-sized PDF")
    return data

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
 for url in [ARCHIVE_URL]+URLS:
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 (scientific supplementary data audit)'})
   with urllib.request.urlopen(req,timeout=35) as resp:
    limit=MAX_ARCHIVE if url==ARCHIVE_URL else MAX_PDF
    content=resp.read(limit+1)
   if len(content)>limit:raise ValueError('source payload exceeds size limit')
   if url==ARCHIVE_URL:
    report['official_archive_sha256']=hashlib.sha256(content).hexdigest()
    content=pdf_from_archive(content)
   if not content.startswith(b'%PDF-') or len(content)<30000:raise ValueError('non-PDF or unexpectedly short file')
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
       'status':'TEXT_EXTRACTED_REQUIRE_TABLE_SCHEMA_CHECK'})
   if len(lines)<30:report['status']='HOLD_TEXT_OCR_OR_PARSE'
  except Exception as e:report['status']='HOLD_PDF_TEXT_EXTRACTION';report['error']=str(e)[:300]
 (dest/'supplement_audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
 print(json.dumps({'status':report['status'],'fetched':report['fetched']}))
if __name__=='__main__':main()