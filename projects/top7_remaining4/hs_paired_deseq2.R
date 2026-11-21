#!/usr/bin/env Rscript
# Raw Entrez gene count paired lesion analysis: adalimumab visit12 vs baseline.
# Treatment effect is observational; no untreated controls or secukinumab patients.
suppressPackageStartupMessages(library(DESeq2))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=1)stop('Usage: Rscript hs_paired_deseq2.R output_folder')
out<-a[[1]]
ctfile<-file.path(out,'hs_therapy_paired_counts.csv')
metafile<-file.path(out,'hs_therapy_paired_metadata.csv')
if(!file.exists(ctfile)||!file.exists(metafile))stop('HOLD: count or metadata missing: no DE performed')
ct<-read.csv(ctfile,row.names=1,check.names=FALSE)
meta<-read.csv(metafile,row.names=1,check.names=FALSE)
stopifnot(ncol(ct)==28, nrow(meta)==28, identical(colnames(ct),rownames(meta)))
stopifnot(length(unique(meta$patient))==14,
          identical(sort(unique(meta$visit)), c('V1','V12')),
          all(table(meta$patient,meta$visit)==1))
mat<-as.matrix(ct)
stopifnot(all(is.finite(mat)), all(mat>=0), max(abs(mat-round(mat)))<1e-7)
storage.mode(mat)<-'integer'
keep<-rowSums(mat>=10)>=4
if(sum(keep)<7000)stop('Insufficient expressed genes for paired DESeq2')
meta$patient<-factor(meta$patient)
meta$visit<-factor(meta$visit,levels=c('V1','V12'))
dds<-DESeqDataSetFromMatrix(mat[keep,,drop=FALSE],meta,design=~patient+visit)
dds<-DESeq(dds,quiet=TRUE)
res<-results(dds,contrast=c('visit','V12','V1'),alpha=.05)
r<-as.data.frame(res)
r$entrez_gene_id<-rownames(r)
r<-r[,c('entrez_gene_id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
r<-r[order(r$padj,r$pvalue,na.last=TRUE),]
write.csv(r,file.path(out,'hs_therapy_paired_DESeq2_gene_results.csv'),row.names=FALSE)
sig <- (!is.na(r$padj)&r$padj<.05)
summary<-data.frame(metric=c('patients_with_paired_lesional_samples','paired_lesional_samples','genes_total','genes_filtered','FDR_lt_0.05','FDR_lt_0.05_absLFC_gt_1','FDR_up','FDR_down'),value=c(14,28,nrow(mat),sum(keep),sum(sig),sum(sig&abs(r$log2FoldChange)>1),sum(sig&r$log2FoldChange>0),sum(sig&r$log2FoldChange<0)))
write.csv(summary,file.path(out,'hs_therapy_paired_DESeq2_summary.csv'),row.names=FALSE)
# Explore nine pre-selected Entrez IDs, do not claim drug efficacy.
ids<-c(IL17A='3605',IL1B='3553',TNF='7124',CXCL8='3576',CXCL13='10563',S100A8='6279',S100A9='6280',KRT17='3872',CCL20='6364')
g<-r[match(unname(ids),r$entrez_gene_id),,drop=FALSE]
g$symbol<-names(ids)
g<-g[,c('symbol','entrez_gene_id','baseMean','log2FoldChange','pvalue','padj')]
write.csv(g,file.path(out,'hs_therapy_paired_candidate_DESeq2.csv'),row.names=FALSE)
writeLines(c('# GSE213761 real-count patient-paired lesion DESeq2',
  paste('DESeq2 version',packageVersion('DESeq2')),
  'Drug: Adalimumab / anti-TNF; NOT Secukinumab / anti-IL17.',
  'Study: 14 distinct paired patients x 2 visits, lesional biopsy only.',
  'Design: ~patient + visit, Wald contrast V12 vs V1.',
  'There is no untreated parallel cohort: statistical time association is not a causal drug effect.',
  'Effect estimates and BH correction are exploratory; disease severity, drug response and timing covariates remain unadjusted.',
  '',paste(summary$metric,summary$value,sep=': ',collapse='\n')),
  file.path(out,'HS_PAIRED_DESEQ2_REPORT.md'))
print(summary)