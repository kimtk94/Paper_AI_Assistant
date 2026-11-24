#!/usr/bin/env Rscript
# GSE213761: exploratory baseline response-group and within-subject interaction models.
# These are observational anti-TNF samples, NOT secukinumab efficacy or prospective prediction.
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Usage: Rscript response_interaction.R output_dir')
out <- args[[1]]
counts_file <- file.path(out,'hs_therapy_paired_counts.csv')
meta_file <- file.path(out,'hs_therapy_paired_metadata.csv')
if(!file.exists(counts_file) || !file.exists(meta_file)) stop('HOLD_MISSING_VERIFIED_COUNT_AND_METADATA')
mat <- read.csv(counts_file,row.names=1,check.names=FALSE)
meta <- read.csv(meta_file,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(meta)==28,ncol(mat)==28,identical(as.character(meta$column),colnames(mat)),
          length(unique(meta$patient))==14,all(meta$responder %in% c('NO','YES')),
          all(meta$visit %in% c('V1','V12')),all(table(meta$patient,meta$visit)==1))
stopifnot(length(unique(meta$patient[meta$responder=='YES']))==7,
          length(unique(meta$patient[meta$responder=='NO']))==7)
stopifnot(all(vapply(split(meta$responder,meta$patient),function(x)length(unique(x))==1,logical(1))))
mat<-as.matrix(mat)
stopifnot(all(is.finite(mat)),all(mat>=0),max(abs(mat-round(mat)))<1e-7)
storage.mode(mat)<-'integer'
rownames(meta)<-meta$column
meta$patient <- factor(meta$patient)
meta$post <- ifelse(meta$visit=='V12',1,0)
meta$response_post <- meta$post*as.integer(meta$responder=='YES')
meta$responder <- factor(meta$responder,levels=c('NO','YES'))
keep<-rowSums(mat>=10)>=4
if(sum(keep)<5000)stop('HOLD_TOO_FEW_GENES_AFTER_QC')
# A responder main effect is entirely absorbed by patient fixed effects.
# post: change in the NO group; response_post: extra change in YES vs NO.
dds<-DESeqDataSetFromMatrix(mat[keep,,drop=FALSE],colData=meta,design=~patient+post+response_post)
if(qr(model.matrix(design(dds),as.data.frame(colData(dds))))$rank!=ncol(model.matrix(design(dds),as.data.frame(colData(dds)))))stop('HOLD_RANK_DEFICIENT_MODEL')
dds<-DESeq(dds,quiet=TRUE)
res <- results(dds,name='response_post',alpha=.05)
res_post <- results(dds,name='post',alpha=.05)
write_results<-function(r,fn) {
  z<-as.data.frame(r);z$entrez_gene_id<-rownames(z)
  z<-z[,c('entrez_gene_id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
  z<-z[order(z$padj,z$pvalue,na.last=TRUE),]
  write.csv(z,file.path(out,fn),row.names=FALSE)
  z
}
interaction<-write_results(res,'hs_response_time_interaction_deseq2.csv')
noresponse_change<-write_results(res_post,'hs_nonresponders_time_change_deseq2.csv')
# Baseline association with *future* response is not validated prediction.
pre<-meta$visit=='V1'
base_meta<-data.frame(responder=factor(meta$responder[pre],levels=c('NO','YES')),
                      row.names=as.character(meta$column[pre]))
base_mat<-mat[,pre,drop=FALSE]
base_keep<-rowSums(base_mat>=10)>=3
if(sum(base_keep)<5000)stop('HOLD_BASELINE_LOW_GENE_COUNT')
bdds<-DESeqDataSetFromMatrix(base_mat[base_keep,,drop=FALSE],base_meta,design=~responder)
bdds<-DESeq(bdds,quiet=TRUE)
baseline<-write_results(results(bdds,contrast=c('responder','YES','NO'),alpha=.05),
                        'hs_pre_response_group_deseq2.csv')
ids<-c(IL17A='3605',IL1B='3553',TNF='7124',CXCL8='3576',CXCL13='10563',
       S100A8='6279',S100A9='6280',KRT17='3872',CCL20='6364',
       MS4A1='931',CD79A='973',JCHAIN='3512',CFH='3075')
summary_targets <- data.frame(symbol=names(ids),entrez_gene_id=unname(ids),stringsAsFactors=FALSE)
for(name in c('interaction','baseline')){
  z <- if(name=='interaction')interaction else baseline
  m<-z[match(summary_targets$entrez_gene_id,z$entrez_gene_id),c('baseMean','log2FoldChange','pvalue','padj')]
  colnames(m)<-paste0(name,'_',colnames(m));summary_targets<-cbind(summary_targets,m)
}
write.csv(summary_targets,file.path(out,'hs_targeted_interaction_baseline.csv'),row.names=FALSE)
# Only selected predeclared genes in normalized-count table (NOT raw paired matrix in artifact)
norm<-counts(dds,normalized=TRUE)
ix<-match(unname(ids),rownames(norm));ok<-which(!is.na(ix));
tn<-as.data.frame(log2(norm[ix[ok],,drop=FALSE]+1))
tn$symbol<-names(ids)[ok];tn$entrez_gene_id<-unname(ids)[ok]
tn<-tn[,c('symbol','entrez_gene_id',colnames(mat))]
write.csv(tn,file.path(out,'hs_targets_log2_normalized_counts.csv'),row.names=FALSE)
metrics<-data.frame(metric=c('patients','responder_patients','nonresponder_patients','lesional_biopsies',
 'interaction_genes_tested','interaction_FDR_0.05','interaction_FDR_0.05_absLFC_gt1',
 'baseline_genes_tested','baseline_FDR_0.05','baseline_FDR_0.05_absLFC_gt1'),
 value=c(14,7,7,28,nrow(interaction),sum(interaction$padj<.05,na.rm=TRUE),
 sum(interaction$padj<.05&abs(interaction$log2FoldChange)>1,na.rm=TRUE),
 nrow(baseline),sum(baseline$padj<.05,na.rm=TRUE),
 sum(baseline$padj<.05&abs(baseline$log2FoldChange)>1,na.rm=TRUE)))
write.csv(metrics,file.path(out,'hs_response_model_summary.csv'),row.names=FALSE)
writeLines(c('# GSE213761 response-stratified analysis (exploratory)',
  'Design interaction = patient-fixed-effect + numeric(post) + numeric(post x future responder)',
  'This test compares DIFFERENCE IN TIME CHANGES between responder/nonresponder groups.',
  'Baseline = YES vs NO across 14 independent people. No random split or ML prediction.',
  'No untreated contemporaneous controls; responder definition from GEO field not independently adjudicated.',
  'Drug = Adalimumab; cannot infer Secukinumab response. Small n=7/group; confounding possible.',
  'FDR only within each genome-wide model. Targets are exploratory and not validated biomarkers.',
  '',paste(metrics$metric,metrics$value,sep=': ',collapse='\n')),
  file.path(out,'RESPONSE_MODEL_QC.md'))
print(metrics)