#!/usr/bin/env Rscript
# Read count matrix 4 vs 4: sample leave-one-out sensitivity, DESeq2 refit.
# This is nonpaired, observational tissue comparison. No patient-level response variable.
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) stop('Usage Rscript loo_deseq2.R count_matrix.csv outdir')
ct <- read.csv(args[1],row.names=1,check.names=FALSE)
out <- args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
ids <- c('GSM8671652','GSM8671653','GSM8671654','GSM8671655','GSM8671656','GSM8671657','GSM8671658','GSM8671659')
stopifnot(identical(colnames(ct),ids),ncol(ct)==8,nrow(ct)>40000)
ct <- as.matrix(ct);stopifnot(all(is.finite(ct)),all(ct>=0),max(abs(ct-round(ct)))<1e-7)
storage.mode(ct) <- 'integer'
group <- factor(c('untreated','dienogest','untreated','dienogest','dienogest','untreated','dienogest','untreated'),levels=c('untreated','dienogest'))
meta <- data.frame(group=group,row.names=ids)
targets <- c('HLA-DOB','ELOVL2','NSFP1','GSTM1','LRRC15')
base_file <- file.path(out,'deseq2_author_filter_gene_results.csv')
stopifnot(file.exists(base_file))
baseline <- read.csv(base_file,check.names=FALSE)
colnames(baseline)[1] <- 'gene'
baseline_lfc <- setNames(baseline$log2FoldChange,baseline$gene)
rows <- list(); fold_deg <- list()
for (j in seq_along(ids)) {
  id <- ids[j];keep_sample <- ids!=id; subm <- ct[,keep_sample,drop=FALSE];submeta <- meta[keep_sample,,drop=FALSE]
  keep_gene <- rowSums(subm[,submeta$group=='untreated',drop=FALSE]>0)>=1 &
               rowSums(subm[,submeta$group=='dienogest',drop=FALSE]>0)>=1
  stopifnot(sum(keep_gene)>10000)
  dds <- DESeqDataSetFromMatrix(countData=subm[keep_gene,,drop=FALSE],colData=submeta,design=~group)
  dds <- DESeq(dds,quiet=TRUE)
  rr <- results(dds,contrast=c('group','dienogest','untreated'),alpha=.05)
  dat <- as.data.frame(rr);dat$gene<-rownames(dat)
  sig <- dat$gene[!is.na(dat$padj)&dat$padj<.05&abs(dat$log2FoldChange)>1]
  fold_deg[[id]] <- sig
  for (gene in targets) {
    row <- dat[dat$gene==gene,,drop=FALSE]
    rows[[length(rows)+1]] <- data.frame(excluded_gsm=id,excluded_group=as.character(meta[id,'group']),gene=gene,
       n_untreated=sum(submeta$group=='untreated'),n_dienogest=sum(submeta$group=='dienogest'),
       log2fc=if(nrow(row))row$log2FoldChange else NA_real_,
       padj=if(nrow(row))row$padj else NA_real_,
       direction_baseline=if(!is.na(baseline_lfc[gene]))sign(baseline_lfc[gene]) else NA_integer_,
       direction_loo=if(nrow(row))sign(row$log2FoldChange) else NA_integer_,
       significant=if(nrow(row))(!is.na(row$padj)&row$padj<.05&abs(row$log2FoldChange)>1) else FALSE,
       filtered_n=sum(keep_gene),total_deg_n=length(sig))
  }
  write.csv(data.frame(gene=sig),file.path(out,paste0('LOO_DEG_',id,'.csv')),row.names=FALSE)
  cat('LOO',id,'filtered',sum(keep_gene),'sig_n',length(sig),'\n')
}
tab <- do.call(rbind,rows)
write.csv(tab,file.path(out,'loo_target_genes.csv'),row.names=FALSE)
summary <- do.call(rbind,lapply(split(tab,tab$gene),function(s){
  data.frame(gene=s$gene[1],n_loo=length(s$excluded_gsm),
    n_direction_stable=sum(s$direction_baseline==s$direction_loo,na.rm=TRUE),
    n_significant=sum(s$significant,na.rm=TRUE),
    min_abs_lfc=min(abs(s$log2fc),na.rm=TRUE),
    min_lfc=min(s$log2fc,na.rm=TRUE),max_lfc=max(s$log2fc,na.rm=TRUE),
    min_padj=min(s$padj,na.rm=TRUE),max_padj=max(s$padj,na.rm=TRUE))
}))
write.csv(summary,file.path(out,'loo_target_summary.csv'),row.names=FALSE)
cat('Overall DEGs per LOO range:',min(tab$total_deg_n),max(tab$total_deg_n),'\n')
print(summary)
writeLines(c('# GSE283755 DESeq2 Leave-One-Out sensitivity audit',
  paste('DESeq2 version',as.character(packageVersion('DESeq2'))),
  '8 independent patients, group counts 3 vs 4 in each LOO; patient-level independent unit.',
  'Baseline 421 DEGs used author-like count filter; each LOO re-estimates dispersion/size factors.',
  'A direction-stable gene can still fail FDR under n=7; external replication is still required.',
  'This is not a treatment-response biomarker validation.'),file.path(out,'LOO_REPORT.md'))