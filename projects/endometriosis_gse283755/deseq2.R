#!/usr/bin/env Rscript
# GSE283755 exploratory replication using DESeq2 4 vs 4 unpaired ovarian endometrioma tissue.
# Caveat: source table uses gene-symbol aggregation and is not identical to the original raw feature pipeline.
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop('Usage: Rscript deseq2.R raw_counts_by_gene_symbol.csv output_dir')
infile <- args[[1]]
out <- args[[2]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
ct <- read.csv(infile, row.names=1, check.names=FALSE)
expected <- c('GSM8671652','GSM8671653','GSM8671654','GSM8671655','GSM8671656','GSM8671657','GSM8671658','GSM8671659')
stopifnot(identical(colnames(ct),expected),ncol(ct)==8,nrow(ct)>10000)
mat <- as.matrix(ct)
stopifnot(all(is.finite(mat)),all(mat>=0),max(abs(mat-round(mat)))<1e-7)
storage.mode(mat)<-'integer'
meta <- data.frame(group=factor(c('untreated','dienogest','untreated','dienogest','dienogest','untreated','dienogest','untreated'), levels=c('untreated','dienogest')),row.names=expected)
keep <- rowSums(mat>=10) >= 4
if(sum(keep) < 10000) stop('Too few count features after prefilter')
dds <- DESeqDataSetFromMatrix(countData=mat[keep,,drop=FALSE],colData=meta, design=~group)
dds <- DESeq(dds, quiet=TRUE)
res <- results(dds,contrast=c('group','dienogest','untreated'),alpha=0.05)
resout <- as.data.frame(res)
resout$gene <- rownames(resout)
resout <- resout[,c('gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
resout <- resout[order(resout$padj,resout$pvalue,na.last=TRUE),,drop=FALSE]
write.csv(resout,file.path(out,'deseq2_gene_results.csv'),row.names=FALSE)
summary <- data.frame(criterion=c('raw_genes','prefiltered_genes','padj_below_0.05','padj_below_0.05_and_abs_log2FC_above_1','significant_up','significant_down'),value=c(nrow(mat),sum(keep),sum(!is.na(res$padj)&res$padj<0.05),sum(!is.na(res$padj)&res$padj<0.05 & abs(res$log2FoldChange)>1),sum(!is.na(res$padj)&res$padj<0.05&res$log2FoldChange>1),sum(!is.na(res$padj)&res$padj<0.05&res$log2FoldChange< -1)))
write.csv(summary,file.path(out,'deseq2_summary.csv'),row.names=FALSE)
write.csv(data.frame(gsm=rownames(meta),group=as.character(meta$group),size_factor=sizeFactors(dds)),file.path(out,'deseq2_size_factors.csv'),row.names=FALSE)
writeLines(c('# GSE283755 independent DESeq2 audit',paste('Bioconductor DESeq2 version',as.character(packageVersion('DESeq2'))),'Design ~ group, 4 untreated vs 4 dienogest-exposed; separate patients','Input: GEO Read_Count values summed across duplicate gene symbols, prefiltered genes with >=10 counts in >=4 samples','Not guaranteed to reproduce original 406 DEGs: feature summarization, filtering, package versions, and analysis details may differ.','',paste(summary$criterion, summary$value,sep=': ',collapse='\n')),file.path(out,'DESEQ2_REPORT.md'))
print(summary)