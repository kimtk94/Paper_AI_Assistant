#!/usr/bin/env Rscript
# GSE294009 full nuisance-reestimated leave-one-independent-donor-out
# and explicit slide effects. Research sensitivity only, not exact paper DEG.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: source_dir preflight_dir out_dir")
source_dir <- args[[1L]]; preflight_dir <- args[[2L]]; out_dir <- args[[3L]]
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
for(pkg in c("standR","edgeR","limma","S4Vectors","SummarizedExperiment","jsonlite","statmod")){
  if(!requireNamespace(pkg,quietly=TRUE))stop("MISSING_R_PACKAGE:",pkg)
}
set.seed(20261122)
audit <- jsonlite::read_json(file.path(source_dir,"HS_132_FULL_WTA_SOURCE_AUDIT.json"),simplifyVector=TRUE)
pre <- jsonlite::read_json(file.path(preflight_dir,"HS_132_STANDR_MODEL_PREFLIGHT.json"),simplifyVector=TRUE)
if(!identical(audit$status,"FULL_WTA_PUBLIC_132_HS_SOURCE_COUNTS_VERIFIED__NO_DE") ||
   !identical(pre$status,"PRIMARY_FFPE_SOURCE_AND_DESIGN_READY__CONFIRMATORY_STAGE_DE_HOLD"))
  stop("SOURCE_PROVENANCE_STOP")
if(audit$endogenous_gene_rows!=18676 || audit$included_S2_GEO_HS_ROI!=132 ||
   audit$n_independent_HS_donors!=15)stop("SOURCE_DIMENSIONS_STOP")
tab <- read.csv(gzfile(file.path(source_dir,"HS_132_HS_GENE_WTA_COUNTS.csv.gz")),
                check.names=FALSE,stringsAsFactors=FALSE)
meta <- read.csv(file.path(source_dir,"HS_132_AUTHOR_S2_MANIFEST.csv"),
                 check.names=FALSE,stringsAsFactors=FALSE)
if(nrow(tab)!=18676 || nrow(meta)!=132 || anyDuplicated(tab$gene) ||
   !identical(names(tab)[-1L],meta$gsm))stop("SOURCE_JOIN_STOP")
x <- as.matrix(tab[,-1L,drop=FALSE]);rownames(x) <- tab$gene;storage.mode(x)<-"integer"
if(anyNA(x) || any(x<0L))stop("COUNT_CHECK_STOP")
keep <- meta$site=="Yonsei University" &
        meta$preparation=="Formalin-Fixed Paraffin-Embedded slide" &
        meta$donor!="HS12"
md0 <- meta[keep,,drop=FALSE]; x0 <- x[,keep,drop=FALSE]
if(nrow(md0)!=82L || length(unique(md0$donor))!=10L ||
   sum(md0$stage_s2=="Early")!=57L || sum(md0$stage_s2=="Late")!=25L)
  stop("PRIMARY_COHORT_STOP")
md0$Stage <- factor(md0$stage_s2,levels=c("Early","Late"))
md0$Donor <- factor(md0$donor)
md0$SlideName <- factor(md0$slide_label_not_donor)
md0$Group <- factor(paste(md0$stage_s2,md0$compartment,sep="__"))
md0$Compartment <- factor(md0$compartment)
deep_name <- "Inflammatory Aggregate_DeepDermis"
all_donors <- sort(unique(as.character(md0$Donor[md0$compartment==deep_name])))
if(length(all_donors)!=10L)stop("DEEP_DONOR_STOP")
if(!"JCHAIN" %in% rownames(x0))stop("JCHAIN_SOURCE_STOP")

# Fit cohort-specific standR normalization, background controls and RUV4.
# No data from withheld donor is passed to upstream transformations.
fit_ruv <- function(submat,submeta,norm){
  dge <- edgeR::DGEList(counts=submat)
  spe <- standR::readGeoMxFromDGE(dge)
  for(field in c("Stage","Donor","SlideName","Group","Compartment"))
    SummarizedExperiment::colData(spe)[[field]] <- submeta[[field]]
  spe <- standR::geomxNorm(spe,method=norm,log=TRUE)
  if(length(SummarizedExperiment::assays(spe))<2L)stop("NORMALIZED_ASSAY_MISSING")
  spe <- standR::findNCGs(spe,n_assay=2L,batch_name="SlideName",top_n=200L)
  ncs <- S4Vectors::metadata(spe)$NCGs
  if(length(ncs)<100L || !all(ncs %in% rownames(spe)) ||
     "JCHAIN" %in% ncs)stop("NCG_QUALITY_STOP")
  spe <- standR::geomxBatchCorrection(spe,k=2L,factors="Group",
                                     NCGs=ncs,n_assay=2L,method="RUV4")
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  if(!all(c("ruv_W1","ruv_W2") %in% names(cd)))stop("RUV4_FACTORS_MISSING")
  if(any(!is.finite(as.matrix(cd[,c("ruv_W1","ruv_W2"),drop=FALSE]))))
    stop("RUV4_FACTORS_NONFINITE")
  if(nrow(cd)!=nrow(submeta))stop("META_ALIGNMENT_STOP")
  list(coldata=cd,n_controls=length(ncs),
       ruv=as.matrix(cd[,c("ruv_W1","ruv_W2"),drop=FALSE]))
}
fit_voom <- function(submat,submeta,ruv_cd,norm,model){
  sel <- submeta$compartment==deep_name
  m <- droplevels(submeta[sel,,drop=FALSE])
  if(nrow(m)!=sum(sel))stop("MODEL_JOIN_STOP")
  early_n <- length(unique(as.character(m$Donor[m$Stage=="Early"])))
  late_n <- length(unique(as.character(m$Donor[m$Stage=="Late"])))
  if(early_n<3L || late_n<2L)stop("INSUFFICIENT_INDEPENDENT_DONORS")
  cd <- ruv_cd[sel,,drop=FALSE]
  cd$Stage <- factor(as.character(m$Stage),levels=c("Early","Late"))
  cd$SlideName <- droplevels(m$SlideName)
  formula <- if(identical(model,"Stage_RUV"))
                ~0 + Stage + ruv_W1 + ruv_W2
             else if(identical(model,"Stage_Slide_RUV"))
                ~0 + Stage + SlideName + ruv_W1 + ruv_W2
             else stop("UNRECOGNIZED_DESIGN")
  d <- stats::model.matrix(formula,data=cd)
  # limma::makeContrasts requires syntactically valid coefficient names.
  # Source slide labels such as "HS slide 2 plz" are legitimate categories
  # but produce invalid coefficient names. Rename DESIGN COLUMNS only;
  # never mutate the source slide labels, batch grouping, or donor metadata.
  source_design_names <- colnames(d)
  colnames(d) <- make.names(source_design_names, unique=TRUE)
  if (!all(c("StageEarly","StageLate") %in% colnames(d)) ||
      anyDuplicated(colnames(d))) stop("CONTRAST_DESIGN_NAMES_STOP")
  if(qr(d)$rank!=ncol(d))stop("SLIDE_OR_RUV_STAGE_ALIASED")
  if(ncol(d)>=nrow(m)-3L)stop("INSUFFICIENT_RESIDUAL_DF")
  ng <- length(unique(as.character(m$Donor)))
  if(ng!=early_n+late_n)stop("DONOR_STAGE_NESTING_ERROR")
  ec <- edgeR::DGEList(counts=submat[,sel,drop=FALSE])
  ec <- edgeR::calcNormFactors(ec,method=norm)
  genes <- edgeR::filterByExpr(ec,design=d)
  gix <- match("JCHAIN",rownames(ec))
  if(is.na(gix) || !genes[gix])stop("JCHAIN_NOT_DETECTED_AFTER_FILTER")
  if(sum(genes)<500L)stop("EXPRESSED_FEATURES_TOO_FEW")
  ec <- ec[genes,,keep.lib.sizes=FALSE]
  # Re-estimate voom and within-patient correlation from fold-specific input.
  v0 <- limma::voom(ec,design=d,plot=FALSE)
  rho <- limma::duplicateCorrelation(v0,design=d,block=m$Donor)$consensus.correlation
  if(length(rho)!=1L || !is.finite(rho) || abs(rho)>=1)stop("UNSTABLE_DONOR_CORRELATION")
  v1 <- limma::voom(ec,design=d,block=m$Donor,correlation=rho,plot=FALSE)
  fit <- limma::lmFit(v1,design=d,block=m$Donor,correlation=rho)
  cmat <- limma::makeContrasts(StageLate-StageEarly,levels=d)
  z <- limma::contrasts.fit(fit,cmat)
  # Our reported endpoint is the exploratory coefficient and not a p value.
  coef <- as.numeric(z$coefficients["JCHAIN",1L])
  if(length(coef)!=1L || !is.finite(coef))stop("JCHAIN_NONFINITE_LOGFC")
  list(coef=coef,rho=as.numeric(rho),n_gene=sum(genes),
       n_early_donors=early_n,n_late_donors=late_n,
       n_deep_roi=nrow(m),n_slides=length(unique(as.character(m$SlideName))),
       rank=qr(d)$rank,n_cols=ncol(d))
}
result <- list()
model_modes <- c("Stage_RUV","Stage_Slide_RUV")
normalizations <- c("TMM","upperquartile")
for(norm in normalizations){
  norm_out <- list()
  for(held in c("NONE",all_donors)){
    train <- if(held=="NONE") rep(TRUE,nrow(md0)) else as.character(md0$Donor)!=held
    sm <- droplevels(md0[train,,drop=FALSE])
    sx <- x0[,train,drop=FALSE]
    if(held!="NONE" && held %in% as.character(sm$Donor))stop("DONOR_LEAK_IN_FOLD")
    if(!identical(colnames(sx),sm$gsm))stop("FOLD_SAMPLE_ID_MISMATCH")
    # Complete nuisance re-estimation for this 81/74/etc-ROI fold.
    upstream <- tryCatch(fit_ruv(sx,sm,norm),error=function(e)e)
    for(mode in model_modes){
      label <- if(held=="NONE")"FULL" else paste0("WITHOUT_",held)
      if(inherits(upstream,"error")){
        val <- list(status="HOLD_RUV_ERROR",reason=conditionMessage(upstream),
                    donor_excluded=held,normalization=norm,model=mode,
                    n_upstream_roi=nrow(sm),full_nuisance_reestimation=TRUE)
      }else{
        fit <- tryCatch(fit_voom(sx,sm,upstream$coldata,norm,mode),
                      error=function(e)e)
        if(inherits(fit,"error")){
          val <- list(status="HOLD_MODEL_ERROR",reason=conditionMessage(fit),
                      donor_excluded=held,normalization=norm,model=mode,
                      n_upstream_roi=nrow(sm),n_NCG=upstream$n_controls,
                      full_nuisance_reestimation=TRUE)
        }else{
          val <- c(list(status="PASS_MODEL_FIT",donor_excluded=held,
                        normalization=norm,model=mode,
                        n_upstream_roi=nrow(sm),n_NCG=upstream$n_controls,
                        full_nuisance_reestimation=TRUE),fit)
        }
      }
      norm_out[[paste(label,mode,sep="__")]] <- val
      cat("FULL_REFIT",norm,mode,held,val$status,
          if(!is.null(val$coef))format(round(val$coef,6),nsmall=6) else val$reason,
          "\n")
    }
  }
  result[[norm]] <- norm_out
}
summaries <- list()
for(norm in normalizations){
  for(mode in model_modes){
    rows <- unname(Filter(function(t)identical(t$model,mode),result[[norm]]))
    base <- rows[[which(vapply(rows,function(t)t$donor_excluded=="NONE",logical(1)))]]
    drops <- Filter(function(t)t$donor_excluded!="NONE",rows)
    passes <- Filter(function(t)t$status=="PASS_MODEL_FIT",drops)
    v <- vapply(passes,function(t)t$coef,numeric(1))
    all_pass <- length(passes)==10L && identical(base$status,"PASS_MODEL_FIT")
    same <- all_pass && all(sign(v)==sign(base$coef))
    summaries[[paste(norm,mode,sep="__")]] <- list(
      normalization=norm,model=mode,
      base_status=base$status,
      base_effect=if(!is.null(base$coef))base$coef else NULL,
      base_reason=if(!is.null(base$reason))base$reason else NULL,
      n_donors_loo=10L,n_successful_refits=length(passes),
      all_success_and_same_direction=same,
      loo_effect_range=if(length(v))unname(range(v)) else NULL,
      n_HOLD=10L-length(passes),
      note="Every LOO fold reestimates standR normalization, NCG and RUV4 from training donors; not an external validation or holdout prediction")
    cat("REFIT_SUMMARY",norm,mode,"base",base$status,"success",length(passes),
        "/10 all_positive",same,
        "range",if(length(v))paste(round(range(v),5),collapse=",") else "NA","\n")
  }
}
report <- list(
  status="FULL_NUISANCE_REESTIMATED_DONOR_LOO_AND_SLIDE_DESIGN_AUDIT__NOT_PAPER_REPLICATION",
  source_accession="GSE294009",
  source_matrix_sha256=audit$matrix_csv_gz_sha256,
  n_source_HS_roi=132L,n_primary_HS_roi=82L,n_primary_donors=10L,
  n_primary_early_donors=7L,n_primary_late_donors=3L,
  n_primary_deep_roi=24L,
  original_156_ROI_reproduced=FALSE,
  confirmed_DE_biomarker=FALSE,
  paper_equivalent_normalization=FALSE,
  full_reestimation_per_heldout_donor=TRUE,
  slide_fixed_effect_explicitly_tested=TRUE,
  initial_workflow_run_type="RESEARCH_HYPOTHESIS_ONLY",
  summaries=summaries,
  results=result
)
jsonlite::write_json(report,file.path(out_dir,"HS_FULL_REFIT_SLIDE_DONOR_LOO_AUDIT.json"),
                     pretty=TRUE,auto_unbox=TRUE,na="null")
cat("FULL_REFIT_COMPLETE",report$status,"\n")
