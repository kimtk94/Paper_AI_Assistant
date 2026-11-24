#!/usr/bin/env Rscript
# Public GSE294009 132-HS-ROI whole-WTA standR RUV4 feasibility.
# Source and donor gates must already PASS. NOT 156-ROI author reproduction.
# No confirmatory DEG / prognostic / treatment response inference.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: run_standr_ruv4_feasibility.R count_dir preflight_dir out_dir")
source_dir <- args[[1L]]
pref_dir <- args[[2L]]
out_dir <- args[[3L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("standR", "edgeR", "limma", "SpatialExperiment",
              "SummarizedExperiment", "S4Vectors", "jsonlite", "statmod")
unavailable <- required[!vapply(required, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
if (length(unavailable) > 0L) stop("R_PACKAGES_MISSING: ", paste(unavailable, collapse = ","))
options(stringsAsFactors = FALSE)
set.seed(20261121)

audit <- jsonlite::read_json(file.path(source_dir, "HS_132_FULL_WTA_SOURCE_AUDIT.json"),
                            simplifyVector = TRUE)
pref <- jsonlite::read_json(file.path(pref_dir, "HS_132_STANDR_MODEL_PREFLIGHT.json"),
                           simplifyVector = TRUE)
if (!identical(audit$status, "FULL_WTA_PUBLIC_132_HS_SOURCE_COUNTS_VERIFIED__NO_DE")) stop("SOURCE_AUDIT_STOP")
if (!identical(pref$status, "PRIMARY_FFPE_SOURCE_AND_DESIGN_READY__CONFIRMATORY_STAGE_DE_HOLD")) stop("DESIGN_PREFLIGHT_STOP")
if (audit$included_S2_GEO_HS_ROI != 132 || audit$n_independent_HS_donors != 15) stop("DONOR_COUNT_STOP")

raw_table <- read.csv(gzfile(file.path(source_dir, "HS_132_HS_GENE_WTA_COUNTS.csv.gz")),
                      check.names=FALSE, stringsAsFactors=FALSE)
meta <- read.csv(file.path(source_dir, "HS_132_AUTHOR_S2_MANIFEST.csv"),
                 check.names=FALSE, stringsAsFactors=FALSE)
if (nrow(meta) != 132 || nrow(raw_table) != 18676) stop("FEATURE_SAMPLE_DIMENSIONS_STOP")
if (!identical(names(raw_table)[-1L], meta$gsm)) stop("GSM_SAMPLE_ORDER_STOP")
if (anyDuplicated(raw_table$gene) || anyDuplicated(meta$gsm)) stop("GENE_OR_GSM_DUPLICATE_STOP")
X <- as.matrix(raw_table[, -1L, drop=FALSE])
rownames(X) <- raw_table$gene
storage.mode(X) <- "integer"
if (anyNA(X) || any(X < 0L)) stop("BAD_NONNEGATIVE_RAW_COUNTS")

idx <- meta$site == "Yonsei University" &
       meta$preparation == "Formalin-Fixed Paraffin-Embedded slide" &
       meta$donor != "HS12"
if (sum(idx) != 82L || length(unique(meta$donor[idx])) != 10L) stop("PRIMARY_INDEPENDENT_DONORS_STOP")
md <- droplevels(meta[idx, , drop=FALSE])
mat <- X[, idx, drop=FALSE]
if (!identical(colnames(mat), md$gsm)) stop("SAME_SITE_GSM_STOP")
md$Stage <- factor(md$stage_s2, levels=c("Early","Late"))
md$Donor <- factor(md$donor)
md$SlideName <- factor(md$slide_label_not_donor)
md$Compartment <- factor(md$compartment)
md$Group <- factor(paste(md$stage_s2, md$compartment, sep="__"))
if (anyNA(md$Stage) || length(unique(md$SlideName)) != 4L) stop("STAGE_OR_SLIDE_STOP")

deep <- md$compartment == "Inflammatory Aggregate_DeepDermis"
deep_md <- droplevels(md[deep,,drop=FALSE])
if (sum(deep) != 24L || length(unique(deep_md$Donor)) != 10L ||
    length(unique(deep_md$Donor[deep_md$Stage == "Early"])) != 7L ||
    length(unique(deep_md$Donor[deep_md$Stage == "Late"])) != 3L) {
  stop("DEEP_INDEPENDENT_PATIENT_STOP")
}
# Stage and donor fixed effects are structurally nonidentifiable (stage nested within donor).
fixed <- model.matrix(~0 + Donor + Stage, data=deep_md)
if (qr(fixed)$rank == ncol(fixed)) stop("EXPECTED_DONOR_FIXED_STAGE_ALIAS_NOT_FOUND")

dge_all <- edgeR::DGEList(counts=mat)
spe_template <- standR::readGeoMxFromDGE(dge_all)
SummarizedExperiment::colData(spe_template)$Stage <- md$Stage
SummarizedExperiment::colData(spe_template)$Donor <- md$Donor
SummarizedExperiment::colData(spe_template)$SlideName <- md$SlideName
SummarizedExperiment::colData(spe_template)$Compartment <- md$Compartment
SummarizedExperiment::colData(spe_template)$Group <- md$Group

results <- list()
for (norm_method in c("TMM","upperquartile")) {
  cat("BEGIN_STANDR", norm_method, "N_ROI", ncol(spe_template),
      "N_GENES", nrow(spe_template), "\n")
  spe_norm <- standR::geomxNorm(spe_template, method=norm_method, log=TRUE)
  if (length(SummarizedExperiment::assays(spe_norm)) < 2L) stop("NORM_ASSAY_ABSENT")
  spe_ncg <- standR::findNCGs(spe_norm, n_assay=2L,
                              batch_name="SlideName", top_n=200L)
  ncgs <- S4Vectors::metadata(spe_ncg)$NCGs
  if (length(ncgs) < 100L || !all(ncgs %in% rownames(spe_ncg))) stop("INVALID_NCG_SELECTION")
  if ("JCHAIN" %in% ncgs) stop("JCHAIN_AS_NEGATIVE_CONTROL_GENE_STOP")
  spe_ruv <- standR::geomxBatchCorrection(spe_ncg, k=2L, factors="Group",
                    NCGs=ncgs, n_assay=2L, method="RUV4")
  cd <- as.data.frame(SummarizedExperiment::colData(spe_ruv))
  required_W <- c("ruv_W1","ruv_W2")
  if (!all(required_W %in% names(cd))) stop("MISSING_STANDR_RUV4_COVARIATES")
  if (any(!is.finite(as.matrix(cd[,required_W,drop=FALSE])))) stop("NONFINITE_RUV_COVARIATES")
  deep_cd <- cd[deep,,drop=FALSE]
  deep_cd$Stage <- factor(deep_cd$Stage, levels=c("Early","Late"))
  design <- stats::model.matrix(~0 + Stage + ruv_W1 + ruv_W2, data=deep_cd)
  if (qr(design)$rank != ncol(design)) stop("RUV_STAGE_DESIGN_RANK_DEFICIENT")

  deep_dge <- edgeR::DGEList(counts=mat[,deep,drop=FALSE])
  deep_dge <- edgeR::calcNormFactors(deep_dge, method=norm_method)
  keep <- edgeR::filterByExpr(deep_dge, design=design)
  if (sum(keep)<500L) stop("TOO_FEW_EXPRESSED_GENES_FOR_MODEL")
  deep_dge <- deep_dge[keep,,keep.lib.sizes=FALSE]

  voom0 <- limma::voom(deep_dge, design=design, plot=FALSE)
  corr <- limma::duplicateCorrelation(voom0, design=design,
                    block=deep_md$Donor)
  rho <- corr$consensus.correlation
  if (length(rho)!=1L || !is.finite(rho) || rho<= -1 || rho>= 1)
    stop("DONOR_CORRELATION_UNSTABLE")
  voom1 <- limma::voom(deep_dge, design=design,
                       block=deep_md$Donor, correlation=rho, plot=FALSE)
  fit <- limma::lmFit(voom1, design=design,
                      block=deep_md$Donor, correlation=rho)
  contrast <- limma::makeContrasts(StageLate - StageEarly, levels=design)
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2, robust=TRUE)
  if (!all(is.finite(fit2$coefficients))) stop("LIMMA_COEFFICIENTS_NONFINITE")
  gene <- "JCHAIN"
  if (!(gene %in% rownames(fit2$coefficients))) stop("JCHAIN_FILTERED_OUT")
  logfc <- unname(fit2$coefficients[gene,1L])

  # CONDITIONAL sensitivity only: drop one *independent donor* at a time
  # from the downstream voom/limma fitting, leaving the upstream 82-ROI
  # standR normalization/NCG/RUV4 nuisance estimates fixed.
  # Since upstream nuisance estimates used all donors, this IS NOT
  # a fully independent end-to-end donor holdout or paper replication.
  excluded <- unique(as.character(deep_md$Donor))
  loo <- lapply(excluded, function(donor_out) {
    train <- as.character(deep_md$Donor) != donor_out
    late_donors <- length(unique(as.character(deep_md$Donor[train & deep_md$Stage == "Late"])))
    early_donors <- length(unique(as.character(deep_md$Donor[train & deep_md$Stage == "Early"])))
    if (min(early_donors,late_donors) < 2L)
      return(list(withheld_donor=donor_out, status="HOLD_LT2_DONORS",logFC=NULL))
    partial <- tryCatch({
      di <- edgeR::DGEList(counts=mat[,deep,drop=FALSE][,train,drop=FALSE])
      di <- edgeR::calcNormFactors(di, method=norm_method)
      design_i <- stats::model.matrix(~0 + Stage + ruv_W1 + ruv_W2,
                                     data=deep_cd[train,,drop=FALSE])
      if (qr(design_i)$rank != ncol(design_i)) stop("MODEL_RANK_DEFICIENT")
      retained <- edgeR::filterByExpr(di,design=design_i)
      gene_index <- match(gene,rownames(di))
      if (is.na(gene_index) || !retained[gene_index] || sum(retained)<500L)
        stop("JCHAIN_FILTERED_OR_FEATURE_COUNT_LOW")
      di <- di[retained,,keep.lib.sizes=FALSE]
      z0 <- limma::voom(di, design_i, plot=FALSE)
      rho_i <- limma::duplicateCorrelation(z0, design=design_i,
                   block=droplevels(deep_md$Donor[train]))$consensus.correlation
      if (!is.finite(rho_i) || abs(rho_i)>=1) stop("DONOR_CORRELATION_UNSTABLE")
      z1 <- limma::voom(di,design_i,block=droplevels(deep_md$Donor[train]),
                        correlation=rho_i, plot=FALSE)
      fi <- limma::lmFit(z1,design=design_i,block=droplevels(deep_md$Donor[train]),
                         correlation=rho_i)
      fi <- limma::contrasts.fit(fi,limma::makeContrasts(StageLate-StageEarly,levels=design_i))
      ef <- unname(fi$coefficients[gene,1L])
      if (!is.finite(ef)) stop("JCHAIN_EFFECT_NOT_FINITE")
      list(withheld_donor=donor_out,status="PASS",
           early_independent_donors=early_donors,
           late_independent_donors=late_donors,
           effect_log2_Late_over_Early=as.numeric(ef),
           donor_correlation=as.numeric(rho_i))
    },error=function(e)list(withheld_donor=donor_out,status="HOLD_MODEL_ERROR",
                             error=conditionMessage(e),logFC=NULL))
    partial
  })
  passed_loo <- vapply(loo,function(z)identical(z$status,"PASS"),logical(1))
  effects_loo <- vapply(loo[passed_loo],function(z)z$effect_log2_Late_over_Early,numeric(1))
  robust_sign <- length(effects_loo)==length(excluded) &&
                 all(sign(effects_loo)==sign(logfc))
  cat("CONDITIONAL_LOO",norm_method,"PASS",sum(passed_loo),"/",length(excluded),
      "same_sign",robust_sign,
      "range",if(length(effects_loo)) paste(round(range(effects_loo),5),collapse=",") else "NONE","\n")
  result <- list(
    norm_method=norm_method,
    standR_version=as.character(utils::packageVersion("standR")),
    limma_version=as.character(utils::packageVersion("limma")),
    n_WTA_gene_rows=nrow(X),
    n_Yonsei_FFPE_HS_ROI=sum(idx),
    n_Yonsei_FFPE_donors=length(unique(md$donor)),
    n_deep_dermis_ROI=sum(deep),
    n_deep_Early_donors=7L,
    n_deep_Late_donors=3L,
    n_RUV_factors=2L,
    n_control_genes=length(ncgs),
    n_genes_pass_filter=sum(keep),
    donor_blocked_duplicateCorrelation=rho,
    full_design_rank=qr(design)$rank,
    full_design_ncol=ncol(design),
    JCHAIN_exploratory_model_logFC=unname(logfc),
    JCHAIN_conditional_LOO=list(
      description="Leave-one-donor-out downstream voom/limma only, with all-donor-estimated standR normalization and RUV nuisance factors fixed; NOT unbiased external/fully refit validation.",
      n_withheld_independent_donors=length(excluded),
      n_model_success=sum(passed_loo),
      all_valid_and_same_sign=robust_sign,
      effect_range=if (length(effects_loo)) unname(range(effects_loo)) else c(NA_real_,NA_real_),
      withheld_donor_results=loo),
    # No nominal p, adjusted p or biomarker call emitted; source
    # is only a 132-ROI subset with unknown original normalization.
    p_FDR_published=FALSE,
    author_equivalent=FALSE,
    early_late_de_confirmed=FALSE,
    status="STANDR_RUV4_LIMMA_MODEL_FIT_ONLY__INFERENTIAL_CLAIMS_HOLD"
  )
  results[[norm_method]] <- result
  cat("STANDR_PILOT_PASS", norm_method,
      "JCHAIN_logFC", round(logfc,5), "donor_rho", round(rho,5),
      "filtered_genes", sum(keep), "\n")
}
model_report <- list(
  status="TWO_NORMALIZATION_ROBUST_FEASIBILITY_COMPLETED__NOT_PAPER_REPLICATION",
  accession="GSE294009",
  raw_mat_sha256=audit$matrix_csv_gz_sha256,
  donor_stage_mismatch_excluded="HS12",
  source_full_original_156_ROI_available=FALSE,
  original_geomxNorm_method_parameter="UNSPECIFIED",
  primary_and_sensitivity_methods=c("TMM","upperquartile"),
  complete_156_author_FDR_reproduction=FALSE,
  confirmatory_stage_biomarker=FALSE,
  donor_repeat_correlation_checked=TRUE,
  results=results
)
jsonlite::write_json(model_report, file.path(out_dir,"HS_STANDR_RUV4_VOOM_FEASIBILITY.json"),
                     pretty=TRUE,auto_unbox=TRUE)
cat("MODEL_FEASIBILITY_DONE", model_report$status, "\n")
