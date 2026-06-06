#!/usr/bin/env Rscript

## Step 1 (UNRELATED samples): fast GLMM-FREE null model for the cross-ancestry rg / h2
## RHE estimator. Thin CLI wrapper around SAIGE::fitNULL_noGRM() -- produces a step1-compatible
## <out>.rda (modglmm) + <out>.varianceRatio.txt WITHOUT fitting a GRM/mixed model. Valid when
## samples are unrelated (tau_g = 0, variance ratio = 1, RHE reduces to ordinary Haseman-Elston).
## The output has the SAME structure step2 (--estimate_cross_anc_rg=TRUE) + step3 consume, so
## the existing validated steps run unchanged. Do NOT use for related samples -- use
## step1_fitNULLGLMM.R there; for related cohorts restrict to an unrelated subset via --sampleFile.

suppressPackageStartupMessages({
  library(optparse)
  library(methods)
  library(SAIGE)
})

opt_list <- list(
  make_option("--phenoFile", type = "character", default = ""),
  make_option("--phenoCol", type = "character", default = ""),
  make_option("--covarColList", type = "character", default = "",
              help = "comma-separated covariates; MUST include ancestry PCs/global-anc"),
  make_option("--sampleIDColinphenoFile", type = "character", default = "IID"),
  make_option("--traitType", type = "character", default = "quantitative",
              help = "quantitative or binary"),
  make_option("--invNormalize", type = "character", default = "FALSE",
              help = "rank-inverse-normal transform a quantitative phenotype (TRUE/FALSE)"),
  make_option("--sampleFile", type = "character", default = "",
              help = "optional: restrict to these unrelated IIDs (one per line)"),
  make_option("--outputPrefix", type = "character", default = "")
)
opt <- parse_args(OptionParser(option_list = opt_list))

invNorm <- toupper(opt$invNormalize) %in% c("TRUE", "T", "1")

fitNULL_noGRM(
  phenoFile = opt$phenoFile,
  phenoCol = opt$phenoCol,
  outputPrefix = opt$outputPrefix,
  covarColList = opt$covarColList,
  sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
  traitType = opt$traitType,
  invNormalize = invNorm,
  sampleFile = opt$sampleFile
)
