#!/usr/bin/env Rscript

## Step 1 (UNRELATED samples): fast GLMM-FREE null model for the cross-ancestry rg / h2
## RHE estimator. Produces a step1-compatible <out>.rda (modglmm) + <out>.varianceRatio.txt
## WITHOUT fitting a GRM/mixed model -- valid when samples are unrelated (tau_g = 0, so the
## working covariance Sigma is diagonal, the variance ratio is 1, and the RHE reduces to
## ordinary Haseman-Elston; see METHODS_RHE_in_SAIGE.tex Prop. 1).
##
## The output .rda has the SAME structure step2 (SAIGE.Admixed --estimate_cross_anc_rg=TRUE)
## expects, so the EXISTING, validated step2 + step3_crossAncestryRgRHE.R run unchanged:
##   quantitative: OLS fit (optional rank-inverse-normal of y); theta=(1/var(res), 0).
##   binary:       logistic GLM fit; theta=(1, 0); RHE uses the PCGC residual internally.
##
## Relatedness: this step assumes the input is ALREADY unrelated. Provide --sampleFile with
## the unrelated IIDs (e.g. plink2 --king-cutoff 0.0884 .in.id), OR pre-restrict the data.
## NOTE: do NOT use this for related samples -- use the standard step1_fitNULLGLMM.R there.

suppressPackageStartupMessages({
  library(optparse)
  library(methods)
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

stopifnot(nzchar(opt$phenoFile), nzchar(opt$phenoCol), nzchar(opt$outputPrefix))
trait <- match.arg(opt$traitType, c("quantitative", "binary"))
invNorm <- toupper(opt$invNormalize) %in% c("TRUE", "T", "1")

## ---- ScoreTest_NULL_Model (GRM-free covariate projection; inlined from SAIGE) ----
ScoreTest_NULL_Model <- function(mu, mu2, y, X) {
  V <- as.vector(mu2)
  res <- as.vector(y - mu)
  XV <- t(X * V)
  XVX <- t(X) %*% t(XV)
  XVX_inv <- solve(XVX)
  XXVX_inv <- X %*% XVX_inv
  XVX_inv_XV <- XXVX_inv * V
  S_a <- colSums(X * res)
  re <- list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv,
             S_a = S_a, XVX_inv_XV = XVX_inv_XV, V = V)
  class(re) <- "SA_NULL"
  re
}

inv_normal <- function(y) {
  qnorm((rank(y) - 0.5) / length(y))
}

## ---- read phenotype + covariates ----
ph <- read.table(opt$phenoFile, header = TRUE, stringsAsFactors = FALSE,
                 sep = "\t", check.names = FALSE)
idcol <- opt$sampleIDColinphenoFile
stopifnot(idcol %in% colnames(ph), opt$phenoCol %in% colnames(ph))
if (nzchar(opt$sampleFile)) {
  keep <- readLines(opt$sampleFile)
  keep <- keep[nzchar(keep)]
  ph <- ph[ph[[idcol]] %in% keep, , drop = FALSE]
  cat(sprintf("[step1-noGRM] restricted to %d unrelated samples from %s\n",
              nrow(ph), opt$sampleFile))
}
covars <- if (nzchar(opt$covarColList))
  strsplit(opt$covarColList, ",")[[1]] else character(0)
covars <- trimws(covars)
miss <- setdiff(covars, colnames(ph))
if (length(miss)) stop("covariates not found in phenoFile: ", paste(miss, collapse = ", "))

ph <- ph[stats::complete.cases(ph[, c(opt$phenoCol, covars, idcol)]), , drop = FALSE]
sampleID <- as.character(ph[[idcol]])
N <- nrow(ph)
y <- as.numeric(ph[[opt$phenoCol]])
X <- cbind(Intercept = 1, as.matrix(ph[, covars, drop = FALSE]))
storage.mode(X) <- "double"
cat(sprintf("[step1-noGRM] N=%d trait=%s covars=[%s] invNorm=%s\n",
            N, trait, paste(covars, collapse = ","), invNorm))

## ---- GRM-free null fit ----
if (trait == "quantitative") {
  if (invNorm) y <- inv_normal(y)
  alpha <- as.vector(solve(t(X) %*% X, t(X) %*% y))
  eta <- as.vector(X %*% alpha)
  mu <- eta
  res <- y - mu
  # theta[1]=1 (tau_g=0): with the variance ratio = 1 this makes the RHE score exactly
  # g~'res (NO tau/vr divisor), the form consistent with q_e = res'Sigma^-1 res = res.res
  # under Sigma^-1 = diag(mu2*tau1) = I. (h2/rg are scale-invariant in res; the ABSOLUTE
  # tau1 is irrelevant -- only tau1*vr = 1 matters, else the ill-conditioned multi-GRM HE
  # collapses sigma2_anc. Verified: tau1*vr != 1 shifts compScoreSum by (tau1*vr)^2.)
  tau <- c(1, 0)
  mu2 <- rep(1, N)                                          # = 1/tau1; Sigma^-1 diag = mu2*tau1 = 1
  obj_cc <- NULL
} else {                                                    # binary
  fit <- stats::glm.fit(X, y, family = stats::binomial())
  alpha <- as.vector(fit$coefficients)
  eta <- as.vector(X %*% alpha)
  mu <- as.vector(stats::plogis(eta))
  res <- y - mu
  mu2 <- mu * (1 - mu)
  tau <- c(1, 0)
  obj_cc <- tryCatch({
    cc <- SKAT::SKAT_Null_Model(y ~ X - 1, out_type = "D", Adjustment = FALSE)
    cc$mu <- mu; cc$res <- res; cc$pi_1 <- mu2; cc
  }, error = function(e) list(mu = mu, res = res, pi_1 = mu2, res.out = NULL))
}

objnoK <- ScoreTest_NULL_Model(mu, mu2, y, X)

## ---- assemble modglmm in the step1 output format (LOCO=FALSE, no GRM) ----
modglmm <- list(
  theta = tau,
  coefficients = alpha,
  linear.predictors = eta,
  fitted.values = mu,
  Y = NULL,
  residuals = res,
  cov = NULL,
  converged = TRUE,
  sampleID = sampleID,
  obj.noK = objnoK,
  y = y,
  X = X,
  traitType = trait,
  isCovariateOffset = FALSE,
  LOCO = FALSE,
  obj.glm.null = NULL,
  offset = rep(0, N),
  useSparseGRMtoFitNULL = FALSE
)
if (!is.null(obj_cc)) modglmm$obj_cc <- obj_cc

save(modglmm, file = paste0(opt$outputPrefix, ".rda"))
## variance ratio = 1 (no GRM): "value <category> <index>" matching Get_Variance_Ratio
writeLines("1\tnull\t1", paste0(opt$outputPrefix, ".varianceRatio.txt"))

cat(sprintf("[step1-noGRM] wrote %s.rda (theta=%.4f,%.0f) and %s.varianceRatio.txt (=1)\n",
            opt$outputPrefix, tau[1], tau[2], opt$outputPrefix))
