#' Fit a GRM-free (covariate-only) null model for cross-ancestry rg / h2 (UNRELATED samples)
#'
#' Standalone null-model fit for the cross-ancestry genetic-correlation / per-ancestry
#' heritability RHE estimator. Produces a step1-compatible \code{<outputPrefix>.rda}
#' (a \code{modglmm} object) plus \code{<outputPrefix>.varianceRatio.txt}, WITHOUT fitting a
#' GRM / mixed model.
#'
#' This is the correct null for rg/h2 on UNRELATED admixed samples: a GRM-based null
#' (\code{fitNULLGLMM}) fits a genetic variance component \eqn{\tau_g>0} on any heritable
#' trait, which makes the step2 cross-ancestry RHE cross-moment collapse and rg/h2 wrong. With
#' no GRM, \eqn{\tau_g=0}, the working covariance \eqn{\Sigma} is diagonal, the variance ratio
#' is 1, and the RHE reduces to ordinary Haseman-Elston (see METHODS_RHE_in_SAIGE.tex Prop. 1).
#'
#' The output has the SAME structure that step2 (\code{--estimate_cross_anc_rg=TRUE}) and
#' \code{step3_crossAncestryRgRHE.R} consume, so the existing, validated steps run unchanged:
#' \itemize{
#'   \item quantitative: OLS fit (optional rank-inverse-normal of y); theta = (1/var-free, 0).
#'   \item binary: logistic GLM fit; theta = (1, 0); the RHE uses the PCGC residual internally.
#' }
#'
#' Do NOT use this for related samples -- use \code{fitNULLGLMM()} there. For related cohorts,
#' restrict to an unrelated subset first (e.g. \code{plink2 --king-cutoff 0.0884}) via
#' \code{sampleFile}. This function never touches GWAS code paths.
#'
#' @param phenoFile tab-delimited phenotype file with a header.
#' @param phenoCol name of the phenotype column.
#' @param outputPrefix output prefix; writes \code{<prefix>.rda} and
#'   \code{<prefix>.varianceRatio.txt}.
#' @param covarColList comma-separated covariate column names. MUST include ancestry PCs /
#'   global-ancestry fractions to de-confound the cross-ancestry signal.
#' @param sampleIDColinphenoFile sample-ID column name (default "IID").
#' @param traitType "quantitative" or "binary".
#' @param invNormalize logical; rank-inverse-normal transform a quantitative phenotype.
#' @param sampleFile optional file of IIDs (one per line) to restrict to (e.g. the unrelated
#'   subset). If empty, all samples in \code{phenoFile} are used.
#' @return (invisibly) the \code{modglmm} list that is also saved to \code{<prefix>.rda}.
#' @export
fitNULL_noGRM <- function(phenoFile, phenoCol, outputPrefix,
                          covarColList = "", sampleIDColinphenoFile = "IID",
                          traitType = c("quantitative", "binary"),
                          invNormalize = FALSE, sampleFile = "") {
  traitType <- match.arg(traitType)
  stopifnot(nzchar(phenoFile), nzchar(phenoCol), nzchar(outputPrefix))

  ## GRM-free covariate projection (inlined from SAIGE's ScoreTest_NULL_Model)
  .scoreTestNullModel <- function(mu, mu2, y, X) {
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
  .invNormal <- function(y) stats::qnorm((rank(y) - 0.5) / length(y))

  ## ---- read phenotype + covariates ----
  ph <- utils::read.table(phenoFile, header = TRUE, stringsAsFactors = FALSE,
                          sep = "\t", check.names = FALSE)
  idcol <- sampleIDColinphenoFile
  stopifnot(idcol %in% colnames(ph), phenoCol %in% colnames(ph))
  if (nzchar(sampleFile)) {
    keep <- readLines(sampleFile)
    keep <- keep[nzchar(keep)]
    ph <- ph[ph[[idcol]] %in% keep, , drop = FALSE]
    cat(sprintf("[fitNULL_noGRM] restricted to %d unrelated samples from %s\n",
                nrow(ph), sampleFile))
  }
  covars <- if (nzchar(covarColList)) strsplit(covarColList, ",")[[1]] else character(0)
  covars <- trimws(covars)
  miss <- setdiff(covars, colnames(ph))
  if (length(miss)) stop("covariates not found in phenoFile: ", paste(miss, collapse = ", "))

  ph <- ph[stats::complete.cases(ph[, c(phenoCol, covars, idcol)]), , drop = FALSE]
  sampleID <- as.character(ph[[idcol]])
  N <- nrow(ph)
  y <- as.numeric(ph[[phenoCol]])
  X <- cbind(Intercept = 1, as.matrix(ph[, covars, drop = FALSE]))
  storage.mode(X) <- "double"
  cat(sprintf("[fitNULL_noGRM] N=%d trait=%s covars=[%s] invNorm=%s\n",
              N, traitType, paste(covars, collapse = ","), isTRUE(invNormalize)))

  ## ---- GRM-free null fit ----
  if (traitType == "quantitative") {
    if (isTRUE(invNormalize)) y <- .invNormal(y)
    alpha <- as.vector(solve(t(X) %*% X, t(X) %*% y))
    eta <- as.vector(X %*% alpha)
    mu <- eta
    res <- y - mu
    # theta[1]=1 (tau_g=0): with variance ratio = 1 the RHE score is exactly g~'res (no
    # tau/vr divisor), consistent with q_e = res'Sigma^-1 res = res.res under Sigma^-1 = I.
    # h2/rg are scale-invariant in res; only tau1*vr = 1 matters (else the ill-conditioned
    # multi-GRM HE collapses sigma2_anc -- verified: tau1*vr != 1 shifts compScoreSum by
    # (tau1*vr)^2).
    tau <- c(1, 0)
    mu2 <- rep(1, N)                                       # 1/tau1; Sigma^-1 diag = mu2*tau1 = 1
    obj_cc <- NULL
  } else {                                                 # binary
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

  objnoK <- .scoreTestNullModel(mu, mu2, y, X)

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
    traitType = traitType,
    isCovariateOffset = FALSE,
    LOCO = FALSE,
    obj.glm.null = NULL,
    offset = rep(0, N),
    useSparseGRMtoFitNULL = FALSE
  )
  if (!is.null(obj_cc)) modglmm$obj_cc <- obj_cc

  save(modglmm, file = paste0(outputPrefix, ".rda"))
  ## variance ratio = 1 (no GRM): "value <category> <index>" matching Get_Variance_Ratio
  writeLines("1\tnull\t1", paste0(outputPrefix, ".varianceRatio.txt"))
  cat(sprintf("[fitNULL_noGRM] wrote %s.rda (theta=%.4f,%.0f) and %s.varianceRatio.txt (=1)\n",
              outputPrefix, tau[1], tau[2], outputPrefix))
  invisible(modglmm)
}
