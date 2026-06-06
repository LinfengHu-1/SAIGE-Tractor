## ============================================================================
## RHE individual-level cross-ancestry rg (map-reduce combine + solve)
## ============================================================================
## Reads per-chunk RHE partial accumulators written by step2 (when run with
## --estimate_cross_anc_rg=TRUE): file <prefix>.rg_partial.bin. Each holds, over
## that chunk's markers and a SHARED set of B +/-1 probes:
##   Ar[c] (n x B)  = sum_j gtilde_{c,j} (gtilde_{c,j} . r_b)          (A-side)
##   Kr[c] (n x B)  = sum_j (Sig^-1 gtilde_{c,j})((Sig^-1 gtilde). r_b)(K-side)
##   for components c = anc_1 .. anc_K, cross; plus q-side score sums and M.
## The combine sums these across chunks (additive: same probes), builds the HE
## trace matrix S and rhs q, and solves S theta = q -> sigma2_anc, r_admix.
## Block jackknife uses each chunk as one leave-one-out block.

#' @keywords internal
.read_rg_partial2 <- function(path) {
  # Stream from a file connection: readBin advances the position itself, so we avoid
  # loading the whole file as a raw vector and the per-field raw[(off+1):(off+k)] index
  # slicing (which allocated a large index vector + raw copy per read -- the memory peak).
  con <- file(path, "rb"); on.exit(close(con))
  rd_i32 <- function(nn = 1L) readBin(con, "integer", n = nn, size = 4L, endian = "little")
  rd_dbl <- function(nn) readBin(con, "double", n = nn, size = 8L, endian = "little")
  rd_f32 <- function(nn) readBin(con, "double", n = nn, size = 4L, endian = "little")  # float32 -> double
  magic <- rawToChar(readBin(con, "raw", n = 8L))
  if (magic != "RHEPART4") stop("bad RHE partial magic in ", path,
       " (got '", magic, "'; expected RHEPART4 -- re-run step2 with the current build)")
  n <- rd_i32(); B <- rd_i32(); K <- rd_i32(); seed <- rd_i32()
  ncomp <- rd_i32(); P <- rd_i32(); isBin <- rd_i32(); storeKr <- rd_i32(); perH2 <- rd_i32()
  M <- rd_dbl(1); caseProp <- rd_dbl(1)
  # component metadata: 4 ints per component (gA, gB, kind, pairId), 0-indexed.
  meta <- if (ncomp > 0) matrix(rd_i32(4L * ncomp), nrow = 4) else matrix(integer(0), 4, 0)
  gA <- meta[1, ]; gB <- meta[2, ]; kind <- meta[3, ]; pairId <- meta[4, ]
  pmeta <- if (P > 0) matrix(rd_i32(2L * P), nrow = 2) else matrix(integer(0), 2, 0)
  pairA <- pmeta[1, ]; pairB <- pmeta[2, ]
  compScoreSum <- rd_dbl(ncomp); compVar <- rd_dbl(ncomp); Mcomp <- rd_dbl(ncomp)
  trSigmaInv <- rd_dbl(1); resDotRes <- rd_dbl(1)
  Ddiag <- if (storeKr == 0L) rd_f32(n) else NULL
  Ar <- vector("list", ncomp)
  for (c in seq_len(ncomp)) Ar[[c]] <- matrix(rd_f32(n * B), n, B)   # column-major float32
  # Non-sparse (storeKr==0): Kr = Ddiag %*% Ar with Ddiag a per-sample vector identical
  # across files, so we do NOT materialize Kr here (it doubled the matrix memory). The
  # solve forms the needed Kr.Ar inner products from Ar + Ddiag on the fly. Sparse path
  # (storeKr!=0) stores Kr explicitly as before.
  Kr <- NULL
  if (storeKr != 0L) {
    Kr <- vector("list", ncomp)
    for (c in seq_len(ncomp)) Kr[[c]] <- matrix(rd_f32(n * B), n, B)
  }
  list(n = n, B = B, K = K, seed = seed, ncomp = ncomp, P = P,
       isBinary = (isBin != 0L), perAncestryH2 = (perH2 != 0L), caseProp = caseProp,
       storeKr = (storeKr != 0L), Ddiag = Ddiag,
       gA = gA, gB = gB, kind = kind, pairId = pairId, pairA = pairA, pairB = pairB,
       compScoreSum = compScoreSum, compVar = compVar, Mcomp = Mcomp, M = M,
       trSigmaInv = trSigmaInv, resDotRes = resDotRes, Ar = Ar, Kr = Kr)
}

#' Solve a Haseman-Elston system over a chosen set of genetic components + e.
#' gidx are 1-based component indices; returns the variance components (aligned
#' with gidx) and sigma2_e. Cross components (gA!=gB) carry a factor 2 in the
#' q / e-trace scalars (K_ab is the symmetric Z_a Z_b' + Z_b Z_a'); the genetic
#' block already encodes it via the accumulated Ar/Kr.
#' @keywords internal
#' All component-pair Kr.Ar inner products in one BLAS crossprod (G[k,l] = sum(Kr_k * Ar_l)).
#' Replaces the per-pair, per-solve interpreted sum(Kr_k * Ar_l) (which was ~76% of step3
#' CPU) with a single (n*B x ncomp)^T (n*B x ncomp) BLAS multiply. Non-sparse: Kr_k = Ar_k *
#' Ddiag (per-sample weight, recycled over the n*B column-major vec), so G = crossprod(Ar*w, Ar)
#' (symmetric). Sparse: G = crossprod(Kr, Ar) using the stored Kr.
#' @keywords internal
.gram_matrix <- function(acc) {
  ncomp <- length(acc$Ar)
  nB <- length(acc$Ar[[1]])
  Am <- matrix(0, nB, ncomp)
  for (c in seq_len(ncomp)) Am[, c] <- acc$Ar[[c]]        # column-major vec of each Ar
  if (is.null(acc$Kr)) {
    w <- rep(acc$Ddiag, length.out = nB)                  # Ddiag per entry (recycled B times)
    crossprod(Am * w, Am)
  } else {
    Km <- matrix(0, nB, ncomp)
    for (c in seq_len(ncomp)) Km[, c] <- acc$Kr[[c]]
    crossprod(Km, Am)
  }
}

.he_solve_components <- function(acc, gidx, B) {
  d <- length(gidx) + 1L
  S <- matrix(0, d, d); q <- numeric(d)
  Mc <- acc$Mcomp
  G <- acc$Gram                                           # precomputed in .solve_all
  for (i in seq_along(gidx)) {
    k <- gidx[i]
    fk <- if (acc$gA[k] != acc$gB[k]) 2 else 1     # cross vs diagonal
    S[d, i] <- S[i, d] <- fk * acc$compVar[k] / Mc[k]      # tr(Sig^-1 K_k)
    q[i] <- fk * acc$compScoreSum[k] / Mc[k]               # y' K_k y
    for (j in seq_along(gidx)) {
      l <- gidx[j]
      S[i, j] <- (G[k, l] + G[l, k]) / (2 * B * Mc[k] * Mc[l])
    }
  }
  S[d, d] <- acc$trSigmaInv; q[d] <- acc$resDotRes
  sv <- svd(S); tol <- max(dim(S)) * .Machine$double.eps * max(sv$d)
  dinv <- ifelse(sv$d > tol, 1 / sv$d, 0)
  theta <- as.numeric(sv$v %*% (dinv * (t(sv$u) %*% q)))
  list(genetic = theta[seq_along(gidx)], sigma2_e = theta[d])
}

#' Observed-scale -> liability-scale heritability factor (Lee et al. 2011, AJHG).
#'
#' h2_liab = h2_obs * K(1-K)^2 ... = h2_obs * [K(1-K)]^2 / (P(1-P) z^2), where
#' K = population prevalence, P = sample case proportion, t = qnorm(1-K) the
#' liability threshold, z = dnorm(t). With P = K (population sample, no
#' ascertainment) this reduces to K(1-K)/z^2.
#' @keywords internal
.rg_liability_factor <- function(prevalence, caseProp) {
  K <- prevalence; P <- caseProp
  if (!is.finite(K) || K <= 0 || K >= 1) return(NA_real_)
  t <- stats::qnorm(1 - K)
  z <- stats::dnorm(t)
  if (z <= 0) return(NA_real_)
  if (!is.finite(P) || P <= 0 || P >= 1) P <- K   # assume population sample
  (K * (1 - K))^2 / (P * (1 - P) * z^2)
}

#' Estimate cross-ancestry rg from RHE partial accumulators (combine step)
#'
#' @param partialFiles Character vector of <prefix>.rg_partial.bin files (one per
#'   step2 chunk/chromosome, all using the same --rg_seed). Supports K>=2
#'   ancestries (all K(K-1)/2 pairwise rg are reported).
#' @param outFile Optional path to write the per-pair rg table. Per-ancestry
#'   heritabilities are written to <outFile>.h2 when outFile is given.
#' @param prevalence Optional population prevalence K for a binary trait. When
#'   given (and the trait was binary), per-ancestry h2 is also reported on the
#'   liability scale (Lee et al. 2011). rg itself is scale-free (no prevalence
#'   needed). If omitted for a binary trait, the sample case proportion is used.
#' @return data.frame of per-pair rg (one row per ancestry pair), invisibly. The
#'   per-ancestry heritability table (with h2_se, h2_pval) is attached as attribute
#'   "h2". Standard errors and p-values are block-jackknife: they require >1 partial
#'   file, obtained either from multiple step2 chunks/chromosomes OR from a single
#'   step2 run with --rg_nJackknifeBlocks>1 (one partial file per block).
#'   rg_pval is a TWO-SIDED Wald test of rg != 0 (rg is interior; unreliable near
#'   rg = +/-1, where a profile-likelihood LRT is preferred). rg_pval_vs1 is the
#'   radmix-style ONE-SIDED test of perfect sharing H0: rg=1 vs H1: rg<1 (reject =>
#'   causal effects DIFFER across ancestries); Wald-based, so unreliable at the rg=1
#'   boundary / large rg_se -- radmix uses a profile-likelihood LRT. h2_pval is the
#'   ONE-SIDED boundary test of h2 > 0 (variance component >= 0; the one-sided
#'   Wald p equals the 1/2:1/2 chi-square mixture LRT p).
#'
#'   The partial files encode which analysis was run in step2:
#'   * default (joint all-K, "A"): rg for every pair on markers shared by all
#'     ancestries; h2 from the same joint model.
#'   * --rg_pairs="1-2,..": rg estimated COHERENTLY per listed pair on M_ab
#'     (markers polymorphic in both) -- numerator and denominator share a set.
#'   * --rg_perAncestryH2: each ancestry's h2 on its OWN markers M_a (h2_ownMarkers).
#' @export
estimateCrossAncestryRgRHE <- function(partialFiles, outFile = "", prevalence = NA_real_) {
  stopifnot(length(partialFiles) >= 1)
  parts <- lapply(partialFiles, .read_rg_partial2)
  h <- parts[[1]]
  n <- h$n; B <- h$B; K <- h$K; ncomp <- h$ncomp; P <- h$P
  for (p in parts) if (p$n != n || p$B != B || p$K != K || p$ncomp != ncomp)
    stop("RHE partials disagree on n/B/K/ncomp (same --rg_seed / --rg_pairs / samples?)")
  isBinary <- h$isBinary; caseProp <- h$caseProp; perH2 <- h$perAncestryH2
  aMode <- any(h$kind %in% c(0L, 1L))            # joint all-K components present

  .combine <- function(ps) {
    acc <- list(gA = h$gA, gB = h$gB, kind = h$kind, pairId = h$pairId,
                trSigmaInv = h$trSigmaInv, resDotRes = h$resDotRes,
                storeKr = h$storeKr, Ddiag = h$Ddiag,   # Ddiag identical across blocks
                compScoreSum = Reduce(`+`, lapply(ps, `[[`, "compScoreSum")),
                compVar = Reduce(`+`, lapply(ps, `[[`, "compVar")),
                Mcomp = Reduce(`+`, lapply(ps, `[[`, "Mcomp")),
                M = sum(vapply(ps, `[[`, numeric(1), "M")))
    acc$Ar <- lapply(seq_len(ncomp), function(c) Reduce(`+`, lapply(ps, function(p) p$Ar[[c]])))
    # Non-sparse: Kr is implied by Ar + Ddiag (see .he_solve_components), so don't sum it.
    acc$Kr <- if (h$storeKr)
      lapply(seq_len(ncomp), function(c) Reduce(`+`, lapply(ps, function(p) p$Kr[[c]]))) else NULL
    acc
  }

  # Leave-one-out accumulator by SUBTRACTING one block from the full combine, instead of
  # re-summing the other G-1 blocks for every block (turns the jackknife from O(G^2) to
  # O(G) matrix ops). All combined fields are additive over blocks; Ddiag/trSigmaInv/
  # resDotRes/metadata are block-invariant (n-level), so they carry over unchanged.
  .combine_drop <- function(acc_full, p) {
    acc <- acc_full
    acc$compScoreSum <- acc_full$compScoreSum - p$compScoreSum
    acc$compVar      <- acc_full$compVar      - p$compVar
    acc$Mcomp        <- acc_full$Mcomp        - p$Mcomp
    acc$M            <- acc_full$M            - p$M
    acc$Ar <- lapply(seq_len(ncomp), function(c) acc_full$Ar[[c]] - p$Ar[[c]])
    if (h$storeKr)
      acc$Kr <- lapply(seq_len(ncomp), function(c) acc_full$Kr[[c]] - p$Kr[[c]])
    acc
  }

  # Solve every reported quantity from a combined accumulator.
  .solve_all <- function(acc) {
    acc$Gram <- .gram_matrix(acc)                  # all Kr.Ar inner products, one BLAS call
    rg <- rep(NA_real_, P); gam <- rep(NA_real_, P)
    s2a <- rep(NA_real_, P); s2b <- rep(NA_real_, P)
    h2_joint <- rep(NA_real_, K); h2_own <- rep(NA_real_, K)
    if (aMode) {
      gd <- which(acc$kind == 0L); gc <- which(acc$kind == 1L)
      sol <- .he_solve_components(acc, c(gd, gc), B)
      s2 <- sol$genetic[seq_along(gd)]; gm <- sol$genetic[length(gd) + seq_along(gc)]
      names(s2) <- acc$gA[gd]                          # sigma2 keyed by ancestry (0-idx)
      tot <- sum(sol$genetic) + sol$sigma2_e
      h2v <- if (tot > 0) s2 / tot else rep(NA_real_, length(s2))
      for (a in seq_len(K)) { key <- as.character(a - 1L)
        if (key %in% names(s2)) h2_joint[a] <- h2v[key] }
      for (i in seq_along(gc)) { cc <- gc[i]; pid <- acc$pairId[cc] + 1L
        a <- acc$gA[cc]; b <- acc$gB[cc]
        va <- s2[as.character(a)]; vb <- s2[as.character(b)]
        gam[pid] <- gm[i]; s2a[pid] <- va; s2b[pid] <- vb
        if (va > 0 && vb > 0) rg[pid] <- gm[i] / sqrt(va * vb) }
    } else {
      for (pid in seq_len(P)) {                        # per-pair coherent on M_ab
        cc <- which(acc$pairId == (pid - 1L))
        a <- h$pairA[pid]; b <- h$pairB[pid]
        ad <- cc[acc$kind[cc] == 3L & acc$gA[cc] == a]
        bd <- cc[acc$kind[cc] == 3L & acc$gA[cc] == b]
        xc <- cc[acc$kind[cc] == 4L]
        sol <- .he_solve_components(acc, c(ad, bd, xc), B)
        va <- sol$genetic[1]; vb <- sol$genetic[2]; gm <- sol$genetic[3]
        s2a[pid] <- va; s2b[pid] <- vb; gam[pid] <- gm
        if (va > 0 && vb > 0) rg[pid] <- gm / sqrt(va * vb)
      }
    }
    if (perH2) for (cc in which(acc$kind == 2L)) {      # single-GRM h2 on M_a
      a <- acc$gA[cc] + 1L
      sol <- .he_solve_components(acc, cc, B)
      s <- sol$genetic[1]; e <- sol$sigma2_e
      if ((s + e) > 0) h2_own[a] <- s / (s + e)
    }
    list(rg = rg, gamma = gam, s2a = s2a, s2b = s2b, h2_joint = h2_joint, h2_own = h2_own)
  }

  acc_full <- .combine(parts)
  full <- .solve_all(acc_full)
  h2_primary <- if (perH2) full$h2_own else full$h2_joint

  .jk_se <- function(mat) apply(mat, 2, function(col) {
    v <- col[is.finite(col)]
    if (length(v) > 1) { g <- length(v); sqrt((g - 1) / g * sum((v - mean(v))^2)) } else NA_real_ })
  rg_se <- rep(NA_real_, max(P, 1L)); h2_se <- rep(NA_real_, K)
  h2joint_se <- rep(NA_real_, K)               # SE of the joint-fit per-ancestry h2
  if (length(parts) > 1) {
    rg_jk <- matrix(NA_real_, length(parts), max(P, 1L))
    h2_jk <- matrix(NA_real_, length(parts), K)
    h2joint_jk <- matrix(NA_real_, length(parts), K)
    for (idx in seq_along(parts)) {
      s <- .solve_all(.combine_drop(acc_full, parts[[idx]]))
      if (P > 0) rg_jk[idx, ] <- s$rg
      h2_jk[idx, ] <- if (perH2) s$h2_own else s$h2_joint
      if (aMode) h2joint_jk[idx, ] <- s$h2_joint
    }
    if (P > 0) rg_se <- .jk_se(rg_jk)
    h2_se <- .jk_se(h2_jk)
    if (aMode) h2joint_se <- .jk_se(h2joint_jk)
  }

  nmark <- sum(vapply(parts, `[[`, numeric(1), "M"))
  .z <- function(est, se) ifelse(is.finite(se) & se > 0, est / se, NA_real_)
  rg_z <- .z(full$rg, rg_se)
  rg_pval <- 2 * stats::pnorm(abs(rg_z), lower.tail = FALSE)
  # radmix-style test of perfect sharing: H0 rg=1 vs H1 rg<1 (one-sided boundary Wald).
  # Reject => causal effects DIFFER across ancestries (rg<1). Wald-based, so unreliable at
  # the rg=1 boundary / large rg_se; radmix uses a profile-likelihood LRT instead.
  rg_z_vs1 <- .z(full$rg - 1, rg_se)
  rg_pval_vs1 <- stats::pnorm(rg_z_vs1, lower.tail = TRUE)
  h2_z <- .z(h2_primary, h2_se)
  h2_pval <- stats::pnorm(h2_z, lower.tail = FALSE)
  # joint-fit per-ancestry h2 gets its own SE/p-value (one-sided H1: h2 > 0)
  h2joint_z <- .z(full$h2_joint, h2joint_se)
  h2joint_pval <- stats::pnorm(h2joint_z, lower.tail = FALSE)

  out <- data.frame(
    anc_a = h$pairA + 1L, anc_b = h$pairB + 1L,
    rg = full$rg, rg_se = rg_se, rg_z = rg_z, rg_pval = rg_pval,
    rg_z_vs1 = rg_z_vs1, rg_pval_vs1 = rg_pval_vs1,
    cov_cross = full$gamma, sigma2_anc_a = full$s2a, sigma2_anc_b = full$s2b,
    n_indiv = n, n_markers = nmark, n_probes = B, stringsAsFactors = FALSE)

  h2tab <- data.frame(ancestry = seq_len(K), stringsAsFactors = FALSE)
  if (aMode) {
    h2tab$h2_joint <- full$h2_joint
    h2tab$h2_joint_se <- h2joint_se
    h2tab$h2_joint_pval <- h2joint_pval
  }
  if (perH2) h2tab$h2_ownMarkers <- full$h2_own
  h2tab$h2 <- h2_primary; h2tab$h2_se <- h2_se; h2tab$h2_pval <- h2_pval
  if (isBinary) {
    K_pop <- if (is.finite(prevalence)) prevalence else caseProp
    fac <- .rg_liability_factor(K_pop, caseProp)
    h2tab$h2_liability <- h2_primary * fac
    h2tab$h2_liability_se <- h2_se * fac
    attr(out, "prevalence") <- K_pop; attr(out, "caseProp") <- caseProp
    if (!is.finite(prevalence))
      message("[rg] binary trait, no --prevalence given: liability h2 uses the sample ",
              "case proportion (", signif(caseProp, 3), ") as the population prevalence.")
  }
  attr(out, "h2") <- h2tab; attr(out, "isBinary") <- isBinary

  mode <- if (aMode) "joint all-K" else "per-pair coherent (M_ab)"
  cat(sprintf("Cross-ancestry rg [%s]: K=%d, %d pair(s), n=%d, M=%g, B=%d%s\n",
              mode, K, P, n, nmark, B, if (isBinary) " [binary]" else ""))
  cat("-- per-pair genetic correlation (rg_pval: 2-sided H1 rg!=0; ",
      "rg_pval_vs1: 1-sided H0 rg=1 vs H1 rg<1, radmix 'effects shared?' test) --\n", sep = "")
  print(out[, c("anc_a", "anc_b", "rg", "rg_se", "rg_pval", "rg_pval_vs1",
                "cov_cross", "sigma2_anc_a", "sigma2_anc_b")], row.names = FALSE)
  cat("-- per-ancestry heritability (h2_pval: one-sided H1 h2 > 0) --\n")
  print(h2tab, row.names = FALSE)
  if (length(parts) <= 1)
    cat("NOTE: SEs/p-values are NA (need >1 jackknife block). Re-run step2 with",
        "--rg_nJackknifeBlocks>1 (e.g. 20) or split across chunks/chromosomes.\n")

  if (nzchar(outFile)) {
    utils::write.table(out, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Written rg table to", outFile, "\n")
    h2file <- paste0(outFile, ".h2")
    utils::write.table(h2tab, h2file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Written h2 table to", h2file, "\n")
  }
  invisible(out)
}

