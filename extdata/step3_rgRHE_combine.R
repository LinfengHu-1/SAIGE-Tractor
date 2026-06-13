#!/usr/bin/env Rscript
## ============================================================================
## Step 3 (rg/h2 RHE) -- GENOME-WIDE COMBINE step. Reads the per-chromosome RHE
## partials written by step2_rgRHE_chrom.R, sums them, and solves the
## Haseman-Elston system with a LEAVE-ONE-CHROMOSOME-OUT (LOCO) jackknife: each
## chromosome is one block, so the SE is LD-aware at chromosome granularity.
## Output is identical in format to step2_rgRHEonly.R: <outFile> (per-pair rg) +
## <outFile>.h2 (per-ancestry h2). All partials MUST share rg_seed / phenotype.
## ============================================================================
suppressPackageStartupMessages({ library(optparse); library(data.table) })

opt <- parse_args(OptionParser(option_list=list(
  make_option("--partialFiles", type="character", default="", help="comma list of .rds partials"),
  make_option("--partialGlob", type="character", default="", help="glob for the cell's .rds partials"),
  make_option("--prevalence", type="double", default=NA_real_),
  make_option("--outFile", type="character", default=""))))
stopifnot(nzchar(opt$outFile))
files <- character(0)
if (nzchar(opt$partialFiles)) files <- trimws(strsplit(opt$partialFiles, ",")[[1]])
if (nzchar(opt$partialGlob))  files <- c(files, Sys.glob(opt$partialGlob))
files <- unique(files[nzchar(files)])
if (length(files) < 2) stop(sprintf("need >=2 chromosome partials for a LOCO jackknife; got %d", length(files)))

## ---- read + validate partials ----
parts <- lapply(files, readRDS)
p1 <- parts[[1]]
N <- p1$N; B <- p1$B; K <- p1$K; aidx <- p1$aidx; pairs <- p1$pairs; npair <- p1$npair
yty <- p1$yty; trait <- p1$trait; caseProp <- p1$caseProp; doOwn <- p1$perAncestryH2
for (p in parts[-1]) {
  if (p$N!=N || p$B!=B || p$K!=K || p$rg_seed!=p1$rg_seed || !identical(p$aidx,aidx))
    stop("partials are inconsistent (N/B/K/rg_seed/ancestry differ) -- re-run the maps with matching options")
}
G <- length(parts)                                              # blocks = chromosomes (LOCO)
## regenerate the SHARED probes (same seed/N/B as every map)
set.seed(p1$rg_seed); P <- matrix(sample(c(-1,1), N*B, replace=TRUE), N, B)

## ---- assemble per-block (per-chromosome) accumulators ----
uA <- lapply(1:G, function(g) parts[[g]]$uA)                    # uA[[g]][[a]]  (N x B)
uX <- if (npair>0) lapply(1:G, function(g) parts[[g]]$uX) else NULL
qd <- t(sapply(parts, function(p) p$qd)); if (K==1) qd <- matrix(qd, ncol=1)   # G x K
qx <- t(sapply(parts, function(p) p$qx)); if (npair<=1) qx <- matrix(qx, ncol=max(npair,1))  # G x npair
MblkJ <- sapply(parts, function(p) p$MblkJ)                     # G
mJoint <- sum(sapply(parts, function(p) p$mJoint))
uH <- if (doOwn) lapply(1:G, function(g) parts[[g]]$uH) else NULL
qH <- t(sapply(parts, function(p) p$qH)); if (K==1) qH <- matrix(qH, ncol=1)
MblkH <- t(sapply(parts, function(p) p$MblkH)); if (K==1) MblkH <- matrix(MblkH, ncol=1)
cat(sprintf("[rgRHE-combine] %d chrom blocks, N=%d K=%d joint_mk=%d (LOCO jackknife)\n", G, N, K, mJoint))
stopifnot(mJoint > 0)

## ---- solvers (identical to step2_rgRHEonly.R, operating over the G chromosome blocks) ----
solve_joint <- function(blocks) {
  M <- sum(MblkJ[blocks]); if (M<=0) return(NULL)
  U <- c(lapply(seq_len(K), function(s) Reduce(`+`, lapply(blocks, function(g) uA[[g]][[s]]))/M),
         if (npair>0) lapply(seq_len(npair), function(p) Reduce(`+`, lapply(blocks, function(g) uX[[g]][[p]]))/M),
         list(P))
  d <- length(U); S <- matrix(0,d,d)
  for (k in 1:d) for (l in k:d) S[k,l] <- S[l,k] <- sum(U[[k]]*U[[l]])/B
  q <- c(colSums(qd[blocks,,drop=FALSE])/M, if (npair>0) 2*colSums(qx[blocks,,drop=FALSE])/M, yty)
  as.vector(solve(S,q))
}
## probe-subset solve: identical to solve_joint but uses ONLY probe columns `pcols` of every U (incl P),
## for the probe-convergence diagnostic (recompute rg from fewer probes; no genotype re-read).
solve_joint_probes <- function(blocks, pcols) {
  M <- sum(MblkJ[blocks]); if (M<=0) return(NULL)
  U <- c(lapply(seq_len(K), function(s) Reduce(`+`, lapply(blocks, function(g) uA[[g]][[s]]))/M),
         if (npair>0) lapply(seq_len(npair), function(p) Reduce(`+`, lapply(blocks, function(g) uX[[g]][[p]]))/M),
         list(P))
  bs <- length(pcols); d <- length(U); S <- matrix(0,d,d)
  for (k in 1:d) for (l in k:d) S[k,l] <- S[l,k] <- sum(U[[k]][,pcols]*U[[l]][,pcols])/bs
  q <- c(colSums(qd[blocks,,drop=FALSE])/M, if (npair>0) 2*colSums(qx[blocks,,drop=FALSE])/M, yty)
  as.vector(solve(S,q))
}
est_joint <- function(th) { s2<-th[1:K]; gam<-if(npair>0) th[(K+1):(K+npair)] else numeric(0)
  tot<-sum(s2)+sum(gam)+th[length(th)]; h2<-if(tot>0) s2/tot else rep(NA,K)
  rg<-if(npair>0) sapply(1:npair,function(p){a<-pairs[p,1];b<-pairs[p,2]
    if(s2[a]>0&&s2[b]>0) gam[p]/sqrt(s2[a]*s2[b]) else NA}) else numeric(0)
  list(s2=s2,h2=h2,rg=rg,gam=gam) }
solve_own <- function(blocks, s) {
  M <- sum(MblkH[blocks,s]); if (M<=0) return(NA)
  u <- Reduce(`+`, lapply(blocks, function(g) uH[[g]][[s]]))/M
  S <- matrix(c(sum(u*u)/B, sum(u*P)/B, sum(u*P)/B, sum(P*P)/B), 2, 2)
  qv <- c(sum(qH[blocks,s])/M, yty); th <- solve(S, qv)
  if ((th[1]+th[2])>0) th[1]/(th[1]+th[2]) else NA
}

fullJ <- est_joint(solve_joint(1:G))
jk_h2 <- matrix(NA,G,K); jk_rg <- matrix(NA,G,max(npair,1)); jk_own <- matrix(NA,G,K)
jk_s2 <- matrix(NA,G,K); jk_gam <- matrix(NA,G,max(npair,1))   # components for the delta-method rg SE
for (g in 1:G) {
  if (MblkJ[g] > 0) { th <- solve_joint(setdiff(1:G,g))
    if (!is.null(th)) { e<-est_joint(th); jk_h2[g,]<-e$h2; jk_s2[g,]<-e$s2
      if (npair>0) { jk_rg[g,]<-e$rg; jk_gam[g,]<-e$gam } } }
  if (doOwn) for (s in seq_len(K)) if (MblkH[g,s] > 0) jk_own[g,s] <- solve_own(setdiff(1:G,g), s) }
jkse <- function(mat) apply(mat,2,function(v){v<-v[is.finite(v)]
  if(length(v)>1) sqrt((length(v)-1)/length(v)*sum((v-mean(v))^2)) else NA})
h2_se <- jkse(jk_h2); rg_se <- if (npair>0) jkse(jk_rg) else numeric(0)
## delta-method jackknife SE for rg: jackknife the smooth/additive COMPONENTS (sigma2_a, sigma2_b,
## gamma) -- which never blow up -- and propagate to rg = gamma/sqrt(s_a s_b) analytically. Robust to
## the ratio exploding when a leave-one-out denominator is near zero (that heavy tail inflates the DIRECT
## rg jackknife at low h2). rg_se stays the (validated, calibrated at moderate h2) direct estimate;
## rg_se_delta is the robust alternative -- compare/choose downstream.
rg_se_delta <- if (npair>0) sapply(1:npair, function(p){
  a<-pairs[p,1]; b<-pairs[p,2]; Mc<-cbind(jk_s2[,a], jk_s2[,b], jk_gam[,p])
  ok<-rowSums(!is.finite(Mc))==0; Mc<-Mc[ok,,drop=FALSE]; nb<-nrow(Mc)
  sa<-fullJ$s2[a]; sb<-fullJ$s2[b]; gm<-fullJ$gam[p]
  if (nb<2 || sa<=0 || sb<=0) return(NA)
  cm<-colMeans(Mc); Cov<-(nb-1)/nb*crossprod(sweep(Mc,2,cm))     # 3x3 jackknife covariance of components
  rgv<-gm/sqrt(sa*sb); grad<-c(-rgv/(2*sa), -rgv/(2*sb), 1/sqrt(sa*sb))
  v<-as.numeric(t(grad)%*%Cov%*%grad); if (is.finite(v) && v>0) sqrt(v) else NA }) else numeric(0)
## probe-convergence diagnostic (mirrors estimateCrossAncestryRgRHE): jackknife over PROBES (leave out
## probe-blocks, recompute rg from the stored per-probe accumulators -- no genotype re-read) to estimate the
## Monte-Carlo SE the finite probe count B adds to rg. rg_converged = that MC SE is <=20% of the LOCO rg_se
## (probe noise negligible vs sampling); NA if there is no rg_se (single chrom) or B<2.
rg_mcse <- rep(NA_real_, max(npair,1)); rg_converged <- rep(NA, max(npair,1))
if (npair>0 && B>=2) {
  nbp <- min(B,10L); grp <- ((seq_len(B)-1L) %% nbp)+1L
  rg_pj <- matrix(NA_real_, nbp, npair)
  for (gg in seq_len(nbp)) { pc <- which(grp!=gg)
    th <- solve_joint_probes(1:G, pc); if (!is.null(th)) rg_pj[gg,] <- est_joint(th)$rg }
  rg_mcse <- jkse(rg_pj)
  rg_converged <- ifelse(is.finite(rg_mcse) & is.finite(rg_se) & rg_se>0, rg_mcse <= 0.2*rg_se, NA)
}
own <- if (doOwn) sapply(seq_len(K), function(s) solve_own(1:G,s)) else rep(NA,K)
own_se <- if (doOwn) jkse(jk_own) else rep(NA,K)

lf <- NA
if (trait=="binary") { Kp <- if (is.na(opt$prevalence)) caseProp else opt$prevalence
  if (Kp>0 && Kp<1) { z<-dnorm(qnorm(1-Kp)); P0<-if(caseProp>0&&caseProp<1) caseProp else Kp
    lf <- (Kp*(1-Kp))^2/(P0*(1-P0)*z*z) } }

if (npair>0) {
  rg_z <- fullJ$rg/rg_se; rg_p <- 2*pnorm(abs(rg_z),lower.tail=FALSE)
  rg_z1 <- (fullJ$rg-1)/rg_se; rg_p1 <- pnorm(rg_z1)
  # rg_constrained: raw rg clipped to the valid [-1,1] range -- the value to REPORT. Raw `rg` stays
  # PRIMARY (rg_se / rg_pval / rg_pval_vs1 on the raw scale). [-1,1] not [0,1] (never folds the rg=0
  # null); big boundary-RMSE win, small known boundary bias, coverage/type-I preserved. See §2.5.
  rgtab <- data.frame(anc_a=aidx[pairs[,1]], anc_b=aidx[pairs[,2]], rg=fullJ$rg,
                      rg_constrained=pmin(pmax(fullJ$rg,-1),1), rg_se=rg_se,
                      rg_se_delta=rg_se_delta, rg_z=rg_z, rg_pval=rg_p, rg_z_vs1=rg_z1, rg_pval_vs1=rg_p1,
                      rg_mcse_probe=rg_mcse, rg_converged=rg_converged,
                      cov_cross=fullJ$gam, sigma2_anc_a=fullJ$s2[pairs[,1]],
                      sigma2_anc_b=fullJ$s2[pairs[,2]], n_indiv=N, n_markers=mJoint, n_chrom=G)
  fwrite(rgtab, opt$outFile, sep="\t")
}
h2_z <- fullJ$h2/h2_se; h2_p <- pnorm(h2_z,lower.tail=FALSE)
# h2_flag (mirrors estimateCrossAncestryRgRHE): precision of the per-ancestry joint h2. "out_of_range" =
# h2<0 or >1; "imprecise" = wide CI (not significantly >0, or CI crosses 1, or SE is NA); else "ok".
# Low-exposure (minority) ancestries are high-variance -> watch for non-"ok"; rg (a ratio) is more robust.
.h2flag <- function(h, se) { if (!is.finite(h)) return("out_of_range")
  if (h<0 || h>1) return("out_of_range"); if (!is.finite(se)) return("imprecise")
  if (1.96*se >= h || h + 1.96*se > 1) return("imprecise"); "ok" }
h2tab <- data.frame(ancestry=aidx, h2_joint=fullJ$h2, h2_joint_se=h2_se, h2_joint_z=h2_z, h2_joint_pval=h2_p,
                    h2_flag=vapply(seq_len(K), function(s) .h2flag(fullJ$h2[s], h2_se[s]), character(1)))
if (doOwn) { ow_z <- own/own_se
  h2tab$h2_ownMarkers <- own; h2tab$h2_ownMarkers_se <- own_se; h2tab$h2_ownMarkers_pval <- pnorm(ow_z, lower.tail=FALSE) }
if (trait=="binary") { h2tab$h2_joint_liability <- fullJ$h2*lf; h2tab$h2_joint_liability_se <- h2_se*lf
  if (doOwn) h2tab$h2_ownMarkers_liability <- own*lf }
fwrite(h2tab, paste0(opt$outFile, ".h2"), sep="\t")

cat("\n--- rg (LOCO jackknife; rg_pval_vs1 = 'effects shared?' H0 rg=1) ---\n")
if (npair>0) for (p in 1:npair) cat(sprintf("  rg(%d,%d) = %7.3f +/- %.3f   p(rg!=0)=%.2g  p(rg=1)=%.2g\n",
  aidx[pairs[p,1]], aidx[pairs[p,2]], fullJ$rg[p], rg_se[p],
  2*pnorm(abs(fullJ$rg[p]/rg_se[p]),lower.tail=FALSE), pnorm((fullJ$rg[p]-1)/rg_se[p])))
cat("--- h2 ---\n")
for (s in seq_len(K)) cat(sprintf("  h2_%d (joint) = %7.3f +/- %.3f%s%s\n", aidx[s], fullJ$h2[s], h2_se[s],
  if (doOwn) sprintf("   own=%.3f", own[s]) else "", if (trait=="binary") sprintf("   liab=%.3f", fullJ$h2[s]*lf) else ""))
cat(sprintf("[rgRHE-combine] wrote %s and %s.h2 (%d chrom LOCO blocks)\n", opt$outFile, opt$outFile, G))
