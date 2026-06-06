#!/usr/bin/env Rscript

## Step 2 (rg/h2 ONLY, UNRELATED samples): standalone cross-ancestry genetic correlation
## and per-ancestry heritability via randomized Haseman-Elston, WITHOUT the association
## scan. GLMM-free: residualize the phenotype by OLS (quant) or PCGC (binary), diagonal
## Sigma (= identity; see Prop.1). Reads per-ancestry Z + MAF mask + marker ids via the
## readZaChunkInCPP block reader; ALL subset/marker logic is done here in R.
##
## OPTIONS (mirror the integrated step2 --estimate_cross_anc_rg):
##   --rg_ancestry_list "2,3"  restrict to a SUBSET of ancestries (1-indexed; default all).
##                             rg = all pairs among the subset, h2 = joint over the subset.
##   --rg_markerFile FILE      markers used for rg (one variant-id or chr:pos per line; default all)
##   --h2_markerFile FILE      markers used for h2 (default all)
##   --rg_perAncestryH2        also report per-ancestry h2 on each ancestry's OWN markers (M_a)
## Output: <outFile> (per-pair rg) + <outFile>.h2 (per-ancestry h2). rg_pval = 2-sided H1 rg!=0;
## rg_pval_vs1 = radmix-style 1-sided H0 rg=1 vs H1 rg<1 ("are causal effects shared?").

suppressPackageStartupMessages({ library(SAIGE); library(optparse); library(data.table) })

op <- list(
  make_option("--tractorHybridPrefix", type="character", default=""),
  make_option("--chrom", type="character", default=""),
  make_option("--phenoFile", type="character", default=""),
  make_option("--phenoCol", type="character", default=""),
  make_option("--sampleIDColinphenoFile", type="character", default="IID"),
  make_option("--covarColList", type="character", default="",
              help="comma-separated covariates; MUST include ancestry PCs/global-anc"),
  make_option("--traitType", type="character", default="quantitative"),
  make_option("--invNormalize", type="character", default="FALSE"),
  make_option("--prevalence", type="double", default=NA_real_),
  make_option("--sampleFile", type="character", default="",
              help="unrelated IID list (one per line); e.g. plink2 --king-cutoff 0.0884 .in.id"),
  make_option("--numberOfAncestry", type="integer", default=0),
  make_option("--rg_ancestry_list", type="character", default="",
              help="subset of ancestries, 1-indexed comma list e.g. '2,3' (default: all)"),
  make_option("--rg_markerFile", type="character", default="",
              help="markers for rg (variant-id or chr:pos per line; default all)"),
  make_option("--h2_markerFile", type="character", default="",
              help="markers for h2 (default all)"),
  make_option("--rg_perAncestryH2", type="logical", default=FALSE,
              help="also report per-ancestry h2 on each ancestry's OWN markers (h2_ownMarkers)"),
  make_option("--rg_nProbes", type="integer", default=60),
  make_option("--rg_nJackknifeBlocks", type="integer", default=20),
  make_option("--maf", type="double", default=0.01),
  make_option("--chunkSize", type="integer", default=2000),
  make_option("--rg_seed", type="integer", default=1),
  make_option("--outFile", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=op))
stopifnot(nzchar(opt$tractorHybridPrefix), nzchar(opt$phenoFile), nzchar(opt$outFile))
trait <- match.arg(opt$traitType, c("quantitative","binary"))
invNorm <- toupper(opt$invNormalize) %in% c("TRUE","T","1")
Kall <- opt$numberOfAncestry
if (Kall <= 0) {
  meta <- readLines(paste0(opt$tractorHybridPrefix, ".meta"))
  Kall <- as.integer(sub(".*\\s", "", grep("^n_ancestries", meta, value=TRUE)[1]))
}
## ancestry subset (1-indexed into the full Kall)
aidx <- if (nzchar(opt$rg_ancestry_list)) {
  as.integer(trimws(strsplit(opt$rg_ancestry_list, ",")[[1]]))
} else seq_len(Kall)
stopifnot(all(aidx >= 1 & aidx <= Kall), length(aidx) >= 1)
K <- length(aidx)                                       # number of ANALYSED ancestries
if (K < 2) opt$rg_perAncestryH2 <- TRUE                 # single ancestry -> only h2 makes sense
loadset <- function(f) if (nzchar(f)) { v <- readLines(f); unique(trimws(v[nzchar(v)])) } else NULL
rgSet <- loadset(opt$rg_markerFile); h2Set <- loadset(opt$h2_markerFile)

## ---- phenotype + covariates ----
ph <- as.data.frame(fread(opt$phenoFile))
idc <- opt$sampleIDColinphenoFile
stopifnot(idc %in% names(ph), opt$phenoCol %in% names(ph))
if (nzchar(opt$sampleFile)) { keep <- readLines(opt$sampleFile); keep <- keep[nzchar(keep)]
  ph <- ph[ph[[idc]] %in% keep, , drop=FALSE] }
covars <- trimws(strsplit(opt$covarColList, ",")[[1]]); covars <- covars[nzchar(covars)]
ph <- ph[stats::complete.cases(ph[, c(opt$phenoCol, covars, idc)]), , drop=FALSE]
sampleID <- as.character(ph[[idc]]); N <- nrow(ph)
y <- as.numeric(ph[[opt$phenoCol]])
X <- cbind(1, as.matrix(ph[, covars, drop=FALSE])); storage.mode(X) <- "double"
caseProp <- if (trait=="binary") mean(y) else NA_real_
if (trait=="quantitative" && invNorm) y <- qnorm((rank(y)-0.5)/N)
beta <- qr.solve(crossprod(X), crossprod(X, y)); yt <- as.vector((y - X %*% beta))
yt <- yt / sd(yt); yty <- sum(yt*yt); XtX_inv <- solve(crossprod(X))
cat(sprintf("[rgRHEonly] N=%d ancestries=[%s] of %d trait=%s B=%d blocks=%d%s%s%s\n",
            N, paste(aidx, collapse=","), Kall, trait, opt$rg_nProbes, opt$rg_nJackknifeBlocks,
            if (!is.null(rgSet)) sprintf(" rgMk=%d", length(rgSet)) else "",
            if (!is.null(h2Set)) sprintf(" h2Mk=%d", length(h2Set)) else "",
            if (opt$rg_perAncestryH2) " +perAncH2" else ""))

## ---- reader setup ----
objGeno <- setGenoInput(tractorHybridPrefix=opt$tractorHybridPrefix, chrom=opt$chrom,
                        sampleInModel=sampleID)
genoType <- objGeno$genoType; maxGIndex <- 5000000L

## ---- accumulators ----
pairs <- if (K >= 2) t(combn(K, 2)) else matrix(integer(0), 0, 2)   # pairs over the SUBSET (1..K)
npair <- nrow(pairs); dJ <- K + npair + 1
G <- opt$rg_nJackknifeBlocks; B <- opt$rg_nProbes
set.seed(opt$rg_seed); P <- matrix(sample(c(-1,1), N*B, replace=TRUE), N, B)
## joint rg/h2 (over subset ancestries, on rg markers all-subset-ok)
uA <- lapply(1:G, function(g) lapply(seq_len(K), function(a) matrix(0,N,B)))
uX <- if (npair>0) lapply(1:G, function(g) lapply(seq_len(npair), function(p) matrix(0,N,B))) else NULL
qd <- matrix(0,G,K); qx <- matrix(0,G,max(npair,1)); MblkJ <- numeric(G); gptr <- 0L; mJoint <- 0L
## per-ancestry own-marker h2 (on h2 markers, a-ok)
doOwn <- isTRUE(opt$rg_perAncestryH2)
uH <- if (doOwn) lapply(1:G, function(g) lapply(seq_len(K), function(a) matrix(0,N,B))) else NULL
qH <- matrix(0,G,K); MblkH <- matrix(0,G,K)
emptyChunks <- 0L

accum_block <- function(uList, qvec, MblkVec, Zc, ytv, g, ai) {
  pj <- crossprod(Zc, P)                                # nc x B
  uList[[g]][[ai]] <- uList[[g]][[ai]] + Zc %*% pj
  list(u=uList, qadd=sum(ytv^2), nc=ncol(Zc))
}

for (s0 in seq(0L, maxGIndex, by=opt$chunkSize)) {
  gi <- as.character(s0:(s0 + opt$chunkSize - 1L))
  res <- readZaChunkInCPP(genoType, gi, Kall, opt$maf)
  m <- res$nRead
  if (m == 0) { emptyChunks <- emptyChunks+1L; if (emptyChunks>=2L && gptr>0L) break else next }
  emptyChunks <- 0L
  ok <- res$ok                                          # Kall x m
  inset <- function(S) if (is.null(S)) rep(TRUE,m) else (res$markerID %in% S) | (res$markerCP %in% S)
  rgKeep <- inset(rgSet); h2Keep <- inset(h2Set)
  ## covariate-adjust subset ancestries
  Zt <- lapply(seq_len(K), function(s) { Za <- res$Za[[ aidx[s] ]]
    Za - X %*% (XtX_inv %*% crossprod(X, Za)) })
  oks <- ok[aidx, , drop=FALSE]                         # K x m (subset)
  ytZ <- lapply(seq_len(K), function(s) as.vector(crossprod(Zt[[s]], yt)))
  blockOf <- ((gptr + seq_len(m) - 1L) %% G) + 1L; gptr <- gptr + m
  ## JOINT rg/h2: markers where ALL subset ancestries pass MAF AND in rg set
  jmask <- rgKeep & (if (K==1) as.logical(oks[1,]) else apply(oks==1,2,all))
  for (g in 1:G) {
    cols <- which(blockOf==g & jmask); if (length(cols)>0) {
      Zc <- lapply(seq_len(K), function(s) Zt[[s]][,cols,drop=FALSE])
      pj <- lapply(seq_len(K), function(s) crossprod(Zc[[s]], P))
      for (s in seq_len(K)) { uA[[g]][[s]] <- uA[[g]][[s]] + Zc[[s]] %*% pj[[s]]
                              qd[g,s] <- qd[g,s] + sum(ytZ[[s]][cols]^2) }
      if (npair>0) for (p in 1:npair) { a<-pairs[p,1]; b<-pairs[p,2]
        uX[[g]][[p]] <- uX[[g]][[p]] + Zc[[a]] %*% pj[[b]] + Zc[[b]] %*% pj[[a]]
        qx[g,p] <- qx[g,p] + sum(ytZ[[a]][cols]*ytZ[[b]][cols]) }
      MblkJ[g] <- MblkJ[g] + length(cols); mJoint <- mJoint + length(cols)
    }
  }
  ## PER-ANCESTRY own-marker h2: each subset ancestry a, markers a-ok AND in h2 set
  if (doOwn) for (s in seq_len(K)) { amask <- h2Keep & (oks[s,]==1)
    for (g in 1:G) { cols <- which(blockOf==g & amask); if (length(cols)>0) {
      Zc <- Zt[[s]][,cols,drop=FALSE]
      uH[[g]][[s]] <- uH[[g]][[s]] + Zc %*% crossprod(Zc, P)
      qH[g,s] <- qH[g,s] + sum(ytZ[[s]][cols]^2); MblkH[g,s] <- MblkH[g,s] + length(cols) } } }
  cat(sprintf("\r[rgRHEonly] gIndex %d  joint_mk=%d", s0+opt$chunkSize-1L, mJoint), file=stderr())
}
cat("\n", file=stderr()); stopifnot(mJoint > 0)

## ---- solvers ----
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
est_joint <- function(th) { s2<-th[1:K]; gam<-if(npair>0) th[(K+1):(K+npair)] else numeric(0)
  tot<-sum(s2)+sum(gam)+th[length(th)]; h2<-if(tot>0) s2/tot else rep(NA,K)
  rg<-if(npair>0) sapply(1:npair,function(p){a<-pairs[p,1];b<-pairs[p,2]
    if(s2[a]>0&&s2[b]>0) gam[p]/sqrt(s2[a]*s2[b]) else NA}) else numeric(0)
  list(s2=s2,h2=h2,rg=rg,gam=gam) }
solve_own <- function(blocks, s) {  # 2-comp (K_s, I) -> h2_own
  M <- sum(MblkH[blocks,s]); if (M<=0) return(NA)
  u <- Reduce(`+`, lapply(blocks, function(g) uH[[g]][[s]]))/M
  S <- matrix(c(sum(u*u)/B, sum(u*P)/B, sum(u*P)/B, sum(P*P)/B), 2, 2)
  qv <- c(sum(qH[blocks,s])/M, yty); th <- solve(S, qv)
  if ((th[1]+th[2])>0) th[1]/(th[1]+th[2]) else NA
}

fullJ <- est_joint(solve_joint(1:G))
jk_h2 <- matrix(NA,G,K); jk_rg <- matrix(NA,G,max(npair,1)); jk_own <- matrix(NA,G,K)
for (g in 1:G) { th <- solve_joint(setdiff(1:G,g))
  if (!is.null(th)) { e<-est_joint(th); jk_h2[g,]<-e$h2; if (npair>0) jk_rg[g,]<-e$rg }
  if (doOwn) for (s in seq_len(K)) jk_own[g,s] <- solve_own(setdiff(1:G,g), s) }
jkse <- function(mat) apply(mat,2,function(v){v<-v[is.finite(v)]
  if(length(v)>1) sqrt((length(v)-1)/length(v)*sum((v-mean(v))^2)) else NA})
h2_se <- jkse(jk_h2); rg_se <- if (npair>0) jkse(jk_rg) else numeric(0)
own <- if (doOwn) sapply(seq_len(K), function(s) solve_own(1:G,s)) else rep(NA,K)
own_se <- if (doOwn) jkse(jk_own) else rep(NA,K)

## liability factor (binary)
lf <- NA
if (trait=="binary") { Kp <- if (is.na(opt$prevalence)) caseProp else opt$prevalence
  if (Kp>0 && Kp<1) { z<-dnorm(qnorm(1-Kp)); P0<-if(caseProp>0&&caseProp<1) caseProp else Kp
    lf <- (Kp*(1-Kp))^2/(P0*(1-P0)*z*z) } }

## ---- rg table (ancestry labels = original aidx) ----
if (npair>0) {
  rg_z <- fullJ$rg/rg_se; rg_p <- 2*pnorm(abs(rg_z),lower.tail=FALSE)
  rg_z1 <- (fullJ$rg-1)/rg_se; rg_p1 <- pnorm(rg_z1)
  rgtab <- data.frame(anc_a=aidx[pairs[,1]], anc_b=aidx[pairs[,2]], rg=fullJ$rg, rg_se=rg_se,
                      rg_z=rg_z, rg_pval=rg_p, rg_z_vs1=rg_z1, rg_pval_vs1=rg_p1,
                      cov_cross=fullJ$gam, sigma2_anc_a=fullJ$s2[pairs[,1]],
                      sigma2_anc_b=fullJ$s2[pairs[,2]], n_indiv=N, n_markers=mJoint)
  fwrite(rgtab, opt$outFile, sep="\t")
}
## ---- h2 table ----
h2_z <- fullJ$h2/h2_se; h2_p <- pnorm(h2_z,lower.tail=FALSE)
h2tab <- data.frame(ancestry=aidx, h2_joint=fullJ$h2, h2_joint_se=h2_se,
                    h2_joint_z=h2_z, h2_joint_pval=h2_p)
if (doOwn) { ow_z <- own/own_se
  h2tab$h2_ownMarkers <- own; h2tab$h2_ownMarkers_se <- own_se
  h2tab$h2_ownMarkers_pval <- pnorm(ow_z, lower.tail=FALSE) }
if (trait=="binary") { h2tab$h2_joint_liability <- fullJ$h2*lf
  h2tab$h2_joint_liability_se <- h2_se*lf
  if (doOwn) { h2tab$h2_ownMarkers_liability <- own*lf } }
fwrite(h2tab, paste0(opt$outFile, ".h2"), sep="\t")

cat("\n--- rg (rg_pval_vs1 = radmix 'effects shared?' H0 rg=1) ---\n")
if (npair>0) for (p in 1:npair) cat(sprintf("  rg(%d,%d) = %7.3f +/- %.3f   p(rg!=0)=%.2g  p(rg=1)=%.2g\n",
  aidx[pairs[p,1]], aidx[pairs[p,2]], fullJ$rg[p], rg_se[p],
  2*pnorm(abs(fullJ$rg[p]/rg_se[p]),lower.tail=FALSE), pnorm((fullJ$rg[p]-1)/rg_se[p])))
cat("--- h2 ---\n")
for (s in seq_len(K)) cat(sprintf("  h2_%d (joint) = %7.3f +/- %.3f%s%s\n", aidx[s], fullJ$h2[s], h2_se[s],
  if (doOwn) sprintf("   own=%.3f", own[s]) else "",
  if (trait=="binary") sprintf("   liab=%.3f", fullJ$h2[s]*lf) else ""))
cat(sprintf("[rgRHEonly] wrote %s%s and %s.h2\n", if(npair>0) opt$outFile else "(no rg: <2 ancestries)",
            "", opt$outFile))
