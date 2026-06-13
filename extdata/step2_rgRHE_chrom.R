#!/usr/bin/env Rscript
## ============================================================================
## Step 2 (rg/h2 RHE) -- GENOME-WIDE MAP step: accumulate ONE chromosome's RHE
## partial into a single block and write it. The COMBINE step
## (step3_rgRHE_combine.R) sums these per-chromosome partials and solves the
## Haseman-Elston system with a LEAVE-ONE-CHROMOSOME-OUT jackknife (each
## chromosome = one block). This map-reduce parallelises the genome-wide read.
##
## All chromosomes of a phenotype MUST be run with the SAME --rg_seed (so the B
## random probes are identical and the partials are additive) and the same
## --phenoFile / --covarColList / --maf / --rg_nProbes.
## Output: <outFile> (an .rds holding the single-block accumulators + metadata).
## ============================================================================
suppressPackageStartupMessages({ library(SAIGE); library(optparse) })

opt <- parse_args(OptionParser(option_list=list(
  make_option("--tractorHybridPrefix", type="character", default=""),
  make_option("--chrom", type="character", default=""),
  make_option("--chromIndex", type="integer", default=0L, help="1-based chromosome ordinal (block id / metadata)"),
  make_option("--numberOfAncestry", type="integer", default=0L),
  make_option("--rg_ancestry_list", type="character", default=""),
  make_option("--rg_markerFile", type="character", default=""),
  make_option("--h2_markerFile", type="character", default=""),
  make_option("--rg_perAncestryH2", type="logical", default=FALSE),
  make_option("--rg_nProbes", type="integer", default=30),
  make_option("--maf", type="double", default=0.01),
  make_option("--chunkSize", type="integer", default=2000),
  make_option("--rg_seed", type="integer", default=1),
  make_option("--phenoFile", type="character", default=""),
  make_option("--phenoCol", type="character", default=""),
  make_option("--sampleIDColinphenoFile", type="character", default="IID"),
  make_option("--covarColList", type="character", default=""),
  make_option("--traitType", type="character", default="quantitative"),
  make_option("--invNormalize", type="character", default="FALSE"),
  make_option("--sampleFile", type="character", default=""),
  make_option("--prevalence", type="double", default=NA_real_),
  make_option("--outFile", type="character", default=""))))
suppressPackageStartupMessages(library(data.table))
stopifnot(nzchar(opt$tractorHybridPrefix), nzchar(opt$phenoFile), nzchar(opt$outFile))
trait <- match.arg(opt$traitType, c("quantitative","binary"))
invNorm <- toupper(opt$invNormalize) %in% c("TRUE","T","1")
Kall <- opt$numberOfAncestry
if (Kall <= 0) {
  meta <- readLines(paste0(opt$tractorHybridPrefix, ".meta"))
  Kall <- as.integer(sub(".*\\s", "", grep("^n_ancestries", meta, value=TRUE)[1]))
}
aidx <- if (nzchar(opt$rg_ancestry_list)) as.integer(trimws(strsplit(opt$rg_ancestry_list, ",")[[1]])) else seq_len(Kall)
stopifnot(all(aidx >= 1 & aidx <= Kall), length(aidx) >= 1)
K <- length(aidx)
if (K < 2) opt$rg_perAncestryH2 <- TRUE
loadset <- function(f) if (nzchar(f)) { v <- readLines(f); unique(trimws(v[nzchar(v)])) } else NULL
rgSet <- loadset(opt$rg_markerFile); h2Set <- loadset(opt$h2_markerFile)

## ---- phenotype + covariates (identical to step2_rgRHEonly.R) ----
ph <- as.data.frame(fread(opt$phenoFile)); idc <- opt$sampleIDColinphenoFile
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

## ---- reader + single-block accumulators ----
objGeno <- setGenoInput(tractorHybridPrefix=opt$tractorHybridPrefix, chrom=opt$chrom, sampleInModel=sampleID)
genoType <- objGeno$genoType; maxGIndex <- 5000000L
pairs <- if (K >= 2) t(combn(K, 2)) else matrix(integer(0), 0, 2)
npair <- nrow(pairs); B <- opt$rg_nProbes
set.seed(opt$rg_seed); P <- matrix(sample(c(-1,1), N*B, replace=TRUE), N, B)
uA <- lapply(seq_len(K), function(a) matrix(0,N,B))               # ONE block (this chromosome)
uX <- if (npair>0) lapply(seq_len(npair), function(p) matrix(0,N,B)) else NULL
qd <- numeric(K); qx <- numeric(max(npair,1)); MblkJ <- 0; mJoint <- 0L
doOwn <- isTRUE(opt$rg_perAncestryH2)
uH <- if (doOwn) lapply(seq_len(K), function(a) matrix(0,N,B)) else NULL
qH <- numeric(K); MblkH <- numeric(K); emptyChunks <- 0L; seen <- 0L

for (s0 in seq(0L, maxGIndex, by=opt$chunkSize)) {
  gi <- as.character(s0:(s0 + opt$chunkSize - 1L))
  res <- readZaChunkInCPP(genoType, gi, Kall, opt$maf)
  m <- res$nRead
  if (m == 0) { emptyChunks <- emptyChunks+1L; if (emptyChunks>=2L && seen>0L) break else next }
  emptyChunks <- 0L; seen <- seen + m
  ok <- res$ok
  inset <- function(S) if (is.null(S)) rep(TRUE,m) else (res$markerID %in% S) | (res$markerCP %in% S)
  rgKeep <- inset(rgSet); h2Keep <- inset(h2Set)
  Zt <- lapply(seq_len(K), function(s) { Za <- res$Za[[ aidx[s] ]]; Za - X %*% (XtX_inv %*% crossprod(X, Za)) })
  oks <- ok[aidx, , drop=FALSE]
  ytZ <- lapply(seq_len(K), function(s) as.vector(crossprod(Zt[[s]], yt)))
  ## JOINT: markers MAF-ok in ALL subset ancestries AND in rg set -> the single (chromosome) block
  jmask <- rgKeep & (if (K==1) as.logical(oks[1,]) else apply(oks==1,2,all))
  cols <- which(jmask)
  if (length(cols)>0) {
    Zc <- lapply(seq_len(K), function(s) Zt[[s]][,cols,drop=FALSE])
    pj <- lapply(seq_len(K), function(s) crossprod(Zc[[s]], P))
    for (s in seq_len(K)) { uA[[s]] <- uA[[s]] + Zc[[s]] %*% pj[[s]]
                            qd[s] <- qd[s] + sum(ytZ[[s]][cols]^2) }
    if (npair>0) for (p in 1:npair) { a<-pairs[p,1]; b<-pairs[p,2]
      uX[[p]] <- uX[[p]] + Zc[[a]] %*% pj[[b]] + Zc[[b]] %*% pj[[a]]
      qx[p] <- qx[p] + sum(ytZ[[a]][cols]*ytZ[[b]][cols]) }
    MblkJ <- MblkJ + length(cols); mJoint <- mJoint + length(cols)
  }
  if (doOwn) for (s in seq_len(K)) { amask <- h2Keep & (oks[s,]==1); c2 <- which(amask)
    if (length(c2)>0) { Zc <- Zt[[s]][,c2,drop=FALSE]
      uH[[s]] <- uH[[s]] + Zc %*% crossprod(Zc, P)
      qH[s] <- qH[s] + sum(ytZ[[s]][c2]^2); MblkH[s] <- MblkH[s] + length(c2) } }
  cat(sprintf("\r[rgRHE-map] chr%d gIndex %d  joint_mk=%d", opt$chromIndex, s0+opt$chunkSize-1L, mJoint), file=stderr())
}
cat("\n", file=stderr())
if (mJoint == 0) cat(sprintf("[rgRHE-map] WARNING: chr%d contributed 0 joint markers\n", opt$chromIndex))

saveRDS(list(uA=uA, uX=uX, qd=qd, qx=qx, MblkJ=MblkJ, mJoint=mJoint,
             uH=uH, qH=qH, MblkH=MblkH, N=N, B=B, K=K, Kall=Kall, aidx=aidx, pairs=pairs,
             npair=npair, yty=yty, trait=trait, caseProp=caseProp, rg_seed=opt$rg_seed,
             chromIndex=opt$chromIndex, perAncestryH2=doOwn),
        opt$outFile)
cat(sprintf("[rgRHE-map] chr%d N=%d K=%d joint_mk=%d -> %s\n", opt$chromIndex, N, K, mJoint, opt$outFile))
