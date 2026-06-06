#!/usr/bin/env Rscript

## Step 0 (alternative): extract a maximal UNRELATED sample set from a sparse GRM.
##
## Biobank GWAS pipelines already build a sparse GRM (e.g. SAIGE createSparseGRM);
## this derives the unrelated keep-list from it directly, so you don't need to re-run
## plink2/KING. Output is a one-IID-per-line file to feed step1/step2 --sampleFile.
##
## Algorithm = the plink2 --king-cutoff greedy rule: build the graph of sample pairs
## with kinship > cutoff, then repeatedly remove the highest-remaining-degree sample
## until no pair exceeds the cutoff. The kept samples are mutually unrelated.
##
## INPUT: a SAIGE-style sparse GRM in MatrixMarket coordinate format (.mtx) + its
## sample-ID file (one IID per line, in matrix order).
##
## *** IMPORTANT for ADMIXED cohorts ***
## The kinship in the GRM must be ANCESTRY-ROBUST (KING-robust style). A standard
## GCTA-style GRM over-estimates relatedness in admixed samples because shared
## continental ancestry inflates apparent kinship -> it will over-remove individuals.
## If your sparse GRM is NOT ancestry-robust, prefer `plink2 --king-cutoff 0.0884`
## (KING-robust) instead. Also: the GRM must have been built at a kinship cutoff
## <= --relatednessCutoff, or pairs just below the build cutoff won't be in the matrix.

suppressPackageStartupMessages({ library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sparseGRMFile", type = "character", default = "",
    help = "Sparse GRM in MatrixMarket coordinate format (.mtx)."),
  make_option("--sparseGRMSampleIDFile", type = "character", default = "",
    help = "Sample-ID file for the GRM (one IID per line, in matrix order)."),
  make_option("--relatednessCutoff", type = "double", default = 0.0884,
    help = "Kinship cutoff; pairs above this are 'related' (default 0.0884 = <3rd-degree, as radmix)."),
  make_option("--outFile", type = "character", default = "",
    help = "Output: unrelated IIDs to KEEP, one per line.")
)))
stopifnot(nzchar(opt$sparseGRMFile), nzchar(opt$sparseGRMSampleIDFile), nzchar(opt$outFile))

ids <- readLines(opt$sparseGRMSampleIDFile)
ids <- ids[nzchar(ids)]
n <- length(ids)

## --- read MatrixMarket .mtx: skip % comments + the dims/nnz header line ---
con <- file(opt$sparseGRMFile, "r")
ii <- integer(0); jj <- integer(0); vv <- numeric(0)
header_seen <- FALSE
chunk <- 1e6
repeat {
  lines <- readLines(con, n = chunk)
  if (length(lines) == 0) break
  lines <- lines[!startsWith(lines, "%")]
  if (!header_seen && length(lines) > 0) { lines <- lines[-1]; header_seen <- TRUE }  # dims line
  if (length(lines) == 0) next
  parts <- strsplit(trimws(lines), "[[:space:]]+")
  m <- do.call(rbind, parts)
  a <- as.integer(m[, 1]); b <- as.integer(m[, 2]); v <- as.numeric(m[, 3])
  off <- a != b & v > opt$relatednessCutoff           # off-diagonal related pairs
  ii <- c(ii, a[off]); jj <- c(jj, b[off]); vv <- c(vv, v[off])
}
close(con)

nrel_pairs <- length(ii)
cat(sprintf("[step0] N=%d samples; %d related pairs (kinship > %.4f)\n",
            n, nrel_pairs, opt$relatednessCutoff))

## --- greedy maximal unrelated set (plink2 --king-cutoff rule) ---
removed <- logical(n)
if (nrel_pairs > 0) {
  ## adjacency as integer lists; degree = number of remaining related partners
  adj <- vector("list", n)
  for (k in seq_len(nrel_pairs)) {
    adj[[ii[k]]] <- c(adj[[ii[k]]], jj[k])
    adj[[jj[k]]] <- c(adj[[jj[k]]], ii[k])
  }
  deg <- lengths(adj)
  repeat {
    mx <- which.max(deg)
    if (deg[mx] == 0) break
    removed[mx] <- TRUE
    nb <- adj[[mx]]
    for (w in nb) if (!removed[w]) deg[w] <- deg[w] - 1L
    deg[mx] <- 0L
  }
}
keep <- ids[!removed]
writeLines(keep, opt$outFile)
cat(sprintf("[step0] kept %d unrelated / removed %d -> %s\n",
            length(keep), sum(removed), opt$outFile))
if (nrel_pairs > 0)
  cat("[step0] NOTE: for ADMIXED cohorts verify the GRM kinship is ancestry-robust ",
      "(KING-robust); a raw GCTA-style GRM over-removes. See header.\n", sep = "")
