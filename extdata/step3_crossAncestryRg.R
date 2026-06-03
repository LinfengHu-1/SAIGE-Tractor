#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

## Step 3 (admixed): estimate cross-ancestry genetic correlation (rho) from the
## per-variant joint score statistics emitted by step2 when run with
## --is_admixed=TRUE --is_output_cross_anc_cov=TRUE.
##
## This is a summary-statistic estimator (multivariate HESS method of moments):
## O(number of variants), independent of sample size. See estimateCrossAncestryRg.

options(stringsAsFactors = FALSE)
library(SAIGE)
library(optparse)
library(data.table)
library(methods)

option_list <- list(
  make_option("--step2File", type = "character", default = "",
    help = "Path to the step2 association output produced with --is_admixed=TRUE --is_output_cross_anc_cov=TRUE. Must contain Tstat_anc1, var_anc1, Tstat_anc2, var_anc2, covT_anc1_anc2."),
  make_option("--outFile", type = "character", default = "",
    help = "Path to write the results table (tab-delimited). If empty, only printed."),
  make_option("--windowSizebp", type = "numeric", default = 0,
    help = "If >0, also run a sliding-window local scan with this window size in bp [default=0, whole-file only]."),
  make_option("--stepSizebp", type = "numeric", default = -1,
    help = "Step in bp between sliding windows [default=windowSizebp/2]."),
  make_option("--minVariantsWindow", type = "numeric", default = 30,
    help = "Minimum number of variants required to report a window [default=30]."),
  make_option("--regionStart", type = "numeric", default = -1,
    help = "Optional bp start of a candidate region; if set with --regionEnd, reports rho inside vs outside."),
  make_option("--regionEnd", type = "numeric", default = -1,
    help = "Optional bp end of a candidate region."),
  make_option("--nJackknifeBlocks", type = "numeric", default = 20,
    help = "Number of block-jackknife blocks for the SE [default=20].")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (!nzchar(opt$step2File)) {
  stop("--step2File is required.")
}

estimateCrossAncestryRg(
  step2File = opt$step2File,
  outFile = opt$outFile,
  windowSizebp = opt$windowSizebp,
  stepSizebp = opt$stepSizebp,
  minVariantsWindow = opt$minVariantsWindow,
  regionStart = opt$regionStart,
  regionEnd = opt$regionEnd,
  nJackknifeBlocks = opt$nJackknifeBlocks
)
