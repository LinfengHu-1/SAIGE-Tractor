#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

## Step 3 (admixed): individual-level cross-ancestry genetic correlation (rg) and
## per-ancestry SNP heritability via RHE, combining the per-chunk partial
## accumulators written by step2 (run with --estimate_cross_anc_rg=TRUE).
##
## Each step2 run (per chromosome / chunk, all with the SAME --rg_seed) writes a
## <SAIGEOutputFile>.rg_partial.bin. This step sums them (additive) and solves the
## Haseman-Elston system for sigma2_anc, r_admix, h2, with a block-jackknife SE
## where each chunk is one block. See estimateCrossAncestryRgRHE.

options(stringsAsFactors = FALSE)
library(SAIGE)
library(optparse)
library(methods)

option_list <- list(
  make_option("--partialFiles", type = "character", default = "",
    help = "Comma-separated list of <prefix>.rg_partial.bin files (one per step2 chunk/chromosome, all run with the same --rg_seed). Alternatively use --partialDir + --partialPattern."),
  make_option("--partialDir", type = "character", default = "",
    help = "Directory to glob for partial files (used if --partialFiles is empty)."),
  make_option("--partialPattern", type = "character", default = "*.rg_partial.bin",
    help = "Glob pattern within --partialDir [default=*.rg_partial.bin]."),
  make_option("--outFile", type = "character", default = "",
    help = "Path to write the results table (tab-delimited). Per-ancestry h2 goes to <outFile>.h2."),
  make_option("--prevalence", type = "double", default = NA_real_,
    help = "Population prevalence K for a BINARY trait. Enables liability-scale per-ancestry h2 (Lee et al. 2011). rg is scale-free and needs no prevalence. If omitted for a binary trait, the sample case proportion is used (assumes no ascertainment).")
)

opt <- parse_args(OptionParser(option_list = option_list))

files <- character(0)
if (nzchar(opt$partialFiles)) {
  files <- trimws(strsplit(opt$partialFiles, ",")[[1]])
} else if (nzchar(opt$partialDir)) {
  files <- Sys.glob(file.path(opt$partialDir, opt$partialPattern))
}
files <- files[nzchar(files)]
if (length(files) == 0) {
  stop("No partial files. Provide --partialFiles (comma-separated) or --partialDir.")
}
missing <- files[!file.exists(files)]
if (length(missing) > 0) stop("Partial files not found: ", paste(missing, collapse = ", "))

cat(sprintf("Combining %d RHE partial file(s):\n  %s\n",
            length(files), paste(files, collapse = "\n  ")))

estimateCrossAncestryRgRHE(partialFiles = files, outFile = opt$outFile,
                           prevalence = opt$prevalence)
