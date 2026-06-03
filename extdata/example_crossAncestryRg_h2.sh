#!/usr/bin/env bash
# =============================================================================
# Example: cross-ancestry genetic correlation (rg) and per-ancestry SNP
# heritability (h2) in SAIGE-Tractor, via the RHE (randomized Haseman-Elston)
# estimator integrated into the admixed step2 pass.
#
# WHAT IT DOES
#   The estimator runs INSIDE the admixed step2 score-test loop (one genotype
#   read). step2 streams the ancestry-deconvolved dosages, builds randomized
#   trace statistics for ancestry variance-component GRMs, and writes a small
#   per-chunk "partial" sidecar. step3 combines the partials and solves the
#   Haseman-Elston system -> per-pair rg and per-ancestry h2 (+ jackknife SE,
#   Wald p-values, and a liability-scale h2 for binary traits).
#
#   It never forms an N x N matrix; cost is O(N * M * B) time, O(N * B) memory,
#   so it scales to biobank N. See doc/STATS_cross_ancestry_rg_h2.md (local) for
#   the statistical details and the connection to SAIGE's own machinery.
#
# PREREQUISITES
#   * a fitted step1 null model (.rda + .varianceRatio.txt). For correct
#     de-confounding, include global-ancestry PCs / fractions among the step1
#     covariates (the estimator covariate-adjusts the ancestry components with
#     exactly the step1 covariates).
#   * a packed tractor_hybrid step2 input (the converter output prefix).
# =============================================================================
set -euo pipefail

SAIGE_DIR="${SAIGE_DIR:?set SAIGE_DIR to the SAIGE-Tractor package dir}"
PREFIX="${1:?usage: example_crossAncestryRg_h2.sh <tractorHybridPrefix> <K>}"
K="${2:?number of ancestries}"
STEP1="${PREFIX}.step1"          # step1 output prefix (<>.rda, <>.varianceRatio.txt)
OUT="${PREFIX}.step2.rg"         # step2 association + rg partial output prefix

# Optional marker subsets (recommended at scale). rg is a scale-free ratio and
# tolerates aggressive LD-pruning; h2 is tagging-limited and wants a fuller set.
RG_MARKERS="${RG_MARKERS:-}"     # e.g. an LD-pruned ~200k list (one id or chr:pos/line)
H2_MARKERS="${H2_MARKERS:-}"     # e.g. a HapMap3 ~1M list

# -----------------------------------------------------------------------------
# STEP 2 — association test + RHE accumulation (single genotype pass)
# -----------------------------------------------------------------------------
# Key rg/h2 options:
#   --estimate_cross_anc_rg=TRUE  turn on RHE accumulation (off by default)
#   --rg_nProbes=30               # random probes (30-50 is plenty at scale)
#   --rg_seed=1                   # MUST be identical across all chunks/chroms
#   --rg_nJackknifeBlocks=20      # split this run into 20 jackknife blocks so a
#                                 #   single job yields an SE (else SE is NA)
#   --rg_pairs="1-2,1-3"          # estimate rg COHERENTLY per listed pair on the
#                                 #   markers shared by both ancestries (M_ab).
#                                 #   Omit (default) for the joint all-K analysis.
#   --rg_perAncestryH2=TRUE       # also estimate each ancestry's h2 on its OWN
#                                 #   markers (the proper per-ancestry heritability)
#   --rg_markerFile / --h2_markerFile  restrict rg / h2 to a marker subset
Rscript "${SAIGE_DIR}/extdata/step2_SPAtests.R" \
    --tractorHybridPrefix="${PREFIX}" \
    --chrom=1 --LOCO=FALSE \
    --is_admixed=TRUE --number_of_ancestry="${K}" \
    --GMMATmodelFile="${STEP1}.rda" \
    --varianceRatioFile="${STEP1}.varianceRatio.txt" \
    --SAIGEOutputFile="${OUT}.assoc.txt" \
    --minMAF=0 --minMAC=1 \
    --estimate_cross_anc_rg=TRUE \
    --rg_nProbes=30 \
    --rg_seed=1 \
    --rg_nJackknifeBlocks=20 \
    --rg_perAncestryH2=TRUE \
    $( [ -n "${RG_MARKERS}" ] && echo --rg_markerFile="${RG_MARKERS}" ) \
    $( [ -n "${H2_MARKERS}" ] && echo --h2_markerFile="${H2_MARKERS}" )
    # add --rg_pairs="1-2,1-3,2-3" to switch from joint all-K to per-pair coherent rg

# step2 wrote <OUT>.assoc.txt (association results) and, because
# --rg_nJackknifeBlocks=20, twenty <OUT>.assoc.txt.block{0..19}.rg_partial.bin
# sidecars (one .rg_partial.bin if blocks=1).
#
# BIOBANK / MULTI-CHROM: run the block above once per chromosome with the SAME
# --rg_seed; step3 below pools every chunk's partials (the chunks themselves also
# act as jackknife blocks). Per-process memory stays at O(N * B).

# -----------------------------------------------------------------------------
# STEP 3 — combine the partials and solve for rg + h2
# -----------------------------------------------------------------------------
#   --prevalence=K   (binary only) report per-ancestry h2 on the liability scale
Rscript "${SAIGE_DIR}/extdata/step3_crossAncestryRgRHE.R" \
    --partialDir="$(dirname "${OUT}")" \
    --partialPattern="$(basename "${OUT}").assoc.txt*.rg_partial.bin" \
    --outFile="${OUT}.estimates.txt"
    # for a binary trait add e.g. --prevalence=0.05

# OUTPUTS
#   <OUT>.estimates.txt       per-pair rg table:
#       anc_a anc_b rg rg_se rg_z rg_pval cov_cross sigma2_anc_a sigma2_anc_b ...
#       rg_pval = two-sided Wald test of rg != 0.
#   <OUT>.estimates.txt.h2    per-ancestry heritability table:
#       ancestry [h2_joint] [h2_ownMarkers] h2 h2_se h2_pval [h2_liability ...]
#       h2_pval = one-sided (boundary) test of h2 > 0.
echo "Done. rg table: ${OUT}.estimates.txt ; h2 table: ${OUT}.estimates.txt.h2"
