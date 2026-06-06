# Cross-ancestry genetic correlation (rg) and per-ancestry heritability (h2)

Estimate, in **unrelated admixed individuals**, the cross-ancestry genetic correlation
`rg` of causal effects across local-ancestry backgrounds and the per-ancestry SNP
heritability `h2`, using a randomized Haseman-Elston (RHE) estimator on Tractor
ancestry-deconvolved dosages (the same estimand as `radmix`, generalized to K>=2
ancestries and to binary traits, and scalable to biobank N).

> **Scope:** this workflow is for **unrelated samples**. Remove related individuals first
> (Step 0). For related samples, relatedness is not yet handled in this estimator.

---

## Inputs

| input | description |
|---|---|
| packed `tractor_hybrid` genotypes | `<prefix>.meta/.common.*/.ancblock.*/.samples` (the Step-2 input produced by the tractor_hybrid converter). One `--tractorHybridPrefix` per chromosome (e.g. `chr1`). |
| phenotype file | tab-delimited, header with `IID`, the phenotype column, and covariates. |
| covariates | **must include ancestry PCs** (or per-individual global-ancestry fractions). With no GRM, PCs are what control population stratification. |
| unrelated sample list | one IID per line (Step 0 output). |
| (optional) marker lists | variant-id or `chr:pos` per line, to restrict rg and/or h2 to chosen marker sets. |

---

## Step 0 — Select unrelated samples (do this first)

Same criterion as `radmix`: exclude pairs with kinship above 0.0884 (< 3rd-degree).
Produce a one-IID-per-line "keep" list and feed it to Step 1/2 via `--sampleFile`.
Two ways:

### Option 1 — from genotypes (recommended for admixed data)

```bash
# plink2 (KING-robust -> ancestry-robust, correct for admixed cohorts)
plink2 --bfile mydata --king-cutoff 0.0884
#   -> plink2.king.cutoff.in.id   (UNRELATED samples to KEEP -- use this)
#   -> plink2.king.cutoff.out.id  (related samples removed)

# OR KING
king -b mydata.bed --unrelated --degree 3   # -> kingunrelated.txt
```

### Option 2 — from an existing sparse GRM

If you already built a sparse GRM for biobank GWAS, derive the unrelated set directly
(no need to re-run plink2). Greedy pruning = the same rule as `plink2 --king-cutoff`:

```bash
Rscript step0_extractUnrelatedFromGRM.R \
  --sparseGRMFile=sparseGRM.mtx \
  --sparseGRMSampleIDFile=sparseGRM.mtx.sampleIDs.txt \
  --relatednessCutoff=0.0884 \
  --outFile=unrelated.in.id
```

> **Caveats for the GRM route:** (1) the GRM kinship must be **ancestry-robust**
> (KING-robust). A standard GCTA-style GRM over-estimates relatedness in admixed samples
> (shared continental ancestry inflates apparent kinship) and will over-remove -- in that
> case use Option 1. (2) The GRM must have been built at a kinship cutoff <= 0.0884, or
> pairs just below the build cutoff won't be in the matrix.

Use the resulting `*.in.id` list as `--sampleFile` below.

---

## Workflow A — 3 steps (integrated: also produces the per-variant GWAS scan)

Use this when you want the standard admixed association results **and** rg/h2 from one
genotype pass.

### Step 1 — null model (GLMM-free, no GRM)

```bash
Rscript step1_fitNULL_noGRM.R \
  --phenoFile=pheno.txt --phenoCol=Y \
  --covarColList=age,sex,PC1,PC2,...,PC10 \
  --traitType=quantitative --invNormalize=TRUE \
  --sampleFile=plink2.king.cutoff.in.id \        # <-- UNRELATED samples to include
  --outputPrefix=null
# writes null.rda + null.varianceRatio.txt (variance ratio = 1; no GRM fit)
# binary: --traitType=binary (drop --invNormalize)
```

### Step 2 — association scan + RHE accumulation (per chromosome)

```bash
Rscript step2_SPAtests.R \
  --tractorHybridPrefix=chr1 --chrom=1 --LOCO=FALSE \
  --is_admixed=TRUE --number_of_ancestry=3 \
  --GMMATmodelFile=null.rda --varianceRatioFile=null.varianceRatio.txt \
  --SAIGEOutputFile=chr1.assoc.txt --minMAF=0 --minMAC=1 \
  --estimate_cross_anc_rg=TRUE \
  --rg_nProbes=30 --rg_seed=1 --rg_nJackknifeBlocks=20 \
  --rg_perAncestryH2=TRUE \
  --rg_markerFile=ld_pruned.snps \               # <-- MARKER LIST for rg (optional)
  --h2_markerFile=hapmap3.snps                    # <-- MARKER LIST for h2 (optional)
# writes chr1.assoc.txt (GWAS) + chr1.assoc.txt*.rg_partial.bin (RHE partials)
# samples are restricted to the unrelated set automatically (from null.rda).
# Repeat per chromosome with the SAME --rg_seed so partials are additive.
```

### Step 3 — combine partials, solve, report

```bash
Rscript step3_crossAncestryRgRHE.R \
  --partialDir=. --partialPattern="*.rg_partial.bin" \
  --outFile=out.rg.txt
  # binary: add --prevalence=0.01  (population prevalence; for liability-scale h2)
# writes out.rg.txt (per-pair rg) + out.rg.txt.h2 (per-ancestry h2)
```

---

## Workflow B — 1 script (standalone: rg/h2 only, ~3-4x faster)

Use this when you only need rg/h2 (no GWAS scan). It does the null model + RHE in one
pass; no separate Step 1/3.

```bash
Rscript step2_rgRHEonly.R \
  --tractorHybridPrefix=chr1 --chrom=1 \
  --phenoFile=pheno.txt --phenoCol=Y \
  --covarColList=age,sex,PC1,PC2,...,PC10 \
  --traitType=quantitative --invNormalize=TRUE \
  --sampleFile=plink2.king.cutoff.in.id \        # <-- UNRELATED samples
  --numberOfAncestry=3 \
  --rg_markerFile=ld_pruned.snps \               # <-- MARKER LIST for rg (optional)
  --h2_markerFile=hapmap3.snps \                  # <-- MARKER LIST for h2 (optional)
  --rg_perAncestryH2=TRUE \
  --rg_nProbes=60 --rg_nJackknifeBlocks=20 --rg_seed=1 \
  --outFile=out.rg.txt
# writes out.rg.txt + out.rg.txt.h2  (binary: --traitType=binary --prevalence=0.01)
```

---

## Output

`out.rg.txt` (one row per ancestry pair):

| column | meaning |
|---|---|
| `anc_a`, `anc_b` | ancestry indices (1-based) |
| `rg`, `rg_se` | genetic correlation + block-jackknife SE |
| `rg_pval` | two-sided test **H0: rg = 0** (is there *any* cross-ancestry correlation?) |
| `rg_pval_vs1` | one-sided test **H0: rg = 1** vs H1: rg < 1 (radmix "are causal effects *shared*?"; small p => effects **differ** across ancestries) |
| `cov_cross` | cross-ancestry genetic covariance (gamma) |
| `sigma2_anc_a/b` | per-ancestry genetic variance components |

`out.rg.txt.h2` (one row per ancestry):

| column | meaning |
|---|---|
| `h2_joint`, `h2_joint_se`, `h2_joint_pval` | per-ancestry h2 from the joint K-ancestry fit; p is one-sided H1: h2 > 0 |
| `h2_ownMarkers` | (if `--rg_perAncestryH2`) h2 on that ancestry's own marker set |
| `h2_*_liability` | (binary, with `--prevalence`) liability-scale h2 (Lee et al. 2011) |

---

## Key options

| option | step | meaning |
|---|---|---|
| `--sampleFile` | 1 / B | unrelated IID list (one per line) |
| `--covarColList` | 1 / B | covariates; **include ancestry PCs** |
| `--rg_markerFile` | 2 / B | markers used for **rg** (variant-id or `chr:pos`/line). rg is a scale-free ratio and tolerates aggressive **LD-pruning** (e.g. ~200k) -- cheaper. |
| `--h2_markerFile` | 2 / B | markers used for **h2**. h2 is tagging-limited; use a **fuller set** (e.g. HapMap3 ~1M). |
| `--rg_ancestry_list` | B | restrict to a subset of ancestries, e.g. `2,3` (1-based). |
| `--rg_perAncestryH2` | 2 / B | also report per-ancestry h2 on own markers (`h2_ownMarkers`). |
| `--rg_nProbes` | 2 / B | RHE random probes (more = less Monte-Carlo noise; 30-60 typical). |
| `--rg_nJackknifeBlocks` | 2 / B | block-jackknife blocks for the SE (>=2; ~20 for an SE from one run). |
| `--rg_seed` | 2 | probe seed; **identical across chromosomes** so partials are additive. |
| `--prevalence` | 3 / B | population prevalence for liability-scale binary h2. |

---

## Notes

- **Covariates must include ancestry PCs.** Without a GRM, PCs control stratification.
- **Marker sets:** LD-pruned for `rg` (scale-free, tolerates pruning); fuller/HapMap3 for
  `h2` (tagging-limited). Pass separate files to `--rg_markerFile` / `--h2_markerFile`.
- **Quantitative traits** are unbiased (validated against an exact-GRM reference).
- **Binary traits** run (PCGC residual + Lee liability scale) but the cross-ancestry rg
  carries a known liability-threshold attenuation, and per-ancestry h2 on the liability
  scale is approximate; interpret binary rg/h2 conservatively.
- **Multi-chromosome:** run Step 2 per chromosome with the SAME `--rg_seed`, then combine
  all `*.rg_partial.bin` in one Step 3 (each chunk is also a jackknife block).
