# Cross-ancestry rg / per-ancestry h2 in *All of Us* (AoU)

End-to-end recipe for estimating, in admixed AoU cohorts (**Researcher Workbench**, controlled tier), the
cross-ancestry genetic correlation `rg` and per-ancestry SNP heritability `h2`. Covers **two-way**
(African-ancestry admixed, AFR-EUR) and **K>2** (e.g. 3-way Latino, AFR-EUR-AMR). Read
[`README_crossAncestryRg_h2.md`](README_crossAncestryRg_h2.md) first for option-level detail; this doc is
the AoU-specific wrapper.

> **This workflow estimates rg/h2 only — it does NOT run a GWAS.** It uses the standalone RHE path
> (`step2_rgRHE_chrom` → `step3_rgRHE_combine`), which fits its own covariate-only (no-GRM) null
> internally — there is **no step1, no association scan, no `step2_SPAtests`**. (If you separately want
> a Tractor GWAS, run that on the full cohort as usual; keep it independent of this rg/h2 analysis.)

Image: **`wzhou88/saigetractor:rg-h2`**. Pin the digest for a reproducible/citeable run:
`docker pull wzhou88/saigetractor:rg-h2 && docker inspect --format '{{index .RepoDigests 0}}' wzhou88/saigetractor:rg-h2`

---

## TL;DR pipeline (rg/h2 only — no GWAS)

```
phased WGS  +  FLARE local ancestry          (you produce these upstream, per chromosome)
        │
        ▼  flare_subset_to_tractor_hybrid          → packed tractor_hybrid genotypes (per chrom)
        │
        ▼  unrelated subset (KING-robust, kinship<0.0884)   + ancestry PCs as covariates
        │
        ▼  step2_rgRHE_chrom  (MAP: one RHE partial per chrom, no GWAS, no separate null)
        │
        ▼  step3_rgRHE_combine  (REDUCE: sum partials, solve, leave-one-chromosome-out jackknife)
        │
        ▼  out.rg.txt (per-pair rg)  +  out.rg.txt.h2 (per-ancestry h2)
```

---

## 0. Inputs you produce upstream (NOT part of this tool)

| input | how, in AoU |
|---|---|
| **Phased** genotypes, per chrom | AoU srWGS → phase (Eagle/SHAPEIT). Phased `GT` is required (Tractor is haplotype-aware). |
| **Local ancestry** calls, per chrom | **FLARE** (recommended) or RFMix/Tractor, with reference panels (e.g. 1000G/HGDP) for your K continental ancestries. Output = per-haplotype ancestry VCF (`ANC1`/`ANC2` fields). |
| phenotype + covariates table | `IID`, the trait column, age, sex, and **AoU genetic PCs** (PC1..PCk). |
| relatedness / PCs | AoU provides KING relatedness flags and genetic PCs in the workbench (the `ancestry_preds` / relatedness resources). |

> **K must match your FLARE run.** Use `K=2` for AFR-EUR, `K=3` for AFR-EUR-AMR. The tool reports
> ancestry **indices** (1,2,…) in FLARE's component order — record which continental ancestry each
> index is.

---

## 1. Convert phased + FLARE → packed `tractor_hybrid` (per chromosome)

```bash
# one per chromosome; <K> = number of ancestries; <minMAC> e.g. 1 (or use estimate_mac_threshold)
flare_subset_to_tractor_hybrid \
  chr${C}.phased.vcf.gz  chr${C}.flare.anc.vcf.gz  ${K}  1  chr${C}
#   → chr${C}.meta / chr${C}.common.* / chr${C}.ancblock.* / chr${C}.samples
```

This packed format is the input to the rg/h2 MAP step below, built once per cohort. (The same converter
output also feeds a Tractor GWAS if you run one separately — but that is independent of this rg/h2 analysis.)

## 2. Select unrelated samples (rg/h2 requires unrelated)

The RHE estimator does **not** model relatedness — subset to unrelated first (same rule as `radmix`,
kinship < 0.0884 / 3rd-degree). In AoU:

```bash
# Preferred: KING-robust (ancestry-robust; correct for admixed) on an LD-pruned set
plink2 --bfile aou_pruned --king-cutoff 0.0884      # → plink2.king.cutoff.in.id  (KEEP list)
# OR derive from AoU's relatedness table / an existing sparse GRM:
step0_extractUnrelatedFromGRM.R --sparseGRMFile=... --relatednessCutoff=0.0884 --outFile=unrelated.in.id
```

> Do **not** use a standard GCTA GRM cutoff in admixed AoU data — shared continental ancestry inflates
> apparent kinship and over-removes. Use KING-robust.

## 3. Covariates & marker sets

- **Covariates:** `age,sex,PC1,…,PC10` (AoU genetic PCs). With no GRM, PCs control stratification.
- **Markers:** LD-pruned set (~200–500k) for **rg** (`--rg_markerFile`, scale-free, tolerates pruning);
  a fuller set (HapMap3 ~1M) for **h2** (`--h2_markerFile`, tagging-limited).

## 4. Run the estimator (rg/h2 only — no GWAS)

Two steps: a per-chromosome **MAP** (`step2_rgRHE_chrom`, writes one `.rds` partial per chrom — it fits its
own covariate-only null, so there is no separate step1 and no association scan) and a genome-wide
**REDUCE** (`step3_rgRHE_combine`, sums the partials and solves with a leave-one-chromosome-out jackknife).
**All chromosomes must use the SAME `--rg_seed`** (so the random probes are identical and the partials are
additive) and the same `--rg_nProbes` / `--covarColList` / `--maf`.

### (a) Docker image directly in `dsub` / Cromwell-WDL  — recommended on the AoU cloud

AoU batch (`dsub`, Cromwell) pulls Docker images natively — **no Singularity conversion needed**.

```bash
# --- MAP: one RHE partial per chromosome (scatter over C=1..22). No GWAS, no step1. ---
dsub --provider google-cls-v2 --project "$GOOGLE_PROJECT" --regions us-central1 \
  --image wzhou88/saigetractor:rg-h2 --machine-type n1-highmem-4 \
  --input GENO=gs://$WS/packed/chr${C}.'*' PHENO=gs://$WS/pheno.txt \
          KEEP=gs://$WS/plink2.king.cutoff.in.id \
          RGM=gs://$WS/ld_pruned.snps H2M=gs://$WS/hapmap3.snps \
  --output PART=gs://$WS/rg_partials/chr${C}.rds \
  --command 'step2_rgRHE_chrom.R --tractorHybridPrefix=chr'${C}' --chrom='${C}' --chromIndex='${C}' \
       --numberOfAncestry='${K}' \
       --phenoFile=${PHENO} --phenoCol=Y --covarColList=age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
       --traitType=quantitative --invNormalize=TRUE --sampleFile=${KEEP} \
       --rg_markerFile=${RGM} --h2_markerFile=${H2M} --rg_perAncestryH2=TRUE \
       --rg_nProbes=30 --rg_seed=1 --outFile=${PART}'
       #   binary trait: --traitType=binary (drop --invNormalize)

# --- REDUCE: combine all per-chromosome partials (LOCO jackknife) ---
dsub ... --image wzhou88/saigetractor:rg-h2 \
  --input PARTIALS=gs://$WS/rg_partials/'*.rds' \
  --output OUT=gs://$WS/rg/'*' \
  --command 'step3_rgRHE_combine.R --partialGlob="$(dirname ${PARTIALS})/*.rds" \
       --outFile=$(dirname ${OUT})/out.rg.txt'
       #   binary trait: add --prevalence=<K>  (population prevalence; for liability-scale h2)
```

(Paths/`dsub` flags are a template — adapt to your workspace bucket `$WS` and the AoU `dsub` wrapper.
`--number_of_ancestry`/`-K` must match your FLARE run: 2 for AFR-EUR, 3 for a 3-way cohort.)

### (b) Singularity, on an HPC-style node

```bash
docker save wzhou88/saigetractor:rg-h2 | gzip > st.tar.gz
gunzip st.tar.gz && apptainer build st.sif docker-archive://st.tar
# MAP per chrom, then REDUCE:
apptainer exec --bind "$PWD:/work" st.sif step2_rgRHE_chrom.R --tractorHybridPrefix=chr1 --chrom=1 \
  --chromIndex=1 --numberOfAncestry=$K --phenoFile=pheno.txt --phenoCol=Y \
  --covarColList=age,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleFile=unrelated.in.id \
  --rg_markerFile=ld_pruned.snps --h2_markerFile=hapmap3.snps --rg_perAncestryH2=TRUE \
  --rg_nProbes=30 --rg_seed=1 --outFile=chr1.rds       # ...repeat per chrom
apptainer exec --bind "$PWD:/work" st.sif step3_rgRHE_combine.R \
  --partialGlob="*.rds" --outFile=out.rg.txt
```

> **One-chromosome quick check:** `step2_rgRHEonly.R` runs the whole thing (null + RHE + solve) in a
> single pass for one `--tractorHybridPrefix`, writing `out.rg.txt` directly — handy for a single
> chromosome or a fast sanity check before the genome-wide scatter.

## 5. Read the output

`out.rg.txt` — one row per ancestry pair (`K(K-1)/2` rows):
- **`rg_constrained`** = the value to report (clipped to [-1,1]); `rg` is the raw estimate.
- **`rg_se_delta`** = robust SE — prefer it over `rg_se` for low-h2 / minority-ancestry pairs.
- **`rg_pval_vs1`** = radmix "are causal effects *shared*?" (H0: rg=1); small p ⇒ effects **differ**.
- **`rg_converged`** = probe-convergence check: `TRUE` if the Monte-Carlo noise from the probe count
  (`--rg_nProbes`) is negligible vs the sampling SE. If `FALSE` for a pair, re-run the MAP with a larger
  `--rg_nProbes` (keep it identical across all chromosomes). `rg_mcse_probe` is the probe MC-SE itself.

`out.rg.txt.h2` — one row per ancestry: `h2_joint` (+ se/pval), `h2_flag` (precision), liability-scale h2
for binary (with `--prevalence`).

### 2-way vs K>2 in AoU
- **AFR-EUR (K=2):** one `rg`. Both ancestries well powered at AoU AFR-ancestry N (~50k+). The headline
  numbers (single rg, two h2) are robust.
- **3-way Latino (K=3):** three pairs from one joint fit. The **smallest-exposure ancestry** (often AMR)
  is the one to scrutinize — its per-ancestry `h2` is high-variance (effective N ≈ N × mean local-ancestry
  exposure). Check its `h2_flag` is `ok` (not `imprecise`/`out_of_range`), prefer `rg_se_delta`, and lean on
  `rg` — a ratio, more robust to a thin minority ancestry than that ancestry's `h2`.

## 6. Scale, compute & policy notes

- **Scale:** RHE is `O(N·M·B)`, **linear in N** (no N×N GRM), so AoU N is fine. Memory is modest
  (`O(N·B)` for probes + accumulators).
- **Parallelism:** scatter Step 2 over the 22 autosomes; gather in Step 3. Each chromosome is one
  jackknife block ⇒ a chromosome-level (LD-aware) SE comes for free. Keep `--rg_seed` identical across
  all chromosomes (the per-chrom partials must use the same probes to be additive).
- **Reproducibility:** pin the image by digest for a published analysis —
  `docker pull wzhou88/saigetractor:rg-h2 && docker inspect --format '{{index .RepoDigests 0}}' wzhou88/saigetractor:rg-h2`
  and cite that `@sha256:…`.
- **Egress:** only the summary `out.rg.txt` / `.h2` tables (aggregate estimates) leave the workbench,
  per AoU policy — no individual-level data.
