# Example outputs

These files are the **actual outputs** produced by running the full SAIGE-Tractor pipeline (Step 0 packer → Step 1 null GLMM → Step 2 admixed-mode association test) on the simulated data in `../example_inputs/` inside the `lhu1/saige-tractor:1.4.9-bundled` Docker image.  They are committed verbatim so users can diff their own runs against a fixed reference.

| file | what it is | size |
|---|---|---|
| `step1_example.rda` | fitted null GLMM, binary R serialization (loaded by Step 2 via `--GMMATmodelFile`). Reproducible inside the same Docker image. | ~35 KB |
| `step1_example.varianceRatio.txt` | the variance ratio `r̂`. Single number `0.645358852038115 null 1`. Passed to Step 2 via `--varianceRatioFile`. | < 1 KB |
| `step1_example.log` | full console log of the Step 1 run, including the GLMM iterations and the variance-ratio CV. Useful for diagnosing your own run when it diverges. | ~24 KB |
| `step2_results.tsv` | the actual Step 2 association-test output. **71 columns**, 52 variants (8 of the 60 dropped by per-ancestry MAC ≥ 1 filtering). | ~28 KB |
| `step2_example.log` | full console log of the Step 2 run, including the per-chunk timing breakdown emitted by the `tractor_hybrid` profiling instrumentation. | ~5 KB |

## What's in `step2_results.tsv`

For a full column-by-column description see [section 8 of the tutorial README](../README.md#8-reading-the-output). The headline row to look at is the simulated causal variant:

```
MarkerID    BETA_anc1   p.value_anc1   BETA_c_anc1   p.value_c_anc1   Pvalue_haplo_anc1   P_cct_admixed_c
v31_G_A     +1.527      4.46e-13       +1.527        4.46e-13         1.51e-04             5.82e-12
```

The simulated truth was `β_AFR = 1.4`; the recovered AFR-haplotype effect is `+1.53` (`p = 4.46e-13`), per-ancestry. The joint p-value is `5.82e-12`. `Pvalue_haplo_anc1 = 1.51e-04` exceeds the configured `--pvalcutoff_of_haplotype = 0.000005`, so marker-haplotype conditioning did not trigger and the `_c` columns mirror the unconditional ones — both are correct to report. To see triggered `_c` columns on this toy data, re-run Step 2 with `--pvalcutoff_of_haplotype=0.001`.

## Reproducibility

To verify these outputs from scratch:

```bash
cd ../
python3 simulate_example.py                                     # regenerate inputs
mkdir -p packed output

# Step 0a — bgzip + tabix (in docker; htslib 1.10.2 CLIs bundled in 1.4.9-bundled)
docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work \
    lhu1/saige-tractor:1.4.9-bundled bash -lc "
        bgzip -c example_inputs/genotypes.phased.vcf > example_inputs/genotypes.phased.vcf.gz
        bgzip -c example_inputs/localanc.flare.vcf  > example_inputs/localanc.flare.vcf.gz
        tabix -p vcf example_inputs/genotypes.phased.vcf.gz
        tabix -p vcf example_inputs/localanc.flare.vcf.gz
    "

# Step 0b — pack (in docker)
docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work \
    lhu1/saige-tractor:1.4.9-bundled \
    /scripts/bin/flare_subset_to_tractor_hybrid \
      example_inputs/genotypes.phased.vcf.gz \
      example_inputs/localanc.flare.vcf.gz \
      3 6 packed/example

# Step 1 — null GLMM
docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work \
    lhu1/saige-tractor:1.4.9-bundled \
    step1_fitNULLGLMM.R \
      --plinkFile=example_inputs/example_grm \
      --phenoFile=example_inputs/phenotype.tsv \
      --phenoCol=Y --traitType=quantitative --invNormalize=FALSE \
      --covarColList=age,sex,PC1,PC2,PC3 --qCovarColList=sex \
      --sampleIDColinphenoFile=IID \
      --outputPrefix=output/step1_example \
      --LOCO=FALSE --IsOverwriteVarianceRatioFile=TRUE

# Step 2 — association test
docker run --rm --platform linux/amd64 -v "$PWD":/work -w /work \
    lhu1/saige-tractor:1.4.9-bundled \
    step2_SPAtests.R \
      --tractorHybridPrefix=packed/example \
      --chrom=chr22 \
      --is_admixed=TRUE --number_of_ancestry=3 \
      --pvalcutoff_of_haplotype=0.000005 \
      --GMMATmodelFile=output/step1_example.rda \
      --varianceRatioFile=output/step1_example.varianceRatio.txt \
      --SAIGEOutputFile=output/step2_results.tsv \
      --LOCO=FALSE --minMAF=0 --minMAC=1

diff <(awk '{print $1, $2, $3}' output/step2_results.tsv) \
     <(awk '{print $1, $2, $3}' example_outputs/step2_results.tsv)
```

If your `diff` shows no output on the first three columns (CHR / POS / MarkerID), your run produced the same variants in the same order. Exact numerical values may differ at the 6th-or-later significant digit because of BLAS-level non-determinism across hosts, but the order of magnitude and the pattern of which `_c` columns differ from unconditional ones should match.
