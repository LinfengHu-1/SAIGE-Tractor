# Example inputs

Pre-generated, fully simulated input files for the tutorial. All produced by [`../simulate_example.py`](../simulate_example.py) with a fixed seed; **no real human data**.

| file | format | shape |
|---|---|---|
| `samples.txt`               | one ID per line, no header             | 200 |
| `phenotype.tsv`             | tab-separated, header: `IID, age, sex, PC1, PC2, PC3, Y` | 200 × 7 |
| `genotypes.phased.vcf`      | VCF v4.2, phased, single `GT` field    | 60 variants × 200 samples |
| `localanc.flare.vcf`        | VCF v4.2, FLARE format with `GT:AN1:AN2` | 60 variants × 200 samples |
| `example_grm.{bed,bim,fam}` | PLINK 1 binary                         | 60 variants × 200 samples |

To regenerate everything from scratch:

```bash
cd ..
python3 simulate_example.py
```

The phased VCF and the FLARE VCF use **identical sample order and variant coordinates** — this is required by the `tractor_hybrid` packer. The simulation produces one **AFR-haplotype-specific causal variant** at position `chr22:16150000` (`MarkerID = v31`) with simulated `BETA_AFR = 1.4`; all other variants are null.

The PLINK files contain the same 60 variants and are used only by Step 1 for variance-ratio estimation. For real biobank data you would typically use a larger, separately-genotyped marker set here.
