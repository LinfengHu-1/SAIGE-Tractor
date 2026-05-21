# Tractor Hybrid Utilities

designed by Kai, implemented by codex

This directory contains standalone utilities for the tractor_hybrid Step 2 input.

## Build

Install htslib development headers, then run:

```bash
make -C tools/tractor_hybrid
```

The Docker build installs these tools into `/scripts/bin`, matching the
`kyuan1024/saigetractor:1.4.9-tractor-hybrid.2` runtime layout.

## Tools

- `flare_subset_to_tractor_hybrid`: convert phased genotype VCF/BCF plus FLARE LAI VCF/BCF to tractor_hybrid files.
- `estimate_mac_threshold`: scan a genotype VCF/BCF and estimate a rare/common MAC threshold for packed storage.
