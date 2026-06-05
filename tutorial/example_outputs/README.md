# Example outputs

These files show the **format** of what SAIGE-Tractor produces, so users can verify their own runs by inspection.

| file | what it is |
|---|---|
| `step2_results.tsv` | a 9-variant excerpt of the Step 2 output for the simulated dataset in [`../example_inputs/`](../example_inputs/). The row with `MarkerID = v31` is the simulated causal variant (AFR-haplotype-specific effect), and the other rows are representative non-causal variants. |

The column structure of `step2_results.tsv` is the real SAIGE-Tractor admixed-mode quantitative-trait output (`K = 3` ancestries, no conditioning markers). For a full description of every column, see [section 8 of the tutorial README](../README.md#8-reading-the-output).

To re-create this file yourself, run the full pipeline as documented in the tutorial. The simulation seed is fixed (`SEED = 7` in `simulate_example.py`) so the per-variant effect sizes and p-values should be reproducible up to floating-point.
