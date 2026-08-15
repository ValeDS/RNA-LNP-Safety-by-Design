# RNA–LNP Safety-by-Design framework

Code and synthetic data supporting the RNA–LNP prioritization and biological-consistency analyses described in the associated BMC Bioinformatics manuscript.

## Scope

This repository is a computational proof of concept. The response equations, virtual population, formulation space and Safety-by-Design preferences are modeling assumptions. Random Forest performance measures recovery of the synthetic response surface. BPCI evaluates immune-program organization and does not validate individual formulations experimentally or clinically.

## Contents

- `scripts/core/`: generation, scoring, grouped cross-validation, Pareto prioritization and transcriptomic consistency analyses.
- `scripts/revision/`: null, sensitivity, subgroup, sampling, covariance and balanced-design analyses.
- `data/inputs/`: intermediate inputs needed to rerun the analyses without local paths.
- `data/primary/`: principal outputs from the 10,000-record analysis.
- `data/revision/`: machine-readable outputs from the additional analyses, including complete risk-stratified distribution summaries and formulation rankings.
- `figures/supplementary/`: supplementary figures in PNG and PDF.
- `tables/`: machine-readable supplementary tables and captions.

## Requirements

The analyses were run with R 4.4.1. Package versions are listed in `SESSION_INFO.txt`. The main dependencies are `dplyr`, `tidyr`, `ggplot2`, `pheatmap`, `ranger`, `patchwork` and `scales`.

## Running the analyses

Run commands from the repository root.

```bash
Rscript scripts/run_core_pipeline.R
Rscript scripts/run_revision_analyses.R
```

Outputs are written to `results/core`, `results/revision` and `results/figures`. These directories can be changed without editing the scripts by setting `RNA_LNP_CORE_OUTPUT_DIR`, `RNA_LNP_REVISION_OUTPUT_DIR` or `RNA_LNP_REVISION_FIGURE_DIR`.

The core runner executes modules 01–06. The GSE171964 single-cell processing script is separate because it requires the public GEO source files. Set `GSE171964_DATA_DIR` to the directory containing those files. Derived inputs required for the BPCI analyses are included under `data/inputs`.

## Analysis designs

The primary dataset contains 10,000 independent draws with replacement from 200 virtual patients and 500 formulations. It contains 9,509 unique patient–formulation pairs. The balanced sensitivity analysis evaluates all 100,000 pairs once and uses the primary normalization parameters.

Key checks reproduced by the release are:

- exact BPCI label-permutation test: mean BPCI 0.760, exact upper-tail p = 0.025;
- 3,876 Safety-by-Design weight combinations;
- 21 primary Pareto-optimal candidates;
- balanced versus primary rank correlation of 0.9981;
- LNP_0462 ranked first globally and in all three balanced inflammatory-risk strata;
- complete endpoint and Safety-by-Design score distributions for the Low-, Intermediate- and High-risk strata;
- maximum absolute cross-layer residual correlation of 0.097 after independent noise generation.

## Data provenance

All virtual-patient, formulation and response data are synthetic. GSE171964 is a public human single-cell RNA-sequencing dataset; raw GEO files are not redistributed here. File checksums are recorded in `SHA256SUMS.txt`.

## License and citation

Code is released under the MIT License. Synthetic data, figures and documentation are released under CC BY 4.0. Version 1.0.2 is archived on Zenodo under DOI [10.5281/zenodo.21944177](https://doi.org/10.5281/zenodo.21944177); full citation metadata are provided in `CITATION.cff`.
