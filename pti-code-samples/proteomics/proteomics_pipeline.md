# Proteomics Pipeline Overview

This folder contains a documented skeleton for DIA/TMT quantitative proteomics:

1. **Search/Quant** (DIA-NN, FragPipe/Philosopher, or MaxQuant)
2. **Tidy outputs**: long-format peptide/protein tables with metadata
3. **Stats**: MSstats or limma for LFQ/TMT with contrasts, missingness handling
4. **QC**: ID rates, CVs, sample–sample distances, batch drift
5. **Integration**: join with RNA (gene-level), pathways/networks

Adjust paths via `pipelines/config/config.yaml`.
