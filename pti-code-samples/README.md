# PTI Code Samples — Single-Cell & Proteomics

Well‑annotated examples of my hands‑on work in single‑cell RNA‑seq (Scanpy/scVI) and quantitative proteomics (DIA/TMT with DIA‑NN/FragPipe → MSstats/limma). 
The goal is to show **reproducible**, **pipeline‑driven** analysis with clear documentation and interactive review tools.

## Contents
- `single_cell/`
  - `sc_pipeline.py`: End‑to‑end single‑cell pipeline skeleton with QC, normalization, batch correction, clustering, annotation, DEG/TF analysis, and figure export.
  - `env/environment.yml`: Reproducible environment spec.
- `proteomics/`
  - `diann_workflow.py`: DIA‑NN command orchestration + tidy outputs; hooks for FragPipe/MaxQuant.
  - `msstats_example.R`: MSstats/limma example for LFQ/TMT with contrasts and QC.
  - `proteomics_pipeline.md`: Overview and rationale.
- `pipelines/`
  - `Snakefile`: Snakemake DAG wiring scRNA‑seq and proteomics tasks with config‑driven parameters.
  - `config/config.yaml`: Parameters for both stacks.
- `apps/`
  - `dash_app.py`: Minimal Dash app to interactively browse markers/DE/QC.

## Reproducibility
- Environments are pinned; scripts are pure‑Python/R with docstrings.
- Data paths are parameterized via `pipelines/config/config.yaml`.
- Each step writes versioned outputs + logs.

## License
MIT — feel free to adapt.

— Generated 2025-08-26
