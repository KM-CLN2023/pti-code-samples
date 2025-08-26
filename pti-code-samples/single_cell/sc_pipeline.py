#!/usr/bin/env python3
"""
Single-cell RNA-seq pipeline skeleton using Scanpy + scvi-tools.

This script demonstrates the **structure** and **documentation** of a typical
QC → normalization → batch correction → clustering/annotation → DEG/TF → figures pipeline.
It is parameterized by a YAML config (see pipelines/config/config.yaml).

Notes:
- Replace example paths with your data. For PTI, swap gene modules & markers as needed.
- Uses a 'frozen' UMAP pattern to keep figures stable across re-runs.
- Designed to be called from Snakemake, but can run standalone too.

Author: Krishanu Mukherjee
"""
import os, sys, json, yaml, logging, numpy as np, pandas as pd
import scanpy as sc

# Optional imports (install via environment.yml)
try:
    import scvi
except Exception:
    scvi = None

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')

def load_config(cfg_path: str) -> dict:
    with open(cfg_path) as f:
        return yaml.safe_load(f)

def qc_and_norm(adata, cfg):
    sc.pp.filter_cells(adata, min_genes=cfg['qc']['min_genes'])
    sc.pp.filter_genes(adata, min_cells=cfg['qc']['min_cells'])
    adata.var['mt'] = adata.var_names.str.upper().str.startswith(('MT-', 'MT_'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # Simple filters for demo; real thresholds come from EDA
    adata = adata[adata.obs.n_genes_by_counts < cfg['qc']['max_genes'], :].copy()
    adata = adata[adata.obs.pct_counts_mt < cfg['qc']['max_mt_pct'], :].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=cfg['hvg']['n_top_genes'], subset=True)
    return adata

def batch_correct(adata, cfg):
    method = cfg['batch']['method']
    if method == 'scvi':
        assert scvi is not None, 'scvi-tools not installed'
        adata.layers['counts'] = adata.X.copy()
        scvi.model.SCVI.setup_anndata(adata, batch_key=cfg['batch']['key'])
        model = scvi.model.SCVI(adata, n_latent=cfg['batch']['n_latent'])
        model.train(max_epochs=cfg['batch']['epochs'])
        adata.obsm['X_scVI'] = model.get_latent_representation()
        sc.pp.neighbors(adata, use_rep='X_scVI')
    elif method == 'harmony':
        # Placeholder: in practice use harmonypy or scanpy.external.pp.harmony_integrate
        sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    else:
        sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata, random_state=42)
    adata.obsm['X_umap_frozen'] = adata.obsm['X_umap'].copy()
    return adata

def cluster_annotate(adata, cfg):
    sc.tl.leiden(adata, resolution=cfg['cluster']['resolution'], key_added=cfg['cluster']['key'])
    # Marker-based quick annotation (demo)
    for label, genes in cfg.get('markers', {}).items():
        present = [g for g in genes if g in adata.var_names]
        adata.obs[f'module_{label}'] = adata[:, present].X.mean(axis=1) if present else 0
    return adata

def deg(adata, cfg):
    key = cfg['cluster']['key']
    sc.tl.rank_genes_groups(adata, key, method='wilcoxon')
    return adata

def save_figures(adata, cfg, outdir):
    os.makedirs(outdir, exist_ok=True)
    sc.pl.umap(adata, color=[cfg['cluster']['key']], legend_loc='on data', frameon=False, save=None, show=False)
    fig = sc.pl.umap(adata, color=list(cfg.get('markers', {}).keys()), show=False, return_fig=True)
    fig.savefig(os.path.join(outdir, 'umap_modules.png'), dpi=300, bbox_inches='tight')

def main(cfg_path: str):
    cfg = load_config(cfg_path)
    adata = sc.read(cfg['input']['adata'])
    logging.info('Loaded AnnData: %s', adata)
    adata = qc_and_norm(adata, cfg)
    adata = batch_correct(adata, cfg)
    adata = cluster_annotate(adata, cfg)
    adata = deg(adata, cfg)
    out = cfg['output']['adata']
    adata.write(out)
    save_figures(adata, cfg, cfg['output']['figdir'])
    logging.info('Done. Wrote %s', out)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: sc_pipeline.py pipelines/config/config.yaml')
        sys.exit(1)
    main(sys.argv[1])
