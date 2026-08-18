#!/usr/bin/env python3
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--mode', required=True)
    p.add_argument('--results', required=True)
    p.add_argument('--matrix_dir', default='VMR_matrix')
    p.add_argument('--qc', required=True)
    p.add_argument('--clinical', required=True)
    p.add_argument('--out_subdir', default='methylvi_casebatch_v1')
    p.add_argument('--batch_key', default='case_barcode')
    p.add_argument('--n_latent', type=int, default=10)
    p.add_argument('--max_epochs', type=int, default=250)
    p.add_argument('--batch_size', type=int, default=128)
    p.add_argument('--min_cells_per_feature', type=int, default=10)
    p.add_argument('--max_features', type=int, default=5000)
    p.add_argument('--n_neighbors', type=int, default=30)
    p.add_argument('--min_dist', type=float, default=0.10)
    p.add_argument('--resolution', type=float, default=0.12)
    p.add_argument('--seed', type=int, default=1234)
    p.add_argument('--devices', default='auto')
    p.add_argument('--accelerator', default='auto')
    return p.parse_args()


def clean(x):
    return pd.Series(x, dtype='string').astype(str).str.replace('"', '', regex=False).str.strip().tolist()


def fread_gz(path):
    return pd.read_csv(path, compression='infer')


def orient_table(df, qc_cells):
    first = [str(x).replace('"', '').strip() for x in df.iloc[:, 0].tolist()]
    cn_rest = [str(x).replace('"', '').strip() for x in df.columns[1:]]
    first_matches = sum(x in qc_cells for x in first)
    col_matches = sum(x in qc_cells for x in cn_rest)

    if col_matches >= first_matches:
        feature_ids = pd.Index(first).astype(str)
        mat = df.iloc[:, 1:].copy()
        mat.index = feature_ids
        mat.columns = cn_rest
        mat = mat.T
    else:
        cell_ids = pd.Index(first).astype(str)
        mat = df.iloc[:, 1:].copy()
        mat.index = cell_ids
        mat.columns = cn_rest
    return mat


def build_meta(qc_file, clinical_file, cells_keep):
    qc = pd.read_csv(qc_file, sep='\t')
    qc['cell_barcode'] = clean(qc['cell_barcode'])
    qc['case_barcode'] = clean(qc['case_barcode'])

    clin = pd.read_csv(clinical_file, sep='\t')
    clin['case_barcode'] = clean(clin['case_barcode'])
    clin['idh_status'] = np.where(clin['idh_codel_subtype'].astype(str) == 'IDHwt', 'IDHwt', 'IDHmut')

    meta = (
        qc[['cell_barcode', 'case_barcode']]
        .drop_duplicates('cell_barcode')
        .merge(
            clin[['case_barcode', 'idh_status', 'idh_codel_subtype']].drop_duplicates('case_barcode'),
            on='case_barcode',
            how='left'
        )
    )
    meta = meta[meta['cell_barcode'].isin(cells_keep)].copy()
    meta = meta.rename(columns={'cell_barcode': 'cell'})
    meta.index = meta['cell']
    return meta


def matrix_to_anndata(mc_path, cov_path, qc_file, clinical_file, max_features, min_cells_per_feature):
    qc_tmp = pd.read_csv(qc_file, sep='\t')
    qc_cells = set(clean(qc_tmp['cell_barcode']))

    mc_df = fread_gz(mc_path)
    cov_df = fread_gz(cov_path)

    mc = orient_table(mc_df, qc_cells)
    cov = orient_table(cov_df, qc_cells)

    mc.index = [str(x).replace('"', '').strip() for x in mc.index]
    cov.index = [str(x).replace('"', '').strip() for x in cov.index]
    mc.columns = [str(x).replace('"', '').strip() for x in mc.columns]
    cov.columns = [str(x).replace('"', '').strip() for x in cov.columns]

    common_cells = mc.index.intersection(cov.index)
    common_features = mc.columns.intersection(cov.columns)
    mc = mc.loc[common_cells, common_features].copy()
    cov = cov.loc[common_cells, common_features].copy()

    meta = build_meta(qc_file, clinical_file, list(common_cells))
    common_cells = [c for c in common_cells if c in meta.index]
    mc = mc.loc[common_cells, :].copy()
    cov = cov.loc[common_cells, :].copy()
    meta = meta.loc[common_cells, :].copy()

    mc = mc.fillna(0)
    cov = cov.fillna(0)

    covered_cells = (cov.to_numpy() > 0).sum(axis=0)
    keep = covered_cells >= min_cells_per_feature
    mc = mc.loc[:, keep]
    cov = cov.loc[:, keep]
    covered_cells = covered_cells[keep]

    if mc.shape[1] > max_features:
        order = np.argsort(-covered_cells)[:max_features]
        mc = mc.iloc[:, order]
        cov = cov.iloc[:, order]

    frac = np.divide(mc.to_numpy(), cov.to_numpy(), out=np.zeros_like(mc.to_numpy(), dtype=float), where=cov.to_numpy() > 0)

    import anndata as ad
    adata = ad.AnnData(X=sparse.csr_matrix(frac), obs=meta.copy(), var=pd.DataFrame(index=mc.columns.copy()))
    adata.layers['mc'] = sparse.csr_matrix(mc.to_numpy(dtype=np.float32))
    adata.layers['cov'] = sparse.csr_matrix(cov.to_numpy(dtype=np.float32))
    adata.obs_names = meta.index.astype(str)
    adata.var_names = mc.columns.astype(str)
    return adata


def write_group_files(obs_df, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    for cl in sorted(obs_df['methylvi_cluster'].dropna().astype(int).unique()):
        grp = pd.DataFrame({
            'cell': obs_df.index.astype(str),
            'group': np.where(obs_df['methylvi_cluster'].astype(int) == cl, f'cluster_{cl}', 'rest')
        })
        grp.to_csv(outdir / f'cluster_{cl}_vs_rest.csv', index=False)


def main():
    args = parse_args()
    np.random.seed(args.seed)

    outdir = Path(args.results) / args.out_subdir
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'cell_groups').mkdir(exist_ok=True)

    matrix_base = Path(args.results) / args.matrix_dir
    mc_path = matrix_base / 'methylated_sites.csv.gz'
    cov_path = matrix_base / 'total_sites.csv.gz'

    adata = matrix_to_anndata(str(mc_path), str(cov_path), args.qc, args.clinical, args.max_features, args.min_cells_per_feature)

    import mudata
    import scanpy as sc
    import scvi
    import torch
    from scvi.external import METHYLVI
    import matplotlib.pyplot as plt

    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision('high')

    mdata = mudata.MuData({'mCG': adata})
    METHYLVI.setup_mudata(
        mdata,
        mc_layer='mc',
        cov_layer='cov',
        methylation_contexts=['mCG'],
        batch_key=args.batch_key,
        modalities={'batch_key': 'mCG'}
    )

    model = METHYLVI(mdata, n_latent=args.n_latent)
    model.train(
        max_epochs=args.max_epochs,
        early_stopping=True,
        batch_size=args.batch_size,
        accelerator=args.accelerator,
        devices=args.devices,
    )

    latent = model.get_latent_representation()
    adata = mdata['mCG']
    adata.obsm['X_methylvi'] = latent

    sc.pp.neighbors(adata, use_rep='X_methylvi', n_neighbors=args.n_neighbors)
    sc.tl.umap(adata, min_dist=args.min_dist)
    sc.tl.leiden(adata, resolution=args.resolution, key_added='methylvi_cluster')

    latent_df = pd.DataFrame(latent, index=adata.obs_names)
    latent_df.insert(0, 'cell', latent_df.index)
    latent_df.to_csv(outdir / 'methylvi_latent.tsv', sep='\t', index=False)

    umap_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP1', 'UMAP2'])
    umap_df.insert(0, 'cell', umap_df.index)
    umap_df.to_csv(outdir / 'methylvi_umap.tsv', sep='\t', index=False)

    obs = adata.obs.copy()
    obs['cell'] = obs.index.astype(str)
    obs['methylvi_cluster'] = obs['methylvi_cluster'].astype(str)
    obs.to_csv(outdir / 'cell_cluster_assignments.tsv', sep='\t', index=False)

    write_group_files(obs.set_index('cell'), outdir / 'cell_groups')

    cluster_sizes = obs.groupby('methylvi_cluster').size().reset_index(name='n_cells')
    cluster_sizes.to_csv(outdir / 'cluster_sizes.tsv', sep='\t', index=False)

    cluster_by_case = obs.groupby(['case_barcode', 'methylvi_cluster']).size().reset_index(name='n')
    cluster_by_case.to_csv(outdir / 'cluster_by_case.tsv', sep='\t', index=False)
    cluster_by_case_fraction = cluster_by_case.copy()
    cluster_by_case_fraction['frac'] = cluster_by_case_fraction['n'] / cluster_by_case_fraction.groupby('case_barcode')['n'].transform('sum')
    cluster_by_case_fraction.to_csv(outdir / 'cluster_by_case_fraction.tsv', sep='\t', index=False)

    sc.pl.umap(adata, color='methylvi_cluster', show=False)
    plt.savefig(outdir / 'UMAP_clusters.pdf', bbox_inches='tight')
    plt.close()
    sc.pl.umap(adata, color='case_barcode', show=False)
    plt.savefig(outdir / 'UMAP_cases.pdf', bbox_inches='tight')
    plt.close()
    sc.pl.umap(adata, color='idh_status', show=False)
    plt.savefig(outdir / 'UMAP_IDH.pdf', bbox_inches='tight')
    plt.close()

    model_dir = outdir / 'model'
    model.save(model_dir, overwrite=True, save_anndata=False)
    mdata.write(outdir / 'methylvi_input_and_latent.h5mu')

    with open(outdir / 'run_metadata.json', 'w') as fh:
        json.dump({
            'mode': args.mode,
            'matrix_dir': args.matrix_dir,
            'n_latent': args.n_latent,
            'max_epochs': args.max_epochs,
            'batch_key': args.batch_key,
            'n_neighbors': args.n_neighbors,
            'min_dist': args.min_dist,
            'resolution': args.resolution,
            'max_features': args.max_features,
            'min_cells_per_feature': args.min_cells_per_feature,
        }, fh, indent=2)

    print(f'MethylVI pipeline finished: {outdir}')


if __name__ == '__main__':
    main()
