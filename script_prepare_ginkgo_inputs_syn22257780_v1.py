#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path

import pandas as pd


def clean_barcode(x: str) -> str:
    x = str(x).strip().replace('.', '-')
    x = re.sub(r'[^A-Za-z0-9_\-]', '_', x)
    return x


def pick_col(cols, candidates):
    cmap = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in cmap:
            return cmap[cand.lower()]
    return None


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', compression='infer', low_memory=False)


def main():
    ap = argparse.ArgumentParser(description='Prepare BAM manifest filtered by scRRBS QC and mode.')
    ap.add_argument('--bam-manifest', required=True)
    ap.add_argument('--qc', required=True)
    ap.add_argument('--clinical', required=True)
    ap.add_argument('--mode', required=True, choices=['all', 'IDHwt', 'IDHmut'])
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--min-unique-cpg', type=int, default=40000)
    ap.add_argument('--min-bs', type=float, default=95)
    ap.add_argument('--require-tumor', type=int, default=1)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bam = read_table(args.bam_manifest)
    qc = read_table(args.qc)
    clin = read_table(args.clinical)

    bam_cell = pick_col(bam.columns, ['cell_barcode', 'cell', 'barcode', 'cell_id'])
    bam_path = pick_col(bam.columns, ['bam_path', 'bam', 'path'])
    qc_cell = pick_col(qc.columns, ['cell_barcode', 'cell', 'barcode', 'cell_id'])
    qc_case = pick_col(qc.columns, ['case_barcode', 'case', 'sample'])
    qc_cpg = pick_col(qc.columns, ['cpg_unique', 'unique_cpg', 'unique_cpgs'])
    qc_bs = pick_col(qc.columns, ['bisulfite_conversion_rate', 'bisulfite_conversion'])
    qc_tumor = pick_col(qc.columns, ['tumor_status', 'is_tumor'])
    clin_case = pick_col(clin.columns, ['case_barcode', 'case', 'sample'])
    clin_idh = pick_col(clin.columns, ['idh_codel_subtype', 'idh_status', 'idh'])

    if None in [bam_cell, bam_path, qc_cell, qc_case, qc_cpg, qc_bs, qc_tumor, clin_case, clin_idh]:
        raise SystemExit('Could not identify required columns in BAM manifest / QC / clinical tables.')

    bam[bam_cell] = bam[bam_cell].map(clean_barcode)
    qc[qc_cell] = qc[qc_cell].map(clean_barcode)
    qc[qc_case] = qc[qc_case].map(clean_barcode)
    clin[clin_case] = clin[clin_case].map(clean_barcode)
    clin['idh_status'] = clin[clin_idh].astype(str).str.contains('IDHwt', case=False, na=False).map({True: 'IDHwt', False: 'IDHmut'})

    qc[qc_cpg] = pd.to_numeric(qc[qc_cpg], errors='coerce')
    qc[qc_bs] = pd.to_numeric(qc[qc_bs], errors='coerce')
    qc[qc_tumor] = pd.to_numeric(qc[qc_tumor], errors='coerce')

    keep = qc.loc[
        (qc[qc_cpg] > args.min_unique_cpg) &
        (qc[qc_bs] > args.min_bs) &
        (qc[qc_tumor] == args.require_tumor)
    ].copy()
    keep = keep.merge(clin[[clin_case, 'idh_status']].drop_duplicates(), left_on=qc_case, right_on=clin_case, how='left')

    if args.mode == 'IDHwt':
        keep = keep.loc[keep['idh_status'] == 'IDHwt'].copy()
    elif args.mode == 'IDHmut':
        keep = keep.loc[keep['idh_status'] == 'IDHmut'].copy()

    keep = keep.rename(columns={qc_cell: 'cell_barcode', qc_case: 'case_barcode'})
    keep = keep[['cell_barcode', 'case_barcode', 'idh_status']].drop_duplicates()

    bam_keep = bam.merge(keep, left_on=bam_cell, right_on='cell_barcode', how='inner').copy()
    bam_keep = bam_keep[[bam_cell, bam_path, 'case_barcode', 'idh_status']].copy()
    bam_keep.columns = ['cell_barcode', 'bam_path', 'case_barcode', 'idh_status']
    bam_keep['bam_exists'] = bam_keep['bam_path'].map(os.path.exists)

    bam_keep.to_csv(outdir / 'bam_manifest_filtered.tsv', sep='\t', index=False)
    with open(outdir / 'cells_to_keep.txt', 'w', encoding='utf-8') as fh:
        for cell in sorted(bam_keep['cell_barcode'].unique()):
            fh.write(f'{cell}\n')

    print(f'Filtered BAM manifest written to: {outdir / "bam_manifest_filtered.tsv"}')
    print(f'Cells retained: {bam_keep.shape[0]}')


if __name__ == '__main__':
    main()
