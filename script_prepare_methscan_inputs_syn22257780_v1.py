#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


QC_CANDIDATES = {
    "cell": ["cell_barcode", "cell", "barcode", "cell_id"],
    "case": ["case_barcode", "case", "sample", "sample_id", "patient", "patient_id"],
    "cpg_unique": ["cpg_unique", "unique_cpg", "unique_cpgs", "n_unique_cpg", "num_unique_cpgs"],
    "bs": ["bisulfite_conversion_rate", "bisulfite_conversion", "bs_conversion_rate", "conversion_rate"],
    "tumor": ["tumor_status", "is_tumor", "tumour_status"],
}

CLINICAL_CANDIDATES = {
    "case": ["case_barcode", "case", "sample", "sample_id", "patient", "patient_id"],
    "idh": ["idh_codel_subtype", "idh_status", "idh", "IDH_status"],
}

COV_CANDIDATES = {
    "cell": ["cell_barcode", "cell", "barcode", "cell_id"],
    "chrom": ["chromosome", "chrom", "chr", "seqnames"],
    "start": ["start", "position", "pos", "bp"],
    "end": ["end", "stop"],
    "meth": ["count_methylated", "methylated_count", "meth_count", "n_meth", "meth"],
    "unmeth": ["count_unmethylated", "unmethylated_count", "unmeth_count", "n_unmeth", "unmeth"],
    "total": ["coverage", "total_count", "n_total", "total", "count_total"],
    "pct": ["methylation_percentage", "pct_methylation", "percent_methylation", "beta_value", "beta"],
}

REF_CANDIDATES = {
    "chrom": ["chromosome", "chrom", "chr", "seqnames"],
    "start": ["start", "tss_start", "promoter_start", "gene_start"],
    "end": ["end", "tss_end", "promoter_end", "gene_end"],
    "name": ["gene_name", "gene_symbol", "promoter_id", "gene_id", "feature_id", "id", "name"],
}


def clean_barcode(x: str) -> str:
    x = str(x).strip()
    x = x.replace(".", "-")
    x = re.sub(r"[^A-Za-z0-9_\-]", "_", x)
    return x



def pick_col(columns: Iterable[str], candidates: List[str], explicit: Optional[str] = None) -> Optional[str]:
    cols = list(columns)
    colmap = {c.lower(): c for c in cols}
    if explicit:
        if explicit not in cols:
            raise ValueError(f"Requested column '{explicit}' not found. Available: {cols}")
        return explicit
    for cand in candidates:
        if cand.lower() in colmap:
            return colmap[cand.lower()]
    return None



def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)



def prepare_qc_whitelist(qc_file: str, clinical_file: str, mode: str, outdir: str,
                         min_unique_cpg: int, min_bs: float, require_tumor: int) -> pd.DataFrame:
    qc = read_table(qc_file)
    clin = read_table(clinical_file)

    qc_cell = pick_col(qc.columns, QC_CANDIDATES["cell"])
    qc_case = pick_col(qc.columns, QC_CANDIDATES["case"])
    qc_cpg = pick_col(qc.columns, QC_CANDIDATES["cpg_unique"])
    qc_bs = pick_col(qc.columns, QC_CANDIDATES["bs"])
    qc_tumor = pick_col(qc.columns, QC_CANDIDATES["tumor"])

    clin_case = pick_col(clin.columns, CLINICAL_CANDIDATES["case"])
    clin_idh = pick_col(clin.columns, CLINICAL_CANDIDATES["idh"])

    required = [qc_cell, qc_case, qc_cpg, qc_bs, qc_tumor, clin_case, clin_idh]
    if any(x is None for x in required):
        raise ValueError("Could not identify required QC/clinical columns. Please inspect your input tables.")

    qc = qc.copy()
    clin = clin.copy()
    qc[qc_cell] = qc[qc_cell].map(clean_barcode)
    qc[qc_case] = qc[qc_case].map(clean_barcode)
    clin[clin_case] = clin[clin_case].map(clean_barcode)

    clin["idh_status"] = clin[clin_idh].astype(str).str.contains("IDHwt", case=False, na=False).map({True: "IDHwt", False: "IDHmut"})

    merged = qc.merge(clin[[clin_case, clin_idh, "idh_status"]].drop_duplicates(), left_on=qc_case, right_on=clin_case, how="left")

    merged[qc_cpg] = pd.to_numeric(merged[qc_cpg], errors="coerce")
    merged[qc_bs] = pd.to_numeric(merged[qc_bs], errors="coerce")
    merged[qc_tumor] = pd.to_numeric(merged[qc_tumor], errors="coerce")

    keep = merged.loc[
        (merged[qc_cpg] > min_unique_cpg) &
        (merged[qc_bs] > min_bs) &
        (merged[qc_tumor] == require_tumor)
    ].copy()

    if mode == "IDHwt":
        keep = keep.loc[keep["idh_status"] == "IDHwt"].copy()
    elif mode == "IDHmut":
        keep = keep.loc[keep["idh_status"] == "IDHmut"].copy()
    elif mode != "all":
        raise ValueError("mode must be one of: all, IDHwt, IDHmut")

    keep = keep.rename(columns={qc_cell: "cell_barcode", qc_case: "case_barcode"})
    keep = keep.drop_duplicates(subset=["cell_barcode"]).copy()

    Path(outdir).mkdir(parents=True, exist_ok=True)
    keep.to_csv(os.path.join(outdir, "qc_keep.tsv"), sep="\t", index=False)
    with open(os.path.join(outdir, "qc_keep_cell_names.txt"), "w", encoding="utf-8") as fh:
        for cell in sorted(keep["cell_barcode"].unique()):
            fh.write(f"{cell}\n")
    return keep



def stream_split_bismark_coverage(coverage_file: str, keep_cells: set, outdir: str,
                                  cell_col: Optional[str] = None,
                                  chrom_col: Optional[str] = None,
                                  start_col: Optional[str] = None,
                                  end_col: Optional[str] = None,
                                  meth_col: Optional[str] = None,
                                  unmeth_col: Optional[str] = None,
                                  total_col: Optional[str] = None,
                                  pct_col: Optional[str] = None,
                                  chunksize: int = 500000) -> None:
    cov_dir = Path(outdir) / "cov_by_cell"
    cov_dir.mkdir(parents=True, exist_ok=True)

    handles: Dict[str, gzip.GzipFile] = {}
    cols_resolved: Optional[Dict[str, Optional[str]]] = None

    try:
        for chunk in pd.read_csv(
        coverage_file,
        delim_whitespace=True,
        compression="infer",
        low_memory=False,
        chunksize=chunksize,
    ):
            if cols_resolved is None:
                cols_resolved = {
                    "cell": pick_col(chunk.columns, COV_CANDIDATES["cell"], cell_col),
                    "chrom": pick_col(chunk.columns, COV_CANDIDATES["chrom"], chrom_col),
                    "start": pick_col(chunk.columns, COV_CANDIDATES["start"], start_col),
                    "end": pick_col(chunk.columns, COV_CANDIDATES["end"], end_col),
                    "meth": pick_col(chunk.columns, COV_CANDIDATES["meth"], meth_col),
                    "unmeth": pick_col(chunk.columns, COV_CANDIDATES["unmeth"], unmeth_col),
                    "total": pick_col(chunk.columns, COV_CANDIDATES["total"], total_col),
                    "pct": pick_col(chunk.columns, COV_CANDIDATES["pct"], pct_col),
                }
                req = [cols_resolved["cell"], cols_resolved["chrom"], cols_resolved["start"]]
                if any(x is None for x in req):
                    raise ValueError(f"Could not infer required coverage columns from: {list(chunk.columns)}")
                if cols_resolved["meth"] is None and not (cols_resolved["total"] and cols_resolved["pct"]):
                    raise ValueError("Need methylated+unmethylated columns, or total+percent methylation columns.")

            c = chunk.copy()
            c[cols_resolved["cell"]] = c[cols_resolved["cell"]].map(clean_barcode)
            c = c.loc[c[cols_resolved["cell"]].isin(keep_cells)].copy()
            if c.empty:
                continue

            c[cols_resolved["chrom"]] = c[cols_resolved["chrom"]].astype(str)
            c[cols_resolved["start"]] = pd.to_numeric(c[cols_resolved["start"]], errors="coerce")
            if cols_resolved["end"] is not None:
                c[cols_resolved["end"]] = pd.to_numeric(c[cols_resolved["end"]], errors="coerce")
            else:
                c["end_tmp"] = c[cols_resolved["start"]]
                cols_resolved["end"] = "end_tmp"

            if cols_resolved["meth"] is not None and cols_resolved["unmeth"] is not None:
                c[cols_resolved["meth"]] = pd.to_numeric(c[cols_resolved["meth"]], errors="coerce")
                c[cols_resolved["unmeth"]] = pd.to_numeric(c[cols_resolved["unmeth"]], errors="coerce")
                c["meth_tmp"] = c[cols_resolved["meth"]]
                c["unmeth_tmp"] = c[cols_resolved["unmeth"]]
            else:
                c[cols_resolved["total"]] = pd.to_numeric(c[cols_resolved["total"]], errors="coerce")
                c[cols_resolved["pct"]] = pd.to_numeric(c[cols_resolved["pct"]], errors="coerce")
                c["meth_tmp"] = (c[cols_resolved["total"]] * c[cols_resolved["pct"]] / 100.0).round()
                c["unmeth_tmp"] = c[cols_resolved["total"]] - c["meth_tmp"]

            c = c.dropna(subset=[cols_resolved["chrom"], cols_resolved["start"], cols_resolved["end"], "meth_tmp", "unmeth_tmp"]).copy()
            c["total_tmp"] = c["meth_tmp"] + c["unmeth_tmp"]
            c = c.loc[c["total_tmp"] > 0].copy()
            c["pct_tmp"] = (100.0 * c["meth_tmp"] / c["total_tmp"]).astype(float)

            for cell, grp in c.groupby(cols_resolved["cell"], sort=False):
                if cell not in handles:
                    handles[cell] = gzip.open(cov_dir / f"{cell}.cov.gz", "wt", encoding="utf-8")
                fh = handles[cell]
                for row in grp.itertuples(index=False):
                    rowd = row._asdict()
                    fh.write(
                        f"{rowd[cols_resolved['chrom']]}\t"
                        f"{int(rowd[cols_resolved['start']])}\t"
                        f"{int(rowd[cols_resolved['end']])}\t"
                        f"{float(rowd['pct_tmp']):.6f}\t"
                        f"{int(round(float(rowd['meth_tmp'])))}\t"
                        f"{int(round(float(rowd['unmeth_tmp'])))}\n"
                    )
    finally:
        for fh in handles.values():
            fh.close()



def write_bed_from_reference(ref_file: str, out_bed: str) -> None:
    if not ref_file or not os.path.exists(ref_file):
        return
    df = read_table(ref_file)
    chrom = pick_col(df.columns, REF_CANDIDATES["chrom"])
    start = pick_col(df.columns, REF_CANDIDATES["start"])
    end = pick_col(df.columns, REF_CANDIDATES["end"])
    name = pick_col(df.columns, REF_CANDIDATES["name"])
    if chrom is None or start is None or end is None:
        return
    if name is None:
        df["__name__"] = [f"feature_{i+1}" for i in range(df.shape[0])]
        name = "__name__"
    out = df[[chrom, start, end, name]].copy()
    out.columns = ["chrom", "start", "end", "name"]
    out["chrom"] = out["chrom"].astype(str)
    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"] = pd.to_numeric(out["end"], errors="coerce")
    out = out.dropna(subset=["chrom", "start", "end"]).copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    out = out.sort_values(["chrom", "start", "end", "name"])
    Path(out_bed).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_bed, sep="\t", header=False, index=False)



def main() -> None:
    ap = argparse.ArgumentParser(description="Prepare MethSCAn inputs for syn22257780 scRRBS data.")
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--qc", required=True)
    ap.add_argument("--clinical", required=True)
    ap.add_argument("--ref-promoters", required=False)
    ap.add_argument("--ref-genes", required=False)
    ap.add_argument("--mode", required=True, choices=["all", "IDHwt", "IDHmut"])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-unique-cpg", type=int, default=40000)
    ap.add_argument("--min-bs", type=float, default=95)
    ap.add_argument("--require-tumor", type=int, default=1)
    ap.add_argument("--chunksize", type=int, default=500000)
    ap.add_argument("--coverage-cell-col")
    ap.add_argument("--coverage-chrom-col")
    ap.add_argument("--coverage-start-col")
    ap.add_argument("--coverage-end-col")
    ap.add_argument("--coverage-meth-col")
    ap.add_argument("--coverage-unmeth-col")
    ap.add_argument("--coverage-total-col")
    ap.add_argument("--coverage-pct-col")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    keep = prepare_qc_whitelist(
        qc_file=args.qc,
        clinical_file=args.clinical,
        mode=args.mode,
        outdir=str(outdir),
        min_unique_cpg=args.min_unique_cpg,
        min_bs=args.min_bs,
        require_tumor=args.require_tumor,
    )

    keep_cells = set(keep["cell_barcode"].tolist())
    stream_split_bismark_coverage(
        coverage_file=args.coverage,
        keep_cells=keep_cells,
        outdir=str(outdir),
        cell_col=args.coverage_cell_col,
        chrom_col=args.coverage_chrom_col,
        start_col=args.coverage_start_col,
        end_col=args.coverage_end_col,
        meth_col=args.coverage_meth_col,
        unmeth_col=args.coverage_unmeth_col,
        total_col=args.coverage_total_col,
        pct_col=args.coverage_pct_col,
        chunksize=args.chunksize,
    )

    bed_dir = outdir / "bed"
    bed_dir.mkdir(parents=True, exist_ok=True)
    if args.ref_promoters:
        write_bed_from_reference(args.ref_promoters, str(bed_dir / "promoters.bed"))
    if args.ref_genes:
        write_bed_from_reference(args.ref_genes, str(bed_dir / "genes.bed"))

    print(f"Prepared MethSCAn inputs in: {outdir}")
    print(f"Cells retained: {len(keep_cells)}")


if __name__ == "__main__":
    main()
