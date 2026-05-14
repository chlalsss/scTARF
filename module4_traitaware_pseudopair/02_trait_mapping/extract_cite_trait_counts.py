"""Extract selected CITE-seq GEX genes from h5ad for R-side mapping.

The R alignment script only needs a small trait-gene subset from the CITE
reference. This helper avoids loading/converting the full h5ad in R.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--out-mtx", required=True)
    parser.add_argument("--out-genes", required=True)
    parser.add_argument("--out-cells", required=True)
    parser.add_argument("--out-summary", required=True)
    return parser.parse_args()


def read_lines(path: str | Path) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        item = line.strip()
        if not item or item in seen:
            continue
        out.append(item)
        seen.add(item)
    return out


def main() -> None:
    args = parse_args()
    out_mtx = Path(args.out_mtx)
    out_mtx.parent.mkdir(parents=True, exist_ok=True)

    requested = read_lines(args.genes)
    if not requested:
        raise ValueError("No genes requested")

    print(f"Reading h5ad: {args.h5ad}", flush=True)
    adata = ad.read_h5ad(args.h5ad)
    var_names = pd.Index(adata.var_names.astype(str))

    present = [g for g in requested if g in var_names]
    missing = [g for g in requested if g not in var_names]
    if not present:
        raise ValueError("None of the requested genes are present in h5ad")

    print(
        f"Requested genes: {len(requested)}; present: {len(present)}; missing: {len(missing)}",
        flush=True,
    )

    sub = adata[:, present]
    x = sub.X
    if sparse.issparse(x):
        mat = x.T.tocoo()
    else:
        mat = sparse.coo_matrix(np.asarray(x).T)

    print(f"Writing MatrixMarket gene x cell matrix: {out_mtx}", flush=True)
    mmwrite(out_mtx, mat)

    Path(args.out_genes).write_text("\n".join(present) + "\n", encoding="utf-8")
    Path(args.out_cells).write_text(
        "\n".join(map(str, sub.obs_names)) + "\n", encoding="utf-8"
    )

    summary = pd.DataFrame(
        [
            {
                "h5ad": str(Path(args.h5ad).resolve()),
                "n_cells": sub.n_obs,
                "n_requested_genes": len(requested),
                "n_present_genes": len(present),
                "n_missing_genes": len(missing),
                "matrix_nnz": int(mat.nnz),
                "matrix_density": float(mat.nnz / max(mat.shape[0] * mat.shape[1], 1)),
            }
        ]
    )
    summary.to_csv(args.out_summary, index=False)

    missing_path = Path(args.out_summary).with_name("cite_trait_counts_missing_genes.txt")
    missing_path.write_text("\n".join(missing) + ("\n" if missing else ""), encoding="utf-8")

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
