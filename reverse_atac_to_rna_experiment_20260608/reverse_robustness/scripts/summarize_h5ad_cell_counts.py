"""Write per-cell total counts from an h5ad file.

Used by the reverse robustness driver to keep metadata library-size columns
consistent with the degraded RNA count matrix.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out-csv", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    if sparse.issparse(adata.X):
        total = np.asarray(adata.X.sum(axis=1)).ravel()
        detected = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    else:
        x = np.asarray(adata.X)
        total = x.sum(axis=1)
        detected = (x > 0).sum(axis=1)

    pd.DataFrame(
        {
            "barcode": adata.obs_names.astype(str),
            "dropout_total_counts": total.astype(float),
            "dropout_detected_genes": detected.astype(int),
        }
    ).to_csv(out_csv, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
