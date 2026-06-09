"""Create a read-level RNA dropout h5ad for scDRS robustness.

This helper only perturbs the input count matrix. It does not implement scDRS.
The R driver calls the official scDRS CLI on the output h5ad.
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
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--p-drop", type=float, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--summary-csv", required=True)
    return parser.parse_args()


def dropout_matrix(x, p_keep: float, seed: int):
    rng = np.random.default_rng(seed)
    if sparse.issparse(x):
        out = x.copy().tocsr()
        before_nnz = int(out.nnz)
        before_sum = float(out.data.sum())
        out.data = rng.binomial(out.data.astype(np.int64), p_keep).astype(np.float32)
        out.eliminate_zeros()
        return out, before_nnz, int(out.nnz), before_sum, float(out.data.sum())

    out = np.asarray(x).copy()
    before_nnz = int(np.count_nonzero(out))
    before_sum = float(out.sum())
    out = rng.binomial(out.astype(np.int64), p_keep).astype(np.float32)
    return out, before_nnz, int(np.count_nonzero(out)), before_sum, float(out.sum())


def main() -> int:
    args = parse_args()
    if not (0 <= args.p_drop < 1):
        raise ValueError("--p-drop must be in [0, 1).")
    p_keep = 1.0 - args.p_drop

    output = Path(args.output_h5ad)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary = Path(args.summary_csv)
    summary.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad)
    x_new, before_nnz, after_nnz, before_sum, after_sum = dropout_matrix(
        adata.X, p_keep=p_keep, seed=args.seed
    )
    adata.X = x_new
    adata.write_h5ad(output, compression="gzip")

    pd.DataFrame(
        [
            {
                "input_h5ad": str(args.input_h5ad),
                "output_h5ad": str(output),
                "p_drop": args.p_drop,
                "p_keep": p_keep,
                "seed": args.seed,
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "nnz_before": before_nnz,
                "nnz_after": after_nnz,
                "sum_before": before_sum,
                "sum_after": after_sum,
            }
        ]
    ).to_csv(summary, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
