"""Export ATAC LSI npy + cell CSV to a compact CSV for R.

This is a small adapter helper. It does not change the embedding; it only
attaches cell barcodes to the rows so the R adapter can align cells safely.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--npy", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--out-csv", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    x = np.load(args.npy)
    cells = pd.read_csv(args.cells)
    if "cell" not in cells.columns:
        raise ValueError(f"Expected a 'cell' column in {args.cells}")
    if x.shape[0] != cells.shape[0]:
        raise ValueError(
            f"Row mismatch: npy has {x.shape[0]} rows, cells CSV has {cells.shape[0]}"
        )

    df = pd.DataFrame(x, columns=[f"LSI_{i + 1}" for i in range(x.shape[1])])
    df.insert(0, "barcode", cells["cell"].astype(str).to_numpy())
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv} with shape {df.shape[0]} x {df.shape[1]}", flush=True)


if __name__ == "__main__":
    main()
