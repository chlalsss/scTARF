"""Run official scDRS CLI for NeurIPS 2021 BMMC CITE-seq RNA.

Inputs:
  input/cite_rna_raw_counts_gex_only.h5ad
  cite_3traits.gs

Outputs:
  output/<trait>.score.gz
  output/<trait>.full_score.gz
  output/<trait>.scdrs_group.cell_type
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


ROOT = Path("D:/RStudio/work/module4_traitaware_pseudopair")
THIS_DIR = ROOT / "01_data_prep_external_rna" / "03_cite_scdrs_python"
H5AD = THIS_DIR / "input" / "cite_rna_raw_counts_gex_only.h5ad"
GS_FILE = THIS_DIR / "cite_3traits.gs"
OUT_DIR = THIS_DIR / "output"
PYTHON = Path("D:/python/python.exe")
SCDRS = Path("D:/python/Scripts/scdrs")

TRAITS = [
    "UKB_460K.blood_MONOCYTE_COUNT",
    "UKB_460K.blood_PLATELET_COUNT",
    "UKB_460K.blood_LYMPHOCYTE_COUNT",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-ctrl", type=int, default=1000)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--force-compute", action="store_true")
    parser.add_argument("--force-downstream", action="store_true")
    return parser.parse_args()


def run(cmd: list[str]) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def require_files() -> None:
    missing = [path for path in [H5AD, GS_FILE, PYTHON, SCDRS] if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required file(s): " + ", ".join(map(str, missing)))


def main() -> int:
    args = parse_args()
    os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / ".numba_cache"))
    os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    require_files()

    expected_full = [OUT_DIR / f"{trait}.full_score.gz" for trait in TRAITS]
    expected_group = [OUT_DIR / f"{trait}.scdrs_group.cell_type" for trait in TRAITS]

    print("CITE scDRS run configuration", flush=True)
    print(f"  h5ad:   {H5AD}", flush=True)
    print(f"  gs:     {GS_FILE}", flush=True)
    print(f"  output: {OUT_DIR}", flush=True)
    print(f"  n_ctrl: {args.n_ctrl}", flush=True)

    if args.check_only:
        print("check-only OK. scDRS was not run.", flush=True)
        return 0

    print("\n" + "=" * 78, flush=True)
    print("[1/2] scDRS compute-score", flush=True)
    print("=" * 78, flush=True)
    if (not args.force_compute) and all(path.exists() for path in expected_full):
        print("all .full_score.gz files exist, skip compute-score", flush=True)
        for path in expected_full:
            print(f"  {path} ({path.stat().st_size / 1e6:.1f} MB)", flush=True)
    else:
        cmd = [
            str(PYTHON),
            str(SCDRS),
            "compute-score",
            "--h5ad-file",
            str(H5AD),
            "--h5ad-species",
            "human",
            "--gs-file",
            str(GS_FILE),
            "--gs-species",
            "human",
            "--out-folder",
            str(OUT_DIR),
            "--n-ctrl",
            str(args.n_ctrl),
        ]
        run(cmd)

    print("\n" + "=" * 78, flush=True)
    print("[2/2] scDRS perform-downstream by cell_type", flush=True)
    print("=" * 78, flush=True)
    if (not args.force_downstream) and all(path.exists() for path in expected_group):
        print("all .scdrs_group.cell_type files exist, skip perform-downstream", flush=True)
        for path in expected_group:
            print(f"  {path} ({path.stat().st_size / 1e3:.1f} KB)", flush=True)
    else:
        cmd = [
            str(PYTHON),
            str(SCDRS),
            "perform-downstream",
            "--h5ad-file",
            str(H5AD),
            "--score-file",
            str(OUT_DIR / "@.full_score.gz"),
            "--out-folder",
            str(OUT_DIR),
            "--group-analysis",
            "cell_type",
        ]
        run(cmd)

    print("\nrun_scdrs_cite.py OK", flush=True)
    print(f"Outputs in {OUT_DIR}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
