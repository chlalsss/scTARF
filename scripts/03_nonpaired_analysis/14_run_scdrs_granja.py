"""Run official scDRS CLI for Granja RNA.

Outputs:
  output/<trait>.score.gz
  output/<trait>.full_score.gz
  output/<trait>.scdrs_group.cell_type
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(os.environ.get("SCTARF_PROJECT_ROOT", Path.cwd()))
OUTPUT_ROOT = Path(os.environ.get("SCTARF_OUTPUT_ROOT", PROJECT_ROOT / "outputs"))
ROOT = Path(os.environ.get("SCTARF_NONPAIRED_ROOT", OUTPUT_ROOT / "nonpaired"))
THIS_DIR = ROOT / "01_data_prep" / "03_rna_scdrs_python"
H5AD = THIS_DIR / "input" / "granja_rna_raw_counts.h5ad"
GS_FILE = THIS_DIR / "granja_3traits.gs"
OUT_DIR = THIS_DIR / "output"
PYTHON_NAME = os.environ.get("SCTARF_PYTHON", sys.executable)
PYTHON = Path(shutil.which(PYTHON_NAME) or PYTHON_NAME)
SCDRS = os.environ.get("SCTARF_SCDRS") or shutil.which("scdrs")
if SCDRS is None:
    SCDRS = str(Path(sys.executable).resolve().parent / "Scripts" / "scdrs")
SCDRS = Path(shutil.which(SCDRS) or SCDRS)

TRAITS = [
    "UKB_460K.blood_MONOCYTE_COUNT",
    "UKB_460K.blood_PLATELET_COUNT",
    "UKB_460K.blood_LYMPHOCYTE_COUNT",
]


def run(cmd: list[str]) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> int:
    os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / ".numba_cache"))
    os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for path in [H5AD, GS_FILE, PYTHON, SCDRS]:
        if not path.exists():
            raise FileNotFoundError(path)

    expected_full = [OUT_DIR / f"{trait}.full_score.gz" for trait in TRAITS]

    print("=" * 78, flush=True)
    print("[1/2] scDRS compute-score", flush=True)
    print("=" * 78, flush=True)
    if all(path.exists() for path in expected_full):
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
            "1000",
        ]
        run(cmd)

    print("\n" + "=" * 78, flush=True)
    print("[2/2] scDRS perform-downstream by cell_type", flush=True)
    print("=" * 78, flush=True)
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

    print("\nrun_scdrs_granja.py OK", flush=True)
    print(f"Outputs in {OUT_DIR}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
