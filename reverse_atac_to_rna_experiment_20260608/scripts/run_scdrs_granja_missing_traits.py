"""Run official scDRS CLI for Granja RNA missing blood traits.

This is a thin wrapper only. It does not implement scDRS scoring logic.

Inputs:
  Granja RNA h5ad:
    module4_traitaware_pseudopair/01_data_prep/03_rna_scdrs_python/input/granja_rna_raw_counts.h5ad
  Existing official MAGMA gene-set subset:
    module4_traitaware_pseudopair/01_data_prep_external_rna/03_cite_scdrs_python/cite_trait_screening.gs

Outputs:
  module4_traitaware_pseudopair/01_data_prep/03_rna_scdrs_python/output_missing_traits/
    <trait>.score.gz
    <trait>.full_score.gz
    <trait>.scdrs_group.cell_type
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


WORK_ROOT = Path("D:/RStudio/work")
MODULE4_ROOT = WORK_ROOT / "module4_traitaware_pseudopair"
GRANJA_SCDRS_DIR = MODULE4_ROOT / "01_data_prep" / "03_rna_scdrs_python"
H5AD = GRANJA_SCDRS_DIR / "input" / "granja_rna_raw_counts.h5ad"
GS_FILE = (
    MODULE4_ROOT
    / "01_data_prep_external_rna"
    / "03_cite_scdrs_python"
    / "cite_trait_screening.gs"
)
OUT_DIR = GRANJA_SCDRS_DIR / "output_missing_traits"
PYTHON = Path("D:/python/python.exe")
SCDRS = Path("D:/python/Scripts/scdrs")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-ctrl", type=int, default=1000)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--force-compute", action="store_true")
    parser.add_argument("--force-downstream", action="store_true")
    parser.add_argument("--skip-downstream", action="store_true")
    return parser.parse_args()


def read_traits(gs_file: Path) -> list[str]:
    with gs_file.open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if len(header) < 2 or header[0] != "TRAIT":
            raise ValueError(f"Unexpected .gs header in {gs_file}: {header}")
        traits = []
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            traits.append(line.split("\t", 1)[0])
    if not traits:
        raise ValueError(f"No traits found in {gs_file}")
    return traits


def require_files() -> None:
    missing = [path for path in [H5AD, GS_FILE, PYTHON, SCDRS] if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required file(s): " + ", ".join(map(str, missing)))


def run(cmd: list[str]) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> int:
    args = parse_args()
    os.environ.setdefault("NUMBA_CACHE_DIR", str(MODULE4_ROOT / ".numba_cache"))
    os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    require_files()
    traits = read_traits(GS_FILE)

    expected_full = [OUT_DIR / f"{trait}.full_score.gz" for trait in traits]
    expected_group = [OUT_DIR / f"{trait}.scdrs_group.cell_type" for trait in traits]

    print("Granja missing-traits official scDRS run configuration", flush=True)
    print(f"  h5ad:   {H5AD}", flush=True)
    print(f"  gs:     {GS_FILE}", flush=True)
    print(f"  output: {OUT_DIR}", flush=True)
    print(f"  n_ctrl: {args.n_ctrl}", flush=True)
    print("  traits:", flush=True)
    for trait in traits:
        print(f"    {trait}", flush=True)

    if args.check_only:
        print("check-only OK. scDRS was not run.", flush=True)
        return 0

    print("\n" + "=" * 78, flush=True)
    print("[1/2] official scDRS compute-score", flush=True)
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

    if args.skip_downstream:
        print("\nskip-downstream set. Post hoc cell_type group analysis was not run.", flush=True)
        return 0

    print("\n" + "=" * 78, flush=True)
    print("[2/2] official scDRS perform-downstream by cell_type", flush=True)
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

    print("\nrun_scdrs_granja_missing_traits.py OK", flush=True)
    print(f"Outputs in {OUT_DIR}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
