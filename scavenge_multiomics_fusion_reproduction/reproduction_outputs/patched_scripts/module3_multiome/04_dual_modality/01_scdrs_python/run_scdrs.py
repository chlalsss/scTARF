"""
run_scdrs.py — thin wrapper around official SCDRS 1.0 CLI.

不做任何 SCDRS 自定义逻辑, 只是把两条 CLI 命令串起来, 一行跑完.

调用的两条 SCDRS 官方 CLI:
    1. scdrs compute-score
        -> output/<trait>.score.gz
        -> output/<trait>.full_score.gz   (含 ctrl_norm_score_*, 给 perform-downstream 用)
    2. scdrs perform-downstream --group-analysis cell_type
        -> output/<trait>.scdrs_group.cell_type   (cell-type assoc + hetero MC tests)

用法 (PowerShell):
    cd D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/patched_scripts/module3_multiome\\04_dual_modality\\01_scdrs_python
    python run_scdrs.py

或者两条命令手动跑 (完全等价):
    scdrs compute-score      --h5ad-file ... --gs-file ... --out-folder output ...
    scdrs perform-downstream --h5ad-file ... --score-file output/@.full_score.gz ...
"""
import os
import subprocess
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
H5AD = "D:/RStudio/work/scavenge_multiomics_fusion_reproduction/reproduction_outputs/BMMC/04_dual_modality/BMMC_rna.h5ad"
GS_FILE = os.path.join(THIS_DIR, "bmmc_3traits.gs")
OUT_DIR = os.path.join(THIS_DIR, "output")
os.makedirs(OUT_DIR, exist_ok=True)

SCDRS = os.path.join(os.path.dirname(sys.executable), "Scripts", "scdrs")
if not os.path.exists(SCDRS):
    sys.exit(f"FAIL: scdrs CLI not found at {SCDRS}")

# 1. compute-score  (~45-90 min, single thread numpy)
EXPECTED_TRAITS = [
    "UKB_460K.blood_LYMPHOCYTE_COUNT",
    "UKB_460K.blood_MONOCYTE_COUNT",
    "UKB_460K.blood_PLATELET_COUNT",
]
expected_full = [os.path.join(OUT_DIR, f"{t}.full_score.gz") for t in EXPECTED_TRAITS]

print("=" * 78)
print("[1/2] scdrs compute-score")
print("=" * 78)
if all(os.path.exists(f) for f in expected_full):
    print("all .full_score.gz already exist, skip compute-score")
    for f in expected_full:
        print(f"  {f}  ({os.path.getsize(f)/1e6:.1f} MB)")
else:
    cmd1 = [
        sys.executable, SCDRS, "compute-score",
        "--h5ad-file", H5AD,
        "--h5ad-species", "human",
        "--gs-file", GS_FILE,
        "--gs-species", "human",
        "--out-folder", OUT_DIR,
        "--n-ctrl", "1000",
    ]
    print(" ".join(cmd1) + "\n")
    subprocess.run(cmd1, check=True)

# 2. perform-downstream (cell-type group analysis, ~3-5 min)
cmd2 = [
    sys.executable, SCDRS, "perform-downstream",
    "--h5ad-file", H5AD,
    "--score-file", os.path.join(OUT_DIR, "@.full_score.gz"),
    "--out-folder", OUT_DIR,
    "--group-analysis", "cell_type",
]
print("\n" + "=" * 78)
print("[2/2] scdrs perform-downstream  (cell-type group analysis)")
print("=" * 78)
print(" ".join(cmd2) + "\n")
subprocess.run(cmd2, check=True)

print(f"\nrun_scdrs.py OK. Outputs in {OUT_DIR}/")
print("\nSCDRS 标准输出:")
print(f"  {OUT_DIR}/<trait>.score.gz                cell-level scores (主结果)")
print(f"  {OUT_DIR}/<trait>.full_score.gz           cell-level + ctrl 列 (下游用)")
print(f"  {OUT_DIR}/<trait>.scdrs_group.cell_type   cell-type assoc + hetero MC")
