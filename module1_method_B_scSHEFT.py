"""模块 1 · 方法 B：官方 scSHEFT 代码适配脚手架。

结论：CSUBioGroup/scSHEFT 可以作为方法 B 的算法本体使用。
本脚本不重写 scSHEFT，只负责：
1. 检查本项目已有输入和官方仓库接口是否匹配
2. 生成官方代码运行模板
3. 记录还缺哪些可执行环境/中间文件

官方仓库：https://github.com/CSUBioGroup/scSHEFT
官方输入：raw ATAC h5ad、Gene Activity h5ad、RNA h5ad、ATAC LSI npy、
ATAC structural PCA npy。
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


DATA_DIR = Path("D:/SCAVENGE-reproducibility-main/blood_data")
IN_DIR = DATA_DIR / "multiomics_output_v2"
GLUE_INPUT_DIR = IN_DIR / "glue_input"
OUT_DIR = DATA_DIR / "multiomics_module1_pairing" / "method_B_scSHEFT"
SCSHEFT_INPUT_DIR = OUT_DIR / "scsheft_input"
OFFICIAL_FILES = ("scSHEFT.py", "models.py", "utils.py", "dataset.py", "demo.py")


def read_mtx_shape(path: Path) -> dict[str, int]:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(f"无法解析 MatrixMarket 维度行: {path}")
            n_rows, n_cols, nnz = map(int, parts)
            return {"n_rows": n_rows, "n_cols": n_cols, "nnz": nnz}
    raise ValueError(f"未在文件中找到 MatrixMarket 维度: {path}")


def count_lines(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def read_csv_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        return next(reader)


def file_info(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else None,
    }


def check_python_packages() -> dict[str, bool]:
    import importlib.util

    packages = ["numpy", "scipy", "pandas", "scanpy", "anndata", "sklearn", "torch", "annoy"]
    return {pkg: importlib.util.find_spec(pkg) is not None for pkg in packages}


def check_official_repo(repo_dir: Path | None) -> dict[str, Any]:
    if repo_dir is None:
        return {"path": None, "exists": False, "missing_files": list(OFFICIAL_FILES)}

    missing = [name for name in OFFICIAL_FILES if not (repo_dir / name).exists()]
    return {
        "path": str(repo_dir),
        "exists": repo_dir.exists(),
        "missing_files": missing,
        "usable": repo_dir.exists() and not missing,
    }


def build_manifest(repo_dir: Path | None = None) -> dict[str, Any]:
    required = {
        "gene_activity_rds": IN_DIR / "gene_activity_matrix.rds",
        "rna_seurat_rds": IN_DIR / "seurat_rna_normalized.rds",
        "atac_metadata_rds": IN_DIR / "atac_metadata.rds",
        "atac_counts_mtx": GLUE_INPUT_DIR / "atac_counts.mtx",
        "atac_peaks_csv": GLUE_INPUT_DIR / "atac_peaks.csv",
        "atac_barcodes_txt": GLUE_INPUT_DIR / "atac_barcodes.txt",
        "atac_cells_csv": GLUE_INPUT_DIR / "atac_cells.csv",
        "rna_counts_mtx": GLUE_INPUT_DIR / "rna_counts.mtx",
        "rna_genes_csv": GLUE_INPUT_DIR / "rna_genes.csv",
        "rna_barcodes_txt": GLUE_INPUT_DIR / "rna_barcodes.txt",
        "rna_cells_csv": GLUE_INPUT_DIR / "rna_cells.csv",
    }
    exported = {
        "rna_counts_shared_genes_mtx": SCSHEFT_INPUT_DIR / "rna_counts_shared_genes.mtx",
        "gene_activity_counts_shared_genes_mtx": SCSHEFT_INPUT_DIR / "gene_activity_counts_shared_genes.mtx",
        "shared_genes_csv": SCSHEFT_INPUT_DIR / "shared_genes.csv",
        "scsheft_rna_h5ad": SCSHEFT_INPUT_DIR / "blood-GEM.h5ad",
        "scsheft_gene_activity_h5ad": SCSHEFT_INPUT_DIR / "blood-GAM.h5ad",
        "scsheft_raw_atac_h5ad": SCSHEFT_INPUT_DIR / "blood_multiome_atac.h5ad",
        "scsheft_lsi_npy": SCSHEFT_INPUT_DIR / "blood_lsi.npy",
        "scsheft_struct_pca_npy": SCSHEFT_INPUT_DIR / "blood_struct_pca.npy",
    }

    missing = [name for name, path in required.items() if not path.exists()]
    missing_exported = [name for name, path in exported.items() if not path.exists()]
    repo_status = check_official_repo(repo_dir)
    manifest: dict[str, Any] = {
        "method": "scSHEFT",
        "official_repo": "https://github.com/CSUBioGroup/scSHEFT",
        "status": "ready_for_h5ad_export" if not missing else "missing_project_inputs",
        "input_files": {name: file_info(path) for name, path in required.items()},
        "scsheft_prepared_files": {name: file_info(path) for name, path in exported.items()},
        "missing_inputs": missing,
        "missing_prepared_files": missing_exported,
        "python_package_check": check_python_packages(),
        "official_repo_check": repo_status,
        "decision": (
            "官方代码可用；本项目需要先导出 h5ad/npy，再调用官方 scSHEFT 类，"
            "最后把模型输出转成统一接口契约。"
        ),
    }

    if not missing:
        atac_shape = read_mtx_shape(required["atac_counts_mtx"])
        rna_shape = read_mtx_shape(required["rna_counts_mtx"])
        atac_barcodes = count_lines(required["atac_barcodes_txt"])
        rna_barcodes = count_lines(required["rna_barcodes_txt"])
        atac_peak_header = read_csv_header(required["atac_peaks_csv"])
        rna_gene_header = read_csv_header(required["rna_genes_csv"])
        atac_cell_header = read_csv_header(required["atac_cells_csv"])
        rna_cell_header = read_csv_header(required["rna_cells_csv"])

        manifest["shapes"] = {
            "atac_counts_mtx": atac_shape,
            "rna_counts_mtx": rna_shape,
            "atac_barcodes_n": atac_barcodes,
            "rna_barcodes_n": rna_barcodes,
            "atac_peaks_header": atac_peak_header,
            "rna_genes_header": rna_gene_header,
            "atac_cells_header": atac_cell_header,
            "rna_cells_header": rna_cell_header,
        }
        manifest["consistency_checks"] = {
            "atac_mtx_cols_match_barcodes": atac_shape["n_cols"] == atac_barcodes,
            "rna_mtx_cols_match_barcodes": rna_shape["n_cols"] == rna_barcodes,
        }

    return manifest


def write_readme(path: Path) -> None:
    text = """# method_B_scSHEFT

本目录用于模块 1 的方法 B（scSHEFT）。

当前已完成：
- 复用 `multiomics_output_v2/glue_input/` 下的 RNA / ATAC 原始输入
- 生成 `input_manifest.json`
- 明确官方仓库 `CSUBioGroup/scSHEFT` 可作为方法 B 算法本体

下一步执行顺序：
1. 运行 `module1_method_B_export_scSHEFT_input.R`，导出 RNA / Gene Activity 共享基因矩阵
2. 在 Python 3.9/3.10 环境中安装官方依赖
3. 将 MatrixMarket 输入转换成官方要求的 h5ad/npy
4. 运行 `run_official_scSHEFT_template.py`
5. 产出统一接口文件：
   - `shared_embedding.rds`
   - `atac_topk_rna_neighbors.rds`
   - `atac_label_transfer.csv`
   - `timing.json`
"""
    path.write_text(text, encoding="utf-8")


def write_official_notes(path: Path) -> None:
    text = """# 官方 scSHEFT 仓库可用性判断

仓库：<https://github.com/CSUBioGroup/scSHEFT>

结论：可以用，作为方法 B 的算法本体是合适的。

需要注意：
- `demo.py` 是路径占位模板，不能直接原样运行。
- 官方代码要求 `.h5ad` / `.npy` 输入；本项目目前已有 peak MatrixMarket 和 RDS，需要转换。
- 默认 Python 3.13 不适合安装 `torch==1.13.0`，建议新建 Python 3.9/3.10 环境。
- `scSHEFT.preprocess(..., stc_emb_inputs=True)` 会从 `./data/{exp_id}_struct_pca.npy` 读取结构 PCA，所以运行目录和文件名要对齐，或者在模板中用 `stc_emb_inputs=False` 让源码现场计算。

本项目对接路线：
1. `module1_method_B_export_scSHEFT_input.R` 导出共享基因 RNA / Gene Activity 矩阵。
2. Python 环境把 RNA、Gene Activity、raw ATAC peak 写成 h5ad，并计算 ATAC LSI / structural PCA。
3. 调用官方 `scSHEFT` 类训练。
4. 从 `model.source_features`、`model.target_features`、`model.target_predictions` 导出统一接口。
"""
    path.write_text(text, encoding="utf-8")


def write_runner_template(path: Path) -> None:
    text = '''"""项目方法 B：官方 scSHEFT 运行模板。

把本文件复制或直接在 scSHEFT Python 环境中运行。运行前需：
1. 克隆 https://github.com/CSUBioGroup/scSHEFT
2. 准备 scsheft_input/blood-*.h5ad 和 blood_*.npy
3. 把 SCSHEFT_REPO 改为官方仓库本地路径
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SCSHEFT_REPO = Path("D:/tools/scSHEFT")
DATA_DIR = Path("D:/SCAVENGE-reproducibility-main/blood_data")
IN_DIR = DATA_DIR / "multiomics_module1_pairing" / "method_B_scSHEFT" / "scsheft_input"
OUT_DIR = DATA_DIR / "multiomics_module1_pairing" / "method_B_scSHEFT"
EXP_ID = "blood"

sys.path.insert(0, str(SCSHEFT_REPO))
from scSHEFT import scSHEFT  # noqa: E402


def main() -> int:
    start = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    adata_raw_target = sc.read_h5ad(IN_DIR / "blood_multiome_atac.h5ad")
    adata_target = sc.read_h5ad(IN_DIR / "blood-GAM.h5ad")
    adata_source = sc.read_h5ad(IN_DIR / "blood-GEM.h5ad")
    target_raw_emb = np.load(IN_DIR / "blood_lsi.npy")

    ppd = {
        "binz": False,
        "hvg_num": None,
        "lognorm": [False, True],
        "scale_per_batch": True,
        "type_label": "cell_type",
        "knn": 10,
    }

    model = scSHEFT(
        encoder_type="linear",
        use_struct=True,
        stc_w=0.5,
        stc_cutoff=100,
        n_latent=128,
        bn=False,
        dr=0.2,
        cont_w=0.1,
        cont_tau=0.1,
        align_w=0.2,
        align_p=0.8,
        align_cutoff=100,
        center_w=5,
        anchor_w=0.5,
    )

    model.preprocess(
        EXP_ID,
        [adata_source, adata_target, adata_raw_target],
        target_raw_emb,
        ppd,
        adata_adt_inputs=None,
        stc_emb_inputs=False,
    )
    model.train(
        batch_size=256,
        training_steps=1000,
        lr=0.01,
        weight_decay=0,
        log_step=50,
        eval_target=True,
        eval_top_k=1,
    )
    model.eval(inplace=True)

    target_pred_type = model.annotate()
    target_prob = model.target_predictions.max(axis=1)
    label_df = pd.DataFrame({
        "atac_cell": adata_target.obs_names.to_numpy(),
        "original.label": adata_target.obs["cell_type"].astype(str).to_numpy(),
        "predicted.id": target_pred_type,
        "score": target_prob,
    })
    label_df.to_csv(OUT_DIR / "atac_label_transfer.csv", index=False)

    shared = np.vstack([model.source_features, model.target_features])
    shared_df = pd.DataFrame(shared, index=list(adata_source.obs_names) + list(adata_target.obs_names))
    shared_df.insert(0, "modality", ["rna"] * adata_source.n_obs + ["atac"] * adata_target.n_obs)
    shared_df.insert(0, "cell", shared_df.index)
    shared_df.to_csv(OUT_DIR / "shared_embedding.csv", index=False)

    timing = {
        "method": "scSHEFT official CSUBioGroup implementation",
        "official_repo": str(SCSHEFT_REPO),
        "total": time.time() - start,
        "note": "shared_embedding.csv should be converted to shared_embedding.rds for the final R interface contract.",
    }
    (OUT_DIR / "timing.json").write_text(json.dumps(timing, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''
    path.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="仅做输入预检并生成 manifest（默认行为）",
    )
    parser.add_argument(
        "--repo-dir",
        type=Path,
        default=None,
        help="官方 scSHEFT 仓库本地路径，用于检查是否已克隆",
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SCSHEFT_INPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest(args.repo_dir)
    (OUT_DIR / "input_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    write_readme(OUT_DIR / "README.md")
    write_official_notes(OUT_DIR / "OFFICIAL_REPO_NOTES.md")
    write_runner_template(OUT_DIR / "run_official_scSHEFT_template.py")

    print("========================================")
    print("方法 B · scSHEFT 输入预检完成")
    print("========================================")
    print("输出目录:", OUT_DIR)
    print("manifest :", OUT_DIR / "input_manifest.json")
    print("状态     :", manifest["status"])
    print("官方代码 :", "可用，但需本地克隆和 Python 3.9/3.10 环境")
    if manifest["missing_inputs"]:
        print("缺失输入 :", ", ".join(manifest["missing_inputs"]))
        return 1

    if manifest["missing_prepared_files"]:
        print("待生成 scSHEFT 输入:", ", ".join(manifest["missing_prepared_files"]))
        print("先运行: source('D:/RStudio/work/module1_method_B_export_scSHEFT_input.R')")

    checks = manifest.get("consistency_checks", {})
    for key, value in checks.items():
        print(f"{key}: {value}")

    if not args.prepare_only:
        print("已生成官方 scSHEFT 运行模板；正式训练需在 scSHEFT Python 环境中执行。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
