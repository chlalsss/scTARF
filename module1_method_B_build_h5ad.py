"""Build official scSHEFT h5ad/npy inputs from project MatrixMarket files.

Run this in a Python 3.9/3.10 environment with:
  numpy scipy pandas anndata scikit-learn

Inputs are produced by:
  source("D:/RStudio/work/module1_method_B_export_scSHEFT_input.R")

Outputs:
  blood-GEM.h5ad
  blood-GAM.h5ad
  blood_multiome_atac.h5ad
  blood_lsi.npy
  blood_struct_pca.npy
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize

try:
    import anndata as ad
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "缺少 anndata。请在 Python 3.9/3.10 环境中安装: "
        "pip install anndata pandas scipy scikit-learn numpy"
    ) from exc

try:
    import scanpy as sc
except ImportError:  # pragma: no cover
    sc = None


DATA_DIR = Path("D:/SCAVENGE-reproducibility-main/blood_data")
PROJECT_INPUT_DIR = DATA_DIR / "multiomics_output_v2"
GLUE_INPUT_DIR = PROJECT_INPUT_DIR / "glue_input"
SCSHEFT_INPUT_DIR = (
    DATA_DIR / "multiomics_module1_pairing" / "method_B_scSHEFT" / "scsheft_input"
)


def read_mtx_cells_by_features(path: Path) -> sparse.csr_matrix:
    """Read genes/peaks x cells MatrixMarket and return cells x features CSR."""
    print(f"  读取矩阵: {path}")
    mat = mmread(path)
    if not sparse.issparse(mat):
        mat = sparse.coo_matrix(mat)
    mat = mat.tocsr()
    print(f"    原始维度: {mat.shape[0]} features x {mat.shape[1]} cells; nnz={mat.nnz}")
    out = mat.T.tocsr()
    print(f"    转置后: {out.shape[0]} cells x {out.shape[1]} features")
    return out


def load_single_column(path: Path, column: str) -> pd.Index:
    df = pd.read_csv(path)
    if column not in df.columns:
        raise ValueError(f"{path} 缺少列: {column}")
    return pd.Index(df[column].astype(str), name=column)


def load_txt_lines(path: Path) -> pd.Index:
    with path.open("r", encoding="utf-8") as handle:
        vals = [line.strip() for line in handle if line.strip()]
    return pd.Index(vals)


def build_obs(path: Path) -> pd.DataFrame:
    obs = pd.read_csv(path)
    if "cell" not in obs.columns:
        raise ValueError(f"{path} 缺少 cell 列")
    if "cell_type" not in obs.columns:
        raise ValueError(f"{path} 缺少 cell_type 列")
    obs["cell"] = obs["cell"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs = obs.set_index("cell", drop=False)
    return obs


def align_obs(obs: pd.DataFrame, expected_cells: pd.Index, label: str) -> pd.DataFrame:
    expected_cells = pd.Index(expected_cells.astype(str))
    if obs.index.equals(expected_cells):
        return obs
    if set(obs.index) != set(expected_cells):
        missing_in_obs = list(expected_cells.difference(obs.index)[:5])
        missing_in_expected = list(obs.index.difference(expected_cells)[:5])
        raise ValueError(
            f"{label} 细胞集合不一致。missing_in_obs={missing_in_obs}, "
            f"missing_in_expected={missing_in_expected}"
        )
    return obs.loc[expected_cells].copy()


def write_h5ad(matrix: sparse.csr_matrix, obs: pd.DataFrame, var_names: pd.Index, path: Path) -> None:
    if matrix.shape[0] != obs.shape[0]:
        raise ValueError(f"{path.name}: 细胞数不匹配 matrix={matrix.shape[0]}, obs={obs.shape[0]}")
    if matrix.shape[1] != len(var_names):
        raise ValueError(f"{path.name}: feature 数不匹配 matrix={matrix.shape[1]}, var={len(var_names)}")

    var = pd.DataFrame(index=var_names.astype(str))
    var.index.name = None
    adata = ad.AnnData(X=matrix.astype(np.float32), obs=obs.copy(), var=var)
    adata.write_h5ad(path, compression="gzip")
    print(f"  已保存: {path} ({adata.n_obs} cells x {adata.n_vars} features)")


def compute_lsi(atac_counts: sparse.csr_matrix, n_components: int, random_state: int) -> np.ndarray:
    """Compute a standard sparse TF-IDF + SVD LSI embedding."""
    print(f"  计算 LSI: n_components={n_components}")
    x = atac_counts.copy().astype(np.float32)

    cell_sums = np.asarray(x.sum(axis=1)).ravel()
    cell_sums[cell_sums == 0] = 1.0
    tf = sparse.diags(1.0 / cell_sums, format="csr") @ x

    n_cells = x.shape[0]
    df = np.asarray((x > 0).sum(axis=0)).ravel()
    idf = np.log1p(n_cells / np.maximum(df, 1))
    tfidf = tf @ sparse.diags(idf.astype(np.float32), format="csr")
    tfidf = tfidf * 1e4
    tfidf.data = np.log1p(tfidf.data)

    svd = TruncatedSVD(n_components=n_components, random_state=random_state)
    emb = svd.fit_transform(tfidf)
    emb = normalize(emb, norm="l2")
    return emb.astype(np.float32)


def compute_struct_embedding(
    matrix: sparse.csr_matrix,
    obs: pd.DataFrame,
    peak_var: pd.DataFrame,
    n_components: int,
    random_state: int,
    method: str,
) -> np.ndarray:
    print(f"  计算 structural embedding: method={method}, n_components={n_components}")
    x = matrix.copy().astype(np.float32)

    if method == "scanpy_pca":
        if sc is None:
            raise SystemExit("缺少 scanpy，无法使用 scanpy_pca 计算 structural PCA。")
        adata = ad.AnnData(X=x, obs=obs.copy(), var=peak_var.copy())
        sc.settings.verbosity = 0
        sc.tl.pca(adata, n_comps=n_components, svd_solver="auto")
        return np.asarray(adata.obsm["X_pca"], dtype=np.float32)

    if method == "svd_fallback":
        svd = TruncatedSVD(n_components=n_components, random_state=random_state)
        return svd.fit_transform(x).astype(np.float32)

    raise ValueError(f"未知 structural embedding 方法: {method}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lsi-components", type=int, default=50)
    parser.add_argument("--struct-components", type=int, default=500)
    parser.add_argument("--random-state", type=int, default=9527)
    parser.add_argument(
        "--struct-method",
        choices=("scanpy_pca", "svd_fallback"),
        default="scanpy_pca",
        help="Official scSHEFT computes structural PCA from raw ATAC h5ad via scanpy.tl.pca.",
    )
    parser.add_argument(
        "--skip-raw-atac-h5ad",
        action="store_true",
        help="只生成 RNA/GAM h5ad 和 npy；raw ATAC h5ad 很大时可先跳过。",
    )
    args = parser.parse_args()

    start = time.time()
    SCSHEFT_INPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("========== 构建 scSHEFT h5ad/npy 输入 ==========")
    print("输出目录:", SCSHEFT_INPUT_DIR)

    genes = load_single_column(SCSHEFT_INPUT_DIR / "shared_genes.csv", "gene")
    rna_obs = build_obs(SCSHEFT_INPUT_DIR / "rna_cells.csv")
    atac_obs = build_obs(SCSHEFT_INPUT_DIR / "atac_cells.csv")
    rna_barcodes = load_txt_lines(GLUE_INPUT_DIR / "rna_barcodes.txt")
    atac_barcodes = load_txt_lines(GLUE_INPUT_DIR / "atac_barcodes.txt")

    print("\n========== RNA h5ad ==========")
    rna = read_mtx_cells_by_features(SCSHEFT_INPUT_DIR / "rna_counts_shared_genes.mtx")
    rna_obs = align_obs(rna_obs, rna_barcodes, "RNA")
    write_h5ad(rna, rna_obs, genes, SCSHEFT_INPUT_DIR / "blood-GEM.h5ad")
    del rna

    print("\n========== Gene Activity h5ad ==========")
    gam = read_mtx_cells_by_features(SCSHEFT_INPUT_DIR / "gene_activity_counts_shared_genes.mtx")
    atac_obs = align_obs(atac_obs, atac_barcodes, "ATAC")
    write_h5ad(gam, atac_obs, genes, SCSHEFT_INPUT_DIR / "blood-GAM.h5ad")
    del gam

    print("\n========== Raw ATAC peak h5ad + embeddings ==========")
    atac_peaks = pd.read_csv(GLUE_INPUT_DIR / "atac_peaks.csv")
    if "peak_id" not in atac_peaks.columns:
        raise ValueError("atac_peaks.csv 缺少 peak_id 列")
    peak_names = pd.Index(atac_peaks["peak_id"].astype(str), name="peak")

    raw_atac = read_mtx_cells_by_features(GLUE_INPUT_DIR / "atac_counts.mtx")
    if raw_atac.shape[0] != atac_obs.shape[0]:
        raise ValueError(f"raw ATAC 细胞数不匹配: {raw_atac.shape[0]} vs {atac_obs.shape[0]}")

    lsi = compute_lsi(raw_atac, args.lsi_components, args.random_state)
    np.save(SCSHEFT_INPUT_DIR / "blood_lsi.npy", lsi)
    print(f"  已保存: {SCSHEFT_INPUT_DIR / 'blood_lsi.npy'} {lsi.shape}")

    peak_var = pd.DataFrame(index=peak_names.astype(str))
    for col in ("chrom", "chromStart", "chromEnd"):
        if col in atac_peaks.columns:
            peak_var[col] = atac_peaks[col].to_numpy()

    struct = compute_struct_embedding(
        raw_atac,
        atac_obs,
        peak_var,
        args.struct_components,
        args.random_state,
        args.struct_method,
    )
    np.save(SCSHEFT_INPUT_DIR / "blood_struct_pca.npy", struct)
    print(f"  已保存: {SCSHEFT_INPUT_DIR / 'blood_struct_pca.npy'} {struct.shape}")

    if not args.skip_raw_atac_h5ad:
        adata_atac = ad.AnnData(X=raw_atac.astype(np.float32), obs=atac_obs.copy(), var=peak_var)
        adata_atac.write_h5ad(SCSHEFT_INPUT_DIR / "blood_multiome_atac.h5ad", compression="gzip")
        print(
            f"  已保存: {SCSHEFT_INPUT_DIR / 'blood_multiome_atac.h5ad'} "
            f"({adata_atac.n_obs} cells x {adata_atac.n_vars} peaks)"
        )

    manifest = {
        "script": "module1_method_B_build_h5ad.py",
        "outputs": {
            "rna_h5ad": str(SCSHEFT_INPUT_DIR / "blood-GEM.h5ad"),
            "gene_activity_h5ad": str(SCSHEFT_INPUT_DIR / "blood-GAM.h5ad"),
            "raw_atac_h5ad": str(SCSHEFT_INPUT_DIR / "blood_multiome_atac.h5ad"),
            "lsi_npy": str(SCSHEFT_INPUT_DIR / "blood_lsi.npy"),
            "struct_pca_npy": str(SCSHEFT_INPUT_DIR / "blood_struct_pca.npy"),
        },
        "n_genes": len(genes),
        "n_rna_cells": len(rna_obs),
        "n_atac_cells": len(atac_obs),
        "n_peaks": len(peak_names),
        "lsi_components": args.lsi_components,
        "struct_components": args.struct_components,
        "struct_method": args.struct_method,
        "total_seconds": time.time() - start,
    }
    (SCSHEFT_INPUT_DIR / "h5ad_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("\n======================================================")
    print("  scSHEFT h5ad/npy 输入构建完成")
    print("======================================================")
    print("  RNA cells:", len(rna_obs))
    print("  ATAC cells:", len(atac_obs))
    print("  genes:", len(genes))
    print("  peaks:", len(peak_names))
    print("  total min:", round((time.time() - start) / 60, 1))
    print("======================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
