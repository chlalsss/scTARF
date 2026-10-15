"""Prepare NeurIPS 2021 CITE-seq GEX as a clean scDRS input h5ad.

Input:
  - GEO/Open Problems processed CITE h5ad with GEX + ADT features.

Output:
  - GEX-only AnnData.
  - X is raw counts copied from layers["counts"].
  - var_names are unique gene symbols.
  - obs keeps cell_type / donor / site / batch and related metadata.

This script intentionally does not run scDRS.
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
    parser.add_argument("--out-h5ad", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--granja-ga-genes", required=True)
    parser.add_argument("--trait-genes", required=True)
    return parser.parse_args()


def read_gene_list(path: str | Path) -> set[str]:
    return {
        line.strip()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    }


def write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)


def sparse_row_sum(mat) -> np.ndarray:
    if sparse.issparse(mat):
        return np.asarray(mat.sum(axis=1)).ravel()
    return np.asarray(mat).sum(axis=1)


def sparse_col_sum(mat) -> np.ndarray:
    if sparse.issparse(mat):
        return np.asarray(mat.sum(axis=0)).ravel()
    return np.asarray(mat).sum(axis=0)


def sparse_row_nnz(mat) -> np.ndarray:
    if sparse.issparse(mat):
        return np.asarray((mat > 0).sum(axis=1)).ravel()
    return np.count_nonzero(np.asarray(mat) > 0, axis=1)


def matrix_nonzero_sample_stats(mat, max_nnz: int = 1_000_000) -> dict[str, float]:
    if sparse.issparse(mat):
        data = mat.data
    else:
        arr = np.asarray(mat)
        data = arr[arr != 0]
    if data.size == 0:
        return {
            "nonzero_median": 0.0,
            "nonzero_p95": 0.0,
            "nonzero_max": 0.0,
            "nonzero_intlike": True,
        }
    if data.size > max_nnz:
        rng = np.random.default_rng(20260514)
        idx = rng.choice(data.size, size=max_nnz, replace=False)
        data = data[idx]
    data = np.asarray(data, dtype=np.float64)
    return {
        "nonzero_median": float(np.median(data)),
        "nonzero_p95": float(np.quantile(data, 0.95)),
        "nonzero_max": float(np.max(data)),
        "nonzero_intlike": bool(np.allclose(data[: min(100_000, data.size)], np.round(data[: min(100_000, data.size)]))),
    }


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_h5ad = Path(args.out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)

    print("Reading source CITE h5ad:", args.input_h5ad, flush=True)
    src = ad.read_h5ad(args.input_h5ad)
    print(f"  source shape: {src.n_obs} cells x {src.n_vars} features", flush=True)

    if "feature_types" not in src.var.columns:
        raise ValueError("source h5ad is missing var['feature_types']")
    if "counts" not in src.layers:
        raise ValueError("source h5ad is missing layers['counts']; raw counts are required for scDRS input")

    feature_types = src.var["feature_types"].astype(str)
    gex_mask = feature_types.eq("GEX").to_numpy()
    n_gex_before_zero_filter = int(gex_mask.sum())
    if n_gex_before_zero_filter == 0:
        raise ValueError("No GEX features found")

    gex_gene_names = pd.Index(src.var_names[gex_mask].astype(str))
    if not gex_gene_names.is_unique:
        dup = sorted(gex_gene_names[gex_gene_names.duplicated()].unique())[:20]
        raise ValueError(f"GEX gene names are not unique after ADT removal. Examples: {dup}")

    obs = src.obs.copy()
    obs["barcode"] = src.obs_names.astype(str)
    obs.index = pd.Index(src.obs_names.astype(str), name="cell")

    var = src.var.loc[gex_mask].copy()
    var["gene_symbol"] = gex_gene_names.to_numpy()
    var.index = pd.Index(gex_gene_names, name="gene_symbol")

    counts = src.layers["counts"][:, gex_mask]
    if sparse.issparse(counts):
        counts = counts.tocsr().astype(np.float32)
    else:
        counts = sparse.csr_matrix(np.asarray(counts, dtype=np.float32))

    gene_sum = sparse_col_sum(counts)
    keep_gene = gene_sum > 0
    n_zero_genes = int((~keep_gene).sum())
    if n_zero_genes > 0:
        counts = counts[:, keep_gene].tocsr()
        var = var.iloc[np.where(keep_gene)[0]].copy()

    cell_sum = sparse_row_sum(counts)
    cell_n_genes = sparse_row_nnz(counts)
    obs["scdrs_input_total_counts"] = cell_sum.astype(np.float64)
    obs["scdrs_input_n_genes"] = cell_n_genes.astype(np.int32)

    out = ad.AnnData(X=counts, obs=obs, var=var)
    out.uns["module4_preprocessing"] = {
        "source": "GEO GSE194122 Open Problems NeurIPS 2021 CITE BMMC",
        "feature_filter": "feature_types == GEX",
        "x_source": "layers['counts']",
        "scdrs_was_run": False,
    }

    if not out.var_names.is_unique:
        raise ValueError("Output var_names are not unique")
    if "cell_type" not in out.obs.columns:
        raise ValueError("Output obs is missing cell_type column")

    print("Writing:", out_h5ad, flush=True)
    out.write_h5ad(out_h5ad, compression="gzip")

    granja_ga_genes = read_gene_list(args.granja_ga_genes)
    out_genes = set(map(str, out.var_names))
    trait_genes = pd.read_csv(args.trait_genes)

    gene_overlap = pd.DataFrame(
        [
            {
                "comparison": "CITE_scDRS_GEX_vs_Granja_GA",
                "n_cite_scdrs_gex_genes": len(out_genes),
                "n_granja_ga_genes": len(granja_ga_genes),
                "n_overlap": len(out_genes & granja_ga_genes),
                "pct_cite_scdrs_covered": len(out_genes & granja_ga_genes) / max(len(out_genes), 1),
                "pct_granja_ga_covered": len(out_genes & granja_ga_genes) / max(len(granja_ga_genes), 1),
            }
        ]
    )
    write_csv(gene_overlap, out_dir / "cite_scdrs_gene_overlap.csv")

    trait_rows = []
    for trait, sub in trait_genes.groupby("trait", dropna=False):
        genes = set(sub["gene"].astype(str))
        trait_rows.append(
            {
                "trait": trait,
                "n_genes": len(genes),
                "n_in_cite_scdrs_gex": len(genes & out_genes),
                "n_in_granja_ga": len(genes & granja_ga_genes),
                "n_in_both_cite_gex_and_granja_ga": len(genes & out_genes & granja_ga_genes),
            }
        )
    trait_overlap = pd.DataFrame(trait_rows)
    write_csv(trait_overlap, out_dir / "cite_scdrs_trait_gene_overlap.csv")

    cell_counts = (
        out.obs["cell_type"]
        .astype(str)
        .value_counts()
        .rename_axis("cell_type")
        .reset_index(name="n_cells")
    )
    cell_counts["fraction"] = cell_counts["n_cells"] / out.n_obs
    write_csv(cell_counts, out_dir / "cite_scdrs_celltype_counts.csv")

    preferred_meta_cols = [
        "barcode",
        "cell_type",
        "DonorNumber",
        "Site",
        "batch",
        "Samplename",
        "Modality",
        "DonorID",
        "DonorAge",
        "DonorGender",
        "DonorRace",
        "Ethnicity",
        "scdrs_input_total_counts",
        "scdrs_input_n_genes",
    ]
    meta_cols = [col for col in preferred_meta_cols if col in out.obs.columns]
    meta_out = out.obs[meta_cols].copy()
    meta_out.insert(0, "cell", out.obs_names.astype(str))
    meta_out.to_csv(out_dir / "cite_scdrs_cell_meta.csv.gz", index=False)

    stats = matrix_nonzero_sample_stats(out.X)
    summary_rows = [
        ("source_h5ad", str(Path(args.input_h5ad).resolve())),
        ("out_h5ad", str(out_h5ad.resolve())),
        ("n_cells", out.n_obs),
        ("n_source_features", src.n_vars),
        ("n_source_gex_features", n_gex_before_zero_filter),
        ("n_output_genes", out.n_vars),
        ("n_zero_count_gex_genes_removed", n_zero_genes),
        ("x_source", "layers['counts']"),
        ("var_names_unique", out.var_names.is_unique),
        ("cell_type_column_present", "cell_type" in out.obs.columns),
        ("n_cell_types", int(out.obs["cell_type"].astype(str).nunique())),
        ("n_obs_columns", len(out.obs.columns)),
        ("n_var_columns", len(out.var.columns)),
        ("matrix_nnz", int(out.X.nnz) if sparse.issparse(out.X) else int(np.count_nonzero(out.X))),
        ("matrix_density", float((out.X.nnz if sparse.issparse(out.X) else np.count_nonzero(out.X)) / (out.n_obs * out.n_vars))),
        ("cell_total_counts_median", float(np.median(cell_sum))),
        ("cell_total_counts_p05", float(np.quantile(cell_sum, 0.05))),
        ("cell_total_counts_p95", float(np.quantile(cell_sum, 0.95))),
        ("cell_n_genes_median", float(np.median(cell_n_genes))),
        ("cell_n_genes_p05", float(np.quantile(cell_n_genes, 0.05))),
        ("cell_n_genes_p95", float(np.quantile(cell_n_genes, 0.95))),
        ("gene_overlap_with_granja_ga", int(gene_overlap["n_overlap"].iloc[0])),
        ("cell_meta_csv", str((out_dir / "cite_scdrs_cell_meta.csv.gz").resolve())),
    ]
    for key, value in stats.items():
        summary_rows.append((key, value))
    summary = pd.DataFrame(summary_rows, columns=["metric", "value"])
    write_csv(summary, out_dir / "cite_scdrs_input_summary.csv")

    report = [
        "CITE RNA pre-scDRS preprocessing report",
        "",
        f"Output h5ad: {out_h5ad.resolve()}",
        f"Shape: {out.n_obs:,} cells x {out.n_vars:,} GEX genes",
        "X source: raw counts from source layers['counts']",
        f"Gene names unique: {out.var_names.is_unique}",
        f"Zero-count GEX genes removed: {n_zero_genes:,}",
        f"Cell types: {out.obs['cell_type'].astype(str).nunique():,}",
        f"Gene overlap with Granja GA: {int(gene_overlap['n_overlap'].iloc[0]):,}",
        "",
        "Trait gene overlap in both CITE GEX and Granja GA:",
    ]
    for row in trait_overlap.itertuples(index=False):
        report.append(f"  {row.trait}: {row.n_in_both_cite_gex_and_granja_ga:,} / {row.n_genes:,}")
    report.extend(
        [
            "",
            "Next step:",
            "  Run scDRS on this h5ad, then check whether PLT enriches MK/E or erythroid-megakaryocyte cell states.",
        ]
    )
    (out_dir / "cite_scdrs_ready_report.txt").write_text("\n".join(report), encoding="utf-8")

    print(f"  output shape: {out.n_obs} cells x {out.n_vars} genes", flush=True)
    print(f"  output file size: {out_h5ad.stat().st_size / 1024**2:.1f} MiB", flush=True)
    print("Pre-scDRS input preparation complete.", flush=True)


if __name__ == "__main__":
    main()
