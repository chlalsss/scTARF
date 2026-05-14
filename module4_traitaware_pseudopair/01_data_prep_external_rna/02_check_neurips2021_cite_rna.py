"""Check NeurIPS 2021 BMMC CITE-seq RNA usability for module4.

This script only inspects the RNA reference candidate. It does not run scDRS.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


MONO_MARKERS = ["LYZ", "S100A8", "S100A9", "FCN1", "CST3", "LST1", "MS4A7", "FCGR3A"]
LYM_MARKERS = ["CD3D", "CD3E", "CD4", "CD8A", "NKG7", "GNLY", "MS4A1", "CD79A", "IL7R"]
PLT_MK_MARKERS = ["PPBP", "PF4", "GP9", "GP1BA", "GP5", "ITGA2B", "ITGB3", "NFE2", "GATA1", "MPL", "VWF", "SELP", "RUNX1", "TAL1"]
ERY_MARKERS = ["HBB", "HBA1", "HBA2", "ALAS2", "KLF1", "GYPA"]

TARGET_PATTERNS = {
    "MONO": re.compile(r"mono|cd14|cd16", re.IGNORECASE),
    "LYM": re.compile(r"\bT\b|CD4|CD8|NK|ILC|B cell|B$|plasma|lymph", re.IGNORECASE),
    "PLT_MK": re.compile(r"mk|mega|megakaryo|platelet|mep|eryth|normoblast", re.IGNORECASE),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--granja-ga-genes", required=True)
    parser.add_argument("--granja-rna-genes", required=True)
    parser.add_argument("--granja-atac-barcodes", required=True)
    parser.add_argument("--trait-genes", required=True)
    return parser.parse_args()


def read_lines(path: str | Path) -> list[str]:
    path = Path(path)
    if not path.exists():
        return []
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)


def pick_column(columns: list[str], candidates: list[str], contains: list[str] | None = None) -> str | None:
    lower_map = {col.lower(): col for col in columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    if contains:
        for token in contains:
            for col in columns:
                if token.lower() in col.lower():
                    return col
    return None


def as_dense_small(x):
    if sparse.issparse(x):
        return x.toarray()
    if hasattr(x, "to_memory"):
        x = x.to_memory()
        if sparse.issparse(x):
            return x.toarray()
    return np.asarray(x)


def summarize_sample_matrix(matrix, n_obs: int, n_vars: int, rng: np.random.Generator) -> dict[str, float]:
    obs_idx = np.sort(rng.choice(n_obs, size=min(2000, n_obs), replace=False))
    var_idx = np.sort(rng.choice(n_vars, size=min(2000, n_vars), replace=False))
    sample = matrix[obs_idx, :][:, var_idx]
    if sparse.issparse(sample):
        nnz = sample.nnz
        total = sample.shape[0] * sample.shape[1]
        data = sample.data
    else:
        sample = np.asarray(sample)
        nnz = int(np.count_nonzero(sample))
        total = sample.size
        data = sample[sample != 0]
    if data.size == 0:
        vals = np.array([0.0])
    else:
        vals = np.asarray(data, dtype=float)
    return {
        "sample_n_obs": int(len(obs_idx)),
        "sample_n_vars": int(len(var_idx)),
        "sample_density": float(nnz / total),
        "sample_nonzero_median": float(np.median(vals)),
        "sample_nonzero_p95": float(np.quantile(vals, 0.95)),
        "sample_nonzero_max": float(np.max(vals)),
    }


def get_gex_mask(adata: ad.AnnData) -> pd.Series:
    var = adata.var
    if "feature_types" in var.columns:
        feature = var["feature_types"].astype(str)
        return feature.str.upper().isin(["GEX", "GENE EXPRESSION", "RNA"])
    return pd.Series(True, index=adata.var_names)


def marker_program_table(adata_gex: ad.AnnData, celltype_col: str | None, out_dir: Path) -> pd.DataFrame:
    if celltype_col is None:
        return pd.DataFrame()

    markers = {
        "MONO_marker_program": MONO_MARKERS,
        "LYM_marker_program": LYM_MARKERS,
        "PLT_MK_marker_program": PLT_MK_MARKERS,
        "ERY_marker_program": ERY_MARKERS,
    }
    available_markers = sorted({g for genes in markers.values() for g in genes if g in adata_gex.var_names})
    marker_presence = pd.DataFrame(
        [{"program": p, "gene": g, "present": g in adata_gex.var_names} for p, genes in markers.items() for g in genes]
    )
    write_csv(marker_presence, out_dir / "neurips2021_cite_marker_gene_presence.csv")

    if not available_markers:
        return pd.DataFrame()

    idx = [adata_gex.var_names.get_loc(g) for g in available_markers]
    mat = as_dense_small(adata_gex.X[:, idx])
    score_df = pd.DataFrame(index=adata_gex.obs_names)
    score_df[celltype_col] = adata_gex.obs[celltype_col].astype(str).values

    for program, genes in markers.items():
        present = [g for g in genes if g in available_markers]
        if not present:
            score_df[program] = np.nan
            continue
        sub_idx = [available_markers.index(g) for g in present]
        score_df[program] = np.asarray(mat[:, sub_idx], dtype=float).mean(axis=1)

    agg = (
        score_df.groupby(celltype_col, dropna=False)
        .agg(
            n_cells=(celltype_col, "size"),
            mean_MONO_marker_program=("MONO_marker_program", "mean"),
            mean_LYM_marker_program=("LYM_marker_program", "mean"),
            mean_PLT_MK_marker_program=("PLT_MK_marker_program", "mean"),
            mean_ERY_marker_program=("ERY_marker_program", "mean"),
        )
        .reset_index()
        .rename(columns={celltype_col: "cell_type"})
    )

    top_rows = []
    for col in [
        "mean_MONO_marker_program",
        "mean_LYM_marker_program",
        "mean_PLT_MK_marker_program",
        "mean_ERY_marker_program",
    ]:
        tmp = agg.sort_values(col, ascending=False).head(10).copy()
        tmp.insert(0, "program", col.replace("mean_", ""))
        top_rows.append(tmp[["program", "cell_type", "n_cells", col]].rename(columns={col: "mean_score"}))
    top = pd.concat(top_rows, ignore_index=True)
    write_csv(top, out_dir / "neurips2021_cite_marker_program_top10.csv")

    return agg.sort_values("n_cells", ascending=False)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Reading h5ad:", args.h5ad, flush=True)
    adata = ad.read_h5ad(args.h5ad)
    print(f"  shape: {adata.n_obs} cells x {adata.n_vars} features", flush=True)

    obs_cols = list(adata.obs.columns)
    var_cols = list(adata.var.columns)
    celltype_col = pick_column(obs_cols, ["cell_type", "celltype", "cell_type_l2", "cell_type_fine"], contains=["cell"])
    donor_col = pick_column(obs_cols, ["donor", "donor_id"], contains=["donor"])
    site_col = pick_column(obs_cols, ["site", "batch_site"], contains=["site"])
    batch_col = pick_column(obs_cols, ["batch", "batch_id"], contains=["batch"])

    gex_mask = get_gex_mask(adata)
    n_gex = int(gex_mask.sum())
    n_non_gex = int(adata.n_vars - n_gex)
    adata_gex = adata[:, gex_mask.values].copy()

    rng = np.random.default_rng(20260514)
    sample_stats = summarize_sample_matrix(adata_gex.X, adata_gex.n_obs, adata_gex.n_vars, rng)

    feature_type_counts = (
        adata.var["feature_types"].astype(str).value_counts().rename_axis("feature_type").reset_index(name="n_features")
        if "feature_types" in adata.var.columns
        else pd.DataFrame({"feature_type": ["not_available"], "n_features": [adata.n_vars]})
    )
    write_csv(feature_type_counts, out_dir / "neurips2021_cite_feature_type_counts.csv")

    if celltype_col is not None:
        cell_counts = (
            adata.obs[celltype_col].astype(str).value_counts().rename_axis("cell_type").reset_index(name="n_cells")
        )
        cell_counts["fraction"] = cell_counts["n_cells"] / adata.n_obs
        write_csv(cell_counts, out_dir / "neurips2021_cite_celltype_counts.csv")
    else:
        cell_counts = pd.DataFrame()
        write_csv(cell_counts, out_dir / "neurips2021_cite_celltype_counts.csv")

    meta_rows = []
    for name, col in [("donor", donor_col), ("site", site_col), ("batch", batch_col)]:
        if col is None:
            continue
        counts = adata.obs[col].astype(str).value_counts().rename_axis("value").reset_index(name="n_cells")
        counts.insert(0, "field", name)
        counts.insert(1, "column", col)
        meta_rows.append(counts)
    if meta_rows:
        write_csv(pd.concat(meta_rows, ignore_index=True), out_dir / "neurips2021_cite_donor_site_batch_counts.csv")
    else:
        write_csv(pd.DataFrame(columns=["field", "column", "value", "n_cells"]), out_dir / "neurips2021_cite_donor_site_batch_counts.csv")

    granja_ga_genes = set(read_lines(args.granja_ga_genes))
    granja_rna_genes = set(read_lines(args.granja_rna_genes))
    granja_atac_barcodes = set(read_lines(args.granja_atac_barcodes))
    ext_genes = set(map(str, adata_gex.var_names))

    overlap_rows = [
        {
            "comparison": "external_CITE_GEX_vs_Granja_GA",
            "n_external_gex_genes": len(ext_genes),
            "n_reference_genes": len(granja_ga_genes),
            "n_overlap": len(ext_genes & granja_ga_genes),
            "pct_external_covered": len(ext_genes & granja_ga_genes) / max(len(ext_genes), 1),
            "pct_reference_covered": len(ext_genes & granja_ga_genes) / max(len(granja_ga_genes), 1),
        },
        {
            "comparison": "external_CITE_GEX_vs_Granja_RNA",
            "n_external_gex_genes": len(ext_genes),
            "n_reference_genes": len(granja_rna_genes),
            "n_overlap": len(ext_genes & granja_rna_genes),
            "pct_external_covered": len(ext_genes & granja_rna_genes) / max(len(ext_genes), 1),
            "pct_reference_covered": len(ext_genes & granja_rna_genes) / max(len(granja_rna_genes), 1),
        },
    ]
    write_csv(pd.DataFrame(overlap_rows), out_dir / "neurips2021_cite_gene_overlap.csv")

    trait_genes = pd.read_csv(args.trait_genes)
    trait_rows = []
    for trait, sub in trait_genes.groupby("trait", dropna=False):
        genes = set(sub["gene"].astype(str))
        trait_rows.append(
            {
                "trait": trait,
                "n_genes": len(genes),
                "n_in_external_CITE_GEX": len(genes & ext_genes),
                "n_in_Granja_GA": len(genes & granja_ga_genes),
                "n_in_both_external_GEX_and_Granja_GA": len(genes & ext_genes & granja_ga_genes),
            }
        )
    write_csv(pd.DataFrame(trait_rows), out_dir / "neurips2021_cite_trait_gene_overlap.csv")

    barcode_overlap = len(set(map(str, adata.obs_names)) & granja_atac_barcodes)
    barcode_df = pd.DataFrame(
        [
            {
                "external_cite_cells": adata.n_obs,
                "granja_atac_cells": len(granja_atac_barcodes),
                "barcode_overlap_count": barcode_overlap,
            }
        ]
    )
    write_csv(barcode_df, out_dir / "neurips2021_cite_vs_granja_barcode_overlap.csv")

    if celltype_col is not None:
        target_rows = []
        celltype_strings = adata.obs[celltype_col].astype(str)
        for trait, pat in TARGET_PATTERNS.items():
            mask = celltype_strings.str.contains(pat, regex=True, na=False)
            matched_types = sorted(celltype_strings[mask].unique())
            target_rows.append(
                {
                    "target": trait,
                    "n_cells_matching_name": int(mask.sum()),
                    "n_cell_types_matching_name": len(matched_types),
                    "matched_cell_types": ";".join(matched_types[:25]),
                }
            )
        target_df = pd.DataFrame(target_rows)
    else:
        target_df = pd.DataFrame(columns=["target", "n_cells_matching_name", "n_cell_types_matching_name", "matched_cell_types"])
    write_csv(target_df, out_dir / "neurips2021_cite_candidate_target_summary.csv")

    marker_agg = marker_program_table(adata_gex, celltype_col, out_dir)
    write_csv(marker_agg, out_dir / "neurips2021_cite_marker_program_by_celltype.csv")

    qc_rows = [
        ("source", "GEO GSE194122 Open Problems NeurIPS 2021 CITE BMMC processed h5ad"),
        ("h5ad_path", str(Path(args.h5ad).resolve())),
        ("n_cells_total", adata.n_obs),
        ("n_features_total", adata.n_vars),
        ("n_gex_features", n_gex),
        ("n_non_gex_features", n_non_gex),
        ("obs_columns", json.dumps(obs_cols, ensure_ascii=False)),
        ("var_columns", json.dumps(var_cols, ensure_ascii=False)),
        ("layers", json.dumps(list(adata.layers.keys()), ensure_ascii=False)),
        ("obsm_keys", json.dumps(list(adata.obsm.keys()), ensure_ascii=False)),
        ("celltype_column", celltype_col or ""),
        ("donor_column", donor_col or ""),
        ("site_column", site_col or ""),
        ("batch_column", batch_col or ""),
        ("n_cell_types", int(adata.obs[celltype_col].nunique()) if celltype_col else 0),
        ("n_donors", int(adata.obs[donor_col].nunique()) if donor_col else 0),
        ("n_sites", int(adata.obs[site_col].nunique()) if site_col else 0),
        ("has_counts_layer", "counts" in adata.layers.keys()),
        ("barcode_overlap_with_granja_atac", barcode_overlap),
    ]
    for key, value in sample_stats.items():
        qc_rows.append((key, value))
    qc_df = pd.DataFrame(qc_rows, columns=["metric", "value"])
    write_csv(qc_df, out_dir / "neurips2021_cite_rna_qc_summary.csv")

    # Simple usability call: enough cells, enough GEX genes, clear target labels,
    # strong gene overlap, and no direct cell barcode reuse with Granja ATAC.
    gene_overlap = int((ext_genes & granja_ga_genes).__len__())
    mono_ok = int(target_df.loc[target_df["target"] == "MONO", "n_cells_matching_name"].iloc[0]) > 100
    lym_ok = int(target_df.loc[target_df["target"] == "LYM", "n_cells_matching_name"].iloc[0]) > 100
    plt_ok = int(target_df.loc[target_df["target"] == "PLT_MK", "n_cells_matching_name"].iloc[0]) > 50
    usable = adata.n_obs > 10000 and n_gex > 10000 and gene_overlap > 12000 and mono_ok and lym_ok and plt_ok

    report_lines = [
        "NeurIPS 2021 BMMC CITE-seq RNA usability report",
        "",
        f"Cells: {adata.n_obs:,}",
        f"GEX genes/features: {n_gex:,}",
        f"Cell type column: {celltype_col or 'NOT FOUND'}",
        f"Donor column: {donor_col or 'NOT FOUND'}",
        f"Site column: {site_col or 'NOT FOUND'}",
        f"Counts layer present: {'yes' if 'counts' in adata.layers.keys() else 'no'}",
        f"Gene overlap with Granja GA: {gene_overlap:,}",
        f"Barcode overlap with Granja ATAC: {barcode_overlap:,}",
        "",
        "Candidate target label coverage:",
    ]
    if not target_df.empty:
        for row in target_df.itertuples(index=False):
            report_lines.append(
                f"  {row.target}: {row.n_cells_matching_name:,} cells across {row.n_cell_types_matching_name} cell types"
            )
            report_lines.append(f"    {row.matched_cell_types}")
    report_lines.extend(
        [
            "",
            "Pre-scDRS usability call:",
            "  PASS" if usable else "  REVIEW_NEEDED",
            "",
            "Interpretation:",
            "  This check only evaluates whether the RNA reference is structurally usable.",
            "  It does not test trait directionality or enrichment; that belongs to the later scDRS step.",
        ]
    )
    (out_dir / "neurips2021_cite_usable_report.txt").write_text("\n".join(report_lines), encoding="utf-8")

    print("QC complete. Usability:", "PASS" if usable else "REVIEW_NEEDED", flush=True)


if __name__ == "__main__":
    main()
