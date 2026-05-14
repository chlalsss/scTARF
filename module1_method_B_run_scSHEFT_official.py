"""Run official CSUBioGroup/scSHEFT on this project's method B inputs.

Run in the scsheft conda environment:
  python module1_method_B_run_scSHEFT_official.py --repo-dir D:/tools/scSHEFT --steps 1000

This script uses:
  blood-GEM.h5ad
  blood-GAM.h5ad
  blood_lsi.npy

And prefers the official structural-information path when available:
  blood_multiome_atac.h5ad + stc_emb_inputs=False

Fallback:
  blood_struct_pca.npy + stc_emb_inputs=True
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sps
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from annoy import AnnoyIndex
from sklearn.neighbors import NearestNeighbors


DATA_DIR = Path("D:/SCAVENGE-reproducibility-main/blood_data")
OUT_DIR = DATA_DIR / "multiomics_module1_pairing" / "method_B_scSHEFT"
SCSHEFT_INPUT_DIR = OUT_DIR / "scsheft_input"
EXP_ID = "blood"


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"缺少文件: {path}")


def looks_like_raw_counts(matrix: sps.spmatrix | np.ndarray, sample_size: int = 100000) -> bool:
    """Heuristic: exported raw counts are non-negative and almost entirely integer-valued."""
    values = matrix.data if sps.issparse(matrix) else np.asarray(matrix).ravel()
    values = np.asarray(values)
    if values.size == 0:
        return False
    if values.size > sample_size:
        rng = np.random.default_rng(9527)
        values = values[rng.choice(values.size, size=sample_size, replace=False)]
    values = values[np.isfinite(values)]
    if values.size == 0 or np.any(values < 0):
        return False
    integer_like_ratio = np.mean(np.abs(values - np.round(values)) < 1e-6)
    return bool(integer_like_ratio >= 0.99 and np.max(values) >= 20)


def parse_bool_mode(mode: str, detected_default: bool | None = None) -> bool:
    mode = mode.lower()
    if mode == "true":
        return True
    if mode == "false":
        return False
    if mode == "auto":
        if detected_default is None:
            raise ValueError("auto mode requires detected_default")
        return detected_default
    raise ValueError(f"Unsupported mode: {mode}")


def write_topk_neighbors(
    source_features: np.ndarray,
    target_features: np.ndarray,
    rna_cells: np.ndarray,
    atac_cells: np.ndarray,
    out_csv: Path,
    k: int = 30,
) -> None:
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nn.fit(source_features)
    dists, idx = nn.kneighbors(target_features)

    rows = []
    for i, atac_cell in enumerate(atac_cells):
        for rank in range(k):
            rows.append(
                {
                    "atac_cell": atac_cell,
                    "rank": rank + 1,
                    "rna_cell": rna_cells[idx[i, rank]],
                    "distance": float(dists[i, rank]),
                }
            )
    pd.DataFrame(rows).to_csv(out_csv, index=False)


def patch_eval_single_worker(model: object, cls_dataset: object) -> None:
    """Patch official model.eval to avoid Windows multiprocessing pickle errors."""

    def eval_single_worker(batch_size: int = 500, inplace: bool = False):
        src_ds = cls_dataset(model.source_data, model.source_label_ids, binz=False, train=False)
        tgt_ds = cls_dataset(
            model.target_data,
            np.ones(model.n_target, dtype="int32"),
            binz=False,
            train=False,
        )

        src_dl = DataLoader(src_ds, batch_size=batch_size, shuffle=False, num_workers=0, drop_last=False)
        tgt_dl = DataLoader(tgt_ds, batch_size=batch_size, shuffle=False, num_workers=0, drop_last=False)

        model.encoder.eval()
        model.head.eval()

        def evaluate_domain(dataloader):
            features, predictions = [], []
            with torch.no_grad():
                for x, _y in dataloader:
                    x = x.to(model.device)
                    z = model.encoder(x)
                    h = nn.Softmax(dim=1)(model.head(z))
                    features.append(z.detach().cpu().numpy())
                    predictions.append(h.detach().cpu().numpy())
            return np.vstack(features), np.vstack(predictions)

        source_features_eval, source_predictions_eval = evaluate_domain(src_dl)
        target_features_eval, target_predictions_eval = evaluate_domain(tgt_dl)

        if inplace:
            model.source_features = source_features_eval
            model.target_features = target_features_eval
            model.source_predictions = source_predictions_eval
            model.target_predictions = target_predictions_eval
            model.features_combined = np.vstack([source_features_eval, target_features_eval])
            model.predictions_combined = np.vstack([source_predictions_eval, target_predictions_eval])
            return None
        return source_features_eval, target_features_eval, source_predictions_eval, target_predictions_eval

    model.eval = eval_single_worker


def make_annoy_anchor_finder(n_trees: int = 20, search_k: int = -1):
    """Return a memory-friendly mutual nearest-neighbor anchor finder."""

    def find_anchor_points_annoy(f_A: np.ndarray, f_B: np.ndarray, k: int = 5):
        f_A = np.asarray(f_A, dtype=np.float32)
        f_B = np.asarray(f_B, dtype=np.float32)
        dim = f_A.shape[1]
        if f_B.shape[1] != dim:
            raise ValueError(f"anchor dim mismatch: {f_A.shape} vs {f_B.shape}")

        print(
            f"Finding anchors with Annoy mutual NN "
            f"(n_A={f_A.shape[0]}, n_B={f_B.shape[0]}, dim={dim}, k={k}, trees={n_trees})"
        )

        b_index = AnnoyIndex(dim, "angular")
        for j in range(f_B.shape[0]):
            b_index.add_item(j, f_B[j])
        b_index.build(n_trees)

        a_to_b = []
        for i in range(f_A.shape[0]):
            a_to_b.append(b_index.get_nns_by_vector(f_A[i], k, search_k=search_k))

        a_index = AnnoyIndex(dim, "angular")
        for i in range(f_A.shape[0]):
            a_index.add_item(i, f_A[i])
        a_index.build(n_trees)

        b_to_a_sets = []
        for j in range(f_B.shape[0]):
            b_to_a_sets.append(set(a_index.get_nns_by_vector(f_B[j], k, search_k=search_k)))

        result = set()
        filtered_A_B_inds = [[] for _ in range(f_A.shape[0])]
        filtered_B_A_inds = [[] for _ in range(f_B.shape[0])]

        for i, b_neighbors in enumerate(a_to_b):
            for j in b_neighbors:
                if i in b_to_a_sets[j]:
                    result.add((i, j))
                    filtered_A_B_inds[i].append(j)
                    filtered_B_A_inds[j].append(i)

        print(f"Finding {len(result)} anchors")
        return filtered_A_B_inds, filtered_B_A_inds, result

    return find_anchor_points_annoy


def stratified_indices(labels: pd.Series, n: int, seed: int) -> np.ndarray:
    """Sample up to n rows with approximate class proportions."""
    rng = np.random.default_rng(seed)
    labels = labels.astype(str).reset_index(drop=True)
    total = len(labels)
    if n <= 0 or n >= total:
        return np.arange(total)

    counts = labels.value_counts()
    raw_quota = counts / total * n
    quota = np.floor(raw_quota).astype(int)
    quota[counts > 0] = np.maximum(quota[counts > 0], 1)

    while quota.sum() > n:
        candidates = quota[quota > 1].sort_values(ascending=False)
        if len(candidates) == 0:
            break
        quota.loc[candidates.index[0]] -= 1

    while quota.sum() < n:
        room = counts - quota
        candidates = room[room > 0].sort_values(ascending=False)
        if len(candidates) == 0:
            break
        quota.loc[candidates.index[0]] += 1

    selected: list[np.ndarray] = []
    for label, k in quota.items():
        pool = np.flatnonzero(labels.to_numpy() == label)
        if k > 0:
            selected.append(rng.choice(pool, size=min(k, len(pool)), replace=False))
    idx = np.concatenate(selected) if selected else np.array([], dtype=int)
    rng.shuffle(idx)
    return np.sort(idx)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-dir", type=Path, default=Path("D:/tools/scSHEFT"))
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--latent-dim", type=int, default=128)
    parser.add_argument("--log-step", type=int, default=50)
    parser.add_argument("--topk", type=int, default=30)
    parser.add_argument(
        "--struct-source",
        choices=("official_raw_atac", "precomputed"),
        default="official_raw_atac",
        help="Prefer official raw-ATAC PCA when blood_multiome_atac.h5ad is available.",
    )
    parser.add_argument(
        "--subsample-cells",
        type=int,
        default=0,
        help="Stratified cells per modality. 0 means full data. Official algorithm remains unchanged.",
    )
    parser.add_argument("--subsample-seed", type=int, default=9527)
    parser.add_argument(
        "--anchor-weight",
        type=float,
        default=0.5,
        help="Official anchor_w. Use 0 for smoke tests or memory-limited CPU runs.",
    )
    parser.add_argument(
        "--anchor-search",
        choices=("annoy", "official"),
        default="official",
        help="Anchor search implementation. Default keeps the official dense cdist version.",
    )
    parser.add_argument("--anchor-trees", type=int, default=20)
    parser.add_argument(
        "--rna-lognorm",
        choices=("auto", "true", "false"),
        default="auto",
        help=(
            "Whether to log-normalize RNA inside scSHEFT preprocess. "
            "Use auto to normalize when blood-GEM.h5ad looks like raw counts."
        ),
    )
    parser.add_argument(
        "--atac-lognorm",
        choices=("auto", "true", "false"),
        default="true",
        help="Whether to log-normalize Gene Activity inside scSHEFT preprocess.",
    )
    args = parser.parse_args()

    start = time.time()
    repo_dir = args.repo_dir.resolve()
    if not repo_dir.exists():
        raise FileNotFoundError(f"官方 scSHEFT 仓库不存在: {repo_dir}")

    for name in ("scSHEFT.py", "models.py", "utils.py", "dataset.py"):
        require_file(repo_dir / name)

    rna_h5ad = SCSHEFT_INPUT_DIR / "blood-GEM.h5ad"
    gam_h5ad = SCSHEFT_INPUT_DIR / "blood-GAM.h5ad"
    lsi_npy = SCSHEFT_INPUT_DIR / "blood_lsi.npy"
    struct_npy = SCSHEFT_INPUT_DIR / "blood_struct_pca.npy"
    raw_atac_h5ad = SCSHEFT_INPUT_DIR / "blood_multiome_atac.h5ad"
    for path in (rna_h5ad, gam_h5ad, lsi_npy):
        require_file(path)
    if args.struct_source == "precomputed":
        require_file(struct_npy)

    compat_dir = Path(__file__).resolve().parent / "scsheft_compat"
    sys.path.insert(0, str(compat_dir))
    sys.path.insert(0, str(repo_dir))
    os.chdir(repo_dir)

    import scSHEFT as scsheft_module  # noqa: E402
    from dataset import ClsDataset  # noqa: E402
    from scSHEFT import scSHEFT  # noqa: E402

    if args.anchor_weight != 0 and args.anchor_search == "annoy":
        scsheft_module.find_anchor_points = make_annoy_anchor_finder(
            n_trees=args.anchor_trees,
            search_k=-1,
        )
        print(f"  使用 Annoy anchor 搜索替换官方 dense cdist；trees={args.anchor_trees}")

    print("========== 加载 scSHEFT 输入 ==========")
    adata_source = sc.read_h5ad(rna_h5ad)
    adata_target = sc.read_h5ad(gam_h5ad)
    target_raw_emb = np.load(lsi_npy)
    rna_detected_raw = looks_like_raw_counts(adata_source.X)
    atac_detected_raw = looks_like_raw_counts(adata_target.X)
    rna_lognorm = parse_bool_mode(args.rna_lognorm, detected_default=rna_detected_raw)
    atac_lognorm = parse_bool_mode(args.atac_lognorm, detected_default=atac_detected_raw)

    if not adata_target.obs_names.equals(pd.Index(adata_target.obs["cell"].astype(str))):
        adata_target = adata_target[adata_target.obs["cell"].astype(str).to_numpy(), :].copy()

    target_cells = pd.Index(adata_target.obs_names.astype(str))

    adata_raw_target = None
    struct_emb = None
    stc_emb_inputs = False
    if args.struct_source == "official_raw_atac" and raw_atac_h5ad.exists():
        adata_raw_target = sc.read_h5ad(raw_atac_h5ad)
        raw_cells = pd.Index(adata_raw_target.obs_names.astype(str))
        if set(raw_cells) != set(target_cells):
            raise ValueError("raw ATAC h5ad 与 GAM h5ad 的细胞集合不一致")
        if not raw_cells.equals(target_cells):
            adata_raw_target = adata_raw_target[target_cells, :].copy()
    else:
        require_file(struct_npy)
        struct_emb = np.load(struct_npy)
        stc_emb_inputs = True

    if args.subsample_cells > 0:
        print(f"  子采样: 每个模态最多 {args.subsample_cells} cells，按 cell_type 分层")
        src_idx = stratified_indices(
            adata_source.obs["cell_type"],
            args.subsample_cells,
            args.subsample_seed,
        )
        tgt_idx = stratified_indices(
            adata_target.obs["cell_type"],
            args.subsample_cells,
            args.subsample_seed + 1,
        )
        adata_source = adata_source[src_idx, :].copy()
        adata_target = adata_target[tgt_idx, :].copy()
        target_raw_emb = target_raw_emb[tgt_idx, :]
        if adata_raw_target is not None:
            adata_raw_target = adata_raw_target[tgt_idx, :].copy()
        else:
            struct_emb = struct_emb[tgt_idx, :]

    struct_dst = None
    if stc_emb_inputs:
        # Official scSHEFT loads precomputed structural PCA from ./data/{exp_id}_struct_pca.npy.
        repo_data_dir = repo_dir / "data"
        repo_data_dir.mkdir(exist_ok=True)
        struct_dst = repo_data_dir / f"{EXP_ID}_struct_pca.npy"
        np.save(struct_dst, struct_emb)

    print("  RNA:", adata_source.shape)
    print("  Gene Activity:", adata_target.shape)
    print("  LSI:", target_raw_emb.shape)
    print(f"  RNA looks like raw counts: {rna_detected_raw} -> lognorm={rna_lognorm}")
    print(f"  GAM looks like raw counts: {atac_detected_raw} -> lognorm={atac_lognorm}")
    if stc_emb_inputs:
        print("  struct PCA (precomputed):", struct_emb.shape)
    else:
        print("  raw ATAC for official structural PCA:", adata_raw_target.shape)

    if not hasattr(sps.csr_matrix, "A"):
        sps.csr_matrix.A = property(lambda self: self.toarray())
    if not hasattr(sps.csc_matrix, "A"):
        sps.csc_matrix.A = property(lambda self: self.toarray())

    ppd = {
        "binz": False,
        "hvg_num": None,
        "lognorm": [rna_lognorm, atac_lognorm],
        "scale_per_batch": True,
        "type_label": "cell_type",
        "knn": 10,
    }

    model = scSHEFT(
        encoder_type="linear",
        use_struct=True,
        stc_w=0.5,
        stc_cutoff=100,
        n_latent=args.latent_dim,
        bn=False,
        dr=0.2,
        cont_w=0.1,
        cont_tau=0.1,
        align_w=0.2,
        align_p=0.8,
        align_cutoff=100,
        center_w=5,
        anchor_w=args.anchor_weight,
    )
    patch_eval_single_worker(model, ClsDataset)

    print("\n========== scSHEFT preprocess ==========")
    t0 = time.time()
    model.preprocess(
        EXP_ID,
        [adata_source, adata_target, adata_raw_target],
        target_raw_emb,
        ppd,
        adata_adt_inputs=None,
        stc_emb_inputs=stc_emb_inputs,
    )
    preprocess_seconds = time.time() - t0
    print(f"  preprocess seconds: {preprocess_seconds:.1f}")

    print("\n========== scSHEFT train ==========")
    t0 = time.time()
    model.train(
        batch_size=args.batch_size,
        training_steps=args.steps,
        lr=0.01,
        weight_decay=0,
        log_step=args.log_step,
        eval_target=True,
        eval_top_k=1,
    )
    train_seconds = time.time() - t0

    print("\n========== 导出结果 ==========")
    model.eval(inplace=True)
    target_pred_type = model.annotate()
    target_score = model.target_predictions.max(axis=1)

    label_df = pd.DataFrame(
        {
            "atac_cell": adata_target.obs_names.astype(str),
            "original.label": adata_target.obs["cell_type"].astype(str).to_numpy(),
            "predicted.id": target_pred_type,
            "score": target_score,
        }
    )
    label_df.to_csv(OUT_DIR / "atac_label_transfer.csv", index=False)
    print("  已保存:", OUT_DIR / "atac_label_transfer.csv")

    shared = np.vstack([model.source_features, model.target_features])
    shared_df = pd.DataFrame(shared, columns=[f"latent_{i + 1}" for i in range(shared.shape[1])])
    shared_df.insert(0, "modality", ["rna"] * adata_source.n_obs + ["atac"] * adata_target.n_obs)
    shared_df.insert(
        0,
        "cell",
        list(adata_source.obs_names.astype(str)) + list(adata_target.obs_names.astype(str)),
    )
    shared_df.to_csv(OUT_DIR / "shared_embedding.csv", index=False)
    print("  已保存:", OUT_DIR / "shared_embedding.csv")

    if args.topk > 0:
        write_topk_neighbors(
            model.source_features,
            model.target_features,
            adata_source.obs_names.astype(str).to_numpy(),
            adata_target.obs_names.astype(str).to_numpy(),
            OUT_DIR / "atac_topk_rna_neighbors.csv",
            k=args.topk,
        )
        print("  已保存:", OUT_DIR / "atac_topk_rna_neighbors.csv")
    else:
        print("  跳过 atac_topk_rna_neighbors.csv (topk=0)")

    ok = label_df.dropna(subset=["original.label", "predicted.id"])
    agreement = float((ok["original.label"] == ok["predicted.id"]).mean())
    timing = {
        "method": "scSHEFT official CSUBioGroup implementation",
        "repo_dir": str(repo_dir),
        "steps": args.steps,
        "batch_size": args.batch_size,
        "latent_dim": args.latent_dim,
        "anchor_weight": args.anchor_weight,
        "anchor_search": args.anchor_search,
        "anchor_trees": args.anchor_trees,
        "rna_lognorm_mode": args.rna_lognorm,
        "rna_detected_raw_counts": rna_detected_raw,
        "rna_lognorm_applied": rna_lognorm,
        "atac_lognorm_mode": args.atac_lognorm,
        "atac_detected_raw_counts": atac_detected_raw,
        "atac_lognorm_applied": atac_lognorm,
        "subsample_cells": args.subsample_cells,
        "subsample_seed": args.subsample_seed,
        "struct_source": args.struct_source,
        "preprocess_seconds": preprocess_seconds,
        "train_seconds": train_seconds,
        "total_seconds": time.time() - start,
        "overall_agreement": agreement,
        "note": "Primary outputs are CSV; convert to RDS in a follow-up R bridge if strict interface files are needed.",
    }
    (OUT_DIR / "timing.json").write_text(
        json.dumps(timing, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print("  已保存:", OUT_DIR / "timing.json")

    print("\n======================================================")
    print("  方法 B · scSHEFT 官方代码运行完成")
    print("======================================================")
    print("  overall agreement:", round(agreement * 100, 2), "%")
    print("  total min:", round((time.time() - start) / 60, 1))
    print("  输出目录:", OUT_DIR)
    print("======================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
