"""Build the historical 50-component ATAC LSI consumed by non-paired workflows.

compute_lsi is copied unchanged from module1_method_B_build_h5ad.py.
Run 02_export_granja_lsi_counts.R first. This does not train scSHEFT.
"""
import os
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize

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

if __name__ == "__main__":
    project = Path(os.environ.get("SCTARF_PROJECT_ROOT", "."))
    data = Path(os.environ.get("SCTARF_DATA_ROOT", str(project / "local_data")))
    blood = Path(os.environ.get("SCTARF_BLOOD_ROOT", str(data / "Granja")))
    out = blood / "multiomics_module1_pairing/method_B_scSHEFT/scsheft_input"
    output = out / "blood_lsi.npy"
    if output.exists():
        raise SystemExit("blood_lsi.npy already exists; reuse it or choose a new data root.")
    counts = sparse.csr_matrix(mmread(out / "atac_counts.mtx").T)
    import csv
    with (out / "atac_cells.csv").open(encoding="utf-8-sig", newline="") as handle:
        cells = [row["cell"] for row in csv.DictReader(handle)]
    if len(cells) != counts.shape[0] or len(set(cells)) != len(cells):
        raise ValueError("ATAC count rows and unique cell barcodes must match")
    result = compute_lsi(counts, n_components=50, random_state=9527)
    np.save(output, result)
    print("Wrote", output, result.shape)
