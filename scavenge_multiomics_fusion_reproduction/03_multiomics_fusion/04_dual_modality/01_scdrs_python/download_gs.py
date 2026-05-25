"""
download_gs.py
拉 SCDRS 官方 UKB 460K MAGMA .gs 文件, 抽 MONO/PLT/LYM 三 trait.

用法 (在 scdrs conda env 里):
    conda activate scdrs
    cd D:/RStudio/work/module3_multiome/04_dual_modality/01_scdrs_python
    python download_gs.py

输出 (同目录下):
    magma_10kb_top1000_zscore.74_traits.rv1.gs   原始 74 trait
    bmmc_3traits.gs                              抽取后的三 trait
    scdrs_top_genes.csv                          关卡 1 用 (top 20 gene per trait)
"""
import os
import sys
import zipfile
import requests
import pandas as pd

ARTICLE_ID = "19312583"
ZIP_FILENAME = "gs_file.zip"
TARGET_GS = "magma_10kb_top1000_zscore.74_traits.rv1.gs"
TARGET_TRAITS = [
    "UKB_460K.blood_MONOCYTE_COUNT",
    "UKB_460K.blood_PLATELET_COUNT",
    "UKB_460K.blood_LYMPHOCYTE_COUNT",
]
OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    # ---- 1. figshare API 查文件 ----
    print(f"[1/4] querying figshare article {ARTICLE_ID} ...")
    api = f"https://api.figshare.com/MultiOmicsFusion/articles/{ARTICLE_ID}/files"
    resp = requests.get(api, timeout=30)
    resp.raise_for_status()
    files = resp.json()
    match = [f for f in files if f["name"] == ZIP_FILENAME]
    if not match:
        available = sorted({f["name"] for f in files})
        print(f"FAIL: '{ZIP_FILENAME}' not in article. Available files:")
        for n in available:
            print(f"  {n}")
        sys.exit(1)
    finfo = match[0]
    print(f"  found: {finfo['name']} ({finfo['size']/1e6:.1f} MB)")

    # ---- 2. 下载 + 解压 ----
    zip_path = os.path.join(OUT_DIR, ZIP_FILENAME)
    if os.path.exists(zip_path) and os.path.getsize(zip_path) == finfo["size"]:
        print(f"[2/4] zip already downloaded, skip: {zip_path}")
    else:
        print(f"[2/4] downloading from {finfo['download_url']} ...")
        with requests.get(finfo["download_url"], stream=True, timeout=600) as r:
            r.raise_for_status()
            total = int(r.headers.get("Content-Length", 0))
            done = 0
            with open(zip_path, "wb") as fp:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    fp.write(chunk)
                    done += len(chunk)
                    if total:
                        pct = done / total * 100
                        sys.stdout.write(
                            f"\r  {done/1e6:6.1f} / {total/1e6:.1f} MB ({pct:5.1f}%)"
                        )
                        sys.stdout.flush()
                print()
        print(f"  saved: {zip_path}")

    # 在 zip 里找目标 .gs (不全部解压, 只抽一份)
    print(f"  extracting '{TARGET_GS}' from zip ...")
    raw_path = os.path.join(OUT_DIR, TARGET_GS)
    with zipfile.ZipFile(zip_path) as zf:
        members = zf.namelist()
        match_in_zip = [m for m in members if m.endswith(TARGET_GS)]
        if not match_in_zip:
            print(f"FAIL: '{TARGET_GS}' not in zip. Zip contents (first 20):")
            for m in members[:20]:
                print(f"  {m}")
            print(f"  ... total {len(members)} entries")
            sys.exit(1)
        member = match_in_zip[0]
        with zf.open(member) as src, open(raw_path, "wb") as dst:
            dst.write(src.read())
        print(f"  extracted: {raw_path}  (from zip entry: {member})")

    # ---- 3. 抽 3 trait ----
    print(f"[3/4] reading {raw_path} ...")
    df = pd.read_csv(raw_path, sep="\t")
    print(f"  total {len(df)} traits; columns: {list(df.columns)}")

    have = set(df["TRAIT"])
    missing = [t for t in TARGET_TRAITS if t not in have]
    if missing:
        print(f"FAIL: traits not found: {missing}")
        print("hint: traits with 'blood_' in name:")
        for t in sorted(df["TRAIT"]):
            if "blood_" in t:
                print(f"  {t}")
        sys.exit(1)

    sub = df[df["TRAIT"].isin(TARGET_TRAITS)].copy()
    out_gs = os.path.join(OUT_DIR, "bmmc_3traits.gs")
    sub.to_csv(out_gs, sep="\t", index=False)
    print(f"  saved 3-trait gs: {out_gs}")

    # ---- 4. top 20 genes 摘要 (关卡 1) ----
    print("[4/4] extracting top 20 genes per trait ...")
    rows = []
    for _, r in sub.iterrows():
        trait = r["TRAIT"]
        for rank, gw in enumerate(r["GENESET"].split(","), 1):
            if rank > 20:
                break
            parts = gw.split(":")
            rows.append(
                {
                    "TRAIT": trait,
                    "rank": rank,
                    "GENE": parts[0],
                    "WEIGHT": float(parts[1]) if len(parts) == 2 else None,
                }
            )
    top = pd.DataFrame(rows)
    out_top = os.path.join(OUT_DIR, "scdrs_top_genes.csv")
    top.to_csv(out_top, index=False)
    print(f"  saved: {out_top}")
    print()

    # 屏幕打印 — 关卡 1 直接看这个
    print("=" * 60)
    print("CHECKPOINT 1 INPUT — top 20 genes per trait")
    print("=" * 60)
    for trait in TARGET_TRAITS:
        sub_t = top[top["TRAIT"] == trait]
        print(f"\n[{trait}]")
        print("  " + ", ".join(sub_t["GENE"].tolist()))
    print()


if __name__ == "__main__":
    main()
