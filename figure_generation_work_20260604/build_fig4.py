# -*- coding: utf-8 -*-
"""Rebuild fig4: 测试实验汇总图 (test-experiment overview).
两段式：
  一、主结果对比 (scTARF 各版本 vs Baseline / scDRS-only)
       M1 各性状 AUROC 分组柱、M2 相对 baseline 的 Δ、M3 版本×性状热图
  二、方法学测试
       A 模块消融、B 参数权重扫描、C 计算效率
目标性状：paired = MONO/PLT/LYM；non-paired = MONO/RDW/LYM。
版本命名：scTARF x.x-关键词；做法留到流程图，此处只给名称与结果。
"""
import base64
import struct

DELIV = (r"D:/RStudio/work/multiomics_fusion_final_delivery_plan_20260528/"
         r"final_outputs")
FIG4ASSET = r"D:/RStudio/work/figure_generation_work_20260604/fig4_main_assets"
ASSETS = {
    # 主结果（由 scripts/02_plot_main_results.R 按现行标准重新生成）
    "M1": FIG4ASSET + "/fig_main_auroc_by_trait.png",
    "M2": FIG4ASSET + "/fig_main_delta_vs_baseline.png",
    "M3": FIG4ASSET + "/fig_main_method_heatmap.png",
    # 方法学测试
    "A": DELIV + "/03_ablation_parameter/figures/fig1_paired_ablation.png",
    "B": DELIV + "/03_ablation_parameter/figures/fig3_paired_weight_sweep.png",
    "C": DELIV + "/04_efficiency/fig_efficiency_final_overview.png",
}
OUT = (r"D:/RStudio/work/figure_generation_work_20260604/"
       r"drawio_xml/fig4_test_experiments_overview.drawio")


def png_size(path):
    d = open(path, "rb").read(33)
    return struct.unpack(">II", d[16:24])


def load_uri(path):
    with open(path, "rb") as f:
        return "data:image/png," + base64.b64encode(f.read()).decode("ascii")


URIS = {k: load_uri(p) for k, p in ASSETS.items()}
PXSIZE = {k: png_size(p) for k, p in ASSETS.items()}


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;"))


cells = []


def box(id_, value, style, x, y, w, h):
    cells.append(f'<mxCell id="{id_}" value="{esc(value)}" style="{style}" '
                 f'vertex="1" parent="1"><mxGeometry x="{x}" y="{y}" '
                 f'width="{w}" height="{h}" as="geometry"/></mxCell>')


def image(id_, uri, x, y, w, h):
    cells.append(f'<mxCell id="{id_}" value="" style="shape=image;html=1;'
                 f'imageAspect=0;aspect=fixed;image={uri};" vertex="1" '
                 f'parent="1"><mxGeometry x="{x}" y="{y}" width="{w}" '
                 f'height="{h}" as="geometry"/></mxCell>')


def text(id_, value, x, y, w, h, fs=14, bold=False, color="#222222",
         align="left", va="middle"):
    st = (f"text;html=1;strokeColor=none;fillColor=none;fontFamily=Arial;"
          f"fontSize={fs};fontColor={color};align={align};verticalAlign={va};"
          + ("fontStyle=1;" if bold else ""))
    cells.append(f'<mxCell id="{id_}" value="{esc(value)}" style="{st}" '
                 f'vertex="1" parent="1"><mxGeometry x="{x}" y="{y}" '
                 f'width="{w}" height="{h}" as="geometry"/></mxCell>')


def fit(key, bx, by, bw, bh):
    pw, ph = PXSIZE[key]
    r = pw / ph
    w = bw
    h = w / r
    if h > bh:
        h = bh
        w = h * r
    return round(bx + (bw - w) / 2), round(by + (bh - h) / 2), round(w), round(h)


# palette
WHITE, BORDER = "#FFFFFF", "#DDDDDD"
GREEN_S, BLUE_S, AMBER_S = "#59A14F", "#4E79A7", "#E8A33D"
PURPLE_S, TEAL_S, RED_S = "#8E63CE", "#4E9C8F", "#C44E4E"
B = "rounded=1;whiteSpace=wrap;html=1;fontFamily=Arial;"

W = 1700           # content width
X0 = 40

# ---- title ----
text("title", "Figure 4. scTARF 各版本主结果对比与方法学测试汇总",
     X0, 18, W, 40, fs=26, bold=True, color="#1A1A1A")

# ---- overview bar ----
box("ovf", "", B + f"fillColor=#F7F9FC;strokeColor={BORDER};arcSize=3;",
    X0, 70, W, 92)
text("ov1", "测试实验总览：上半区 = 主结果对比（scTARF 各版本 vs Baseline / "
     "scDRS-only），下半区 = 方法学测试（消融 / 参数 / 效率）",
     X0 + 18, 80, W - 36, 26, fs=15, bold=True, color="#2A2A2A")
text("ov2", "版本命名 scTARF x.x-关键词；Baseline = ATAC-only SCAVENGE，"
     "scDRS-only = RNA 基线。目标性状：paired = MONO / PLT / LYM；"
     "non-paired = MONO / RDW / LYM。各历史版本的具体做法见后续方法流程图。",
     X0 + 18, 110, W - 36, 44, fs=13, color="#444444", va="top")


def section(id_, y, label, color):
    box(id_, "", f"rounded=0;whiteSpace=wrap;html=1;fontFamily=Arial;"
        f"fillColor={color};strokeColor={color};", X0, y, W, 34)
    text(id_ + "t", label, X0 + 16, y + 3, W - 32, 28, fs=16, bold=True,
         color="#FFFFFF")


def panel(fid, key, x, y, w, h, letter, title, color, result):
    box(fid, "", B + f"fillColor={WHITE};strokeColor={BORDER};arcSize=2;",
        x, y, w, h)
    text("L" + fid, letter, x + 14, y + 9, 60, 28, fs=20, bold=True,
         color="#1A1A1A")
    text("T" + fid, title, x + 60, y + 10, w - 74, 26, fs=16, bold=True,
         color=color)
    ix, iy, iw, ih = fit(key, x + 16, y + 46, w - 32, h - 100)
    image("I" + fid, URIS[key], ix, iy, iw, ih)
    text("R" + fid, result, x + 16, y + h - 50, w - 32, 42, fs=12,
         color="#333333", align="center")


# ============ 一、主结果对比 ============
section("sec1", 174, "一、主结果对比：scTARF 各版本 vs Baseline / scDRS-only",
        "#3B6CA8")

# Row1: M1 (auroc) + M2 (delta)
HALF = 830
panel("pM1", "M1", X0, 218, HALF, 470, "M1",
      "各性状 AUROC 分组柱", PURPLE_S,
      "paired：scTARF1.0-paired 三性状全面超 Baseline 与 scDRS-only；"
      "non-paired：scTARF3.3-sparse 在 MONO/RDW/LYM 均 ≥ Baseline")
panel("pM2", "M2", X0 + HALF + 40, 218, HALF, 470, "M2",
      "相对 ATAC-only baseline 的 Δ", PURPLE_S,
      "paired 增益最大（MONO Δ≈+0.44）；non-paired 最终版 scTARF3.3-sparse "
      "三性状均为正向，中间版本在部分性状为负")

# Row2: M3 heatmap (full width)
panel("pM3", "M3", X0, 700, W, 612, "M3",
      "版本 × 性状 AUROC 热图", PURPLE_S,
      "横轴为各 scTARF 版本与两条参考基线，纵轴为 paired / non-paired 目标性状；"
      "灰格表示该版本未在对应性状上运行")

# ============ 二、方法学测试 ============
section("sec2", 1330, "二、方法学测试：模块消融 / 参数权重扫描 / 计算效率"
        "（仅展示配对正效果性状 MONO / PLT / LYM）", "#7A6B3A")

# Row3: A ablation + B weight
panel("pA", "A", X0, 1374, HALF, 600, "A",
      "模块消融：method design", GREEN_S,
      "BMMC paired：MONO 0.53→0.98、PLT 0.72→0.87、LYM 0.70→0.86"
      "（ATAC-only → full fusion）")
panel("pB", "B", X0 + HALF + 40, 1374, HALF, 600, "B",
      "参数权重扫描：weight sensitivity", BLUE_S,
      "atac_w 0.5–1.0 三性状全程高于各自 ATAC-only baseline（虚线）；"
      "曲线平缓，推荐 atac_w = 0.7")

# Row4: C efficiency (full width)
panel("pC", "C", X0, 1986, W, 700, "C",
      "计算效率：runtime / memory scaling", AMBER_S,
      "MONO 基准，1K→100K：runtime 0.5→25 s、peak memory ≤ 0.31 GB，"
      "AUROC/ACC 跨规模基本稳定；基于 subsampling / bootstrap scaling，非端到端实测")

# ---- caption ----
CAPY = 2710
box("capf", "", B + f"fillColor=#FBFBFD;strokeColor={BORDER};arcSize=2;",
    X0, CAPY, W, 250)
caption = (
    "图 4. scTARF 各版本主结果对比与方法学测试汇总。版本统一命名为 "
    "scTARF x.x-关键词，Baseline（ATAC-only SCAVENGE）与 scDRS-only（RNA）为两条"
    "参考基线；目标性状在 paired BMMC multiome 取 MONO、PLT、LYM，在 non-paired "
    "Granja ATAC + RNA reference 取 MONO、RDW、LYM。\n"
    "【主结果】(M1) 各性状 AUROC 分组柱：paired 场景 scTARF1.0-paired 在三性状均"
    "明显高于 Baseline 与 scDRS-only；non-paired 场景最终版 scTARF3.3-sparse 在 "
    "MONO/RDW/LYM 均不低于 Baseline。(M2) 相对 ATAC-only baseline 的 Δ AUROC："
    "paired 提升最大（MONO Δ≈+0.44、PLT/LYM Δ≈+0.15），non-paired 最终版三性状"
    "均为正向，部分中间版本（如 scTARF2.1-pseudopair）在个别性状为负，说明版本"
    "迭代的必要性。(M3) 版本 × 性状 AUROC 热图，灰格表示该版本未在对应性状上运行。\n"
    "【方法学测试】(A) 模块消融：full fusion 优于任一单模态与早期版本；(B) 参数"
    "权重扫描：atac_w 全程优于基线、推荐 w = 0.7；(C) 计算效率：runtime 随规模上升"
    "至约 25 s、peak memory ≤ 0.31 GB。各历史版本的具体做法将在方法流程图中展开，"
    "此处只给名称与结果。")
text("cap", caption, X0 + 18, CAPY + 14, W - 36, 222, fs=13, color="#222222",
     align="left", va="top")

# ---- assemble ----
PAGE_H = CAPY + 250 + 40
inner = "\n        ".join(cells)
xml = f'''<mxfile host="app.diagrams.net">
  <diagram id="fig4" name="Figure 4">
    <mxGraphModel dx="1400" dy="900" grid="0" gridSize="10" guides="1" \
tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" \
pageWidth="1780" pageHeight="{PAGE_H}" math="0" shadow="0">
      <root>
        <mxCell id="0"/>
        <mxCell id="1" parent="0"/>
        {inner}
      </root>
    </mxGraphModel>
  </diagram>
</mxfile>'''

with open(OUT, "w", encoding="utf-8") as f:
    f.write(xml)
print("wrote", OUT, "len", len(xml), "pageH", PAGE_H)
