# -*- coding: utf-8 -*-
"""non-paired scTARF 历史版本 vs 最终方法 (scTARF3.3-sparse) 流程差异图。
单色主轴 = 最终方法 v3.3 流程；彩色分支 = 各历史版本在对应步骤上的不同做法。
格式参照 fig1x_overview_nonpaired_three_routes.drawio (圆角框 / callout / block 箭头 / Arial)。
"""
OUT = (r"D:/RStudio/work/figure_generation_work_20260604/"
       r"drawio_xml/fig1x_overview_nonpaired_three_routes.drawio")


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;").replace("\n", "&#10;"))


cells = []


def node(id_, value, x, y, w, h, fill, stroke, fs=15, bold=False, rounded=True):
    shape = "rounded=1;" if rounded else "rounded=0;"
    st = (f"{shape}whiteSpace=wrap;html=1;fillColor={fill};strokeColor={stroke};"
          f"fontFamily=Arial;fontSize={fs};fontColor=#1A1A1A;"
          + ("fontStyle=1;" if bold else ""))
    cells.append(f'<mxCell id="{id_}" value="{esc(value)}" style="{st}" '
                 f'vertex="1" parent="1"><mxGeometry x="{x}" y="{y}" '
                 f'width="{w}" height="{h}" as="geometry"/></mxCell>')


def callout(id_, value, x, y, w, h):
    st = ("shape=callout;whiteSpace=wrap;html=1;perimeter=calloutPerimeter;"
          "fillColor=#FFF4D9;strokeColor=#C49A00;fontFamily=Arial;fontSize=13;")
    cells.append(f'<mxCell id="{id_}" value="{esc(value)}" style="{st}" '
                 f'vertex="1" parent="1"><mxGeometry x="{x}" y="{y}" '
                 f'width="{w}" height="{h}" as="geometry"/></mxCell>')


def text(id_, value, x, y, w, h, fs=14, bold=False, color="#222222",
         align="left", va="middle"):
    st = (f"text;html=1;strokeColor=none;fillColor=none;fontFamily=Arial;"
          f"fontSize={fs};fontColor={color};align={align};verticalAlign={va};"
          + ("fontStyle=1;" if bold else ""))
    cells.append(f'<mxCell id="{id_}" value="{esc(value)}" style="{st}" '
                 f'vertex="1" parent="1"><mxGeometry x="{x}" y="{y}" '
                 f'width="{w}" height="{h}" as="geometry"/></mxCell>')


def edge(id_, src, tgt, color="#555555", dashed=False, label=""):
    st = (f"endArrow=block;html=1;strokeColor={color};rounded=0;"
          f"strokeWidth={'2' if dashed else '1.5'};"
          + ("dashed=1;" if dashed else ""))
    lab = f' value="{esc(label)}"' if label else ' value=""'
    cells.append(f'<mxCell id="{id_}"{lab} style="{st}" edge="1" parent="1" '
                 f'source="{src}" target="{tgt}"><mxGeometry relative="1" '
                 f'as="geometry"/></mxCell>')


# ---- palette ----
MAIN_F, MAIN_S = "#ECEFF1", "#455A64"      # 主轴 slate (最终 v3.3)
ACC_F, ACC_S = "#FBE3E3", "#C44E4E"        # v3.3 关键步骤强调 (红)
IN_F, IN_S = "#F4F6F8", "#9AA3AD"          # 输入/输出
TEAL_F, TEAL_S = "#DDF0EC", "#4E9C8F"      # v2.1
BLUE_F, BLUE_S = "#DCEBFA", "#4C78A8"      # v2.2
LBLU_F, LBLU_S = "#E6F2FB", "#6FB1E0"      # v3.1
PINK_F, PINK_S = "#F7E3EF", "#D17DB0"      # v3.2

# main column geometry
CX, CW = 600, 380          # main boxes x / width
LX, LW = 60, 470           # left branch column
RX, RW = 1080, 480         # right branch column

# ---- title + note ----
text("title", "Figure 1x. non-paired scTARF 历史版本与最终方法 "
     "(scTARF3.3-sparse) 的流程差异",
     40, 26, 1500, 36, fs=26, bold=True, color="#1A1A1A")
text("note", "单色主轴 = 最终方法 scTARF3.3-sparse 的流程；彩色分支 = 各历史版本"
     "在对应步骤上的不同做法（只标差异，不重复主轴步骤）。",
     40, 66, 1500, 24, fs=15, color="#555555")

# ---- legend (top-right) ----
leg = [("最终 scTARF3.3-sparse（主轴）", MAIN_F, MAIN_S),
       ("scTARF3.1-top1", LBLU_F, LBLU_S),
       ("scTARF2.1-pseudopair", TEAL_F, TEAL_S),
       ("scTARF2.2-nbFDR", BLUE_F, BLUE_S),
       ("scTARF3.2-top5sig", PINK_F, PINK_S)]
node("legbox", "", 1180, 100, 470, 168, "#FFFFFF", "#DDDDDD", rounded=True)
text("legt", "图例", 1196, 106, 200, 22, fs=14, bold=True, color="#333333")
for i, (lab, f, s) in enumerate(leg):
    yy = 134 + i * 28
    node(f"lg{i}", "", 1196, yy, 26, 18, f, s)
    text(f"lt{i}", lab, 1232, yy - 3, 410, 24, fs=13, color="#333333")

# ============ MAIN AXIS (v3.3) ============
m = [
    ("m0", "输入：Granja ATAC query · RNA reference · trait-aware 基因集", 118, 54,
     IN_F, IN_S, 14, False),
    ("m1", "scDRS 在 RNA reference 上打分 → 每个 RNA 细胞 trait score", 200, 58,
     MAIN_F, MAIN_S, 15, True),
    ("m2", "ATAC → RNA 匹配：每个 ATAC 细胞取 top-5 RNA 邻居", 292, 58,
     MAIN_F, MAIN_S, 15, True),
    ("m3", "Top-5 显著性补全 + sparse 映射\n(random-top5 背景，empirical p)", 384, 76,
     ACC_F, ACC_S, 15, True),
    ("m4", "ATAC / RNA 加权融合向量 → fused KNN", 488, 58,
     MAIN_F, MAIN_S, 15, True),
    ("m5", "Neighborhood-FDR 种子选择", 580, 58,
     ACC_F, ACC_S, 15, True),
    ("m6", "SCAVENGE / RWR → 每个细胞 TRS", 672, 58,
     MAIN_F, MAIN_S, 15, True),
    ("m7", "按 target cell type 评估 AUROC / ACC / F1", 764, 58,
     IN_F, IN_S, 14, False),
]
for id_, val, y, h, f, s, fs, bold in m:
    node(id_, val, CX, y, CW, h, f, s, fs=fs, bold=bold)
for a, b in [("m0", "m1"), ("m1", "m2"), ("m2", "m3"), ("m3", "m4"),
             ("m4", "m5"), ("m5", "m6"), ("m6", "m7")]:
    edge("e_" + a + b, a, b)

# main-axis defining callouts
callout("ca_top5", "top-5 邻居", 990, 286, 120, 44)
callout("ca_sparse", "sparse 稀疏特征\n逐基因背景过滤", 990, 372, 170, 60)

# 主轴竖向标注
text("axislab", "最终方法主轴", CX - 150, 470, 130, 40, fs=14, bold=True,
     color=MAIN_S, align="right")

# ============ LEFT BRANCHES (匹配 / 补全 / 映射密度 的历史做法) ============
text("lhdr", "RNA 匹配与补全：历史做法", LX, 250, LW, 24, fs=15, bold=True,
     color="#444444")
# v3.1-top1  (匹配粒度) —— 对齐 m2
node("bv31", "scTARF3.1-top1\n每个 ATAC 仅取 top-1 RNA 邻居（top5 → top1 消融，"
     "无多参考平均/显著性）", LX, 292, LW, 58, LBLU_F, LBLU_S, fs=13)
edge("eb31", "bv31", "m2", color=LBLU_S, dashed=True)
# v2.1-pseudopair (dense 补全 + early 种子)
node("bv21", "scTARF2.1-pseudopair\ndense 伪配对：转移全部映射 RNA，无 top5 显著性"
     "筛选；种子用 early/basic 规则", LX, 372, LW, 62, TEAL_F, TEAL_S, fs=13)
edge("eb21", "bv21", "m3", color=TEAL_S, dashed=True)
# v2.2-nbFDR (dense 补全, 种子已升级)
node("bv22", "scTARF2.2-nbFDR\ndense completion（同 v2.1 补全），仅把种子升级为 "
     "neighborhood-FDR", LX, 452, LW, 58, BLUE_F, BLUE_S, fs=13)
edge("eb22", "bv22", "m3", color=BLUE_S, dashed=True)
# v3.2-top5sig (dense 映射, 无 sparse)
node("bv32", "scTARF3.2-top5sig\ntop5 显著性补全，但用 dense 映射（无 sparse 稀疏化），"
     "其余同 v3.3", LX, 528, LW, 58, PINK_F, PINK_S, fs=13)
edge("eb32", "bv32", "m3", color=PINK_S, dashed=True)

# ============ RIGHT BRANCH (种子选择 的历史做法) ============
text("rhdr", "种子选择：历史做法", RX, 540, RW, 24, fs=15, bold=True,
     color="#444444")
node("rv21", "scTARF2.1-pseudopair\nearly/basic 种子（scDRS 阈值法），"
     "非 neighborhood-FDR", RX, 580, RW, 58, TEAL_F, TEAL_S, fs=13)
edge("erv21", "rv21", "m5", color=TEAL_S, dashed=True)

# ---- bottom one-line takeaways ----
text("foot", "迭代主线：v2.1（dense + early 种子，最朴素）→ v2.2（升级 "
     "neighborhood-FDR 种子）→ v3.1/v3.2（top-k 与显著性补全消融）→ "
     "v3.3（top5 显著性 + sparse 映射，最终方法）。各版本只在上述彩色步骤与主轴不同，"
     "其余步骤完全一致。",
     60, 858, 1590, 56, fs=14, color="#333333", va="top")

# ---- assemble ----
inner = "\n        ".join(cells)
xml = f'''<mxfile host="app.diagrams.net" agent="Claude" version="24.0.0" type="device">
  <diagram id="fig1x_version_diff" name="Figure 1x - 版本流程差异">
    <mxGraphModel dx="1600" dy="1000" grid="1" gridSize="10" guides="1" \
tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" \
pageWidth="1700" pageHeight="950" math="0" shadow="0">
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
print("wrote", OUT, "len", len(xml))
