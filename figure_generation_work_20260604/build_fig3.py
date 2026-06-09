# -*- coding: utf-8 -*-
"""Rebuild fig3: Panel A vector flowchart (label-free two-gate process),
Panel B embedded suitability bar chart, Panel C native decision table.
All content driven by 01_data_prep_rna_reference_screening.
Value text uses raw chars (<br>, <b>, >, <) and is entity-escaped on emit."""
import base64

# bar chart lives with the screening tool's outputs; drawio is written to drawio_xml/
ASSET = (r"D:/RStudio/work/module4_traitaware_pseudopair/"
         r"01_data_prep_rna_reference_screening/output/"
         r"rna_suitability_scores_barplot.png")
OUT = (r"D:/RStudio/work/figure_generation_work_20260604/"
       r"drawio_xml/fig3_rna_reference_suitability.drawio")

with open(ASSET, "rb") as f:
    img_uri = "data:image/png," + base64.b64encode(f.read()).decode("ascii")

def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;"))

cells = []
def box(id_, value, style, x, y, w, h):
    cells.append(f'<mxCell id="{id_}" value="{esc(value)}" style="{style}" vertex="1" parent="1">'
                 f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>')
def image(id_, uri, x, y, w, h):
    cells.append(f'<mxCell id="{id_}" value="" style="shape=image;html=1;imageAspect=0;aspect=fixed;'
                 f'image={uri};" vertex="1" parent="1">'
                 f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>')
def edge(id_, src, tgt, style="endArrow=block;html=1;strokeColor=#555555;rounded=0;"):
    cells.append(f'<mxCell id="{id_}" style="{style}" edge="1" parent="1" source="{src}" target="{tgt}">'
                 f'<mxGeometry relative="1" as="geometry"/></mxCell>')
def text(id_, value, x, y, w, h, fs=14, bold=False, color="#222222", align="left", va="middle"):
    st = (f"text;html=1;strokeColor=none;fillColor=none;fontFamily=Arial;fontSize={fs};"
          f"fontColor={color};align={align};verticalAlign={va};" + ("fontStyle=1;" if bold else ""))
    cells.append(f'<mxCell id="{id_}" value="{esc(value)}" style="{st}" vertex="1" parent="1">'
                 f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>')

# palette
GREEN_F, GREEN_S = "#E4F4E2", "#59A14F"
AMBER_F, AMBER_S = "#FFF1D6", "#E8A33D"
RED_F, RED_S = "#F8E0E0", "#C44E4E"
GRAY_F, GRAY_S = "#F2F4F7", "#9AA3AD"
BLUE_F, BLUE_S = "#E8F0FB", "#4E79A7"
WHITE, BORDER = "#FFFFFF", "#DDDDDD"
B = "rounded=1;whiteSpace=wrap;html=1;fontFamily=Arial;"

# ---- title ----
text("title", "Figure 3. 非配对 scTARF 的 RNA reference 适合性筛选 (label-free)",
     40, 22, 1200, 40, fs=26, bold=True, color="#1A1A1A")

# ---- panel frames ----
box("pa", "", B + f"fillColor={WHITE};strokeColor={BORDER};arcSize=2;shadow=0;", 40, 78, 1480, 322)
box("pb", "", B + f"fillColor={WHITE};strokeColor={BORDER};arcSize=2;shadow=0;", 40, 416, 776, 454)
box("pc", "", B + f"fillColor={WHITE};strokeColor={BORDER};arcSize=2;shadow=0;", 836, 416, 684, 454)

# ============ PANEL A : flowchart ============
text("la", "A", 54, 90, 30, 28, fs=20, bold=True, color="#1A1A1A")
text("ta", "适合性门控流程：两道 label-free gate → 打分 → 判定 → 选择",
     90, 91, 740, 26, fs=16, bold=True, color="#2A2A2A")

y, h = 138, 66
box("a_in", "输入<br><b>RNA reference 基因</b><br>+ ATAC gene-activity 基因<br>+ trait gene set",
    B + f"fillColor={GREEN_F};strokeColor={GREEN_S};fontSize=12;", 62, y, 188, h)
box("a_g1", "<b>Gate 1 · 基因覆盖</b><br>gene_coverage_score",
    B + f"fillColor={AMBER_F};strokeColor={AMBER_S};fontSize=13;", 286, y, 178, h)
box("a_g2", "<b>Gate 2 · RNA 内在信号</b><br>rna_signal_score (scDRS)",
    B + f"fillColor={AMBER_F};strokeColor={AMBER_S};fontSize=13;", 500, y, 184, h)
box("a_sc", "<b>suitability score</b><br>0.45 · coverage<br>+ 0.55 · signal",
    B + f"fillColor={BLUE_F};strokeColor={BLUE_S};fontSize=12;", 720, y, 168, h)
box("a_dec", "<b>final_gate</b><br>判定",
    B + f"fillColor={GRAY_F};strokeColor={GRAY_S};fontSize=13;", 924, y, 96, h)
edge("ae1", "a_in", "a_g1"); edge("ae2", "a_g1", "a_g2")
edge("ae3", "a_g2", "a_sc"); edge("ae4", "a_sc", "a_dec")

ann = ("rounded=1;whiteSpace=wrap;html=1;fillColor=#FFFDF2;dashed=1;"
       f"strokeColor={AMBER_S};fontFamily=Arial;fontSize=10;align=left;"
       "spacingLeft=6;spacingRight=4;verticalAlign=top;spacingTop=4;")
box("c1", "<b>阈值</b><br>PASS: rna_w>=0.70 且 shared_w>=0.65 且 shared基因>=min(500,0.60·n)<br>"
    "RISK: rna_w>=0.55 且 shared_w>=0.50 且 shared基因>=min(300,0.40·n)",
    ann, 280, 218, 200, 92)
box("c2", "<b>阈值</b><br>PASS: p<0.01 富集>=3.0 且 frac(p<=0.05)>=0.08 且 n(FDR<=0.10,+)>=max(20,0.001·n)<br>"
    "RISK: 富集>=1.8 且 frac>=0.06 且 n>=1；无 score = NO_SCORE_FOR_TRAIT",
    ann, 494, 218, 230, 92)
edge("c1e", "a_g1", "c1", "endArrow=none;html=1;strokeColor=#C9B26A;dashed=1;")
edge("c2e", "a_g2", "c2", "endArrow=none;html=1;strokeColor=#C9B26A;dashed=1;")

box("o_pass", "PASS_PRIOR<br><b>→ USE_FOR_FUSION_CANDIDATE</b>",
    B + f"fillColor={GREEN_F};strokeColor={GREEN_S};fontSize=12;", 1064, 110, 372, 50)
box("o_risk", "RISK_PRIOR<br><b>→ REVIEW_BEFORE_FUSION</b>",
    B + f"fillColor={AMBER_F};strokeColor={AMBER_S};fontSize=12;", 1064, 168, 372, 50)
box("o_fail", "FAIL_* / NO_SCORE_FOR_TRAIT<br><b>→ DO_NOT_USE_FOR_FUSION</b>",
    B + f"fillColor={RED_F};strokeColor={RED_S};fontSize=12;", 1064, 226, 372, 50)
edge("oe1", "a_dec", "o_pass"); edge("oe2", "a_dec", "o_risk"); edge("oe3", "a_dec", "o_fail")

box("a_sel", "<b>选择规则</b>　每个 trait 在 Granja RNA 与 CITE RNA 两套 reference 中，取门控合格（PASS_PRIOR 或 RISK_PRIOR）"
    "且 suitability score 最高者作为 recommended_reference；若两套都 FAIL，则 recommended_reference = none（PLT 即此情况）。"
    "筛选只用基因覆盖与 scDRS 分布，不使用 celltype / annotation / mapping / AUROC。",
    B + f"fillColor=#F7F9FB;strokeColor={BORDER};fontSize=11;align=left;spacingLeft=10;spacingRight=10;",
    62, 320, 1374, 64)

# ============ PANEL B : bar chart ============
text("lb", "B", 54, 426, 30, 28, fs=20, bold=True, color="#1A1A1A")
text("tb", "两套 RNA reference 的 suitability score 对比（按门控着色，★ = 推荐）",
     90, 427, 720, 26, fs=16, bold=True, color="#2A2A2A")
image("img_bar", img_uri, 70, 455, 683, 400)

# ============ PANEL C : decision table ============
text("lc", "C", 850, 426, 30, 28, fs=20, bold=True, color="#1A1A1A")
text("tc", "每个 trait 的推荐 reference 与判定", 886, 427, 500, 26, fs=16, bold=True, color="#2A2A2A")

tx, ty = 850, 470
cw = [62, 116, 150, 150, 178]   # sum = 656
rh_head, rh = 38, 38
cols = ["trait", "推荐 reference", "Granja 分(门控)", "CITE 分(门控)", "判定"]

def tcell(id_, txt, cx, cy, w, hh, fill, stroke="#CCCCCC", bold=False, fs=11, fcolor="#222222"):
    st = (f"rounded=0;whiteSpace=wrap;html=1;fontFamily=Arial;fontSize={fs};fillColor={fill};"
          f"strokeColor={stroke};fontColor={fcolor};align=center;verticalAlign=middle;"
          + ("fontStyle=1;" if bold else ""))
    cells.append(f'<mxCell id="{id_}" value="{esc(txt)}" style="{st}" vertex="1" parent="1">'
                 f'<mxGeometry x="{cx}" y="{cy}" width="{w}" height="{hh}" as="geometry"/></mxCell>')

cx = tx
for j, c in enumerate(cols):
    tcell(f"th{j}", c, cx, ty, cw[j], rh_head, "#3B4A5A", stroke="#3B4A5A", bold=True, fcolor="#FFFFFF")
    cx += cw[j]

PASS, RISK, FAIL, NONE = (GREEN_F, GREEN_S), (AMBER_F, AMBER_S), (RED_F, RED_S), ("#F2F4F7", "#9AA3AD")
def gcell(score, gate):
    if gate.startswith("PASS"): return f"{score} (PASS)", PASS
    if gate.startswith("RISK"): return f"{score} (RISK)", RISK
    if gate.startswith("FAIL"): return f"{score} (FAIL)", FAIL
    return "— (无分)", NONE
DEC = {"USE": ("USE_FOR_FUSION", GREEN_F, GREEN_S),
       "REVIEW": ("REVIEW_BEFORE_FUSION", AMBER_F, AMBER_S),
       "NO": ("DO_NOT_USE", RED_F, RED_S)}

rows = [
    ("MONO", "granja_rna", ("88.65","PASS"), ("79.77","PASS"), "USE"),
    ("LYM",  "granja_rna", ("88.10","PASS"), ("50.23","RISK"), "USE"),
    ("MCH",  "cite_rna",   (None,"NO_SCORE"), ("87.67","PASS"), "USE"),
    ("RDW",  "cite_rna",   (None,"NO_SCORE"), ("89.58","PASS"), "USE"),
    ("RBC",  "cite_rna",   (None,"NO_SCORE"), ("83.01","PASS"), "USE"),
    ("EOS",  "cite_rna",   (None,"NO_SCORE"), ("46.51","RISK"), "REVIEW"),
    ("WBC",  "cite_rna",   (None,"NO_SCORE"), ("49.27","RISK"), "REVIEW"),
    ("PLT",  "none",       ("51.69","FAIL"), ("45.84","FAIL"),  "NO"),
]
cy = ty + rh_head
for i, (tr, rec, gj, ct, dec) in enumerate(rows):
    r = f"r{i}"
    tcell(f"{r}c0", tr, tx, cy, cw[0], rh, "#FBFCFD", bold=True, fs=12)
    tcell(f"{r}c1", rec, tx+cw[0], cy, cw[1], rh,
          ("#FFFFFF" if rec == "none" else "#EFF6EE"), bold=(rec != "none"),
          fcolor=("#999999" if rec == "none" else "#1A1A1A"))
    gtxt, (gf, gs) = gcell(gj[0], gj[1])
    tcell(f"{r}c2", gtxt, tx+cw[0]+cw[1], cy, cw[2], rh, gf, stroke=gs)
    ctxt, (cf, cs) = gcell(ct[0], ct[1])
    tcell(f"{r}c3", ctxt, tx+cw[0]+cw[1]+cw[2], cy, cw[3], rh, cf, stroke=cs)
    dtxt, df, ds = DEC[dec]
    tcell(f"{r}c4", dtxt, tx+cw[0]+cw[1]+cw[2]+cw[3], cy, cw[4], rh, df, stroke=ds, bold=True, fs=10)
    cy += rh

text("tc_fn", "无分 = 该 reference 无此 trait 的 scDRS score 文件，不能作为候选。",
     tx, cy+6, 640, 22, fs=10, color="#777777")

# ---- caption ----
cap = ("图 3. RNA reference 选择对非配对 scTARF 的影响（label-free 先验筛选）。"
       "(A) 适合性门控流程：Gate 1 用 trait 基因在 RNA reference 与 ATAC gene-activity 上的覆盖给出 gene_coverage_score，"
       "Gate 2 用已有 RNA scDRS 分数的 p 值富集、FDR 阳性细胞数与离散度给出 rna_signal_score，"
       "二者按 0.45 / 0.55 加权得到 0–100 的 suitability score，并据门控判为 PASS / RISK / FAIL，再映射到融合候选判定。"
       "(B) 八个性状在 Granja RNA 与 CITE RNA 两套 reference 上的 suitability score，颜色编码门控状态，黑边加 ★ 为推荐 reference；"
       "Granja 在 MCH / RDW / RBC / EOS / WBC 上无 scDRS 分，故不参与这些性状的候选。"
       "(C) 各性状最终推荐的 reference 与判定。"
       "MONO / LYM 选 Granja RNA，MCH / RDW / RBC 选 CITE RNA 且通过门控，EOS / WBC 仅达 RISK 需复核，"
       "PLT 在两套 reference 下均 FAIL（RNA 内在信号过弱），不进入融合，须依赖 ATAC 端。")
text("caption", cap, 40, 884, 1480, 170, fs=14, color="#222222", align="left", va="top")

# ---- assemble ----
body = "\n        ".join(cells)
xml = ('<mxfile host="app.diagrams.net" modified="2026-06-05T00:00:00.000Z" agent="Claude" version="24.0.0" type="device">\n'
       '  <diagram id="fig3" name="Figure 3 - RNA reference suitability">\n'
       '    <mxGraphModel dx="1600" dy="1000" grid="1" gridSize="10" guides="1" tooltips="1" connect="1" '
       'arrows="1" fold="1" page="1" pageScale="1" pageWidth="1600" pageHeight="1080" math="0" shadow="0">\n'
       '      <root>\n        <mxCell id="0"/>\n        <mxCell id="1" parent="0"/>\n'
       f'        {body}\n      </root>\n    </mxGraphModel>\n  </diagram>\n</mxfile>\n')
with open(OUT, "w", encoding="utf-8") as f:
    f.write(xml)
print("wrote", OUT, "bytes", len(xml))
