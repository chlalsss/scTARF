# Prepare final output folders and method-flow source files.
# Run after the final analysis outputs are available.

source(file.path("scripts", "07_figures_tables", "00_final_delivery_common.R"))
ensure_final_dirs()

out_dir <- file.path(FINAL_ROOT, "00_method_flow")

steps <- data.frame(
  branch = c(
    rep("paired_multiome", 7),
    rep("non_paired_reference_assisted", 11),
    rep("shared_scavenge_rwr", 4)
  ),
  step_id = c(
    "P1", "P2", "P3", "P4", "P5", "P6", "P7",
    "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10", "N11",
    "S1", "S2", "S3", "S4"
  ),
  step_name = c(
    "Select paired BMMC multiome cells",
    "Build ATAC LSI and RNA latent/scDRS signals",
    "Construct paired ATAC/RNA feature blocks",
    "Build weighted fused KNN",
    "Select trait-aware seeds",
    "Run RWR/TRS",
    "Evaluate against paired cell-type targets",
    "Select non-paired ATAC query and RNA reference",
    "Run scDRS on RNA reference",
    "Build ATAC gene activity and ATAC z block",
    "Find top-k RNA reference cells for each ATAC cell",
    "Average top5 RNA reference profiles",
    "Filter RNA-completed genes by random top5 background",
    "Construct ATAC/RNA fused representation",
    "Build fused KNN",
    "Select neighborhood-FDR trait seeds",
    "Run RWR/TRS",
    "Evaluate against post-hoc target cell types",
    "Collect per-cell TRS",
    "Compute AUROC/ACC/F1/AUPRC",
    "Summarize ablations and reference dependence",
    "Prepare case-study interpretation"
  ),
  current_source = c(
    rep("module3_multiome/ and ppt_dataset_figures/", 7),
    rep("final non-paired sparse-mapping workflow", 11),
    rep("final delivery scripts", 4)
  ),
  stringsAsFactors = FALSE
)

write_csv(steps, file.path(out_dir, "method_steps_table.csv"))
version_defs <- version_definition_table()
write_csv(version_defs, file.path(out_dir, "method_version_definitions.csv"))

version_md <- c(
  "# Method Version Definitions",
  "",
  "Final reporting uses method-level version names. If the same algorithm was run on different traits in separate batches, those runs are merged into the same version.",
  "",
  "| Version | Experiment | Raw results | Data setting | Traits | Main method |",
  "|---|---|---|---|---|---|",
  apply(version_defs, 1, function(r) {
    paste0(
      "| ", r[["version"]],
      " | ", r[["experiment"]],
      " | `", r[["raw_results"]],
      "` | ", r[["data_setting"]],
      " | ", r[["traits"]],
      " | ", r[["main_method"]],
      " |"
    )
  }),
  "",
  "The final non-paired runs use the same sparse-mapping workflow across reported traits: RNA reference evidence is mapped onto ATAC cells, significant completed RNA features are fused with ATAC features, trait-aware seeds are selected, and the shared SCAVENGE/RWR/TRS code is reused for scoring."
)
writeLines(version_md, file.path(out_dir, "method_version_definitions.md"), useBytes = TRUE)

mermaid <- c(
  "flowchart TD",
  "  A1[Paired BMMC multiome] --> A2[ATAC LSI + RNA trait signals]",
  "  A2 --> A3[Paired fused feature blocks]",
  "  A3 --> A4[Weighted fused KNN]",
  "  B1[Non-paired ATAC query] --> B3[ATAC gene activity / z block]",
  "  B2[scRNA reference] --> B4[scDRS RNA trait scores]",
  "  B3 --> B5[Top-k ATAC-to-RNA matching]",
  "  B4 --> B5",
  "  B5 --> B6[Top5 RNA completion]",
  "  B6 --> B7[Random-top5 significance filter]",
  "  B7 --> B8[ATAC/RNA fused representation]",
  "  B8 --> B9[Fused KNN]",
  "  A4 --> C1[Trait-aware seed selection]",
  "  B9 --> C1",
  "  C1 --> C2[SCAVENGE/RWR propagation]",
  "  C2 --> C3[Cell-level TRS]",
  "  C3 --> C4[AUROC / ACC / F1 / case study]"
)
writeLines(mermaid, file.path(out_dir, "method_flowchart.mmd"), useBytes = TRUE)

protocol <- c(
  "# Final Method Protocol",
  "",
  "Method name: Trait-aware RNA-ATAC fusion SCAVENGE.",
  "",
  "Paired branch: use paired BMMC multiome data to build ATAC and RNA feature blocks for the same cells, construct a fused KNN, select trait-aware seeds, and run SCAVENGE/RWR to produce TRS.",
  "",
  "Non-paired branch: use scATAC as query and scRNA as reference. RNA reference cells are scored by scDRS. Each ATAC cell is matched to top-k RNA reference cells, top5 RNA profiles are averaged, and RNA-completed signals are filtered against random top5 RNA backgrounds before fusion.",
  "",
  "Shared downstream step: the fused graph and selected trait seeds enter the existing SCAVENGE/RWR/TRS workflow. Evaluation uses post-hoc target cell-type definitions and reports AUROC, ACC, F1, and related metrics where available.",
  "",
  "Source modules: paired analyses reuse the paired multiome workflow; non-paired analyses use the final sparse-mapping workflow."
)
writeLines(protocol, file.path(out_dir, "method_protocol.md"), useBytes = TRUE)

append_manifest(data.frame(
  final_file = c(
    file.path(out_dir, "method_steps_table.csv"),
    file.path(out_dir, "method_version_definitions.csv"),
    file.path(out_dir, "method_version_definitions.md"),
    file.path(out_dir, "method_flowchart.mmd"),
    file.path(out_dir, "method_protocol.md")
  ),
  source_file = c(
    "module3_multiome; module4_traitaware_pseudopair",
    "README.md; final version classification",
    "README.md; final version classification",
    "manual method synthesis",
    "manual method synthesis"
  ),
  script = "scripts/07_figures_tables/01_prepare_delivery_dirs_and_method_flow.R",
  dataset = "paired and non-paired",
  trait = "all",
  method = "final workflow",
  status = "planned_output",
  note = "Generated when this script is run.",
  stringsAsFactors = FALSE
))

cat("Prepared method-flow files under: ", out_dir, "\n", sep = "")
