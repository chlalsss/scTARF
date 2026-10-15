# Common-cell evaluation diagnostic for paired BMMC

This supplement documents the existing 2026-10-09 diagnostic of the evaluation-cell differences between SCAVENGE and scTARF. It evaluates **saved scores** on shared barcodes; it does not rerun seed selection, graph construction, or propagation, and does not replace manuscript figures or tables.

## Scope and results

For each trait, the shared set is the intersection of annotated cells with finite scores in three methods: SCAVENGE baseline, RNA-guided scTARF with an ATAC-only graph (`w_ATAC=1.0`), and RNA-guided scTARF with representation fusion (`w_ATAC=0.7`). The proxy definitions match `ppt_dataset_figures/03_plot_paper_roc_curves.R`; LYM uses the broad `LYM_wide` population. RNA-only scDRS, intersection-seed ablation, non-paired datasets, and downsampling are outside this diagnostic.

| Trait | Baseline evaluated cells | Each scTARF setting | Shared cells | Original fusion minus baseline AUROC | Shared-set difference |
|---|---:|---:|---:|---:|---:|
| MONO | 69,212 | 68,303 | 68,267 | 0.443811 | 0.443442 |
| PLT | 69,204 | 68,303 | 68,259 | 0.146653 | 0.146267 |
| LYM | 69,200 | 68,303 | 68,256 | 0.153583 | 0.153302 |

Using the historical pROC automatic directions, the maximum absolute difference between original-set and shared-set AUROCs across the nine method/trait results is **0.0013469489066794**. The method ordering is unchanged. This indicates a small numerical effect of restricting evaluation to shared cells in these saved-score comparisons. It does not establish fairness under identical preprocessing or show that all possible benchmark differences are negligible. The two scTARF settings already use the same cells.

## AUROC direction is a separate issue

The CSV reports both the historical pROC automatic direction and a fixed higher-score-is-positive direction. For the LYM SCAVENGE baseline, pROC selects `>` (lower scores rank proxy-positive cells), whereas both scTARF settings select `<` (higher scores rank proxy-positive cells). The original LYM baseline AUROC is 0.703240 under the historical convention and 0.296760 under the fixed higher-score convention; shared-set values are 0.703795 and 0.296205, respectively. Restricting cell sets does **not** resolve that direction difference. Conclusions about the cell-set effect are conditional on keeping the scoring convention unchanged.

## Files and reproduction

- `historical_diagnostic_20261009.json`: archived local diagnostic, preserving its original fields and scope.
- `common_cell_auroc.csv`: nine rows reproduced from existing saved inputs, including cell counts and both AUROC directions.
- `recompute_common_cells.R`: read-only input evaluation script. Requires R and pROC.
- `input_checksums.csv`: MD5 identifiers of the saved annotation and score inputs, with paths relative to the BMMC root.
- `session_info.txt`: environment of the 2026-10-11 diagnostic reproduction.

Required local inputs under a BMMC data root:

```text
01_data_qc/BMMC_cell_meta.rds
03_scavenge_baseline/TRS_matrix_baseline.rds
04_dual_modality/TRS_v2_{MONO,PLT,LYM}_w{0.7,1.0}.tsv.gz
```

Run from any working directory, providing the data root and a separate output directory:

```sh
Rscript validation/common_cells/recompute_common_cells.R /path/to/BMMC /path/to/diagnostic_output
```

Cell-level third-party input files are not included. This command recomputes the diagnostic from saved scores; it is not an end-to-end reproduction of the models. The archived diagnostic and the recomputed results agree to numerical precision.
