# Replay the manuscript's saved analysis

This entry point reproduces the stored manuscript computation using the actual
saved graphs, seed vectors, baseline scores, scDRS scores and cell annotations.
It preserves the original inputs, method combinations and evaluation behavior.
It does not perform coordinate conversion or replace the historical analysis.

## Run

Required R packages: Matrix, BiocNeighbors, data.table and pROC. The checked
environment used R 4.5.2, Matrix 1.7-4, BiocNeighbors 2.4.0, data.table 1.17.8
and pROC 1.19.0.1. These identify the tested environment, not minimum versions.
Use the existing repository environment instructions before running.

From the repository root, choose an output directory that does not exist:

```sh
Rscript --vanilla paper_reproduction/run_paper.R paper_reproduction/output/run_01
```

The code reads three unchanged source files from the repository:

- `examples/run_core_demo.R` supplies the checked original-function loader.
- `module3_multiome/04_dual_modality/03_build_weighted_knn.R` supplies definitions
  required by that loader; graphs in this replay are already saved inputs.
- `module3_multiome/04_dual_modality/05_rwr_trs.R` supplies original RWR and TRS
  constants/functions. Its data-loading and writing loops are not executed.

Only output location and path resolution are adapted. The runner checks frozen
input bytes and original source content before calculation. Source CRLF/LF
line-ending differences are normalized for checking. The `.gitattributes` file
preserves frozen input bytes on checkout. Existing output directories are
rejected so a prior result cannot silently stand in for a new calculation.

For the author machine's inherited `LC_ALL=C.UTF-8` issue on Windows, temporarily
use a valid child locale, restoring the shell value afterwards:

```powershell
$previousLocale = $env:LC_ALL
try {
    $env:LC_ALL = 'Chinese (Simplified)_China.utf8'
    Rscript --vanilla paper_reproduction/run_paper.R paper_reproduction/output/run_01
} finally {
    $env:LC_ALL = $previousLocale
}
```

## Inputs and outputs

The selected data/reference set has 87 files, approximately 295 MB.
`replay_manifest.tsv` identifies 30 graph/seed/reference cases:

| Setting | Cases | Cells per case |
| --- | ---: | ---: |
| Paired MONO/PLT/LYM, six weights 0.5 through 1.0 | 18 | 68,303 |
| Non-paired MONO/RDW/LYM, weights 0.7 and 1.0 | 6 | 35,000 |
| Reverse PLT/MCH/RBC, weights 0.7 and 1.0 | 6 | 35,582 |

The runner recomputes RWR/TRS for every case and writes cell-level scores into
`paired/`, `nonpaired/` and `reverse/` under the new output directory.
`replay_comparison.csv` checks each result against stored cell scores with
absolute tolerance 1e-12. Historical baseline and scDRS cell scores are fixed
inputs; their upstream enrichment/scDRS computation is not rerun here.

`auroc_comparison.csv` recalculates 61 AUROCs: 40 paired baseline/ablation/sweep
entries, three RNA-only entries, nine non-paired entries and nine reverse
entries. It records cell counts, positive counts, actual ROC direction, source
value and difference. The original paired automatic ROC direction and the
other workflows' fixed direction are retained. The earlier `v1_paired`
evaluator's first-two-column fallback is also retained when its expected barcode
column is absent. These are historical numerical-replay conventions.

`session_info.txt` records the execution environment. All 30 score comparisons
and all 61 AUROC comparisons must pass; otherwise the runner exits with an error.
No saved AUROC is substituted for a newly calculated value to make a check pass.

## Scope of a successful run

The validated start is frozen real graphs/seeds and historical baseline/scDRS
scores. Outputs are new cell scores and recomputed evaluation numbers.
This entry does not rebuild raw-count preprocessing, GWAS enrichment, RNA
mapping, graph construction, dropout experiments, efficiency timings or final
figure layouts. Those materials have their own saved sources and are not
claimed as recomputed by this command.

Comparisons target the identified calculation tables. Previously recorded
manuscript display/source differences remain unchanged: for example ablation
MONO/PLT display 0.500/0.537 versus table 0.502139/0.538324, and reverse MCH
baseline displays 0.482 versus table 0.481450336614053. Replaying the calculation
does not silently force a printed constant into its output. See the repository's
reproduction-status documentation for the complete recorded display differences.

Data/reference checksums are in `input_checksums.tsv`; required source checks
are in `required_sources.tsv`. Raw sequencing/count files and unused research
iterations are not included in this selected set. License and final citation
metadata remain pending author review; this local bundle has not been published.
