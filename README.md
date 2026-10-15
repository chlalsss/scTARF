# scTARF

Code and data accompanying *scTARF: Transcriptome-Assisted Trait Relevance
Mapping in Single-Cell Chromatin Accessibility Data*.

This selected snapshot recomputes the manuscript's stored cell scores and
evaluation results from its frozen real graphs, seeds, baseline/scDRS scores
and cell annotations. It preserves the historical inputs and evaluation rules.

## Install and run

The local validation used R 4.5.2, Matrix 1.7-4, BiocNeighbors 2.4.0,
data.table 1.17.8 and pROC 1.19.0.1. In a separate R environment, install
the dependencies from their package sources:

```r
install.packages(c("Matrix", "data.table", "pROC", "BiocManager"),
                 repos = "https://cloud.r-project.org")
BiocManager::install("BiocNeighbors", update = FALSE, ask = FALSE)
```

The installation commands have not been validated in a fresh library; use
a Bioconductor release compatible with your R installation. Then run from
the repository root, choosing an output directory that does not exist:

```sh
git clone https://github.com/chlalsss/scTARF.git
cd scTARF
Rscript --vanilla paper_reproduction/run_paper.R output/run_01
```

The run checks the selected input/source hashes, recomputes 30 real-data
RWR/TRS cases, and recalculates 61 AUROCs against the saved calculation tables.
It exits with an error if a comparison fails. Generated results include
cell-level TSVs, score/AUROC comparison CSVs and R session information.

## Repository contents

- `paper_reproduction/`: the reproduction entry point, evaluation code, case
  definitions, input checksums and frozen data/reference files.
- `module3_multiome/04_dual_modality/`: original graph and RWR/TRS function
  definitions used by the reproduction entry point.
- `examples/run_core_demo.R`: the original-function loader and synthetic example.

Generated outputs are excluded from version control. Earlier development
iterations remain accessible through Git history.

## Data and validation scope

The 87 selected data/reference files total approximately 295 MB. Details and
Windows locale instructions are in the [replay guide](paper_reproduction/README.md).
See [case definitions](paper_reproduction/replay_manifest.tsv),
[data checksums](paper_reproduction/input_checksums.tsv) and
[required source checks](paper_reproduction/required_sources.tsv).

This entry starts from frozen intermediate inputs. Raw-count preprocessing,
GWAS enrichment, scDRS scoring, mapping, graph reconstruction, dropout runs,
efficiency measurement and final figure layouts are outside this command.

## Citation and license

Yushan Ji, Bingbo Wang, and Xianan Dong. *scTARF: Transcriptome-Assisted Trait
Relevance Mapping in Single-Cell Chromatin Accessibility Data.* Manuscript,
2026. Citation details and the repository license remain pending author review.
Third-party packages and source datasets retain their respective terms.
