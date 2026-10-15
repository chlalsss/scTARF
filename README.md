# scTARF

Code and data accompanying *scTARF: Transcriptome-Assisted Trait Relevance
Mapping in Single-Cell Chromatin Accessibility Data*.

This repository provides two complementary paths:

- **Experimental workflows:** preprocessing, scDRS scoring and seed generation,
  reference-assisted non-paired mapping, sparsification, ablation, efficiency,
  reverse analysis, and figure/table scripts under `scripts/`. Start with the
  [workflow guide](docs/reproduction.md) and its ordered execution manifest.
  Preview a stage with `Rscript --vanilla scripts/run_workflow.R nonpaired --dry-run`.
- **Frozen-input result replay:** `paper_reproduction/` recomputes stored cell
  scores and evaluation results from frozen graphs, seeds, baseline/scDRS
  scores and annotations, preserving the historical evaluation rules.

The workflow sources and the replay package are retained together. See
[validation scope and provenance](docs/workflow-restoration.md) before running
the experimental workflows; restoring sources does not establish that every
workflow reproduces the final manuscript from raw data.

## Install and run the frozen-input replay

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

- `scripts/`: experimental workflows, including non-paired mapping and actual
  seed generation; [execution manifest](scripts/final_reproduction_manifest.tsv).
- `configs/`, `environment.yml`, and `data/`: workflow path templates,
  environment specification, and public input-data manifest.
- `results/expected_metrics.tsv`: manuscript reference values, not execution
  evidence or substitutes for computed outputs.

Generated outputs are excluded from version control. Earlier development
iterations remain accessible through Git history.

## Data and validation scope

The 87 selected data/reference files total approximately 295 MB. Details and
Windows locale instructions are in the [replay guide](paper_reproduction/README.md).
See [case definitions](paper_reproduction/replay_manifest.tsv),
[data checksums](paper_reproduction/input_checksums.tsv) and
[required source checks](paper_reproduction/required_sources.tsv).

The replay command starts from frozen intermediate inputs. The additional
experimental workflows provide upstream and supplementary analysis code with
separate dependencies and input requirements. They are not executed by this
command. Full raw-data reproduction and final manuscript figure layout
equivalence remain unverified. The exact author-selected GWAS table and some
historical ablation generators remain outside this release; see the explicit
[remaining requirements](docs/workflow-restoration.md#validation-boundary).

## Citation and license

Yushan Ji, Bingbo Wang, and Xianan Dong. *scTARF: Transcriptome-Assisted Trait
Relevance Mapping in Single-Cell Chromatin Accessibility Data.* Manuscript,
2026. Citation details and the repository license remain pending author review.
Third-party packages and source datasets retain their respective terms.
