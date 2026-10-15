# Experimental workflow instructions

This guide runs the experimental generators. For the separately validated result
replay, use [paper_reproduction](../paper_reproduction/README.md). Read the
[validation boundary](workflow-restoration.md) before interpreting a run as
reproduction of the final manuscript.

## Environment and configuration

Run all commands from the repository root. [environment.yml](../environment.yml)
lists R/Bioconductor and Python dependencies, including Seurat, Signac, the hg19
annotation, scDRS and scikit-learn. It is a dependency specification, not a tested
lockfile. A fresh installation has not been validated.

```sh
conda env create -f environment.yml
conda activate sctarf-workflows
```

Install SCAVENGE in that R library, following its
[upstream installation instructions](https://github.com/sankaranlab/SCAVENGE#installation).
The dependency helper can install it when the repository is supplied explicitly:

```r
Sys.setenv(SCTARF_SCAVENGE_GITHUB = "sankaranlab/SCAVENGE")
source("scripts/01_paired_analysis/00_install_deps.R")
```

Restart R after installation. Copy [configs/paths.example.R](../configs/paths.example.R)
to `paths.local.R` in the repository root and edit paths. This local file is
ignored by Git. The runner loads it for each R step; Python inherits environment
variables. Alternatively, set environment variables in the launching shell.
`SCTARF_PYTHON` accepts a PATH command or executable path. `SCTARF_SCDRS` identifies
the Python scDRS entry script, which the historical wrappers invoke through Python
(on Windows, use the script rather than a binary `.exe` launcher).

The YAML files in `configs/` are historical parameter records; scripts do not
parse them. Use documented R overrides in `paths.local.R`. Use fresh data/output
roots when changing inputs, software or parameters: several original scripts
reuse graph/score caches without cache-key validation. Baseline and LSI generation
refuse to replace existing outputs.

## Required inputs

The [input manifest](../data/data_manifest.tsv) records paths and generators.
Cell counts in it are historical sizes, not checks on a new download.

- BMMC paired multiome: GSE194122 processed H5AD. Paired scripts obtain/extract
  metadata, peaks, matrices, background peaks and LSI.
- Granja ATAC/RNA: obtain the two `Healthy-Hematopoiesis-191120.rds` objects from
  the original [GreenleafLab resource page](https://github.com/GreenleafLab/MPAL-Single-Cell-2019#healthy-hematopoiesis)
  under `SCTARF_BLOOD_ROOT`. Use the peak-count ATAC object; the included generator
  computes this study's gene activity. The public Cicero gene-activity object
  does not implement the same recipe.
- CITE RNA: step 18 downloads the GSE194122 CITE H5AD; step 20 prepares its RNA.
- GWAS: supply `CS.trans_all.PP01.strong_atac_overlap.tsv`, containing `var`,
  `trait`, `PP` (`var` format such as `1:2164935_G_T`). Set `M3_GWAS_TSV` for paired
  analysis; the Granja baseline reads a copy under `SCTARF_BLOOD_ROOT`. The exact
  author-selected table has no verified download or generation/filtering pipeline
  in this release. Changing its variants or coordinate conventions changes the
  analysis; scientific reconciliation is outside this source-restoration release.

Gene activity, metadata, baseline TRS and LSI now have included generators.
The LSI adapter preserves the original TF-IDF, log transform, 50-component SVD,
random state 9527 and row normalization, with barcodes in count-matrix order.
Fresh LSI equivalence to the stored historical object remains unverified.

## Execution

Inspect a stage without loading data or creating outputs:

```sh
Rscript --vanilla scripts/run_workflow.R nonpaired --dry-run
```

The runner follows the [manifest](../scripts/final_reproduction_manifest.tsv),
launches each step in a fresh process, records logs under `outputs/run_logs/`,
records R session information, and stops on failure. It does not install packages
or automatically run prerequisite stages.

```sh
Rscript --vanilla scripts/run_workflow.R paired
Rscript --vanilla scripts/run_workflow.R nonpaired
Rscript --vanilla scripts/run_workflow.R sparsification
Rscript --vanilla scripts/run_workflow.R ablation
Rscript --vanilla scripts/run_workflow.R efficiency
Rscript --vanilla scripts/run_workflow.R reverse
Rscript --vanilla scripts/run_workflow.R figures_tables
```

`sparsification` requires paired preprocessing/scDRS results. `ablation`,
`efficiency`, and `figures_tables` consume paired/nonpaired results. `reverse`
consumes nonpaired inputs/results. The optional `nonpaired_dropout` stage runs
the original Granja entry-dropout baseline, distinct from paired read thinning.

After fixing a failed step, resume the remaining manifest entries individually
in fresh sessions. Example from R (choose the first unfinished row from logs):

```r
source("scripts/run_workflow.R")
p <- workflow_plan("nonpaired")
if (file.exists("paths.local.R")) source("paths.local.R")
dir.create("outputs/resume_logs", recursive = TRUE, showWarnings = FALSE)
for (i in 7:nrow(p)) { # example: resume at row 7
  workflow_run_step(p$script[i], sprintf("outputs/resume_logs/%02d.log", i))
}
```

## Nonpaired dependency chain

1. Generate Granja gene activity/metadata, run original SCAVENGE baseline,
   standardize with step 54, evaluate with step 55, and generate LSI.
2. Prepare and score Granja RNA and CITE RNA using scDRS.
3. Run `32_build_legacy_dense_mapping.R`, then
   `64_diagnose_scdrs_neighborhood_fdr.R` to produce historical RNA neighborhoods
   and MONO seed candidates. Step 29 reuses membership for additional traits.
4. Run final drivers 03/04/05 for MONO/LYM/RDW. Their helpers generate sparse
   mapping, mapped-expression completion, seeds, graphs, TRS and evaluations.

The MONO completion uses sparse mapping; its historical seed candidates come
from the earlier neighborhood chain. Both paths are retained. Substituting new
sparse-route neighborhoods would change the historical analysis.

## Dropout, sensitivity and efficiency

Paired dropout defaults: depths `seq(0, 0.9, 0.1)`, seeds `1, 2, 3`, ATAC weights
`0.7, 1.0`. Override `P_DROP_LIST`, `SEED_LIST`, `ATAC_W_LIST` in `paths.local.R`
for smaller runs and record the changes. The paired weight-sweep generator uses
`0.5, 0.6, 0.7, 0.8, 0.9` from `paired_config.yaml`, plus control `1.0`, through
the original graph/RWR/evaluation code.

The ablation stage builds master tables before summaries. Historical v1
intersection-seed variants and optional `rna_source_diagnostics_summary.csv`
are not regenerated by this stage. Summaries of available variants must not be
interpreted as a complete historical ablation comparison.

Set `EFFICIENCY_FORCE_RERUN=true` for new measurements; otherwise existing
summary files can be reused. Defaults are MONO, 1k/10k/50k/100k sampled cells,
ATAC weight 0.7. Above the available cell count sampling uses replacement.
Memory estimates use R heap/object sizes, not OS peak resident memory. Record
hardware, package versions, threads and cache settings with timings.
Final manuscript figure layout equivalence is unverified. Values in
`results/expected_metrics.tsv` are references, not computed evidence.

## Lightweight checks

```sh
Rscript --vanilla tests/test_workflow_configuration.R
Rscript --vanilla tests/test_workflow_paths.R
Rscript --vanilla tests/test_workflow_runner.R
Rscript --vanilla tests/test_result_weights.R
python -m unittest discover -s tests -p "test_*.py"
```

Path/Python checks require NumPy/SciPy/scikit-learn. They exercise configuration,
process execution, input ordering and failure propagation, not full experiments.
