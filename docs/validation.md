# Release validation — 2026-10-15

Checks below were run on the integrated release candidate in Windows, using
R 4.5.2 and Python 3.13.12. NumPy 2.4.4, SciPy
1.17.1, scikit-learn 1.8.0. R replay packages:
Matrix 1.7-4, data.table 1.17.8, pROC 1.19.0.1. The local pROC build reports a
R 4.5.3 build-version warning; the checks nevertheless exit successfully.

| Check | Result |
| --- | --- |
| Original formal replay/core/data files | 98/98 byte-identical |
| Original anonymous workflow files | 95/95 non-conflicting paths retained; adaptations recorded in provenance |
| R syntax | 95 files parsed successfully |
| Python syntax | 13 files parsed successfully |
| Execution manifest | 62 entries reference included scripts |
| Configuration | Full dropout defaults and caller overrides pass |
| Process paths | PATH resolution and paths with spaces pass; real Python LSI export preserves values/barcodes |
| Stage runner | Prerequisite ordering, fresh subprocess execution and failure propagation pass |
| Result table | All six paired graph weights survive into master output |
| Python tests | 2 pass: wrapper executable resolution; actual synthetic LSI CLI, finite normalized output and overwrite guard |
| Frozen real-data replay | 30/30 score cases and 61/61 AUROC checks pass |
| Maximum score absolute difference | 6.66134e-16 |
| Maximum AUROC absolute difference | 5.55112e-16 |

The LSI function's syntax tree matches the original local implementation.
No full scientific experiment was run by these small tests. The real-data
replay starts from frozen intermediate inputs and historical baseline scores.
Its success does not validate raw-input regeneration, fresh environment setup,
GWAS preprocessing, complete historical ablations, or manuscript figure layouts.
See [remaining requirements](workflow-restoration.md#validation-boundary).
