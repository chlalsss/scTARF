# Source provenance and validation scope

Experimental generators and frozen-input replay are distributed together.
Replay does not replace preprocessing, seed selection, reference mapping,
dropout or efficiency experiments.

## Retained sources

- Replay: `chlalsss/scTARF` commit `7ab5fa6ecdc67eff54ff9ad7a696c1a7b80a4c78`.
  All 98 files except the root README/ignore rules remain byte-identical,
  including frozen data, checksums, reference tables and original-function code.
- Workflow base: `artifact-review-4621/artifact-4621` commit
  `23b831fee7b7efc86465a53f816bc47b49f327be`. All 95 non-conflicting source files
  are retained. Selected files receive path/configuration/dependency repairs;
  they are no longer all byte-identical to the source snapshot.
- Additional original local generators supply Granja gene activity, baseline,
  legacy dense mapping and optional baseline dropout. The LSI adapter copies
  the original `compute_lsi` function unchanged.

The [file-level provenance table](source-provenance.tsv) records original paths,
source hashes and integrated hashes; new wrappers/tests are marked as such.
Local-source SHA-256 identifies inspected bytes, not a claimed Git revision.

## Integration repairs

The release restores upstream sources, adds an ordered stage runner, fills
Granja preprocessing and historical MONO seed dependencies, resolves executable
names through PATH, and quotes process arguments. Paired dropout exposes the
full depth grid and configurable seeds/weights. A weight-sweep entry reuses
existing functions. Master-table generation precedes ablation consumers and
preserves numeric graph weights. Reverse analysis reads the configured
nonpaired root rather than the former local directory layout.

Scientific input files, score formulas, historical seed rules and frozen replay
references are preserved. Additional dropout/sweep runs have not been substituted
for historical results.

## Validation boundary

See [release validation](validation.md) for checks and environment. Source syntax,
manifest order/files, replay-file hashes, documentation links, configuration,
process handling (including spaces), and small synthetic LSI inputs are checked.
Real-data frozen-input replay separately compares 30 RWR/TRS cases and 61 AUROCs.

Full raw-input experiments and fresh dependency installation are unverified.
The local R environment runs the replay, but loading the full SCAVENGE/gchromVAR
stack crashed in the Windows diagnostic. Syntax and small checks therefore do
not establish full experiment success.

Remaining requirements for independent manuscript regeneration:

1. Supply/version the exact author-selected GWAS table and upstream filtering
   provenance; public raw data do not determine that table.
2. Run the restored preprocessing, seed/mapping, dropout and efficiency stages
   in a clean environment with original inputs. Confirm fresh LSI/neighborhoods
   against stored historical objects.
3. Restore additional historical v1/intersection-ablation generation if that
   comparison is required. Summary scripts currently consume available variants;
   the optional RNA-source diagnostic summary is not built.
4. Verify manuscript figure layouts and timing/memory measurement conditions.
   Efficiency memory estimates are not OS peak RSS.

Use the [workflow guide](reproduction.md) for commands and cache behavior. Source
restoration makes the method inspectable and provides generators with specified
inputs; it does not certify all manuscript results from raw sequencing data.
