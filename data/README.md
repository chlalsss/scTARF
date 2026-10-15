# Data

Raw third-party datasets are not included in this anonymized review repository.

The analyses use public datasets and required preprocessed inputs listed in
`data_manifest.tsv`. The manifest gives expected filenames, local relative
paths, and preparation scripts for the minimal review workflow. Download the
corresponding source files and place them under `local_data/`, restore the
listed preprocessed objects, or set the matching environment variables in
`configs/paths.example.R`.

Some non-paired inputs are derived SCAVENGE/Granja objects rather than raw
third-party files. These are marked as "Required preprocessed input; not
redistributed." in the manifest and must be restored locally before running the
non-paired `run` entries.

The manuscript reports analyses on the BMMC paired multiome dataset, Granja
blood ATAC data, Granja RNA reference data, and CITE-seq RNA reference data.
