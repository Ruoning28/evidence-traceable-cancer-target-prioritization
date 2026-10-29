# Evidence-traceable cancer target prioritization framework

Analysis code and reproducibility manifests for:

**An evidence-traceable framework for cancer therapeutic target prioritization using heterogeneous biomedical evidence**

## Final revised manuscript

The current analysis package is **[additional_file_3/](additional_file_3/README.md)**. Its 80 files are byte-identical to the final submitted Additional file 3 and the downloadable **[Additional_file_3.zip](releases/Additional_file_3.zip)**. The [archive checksum](releases/Additional_file_3.zip.sha256) and internal 79-entry [file manifest](additional_file_3/sha256.json) verify those exact bytes.

The matching reference manuscript PDF has SHA-256 `14c6e6d2bbc368c1ad345505759a8f336a42eeeddf475c7613b4b64c7328a7ff`. [Submission artifact identities](supplementary_manifest/submission_artifacts.tsv) record the LaTeX archive, figures, tables and supplementary files. [Comparison with the original repository](docs/ADDITIONAL_FILE_3_COMPARISON.md) describes the earlier analysis changes and the later packaging-only synchronization.

## Data and supplementary files

| File | Contents |
|---|---|
| Additional file 1 | Supplementary methods and results S1–S19 |
| Additional file 2 | Frozen model inputs, configurations, ranking/review results, observation and source metadata, ledgers and the combined data checksum manifest |
| Additional file 3 | Analysis and figure-reproduction code, environment specifications and running instructions |
| Additional file 4 | PDAC and LUAD source-evidence matrices |

Obtain Additional files 2 and 4 with the manuscript and extract both into the repository root, preserving their paths. Together they restore `additional_file_2/`, adjacent to `additional_file_3/`. Verify all 179 input-file entries in `additional_file_2/sha256.json`; it includes the two source matrices from Additional file 4.

See the [data and supplementary manifests](supplementary_manifest/README.md), [source-resource overview](supplementary_manifest/source_data_manifest.tsv), [resource acquisition identities](supplementary_manifest/source_metadata/resource_acquisition_register.tsv) and [source-role records](supplementary_manifest/source_roles/). PDAC GEO cross-cohort support and LUAD internal TCGA/DepMap support have different constructions; the latter is not independent external replication.

## Run the current analysis

Create the analysis environment described in [the package README](additional_file_3/README.md), using Python 3.10 and the pinned dependencies. From the repository root:

```text
cd additional_file_3
python -m pip install -r requirements.txt
python run_experiments.py --data ../additional_file_2 --output results/strict --stage prepare
python run_experiments.py --data ../additional_file_2 --output results/strict --stage primary
python run_review.py --data ../additional_file_2 --output results/review
```

[ANALYSES.md](additional_file_3/ANALYSES.md) maps primary, comparator, permutation, holdout, source-out, measurement and sensitivity stages to the manuscript. Full permutation/holdout refits can be expensive; output directories are explicit and frozen inputs are not modified.

## Figures and tables

Current [six figure assets and captions](figures/README.md), [machine-readable figure source tables](figures/source/) and [nine editable main tables](tables/README.md) match the final submission. Figure 1 is static artwork. For Figures 2–6 use the [locked figure environments and per-figure commands](additional_file_3/figure_code/README.md). The Figure 2 wrapper reads the fully restored data directory; the split Additional_file_2.zip alone is insufficient. Plotting environments are separate from the analysis environment.

## Original code and validation scope

The root `src/`, `scripts/`, `config/`, `requirements.txt` and the original project documentation preserve the historical source layout. Use the current `additional_file_3` instructions for this revision; use the baseline in [SOURCE_ORIGIN.json](additional_file_3/SOURCE_ORIGIN.json) to retrieve a complete historical checkout. Original figure source tables are retained in `figures/legacy_source/`.

Validation records distinguish previous model checks, figure reconstruction and submission packaging. This repository synchronization verifies file identities, data manifests and entry-point compatibility; it does not claim a new full model fit, refitted permutation suite or complete upstream raw-resource reconstruction. [Current submission validation](additional_file_3/verification/submission_validation.json) and earlier release records retain their original scope and fingerprints.

## License and citation

Code is released under the MIT License; source datasets retain their respective licences. Please cite the associated manuscript:

Li R, Wang D, Lan W. An evidence-traceable framework for cancer therapeutic target prioritization using heterogeneous biomedical evidence.
