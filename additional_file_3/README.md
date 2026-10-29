# PDAC/LUAD evidence-traceable target prioritization

This project implements the dual-axis sampled-background model, benchmark assessment, reviewer-requested diagnostics, operational review and figure generation for the manuscript.

## Analysis setup

Use Python 3.10 and install the pinned analysis dependencies. These describe the analysis environment, not a pixel-identical plotting lock:

```text
python -m pip install -r requirements.txt
```

Excel inputs use openpyxl; unreadable TTD or literature workbooks raise an error identifying the affected file.

Extract both Additional_file_2.zip and Additional_file_4.zip into the same parent directory; together they restore `../additional_file_2`, alongside this repository. Run the following commands from the repository root. Analysis inputs are distributed separately across Additional files 2 and 4.

```text
python run_experiments.py --data ../additional_file_2 --output results/strict --stage prepare
python run_experiments.py --data ../additional_file_2 --output results/strict --stage primary
python run_experiments.py --data ../additional_file_2 --output results/comparators --stage comparators
python run_review.py --data ../additional_file_2 --output results/review
```

`primary` writes the complete score vectors, ranks, coefficient summaries, fit/convergence audit, benchmark recovery, reference classes and evidence cards. `run_primary.py` is a shorter aggregate-recovery check. Both use explicit elastic-net regularization. The supplied runtime reproduces the original strict score vectors to floating-point precision; the original analysis runtime is recorded in Additional file 2 metadata.

## Diagnostics

`run_experiments.py` provides `oof`, `source-out`, `ablation`, `perturbation`, `permutation`, `background`, `regularization`, `thresholds`, `support` and `legacy-comparators` stages. `ANALYSES.md` maps them to the manuscript. Each invocation accepts explicit data and output directories, cancer and solver profile. Primary strict scores supplied in Additional file 2 provide the reference for diagnostics run in a new output directory. The reference-tolerance profile requires its own primary run first.

```text
python run_experiments.py --data ../additional_file_2 --output results/strict --stage source-out
python run_experiments.py --data ../additional_file_2 --output results/strict --stage permutation --jobs 4
python summarize_permutations.py --results results/strict --output results/strict/permutation_summary.tsv
python revision_checks/identifier_source_audit/run_audits.py --data ../additional_file_2 --output results/identifier_source_audit
python revision_checks/quantify_overlap.py ../additional_file_2 results/overlap
```

Source-role-out removes the documented label groups and source-specific feature lists. Fixed-benchmark ablation uses its separately recorded feature lists and strata changes. `gate_engine.py` implements the ordered review rules and 31 policies; `measurement_sources.py` supplies the raw constraint observation and interval functions.

## Feature and figure code

`build_features.py` rebuilds model inputs and feature manifests from the packaged source matrices. `scientific_source/reference` retains the upstream PDAC and LUAD source-matrix constructors and their local import dependencies. `source_config/PDAC` provides the corresponding raw-source configuration, while `config/PDAC` and `config/LUAD` provide the frozen model configurations. Historical resource identities and acquisition dates are in Additional file 2.

`build_source_matrix.py` exposes the source-matrix constructors for an existing raw-source layout. For PDAC, `--source-root` points to the original PDAC project containing `config`, `data/raw` and `data/processed`; retained preprocessing helpers can write their intermediate tables there. For LUAD, it points to the LUAD module containing `data/raw`, and `--shared-pdac-root` identifies the PDAC project supplying shared annotations and DepMap files. `--output` selects the target-matrix destination. Model reproduction uses the packaged frozen inputs.

```text
python build_features.py --data ../additional_file_2 --output results/features
```

For Figures 2–6, follow [figure_code/README.md](figure_code/README.md). It provides the tested Windows plotting lock, separate Figure 2 wheel, font checks and commands for every figure. Do not install the plotting wheel into the analysis environment. `figure_code/reproduce_figures.py` validates the selected runtimes and renders from the actual extracted Additional file 2 inputs; it does not retrain models. Figure 1 is a static illustration distributed with the manuscript/LaTeX materials.

`scientific_source/core/src/framework/pipeline.py` contains the reusable feature, model, class and evidence-card functions. `scientific_source/core/scripts/fast_model.py` implements the converged ensemble, and `run_revision.py` provides frozen-input loading for the diagnostics. `SOURCE_ORIGIN.json` identifies the original repository baseline. `sha256.json` lists file hashes. Code is distributed under the MIT License; source datasets retain their respective licences.

## Source-construction checks

```text
python -m unittest discover -s tests -v
```

These checks cover the PDAC table-saving entry point, TTD Excel parsing and workbook failure reporting using temporary input files.

## Code archive and verification

`sha256.json` covers every distributed code/document/table file except itself. After reviewed changes, refresh it and create a new archive outside this directory:

```text
python tools/package_code.py --refresh-manifest --output ../Additional_file_3.zip
```

The command refuses to overwrite an existing archive, excludes runtime environments and generated results, and verifies the ZIP CRC, entry list and SHA-256 of every extracted file against this directory. It writes an archive checksum and verification report beside the ZIP. Analysis inputs remain in Additional file 2.

## Validation scope

`verification/release_validation.json` records the environment rebuild, frozen-table figure generation and reference comparisons performed for this code release. The retained primary-model consistency files under `reproduction_checks/` are earlier supplied run records. Rendering the figures does not repeat all primary fits, out-of-fold analyses or refitted permutations; no new full-model reproduction claim is made by this release.

The current submission packaging checks are in `verification/submission_validation.json`. The earlier `release_validation.json` retains the prepackaging manuscript hash as historical provenance. The current figure reference hash points to the submission PDF; drawings and all frozen numerical inputs are unchanged.
