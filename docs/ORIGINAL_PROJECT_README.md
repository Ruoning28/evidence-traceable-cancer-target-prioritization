# Evidence-traceable PDAC/LUDA reproducibility bundle

This repository reproduces the analyses, source tables, validation records, and figures reported in the accompanying manuscript. A shared framework is applied to pancreatic ductal adenocarcinoma (PDAC) and lung adenocarcinoma (LUAD). The directory is named `LUDA` only to preserve the original package layout; all scientific outputs use the cancer code `LUAD`.

## Analysis design

The framework produces separate translational-precedence and functional-discovery rankings. For each axis, evidence records are assigned before model fitting to one of four roles:

- `F`: feature construction.
- `L`: development positive labels.
- `V`: independent validation benchmark.
- `A`: interpretation-only annotation.

Benchmark genes are embargoed from same-axis development labels and from the sampled unlabeled background. Development positives receive out-of-fold (OOF) scores. Held-out benchmark targets are used for top-K recovery, target-level bootstrap confidence intervals, comparator evaluation, and 500 label permutations. Full-source-out analyses remove both source-linked labels and source-linked feature columns.

## Project structure

- `main.py`: single command-line entry point.
- `framework/pipeline.py`: shared feature construction, L/V role assignment, leakage checks, sampled-background positive-unlabeled-style modelling, OOF and independent validation, full-source-out analysis, comparator models, sensitivity analysis, review classes, and evidence cards.
- `framework/experiment.py`: coordinates module analyses and figure/table generation.
- `framework/module_figures.py`: creates the three figures stored with each cancer module.
- `framework/paper_figures.py`: creates manuscript-level figures and their source CSV files.
- `framework/paper_tables.py`: creates manuscript LaTeX table fragments.
- `framework/verify.py`: checks reproduced values and output integrity.
- `tests/test_method_contract.py`: automated L/V separation, feature-exclusion, and resampling-configuration tests.
- `PDAC/`: PDAC configuration, processed inputs, data-download code, benchmark literature, extracted benchmark records, and results.
- `LUDA/`: LUAD configuration, processed inputs, data-download/build code, benchmark literature, extracted benchmark records, and results.
- `result/`: cross-cancer manuscript figures, tables, and verification logs.

Each cancer module contains:

- `config/framework_config.json`: all predeclared model, split, validation, and review settings.
- `data/processed/`: analysis-ready evidence and positive-source tables included for direct reproduction.
- `datasources/download_data.py`: public-data download entry point.
- `PDAC/datasources/audit_legacy_formula_scores.py`: recomputes and audits the legacy PDAC fixed-formula comparator columns shipped in the processed PDAC evidence matrix.
- `datasources/benchmark/papers/`: archived benchmark article records.
- `datasources/benchmark/extracted/`: target records extracted from benchmark sources.
- `result/tables/`: raw and summarized analysis output in TSV format.
- `result/figures/`: module-level PDF and PNG figures.

TSV means tab-separated values. It is a standard plain-text table format, not a Codex-specific format. It can be opened with spreadsheet software or read with `pandas.read_csv(path, sep="\t")`.

## Installation

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

## Reproduce the analyses

Run both cancer modules, regenerate all figures and tables, and verify the outputs:

```bash
python main.py
```

Run one module:

```bash
python main.py --module PDAC
python main.py --module LUDA
```

Reuse existing analysis tables while regenerating figures and verification reports:

```bash
python main.py --skip-analysis
```

Run verification only:

```bash
python main.py --verify-only
```

Run the method-contract tests:

```bash
python -m unittest discover -s tests -v
```

The formal configuration uses 30 sampled-background repeats, 500 label permutations, and 30 model repeats within every permutation. Runtime therefore depends on CPU performance.

## Download public raw data

The processed inputs required for exact reproduction are included. Raw public datasets are much larger and are downloaded only when requested:

```bash
python main.py --download-data --module PDAC
python main.py --download-data --module LUDA
```

Existing files are retained unless `--force-download` is added. To rebuild LUAD processed inputs after downloading the required LUAD and shared annotations:

```bash
python main.py --build-luda-inputs --module LUDA --skip-analysis --skip-figures
```

The original PDAC preprocessing workflow used a substantially larger raw-data collection. This bundle therefore includes the exact processed PDAC inputs used by the analysis, together with the migrated download scripts and source manifests needed to trace those inputs.

The PDAC processed matrix also contains two legacy fixed-formula comparator columns, `clinical_actionability_score` and `discovery_potential_score`, imported from the earlier PDAC target-prioritization workflow. Their provenance can be audited without the older project:

```bash
python -m PDAC.datasources.audit_legacy_formula_scores
```

This command recomputes the PDAC disease-relevance, safety-factor, clinical-actionability, and discovery-potential formula columns from the component scores in `PDAC/data/processed/pdac_target_evidence_matrix.tsv`, writes `PDAC/result/tables/pdac_legacy_formula_score_audit.tsv`, and exits with an error if any audited column differs beyond numerical tolerance.

## Result files

The most important raw review and validation tables in each module's `result/tables/` directory are:

- `framework_record_role_manifest.tsv`: record-level F/L/V/A assignments and split reasons.
- `framework_development_positive_sets.tsv`: development labels used for model fitting.
- `framework_validation_benchmarks.tsv`: embargoed held-out benchmark records.
- `framework_leakage_audit.tsv`: L/V overlap and role-conflict checks.
- `framework_feature_manifest.tsv`: feature source, transformation, missing-value rule, and axis inclusion.
- `framework_sampling_and_model_specification.tsv`: executable modelling settings.
- `framework_review_rule_specification.tsv`: ordered review-class rules.
- `framework_priority_scores.tsv`: final scores, ranks, uncertainty, and top-K frequencies for all genes.
- `framework_oof_positive_scores.tsv`: OOF scores and ranks for development positives.
- `framework_oof_validation_summary.tsv`: OOF recovery summaries.
- `framework_independent_benchmark_validation.tsv`: held-out benchmark recovery, enrichment, rank percentiles, and bootstrap intervals.
- `framework_leave_source_out_validation.tsv`: full-source-out transfer results.
- `framework_label_permutation_test.tsv`: 500-permutation null comparison.
- `framework_comparator_summary.tsv`: primary and comparator methods evaluated on the same held-out benchmark.
- `pdac_legacy_formula_score_audit.tsv`: PDAC-only audit showing that shipped legacy fixed-formula comparator columns can be recomputed from documented component scores.
- `framework_review_classes.tsv`: review class assigned to every candidate gene.
- `framework_evidence_cards.tsv`: target-level evidence and limiting-factor records.

Cross-cancer manuscript outputs are written to:

- `result/paper_figures/`: PDF/PNG figures.
- `result/paper_figures/source/`: figure source CSV files.
- `result/paper_tables/`: LaTeX table fragments.
- `result/logs/result_verification.tsv`: numeric and structural verification report.

