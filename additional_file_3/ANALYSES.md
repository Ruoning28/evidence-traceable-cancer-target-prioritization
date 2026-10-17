# Analysis entry points

The analysis scripts read Additional file 2 without modifying its frozen inputs. All outputs are written under the directory supplied with `--output`.

| Manuscript analysis | Entry point |
|---|---|
| Feature construction and axis permissions; S7 | `build_features.py`; `run_experiments.py --stage prepare` |
| Primary dual-axis scoring | `run_experiments.py --stage primary` |
| Benchmark recovery and target bootstrap; S1 | `run_experiments.py --stage enrichment` |
| Uniform and coverage-conditioned label permutations; S2 | `run_experiments.py --stage permutation`; `summarize_permutations.py` |
| Development-positive holdout; S3 and S12 | `run_experiments.py --stage oof` |
| Fixed-benchmark feature-family ablation; S4 | `run_experiments.py --stage ablation` |
| Reference classification, thresholds and case cards; S5–S6 and S10 | `run_experiments.py --stage primary`; `--stage thresholds` |
| Source roles and derivation paths; S8 and S14 | `revision_checks/identifier_source_audit/run_audits.py` |
| Matched-information and native comparators; S9 | `run_experiments.py --stage comparators`; `run_review.py` |
| Source-specific label-group transfer; S13 | `run_experiments.py --stage source-out` |
| Class VI and coverage; S15 | `run_experiments.py --stage support`; identifier/review-denominator audits |
| Measurement-conditional review and policy sensitivity; S16–S17 | `run_review.py`; `gate_engine.py`; `measurement_sources.py` |
| Benchmark origin subsets; S18 | `revision_checks/quantify_overlap.py` |
| Background and regularization sensitivity; S19 | `run_experiments.py --stage background`; `--stage regularization` |
| Background resampling and positive-label thinning | `run_experiments.py --stage perturbation` |
| Figures 2–6 | `figure_code/generate_figure2.py`; `figure_code/generate_figures.py` |

`legacy-comparators` regenerates the original comparator panel retained in the machine-readable tables, including the merged-axis and fixed-formula routes. The current matched-information comparison uses `comparators`, whose percentile mean excludes missingness indicators and reverses the risk predictor. The reference policy uses the complete frozen candidate universe for quantiles; the measurement-review population is defined in `metadata/review_population_exclusions.tsv` in Additional file 2.

Primary reported scores use `--profile strict` (the default). `--profile reference` retains the archived reference solver tolerance; run its primary model first in a separate output directory before its diagnostics. The background and regularization summaries in S19 are the archived reference-configuration diagnostics. Newly generated diagnostics record the selected profile in their output. A `--stage all` run includes 500 refitted permutations per axis and null and can take substantial time.
