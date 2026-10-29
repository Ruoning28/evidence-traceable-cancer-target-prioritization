# Revised manuscript code compared with the original repository

Baseline: `ab73688532fb96256041d7cd44360d1f8b2f1429` (original `main`). Release: the audited Additional file 3 package dated 17 October 2026.

## Current entry point

Use [`additional_file_3/README.md`](../additional_file_3/README.md) for the revised manuscript. Its 80 files are byte-for-byte identical to [`releases/Additional_file_3.zip`](../releases/Additional_file_3.zip). The archive SHA-256 is `f12e8ed927cecab2003ca7c80731a576a668e49c30f1c1a7f42c66a27b35e341`. `additional_file_3/sha256.json` covers 79 files and excludes itself.

The original `src/`, `scripts/` and `config/` remain as the baseline. Current supplementary manifests and previews now match the final submission; original versions remain in Git history. Their existence does not identify them as the current revised analysis or current figure previews; use the released entry points and packaged frozen inputs rather than mixing layouts.

## Substantive differences

| Area | Original repository | Revised package |
|---|---|---|
| Input/output handling | Module directories and large inputs outside the repository | Explicit `--data` and `--output`; frozen Additional files 2 and 4, extracted together into one additional_file_2 layout |
| Solver configuration | No explicit `model.tol`/initial weight-normalization fields | `tol=1e-6` and `normalize_class_weights_initially=true`; same-sample convergence retries through `fast_model` |
| Sensitivity repeats | `sensitivity.n_repeats=10` | `sensitivity.n_repeats=30` in both cancer configurations |
| Background strata | Dependency/expression fields use the original PDAC-oriented rules | Cancer-context dependency/expression columns; ties kept together; unavailable/constant expression handled explicitly |
| Feature metadata | Generic external-evidence labels, including an inaccurate LUAD description | Cancer-specific derivation descriptions; nonfinite model inputs handled explicitly; LUAD support identified as internal TCGA/DepMap evidence |
| Review/measurement and diagnostics | Original workflow | Separate operational review, observation/interval helpers, current comparator and revision-analysis entry points |
| Figure reproduction | Original plotting module and previews | Current Figures 2–6 source tables, isolated locked Windows rendering profiles, font preflight and individual figure commands |
| Distribution | Original source layout and unpinned top-level dependencies | Self-contained code archive with analysis requirements, exact plotting build locks, file manifest, packaging verifier and validation-scope records |

The two configuration objects retain the same seeds, main background ratio, `C`, `l1_ratio`, label/benchmark definitions, feature exclusions and review settings; their six changed flattened fields are listed in [`additional_file_3_configuration_changes.json`](additional_file_3_configuration_changes.json). Differences in JSON formatting are not additional parameter changes.

In the relocated core pipeline, 28 function ASTs are identical and eight differ; no pipeline functions were added or removed. The eight include documentation/output labels as well as real background-stratification, feature-input and model-compatibility changes. See [`additional_file_3_pipeline_functions.json`](additional_file_3_pipeline_functions.json). This is not a claim that all differences from the original repository are environment-only.

The final packaging/environment repair itself left the earlier revised analysis scripts and native plotting functions unchanged; its figure-selector argument does not alter drawing formulas. That narrower release-to-release statement is separate from this comparison against the older GitHub baseline.

## Mapping and verification

[`additional_file_3_vs_original.csv`](additional_file_3_vs_original.csv) maps each released file to a direct same-named/relocated baseline file where one exists. Nine have direct mapped counterparts: two are byte-identical and seven differ; 71 have no direct file counterpart in the original layout. Baseline hashes now use exact committed Git blobs (LF), not a CRLF-converted Windows checkout. Earlier all-nine-different counts included newline-only differences. These counts are not counts of changed algorithms and do not imply that original files were removed.

Before this submission, the Windows plotting environment was recreated from its 72-package exact Conda closure and the SHA-256-locked Figure 2 wheel. All five Figures 2–6 matched the referenced manuscript at 300 DPI under the same PDF transparency-group conditions, with zero differing pixels. Each single-figure command and the wrong-runtime preflight were checked. Reference hashes and the recorded scope are in [`additional_file_3/verification/release_validation.json`](../additional_file_3/verification/release_validation.json).

The Git submission separately checks that the archive and the extracted release tree match all 80 source files and that Git stores their original bytes. `.gitattributes` preserves the audited release tree and exact copied figure/table/data-manifest bytes against newline conversion. Runtime environments, cached dependencies, fonts, manuscript files and Additional file 2 inputs are not added to Git. Figure 1 remains a static manuscript illustration.

No new complete model ensemble, OOF, refitted permutation suite or raw-data reconstruction was run for the Git submission. Existing model-consistency files are retained supplied records, not new training evidence. This change does not modify the original branch's scientific files in place.

## Final submission synchronization (29 October 2026)

The final code tree has 80 files and 79 internal checksum entries. The 31 Python files remain byte-identical to the earlier PR release; five documents/manifests were updated and `verification/submission_validation.json` was added. Data come from Additional files 2 and 4 extracted together. Current previews, figure sources and main tables now match the final manuscript; attachment/source manifests record the actual supplied file identities. Historical validation records are retained without rewriting their original manuscript fingerprints.
