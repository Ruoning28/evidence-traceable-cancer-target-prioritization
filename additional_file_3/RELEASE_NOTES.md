# Code release, 17 October 2026

This release separates analysis and pixel-verified plotting environments. It adds exact Conda builds, a hashed isolated Figure 2 wheel, runtime/font checks, commands for each figure and archive verification. The native Figure 3–6 generator gains only an optional figure-selection argument; drawing functions and analysis algorithms are unchanged.

A new environment was recreated from the lock and all Figures 2–6 were generated with zero pixel differences against the stated manuscript at 300 DPI under matching PDF transparency-group conditions. Each individual figure command and the runtime mismatch guard were also checked. See `verification/release_validation.json` and `figure_code/README.md` for scope and requirements.

The Windows Arial font is required for the verified pixel profile but is not redistributed. Analysis inputs are supplied separately in Additional file 2. Figure 1 remains a static manuscript illustration. This release does not claim a new full-model refit, raw-data reconstruction or remote GitHub verification.

## BMC submission packaging

The manuscript now lists data split across Additional files 2 and 4 to keep each standard ZIP below 20 MB. Extract them together. Figure 2 wrapper instructions were clarified for that split. The current manuscript fingerprint and submission-validation record were added; scientific Python, numerical inputs and drawing functions were not changed.
