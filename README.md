# Evidence-traceable cancer target prioritization framework

This repository contains the analysis code for the manuscript:

**An evidence-traceable framework for cancer therapeutic target assessment and prioritization from public multi-source data**

## Overview

This repository provides code for a source-role-based, evidence-traceable framework for cancer therapeutic target assessment and prioritization using public multi-source data.

The workflow includes:

- cancer-specific evidence harmonization
- target-level feature construction
- sampled-background positive-unlabeled-style prioritization
- translational-precedence and functional-discovery scoring
- review-class assignment
- held-out benchmark recovery analysis
- label-permutation validation
- out-of-fold and source-out validation
- comparator analyses
- evidence-card generation
- figure and table generation

## Repository structure

```text
src/                    Core functions and reusable modules
scripts/                Analysis scripts
config/                 Non-sensitive configuration files
docs/                   Documentation and workflow notes
figures/                Figure-generation scripts, source data, and small preview figures
examples/               Small manifest-style example inputs
supplementary_manifest/ Data and supplementary-file manifests
```

The original analysis package used a module directory named `LUDA` for the LUAD cancer analysis. The name is preserved in `src/LUDA/` to keep imports reproducible.

## Requirements

The analysis was performed using Python.

Python dependencies are listed in `requirements.txt` and can be installed with:

```bash
pip install -r requirements.txt
```

For the source-layout repository, run commands from the repository root and expose `src/` on the Python path when needed:

```bash
set PYTHONPATH=src
```

On PowerShell:

```powershell
$env:PYTHONPATH = "src"
```

## Usage

The scripts should be run in the order described below. Some large processed matrices and public-source downloads are not included in this GitHub repository; see `supplementary_manifest/` for data availability and output-file descriptions.

```bash
# Run the complete workflow if the required input data are available
python scripts/main.py

# Cancer-specific analyses
python scripts/run_pdac_analysis.py
python scripts/run_luad_analysis.py

# Validation and verification checks
python scripts/run_validation.py

# Generate paper figures and tables from available result files
python scripts/generate_figures.py
python scripts/generate_tables.py
```

Additional options from the original workflow include:

```bash
python scripts/main.py --module PDAC
python scripts/main.py --module LUDA
python scripts/main.py --module LUAD
python scripts/main.py --skip-analysis
python scripts/main.py --verify-only
python scripts/main.py --download-data --module PDAC
python scripts/main.py --download-data --module LUDA
python scripts/main.py --build-luda-inputs --module LUDA --skip-analysis --skip-figures
```

## Data availability

This repository contains analysis code and data manifests.

Large processed evidence matrices, label and benchmark manifests, source-trace files, validation summaries, and machine-readable supplementary tables are provided with the manuscript and/or in the associated public data repository.

Raw public data sources are cited in the manuscript.

## Reproducibility notes

The analyses use public multi-source data and cancer-specific processed evidence matrices. Some large input and output files are not stored directly in this GitHub repository. See `supplementary_manifest/` for file descriptions and data-source manifests.

Users should update local input and output paths in the example configuration files before running the analysis. No local absolute paths are required by default.

The `figures/source/` directory contains small machine-readable source tables for reproducing selected paper figures. The `figures/paper_preview/` directory contains small preview copies of final paper figures for orientation only.

## License

This code is released under the MIT License.

## Citation

If you use this code, please cite the associated manuscript:

Li R, Wang D, Lan W. An evidence-traceable framework for cancer therapeutic target assessment and prioritization from public multi-source data.

