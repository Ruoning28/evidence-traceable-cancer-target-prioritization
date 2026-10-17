# Recreating Figures 2–6

The plotting code and frozen source tables reproduce the analysis figures. `figure_source/` and `display_labels.json` are included here. Figure 1 is a static illustration distributed with the manuscript/LaTeX materials.

## Rebuild the plotting runtimes

The pixel-verified profile is **Windows x86-64**, separate from the Python 3.10 analysis environment in the root README. From the code archive root, using a Conda installation:

```text
conda create --prefix .runtime/figures --file environments/figures-win-64.explicit.txt --yes
conda run --prefix .runtime/figures python -m pip install --no-deps --only-binary=:all: --require-hashes --target .runtime/figure2 -r environments/figure2-pip-win64.lock.txt
```

The first command restores exact Conda builds/URLs and package hashes, including Python 3.13.9, NumPy 2.3.5, pandas 2.3.3, Matplotlib 3.10.6 and FreeType 2.13.3 for Figures 3–6. The second installs the hashed Matplotlib 3.10.0 Windows wheel into a separate target for Figure 2, with FreeType 2.6.1. It does not overwrite the Conda plotting library or the analysis environment.

Verified typography uses Windows Arial, SHA-256 `c9b76220a5be42ead4733611e417cd65c5fd8aeaa33eb56576ac378a37d130a1`. Use the matching font from a licensed Windows installation; fonts are not redistributed. The wrapper checks its hash, Python/library versions and FreeType build before drawing. It stops on a mismatch rather than claiming pixel equivalence. Other platforms or fallback fonts are outside this pixel-verified profile.

## Inputs and commands

Extract Additional file 2 into a sibling `additional_file_2` directory, keeping its `frozen_inputs/`, `review_inputs/`, `historical_tables/` and `context/` together. Outputs/logs go under `--output`; inputs are not changed.

Generate all five figures:

```text
conda run --prefix .runtime/figures python figure_code/reproduce_figures.py --data ../additional_file_2 --output results/figures
```

Generate each figure separately:

```text
conda run --prefix .runtime/figures python figure_code/reproduce_figures.py --data ../additional_file_2 --output results/fig2 --figures 2
conda run --prefix .runtime/figures python figure_code/reproduce_figures.py --data ../additional_file_2 --output results/fig3 --figures 3
conda run --prefix .runtime/figures python figure_code/reproduce_figures.py --data ../additional_file_2 --output results/fig4 --figures 4
conda run --prefix .runtime/figures python figure_code/reproduce_figures.py --data ../additional_file_2 --output results/fig5 --figures 5
conda run --prefix .runtime/figures python figure_code/reproduce_figures.py --data ../additional_file_2 --output results/fig6 --figures 6
```

If the Figure 2 wheel target is elsewhere, pass `--figure2-runtime <target-directory>`. The wrapper starts each profile in a fresh subprocess and records versions/font hash in `figure_runtime_record.json`. It temporarily copies six real frozen tables into a ZIP under the output directory to use the unchanged Figure 2 archive reader, records their SHA-256 and removes the container afterward. This is not a replacement data set or a submission attachment.

The underlying entry points remain available for content-oriented rendering:

```text
python figure_code/generate_figure2.py --data-archive ../Additional_file_2.zip --output results/fig2
python figure_code/generate_figures.py --data ../additional_file_2 --output results/fig3 --figures 3
```

These direct commands do not check the pixel profile. Version numbers alone do not specify font/FreeType builds; use the wrapper for the verified profile.

## Verification scope

`environments/figure-reference.json` records native drawing-stream hashes and 300-DPI RGB reference hashes. Pixel validation uses the actual figure objects in the stated manuscript PDF under matching PDF transparency-group conditions; images are not resized or registered to hide differences. These commands draw frozen/archived results and do not refit the upstream models or all permutations.
