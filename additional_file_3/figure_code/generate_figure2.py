import os
from pathlib import Path
import math
import textwrap

import sys
if os.name == 'nt':
    _lib = Path(sys.prefix) / 'Library/bin'
    os.environ['PATH'] = str(_lib) + os.pathsep + sys.prefix + os.pathsep + os.environ.get('PATH', '')
    _dll_handle = os.add_dll_directory(str(_lib)) if _lib.exists() else None
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parent
MANUSCRIPT = Path(os.environ.get("MANUSCRIPT_DIR", REPO_ROOT / "Evidence_traceable_BioDataMining"))
FIG_DIR = MANUSCRIPT / "figures"
SOURCE_DIR = FIG_DIR / "source"

MODULES = {
    "PDAC": PROJECT_ROOT / "modules" / "PDAC" / "results" / "tables",
    "LUAD": PROJECT_ROOT / "modules" / "LUAD" / "results" / "tables",
}

CLASS_ORDER = [
    "Class I",
    "Class II",
    "Class III",
    "Class IV",
    "Class V",
    "Class VI",
]

CLASS_COLORS = {
    "Class I": "#0F4D92",
    "Class II": "#E28E2C",
    "Class III": "#4D8A64",
    "Class IV": "#B64342",
    "Class V": "#8069A8",
    "Class VI": "#8F8F8F",
}

AXIS_COLORS = {
    "translational": "#0F4D92",
    "functional": "#E28E2C",
    "baseline": "#767676",
}

NATURE = {
    "ink": "#272727",
    "muted": "#606060",
    "line": "#D8D8D8",
    "panel": "#F7F8FA",
    "blue_soft": "#DCEBFA",
    "orange_soft": "#F0E0D0",
    "green_soft": "#DDF3DE",
    "red_soft": "#F6CFCB",
    "teal": "#42949E",
}

BLOCK_COLUMNS = {
    "Patient molecular": ["patient_alteration_score", "tcga_gtex_expression_score", "tcga_mutation_score", "tcga_cnv_score", "tcga_survival_score"],
    "Protein-level": ["proteomic_support_score", "cptac_mrna_score", "cptac_protein_score", "phosphoproteomic_support_score"],
    "Functional genomics": ["dependency_score"],
    "Drug knowledge": ["chemical_modulation_score", "clinical_precedence_score"],
    "Tractability": ["target_tractability_score", "druggability_score"],
    "Safety proxy": ["safety_factor", "safety_risk_score"],
    "Reproducibility support": ["external_reproducibility_score", "geo_external_reproducibility_score", "opentargets_overall_score"],
}

matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 8.2,
        "axes.titlesize": 9.0,
        "axes.labelsize": 8.2,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.4,
        "figure.dpi": 300,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.75,
        "axes.edgecolor": "#303030",
        "xtick.color": "#303030",
        "ytick.color": "#303030",
        "axes.labelcolor": "#303030",
    }
)


def ensure_dirs() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)


def read_table(cancer: str, name: str) -> pd.DataFrame:
    return pd.read_csv(MODULES[cancer] / name, sep="\t")


def class_short(value: object) -> str:
    text = str(value)
    return text.split(":", 1)[0] if ":" in text else text


def save_source(df: pd.DataFrame, name: str) -> None:
    df.to_csv(SOURCE_DIR / name, index=False)


def save_fig(fig: plt.Figure, name: str) -> None:
    fig.savefig(FIG_DIR / f"{name}.svg")
    fig.savefig(FIG_DIR / f"{name}.pdf")
    fig.savefig(FIG_DIR / f"{name}.png", dpi=600)
    plt.close(fig)


def observed_feature_fraction(matrix: pd.DataFrame, cols: list[str]) -> float:
    """Return the fraction of rows with processed evidence coverage.

    Missingness indicators are preferred when available. If no indicator is present,
    positive numeric values are used as a conservative proxy for observed support
    rather than treating imputed zeros as observed records.
    """
    observed_parts = []
    for col in cols:
        if col not in matrix.columns:
            continue
        missing_col = f"{col}_missing"
        if missing_col in matrix.columns:
            observed_parts.append(pd.to_numeric(matrix[missing_col], errors="coerce").fillna(1).eq(0))
        else:
            observed_parts.append(pd.to_numeric(matrix[col], errors="coerce").fillna(0).gt(0))
    if not observed_parts:
        return 0.0
    observed = pd.concat(observed_parts, axis=1).any(axis=1)
    return float(observed.mean())


def panel(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.075, 1.035, letter, transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom", color=NATURE["ink"])


def add_grid(ax: plt.Axes, axis: str = "both") -> None:
    ax.grid(axis=axis, color=NATURE["line"], linewidth=0.45, alpha=0.78)
    ax.set_axisbelow(True)


def point_value_label(ax: plt.Axes, x: float, y: float, text: str, fontsize: float = 7.4) -> None:
    ax.annotate(
        text,
        (x, y),
        xytext=(0, 7),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=fontsize,
        bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
        zorder=5,
    )


def score_class_data(cancer: str) -> pd.DataFrame:
    scores = read_table(cancer, "framework_priority_scores.tsv")
    review = read_table(cancer, "framework_review_classes.tsv")
    df = scores.merge(review[["gene_symbol", "review_class"]], on="gene_symbol", how="left")
    df["class_short"] = df["review_class"].map(class_short)
    return df


def figure_2() -> None:
    summaries = []
    coverage_rows = []
    positive_rows = []
    for cancer in ["PDAC", "LUAD"]:
        summary = read_table(cancer, "framework_results_summary.tsv").set_index("metric")["value"]
        summaries.extend(
            [
                {"cancer": cancer, "quantity": "Candidate universe", "count": int(summary["candidate_universe"])},
                {"cancer": cancer, "quantity": "Translational L", "count": int(summary["translational_positives"])},
                {"cancer": cancer, "quantity": "Translational V", "count": int(summary["translational_benchmark_targets"])},
                {"cancer": cancer, "quantity": "Functional L", "count": int(summary["functional_positives"])},
                {"cancer": cancer, "quantity": "Functional V", "count": int(summary["functional_benchmark_targets"])},
            ]
        )
        matrix_name = "source_pdac_target_evidence_matrix.tsv" if cancer == "PDAC" else "source_luad_target_evidence_matrix.tsv"
        matrix = read_table(cancer, matrix_name)
        for block, cols in BLOCK_COLUMNS.items():
            present = [c for c in cols if c in matrix.columns]
            coverage = observed_feature_fraction(matrix, present)
            coverage_rows.append({"cancer": cancer, "evidence_block": block, "coverage": coverage})
        labels = read_table(cancer, "framework_record_role_manifest.tsv")
        positive = labels[labels["record_role"].isin(["L", "V"])]
        grouped = positive.groupby(["axis", "record_role", "source_group"]).size().reset_index(name="count")
        grouped["cancer"] = cancer
        positive_rows.append(grouped)

    counts = pd.DataFrame(summaries)
    coverage = pd.DataFrame(coverage_rows)
    positives = pd.concat(positive_rows, ignore_index=True)
    save_source(counts, "fig2a_candidate_and_label_counts.csv")
    save_source(coverage, "fig2b_evidence_coverage.csv")
    save_source(positives, "fig2c_positive_label_sources.csv")

    fig = plt.figure(figsize=(7.4, 5.25), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], width_ratios=[1.0, 1.05])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]

    ax = axes[0]
    ylabels = []
    ypos = []
    xpos = []
    colors = []
    for i, cancer in enumerate(["PDAC", "LUAD"]):
        sub = counts[counts["cancer"].eq(cancer)]
        quantities = ["Candidate universe", "Translational L", "Translational V", "Functional L", "Functional V"]
        for j, quantity in enumerate(quantities):
            row = sub[sub["quantity"].eq(quantity)].iloc[0]
            ylabels.append(f"{cancer} {quantity}")
            ypos.append(i * 6 + j)
            xpos.append(row["count"])
            colors.append(["#4E79A7", "#F28E2B", "#E15759", "#59A14F", "#B07AA1"][j])
    ax.hlines(ypos, [1] * len(ypos), xpos, color="#BBBBBB", linewidth=1.2)
    ax.scatter(xpos, ypos, color=colors, s=26, zorder=3)
    for x, y in zip(xpos, ypos):
        point_value_label(ax, x, y, f"{int(x):,}", fontsize=7.4)
    ax.set_xscale("log")
    ax.set_yticks(ypos, ylabels)
    ax.set_ylim(max(ypos) + 0.7, -0.7)
    ax.set_xlabel("Count, log scale")
    ax.set_title("Input, development-label and benchmark scale", loc="left")
    add_grid(ax, "x")
    panel(ax, "A")

    ax = axes[1]
    heat = coverage.pivot(index="evidence_block", columns="cancer", values="coverage").loc[list(BLOCK_COLUMNS.keys()), ["PDAC", "LUAD"]]
    im = ax.imshow(heat.values, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(2), heat.columns)
    ax.set_yticks(np.arange(len(heat.index)), heat.index)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            value = heat.iloc[i, j]
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=7.2, color="white" if value >= 0.55 else "#222222")
    ax.set_title("Processed evidence coverage", loc="left")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.ax.set_ylabel("Target fraction", rotation=270, labelpad=10)
    panel(ax, "B")

    ax = axes[2]
    pivot = positives.pivot_table(index=["cancer", "axis", "record_role"], columns="source_group", values="count", aggfunc="sum", fill_value=0)
    source_order = list(pivot.sum(axis=0).sort_values(ascending=False).index)
    y = np.arange(len(pivot.index))
    left = np.zeros(len(y))
    palette = ["#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1", "#9D9D9D", "#76B7B2", "#EDC948"]
    for k, source in enumerate(source_order):
        vals = pivot[source].values
        ax.barh(y, vals, left=left, color=palette[k % len(palette)], label=source, height=0.62)
        left += vals
    ax.set_yticks(y, [f"{idx[0]} {idx[1][:5]}. {idx[2]}" for idx in pivot.index])
    ax.set_xlabel("Role-assigned records")
    ax.set_title("Development (L) and benchmark (V) record composition", loc="left")
    ax.legend(frameon=False, fontsize=7.2, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=4, columnspacing=1.0, handlelength=1.2)
    add_grid(ax, "x")
    panel(ax, "C")
    save_fig(fig, "fig2_evidence_landscape")



# Read the frozen inputs shipped with the manuscript.
import io,zipfile,argparse
def frozen_table(cancer,name):
    prefix=f'additional_file_2/frozen_inputs/{cancer}/'
    with zipfile.ZipFile(DATA_ARCHIVE) as z:
        if name=='framework_results_summary.tsv':
            f=pd.read_csv(io.BytesIO(z.read(prefix+'framework_feature_matrix.tsv')),sep='\t')
            role=pd.read_csv(io.BytesIO(z.read(prefix+'framework_record_role_manifest.tsv')),sep='\t')
            vals={'candidate_universe':len(f)}
            for axis in ['translational','functional']:
                for code,suffix in [('L','positives'),('V','benchmark_targets')]:
                    vals[axis+'_'+suffix]=role.loc[(role.axis==axis)&(role.record_role==code),'gene_symbol'].nunique()
            return pd.DataFrame({'metric':list(vals),'value':list(vals.values())})
        if name.startswith('source_'):name='source_target_evidence_matrix.tsv'
        return pd.read_csv(io.BytesIO(z.read(prefix+name)),sep='\t')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--data-archive',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();DATA_ARCHIVE=a.data_archive
    FIG_DIR=a.output;SOURCE_DIR=FIG_DIR/'figure_source';ensure_dirs()
    read_table=frozen_table
    figure_2()

# Display name: Reproducibility support; historical external_* fields remain unchanged.
# PDAC: GEO cross-cohort support. LUAD: internal TCGA/DepMap composite; not independent external replication.
