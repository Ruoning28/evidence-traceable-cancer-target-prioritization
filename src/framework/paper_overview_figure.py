from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 600,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parent
MANUSCRIPT_FIGURES = REPO_ROOT / "Evidence_traceable_BioDataMining" / "figures"
SUBMISSION_FIGURES = PROJECT_ROOT / "result" / "figures"

MODULES = {
    "PDAC": PROJECT_ROOT / "PDAC",
    "LUAD": PROJECT_ROOT / "LUDA",
}

CLASS_KEYS = ["I", "II", "III", "IV", "V", "VI"]
CLASS_LABELS = [
    "I: translational\nprecedence",
    "II: functional\ndiscovery",
    "III: dependency\nlimited",
    "IV: safety-risk",
    "V: biomarker-like\nlow dependency",
    "VI: insufficient\nor mixed",
]


def read_table(cancer: str, name: str) -> pd.DataFrame:
    """Return a result TSV for ``cancer`` as a pandas DataFrame.

    Parameters
    ----------
    cancer:
        Cancer key, either ``PDAC`` or ``LUAD``.
    name:
        TSV file name under the cancer module's ``result/tables`` directory.

    Returns
    -------
    pandas.DataFrame
        Parsed table contents.
    """
    return pd.read_csv(MODULES[cancer] / "result" / "tables" / name, sep="\t")


def summary(cancer: str) -> pd.Series:
    """Return the metric summary vector for ``cancer``.

    Parameters
    ----------
    cancer:
        Cancer key, either ``PDAC`` or ``LUAD``.

    Returns
    -------
    pandas.Series
        Series indexed by metric name, containing numeric values.
    """
    return read_table(cancer, "framework_results_summary.tsv").set_index("metric")["value"].astype(float)


def permutation_summary() -> pd.DataFrame:
    """Return the combined label-permutation validation summary.

    Returns
    -------
    pandas.DataFrame
        One row per cancer and axis, restricted to the top-50 validation metric.
    """
    frames = []
    for cancer in MODULES:
        df = read_table(cancer, "framework_label_permutation_test.tsv")
        df = df[(df["metric"].eq("top50_recovery")) & (df["evaluation_set"].eq("independent_benchmark"))].copy()
        df["cancer"] = cancer
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def add_box(ax, x: float, y: float, w: float, h: float, text: str, fc: str, ec: str = "#2f3b46", size: int = 9) -> None:
    """Draw a rounded text box and return no value.

    Parameters
    ----------
    ax:
        Matplotlib axes receiving the box.
    x, y, w, h:
        Box position and size in axes coordinates.
    text:
        Box label.
    fc:
        Fill color.
    ec:
        Edge color.
    size:
        Font size for the box label.
    """
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.1,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, color="#1d2730")


def fmt_int(value: float) -> str:
    """Return ``value`` as a comma-grouped integer string."""
    return f"{int(round(float(value))):,}"


def draw_overview() -> plt.Figure:
    """Create the Figure 1 overview from current result tables.

    Returns
    -------
    matplotlib.figure.Figure
        Fully rendered overview figure.
    """
    pdac = summary("PDAC")
    luad = summary("LUAD")
    perm = permutation_summary()

    fig = plt.figure(figsize=(14, 8), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax.text(0.04, 0.955, "Evidence-traceable target prioritization workflow", fontsize=19, weight="bold", color="#1d2730")
    ax.text(0.04, 0.915, "Role-separated inputs, dual-axis sampled-background PU-style ranking, fixed review gates, and auditable validation outputs", fontsize=10.5, color="#4d5b66")

    add_box(ax, 0.045, 0.745, 0.22, 0.115, f"PDAC module\n{fmt_int(pdac['candidate_universe'])} genes\n7 T-dev / 21 F-dev\n11 T-val / 30 F-val", "#d9ebf7", size=9.5)
    add_box(ax, 0.045, 0.585, 0.22, 0.115, f"LUAD module\n{fmt_int(luad['candidate_universe'])} genes\n6 T-dev / 26 F-dev\n7 T-val / 36 F-val", "#e7f1df", size=9.5)
    add_box(ax, 0.315, 0.665, 0.22, 0.12, "Source-role assignment\nF: features\nL: development labels\nV: validation benchmarks\nA: annotation/review", "#f6edd8", size=9.2)
    add_box(ax, 0.585, 0.665, 0.22, 0.12, "Dual-axis modelling\nTranslational-precedence\nFunctional-discovery\n30 sampled backgrounds", "#eadff1", size=9.2)

    for x1, y1, x2, y2 in [(0.265, 0.802, 0.315, 0.725), (0.265, 0.642, 0.315, 0.725), (0.535, 0.725, 0.585, 0.725)]:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=1.2, color="#3e4a54"))

    add_box(ax, 0.835, 0.655, 0.12, 0.14, "Fixed review gates\nDependency\nCancer relevance\nTractability\nSafety proxies", "#f2dfdf", size=8.7)
    ax.annotate("", xy=(0.835, 0.725), xytext=(0.805, 0.725), arrowprops=dict(arrowstyle="->", lw=1.2, color="#3e4a54"))

    ax.text(0.05, 0.50, "Review-class output", fontsize=13, weight="bold", color="#1d2730")
    y0 = 0.445
    colors = ["#0F4D92", "#E28E2C", "#4D8A64", "#B64342", "#8069A8", "#8F8F8F"]
    max_count = max([pdac[f"class_{k}"] for k in CLASS_KEYS] + [luad[f"class_{k}"] for k in CLASS_KEYS])
    for idx, (key, label) in enumerate(zip(CLASS_KEYS, CLASS_LABELS)):
        y = y0 - idx * 0.055
        p = pdac[f"class_{key}"]
        l = luad[f"class_{key}"]
        ax.text(0.05, y + 0.011, f"Class {key}", fontsize=9, weight="bold", color="#1d2730")
        ax.text(0.115, y + 0.011, label, fontsize=7.4, color="#4d5b66", va="center")
        ax.barh(y + 0.016, 0.22 * p / max_count, height=0.014, left=0.255, color=colors[idx], alpha=0.85)
        ax.barh(y - 0.002, 0.22 * l / max_count, height=0.014, left=0.255, color=colors[idx], alpha=0.45)
        ax.text(0.49, y + 0.016, f"PDAC {fmt_int(p)}", fontsize=7.7, va="center", color="#1d2730")
        ax.text(0.49, y - 0.002, f"LUAD {fmt_int(l)}", fontsize=7.7, va="center", color="#1d2730")

    ax.text(0.61, 0.50, "Independent validation and comparator interpretation", fontsize=13, weight="bold", color="#1d2730")
    validation_lines = []
    for cancer, axis in [("LUAD", "translational"), ("PDAC", "translational"), ("PDAC", "functional"), ("LUAD", "functional")]:
        row = perm[(perm["cancer"].eq(cancer)) & (perm["axis"].eq(axis))].iloc[0]
        pval = float(row["empirical_p_value"])
        validation_lines.append(
            f"{cancer} {axis}: top-50 {float(row['observed_value']):.3f} vs permutation mean {float(row['permutation_mean']):.4f}; P={pval:.4f}"
        )
    ax.text(0.61, 0.455, "\n".join(validation_lines), fontsize=9.0, color="#1d2730", linespacing=1.5, va="top")

    add_box(
        ax,
        0.61,
        0.185,
        0.345,
        0.13,
        "Comparator result summary\nRank aggregation recovered more translational benchmarks in top-50,\nwhereas dual-axis scoring preserved objective-specific interpretation\nand linked rankings to post-model evidence-card review.",
        "#f4f4f4",
        size=8.2,
    )

    ax.text(0.05, 0.085, "Figure 1 values are generated from the current result tables; T-dev/F-dev denote development positives and T-val/F-val denote held-out benchmark targets.", fontsize=8.4, color="#5a6670")
    return fig


def main() -> None:
    """Write Figure 1 to manuscript and submission result directories."""
    MANUSCRIPT_FIGURES.mkdir(parents=True, exist_ok=True)
    SUBMISSION_FIGURES.mkdir(parents=True, exist_ok=True)
    fig = draw_overview()
    for directory in [MANUSCRIPT_FIGURES, SUBMISSION_FIGURES]:
        for suffix in ["pdf", "png", "svg"]:
            fig.savefig(directory / f"fig1_evidence_traceable_framework.{suffix}", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote Figure 1 to {MANUSCRIPT_FIGURES} and {SUBMISSION_FIGURES}")


if __name__ == "__main__":
    main()

