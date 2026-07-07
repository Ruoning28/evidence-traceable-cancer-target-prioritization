from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CLASS_ORDER = ["Class I", "Class II", "Class III", "Class IV", "Class V", "Class VI"]
CLASS_COLORS = {
    "Class I": "#1B9E77",
    "Class II": "#2C7FB8",
    "Class III": "#A6761D",
    "Class IV": "#D95F02",
    "Class V": "#7570B3",
    "Class VI": "#8F8F8F",
}


def class_short(value: str) -> str:
    """Return the leading review-class label from a full review-class string."""
    return str(value).split(":", 1)[0]


def setup_style() -> None:
    """Apply a stable Matplotlib style for regenerated module figures."""
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def module_dir(module: str) -> Path:
    """Return the package directory for a requested module name.

    Parameters
    ----------
    module:
        ``PDAC``, ``LUAD`` or ``LUDA``. The code accepts both LUAD, the cancer
        abbreviation used in the paper, and LUDA, the package name requested for
        this submission bundle.

    Returns
    -------
    Path
        Absolute path to the module package directory.
    """
    return PROJECT_ROOT / ("LUDA" if module.upper() in {"LUAD", "LUDA"} else "PDAC")


def module_label(module: str) -> str:
    """Return the biological label printed in figures for a module name."""
    return "LUAD" if module.upper() in {"LUAD", "LUDA"} else "PDAC"


def read_tables(module: str) -> tuple[Path, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load all framework result tables required for module-level figures.

    Parameters
    ----------
    module:
        Module selector accepted by ``module_dir``.

    Returns
    -------
    tuple
        Table directory plus priority-score, review-class, feature-matrix,
        comparator-summary and permutation-test data frames.
    """
    tables = module_dir(module) / "result" / "tables"
    scores = pd.read_csv(tables / "framework_priority_scores.tsv", sep="\t")
    review = pd.read_csv(tables / "framework_review_classes.tsv", sep="\t")
    features = pd.read_csv(tables / "framework_feature_matrix.tsv", sep="\t", low_memory=False)
    comparator = pd.read_csv(tables / "framework_comparator_summary.tsv", sep="\t")
    perm = pd.read_csv(tables / "framework_label_permutation_test.tsv", sep="\t")
    return tables, scores, review, features, comparator, perm


def save(fig: plt.Figure, out_dir: Path, name: str) -> None:
    """Save a figure as PNG and PDF.

    Parameters
    ----------
    fig:
        Matplotlib figure object to write.
    out_dir:
        Directory receiving the figure files.
    name:
        File stem used for both output formats.

    Returns
    -------
    None
        Files are written to disk.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_triage(module: str, out_dir: Path, scores: pd.DataFrame, review: pd.DataFrame) -> None:
    """Create the dual-axis score landscape and review-class count figure."""
    df = scores.merge(review[["gene_symbol", "review_class"]], on="gene_symbol", how="left")
    df["class_short"] = df["review_class"].map(class_short)
    plot_df = df.sample(min(len(df), 12000), random_state=7)
    hi = df[df["class_short"].isin(["Class I", "Class II", "Class IV"])]
    plot_df = pd.concat([plot_df, hi]).drop_duplicates("gene_symbol")

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), gridspec_kw={"width_ratios": [1.45, 0.75]})
    ax = axes[0]
    for cls in CLASS_ORDER[::-1]:
        sub = plot_df[plot_df["class_short"].eq(cls)]
        ax.scatter(
            sub["translational_precedence_score"],
            sub["functional_discovery_score"],
            s=8 if cls in {"Class V", "Class VI"} else 22,
            alpha=0.13 if cls in {"Class V", "Class VI"} else 0.76,
            color=CLASS_COLORS.get(cls, "#777777"),
            linewidths=0,
            label=cls,
        )
    ax.set_xlabel("Translational-precedence score")
    ax.set_ylabel("Functional-discovery score")
    ax.set_title(f"{module} dual-axis target triage", loc="left", weight="bold")
    ax.grid(color="#E5E7EB", linewidth=0.5)
    ax.legend(frameon=False, ncol=2, loc="lower right")

    ax = axes[1]
    counts = df["class_short"].value_counts().reindex(CLASS_ORDER).fillna(0)
    y = np.arange(len(CLASS_ORDER))
    ax.barh(y, counts.values, color=[CLASS_COLORS[c] for c in CLASS_ORDER])
    ax.set_yticks(y, CLASS_ORDER)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Targets, log scale")
    ax.set_title("Review classes", loc="left", weight="bold")
    for yi, v in zip(y, counts.values):
        ax.text(max(v, 1) * 1.08, yi, f"{int(v):,}", va="center", fontsize=9)
    save(fig, out_dir, f"{module}_Figure_1_dual_axis_triage")


def figure_evidence(module: str, out_dir: Path, features: pd.DataFrame, review: pd.DataFrame) -> None:
    """Create the representative evidence heatmap for top-ranked targets."""
    df = features.merge(review[["gene_symbol", "review_class", "translational_precedence_rank", "functional_discovery_rank"]], on="gene_symbol", how="left")
    keep = (
        df.nsmallest(8, "translational_precedence_rank")["gene_symbol"].tolist()
        + [g for g in df.nsmallest(8, "functional_discovery_rank")["gene_symbol"].tolist()]
    )
    keep = list(dict.fromkeys(keep))[:14]
    cols = [
        "patient_alteration_score",
        "proteomic_support_score",
        "external_reproducibility_score",
        "dependency_score",
        "druggability_score",
        "safety_factor",
    ]
    labels = ["Patient", "Protein", "Reprod.", "Dependency", "Drug.", "Safety"]
    sub = df[df["gene_symbol"].isin(keep)].set_index("gene_symbol").reindex(keep)
    mat = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0).clip(0, 1).to_numpy()
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    im = ax.imshow(mat, vmin=0, vmax=1, cmap="YlGnBu", aspect="auto")
    ax.set_xticks(np.arange(len(cols)), labels, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(keep)), keep)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8, color="white" if mat[i, j] > 0.72 else "#111111")
    ax.set_title(f"{module} representative evidence cards", loc="left", weight="bold")
    fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02, label="Feature score")
    save(fig, out_dir, f"{module}_Figure_2_evidence_heatmap")


def figure_validation(module: str, out_dir: Path, comparator: pd.DataFrame, perm: pd.DataFrame) -> None:
    """Create comparator recovery and permutation-validation panels."""
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))
    comp = comparator[comparator["top_k"].eq(50)].copy()
    comp = comp[comp["method"].isin([
        "Dual-axis sampled-background PU",
        "Combined single-axis sampled-background PU",
        "Pseudo-negative L2 logistic",
        "Repeated sampled-background L2",
        "Rank aggregation",
        "Dependency-only",
        "Druggability-only",
    ])]
    comp["label"] = comp["method"].str.replace(" sampled-background PU", "", regex=False).str.replace(" elastic-net logistic", "", regex=False)
    ax = axes[0]
    for i, axis_name in enumerate(["translational", "functional"]):
        sub = comp[comp["axis"].eq(axis_name)].sort_values("positive_recovery")
        ax.barh(
            np.arange(len(sub)) + i * 0.38,
            sub["positive_recovery"],
            height=0.35,
            label=axis_name,
            alpha=0.86,
        )
        if i == 0:
            ax.set_yticks(np.arange(len(sub)) + 0.19, sub["label"])
    ax.set_xlabel("Held-out benchmark recovery at top 50")
    ax.set_title("Comparator benchmark recovery", loc="left", weight="bold")
    ax.legend(frameon=False)

    ax = axes[1]
    x = np.arange(len(perm))
    ax.bar(x - 0.18, perm["observed_value"], width=0.35, label="Observed")
    ax.bar(x + 0.18, perm["permutation_mean"], width=0.35, label="Permutation mean")
    ax.set_xticks(x, perm["axis"])
    ax.set_ylabel("Top-50 recovery")
    ax.set_title("Held-out benchmark permutation check", loc="left", weight="bold")
    for xi, row in perm.iterrows():
        ax.text(xi, max(row["observed_value"], row["permutation_mean"]) + 0.02, f"p={row['empirical_p_value']:.3f}", ha="center", fontsize=9)
    ax.legend(frameon=False)
    save(fig, out_dir, f"{module}_Figure_3_validation")


def main() -> None:
    """Command-line entry point for regenerating one module's figures."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", required=True, choices=["PDAC", "LUAD", "LUDA"])
    args = parser.parse_args()
    setup_style()
    _, scores, review, features, comparator, perm = read_tables(args.module)
    label = module_label(args.module)
    out_dir = module_dir(args.module) / "result" / "figures"
    figure_triage(label, out_dir, scores, review)
    figure_evidence(label, out_dir, features, review)
    figure_validation(label, out_dir, comparator, perm)
    print(f"Wrote figures to {out_dir}")


if __name__ == "__main__":
    main()

