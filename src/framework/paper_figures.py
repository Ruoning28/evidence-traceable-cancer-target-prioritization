from pathlib import Path
import math
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = PROJECT_ROOT / "result" / "paper_figures"
SOURCE_DIR = FIG_DIR / "source"

MODULES = {
    "PDAC": PROJECT_ROOT / "PDAC" / "result" / "tables",
    "LUAD": PROJECT_ROOT / "LUDA" / "result" / "tables",
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
    "External reproducibility": ["external_reproducibility_score", "geo_external_reproducibility_score", "opentargets_overall_score"],
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
    """Create manuscript figure and source-data directories."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)


def read_table(cancer: str, name: str) -> pd.DataFrame:
    """Return TSV ``name`` from the selected cancer module."""
    return pd.read_csv(MODULES[cancer] / name, sep="\t")


def class_short(value: object) -> str:
    """Return the review-class identifier contained in ``value``."""
    text = str(value)
    return text.split(":", 1)[0] if ":" in text else text


def save_source(df: pd.DataFrame, name: str) -> None:
    """Write figure source ``df`` to CSV file ``name``."""
    df.to_csv(SOURCE_DIR / name, index=False)


def save_fig(fig: plt.Figure, name: str) -> None:
    """Save ``fig`` as PDF and high-resolution PNG under ``name``."""
    fig.savefig(FIG_DIR / f"{name}.svg")
    fig.savefig(FIG_DIR / f"{name}.pdf")
    fig.savefig(FIG_DIR / f"{name}.png", dpi=600)
    plt.close(fig)


def panel(ax: plt.Axes, letter: str) -> None:
    """Place panel ``letter`` relative to matplotlib axis ``ax``."""
    ax.text(-0.075, 1.035, letter, transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom", color=NATURE["ink"])


def add_grid(ax: plt.Axes, axis: str = "both") -> None:
    """Add a restrained grid on the requested axis of ``ax``."""
    ax.grid(axis=axis, color=NATURE["line"], linewidth=0.45, alpha=0.78)
    ax.set_axisbelow(True)


def point_value_label(ax: plt.Axes, x: float, y: float, text: str, fontsize: float = 7.4) -> None:
    """Annotate point ``(x, y)`` with ``text`` and return no value."""
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
    """Return merged score and review-class records for ``cancer``."""
    scores = read_table(cancer, "framework_priority_scores.tsv")
    review = read_table(cancer, "framework_review_classes.tsv")
    df = scores.merge(review[["gene_symbol", "review_class"]], on="gene_symbol", how="left")
    df["class_short"] = df["review_class"].map(class_short)
    return df


def observed_feature_fraction(matrix: pd.DataFrame, cols: list[str]) -> float:
    """Return the fraction of rows with source-observed evidence.

    Parameters
    ----------
    matrix:
        Source evidence matrix before manuscript-level plotting.
    cols:
        Feature columns belonging to one evidence block.

    Returns
    -------
    float
        Fraction of candidate records with at least one observed feature in the
        requested block. Missingness indicators are preferred because many
        feature columns are zero-imputed before modelling.
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


def figure_2() -> None:
    """Generate the input-scale, coverage and L/V composition figure."""
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
    ax.set_title("Observed evidence coverage", loc="left")
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


def figure_3() -> None:
    """Generate the LUAD development/benchmark anchor validation figure."""
    df = score_class_data("LUAD")
    cards = read_table("LUAD", "framework_evidence_cards.tsv")
    role_manifest = read_table("LUAD", "framework_record_role_manifest.tsv")
    anchor_roles = role_manifest[(role_manifest["axis"].eq("translational")) & (role_manifest["record_role"].isin(["L", "V"]))][["gene_symbol", "record_role"]].drop_duplicates()
    anchors = anchor_roles["gene_symbol"].tolist()
    benchmark_anchors = anchor_roles[anchor_roles["record_role"].eq("V")]["gene_symbol"].tolist()
    perm = read_table("LUAD", "framework_label_permutation_test.tsv")
    comp = read_table("LUAD", "framework_comparator_summary.tsv")
    anchor_cards = cards[cards["gene_symbol"].isin(anchors)].merge(anchor_roles, on="gene_symbol", how="left")
    save_source(anchor_cards, "fig3_luad_anchor_ranks.csv")

    fig = plt.figure(figsize=(7.4, 5.35), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.28], width_ratios=[1.18, 1.0])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]

    ax = axes[0]
    anchor_cards = anchor_cards.sort_values(["record_role", "translational_precedence_rank"], ascending=[False, True])
    y = np.arange(len(anchor_cards))
    ax.hlines(y, anchor_cards["translational_precedence_rank"], anchor_cards["functional_discovery_rank"], color="#C8C8C8", linewidth=1.3)
    ax.scatter(anchor_cards["translational_precedence_rank"], y, color=AXIS_COLORS["translational"], s=28, label="Translational rank", zorder=3)
    ax.scatter(anchor_cards["functional_discovery_rank"], y, color=AXIS_COLORS["functional"], s=28, label="Functional rank", zorder=3)
    ax.set_xscale("log")
    ax.set_yticks(y, anchor_cards["gene_symbol"])
    ax.invert_yaxis()
    ax.set_xlabel("Rank, lower is better")
    ax.set_title("Development + held-out anchor ranks", loc="left")
    add_grid(ax, "x")
    ax.legend(frameon=False, loc="upper right", fontsize=7.4)
    panel(ax, "A")

    ax = axes[1]
    dual = comp[(comp["method"].eq("Dual-axis sampled-background PU")) & (comp["axis"].eq("translational"))]
    ax.plot(dual["top_k"], dual["positive_recovery"], marker="o", color=AXIS_COLORS["translational"], linewidth=1.7, label="Held-out anchor recovery")
    mean = float(perm[perm["axis"].eq("translational")]["permutation_mean"].iloc[0])
    ax.axhline(mean, color=AXIS_COLORS["baseline"], linestyle="--", linewidth=1.0, label="Permutation mean")
    ax.set_ylim(-0.02, 1.08)
    ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
    ax.set_xlabel("Top-K cutoff")
    ax.set_ylabel("Anchor recovery")
    ax.set_title("Held-out top-K anchor recovery", loc="left")
    row50 = dual[dual["top_k"].eq(50)].iloc[0]
    recovered = int(round(row50["positive_recovery"] * len(benchmark_anchors)))
    ax.annotate(f"{recovered}/{len(benchmark_anchors)} held-out anchors", xy=(50, row50["positive_recovery"]), xytext=(58, min(1.0, row50["positive_recovery"] + 0.18)), textcoords="data", fontsize=7.2, arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.6})
    add_grid(ax)
    ax.legend(frameon=False, loc="lower right", fontsize=7.0)
    panel(ax, "B")

    ax = axes[2]
    bg = df.sample(min(7000, len(df)), random_state=7)
    ax.scatter(bg["translational_precedence_score"], bg["functional_discovery_score"], s=4, color="#CFCFCF", alpha=0.22, linewidths=0, label="Background genes")
    for cls in ["Class I", "Class II", "Class IV"]:
        sub = df[df["class_short"].eq(cls)]
        ax.scatter(sub["translational_precedence_score"], sub["functional_discovery_score"], s=9, color=CLASS_COLORS[cls], alpha=0.42, linewidths=0, label=cls)
    a = df[df["gene_symbol"].isin(anchors)].merge(anchor_roles, on="gene_symbol", how="left")
    for role, marker, label in [("L", "o", "Development anchors"), ("V", "s", "Held-out anchors")]:
        sub = a[a["record_role"].eq(role)]
        ax.scatter(sub["translational_precedence_score"], sub["functional_discovery_score"], s=36, marker=marker, facecolor="white", edgecolor="#1F1F1F", linewidth=0.85, zorder=5, label=label)
    label_pos = {
        "EGFR": (0.85, 1.16),
        "MET": (0.80, 1.06),
        "BRAF": (1.06, 1.16),
        "ALK": (1.13, 1.06),
        "ERBB2": (0.81, 0.83),
        "KRAS": (0.88, 0.76),
        "NTRK1": (1.05, 0.86),
        "RET": (1.13, 0.84),
        "ROS1": (1.13, 0.71),
        "NTRK3": (1.00, 0.50),
        "NTRK2": (0.98, 0.22),
    }
    for _, row in a.iterrows():
        tx, ty = label_pos.get(
            row["gene_symbol"],
            (row["translational_precedence_score"] + 0.02, row["functional_discovery_score"] + 0.02),
        )
        ax.annotate(
            row["gene_symbol"],
            (row["translational_precedence_score"], row["functional_discovery_score"]),
            xytext=(tx, ty),
            textcoords="data",
            fontsize=6.8,
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
            arrowprops={"arrowstyle": "-", "lw": 0.45, "color": "#555555"},
            zorder=6,
            annotation_clip=False,
        )
    ax.set_xlim(-0.05, 1.20)
    ax.set_ylim(-0.05, 1.24)
    ax.set_xlabel("$s_{trans}$")
    ax.set_ylabel("$s_{func}$")
    ax.set_title("LUAD score landscape", loc="left")
    add_grid(ax)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=5, fontsize=7.0, columnspacing=0.9, handletextpad=0.35)
    panel(ax, "C")
    save_fig(fig, "fig3_luad_positive_control")


def figure_4() -> None:
    """Generate the PDAC triage and representative-target figure."""
    df = score_class_data("PDAC")
    review = read_table("PDAC", "framework_review_classes.tsv")
    cards = read_table("PDAC", "framework_evidence_cards.tsv")
    # The first three examples illustrate prespecified interpretive roles, whereas
    # the safety-constrained example is selected from the current Class IV evidence
    # cards so that it follows the experimental output rather than a fixed gene name.
    selected = pdac_representative_genes(cards)
    mat_cols = ["patient_alteration", "proteomic_support", "external_reproducibility", "dependency", "druggability", "safety_factor", "safety_risk", "open_targets"]
    selected_cards = cards[cards["gene_symbol"].isin(selected)].copy()
    missing = sorted(set(selected) - set(selected_cards["gene_symbol"].astype(str)))
    if missing:
        raise ValueError(f"Missing PDAC representative evidence-card rows: {', '.join(missing)}")
    selected_cards["class_short"] = selected_cards["review_class"].map(class_short)
    save_source(selected_cards, "fig4_pdac_representative_cards.csv")

    fig = plt.figure(figsize=(7.4, 5.45), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.0, 1.0])
    axes = [fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])]

    ax = axes[0]
    selected_cards = selected_cards.set_index("gene_symbol").reindex(selected).reset_index()
    mat = selected_cards[mat_cols].astype(float).values
    im = ax.imshow(mat, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(mat_cols)), [c.replace("_", " ") for c in mat_cols], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(selected)), selected_cards["gene_symbol"])
    for i in range(mat.shape[0]):
        cls = selected_cards.loc[i, "class_short"]
        cls = "Class VI" if pd.isna(cls) else str(cls)
        ax.text(len(mat_cols) + 0.08, i, cls.replace("Class ", ""), va="center", fontsize=7.2, color=CLASS_COLORS.get(cls, "#333333"))
    ax.set_xlim(-0.5, len(mat_cols) + 1.0)
    ax.set_title("Representative evidence matrix", loc="left")
    cbar = fig.colorbar(im, ax=ax, fraction=0.036, pad=0.015)
    cbar.ax.set_ylabel("Scaled evidence", rotation=270, labelpad=8)
    panel(ax, "A")

    ax = axes[1]
    x = df["translational_precedence_score"].clip(lower=0)
    y = df["functional_discovery_score"].clip(lower=0)
    hb = ax.hexbin(x, y, gridsize=34, mincnt=1, bins="log", cmap="Greys", linewidths=0, alpha=0.75)
    hi = df[df["gene_symbol"].isin(selected)]
    ax.scatter(hi["translational_precedence_score"], hi["functional_discovery_score"], s=25, color="#E15759", edgecolor="white", linewidth=0.45, zorder=4)
    label_pos = {
        "KIF23": (0.24, 1.10),
        "KRAS": (0.70, 1.16),
        "ARID1A": (1.10, 0.46),
    }
    if len(selected) >= 4:
        label_pos.setdefault(selected[3], (1.10, 0.97))
    for _, row in hi.iterrows():
        tx, ty = label_pos.get(
            row["gene_symbol"],
            (row["translational_precedence_score"] + 0.02, row["functional_discovery_score"] + 0.02),
        )
        ax.annotate(
            row["gene_symbol"],
            (row["translational_precedence_score"], row["functional_discovery_score"]),
            xytext=(tx, ty),
            textcoords="data",
            fontsize=6.8,
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
            arrowprops={"arrowstyle": "-", "lw": 0.45, "color": "#555555"},
            zorder=5,
            annotation_clip=False,
        )
    ax.set_xlim(-0.05, 1.12)
    ax.set_ylim(0.0, 1.22)
    ax.set_xlabel("$s_{trans}$")
    ax.set_ylabel("$s_{func}$")
    ax.set_title("PDAC score density", loc="left")
    panel(ax, "B")

    ax = axes[2]
    counts = review["review_class"].map(class_short).value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)
    total = counts.sum()
    y = np.arange(len(CLASS_ORDER))
    ax.scatter(counts.values, y, color=[CLASS_COLORS[c] for c in CLASS_ORDER], s=28)
    ax.hlines(y, 0, counts.values, color="#CFCFCF", linewidth=1.0)
    for yi, val in zip(y, counts.values):
        point_value_label(ax, val, yi, f"{val:,} ({val / total:.1%})", fontsize=7.0)
    ax.set_xscale("log")
    ax.set_yticks(y, CLASS_ORDER)
    ax.set_ylim(max(y) + 0.7, -0.7)
    ax.set_xlabel("Targets, log scale")
    ax.set_title("Review classes", loc="left")
    add_grid(ax, "x")
    panel(ax, "C")
    save_fig(fig, "fig4_pdac_triage")


def figure_5() -> None:
    """Generate permutation, OOF, source-out and sensitivity panels."""
    perm_rows = []
    recovery_rows = []
    oof_rows = []
    source_out_rows = []
    stability_rows = []
    for cancer in ["PDAC", "LUAD"]:
        perm = read_table(cancer, "framework_label_permutation_test.tsv")
        perm["cancer"] = cancer
        perm_rows.append(perm)
        comp = read_table(cancer, "framework_comparator_summary.tsv")
        dual = comp[comp["method"].eq("Dual-axis sampled-background PU")].copy()
        dual["cancer"] = cancer
        recovery_rows.append(dual)
        oof = read_table(cancer, "framework_oof_validation_summary.tsv")
        oof["cancer"] = cancer
        oof_rows.append(oof)
        source_out = read_table(cancer, "framework_leave_source_out_validation.tsv")
        source_out["cancer"] = cancer
        source_out_rows.append(source_out)
        sens = read_table(cancer, "framework_model_sensitivity.tsv")
        sens["cancer"] = cancer
        stability_rows.append(sens)
    perm_all = pd.concat(perm_rows, ignore_index=True)
    recovery = pd.concat(recovery_rows, ignore_index=True)
    oof = pd.concat(oof_rows, ignore_index=True)
    source_out = pd.concat(source_out_rows, ignore_index=True)
    stability = pd.concat(stability_rows, ignore_index=True)
    ordered = perm_all.sort_values("empirical_p_value").copy()
    adjusted = np.maximum.accumulate(
        np.minimum(1.0, ordered["empirical_p_value"].to_numpy() * np.arange(len(ordered), 0, -1))
    )
    ordered["holm_p_value"] = adjusted
    perm_all = perm_all.merge(ordered[["cancer", "axis", "holm_p_value"]], on=["cancer", "axis"], how="left")
    save_source(perm_all, "fig5ab_label_permutation_summary.csv")
    save_source(source_out, "fig5c_full_source_out_transfer.csv")
    save_source(recovery, "fig5d_independent_topk_recovery.csv")
    save_source(oof, "fig5e_oof_recovery.csv")
    save_source(stability, "fig5f_model_sensitivity.csv")

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(7.4, 6.1),
        gridspec_kw={"height_ratios": [1.05, 1.0], "width_ratios": [1.0, 1.0, 1.12]},
    )
    for ax, cancer, letter in [(axes[0, 0], "LUAD", "A"), (axes[0, 1], "PDAC", "B")]:
        sub = perm_all[perm_all["cancer"].eq(cancer)].copy()
        y = np.arange(len(sub))
        ax.errorbar(sub["permutation_mean"], y, xerr=1.96 * sub["permutation_sd"], fmt="o", color=AXIS_COLORS["baseline"], capsize=2.5, markersize=5, label="Permutation +/-1.96 SD")
        colors = [AXIS_COLORS[a] for a in sub["axis"]]
        ax.scatter(sub["observed_value"], y, color=colors, s=38, label="Observed", zorder=3)
        for i, row in sub.iterrows():
            # Keep the exact and multiplicity-adjusted p-values in a fixed
            # annotation column so they cannot collide with observed points.
            ax.text(
                0.98,
                y[list(sub.index).index(i)],
                f"p={row['empirical_p_value']:.3f}\nHolm={row['holm_p_value']:.3f}",
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=6.2,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.3},
            )
        ax.set_yticks(y, sub["axis"].str.title())
        ax.set_ylim(1.35, -0.35)
        ax.set_xlim(0, max(1.05, sub["observed_value"].max() + 0.2))
        ax.set_xlabel("Top50 recovery")
        ax.set_title(f"{cancer} permutation validation", x=0.16, ha="left")
        add_grid(ax, "x")
        panel(ax, letter)
    ax = axes[0, 2]
    estimated = source_out[source_out["status"].eq("ESTIMATED")].copy().sort_values(["cancer", "axis", "heldout_source"])
    y = np.arange(len(estimated))
    ax.scatter(estimated["median_rank_percentile"], y, color=[AXIS_COLORS[a] for a in estimated["axis"]], s=24)
    source_labels = {
        "known_direct_drug": "drug precedence",
        "known_biomarker": "biomarker",
        "known_pdac_driver": "driver",
        "known_dependency": "dependency",
        "luad_actionable_anchor": "actionable anchors",
        "tcga_luad_nature_2014": "TCGA",
        "depmap_luad_high_dependency": "DepMap",
    }
    labels = []
    for cancer, source in zip(estimated["cancer"], estimated["heldout_source"]):
        source = str(source)
        short_source = source_labels.get(source, source.replace("curated_literature:", ""))
        labels.append(f"{cancer} {short_source}")
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Median rank percentile")
    ax.set_title("Full-source-out transfer", x=0.16, ha="left")
    add_grid(ax, "x")
    panel(ax, "C")

    ax = axes[1, 0]
    for (cancer, axis), sub in recovery.groupby(["cancer", "axis"]):
        color = AXIS_COLORS[axis]
        linestyle = "-" if cancer == "LUAD" else "--"
        label = f"{cancer} {axis[:5]}."
        ax.plot(sub["top_k"], sub["positive_recovery"], marker="o", color=color, linestyle=linestyle, label=label)
    ax.set_xlabel("Top-K cutoff")
    ax.set_ylabel("Benchmark recovery")
    ax.set_title("Held-out top-K recovery", x=0.16, ha="left")
    add_grid(ax)
    shared_handles, shared_labels = ax.get_legend_handles_labels()
    panel(ax, "D")

    ax = axes[1, 1]
    oof_long = oof.melt(
        id_vars=["cancer", "axis"],
        value_vars=["top20_recovery", "top50_recovery", "top100_recovery", "top200_recovery"],
        var_name="cutoff", value_name="recovery",
    )
    oof_long["top_k"] = oof_long["cutoff"].str.extract(r"top(\d+)")[0].astype(int)
    for (cancer, axis), sub in oof_long.groupby(["cancer", "axis"]):
        ax.plot(sub["top_k"], sub["recovery"], marker="o", color=AXIS_COLORS[axis], linestyle="-" if cancer == "LUAD" else "--", label=f"{cancer} {axis[:5]}.")
    ax.set_xlabel("Top-K cutoff")
    ax.set_ylabel("OOF recovery")
    ax.set_title("Held-out development positives", x=0.16, ha="left")
    add_grid(ax)
    panel(ax, "E")

    ax = axes[1, 2]
    keep = stability[
        stability["variant"].ne("primary")
        & ~stability["variant"].str.startswith("matched_background_", na=False)
    ].copy()
    order = ["PDAC translational", "PDAC functional", "LUAD translational", "LUAD functional"]
    vals = []
    for item in order:
        cancer, axis = item.split()
        vals.append(keep[(keep["cancer"].eq(cancer)) & (keep["axis"].eq(axis))]["spearman_with_primary"].min())
    y = np.arange(len(order))
    bar_colors = ["#4E79A7", "#F28E2B", "#4E79A7", "#F28E2B"]
    ax.barh(y, vals, color=bar_colors, alpha=0.22, height=0.52)
    ax.scatter(vals, y, color=bar_colors, s=42, zorder=3)
    for xi, yi in zip(vals, y):
        ax.text(min(xi + 0.035, 0.98), yi, f"{xi:.2f}", va="center", fontsize=7.0)
    ax.set_yticks(y, order)
    ax.set_ylim(len(order) - 0.35, -0.65)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Minimum Spearman rho")
    ax.set_title("Model-parameter stability", x=0.16, ha="left")
    add_grid(ax, "x")
    panel(ax, "F")
    fig.legend(
        shared_handles,
        shared_labels,
        frameon=False,
        ncol=4,
        fontsize=6.8,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        columnspacing=0.9,
        handlelength=1.4,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1), h_pad=1.5, w_pad=0.9)
    save_fig(fig, "fig5_ranking_validation")


def short_class(text: str) -> str:
    """Return a compact review-class label for ``text``."""
    return class_short(text).replace("Class ", "C")


def wrap(text: object, width: int = 23) -> str:
    """Return ``text`` wrapped at ``width`` characters for card display."""
    if pd.isna(text) or str(text).strip() == "":
        return "No dominant limit"
    return "\n".join(textwrap.wrap(str(text), width=width))


def review_route(review_class: str) -> str:
    """Return the review route associated with ``review_class``."""
    cls = class_short(review_class)
    if cls == "Class I":
        return "Promote: translational"
    if cls == "Class II":
        return "Promote: functional"
    if cls == "Class III":
        return "Retain: limited actionability"
    if cls == "Class IV":
        return "Deprioritize: safety"
    if cls == "Class V":
        return "Retain: biomarker-like"
    return "Review: insufficient evidence"


def pdac_representative_genes(cards: pd.DataFrame) -> list[str]:
    """Return synchronized PDAC evidence-card examples for manuscript figures.

    Parameters
    ----------
    cards:
        Evidence-card table generated by the PDAC framework run.

    Returns
    -------
    list[str]
        Four gene symbols used in Fig. 4 and Fig. 6. The safety-constrained
        example is selected automatically from current Class IV evidence-card
        records, preferring safety-constrained benchmark-like records when they
        are available in the exported card table.
    """
    required = ["KRAS", "ARID1A", "KIF23"]
    available = set(cards["gene_symbol"].astype(str))
    missing_required = [gene for gene in required if gene not in available]
    if missing_required:
        raise ValueError(f"Missing required PDAC evidence-card examples: {', '.join(missing_required)}")

    class_iv = cards[cards["review_class"].astype(str).str.startswith("Class IV")].copy()
    if class_iv.empty:
        raise ValueError("No Class IV evidence-card row is available for the PDAC safety-constrained example.")

    # PALB2 is not hard-coded as the representative; it is only part of a small
    # benchmark-like preference set used when multiple Class IV records are present.
    benchmark_like = {"PALB2", "BRCA1", "BRCA2", "ATM", "MLH1", "MSH6"}
    class_iv["benchmark_like"] = class_iv["gene_symbol"].astype(str).isin(benchmark_like)
    class_iv["has_safety_limit"] = class_iv["limiting_evidence"].astype(str).str.contains("safety", case=False, na=False)
    for col in ["safety_risk", "safety_factor", "translational_precedence_score"]:
        class_iv[col] = pd.to_numeric(class_iv[col], errors="coerce")

    class_iv = class_iv.sort_values(
        ["benchmark_like", "has_safety_limit", "safety_risk", "safety_factor", "translational_precedence_score"],
        ascending=[False, False, False, True, False],
    )
    safety_gene = str(class_iv.iloc[0]["gene_symbol"])
    selected = list(dict.fromkeys(required + [safety_gene]))
    if len(selected) != 4:
        raise ValueError(f"Expected four PDAC representative genes, found {len(selected)}: {selected}")
    return selected


def figure_6() -> None:
    """Generate the target-level evidence-card manuscript figure."""
    pdac_cards = read_table("PDAC", "framework_evidence_cards.tsv")
    selections = {
        "PDAC": pdac_representative_genes(pdac_cards),
    }
    rows = []
    for cancer, genes in selections.items():
        cards = pdac_cards if cancer == "PDAC" else read_table(cancer, "framework_evidence_cards.tsv")
        sub = cards[cards["gene_symbol"].isin(genes)].copy()
        missing = sorted(set(genes) - set(sub["gene_symbol"].astype(str)))
        if missing:
            raise ValueError(f"Missing evidence-card rows for {cancer}: {', '.join(missing)}")
        sub["cancer"] = cancer
        rows.append(sub)
    cards = pd.concat(rows, ignore_index=True)
    order = [(cancer, gene) for cancer, genes in selections.items() for gene in genes]
    cards["_order"] = cards.apply(lambda r: order.index((r["cancer"], r["gene_symbol"])), axis=1)
    cards = cards.sort_values("_order")
    save_source(cards.drop(columns=["_order"]), "fig6_evidence_card_examples.csv")

    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.2))
    metrics = [
        ("patient_alteration", "Patient"),
        ("proteomic_support", "Protein"),
        ("external_reproducibility", "Repro."),
        ("dependency", "Depend."),
        ("druggability", "Drug."),
        ("safety_factor", "Safety accept."),
        ("safety_risk", "Safety risk"),
        ("open_targets", "OT"),
    ]
    for ax, (_, row) in zip(axes.ravel(), cards.iterrows()):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        cls = class_short(row["review_class"])
        edge = CLASS_COLORS.get(cls, "#9D9D9D")
        ax.add_patch(Rectangle((0.02, 0.03), 0.96, 0.94, fill=False, edgecolor=edge, linewidth=1.35))
        ax.add_patch(Rectangle((0.04, 0.045), 0.92, 0.075, facecolor=edge, edgecolor="none", alpha=0.16))
        ax.text(0.06, 0.91, f"{row['gene_symbol']} ({row['cancer']})", fontsize=10.2, fontweight="bold", color="#222222")
        ax.text(0.06, 0.825, f"{short_class(row['review_class'])}; trans rank {int(row['translational_precedence_rank'])}; func rank {int(row['functional_discovery_rank'])}", fontsize=7.2)
        y0 = 0.70
        for i, (col, label) in enumerate(metrics):
            y = y0 - i * 0.066
            value = float(row[col])
            ax.text(0.06, y, label, fontsize=6.6, va="center")
            ax.add_patch(Rectangle((0.30, y - 0.018), 0.50, 0.030, facecolor="#EFEFEF", edgecolor="none"))
            ax.add_patch(Rectangle((0.30, y - 0.018), 0.50 * max(0, min(1, value)), 0.030, facecolor="#4E79A7" if col != "safety_factor" else "#59A14F", edgecolor="none"))
            ax.text(0.83, y, f"{value:.2f}", fontsize=6.6, va="center")
        ax.text(0.06, 0.200, "Limiting evidence:", fontsize=6.8, fontweight="bold")
        ax.text(0.06, 0.155, wrap(row.get("limiting_evidence"), 38), fontsize=6.0, va="top", linespacing=0.92)
        ax.text(0.06, 0.070, review_route(row["review_class"]), fontsize=6.7, fontweight="bold", color="#222222", va="center")
    fig.tight_layout(w_pad=0.45, h_pad=0.60)
    save_fig(fig, "fig6_evidence_cards")


def figure_7() -> None:
    """Generate the review-threshold robustness supplementary figure.

    Parameters
    ----------
    None

    Returns
    -------
    None
        Writes ``fig7_threshold_robustness`` PDF/PNG and source CSV files.
    """
    rows = []
    for cancer in ["PDAC", "LUAD"]:
        sens = read_table(cancer, "framework_review_threshold_sensitivity.tsv").copy()
        sens["cancer"] = cancer
        sens["class_I_II_count"] = sens["class_I_count"] + sens["class_II_count"]
        rows.append(sens)
    data = pd.concat(rows, ignore_index=True)
    save_source(data, "fig7_threshold_robustness.csv")

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.25), constrained_layout=True)
    for ax, cancer in zip(axes, ["PDAC", "LUAD"]):
        sub = data[data["cancer"].eq(cancer)]
        pivot = sub.pivot(index="safety_factor_min", columns="high_score_percentile", values="class_I_II_count").sort_index(ascending=False)
        im = ax.imshow(pivot.values, cmap="YlGnBu", aspect="auto")
        ax.set_xticks(np.arange(pivot.shape[1]), [f"{100 * c:.0f}%" for c in pivot.columns])
        ax.set_yticks(np.arange(pivot.shape[0]), [f"{v:.1f}" for v in pivot.index])
        ax.set_xlabel("High-score gate")
        ax.set_ylabel("Safety-acceptability min")
        ax.set_title(f"{cancer}: Class I/II candidates", loc="left")
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                value = int(pivot.iloc[i, j])
                ax.text(j, i, f"{value:,}", ha="center", va="center", fontsize=8, color="white" if value > pivot.values.max() * 0.55 else "#222222")
        panel(ax, "A" if cancer == "PDAC" else "B")
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.82, pad=0.02)
    cbar.ax.set_ylabel("Class I/II count", rotation=270, labelpad=12)
    save_fig(fig, "fig7_threshold_robustness")


def main() -> None:
    """Generate all manuscript figures and their source CSV files."""
    ensure_dirs()
    figure_2()
    figure_3()
    figure_4()
    figure_5()
    figure_6()
    figure_7()
    print(f"Wrote figures to {FIG_DIR}")
    print(f"Wrote source CSV files to {SOURCE_DIR}")


if __name__ == "__main__":
    main()

