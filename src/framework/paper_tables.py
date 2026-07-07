from pathlib import Path
import json

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT_ROOT / "result" / "paper_tables"
LOG_DIR = PROJECT_ROOT / "result" / "logs"

MODULES = {
    "PDAC": PROJECT_ROOT / "PDAC",
    "LUAD": PROJECT_ROOT / "LUDA",
}

CLASS_ROWS = [
    ("Class I", "Translational-precedence candidate"),
    ("Class II", "Functional-discovery candidate"),
    ("Class III", "Dependency-supported candidate with limited actionability"),
    ("Class IV", "Safety-risk candidate"),
    ("Class V", "Biomarker-like or low-dependency candidate"),
    ("Class VI", "Insufficient or mixed evidence"),
]


def read_table(cancer: str, name: str) -> pd.DataFrame:
    """Return TSV ``name`` from the selected cancer module as a DataFrame."""
    return pd.read_csv(MODULES[cancer] / "result" / "tables" / name, sep="\t")


def read_config(cancer: str) -> dict:
    """Return the parsed JSON configuration for ``cancer``."""
    return json.loads((MODULES[cancer] / "config" / "framework_config.json").read_text(encoding="utf-8"))


def summary(cancer: str) -> pd.Series:
    """Return the metric-to-value result summary for ``cancer``."""
    return read_table(cancer, "framework_results_summary.tsv").set_index("metric")["value"]


def thresholds(cancer: str) -> dict:
    """Return the named review thresholds for ``cancer`` as a dictionary."""
    return read_table(cancer, "framework_review_thresholds.tsv").set_index("threshold")["value"].to_dict()


def class_short(value: object) -> str:
    """Return the class identifier before the first colon in ``value``."""
    return str(value).split(":", 1)[0]


def fmt_int(value: float) -> str:
    """Return numeric ``value`` as a comma-grouped integer string."""
    return f"{int(round(float(value))):,}"


def fmt_float(value: float, digits: int = 3) -> str:
    """Return ``value`` rounded to ``digits`` decimal places."""
    return f"{float(value):.{digits}f}"


def fmt_pct(value: float, digits: int = 1) -> str:
    """Return fraction ``value`` as a percentage string.

    Parameters
    ----------
    value:
        Fraction on a 0--1 scale.
    digits:
        Number of decimal places to report.

    Returns
    -------
    str
        Percentage-formatted value.
    """
    return f"{100 * float(value):.{digits}f}%"


def fmt_pct_nonzero(value: float, total: float, digits: int = 1) -> str:
    """Return a table percentage without rounding non-zero counts to 0.0%.

    Parameters
    ----------
    value:
        Numerator count.
    total:
        Denominator count.
    digits:
        Number of percentage decimal places for ordinary values.

    Returns
    -------
    str
        Percentage string; positive values below the displayed precision are
        reported as ``<0.1%`` when ``digits`` is 1.
    """
    pct = 100 * float(value) / float(total)
    threshold = 10 ** (-digits)
    if float(value) > 0 and pct < threshold:
        return f"<{threshold:.{digits}f}%"
    return f"{pct:.{digits}f}%"


def latex_escape(text: object) -> str:
    """Return ``text`` with LaTeX special characters escaped."""
    s = "" if pd.isna(text) else str(text)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        s = s.replace(old, new)
    # Render comparison operators in math mode to avoid template/font encoding
    # artifacts such as inverted punctuation in the generated PDF.
    s = s.replace(">=", r"\(\geq\)")
    s = s.replace("<=", r"\(\leq\)")
    s = s.replace("<", r"\(<\)")
    s = s.replace(">", r"\(>\)")
    return s


def average_precision_from_scores(scores: pd.Series, positives: set[str]) -> float:
    """Return benchmark average precision from full-universe model scores.

    Parameters
    ----------
    scores:
        Gene-indexed score series where larger values indicate higher priority.
    positives:
        Benchmark genes treated as validation positives.

    Returns
    -------
    float
        Average precision, equivalent to the area under the precision-recall
        step curve for a ranked candidate list.
    """
    ordered = scores.sort_values(ascending=False)
    hits = ordered.index.astype(str).isin(positives)
    n_pos = int(hits.sum())
    if n_pos == 0:
        return 0.0
    precision_at_hit = hits.cumsum()[hits] / (pd.Series(range(1, len(hits) + 1), index=ordered.index)[hits])
    return float(precision_at_hit.sum() / n_pos)


def write(path: Path, text: str) -> None:
    """Write UTF-8 ``text`` to ``path`` and return no value."""
    path.write_text(text, encoding="utf-8", newline="\n")


def table3() -> str:
    """Return the LaTeX source for the fixed-settings manuscript table."""
    cfg = read_config("PDAC")
    pdac = summary("PDAC")
    luad = summary("LUAD")
    th_p = thresholds("PDAC")
    th_l = thresholds("LUAD")
    rows = [
        ("Input scale", "Candidate universe", f"{fmt_int(pdac['candidate_universe'])} genes", f"{fmt_int(luad['candidate_universe'])} genes", "Cancer-specific evidence matrix"),
        ("Development labels", "Translational-precedence positives", fmt_int(pdac["translational_positives"]), fmt_int(luad["translational_positives"]), "Benchmark genes are embargoed"),
        ("Development labels", "Functional-discovery positives", fmt_int(pdac["functional_positives"]), fmt_int(luad["functional_positives"]), "Benchmark genes are embargoed"),
        ("Held-out benchmark", "Translational targets", fmt_int(pdac["translational_benchmark_targets"]), fmt_int(luad["translational_benchmark_targets"]), "Never used as development labels"),
        ("Held-out benchmark", "Functional targets", fmt_int(pdac["functional_benchmark_targets"]), fmt_int(luad["functional_benchmark_targets"]), "Never used as development labels"),
        ("Sampling", "Sampled-background repeats", str(cfg["background_sampling"]["n_repeats"]), "Shared", "Random seeds 1--30"),
        ("Sampling", "Background ratio", f"{cfg['background_sampling']['background_ratio']}:1 unlabeled-to-positive", "Shared", "Matched background sampling enabled"),
        ("Model", "Elastic-net logistic model", f"C={cfg['model']['C']}; l1_ratio={cfg['model']['l1_ratio']}; solver={cfg['model']['solver']}", "Shared", f"tol=0.01; max_iter={cfg['model']['max_iter']}"),
        ("Validation", "Top-K recovery cutoffs", ", ".join(map(str, cfg["validation"]["top_k"])), "Shared", "Used for recovery and comparator summaries"),
        ("Validation", "Label permutations", str(cfg["validation"]["n_permutations"]), "Shared", f"Each permutation repeats the full {cfg['validation']['permutation_repeats']}-background evaluation"),
        ("Validation", "OOF design", "LOPO or source-stratified 5-fold", "Shared", "Selected before model fitting from positive-set size and source structure"),
        ("Validation", "Leave-source-out repeats", str(cfg["validation"]["leave_source_out_repeats"]), "Shared where estimable", "Held-out-source features are removed"),
        ("Review gate", "High-score percentile", f"{100 * cfg['review']['high_score_percentile']:.0f}%", "Shared", f"Absolute cutoffs: PDAC trans {th_p['translational_high_score']:.4f}/func {th_p['functional_high_score']:.4f}; LUAD trans {th_l['translational_high_score']:.4f}/func {th_l['functional_high_score']:.4f}"),
        ("Review gate", "Dependency and disease-reference thresholds", "Top 10% each", "Shared", f"Absolute cutoffs: PDAC dep {th_p['dependency_reference']:.4f}/disease {th_p['disease_reference']:.4f}; LUAD dep {th_l['dependency_reference']:.4f}/disease {th_l['disease_reference']:.4f}"),
        ("Review gate", "Safety and actionability gates", f"safety-acceptability factor >= {cfg['review']['safety_factor_min']:.2f}; high safety risk >= {cfg['review']['safety_risk_high_threshold']:.2f}; limited actionability if druggability < {cfg['review']['limited_actionability_druggability_max']:.2f}", "Shared", "Applied in the prespecified class-priority order"),
        ("Sensitivity", "Perturbation grid", "l1_ratio 0.1/0.5/0.9; C 0.3/1.0/3.0; background ratio 4/8/12; matched background true/false", "Shared", "10 repeats per model-sensitivity variant"),
    ]
    body = "\n".join("			" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{table*}}[!tp]
	\centering
	\scriptsize
	\caption{{Fixed modelling, validation and post-model review settings used in the paired-cancer empirical analyses.}}
	\label{{tab:fixed-settings}}
	\setlength{{\tabcolsep}}{{3pt}}
	\begin{{tabularx}}{{\textwidth}}{{p{{0.14\textwidth}}p{{0.18\textwidth}}p{{0.22\textwidth}}p{{0.14\textwidth}}Y}}
		\toprule
		\textbf{{Setting group}} & \textbf{{Setting}} & \textbf{{PDAC}} & \textbf{{LUAD}} & \textbf{{Reporting note}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{table*}}
"""


def table4() -> str:
    """Return the LaTeX source for the cancer-specific input table."""
    pdac = summary("PDAC")
    luad = summary("LUAD")
    pdac_benchmark_text = (
        f"{fmt_int(pdac['translational_benchmark_targets'])} translational and "
        f"{fmt_int(pdac['functional_benchmark_targets'])} functional targets are excluded from all same-axis development labels"
    )
    luad_benchmark_text = (
        f"{fmt_int(luad['translational_benchmark_targets'])} translational and "
        f"{fmt_int(luad['functional_benchmark_targets'])} functional targets are excluded from all same-axis development labels"
    )
    rows = [
        ("Patient molecular evidence", "F: disease-relevance features", "TCGA-PAAD/GDC alteration, expression and survival records", "cBioPortal/TCGA LUAD PanCancer Atlas records; excluded from the leakage-controlled functional model", "Target-level bounded features and missingness indicators"),
        ("Protein-level evidence", "Protein and phosphoprotein support", "CPTAC/PDC-derived PDAC proteomic or proteogenomic support mapped to candidate genes", "CPTAC-LUAD/LinkedOmics proteome, phosphoproteome, acetylproteome and RNA-seq resources mapped to candidate genes", "Protein-level and phosphoproteomic features"),
        ("Functional genomics", "F/L/V: record-separated", "DepMap is a feature source; selected dependency records enter development labels and literature records form an embargoed benchmark", "TCGA and DepMap records are deterministically split within source into development and benchmark records without using model scores; source-linked features are excluded", "Development labels, benchmark targets and dependency annotations"),
        ("Drug and target knowledge", "F/L/V: record-separated", "Direct-drug and biomarker genes are deterministically split into disjoint development and benchmark records; direct drug features are excluded from the translational model", "Six LUAD anchors are development labels and seven independent translational benchmark targets are held out; direct drug features are excluded", "Embargoed translational benchmark and source-independent model features"),
        ("Tractability and safety proxies", "Post-model review gates", "UniProt, GTEx, Human Protein Atlas, gnomAD and shared tractability, druggability and safety-proxy annotations", "Same annotation logic applied after LUAD gene harmonization and feature construction", "Review-class assignment and evidence-card limitations"),
        ("Benchmark and validation sources", "V only", pdac_benchmark_text, luad_benchmark_text, "Independent recovery, OOF, full-source-out and permutation summaries"),
    ]
    body = "\n".join("			" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{table*}}[!tp]
	\centering
	\scriptsize
	\caption{{Cancer-specific input sources and positive-label definitions used in the paired analyses.}}
	\label{{tab:cancer-specific-inputs}}
	\setlength{{\tabcolsep}}{{3pt}}
	\begin{{tabularx}}{{\textwidth}}{{p{{0.15\textwidth}}p{{0.18\textwidth}}Y Y p{{0.18\textwidth}}}}
		\toprule
		\textbf{{Evidence group}} & \textbf{{Framework role}} & \textbf{{PDAC implementation}} & \textbf{{LUAD implementation}} & \textbf{{Output representation}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{table*}}
"""


def table9_source_role_controls() -> str:
    """Return the main-text source-role and leakage-control summary table.

    Returns
    -------
    str
        LaTeX source for the compact F/L/V/A source-role table.
    """
    rows = [
        ("Drug-knowledge databases", "F/A; selected records may define L/V", "Direct drug features excluded from same-axis translational model evaluation", "Limits direct label-feature reuse; cross-database fact overlap remains a limitation"),
        ("PDAC biomarker/direct-drug records", "L/V", "Deterministic development/benchmark split; same-axis benchmark genes embargoed from development labels", "Validation-only benchmark, not proof of biological independence"),
        ("LUAD therapeutic-anchor records", "L/V", "Six anchors used for development and seven held out for translational benchmark evaluation", "Tests held-out anchor recovery under small benchmark size"),
        ("Functional-dependency sources", "F/L/V", "Source-linked dependency features excluded where corresponding records define labels or benchmarks", "Functional recovery interpreted with source-family sensitivity"),
        ("Patient/protein omics sources", "F/V or A", "Cancer-context omics transformed into bounded features; LUAD functional source-linked features excluded", "Supports disease relevance, not direct therapeutic actionability"),
        ("Curated benchmark literature", "V", "Benchmark targets stored as validation records and excluded from same-axis development labels", "Protects against direct benchmark reuse while acknowledging literature-database overlap"),
    ]
    body = "\n".join("		" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{table*}}[!tp]
	\centering
	\scriptsize
	\caption{{Main-text summary of source roles and leakage controls. \(F\), feature construction; \(L\), development positive labels; \(V\), independent validation or benchmark records; \(A\), post-model annotation or review. A fuller split-rationale table is provided in Supplementary Table~\ref{{tab:leakage-control-rationale}}.}}
	\label{{tab:source-role-controls}}
	\setlength{{\tabcolsep}}{{2.5pt}}
	\begin{{tabularx}}{{\textwidth}}{{p{{0.18\textwidth}}p{{0.12\textwidth}}Y Y}}
		\toprule
		\textbf{{Source family}} & \textbf{{Role}} & \textbf{{Main embargo or exclusion rule}} & \textbf{{Interpretation}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{table*}}
"""


def table5() -> str:
    """Return the LaTeX source for the review-class size table."""
    pdac = summary("PDAC")
    luad = summary("LUAD")
    total_p = float(pdac["candidate_universe"])
    total_l = float(luad["candidate_universe"])
    rows = []
    for cls, interp in CLASS_ROWS:
        key = "class_" + cls.split()[1]
        p = float(pdac[key])
        l = float(luad[key])
        rows.append((cls, interp, fmt_int(p), fmt_pct_nonzero(p, total_p), fmt_int(l), fmt_pct_nonzero(l, total_l)))
    body = "\n".join("			" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{table*}}[!tp]
	\centering
	\small
	\caption{{Review-class sizes in the paired-cancer analyses. Percentages use each cancer-specific candidate universe as the denominator.}}
	\label{{tab:review-class-sizes}}
	\setlength{{\tabcolsep}}{{4pt}}
	\begin{{tabularx}}{{\textwidth}}{{p{{0.11\textwidth}}Yp{{0.10\textwidth}}p{{0.09\textwidth}}p{{0.10\textwidth}}p{{0.09\textwidth}}}}
		\toprule
		\textbf{{Class}} & \textbf{{Interpretation}} & \textbf{{PDAC n}} & \textbf{{PDAC \%}} & \textbf{{LUAD n}} & \textbf{{LUAD \%}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{table*}}
"""


def comparator_route_label(method: str) -> str:
    """Return a compact manuscript label for a comparator ``method``.

    Parameters
    ----------
    method:
        Exact method name stored in ``framework_comparator_summary.tsv``.

    Returns
    -------
    str
        Short route label used in Table 6.
    """
    labels = {
        "Dual-axis sampled-background PU": "Dual-axis PU-style",
        "Combined single-axis sampled-background PU": "Merged single-axis PU-style",
        "Dependency-only": "Dependency-only ranking",
        "Rank aggregation": "Rank aggregation",
    }
    return labels.get(method, method)


def comparator_metric(cancer: str, method: str, axis: str) -> dict:
    """Return formatted comparator metrics for one cancer, method and axis.

    Parameters
    ----------
    cancer:
        Cancer module key, either ``PDAC`` or ``LUAD``.
    method:
        Exact comparator route name in the source TSV table.
    axis:
        Benchmark axis to report, either ``translational`` or ``functional``.

    Returns
    -------
    dict
        Route label, benchmark size, top-50 hits/recovery and review fractions
        used to assemble the compact main-text comparator table.
    """
    df = read_table(cancer, "framework_comparator_summary.tsv")
    row = df[(df["method"].eq(method)) & (df["axis"].eq(axis)) & (df["top_k"].eq(50))].iloc[0]
    recovery = float(row["positive_recovery"])
    benchmark_n = int(float(summary(cancer)[f"{axis}_benchmark_targets"]))
    hits = int(round(recovery * benchmark_n))
    return {
        "label": comparator_route_label(method),
        "benchmark_n": benchmark_n,
        "recovery": recovery,
        "hits": hits,
        "hit_text": f"{hits}/{benchmark_n} ({recovery:.3f})",
        "median_rank_percentile": float(row["median_rank_percentile"]),
        "class_I_II_fraction": float(row["class_I_II_fraction"]),
        "high_safety_risk_fraction": float(row["high_safety_risk_fraction"]),
    }


def comparator_summary_row(cancer: str, axis: str, comparator_method: str, interpretation: str) -> tuple:
    """Return one compact row comparing the primary route with a key baseline."""
    primary = comparator_metric(cancer, "Dual-axis sampled-background PU", axis)
    comparator = comparator_metric(cancer, comparator_method, axis)
    signal = "trans." if axis == "translational" else "func."
    return (
        f"{cancer} {signal}",
        primary["hit_text"],
        f"{comparator['label']} {comparator['hit_text']}",
        f"{interpretation}; dual-axis review profile: I/II {primary['class_I_II_fraction']:.2f}, high-risk {primary['high_safety_risk_fraction']:.2f}",
    )


def table6() -> str:
    """Return the compact main-text comparator table.

    Returns
    -------
    str
        LaTeX table summarizing the primary dual-axis route against the most
        important comparator for each cancer--axis setting. Full comparator
        outputs remain available in the machine-readable result tables.
    """
    rows = [
        comparator_summary_row("PDAC", "translational", "Rank aggregation", "Rank aggregation recovered one additional benchmark target but lacks axis-specific model and evidence-card linkage"),
        comparator_summary_row("PDAC", "functional", "Combined single-axis sampled-background PU", "Merged single-axis matched top-50 recovery but collapses translational and functional objectives"),
        comparator_summary_row("LUAD", "translational", "Rank aggregation", "Rank aggregation recovered more translational targets; dual-axis route preserves interpretable review linkage"),
        comparator_summary_row("LUAD", "functional", "Dependency-only", "Dependency-only recovery was higher but bypasses tractability, safety and multi-source review gates"),
    ]
    body = "\n".join("			" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{table*}}[!tp]
	\centering
	\scriptsize
	\caption{{Main-text comparator summary for held-out benchmark recovery. All rows use the prespecified top-50 cutoff; benchmark targets were excluded from same-axis development labels. Full route-level comparator records are provided in the machine-readable result tables.}}
	\label{{tab:comparator-results}}
	\setlength{{\tabcolsep}}{{3pt}}
	\begin{{tabularx}}{{\textwidth}}{{>{{\raggedright\arraybackslash}}p{{0.18\textwidth}}>{{\raggedright\arraybackslash}}p{{0.15\textwidth}}>{{\raggedright\arraybackslash}}p{{0.24\textwidth}}Y}}
		\toprule
		\textbf{{Setting}} & \textbf{{Dual-axis top-50}} & \textbf{{Key comparator top-50}} & \textbf{{Interpretation}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{table*}}
"""


def supplementary_leakage_control_table() -> str:
    """Return the LaTeX source for the leakage-control split-rationale table.

    Returns
    -------
    str
        Supplementary table describing deterministic split rationale, label
        embargoes, and source-linked feature exclusions.
    """
    rows = [
        (
            "Deterministic record split",
            "Alphabetical alternation used only to assign records reproducibly when a source contributed both development and benchmark records",
            "The split used gene/source identifiers and did not use model scores, fitted ranks, benchmark recovery, or post-model classes",
            "Prevents outcome-dependent benchmark construction; it is a deterministic separation rule, not evidence of biological independence by itself",
        ),
        (
            "Same-axis label embargo",
            "Held-out benchmark genes were excluded from same-axis development-positive labels and sampled unlabeled background",
            "Applied separately for translational-precedence and functional-discovery axes",
            "Prevents direct reuse of validation targets as training positives or temporary negatives on the same axis",
        ),
        (
            "Direct drug-feature exclusion",
            "Direct drug-target component columns were excluded from translational-precedence modelling where they would duplicate translational benchmark evidence",
            "The lower-resolution druggability composite was retained as a tractability/actionability proxy and post-model review variable",
            "Reduces direct component lookup while acknowledging residual overlap from composite druggability and cross-database target-level facts",
        ),
        (
            "LUAD functional source-linked exclusion",
            "TCGA/DepMap source-linked features were removed from the leakage-controlled LUAD functional model when corresponding source records were used for labels or benchmarks",
            "Applied before model fitting in the LUAD functional setting",
            "Separates functional labels from direct source-linked feature evidence in the held-out validation setting",
        ),
        (
            "Overlapping drug-knowledge resources",
            "Open Targets, DGIdb, ChEMBL, TTD and DrugCentral may contain overlapping therapeutic or target records",
            "Records were harmonized at target level and assigned to feature, label, validation, or annotation roles before modelling",
            "Avoids treating duplicated drug-knowledge records as independent validation evidence for the same analytical axis",
        ),
        (
            "Leave-source-out validation",
            "Held-out-source labels were removed and source-linked feature columns were excluded where estimable",
            "Reported as full-source-out transfer summaries in Figure 5C",
            "Tests whether ranking signal persists beyond a single curated source or source-linked feature block",
        ),
    ]
    body = "\n".join("			" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{sidewaystable*}}[!tp]
	\centering
	\scriptsize
	\caption{{Supplementary leakage-control and split-rationale summary. The table distinguishes deterministic record separation from biological or database independence and lists the main controls used to reduce label--feature reuse.}}
	\label{{tab:leakage-control-rationale}}
	\setlength{{\tabcolsep}}{{3pt}}
	\begin{{tabularx}}{{\textwidth}}{{p{{0.18\textwidth}}Y Y Y}}
		\toprule
		\textbf{{Control item}} & \textbf{{What was done}} & \textbf{{Implementation detail}} & \textbf{{Rationale and limitation addressed}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{sidewaystable*}}
"""


def supplementary_benchmark_uncertainty_table() -> str:
    """Return a compact validation-uncertainty supplementary table.

    Returns
    -------
    str
        LaTeX table reporting raw top-50 hits, bootstrap confidence intervals,
        precision@K, average precision and early-enrichment summaries.
    """
    rows = []
    for cancer in ["PDAC", "LUAD"]:
        scores = read_table(cancer, "framework_priority_scores.tsv").set_index("gene_symbol")
        validation = read_table(cancer, "framework_independent_benchmark_validation.tsv")
        benchmarks = read_table(cancer, "framework_validation_benchmarks.tsv")
        for axis, score_col in [
            ("translational", "translational_precedence_score"),
            ("functional", "functional_discovery_score"),
        ]:
            row = validation[validation["axis"].eq(axis)].iloc[0]
            genes = set(benchmarks[benchmarks["axis"].eq(axis)]["gene_symbol"].astype(str))
            hits50 = int(round(float(row["top50_recovery"]) * len(genes)))
            ap = average_precision_from_scores(scores[score_col], genes)
            rows.append(
                (
                    cancer,
                    axis.capitalize(),
                    f"{hits50}/{len(genes)}",
                    fmt_float(row["top50_recovery"]),
                    f"{fmt_float(row['top50_recovery_ci_low'])}--{fmt_float(row['top50_recovery_ci_high'])}",
                    fmt_float(row["precision_at_50"]),
                    fmt_float(ap),
                    fmt_float(row["enrichment_top50_over_random"], 1),
                    fmt_float(row["median_rank_percentile"]),
                )
            )
    body = "\n".join("			" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{table*}}[!tp]
	\centering
	\scriptsize
	\caption{{Supplementary small-sample validation uncertainty summary. Top-50 recovery is reported with raw hit counts and target-level bootstrap confidence intervals because each benchmark contains a limited number of genes. Average precision was calculated from the full ranked candidate universe against the same held-out benchmark labels.}}
	\label{{tab:benchmark-uncertainty}}
	\setlength{{\tabcolsep}}{{3pt}}
	\begin{{tabularx}}{{\textwidth}}{{p{{0.07\textwidth}}p{{0.09\textwidth}}p{{0.08\textwidth}}p{{0.09\textwidth}}p{{0.12\textwidth}}p{{0.08\textwidth}}p{{0.08\textwidth}}p{{0.09\textwidth}}Y}}
		\toprule
		\textbf{{Cancer}} & \textbf{{Axis}} & \textbf{{Hits}} & \textbf{{Recovery}} & \textbf{{Bootstrap CI}} & \textbf{{P@50}} & \textbf{{AP}} & \textbf{{Enrich.}} & \textbf{{Rank pct.}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{table*}}
"""


def supplementary_benchmark_rank_table() -> str:
    """Return a longtable listing every held-out benchmark gene rank.

    Returns
    -------
    str
        Paginated LaTeX longtable with one row per cancer, axis and benchmark
        gene. Full source data are also written as TSV by ``write_tables``.
    """
    rows = []
    source_rows = []
    for cancer in ["PDAC", "LUAD"]:
        scores = read_table(cancer, "framework_priority_scores.tsv").set_index("gene_symbol")
        benchmarks = read_table(cancer, "framework_validation_benchmarks.tsv")
        n = len(scores)
        for axis, rank_col, score_col in [
            ("translational", "translational_precedence_rank", "translational_precedence_score"),
            ("functional", "functional_discovery_rank", "functional_discovery_score"),
        ]:
            sub = benchmarks[benchmarks["axis"].eq(axis)].copy()
            grouped = (
                sub.groupby("gene_symbol")
                .agg(
                    source=("source", lambda x: "; ".join(sorted(set(map(str, x))))),
                    source_group=("source_group", lambda x: "; ".join(sorted(set(map(str, x))))),
                    notes=("notes", lambda x: "; ".join(sorted(set(str(v) for v in x if str(v) and str(v) != "nan")))),
                )
                .reset_index()
            )
            for _, item in grouped.sort_values("gene_symbol").iterrows():
                gene = str(item["gene_symbol"])
                if gene not in scores.index:
                    rank = ""
                    rank_pct = ""
                    score = ""
                    top50 = "No"
                else:
                    rank_value = int(scores.loc[gene, rank_col])
                    rank = fmt_int(rank_value)
                    rank_pct_value = 1 - (rank_value - 1) / max(n - 1, 1)
                    rank_pct = fmt_float(rank_pct_value)
                    score = fmt_float(scores.loc[gene, score_col])
                    top50 = "Yes" if rank_value <= 50 else "No"
                record = {
                    "cancer": cancer,
                    "axis": axis,
                    "gene_symbol": gene,
                    "rank": rank,
                    "rank_percentile": rank_pct,
                    "score": score,
                    "top50_hit": top50,
                    "source": item["source"],
                    "source_group": item["source_group"],
                    "notes": item["notes"],
                }
                source_rows.append(record)
                rows.append(
                    (
                        cancer,
                        "Trans." if axis == "translational" else "Func.",
                        gene,
                        rank,
                        rank_pct,
                        score,
                        top50,
                        item["source_group"],
                    )
                )
    source = pd.DataFrame(source_rows)
    source.to_csv(TABLE_DIR / "supplementary_benchmark_gene_ranks.tsv", sep="\t", index=False)
    body = "\n".join("		" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begingroup
\scriptsize
\setlength{{\tabcolsep}}{{2pt}}
\begin{{longtable}}{{@{{}}p{{0.07\textwidth}}p{{0.08\textwidth}}p{{0.12\textwidth}}p{{0.08\textwidth}}p{{0.09\textwidth}}p{{0.08\textwidth}}p{{0.06\textwidth}}p{{0.20\textwidth}}@{{}}}}
	\caption{{Supplementary held-out benchmark gene ranks. Rank percentile is scaled so that values closer to 1 indicate higher-ranked genes. A machine-readable source version with source names and notes is provided as \texttt{{supplementary\_benchmark\_gene\_ranks.tsv}}.}}
	\label{{tab:benchmark-gene-ranks}}\\
	\toprule
	\textbf{{Cancer}} & \textbf{{Axis}} & \textbf{{Gene}} & \textbf{{Rank}} & \textbf{{Rank pct.}} & \textbf{{Score}} & \textbf{{Top 50}} & \textbf{{Source group}} \\
	\midrule
	\endfirsthead
	\caption[]{{Supplementary held-out benchmark gene ranks, continued.}}\\
	\toprule
	\textbf{{Cancer}} & \textbf{{Axis}} & \textbf{{Gene}} & \textbf{{Rank}} & \textbf{{Rank pct.}} & \textbf{{Score}} & \textbf{{Top 50}} & \textbf{{Source group}} \\
	\midrule
	\endhead
{body}
	\bottomrule
\end{{longtable}}
\endgroup
"""


def supplementary_threshold_robustness_table() -> str:
    """Return threshold-robustness table for review-class assignment.

    Returns
    -------
    str
        LaTeX table summarizing class I/II candidate counts across high-score
        and safety-acceptability thresholds.
    """
    rows = []
    for cancer in ["PDAC", "LUAD"]:
        sens = read_table(cancer, "framework_review_threshold_sensitivity.tsv")
        for _, row in sens.iterrows():
            class_i = int(row["class_I_count"])
            class_ii = int(row["class_II_count"])
            rows.append(
                (
                    cancer,
                    fmt_pct(row["high_score_percentile"], 0),
                    fmt_float(row["safety_factor_min"], 1),
                    fmt_int(class_i),
                    fmt_int(class_ii),
                    fmt_int(class_i + class_ii),
                    fmt_int(row["class_IV_count"]),
                    fmt_int(row["class_VI_count"]),
                )
            )
    body = "\n".join("			" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{table*}}[!tp]
	\centering
	\scriptsize
	\caption{{Supplementary review-threshold robustness summary. The primary analysis used a top-2\% high-score gate and safety-acceptability factor threshold of 0.40; neighbouring thresholds show how review-class sizes change under prespecified perturbations.}}
	\label{{tab:threshold-robustness}}
	\setlength{{\tabcolsep}}{{4pt}}
	\begin{{tabularx}}{{\textwidth}}{{p{{0.09\textwidth}}p{{0.13\textwidth}}p{{0.14\textwidth}}p{{0.10\textwidth}}p{{0.10\textwidth}}p{{0.12\textwidth}}p{{0.10\textwidth}}Y}}
		\toprule
		\textbf{{Cancer}} & \textbf{{High-score gate}} & \textbf{{Safety-accept. min}} & \textbf{{Class I}} & \textbf{{Class II}} & \textbf{{Class I+II}} & \textbf{{Class IV}} & \textbf{{Class VI}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{table*}}
"""


def supplementary_feature_mapping_table() -> str:
    """Return a core feature mapping table for method reproducibility.

    Returns
    -------
    str
        LaTeX table describing evidence blocks, transformations, missing-value
        handling and axis/review use.
    """
    rows = []
    manifest = read_table("PDAC", "framework_feature_manifest.tsv")
    keep = manifest[~manifest["feature"].str.endswith("_missing", na=False)].copy()
    for _, row in keep.iterrows():
        axis_use = []
        if int(row["entered_translational_model"]):
            axis_use.append("T")
        if int(row["entered_functional_model"]):
            axis_use.append("F")
        if int(row["used_in_review_gate"]):
            axis_use.append("Review")
        rows.append(
            (
                row["source_family"],
                row["feature"].replace("safety_factor", "safety_acceptability_factor"),
                row["transformation"],
                row["missing_value_rule"],
                "/".join(axis_use) if axis_use else "source only",
            )
        )
    body = "\n".join("			" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{sidewaystable*}}[!tp]
	\centering
	\scriptsize
	\caption{{Supplementary core feature mapping for reproducibility. Feature transformations map heterogeneous source values to bounded 0--1 evidence variables before axis-specific exclusions. Use denotes T, translational-precedence model; F, functional-discovery model; and Review, post-model review gates. Direct drug-knowledge component columns were excluded from leakage-controlled translational-precedence modelling, whereas the lower-resolution druggability composite was retained as a tractability/actionability proxy and review variable. Safety-acceptability factor equals \(1-\) safety-risk score, so high values indicate more acceptable safety-proxy evidence whereas high safety-risk score indicates greater risk.}}
	\label{{tab:feature-mapping}}
	\setlength{{\tabcolsep}}{{3pt}}
	\begin{{tabularx}}{{\textwidth}}{{p{{0.16\textwidth}}p{{0.18\textwidth}}Y Y p{{0.13\textwidth}}}}
		\toprule
		\textbf{{Evidence/source family}} & \textbf{{Feature}} & \textbf{{0--1 transformation}} & \textbf{{Missing-value rule}} & \textbf{{Use}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{sidewaystable*}}
"""


def supplementary_source_family_stress_table() -> str:
    """Return a source-family leakage stress-test summary table.

    Returns
    -------
    str
        LaTeX table distinguishing completed family-level exclusions from
        already-embargoed or underpowered stress tests.
    """
    rows = [
        (
            "Drug-knowledge family",
            "Translational",
            "chemical_modulation_score, clinical_precedence_score and Open Targets association features; composite druggability_score retained",
            "Direct component columns were excluded from the primary translational models in both cancers",
            "Independent translational recovery did not use the direct component columns, but residual risk remains from composite druggability, target-level fact overlap and external benchmark literature/database overlap",
        ),
        (
            "DepMap/functional-dependency family",
            "LUAD functional",
            "dependency_score and external reproducibility linked to DepMap development or benchmark records",
            "Already excluded from the leakage-controlled LUAD functional model",
            "LUAD functional recovery should be interpreted conservatively under source-family stress testing despite nominal significance in the primary label-permutation benchmark",
        ),
        (
            "DepMap/functional-dependency family",
            "PDAC functional",
            "known_dependency source-out excluded dependency_score in the leave-source-out analysis",
            "Top-50 recovery was 0.000 for the held-out known_dependency source-out test",
            "This indicates strong dependence on source definition for some functional signals and supports conservative interpretation of functional-discovery results",
        ),
        (
            "Patient-omics family",
            "LUAD functional",
            "TCGA source-out excluded patient_alteration_score, TCGA expression, mutation, CNV and survival features",
            "Top-50 recovery was 0.000 for the held-out tcga_luad_nature_2014 source-out test",
            "The source-family stress test did not support claiming source-independent LUAD functional recovery",
        ),
        (
            "Literature benchmark family",
            "Both axes",
            "Expanded literature records used only as validation records",
            "Same-axis label embargo and source-data manifest record the literature extraction source",
            "Protects against direct label reuse, but does not prove independence from all database-curated knowledge; this is stated as a limitation",
        ),
    ]
    body = "\n".join("			" + " & ".join(map(latex_escape, r)) + r" \\" for r in rows)
    return rf"""\begin{{sidewaystable*}}[!tp]
	\centering
	\scriptsize
	\caption{{Supplementary source-family leakage stress-test summary. The table separates direct feature exclusions already present in the primary models from source-family stress tests that remain underpowered or nonsignificant.}}
	\label{{tab:source-family-stress}}
	\setlength{{\tabcolsep}}{{3pt}}
	\begin{{tabularx}}{{\textwidth}}{{p{{0.16\textwidth}}p{{0.12\textwidth}}Y Y Y}}
		\toprule
		\textbf{{Source family}} & \textbf{{Setting}} & \textbf{{Feature/source exclusion}} & \textbf{{Observed result or design status}} & \textbf{{Interpretation}} \\
		\midrule
{body}
		\bottomrule
	\end{{tabularx}}
\end{{sidewaystable*}}
"""


def write_tables() -> None:
    """Generate all manuscript table fragments and return no value."""
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    write(TABLE_DIR / "table3_fixed_settings.tex", table3())
    write(TABLE_DIR / "table4_cancer_specific_inputs.tex", table4())
    write(TABLE_DIR / "table9_source_role_controls.tex", table9_source_role_controls())
    write(TABLE_DIR / "table5_review_class_sizes.tex", table5())
    write(TABLE_DIR / "table6_comparator_results.tex", table6())
    write(TABLE_DIR / "supplementary_leakage_control_split_rationale.tex", supplementary_leakage_control_table())
    write(TABLE_DIR / "supplementary_benchmark_uncertainty.tex", supplementary_benchmark_uncertainty_table())
    write(TABLE_DIR / "supplementary_benchmark_gene_ranks.tex", supplementary_benchmark_rank_table())
    write(TABLE_DIR / "supplementary_threshold_robustness.tex", supplementary_threshold_robustness_table())
    write(TABLE_DIR / "supplementary_feature_mapping.tex", supplementary_feature_mapping_table())
    write(TABLE_DIR / "supplementary_source_family_stress_tests.tex", supplementary_source_family_stress_table())


def consistency_log() -> None:
    """Validate cross-table counts and write a PASS/NEEDS_REVIEW log."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    lines = []
    ok = True
    for cancer in ["PDAC", "LUAD"]:
        s = summary(cancer)
        review = read_table(cancer, "framework_review_classes.tsv")
        scores = read_table(cancer, "framework_priority_scores.tsv")
        positives = read_table(cancer, "framework_development_positive_sets.tsv")
        leakage = read_table(cancer, "framework_leakage_audit.tsv")
        class_sum = sum(float(s[f"class_{roman}"]) for roman in ["I", "II", "III", "IV", "V", "VI"])
        universe = float(s["candidate_universe"])
        lines.append(f"{cancer}: class-count sum = {class_sum:.0f}; candidate universe = {universe:.0f}.")
        if int(class_sum) != int(universe) or len(review) != int(universe) or len(scores) != int(universe):
            ok = False
            lines.append(f"TODO_NUMERIC_CONFLICT: {cancer} candidate universe, class table and score table lengths are inconsistent.")
        trans = int(positives["clinical_positive"].sum())
        func = int(positives["discovery_positive"].sum())
        lines.append(f"{cancer}: positive-set counts from source_positive_sets.tsv = translational {trans}, functional {func}.")
        if trans != int(float(s["translational_positives"])) or func != int(float(s["functional_positives"])):
            ok = False
            lines.append(f"TODO_NUMERIC_CONFLICT: {cancer} positive counts differ between source_positive_sets.tsv and framework_results_summary.tsv.")
        if not leakage["status"].eq("PASS").all() or not leakage["gene_overlap_count"].eq(0).all():
            ok = False
            lines.append(f"TODO_NUMERIC_CONFLICT: {cancer} leakage audit did not pass.")
        classes = sorted(review["review_class"].map(class_short).unique())
        if classes != CLASS_ROWS_TO_CHECK:
            ok = False
            lines.append(f"TODO_NUMERIC_CONFLICT: {cancer} review classes are {classes}.")
    lines.append("Comparator top-K cutoff checked against config: top_k includes 50 for both cancers.")
    for cancer in ["PDAC", "LUAD"]:
        cfg = read_config(cancer)
        if 50 not in cfg["validation"]["top_k"]:
            ok = False
            lines.append(f"TODO_NUMERIC_CONFLICT: {cancer} config does not include top_k=50.")
    lines.append("Status: PASS" if ok else "Status: NEEDS_REVIEW")
    write(LOG_DIR / "results_consistency_check.txt", "\n".join(lines) + "\n")


CLASS_ROWS_TO_CHECK = [row[0] for row in CLASS_ROWS]


def main() -> None:
    """Generate table fragments and their consistency log."""
    write_tables()
    consistency_log()
    print(f"Wrote tables to {TABLE_DIR}")
    print(f"Wrote consistency log to {LOG_DIR / 'results_consistency_check.txt'}")


if __name__ == "__main__":
    main()

