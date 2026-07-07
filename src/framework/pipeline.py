from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


CLINICAL_SCORE = "translational_precedence_score"
DISCOVERY_SCORE = "functional_discovery_score"


@dataclass(frozen=True)
class AxisSpec:
    name: str
    label_column: str
    score_column: str
    rank_column: str
    top_prefix: str


AXES = [
    AxisSpec("translational", "clinical_positive", CLINICAL_SCORE, "translational_precedence_rank", "translational"),
    AxisSpec("functional", "discovery_positive", DISCOVERY_SCORE, "functional_discovery_rank", "functional"),
]


FEATURE_DETAILS = {
    "patient_alteration_score": ("patient omics", "noisy-OR of bounded cancer-context alteration components"),
    "proteomic_support_score": ("CPTAC/proteomics", "noisy-OR of bounded protein and phosphoprotein components"),
    "external_reproducibility_score": ("external cohorts", "cancer-specific bounded external-evidence composite"),
    "dependency_score": ("DepMap", "bounded lineage dependency score; larger values indicate stronger dependency"),
    "chemical_modulation_score": ("drug knowledge", "min(drug evidence source count / 5, 1)"),
    "clinical_precedence_score": ("drug/clinical knowledge", "mean of bounded approval and clinical-development components"),
    "target_tractability_score": ("target annotations", "fraction of observed tractability indicators that are positive"),
    "druggability_score": ("drug/target annotations", "mean of chemical modulation, clinical precedence and tractability"),
    "safety_risk_score": ("GTEx/HPA/gnomAD/DepMap", "mean of bounded normal-expression, constraint and common-essentiality risks"),
    "safety_factor": ("safety annotations", "1 - safety_risk_score"),
    "therapeutic_window_score": ("safety annotations", "1 - safety_risk_score"),
    "tcga_gtex_expression_score": ("TCGA/GTEx", "bounded tumour-versus-normal expression contrast"),
    "tcga_mutation_score": ("TCGA", "bounded cancer-cohort mutation frequency"),
    "tcga_cnv_score": ("TCGA", "bounded copy-number alteration burden"),
    "tcga_survival_score": ("TCGA", "bounded survival-association evidence"),
    "cptac_mrna_score": ("CPTAC", "bounded tumour-versus-normal mRNA evidence"),
    "cptac_protein_score": ("CPTAC", "bounded tumour-versus-normal protein evidence"),
    "phosphoproteomic_support_score": ("CPTAC", "bounded phosphoproteomic detection/support evidence"),
    "geo_external_reproducibility_score": ("GEO", "bounded cross-cohort direction and significance reproducibility"),
    "opentargets_overall_score": ("Open Targets", "maximum bounded cancer-target association score"),
}


def read_config(project_root: Path, config_path: Path | None = None) -> dict:
    path = config_path or (project_root / "config" / "framework_config.json")
    with path.open("r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    cfg["_config_path"] = str(path.resolve())
    if "legacy_project" in cfg:
        cfg["legacy_project"] = str((project_root / cfg["legacy_project"]).resolve())
    return cfg


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"gene_symbol": str}, low_memory=False)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def _numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce")
    return pd.Series(default, index=df.index, dtype=float)


def _safe_gene_index(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["gene_symbol"] = out["gene_symbol"].astype(str)
    out = out.drop_duplicates("gene_symbol", keep="first")
    return out.set_index("gene_symbol", drop=True)


def copy_framework_inputs(cfg: dict, out_tables: Path, project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cancer_code = str(cfg.get("cancer_code", "pdac")).lower()
    if "input_tables" in cfg:
        inputs = cfg["input_tables"]
        base = Path(inputs.get("base_dir", "."))
        if not base.is_absolute():
            base = (project_root / base).resolve()
        matrix = read_tsv(base / inputs.get("evidence_matrix", f"{cancer_code}_target_evidence_matrix.tsv"))
        label_sources = read_tsv(base / inputs.get("label_sources", "ml_label_sources.tsv"))
        positive_sets = read_tsv(base / inputs.get("positive_sets", "ml_positive_sets.tsv"))
    else:
        legacy_tables = Path(cfg["legacy_project"]) / "results" / "tables"
        matrix = read_tsv(legacy_tables / "pdac_target_evidence_matrix.tsv")
        label_sources = read_tsv(legacy_tables / "ml_label_sources.tsv")
        positive_sets = read_tsv(legacy_tables / "ml_positive_sets.tsv")
    label_sources, positive_sets = apply_positive_label_filters(label_sources, positive_sets, matrix, cfg)
    write_tsv(matrix, out_tables / f"source_{cancer_code}_target_evidence_matrix.tsv")
    write_tsv(label_sources, out_tables / "source_label_sources.tsv")
    write_tsv(positive_sets, out_tables / "source_positive_sets.tsv")
    return matrix, label_sources, positive_sets


def build_validation_contract(
    label_sources: pd.DataFrame,
    matrix: pd.DataFrame,
    cfg: dict,
) -> tuple[pd.DataFrame, dict[str, set[str]], pd.DataFrame, pd.DataFrame]:
    """Split positive-label records into development and validation roles.

    Parameters
    ----------
    label_sources:
        Positive-label records after cancer-specific filtering.
    matrix:
        Candidate evidence matrix used to restrict records to eligible genes.
    cfg:
        Framework configuration containing ``evaluation.axes`` split rules.

    Returns
    -------
    tuple
        Development positive sets, benchmark genes by axis, record-level role
        manifest, and leakage-audit rows.
    """
    evaluation = cfg.get("evaluation", {})
    axis_cfgs = evaluation.get("axes", {})
    if not axis_cfgs:
        raise ValueError("evaluation.axes must define development and benchmark records")

    eligible = set(matrix["gene_symbol"].astype(str))
    labels = label_sources.copy()
    labels["gene_symbol"] = labels["gene_symbol"].astype(str)
    labels = labels[labels["gene_symbol"].isin(eligible)].copy()
    labels["is_positive"] = pd.to_numeric(labels.get("is_positive", 1), errors="coerce").fillna(0).astype(int)
    labels = labels[labels["is_positive"].eq(1)].copy()

    manifest_parts: list[pd.DataFrame] = []
    development: dict[str, set[str]] = {}
    benchmarks: dict[str, set[str]] = {}
    audit_rows: list[dict] = []
    for axis in AXES:
        source_axis = "clinical" if axis.name == "translational" else "discovery"
        spec = axis_cfgs.get(axis.name, {})
        axis_rows = labels[labels["priority_axis"].isin([source_axis, "both"])].copy()
        dev_sources = set(spec.get("development_sources", []))
        benchmark_sources = set(spec.get("benchmark_sources", []))
        dev_gene_filter = set(map(str, spec.get("development_genes", [])))
        benchmark_gene_filter = set(map(str, spec.get("benchmark_genes", [])))

        benchmark_mask = axis_rows["source"].isin(benchmark_sources)
        if benchmark_gene_filter:
            benchmark_mask &= axis_rows["gene_symbol"].isin(benchmark_gene_filter)
        benchmark_genes = set(axis_rows.loc[benchmark_mask, "gene_symbol"])

        development_mask = axis_rows["source"].isin(dev_sources)
        if dev_gene_filter:
            development_mask &= axis_rows["gene_symbol"].isin(dev_gene_filter)
        development_genes = set(axis_rows.loc[development_mask, "gene_symbol"])
        if bool(evaluation.get("benchmark_gene_embargo", True)):
            development_genes -= benchmark_genes

        if len(development_genes) < 2:
            raise ValueError(f"{axis.name} development split has fewer than two positives")
        if len(benchmark_genes) < 2:
            raise ValueError(f"{axis.name} benchmark split has fewer than two targets")

        role_rows = axis_rows.copy()
        role_rows["axis"] = axis.name
        role_rows["record_role"] = "A"
        role_rows["split_reason"] = "interpretation-only or unused positive record"
        is_v = role_rows["source"].isin(benchmark_sources) & role_rows["gene_symbol"].isin(benchmark_genes)
        is_l = role_rows["source"].isin(dev_sources) & role_rows["gene_symbol"].isin(development_genes)
        role_rows.loc[is_v, ["record_role", "split_reason"]] = ["V", "held-out benchmark record"]
        role_rows.loc[is_l, ["record_role", "split_reason"]] = ["L", "development positive-label record"]
        embargoed = role_rows["gene_symbol"].isin(benchmark_genes) & ~is_v
        role_rows.loc[embargoed, "split_reason"] = "benchmark gene embargoed from development labels"
        role_rows["record_id"] = (
            role_rows["source"].astype(str) + "|" + axis.name + "|" + role_rows["gene_symbol"].astype(str)
        )
        manifest_parts.append(role_rows)
        development[axis.name] = development_genes
        benchmarks[axis.name] = benchmark_genes
        overlap = development_genes & benchmark_genes
        audit_rows.append({
            "axis": axis.name,
            "development_positive_count": len(development_genes),
            "benchmark_target_count": len(benchmark_genes),
            "gene_overlap_count": len(overlap),
            "record_role_conflicts": 0,
            "status": "PASS" if not overlap else "FAIL",
        })

    role_manifest = pd.concat(manifest_parts, ignore_index=True)
    conflicts = role_manifest.groupby(["record_id", "axis"])["record_role"].nunique()
    conflict_n = int((conflicts > 1).sum())
    if conflict_n:
        raise ValueError(f"Record-role manifest contains {conflict_n} conflicting assignments")
    for row in audit_rows:
        row["record_role_conflicts"] = conflict_n
        if conflict_n:
            row["status"] = "FAIL"
    if any(row["status"] != "PASS" for row in audit_rows):
        raise ValueError("Validation contract failed leakage audit")

    genes = sorted(development["translational"] | development["functional"])
    positive_sets = pd.DataFrame({"gene_symbol": genes})
    positive_sets["clinical_positive"] = positive_sets["gene_symbol"].isin(development["translational"]).astype(int)
    positive_sets["discovery_positive"] = positive_sets["gene_symbol"].isin(development["functional"]).astype(int)
    positive_sets["positive_source_count"] = positive_sets["gene_symbol"].map(
        role_manifest[role_manifest["record_role"].eq("L")].groupby("gene_symbol")["source"].nunique()
    ).fillna(0).astype(int)
    for axis_name, column in [("translational", "clinical_sources"), ("functional", "discovery_sources")]:
        source_map = (
            role_manifest[(role_manifest["record_role"].eq("L")) & (role_manifest["axis"].eq(axis_name))]
            .groupby("gene_symbol")["source"]
            .agg(lambda x: ";".join(sorted(set(map(str, x)))))
        )
        positive_sets[column] = positive_sets["gene_symbol"].map(source_map).fillna("")
    positive_sets["heldout_eligible_sources"] = positive_sets["clinical_sources"].where(
        positive_sets["discovery_sources"].eq(""),
        positive_sets["clinical_sources"] + ";" + positive_sets["discovery_sources"],
    ).str.strip(";")
    return positive_sets, benchmarks, role_manifest, pd.DataFrame(audit_rows)


def axis_feature_columns(feature_matrix: pd.DataFrame, cfg: dict, axis_name: str) -> list[str]:
    """Return the leakage-controlled feature columns for one model axis."""
    excluded = set(cfg.get("evaluation", {}).get("axes", {}).get(axis_name, {}).get("feature_exclusions", []))
    excluded |= {f"{feature}_missing" for feature in excluded}
    return [column for column in feature_matrix.columns if column != "gene_symbol" and column not in excluded]


def apply_positive_label_filters(
    label_sources: pd.DataFrame,
    positive_sets: pd.DataFrame,
    matrix: pd.DataFrame,
    cfg: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    filters = cfg.get("positive_label_filters", {})
    if not filters:
        return label_sources, positive_sets

    labels = label_sources.copy()
    labels["gene_symbol"] = labels["gene_symbol"].astype(str)
    labels = labels[labels["gene_symbol"].isin(set(matrix["gene_symbol"].astype(str)))].copy()

    clinical = filters.get("clinical", {})
    allowed_sources = set(clinical.get("include_sources", []))
    allowed_groups = set(clinical.get("include_source_groups", []))
    require_symbol_like = bool(clinical.get("require_symbol_like_gene", False))
    is_clinical = labels["priority_axis"].isin(["clinical", "both"])
    keep_clinical = pd.Series(True, index=labels.index)
    if allowed_sources:
        keep_clinical &= labels["source"].isin(allowed_sources)
    if allowed_groups:
        keep_clinical &= labels["source_group"].isin(allowed_groups)
    if require_symbol_like:
        keep_clinical &= labels["gene_symbol"].str.match(r"^[A-Za-z][A-Za-z0-9.-]*$", na=False)
    labels = labels[(~is_clinical) | keep_clinical].drop_duplicates(["gene_symbol", "source", "priority_axis"]).copy()

    if labels.empty:
        raise ValueError("Positive-label filtering removed all label sources")

    rows = []
    for gene, group in labels.groupby("gene_symbol", sort=True):
        clinical_rows = group[group["priority_axis"].isin(["clinical", "both"])]
        discovery_rows = group[group["priority_axis"].isin(["discovery", "both"])]
        rows.append({
            "gene_symbol": gene,
            "clinical_positive": int(not clinical_rows.empty),
            "discovery_positive": int(not discovery_rows.empty),
            "positive_source_count": int(group["source"].nunique()),
            "clinical_sources": ";".join(sorted(clinical_rows["source"].dropna().astype(str).unique())),
            "discovery_sources": ";".join(sorted(discovery_rows["source"].dropna().astype(str).unique())),
            "heldout_eligible_sources": ";".join(sorted(group["source"].dropna().astype(str).unique())),
        })
    return labels, pd.DataFrame(rows)


def build_feature_matrix(matrix: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = pd.DataFrame({"gene_symbol": matrix["gene_symbol"].astype(str)})
    manifest_rows: list[dict] = []
    missing_features = set(cfg.get("missing_indicator_features", []))
    cancer_code = str(cfg.get("cancer_code", "PDAC")).upper()
    gate_features = {"dependency_score", "druggability_score", "safety_risk_score", "safety_factor"}
    trans_excluded = set(cfg.get("evaluation", {}).get("axes", {}).get("translational", {}).get("feature_exclusions", []))
    func_excluded = set(cfg.get("evaluation", {}).get("axes", {}).get("functional", {}).get("feature_exclusions", []))
    for feat in cfg["features"]:
        if feat == "safety_factor" and feat not in matrix.columns:
            values = 1.0 - _numeric(matrix, "safety_risk_score", 1.0)
        else:
            values = _numeric(matrix, feat, np.nan)
        if feat in missing_features:
            out[f"{feat}_missing"] = values.isna().astype(int)
            manifest_rows.append({
                "feature": f"{feat}_missing",
                "source_family": FEATURE_DETAILS.get(feat, ("unspecified", "bounded feature"))[0],
                "record_role": "F",
                "used_in_review_gate": 0,
                "transformation": f"1 if {feat} is unavailable, otherwise 0",
                "missing_value_rule": "indicator is always observed",
                "entered_translational_model": int(feat not in trans_excluded),
                "entered_functional_model": int(feat not in func_excluded),
                "entered_model": int(feat not in trans_excluded or feat not in func_excluded),
                "available": int(feat in matrix.columns),
            })
        out[feat] = values.fillna(0).clip(0, 1)
        transformation = FEATURE_DETAILS.get(feat, ("unspecified", "precomputed bounded score"))[1]
        if feat == "external_reproducibility_score" and cancer_code == "LUAD":
            transformation = "0.5 * proteomic_support_score + 0.5 * dependency_score"
        manifest_rows.append({
            "feature": feat,
            "source_family": FEATURE_DETAILS.get(feat, ("unspecified", transformation))[0],
            "record_role": "F",
            "used_in_review_gate": int(feat in gate_features),
            "transformation": transformation,
            "missing_value_rule": "explicit missing indicator where configured; numeric value set to 0 before scaling",
            "entered_translational_model": int(feat not in trans_excluded),
            "entered_functional_model": int(feat not in func_excluded),
            "entered_model": int(feat not in trans_excluded or feat not in func_excluded),
            "available": int(feat in matrix.columns),
        })
    return out, pd.DataFrame(manifest_rows)


def build_source_role_audit(role_manifest: pd.DataFrame, feature_manifest: pd.DataFrame) -> pd.DataFrame:
    """Summarize explicit record roles and feature-family uses."""
    label_rows = []
    for source, group in role_manifest.groupby("source", sort=True):
        roles = ",".join(sorted(group["record_role"].unique()))
        axes = ",".join(sorted(group["axis"].unique()))
        label_rows.append({
            "source": source,
            "roles": roles,
            "framework_use": f"label/benchmark records on {axes} axis",
            "leakage_control": "V genes embargoed from L; record_id has one role per axis",
        })
    for source, group in feature_manifest.groupby("source_family", sort=True):
        roles = ",".join(sorted(group["record_role"].unique()))
        label_rows.append({
            "source": source,
            "roles": roles,
            "framework_use": ";".join(group["feature"].astype(str)),
            "leakage_control": "axis-specific feature exclusions are recorded in the feature manifest",
        })
    return pd.DataFrame(label_rows)


def build_method_specifications(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build machine-readable sampling, validation and review-rule specifications."""
    bg = cfg["background_sampling"]
    model = cfg["model"]
    validation = cfg["validation"]
    method_rows = [
        ("candidate matching", "database-count bin | DepMap-measured flag | expression-coverage bin | symbol-validity flag"),
        ("matching fallback", "fill unmatched stratum quota from remaining eligible unlabeled genes without replacement"),
        ("background ratio", str(bg["background_ratio"])),
        ("development repeats", str(bg["n_repeats"])),
        ("permutation repeats", str(validation["permutation_repeats"])),
        ("label permutations", str(validation["n_permutations"])),
        ("scaling", "StandardScaler fitted within each sampled training set"),
        ("classifier", "class-weight-balanced elastic-net logistic regression"),
        ("C", str(model["C"])),
        ("l1_ratio", str(model["l1_ratio"])),
        ("solver", str(model["solver"])),
        ("max_iter", str(model["max_iter"])),
        ("convergence tolerance", "0.01 for elastic-net SAGA; 0.0001 for L2 liblinear"),
        ("convergence replacement", "exclude fits reaching max_iter and deterministically replace them until the requested number of converged repeats is obtained"),
    ]
    review = cfg["review"]
    rule_rows = [
        (1, "Class IV", f"safety_risk_score >= {review['safety_risk_high_threshold']}", "safety-risk gate overrides all other classes"),
        (2, "Class I", "translational score in top percentile AND druggability >= 0.1 AND safety factor passes AND dependency/override passes", "translational-precedence candidate"),
        (3, "Class III", f"dependency passes AND druggability < {review['limited_actionability_druggability_max']}", "limited-actionability gate precedes Class II"),
        (4, "Class II", "functional score, dependency, disease relevance and safety factor all pass", "functional-discovery candidate"),
        (5, "Class V", "disease relevance passes AND dependency is below half-threshold", "biomarker-like or low-dependency candidate"),
        (6, "Class VI", "none of the preceding rules pass", "insufficient or mixed evidence"),
    ]
    statistical_rows = [
        ("primary validation endpoint", "held-out benchmark top-50 recovery"),
        ("OOF endpoint", "held-out development-positive rank and top-K recovery"),
        ("empirical p-value", "(1 + number of permuted statistics >= observed) / (1 + number of permutations)"),
        ("multiple testing", "Holm adjustment across the four cancer-by-axis permutation tests"),
        ("benchmark uncertainty", "target-level bootstrap 95% confidence interval for top-K recovery"),
    ]
    return (
        pd.DataFrame(method_rows, columns=["setting", "specification"]),
        pd.DataFrame(rule_rows, columns=["priority", "class", "executable_rule", "interpretation"]),
        pd.DataFrame(statistical_rows, columns=["item", "specification"]),
    )


def build_background_strata(matrix: pd.DataFrame) -> pd.Series:
    idx = matrix["gene_symbol"].astype(str)
    db_cols = [
        "opentargets_overall_score", "drugcentral_interaction_count", "dgidb_interaction_count",
        "ttd_mapping_count", "iuphar_approved_interaction_count", "clinical_precedence_score",
    ]
    database_coverage = pd.Series(0, index=matrix.index, dtype=float)
    for col in db_cols:
        if col in matrix.columns:
            database_coverage += _numeric(matrix, col).fillna(0).gt(0).astype(int)
    database_bin = pd.cut(database_coverage, bins=[-1, 0, 2, 99], labels=["db0", "db1_2", "db3p"])

    depmap_measured = _numeric(matrix, "depmap_pdac_median_gene_effect", np.nan).notna().map({True: "depmap_y", False: "depmap_n"})
    expr_rank = _numeric(matrix, "tcga_paad_median_log2tpm", 0).rank(method="first")
    expr_bin = pd.qcut(expr_rank, q=4, labels=["expr_q1", "expr_q2", "expr_q3", "expr_q4"], duplicates="drop")
    symbol_bin = idx.str.match(r"^[A-Z][A-Z0-9.-]{1,15}$").map({True: "symbol_like", False: "symbol_other"})
    strata = (
        database_bin.astype(str).fillna("db_missing") + "|" +
        depmap_measured.astype(str).fillna("depmap_missing") + "|" +
        expr_bin.astype(str).fillna("expr_missing") + "|" +
        symbol_bin.astype(str).fillna("symbol_missing")
    )
    return pd.Series(strata.to_numpy(), index=idx, name="background_stratum")


def make_elastic_net_model(C: float, l1_ratio: float, max_iter: int, solver: str = "saga") -> Pipeline:
    return Pipeline([
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(
            solver=solver,
            C=C,
            l1_ratio=l1_ratio,
            max_iter=max_iter,
            tol=1e-2,
            class_weight="balanced",
            random_state=0,
        )),
    ])


def make_l2_model(C: float, max_iter: int, solver: str = "liblinear") -> Pipeline:
    return Pipeline([
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(
            solver="liblinear",
            C=C,
            max_iter=max_iter,
            tol=1e-4,
            class_weight="balanced",
            random_state=0,
        )),
    ])


def _matched_background_sample(
    positives: list[str],
    unlabeled: list[str],
    strata: pd.Series,
    ratio: int,
    rng: np.random.Generator,
) -> list[str]:
    selected: list[str] = []
    used: set[str] = set()
    unlabeled_set = set(unlabeled)
    stratum_pools = {
        stratum: [g for g in members.index.astype(str) if g in unlabeled_set]
        for stratum, members in strata.groupby(strata)
    }
    positive_strata = strata.reindex(positives).fillna("missing").astype(str)
    for stratum, n_pos in positive_strata.value_counts().items():
        pool = stratum_pools.get(stratum, [])
        if not pool:
            continue
        draw_n = min(int(n_pos) * ratio, len(pool))
        draw = rng.choice(pool, size=draw_n, replace=False).tolist()
        selected.extend(draw)
        used.update(draw)
    target_n = min(len(unlabeled), max(len(positives) * ratio, len(positives) + 1))
    if len(selected) < target_n:
        remaining = [g for g in unlabeled if g not in used]
        if remaining:
            selected.extend(rng.choice(remaining, size=min(target_n - len(selected), len(remaining)), replace=False).tolist())
    return selected[:target_n]


def fit_sampled_background_model(
    X: pd.DataFrame,
    positive_genes: set[str],
    all_positive_genes: set[str],
    seeds: list[int],
    cfg: dict,
    strata: pd.Series | None,
    matched_background: bool,
    background_ratio: int,
    penalty: str = "elasticnet",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Fit a fixed-size ensemble of converged sampled-background models.

    Parameters
    ----------
    X:
        Gene-by-feature matrix indexed by gene symbol.
    positive_genes:
        Development positives for the fitted axis or resampled null set.
    all_positive_genes:
        Genes excluded from background sampling, including benchmark embargoes.
    seeds:
        Deterministic seeds defining the requested ensemble size.
    cfg:
        Framework configuration containing model hyperparameters.
    strata:
        Optional matched-background stratum for every candidate gene.
    matched_background:
        Whether to sample background genes within positive strata.
    background_ratio:
        Requested number of unlabeled genes per positive gene.
    penalty:
        ``elasticnet`` for the primary model or ``l2`` for a comparator.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame]
        Ensemble scores, accepted-fit coefficients, and an attempt-level
        convergence audit. A fit that reaches ``max_iter`` is excluded and
        deterministically replaced so the ensemble size remains fixed.
    """
    model_cfg = cfg["model"]
    genes = X.index.astype(str)
    positives = [g for g in genes if g in positive_genes]
    unlabeled = [g for g in genes if g not in all_positive_genes]
    if len(positives) < 2 or len(unlabeled) < 10:
        raise ValueError("Not enough positives or unlabeled genes")

    target_repeats = len(seeds)
    seed_queue = list(map(int, seeds))
    next_seed = (max(seed_queue) + 1) if seed_queue else 1
    max_attempts = max(target_repeats * 10, target_repeats + 100)
    pred_runs: list[pd.Series] = []
    coef_runs: list[pd.DataFrame] = []
    train_rows: list[dict] = []
    attempt = 0
    while len(pred_runs) < target_repeats and attempt < max_attempts:
        # Additional seeds are consecutive and deterministic, making numerical
        # replacement exactly reproducible without accepting a failed fit.
        if attempt >= len(seed_queue):
            seed_queue.append(next_seed)
            next_seed += 1
        seed = seed_queue[attempt]
        attempt += 1
        rng = np.random.default_rng(seed)
        if matched_background and strata is not None:
            bg = _matched_background_sample(positives, unlabeled, strata, background_ratio, rng)
        else:
            n_bg = min(len(unlabeled), max(len(positives) * background_ratio, len(positives) + 1))
            bg = rng.choice(unlabeled, size=n_bg, replace=False).tolist()
        train_genes = positives + bg
        y = np.array([1] * len(positives) + [0] * len(bg))
        if penalty == "elasticnet":
            model = make_elastic_net_model(
                C=float(model_cfg.get("C", 1.0)),
                l1_ratio=float(model_cfg.get("l1_ratio", 0.5)),
                max_iter=int(model_cfg.get("max_iter", 3000)),
                solver=str(model_cfg.get("solver", "saga")),
            )
        else:
            model = make_l2_model(C=float(model_cfg.get("C", 1.0)), max_iter=int(model_cfg.get("max_iter", 3000)))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            model.fit(X.loc[train_genes], y)
        iterations = int(model.named_steps["clf"].n_iter_[0])
        converged = iterations < int(model_cfg.get("max_iter", 3000))
        if converged:
            pred_runs.append(pd.Series(model.predict_proba(X)[:, 1], index=genes, name=f"seed_{seed}"))
            coef_runs.append(pd.DataFrame({
                "feature": X.columns,
                "coefficient": model.named_steps["clf"].coef_[0],
                "seed": seed,
            }))
        train_rows.append({
            "seed": seed,
            "attempt": attempt,
            "n_positive": len(positives),
            "n_background": len(bg),
            "background_ratio": background_ratio,
            "matched_background": int(bool(matched_background and strata is not None)),
            "penalty": penalty,
            "C": float(model_cfg.get("C", 1.0)),
            "l1_ratio": float(model_cfg.get("l1_ratio", 0.5)) if penalty == "elasticnet" else np.nan,
            "iterations": iterations,
            "converged": int(converged),
            "used_in_ensemble": int(converged),
        })

    if len(pred_runs) < target_repeats:
        raise RuntimeError(
            f"Only {len(pred_runs)} of {target_repeats} requested fits converged after {attempt} attempts"
        )

    pred_df = pd.concat(pred_runs, axis=1)
    ranks = pred_df.rank(ascending=False, method="min")
    scores = pd.DataFrame({
        "gene_symbol": pred_df.index,
        "score_mean": pred_df.mean(axis=1).to_numpy(),
        "score_sd": pred_df.std(axis=1).fillna(0).to_numpy(),
        "top20_frequency": (ranks <= 20).mean(axis=1).to_numpy(),
        "top50_frequency": (ranks <= 50).mean(axis=1).to_numpy(),
        "top100_frequency": (ranks <= 100).mean(axis=1).to_numpy(),
    })
    scores["rank"] = scores["score_mean"].rank(ascending=False, method="min").astype(int)
    return scores, pd.concat(coef_runs, ignore_index=True), pd.DataFrame(train_rows)


def fit_pseudo_negative_model(
    X: pd.DataFrame,
    positive_genes: set[str],
    all_positive_genes: set[str],
    cfg: dict,
) -> pd.Series:
    genes = X.index.astype(str)
    positives = [g for g in genes if g in positive_genes]
    unlabeled = [g for g in genes if g not in all_positive_genes]
    train_genes = positives + unlabeled
    y = np.array([1] * len(positives) + [0] * len(unlabeled))
    model_cfg = cfg["model"]
    model = make_l2_model(C=float(model_cfg.get("C", 1.0)), max_iter=int(model_cfg.get("max_iter", 3000)))
    model.fit(X.loc[train_genes], y)
    return pd.Series(model.predict_proba(X)[:, 1], index=genes)


def train_dual_axis_models(
    feature_matrix: pd.DataFrame,
    labels: pd.DataFrame,
    benchmark_sets: dict[str, set[str]],
    matrix: pd.DataFrame,
    cfg: dict,
    out_tables: Path,
) -> pd.DataFrame:
    """Fit final development models while keeping benchmark genes out of training."""
    X_all = _safe_gene_index(feature_matrix)
    labels = _safe_gene_index(labels)
    all_benchmarks = set().union(*benchmark_sets.values())
    bg_cfg = cfg["background_sampling"]
    seeds = cfg["random_seeds"][: int(bg_cfg["n_repeats"])]
    strata = build_background_strata(matrix) if bool(bg_cfg.get("matched_background", True)) else None
    score_parts: list[pd.DataFrame] = []
    coef_parts: list[pd.DataFrame] = []
    train_parts: list[pd.DataFrame] = []
    for axis in AXES:
        columns = axis_feature_columns(feature_matrix, cfg, axis.name)
        X = X_all[columns]
        pos = set(labels.index[pd.to_numeric(labels[axis.label_column], errors="coerce").fillna(0).astype(int).eq(1)])
        axis_reserved = pos | all_benchmarks
        score, coef, train = fit_sampled_background_model(
            X,
            pos,
            axis_reserved,
            seeds,
            cfg,
            strata,
            bool(bg_cfg.get("matched_background", True)),
            int(bg_cfg.get("background_ratio", 8)),
            penalty="elasticnet",
        )
        score = score.rename(columns={
            "score_mean": axis.score_column,
            "score_sd": f"{axis.score_column}_sd",
            "rank": axis.rank_column,
            "top20_frequency": f"{axis.top_prefix}_top20_frequency",
            "top50_frequency": f"{axis.top_prefix}_top50_frequency",
            "top100_frequency": f"{axis.top_prefix}_top100_frequency",
        })
        score_parts.append(score)
        coef["axis"] = axis.name
        train["axis"] = axis.name
        train["n_model_features"] = len(columns)
        train["n_reserved_benchmark_genes"] = len(benchmark_sets[axis.name])
        coef_parts.append(coef)
        train_parts.append(train)

    scores = score_parts[0].merge(score_parts[1], on="gene_symbol", how="outer")
    coef_all = pd.concat(coef_parts, ignore_index=True)
    coef_summary = coef_all.groupby(["axis", "feature"], as_index=False).agg(
        coefficient_mean=("coefficient", "mean"),
        coefficient_sd=("coefficient", "std"),
        nonzero_frequency=("coefficient", lambda x: float((x.abs() > 1e-9).mean())),
    )
    write_tsv(scores, out_tables / "framework_priority_scores.tsv")
    write_tsv(coef_summary, out_tables / "framework_feature_coefficients.tsv")
    write_tsv(pd.concat(train_parts, ignore_index=True), out_tables / "framework_training_runs.tsv")
    return scores


def assign_review_classes(matrix: pd.DataFrame, scores: pd.DataFrame, cfg: dict, out_tables: Path) -> pd.DataFrame:
    review = cfg["review"]
    df = _safe_gene_index(matrix).join(_safe_gene_index(scores), how="left")
    clinical_cut = df[CLINICAL_SCORE].quantile(1 - float(review["high_score_percentile"]))
    discovery_cut = df[DISCOVERY_SCORE].quantile(1 - float(review["high_score_percentile"]))
    dep_cut = df["dependency_score"].quantile(1 - float(review["dependency_top_percentile"]))
    disease_source = df["disease_relevance_score"] if "disease_relevance_score" in df.columns else df["patient_alteration_score"]
    disease_cut = disease_source.quantile(1 - float(review["disease_top_percentile"]))
    safety_min = float(review["safety_factor_min"])
    high_safety = float(review["safety_risk_high_threshold"])
    translational_override = bool(review.get("allow_strong_translational_evidence_without_dependency", False))
    clinical_override_min = float(review.get("clinical_precedence_override_min", 0.6))
    opentargets_override_min = float(review.get("opentargets_override_min", 0.5))
    limited_actionability_first = bool(review.get("limited_actionability_before_functional_class", False))
    limited_actionability_drug_max = float(review.get("limited_actionability_druggability_max", 0.15))

    classes: list[str] = []
    for _, row in df.iterrows():
        safety = float(row.get("safety_risk_score", 0) or 0)
        sf = float(row.get("safety_factor", 1 - safety) or 0)
        dep = float(row.get("dependency_score", 0) or 0)
        disease = float(row.get("disease_relevance_score", row.get("patient_alteration_score", 0)) or 0)
        drug = float(row.get("druggability_score", 0) or 0)
        strong_translational_support = (
            float(row.get("clinical_precedence_score", 0) or 0) >= clinical_override_min
            or float(row.get("opentargets_overall_score", 0) or 0) >= opentargets_override_min
        )
        clinical_gate = (
            row[CLINICAL_SCORE] >= clinical_cut
            and drug >= 0.1
            and sf >= safety_min
            and (dep >= dep_cut * 0.5 or (translational_override and strong_translational_support))
        )
        limited_actionability_gate = dep >= dep_cut and drug < limited_actionability_drug_max
        functional_gate = (
            row[DISCOVERY_SCORE] >= discovery_cut
            and dep >= dep_cut * 0.5
            and disease >= disease_cut * 0.5
            and sf >= safety_min
        )
        if safety >= high_safety:
            cls = "Class IV: safety-risk candidate"
        elif clinical_gate:
            cls = "Class I: translational-precedence candidate"
        elif limited_actionability_first and limited_actionability_gate:
            cls = "Class III: dependency-supported candidate with limited actionability"
        elif functional_gate:
            cls = "Class II: functional-discovery candidate"
        elif limited_actionability_gate:
            cls = "Class III: dependency-supported candidate with limited actionability"
        elif disease >= disease_cut and dep < dep_cut * 0.5:
            cls = "Class V: biomarker-like or low-dependency candidate"
        else:
            cls = "Class VI: insufficient or mixed evidence"
        classes.append(cls)

    cols = [
        CLINICAL_SCORE, "translational_precedence_rank", DISCOVERY_SCORE, "functional_discovery_rank",
        "dependency_score", "druggability_score", "disease_relevance_score", "safety_risk_score", "safety_factor",
    ]
    out = df.reset_index()[["gene_symbol"] + [c for c in cols if c in df.columns]].copy()
    out["review_class"] = classes
    write_tsv(out, out_tables / "framework_review_classes.tsv")
    write_tsv(pd.DataFrame([
        {"threshold": "translational_high_score", "value": clinical_cut},
        {"threshold": "functional_high_score", "value": discovery_cut},
        {"threshold": "dependency_reference", "value": dep_cut},
        {"threshold": "disease_reference", "value": disease_cut},
        {"threshold": "safety_factor_min", "value": safety_min},
        {"threshold": "safety_risk_high", "value": high_safety},
    ]), out_tables / "framework_review_thresholds.tsv")
    return out


def _top_set(scores: pd.Series, k: int) -> set[str]:
    return set(scores.sort_values(ascending=False).head(k).index.astype(str))


def _metric_row(method: str, axis: str, scores: pd.Series, positives: set[str], review_classes: pd.DataFrame, k: int, universe_n: int) -> dict:
    top = _top_set(scores, k)
    sub = _safe_gene_index(review_classes).reindex(list(top))
    classes = sub["review_class"].astype(str)
    recovery = len(top & positives) / max(len(positives), 1)
    positive_ranks = scores.rank(ascending=False, method="min").reindex(sorted(positives)).dropna()
    return {
        "method": method,
        "axis": axis,
        "top_k": k,
        "positive_recovery": recovery,
        "precision_at_k": len(top & positives) / max(k, 1),
        "enrichment_over_random": recovery / max(k / universe_n, 1e-12),
        "median_rank_percentile": float(1 - (positive_ranks.median() - 1) / max(universe_n - 1, 1)) if len(positive_ranks) else 0.0,
        "mean_rank_percentile": float(1 - (positive_ranks.mean() - 1) / max(universe_n - 1, 1)) if len(positive_ranks) else 0.0,
        "class_I_II_fraction": float(classes.str.contains("Class I:|Class II:").mean()),
        "high_safety_risk_fraction": float(pd.to_numeric(sub.get("safety_risk_score"), errors="coerce").fillna(0).ge(0.75).mean()),
        "biomarker_like_fraction": float(classes.str.contains("Class V:").mean()),
        "median_safety_risk": float(pd.to_numeric(sub.get("safety_risk_score"), errors="coerce").median()),
        "median_dependency_score": float(pd.to_numeric(sub.get("dependency_score"), errors="coerce").median()),
    }


def _component_score(matrix_idx: pd.DataFrame, cols: Iterable[str]) -> pd.Series:
    existing = [c for c in cols if c in matrix_idx.columns]
    if not existing:
        return pd.Series(0.0, index=matrix_idx.index)
    return matrix_idx[existing].apply(pd.to_numeric, errors="coerce").fillna(0).clip(0, 1).mean(axis=1)


def _rank_aggregation(matrix_idx: pd.DataFrame, cols: Iterable[str]) -> pd.Series:
    existing = [c for c in cols if c in matrix_idx.columns]
    if not existing:
        return pd.Series(0.0, index=matrix_idx.index)
    ranks = matrix_idx[existing].apply(pd.to_numeric, errors="coerce").fillna(0).rank(ascending=False, pct=True)
    return 1 - ranks.mean(axis=1)


def run_comparators(
    feature_matrix: pd.DataFrame,
    matrix: pd.DataFrame,
    labels: pd.DataFrame,
    benchmark_sets: dict[str, set[str]],
    scores: pd.DataFrame,
    review_classes: pd.DataFrame,
    cfg: dict,
    out_tables: Path,
) -> pd.DataFrame:
    """Evaluate all routes against the same embargoed benchmark targets."""
    X = _safe_gene_index(feature_matrix)
    matrix_idx = _safe_gene_index(matrix)
    labels_idx = _safe_gene_index(labels)
    scores_idx = _safe_gene_index(scores)
    all_benchmarks = set().union(*benchmark_sets.values())
    positives = {
        "translational": set(labels_idx.index[pd.to_numeric(labels_idx["clinical_positive"], errors="coerce").fillna(0).astype(int).eq(1)]),
        "functional": set(labels_idx.index[pd.to_numeric(labels_idx["discovery_positive"], errors="coerce").fillna(0).astype(int).eq(1)]),
    }
    combined_pos = positives["translational"] | positives["functional"]
    bg_cfg = cfg["background_sampling"]
    short_cfg = {**cfg, "random_seeds": cfg["random_seeds"][:10]}
    strata = build_background_strata(matrix) if bool(bg_cfg.get("matched_background", True)) else None
    combined_columns = sorted(set(axis_feature_columns(feature_matrix, cfg, "translational")) & set(axis_feature_columns(feature_matrix, cfg, "functional")))
    combined_score, _, _ = fit_sampled_background_model(
        X[combined_columns], combined_pos, combined_pos | all_benchmarks, cfg["random_seeds"][: int(bg_cfg["n_repeats"])], cfg, strata,
        bool(bg_cfg.get("matched_background", True)), int(bg_cfg.get("background_ratio", 8)), penalty="elasticnet",
    )
    axis_X = {axis: X[axis_feature_columns(feature_matrix, cfg, axis)] for axis in positives}
    pseudo_trans = fit_pseudo_negative_model(axis_X["translational"], positives["translational"], positives["translational"] | all_benchmarks, cfg)
    pseudo_func = fit_pseudo_negative_model(axis_X["functional"], positives["functional"], positives["functional"] | all_benchmarks, cfg)
    l2_scores = {}
    for axis_name in positives:
        score, _, _ = fit_sampled_background_model(
            axis_X[axis_name], positives[axis_name], positives[axis_name] | all_benchmarks,
            cfg["random_seeds"][: int(bg_cfg["n_repeats"])], cfg, strata,
            bool(bg_cfg.get("matched_background", True)), int(bg_cfg.get("background_ratio", 8)), penalty="l2",
        )
        l2_scores[axis_name] = score.set_index("gene_symbol")["score_mean"]

    method_scores = {
        ("Dual-axis sampled-background PU", "translational"): scores_idx[CLINICAL_SCORE],
        ("Dual-axis sampled-background PU", "functional"): scores_idx[DISCOVERY_SCORE],
        ("Pseudo-negative L2 logistic", "translational"): pseudo_trans,
        ("Pseudo-negative L2 logistic", "functional"): pseudo_func,
        ("Repeated sampled-background L2", "translational"): l2_scores["translational"],
        ("Repeated sampled-background L2", "functional"): l2_scores["functional"],
        ("Combined single-axis sampled-background PU", "translational"): combined_score.set_index("gene_symbol")["score_mean"],
        ("Combined single-axis sampled-background PU", "functional"): combined_score.set_index("gene_symbol")["score_mean"],
        ("Rank aggregation", "translational"): _rank_aggregation(matrix_idx, ["clinical_precedence_score", "chemical_modulation_score", "druggability_score", "opentargets_overall_score"]),
        ("Rank aggregation", "functional"): _rank_aggregation(matrix_idx, ["dependency_score", "patient_alteration_score", "proteomic_support_score", "external_reproducibility_score"]),
        ("DEG-only", "translational"): _component_score(matrix_idx, ["cptac_mrna_score", "tcga_gtex_expression_score", "geo_external_reproducibility_score"]),
        ("Protein-only", "translational"): _component_score(matrix_idx, ["cptac_protein_score", "phosphoproteomic_support_score"]),
        ("Dependency-only", "functional"): _component_score(matrix_idx, ["dependency_score"]),
        ("Druggability-only", "translational"): _component_score(matrix_idx, ["druggability_score"]),
        ("OpenTargets-only", "translational"): _component_score(matrix_idx, ["opentargets_overall_score"]),
        ("Formula clinical-actionability score", "translational"): _component_score(matrix_idx, ["clinical_actionability_score"]),
        ("Formula discovery-potential score", "functional"): _component_score(matrix_idx, ["discovery_potential_score"]),
    }
    rows = []
    for (method, axis), s in method_scores.items():
        s = s.reindex(matrix_idx.index).fillna(0)
        for k in cfg["validation"]["top_k"]:
            rows.append(_metric_row(method, axis, s, benchmark_sets[axis], review_classes, int(k), len(matrix_idx)))
    out = pd.DataFrame(rows)
    write_tsv(out, out_tables / "framework_comparator_summary.tsv")

    primary_t50 = _top_set(scores_idx[CLINICAL_SCORE], 50)
    primary_f50 = _top_set(scores_idx[DISCOVERY_SCORE], 50)
    combined_t50 = _top_set(combined_score.set_index("gene_symbol")["score_mean"], 50)
    write_tsv(pd.DataFrame([{
        "comparison": "combined_single_axis_vs_dual_axis_top50",
        "combined_overlap_with_translational_top50": len(combined_t50 & primary_t50) / 50,
        "combined_overlap_with_functional_top50": len(combined_t50 & primary_f50) / 50,
        "combined_recovered_translational_benchmarks_top50": len(combined_t50 & benchmark_sets["translational"]),
        "combined_recovered_functional_benchmarks_top50": len(combined_t50 & benchmark_sets["functional"]),
        "dual_functional_recovered_functional_benchmarks_top50": len(primary_f50 & benchmark_sets["functional"]),
    }]), out_tables / "framework_combined_axis_overlap.tsv")
    return out


def _recovery(scores: pd.Series, heldout: set[str], topks: list[int]) -> dict:
    ranks = scores.rank(ascending=False, method="min")
    n = max(1, len(heldout))
    out = {"n_heldout_targets": len(heldout)}
    for k in topks:
        top = _top_set(scores, int(k))
        rec = len(top & heldout) / n
        out[f"top{k}_recovery"] = rec
        out[f"enrichment_top{k}_over_random"] = rec / max(int(k) / max(len(scores), 1), 1e-12)
        out[f"precision_at_{k}"] = len(top & heldout) / max(int(k), 1)
    hr = ranks.loc[list(heldout & set(ranks.index))] if heldout else pd.Series(dtype=float)
    out["median_rank_percentile"] = float(1 - (hr.median() - 1) / max(len(scores) - 1, 1)) if len(hr) else 0.0
    out["mean_rank_percentile"] = float(1 - (hr.mean() - 1) / max(len(scores) - 1, 1)) if len(hr) else 0.0
    return out


def _bootstrap_recovery_ci(scores: pd.Series, targets: set[str], k: int, seed: int = 91) -> tuple[float, float]:
    """Return a target-level bootstrap interval for recovery at ``k``."""
    target_list = sorted(targets & set(scores.index.astype(str)))
    if not target_list:
        return 0.0, 0.0
    top = _top_set(scores, k)
    hits = np.array([int(gene in top) for gene in target_list], dtype=float)
    rng = np.random.default_rng(seed)
    means = np.array([rng.choice(hits, size=len(hits), replace=True).mean() for _ in range(2000)])
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def run_independent_benchmark_validation(
    scores: pd.DataFrame,
    benchmark_sets: dict[str, set[str]],
    cfg: dict,
    out_tables: Path,
) -> pd.DataFrame:
    """Evaluate final development models on genes never used as labels."""
    score_idx = _safe_gene_index(scores)
    rows = []
    for axis in AXES:
        axis_scores = score_idx[axis.score_column]
        metrics = _recovery(axis_scores, benchmark_sets[axis.name], [int(k) for k in cfg["validation"]["top_k"]])
        low, high = _bootstrap_recovery_ci(axis_scores, benchmark_sets[axis.name], 50)
        rows.append({
            "axis": axis.name,
            "evaluation_set": "independent_benchmark",
            **metrics,
            "top50_recovery_ci_low": low,
            "top50_recovery_ci_high": high,
        })
    out = pd.DataFrame(rows)
    write_tsv(out, out_tables / "framework_independent_benchmark_validation.tsv")
    return out


def _positive_folds(genes: set[str], role_manifest: pd.DataFrame, axis_name: str, strategy: str) -> list[list[str]]:
    """Create deterministic held-out positive folds without using outcomes."""
    ordered = sorted(genes)
    if strategy == "leave_one_positive_out" or len(ordered) <= 10:
        return [[gene] for gene in ordered]
    source_rows = role_manifest[(role_manifest["axis"].eq(axis_name)) & (role_manifest["record_role"].eq("L"))]
    source_map = source_rows.groupby("gene_symbol")["source"].agg(lambda x: sorted(set(map(str, x)))[0]).to_dict()
    n_folds = min(5, len(ordered))
    folds: list[list[str]] = [[] for _ in range(n_folds)]
    for _, group in pd.DataFrame({"gene": ordered, "source": [source_map.get(g, "unknown") for g in ordered]}).groupby("source", sort=True):
        for index, gene in enumerate(sorted(group["gene"])):
            folds[index % n_folds].append(gene)
    return [fold for fold in folds if fold]


def run_out_of_fold_validation(
    feature_matrix: pd.DataFrame,
    positive_sets: pd.DataFrame,
    benchmark_sets: dict[str, set[str]],
    role_manifest: pd.DataFrame,
    matrix: pd.DataFrame,
    cfg: dict,
    out_tables: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate held-out scores for every development positive."""
    X_all = _safe_gene_index(feature_matrix)
    labels = _safe_gene_index(positive_sets)
    bg = cfg["background_sampling"]
    seeds = cfg["random_seeds"][: int(bg["n_repeats"])]
    strata = build_background_strata(matrix) if bool(bg.get("matched_background", True)) else None
    detail_rows = []
    summary_rows = []
    all_benchmarks = set().union(*benchmark_sets.values())
    for axis in AXES:
        positives = set(labels.index[pd.to_numeric(labels[axis.label_column], errors="coerce").fillna(0).astype(int).eq(1)])
        strategy = cfg["evaluation"]["axes"][axis.name].get("oof_strategy", "source_stratified_5fold")
        folds = _positive_folds(positives, role_manifest, axis.name, strategy)
        X = X_all[axis_feature_columns(feature_matrix, cfg, axis.name)]
        axis_details = []
        for fold_id, heldout_genes in enumerate(folds, start=1):
            heldout = set(heldout_genes)
            train_pos = positives - heldout
            score, _, _ = fit_sampled_background_model(
                X, train_pos, positives | all_benchmarks, seeds, cfg, strata,
                bool(bg.get("matched_background", True)), int(bg.get("background_ratio", 8)), penalty="elasticnet",
            )
            s = score.set_index("gene_symbol")["score_mean"]
            ranks = s.rank(ascending=False, method="min")
            for gene in sorted(heldout):
                rank = float(ranks.get(gene, np.nan))
                axis_details.append({
                    "axis": axis.name,
                    "fold": fold_id,
                    "gene_symbol": gene,
                    "score": float(s.get(gene, np.nan)),
                    "rank": rank,
                    "rank_percentile": float(1 - (rank - 1) / max(len(s) - 1, 1)) if np.isfinite(rank) else np.nan,
                    "n_training_positives": len(train_pos),
                    "n_repeats": len(seeds),
                })
        detail_rows.extend(axis_details)
        detail = pd.DataFrame(axis_details)
        row = {
            "axis": axis.name,
            "strategy": strategy,
            "n_heldout_targets": len(detail),
            "median_rank_percentile": float(detail["rank_percentile"].median()),
            "mean_rank_percentile": float(detail["rank_percentile"].mean()),
        }
        for k in cfg["validation"]["top_k"]:
            row[f"top{k}_recovery"] = float(detail["rank"].le(int(k)).mean())
        summary_rows.append(row)
    detail_out = pd.DataFrame(detail_rows)
    summary_out = pd.DataFrame(summary_rows)
    write_tsv(detail_out, out_tables / "framework_oof_positive_scores.tsv")
    write_tsv(summary_out, out_tables / "framework_oof_validation_summary.tsv")
    return detail_out, summary_out


def run_leave_source_out(
    feature_matrix: pd.DataFrame,
    role_manifest: pd.DataFrame,
    positive_sets: pd.DataFrame,
    benchmark_sets: dict[str, set[str]],
    matrix: pd.DataFrame,
    cfg: dict,
    out_tables: Path,
) -> pd.DataFrame:
    X_all = _safe_gene_index(feature_matrix)
    labels_idx = _safe_gene_index(positive_sets)
    repeats = int(cfg["validation"]["leave_source_out_repeats"])
    seeds = cfg["random_seeds"][:repeats]
    bg_cfg = cfg["background_sampling"]
    strata = build_background_strata(matrix) if bool(bg_cfg.get("matched_background", True)) else None
    rows = []
    all_benchmarks = set().union(*benchmark_sets.values())
    for axis in AXES:
        axis_records = role_manifest[(role_manifest["axis"].eq(axis.name)) & (role_manifest["record_role"].eq("L"))]
        axis_pos = set(labels_idx.index[pd.to_numeric(labels_idx[axis.label_column], errors="coerce").fillna(0).astype(int).eq(1)])
        for source in sorted(axis_records["source"].dropna().unique()):
            heldout = set(axis_records.loc[axis_records["source"].eq(source), "gene_symbol"].astype(str))
            train_pos = axis_pos - heldout
            if len(heldout) < 2:
                continue
            if len(train_pos) < 2:
                rows.append({"heldout_source": source, "axis": axis.name, "status": "NOT_ESTIMABLE", "reason": "fewer than two positives remain"})
                continue
            excluded = set(cfg.get("evaluation", {}).get("source_feature_exclusions", {}).get(source, []))
            base_columns = axis_feature_columns(feature_matrix, cfg, axis.name)
            columns = [column for column in base_columns if column not in excluded and not any(column == f"{item}_missing" for item in excluded)]
            X = X_all[columns]
            score, _, _ = fit_sampled_background_model(
                X, train_pos, axis_pos | all_benchmarks, seeds, cfg, strata,
                bool(bg_cfg.get("matched_background", True)), int(bg_cfg.get("background_ratio", 8)), penalty="elasticnet",
            )
            s = score.set_index("gene_symbol")["score_mean"]
            row = {"heldout_source": source, "axis": axis.name, "status": "ESTIMATED", "reason": "", "excluded_features": ";".join(sorted(excluded))}
            row.update(_recovery(s, heldout, [int(k) for k in cfg["validation"]["top_k"]]))
            rows.append(row)
    out = pd.DataFrame(rows)
    write_tsv(out, out_tables / "framework_leave_source_out_validation.tsv")
    return out


def run_permutation_test(
    feature_matrix: pd.DataFrame,
    positive_sets: pd.DataFrame,
    benchmark_sets: dict[str, set[str]],
    scores: pd.DataFrame,
    matrix: pd.DataFrame,
    cfg: dict,
    out_tables: Path,
) -> pd.DataFrame:
    X_all = _safe_gene_index(feature_matrix)
    labels_idx = _safe_gene_index(positive_sets)
    scores_idx = _safe_gene_index(scores)
    bg_cfg = cfg["background_sampling"]
    strata = build_background_strata(matrix) if bool(bg_cfg.get("matched_background", True)) else None
    n_perm = int(cfg["validation"]["n_permutations"])
    repeats = int(cfg["validation"]["permutation_repeats"])
    seeds = cfg["random_seeds"][:repeats]
    rng = np.random.default_rng(17)
    rows = []
    for axis in AXES:
        X = X_all[axis_feature_columns(feature_matrix, cfg, axis.name)]
        y = pd.to_numeric(labels_idx[axis.label_column], errors="coerce").fillna(0).reindex(scores_idx.index).fillna(0).astype(int)
        benchmark = benchmark_sets[axis.name]
        all_benchmarks = set().union(*benchmark_sets.values())
        observed = _recovery(scores_idx[axis.score_column], benchmark, [50])["top50_recovery"]
        eligible = np.array(sorted(set(y.index) - all_benchmarks))
        vals = []
        fit_attempts = 0
        converged_fits = 0
        for perm in range(n_perm):
            perm_pos = set(rng.choice(eligible, size=int(y.sum()), replace=False))
            score, _, train_audit = fit_sampled_background_model(
                X, perm_pos, perm_pos | all_benchmarks, [int(seed) + perm * 1000 for seed in seeds], cfg, strata,
                bool(bg_cfg.get("matched_background", True)), int(bg_cfg.get("background_ratio", 8)), penalty="elasticnet",
            )
            fit_attempts += len(train_audit)
            converged_fits += int(train_audit["used_in_ensemble"].sum())
            vals.append(_recovery(score.set_index("gene_symbol")["score_mean"], benchmark, [50])["top50_recovery"])
        vals_arr = np.array(vals)
        rows.append({
            "axis": axis.name,
            "metric": "top50_recovery",
            "evaluation_set": "independent_benchmark",
            "n_development_positives": int(y.sum()),
            "n_benchmark_targets": len(benchmark),
            "n_model_repeats_per_permutation": len(seeds),
            "n_converged_fits": converged_fits,
            "n_fit_attempts": fit_attempts,
            "n_replacement_fits": fit_attempts - converged_fits,
            "fit_convergence_rate": converged_fits / max(fit_attempts, 1),
            "observed_value": observed,
            "permutation_mean": float(vals_arr.mean()),
            "permutation_sd": float(vals_arr.std()),
            "empirical_p_value": float((np.sum(vals_arr >= observed) + 1) / (len(vals_arr) + 1)),
        })
    out = pd.DataFrame(rows)
    write_tsv(out, out_tables / "framework_label_permutation_test.tsv")
    return out


def run_sensitivity(
    feature_matrix: pd.DataFrame,
    matrix: pd.DataFrame,
    positive_sets: pd.DataFrame,
    benchmark_sets: dict[str, set[str]],
    scores: pd.DataFrame,
    cfg: dict,
    out_tables: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    X_all = _safe_gene_index(feature_matrix)
    labels_idx = _safe_gene_index(positive_sets)
    scores_idx = _safe_gene_index(scores)
    positives = {
        "translational": set(labels_idx.index[pd.to_numeric(labels_idx["clinical_positive"], errors="coerce").fillna(0).astype(int).eq(1)]),
        "functional": set(labels_idx.index[pd.to_numeric(labels_idx["discovery_positive"], errors="coerce").fillna(0).astype(int).eq(1)]),
    }
    all_benchmarks = set().union(*benchmark_sets.values())
    sens_cfg = cfg["sensitivity"]
    seeds = cfg["random_seeds"][: int(sens_cfg["n_repeats"])]
    base_bg = cfg["background_sampling"]
    base_model = cfg["model"]
    strata = build_background_strata(matrix)
    variants = []
    for v in sens_cfg["l1_ratio_values"]:
        new_cfg = json.loads(json.dumps(cfg))
        new_cfg["model"]["l1_ratio"] = v
        variants.append((f"l1_ratio_{v}", "l1_ratio", new_cfg, int(base_bg["background_ratio"]), bool(base_bg["matched_background"])))
    for v in sens_cfg["C_values"]:
        new_cfg = json.loads(json.dumps(cfg))
        new_cfg["model"]["C"] = v
        variants.append((f"C_{v}", "C", new_cfg, int(base_bg["background_ratio"]), bool(base_bg["matched_background"])))
    for v in sens_cfg["background_ratio_values"]:
        variants.append((f"background_ratio_{v}", "background_ratio", cfg, int(v), bool(base_bg["matched_background"])))
    for v in sens_cfg["matched_background_values"]:
        variants.append((f"matched_background_{int(bool(v))}", "matched_background", cfg, int(base_bg["background_ratio"]), bool(v)))

    rows = []
    seen = set()
    for variant, parameter, variant_cfg, ratio, matched in variants:
        key = (variant, ratio, matched, variant_cfg["model"]["C"], variant_cfg["model"]["l1_ratio"])
        if key in seen:
            continue
        seen.add(key)
        for axis in AXES:
            pos_key = "translational" if axis.name == "translational" else "functional"
            X = X_all[axis_feature_columns(feature_matrix, cfg, axis.name)]
            score, _, _ = fit_sampled_background_model(
                X, positives[pos_key], positives[pos_key] | all_benchmarks, seeds, variant_cfg, strata if matched else None,
                matched, ratio, penalty="elasticnet",
            )
            s = score.set_index("gene_symbol")["score_mean"]
            primary = scores_idx[axis.score_column].reindex(s.index)
            top50 = _top_set(s, 50)
            primary_top50 = _top_set(primary, 50)
            rows.append({
                "variant": variant,
                "parameter": parameter,
                "axis": axis.name,
                "spearman_with_primary": float(s.corr(primary, method="spearman")),
                "top50_overlap_with_primary": len(top50 & primary_top50) / 50,
                "top50_benchmark_recovery": len(top50 & benchmark_sets[pos_key]) / max(len(benchmark_sets[pos_key]), 1),
                "n_repeats": len(seeds),
                "background_ratio": ratio,
                "matched_background": int(matched),
                "C": variant_cfg["model"].get("C", base_model["C"]),
                "l1_ratio": variant_cfg["model"].get("l1_ratio", base_model["l1_ratio"]),
            })
    model_sens = pd.DataFrame(rows)
    write_tsv(model_sens, out_tables / "framework_model_sensitivity.tsv")

    review_cfg = cfg["review"]
    threshold_rows = []
    for high_pct in sens_cfg["triage_high_score_percentile_values"]:
        for safety_min in sens_cfg["triage_safety_factor_min_values"]:
            temp_cfg = json.loads(json.dumps(cfg))
            temp_cfg["review"]["high_score_percentile"] = high_pct
            temp_cfg["review"]["safety_factor_min"] = safety_min
            temp_out = assign_review_classes(matrix, scores, temp_cfg, out_tables)
            counts = temp_out["review_class"].str.extract(r"^(Class [IVX]+)")[0].value_counts().to_dict()
            threshold_rows.append({
                "high_score_percentile": high_pct,
                "safety_factor_min": safety_min,
                "class_I_count": counts.get("Class I", 0),
                "class_II_count": counts.get("Class II", 0),
                "class_III_count": counts.get("Class III", 0),
                "class_IV_count": counts.get("Class IV", 0),
                "class_V_count": counts.get("Class V", 0),
                "class_VI_count": counts.get("Class VI", 0),
            })
    # Restore primary review-class file after temporary threshold sweeps.
    assign_review_classes(matrix, scores, cfg, out_tables)
    threshold_sens = pd.DataFrame(threshold_rows)
    write_tsv(threshold_sens, out_tables / "framework_review_threshold_sensitivity.tsv")
    return model_sens, threshold_sens


def build_evidence_cards(matrix: pd.DataFrame, scores: pd.DataFrame, review_classes: pd.DataFrame, out_tables: Path) -> pd.DataFrame:
    matrix_idx = _safe_gene_index(matrix)
    scores_idx = _safe_gene_index(scores)
    review_idx = _safe_gene_index(review_classes)
    genes = list(scores_idx.nsmallest(35, "translational_precedence_rank").index)
    genes += [g for g in scores_idx.nsmallest(35, "functional_discovery_rank").index if g not in genes]
    # Keep a small representative-candidate pool available even when one of the
    # automatically selected evidence-card examples falls outside the top-ranked
    # export window after model updates. The final manuscript examples are chosen
    # downstream from the generated evidence-card table rather than fixed here.
    representative_pool = ["KRAS", "ARID1A", "KIF23", "PKMYT1", "EGFR", "ERBB2", "PALB2", "BRCA1", "BRCA2", "ATM"]
    genes += [g for g in representative_pool if g in scores_idx.index and g not in genes]
    rows = []
    for gene in genes:
        row = matrix_idx.loc[gene]
        sc = scores_idx.loc[gene]
        rv = review_idx.loc[gene]
        limiting = []
        review_class = str(rv["review_class"])
        if float(row.get("dependency_score", 0)) < 0.05:
            if review_class.startswith("Class I:"):
                limiting.append("weak functional-discovery support")
            else:
                limiting.append("low dependency")
        if float(row.get("safety_factor", 1)) < 0.4:
            limiting.append("safety constraint")
        if float(row.get("external_reproducibility_score", 0)) == 0:
            limiting.append("limited external reproducibility")
        if float(row.get("druggability_score", 0)) < 0.15:
            limiting.append("limited actionability")
        rows.append({
            "gene_symbol": gene,
            "review_class": review_class,
            "translational_precedence_rank": int(sc["translational_precedence_rank"]),
            "functional_discovery_rank": int(sc["functional_discovery_rank"]),
            "translational_precedence_score": sc[CLINICAL_SCORE],
            "functional_discovery_score": sc[DISCOVERY_SCORE],
            "patient_alteration": row.get("patient_alteration_score", 0),
            "proteomic_support": row.get("proteomic_support_score", 0),
            "external_reproducibility": row.get("external_reproducibility_score", 0),
            "dependency": row.get("dependency_score", 0),
            "druggability": row.get("druggability_score", 0),
            "safety_factor": row.get("safety_factor", 1 - row.get("safety_risk_score", 0)),
            "safety_risk": row.get("safety_risk_score", 0),
            "open_targets": row.get("opentargets_overall_score", 0),
            "limiting_evidence": "; ".join(limiting),
        })
    out = pd.DataFrame(rows)
    write_tsv(out, out_tables / "framework_evidence_cards.tsv")
    return out


def write_results_summary(
    matrix: pd.DataFrame,
    positive_sets: pd.DataFrame,
    benchmark_sets: dict[str, set[str]],
    review_classes: pd.DataFrame,
    comparator: pd.DataFrame,
    leave_source: pd.DataFrame,
    independent_validation: pd.DataFrame,
    oof_summary: pd.DataFrame,
    permutation: pd.DataFrame,
    sensitivity: pd.DataFrame,
    threshold_sens: pd.DataFrame,
    evidence_cards: pd.DataFrame,
    out_tables: Path,
) -> pd.DataFrame:
    class_counts = review_classes["review_class"].str.extract(r"^(Class [IVX]+)")[0].value_counts().to_dict()
    labels = _safe_gene_index(positive_sets)
    summary_rows = [
        ("candidate_universe", len(matrix)),
        ("translational_positives", int(pd.to_numeric(labels["clinical_positive"], errors="coerce").fillna(0).sum())),
        ("functional_positives", int(pd.to_numeric(labels["discovery_positive"], errors="coerce").fillna(0).sum())),
        ("translational_benchmark_targets", len(benchmark_sets["translational"])),
        ("functional_benchmark_targets", len(benchmark_sets["functional"])),
        ("positive_union", len(labels)),
        ("class_I", class_counts.get("Class I", 0)),
        ("class_II", class_counts.get("Class II", 0)),
        ("class_III", class_counts.get("Class III", 0)),
        ("class_IV", class_counts.get("Class IV", 0)),
        ("class_V", class_counts.get("Class V", 0)),
        ("class_VI", class_counts.get("Class VI", 0)),
    ]
    for axis in ["translational", "functional"]:
        independent_row = independent_validation[independent_validation["axis"].eq(axis)].iloc[0]
        oof_row = oof_summary[oof_summary["axis"].eq(axis)].iloc[0]
        summary_rows.extend([
            (f"{axis}_independent_top50_recovery", independent_row["top50_recovery"]),
            (f"{axis}_independent_top50_ci_low", independent_row["top50_recovery_ci_low"]),
            (f"{axis}_independent_top50_ci_high", independent_row["top50_recovery_ci_high"]),
            (f"{axis}_independent_median_rank_percentile", independent_row["median_rank_percentile"]),
            (f"{axis}_oof_top50_recovery", oof_row["top50_recovery"]),
            (f"{axis}_oof_median_rank_percentile", oof_row["median_rank_percentile"]),
        ])
    for axis in ["translational", "functional"]:
        row = permutation[permutation["axis"].eq(axis)].iloc[0]
        summary_rows.extend([
            (f"{axis}_permutation_observed_top50", row["observed_value"]),
            (f"{axis}_permutation_mean_top50", row["permutation_mean"]),
            (f"{axis}_permutation_p", row["empirical_p_value"]),
        ])
    comp50 = comparator[comparator["top_k"].eq(50)]
    for _, row in comp50.iterrows():
        key = row["method"].lower().replace(" ", "_").replace("-", "_").replace("/", "_")
        summary_rows.append((f"{key}_{row['axis']}_recovery_at_50", row["positive_recovery"]))
        summary_rows.append((f"{key}_{row['axis']}_class_I_II_at_50", row["class_I_II_fraction"]))
    for axis in ["translational", "functional"]:
        sub = sensitivity[sensitivity["axis"].eq(axis)]
        summary_rows.append((f"{axis}_min_parameter_spearman", sub[sub["parameter"].ne("matched_background")]["spearman_with_primary"].min()))
        mb = sub[sub["parameter"].eq("matched_background")]
        if not mb.empty:
            summary_rows.append((f"{axis}_matched_to_unmatched_spearman", mb["spearman_with_primary"].min()))
    primary_threshold = threshold_sens[
        threshold_sens["high_score_percentile"].eq(0.02) &
        threshold_sens["safety_factor_min"].eq(0.4)
    ]
    if not primary_threshold.empty:
        for col in primary_threshold.columns:
            if col.endswith("_count"):
                summary_rows.append((f"primary_threshold_{col}", primary_threshold.iloc[0][col]))
    out = pd.DataFrame(summary_rows, columns=["metric", "value"])
    write_tsv(out, out_tables / "framework_results_summary.tsv")
    return out


def run_framework(project_root: Path, config_path: Path | None = None, output_root: Path | None = None) -> None:
    cfg = read_config(project_root, config_path)
    out_root = output_root or project_root / "results"
    out_tables = out_root / "tables"
    cancer_code = str(cfg.get("cancer_code", "PDAC")).upper()
    matrix, label_sources, source_positive_sets = copy_framework_inputs(cfg, out_tables, project_root)
    positive_sets, benchmark_sets, role_manifest, leakage_audit = build_validation_contract(label_sources, matrix, cfg)
    write_tsv(positive_sets, out_tables / "framework_development_positive_sets.tsv")
    benchmark_rows = role_manifest[role_manifest["record_role"].eq("V")][
        ["record_id", "axis", "gene_symbol", "source", "source_group", "label_type", "notes"]
    ].drop_duplicates()
    write_tsv(benchmark_rows, out_tables / "framework_validation_benchmarks.tsv")
    write_tsv(role_manifest, out_tables / "framework_record_role_manifest.tsv")
    write_tsv(leakage_audit, out_tables / "framework_leakage_audit.tsv")
    feature_matrix, manifest = build_feature_matrix(matrix, cfg)
    write_tsv(feature_matrix, out_tables / "framework_feature_matrix.tsv")
    write_tsv(manifest, out_tables / "framework_feature_manifest.tsv")
    write_tsv(build_source_role_audit(role_manifest, manifest), out_tables / "framework_source_role_audit.tsv")
    method_spec, review_spec, statistical_spec = build_method_specifications(cfg)
    write_tsv(method_spec, out_tables / "framework_sampling_and_model_specification.tsv")
    write_tsv(review_spec, out_tables / "framework_review_rule_specification.tsv")
    write_tsv(statistical_spec, out_tables / "framework_statistical_analysis_specification.tsv")
    scores = train_dual_axis_models(feature_matrix, positive_sets, benchmark_sets, matrix, cfg, out_tables)
    independent_validation = run_independent_benchmark_validation(scores, benchmark_sets, cfg, out_tables)
    _, oof_summary = run_out_of_fold_validation(
        feature_matrix, positive_sets, benchmark_sets, role_manifest, matrix, cfg, out_tables
    )
    review_classes = assign_review_classes(matrix, scores, cfg, out_tables)
    evidence_cards = build_evidence_cards(matrix, scores, review_classes, out_tables)
    comparator = run_comparators(feature_matrix, matrix, positive_sets, benchmark_sets, scores, review_classes, cfg, out_tables)
    leave_source = run_leave_source_out(feature_matrix, role_manifest, positive_sets, benchmark_sets, matrix, cfg, out_tables)
    permutation = run_permutation_test(feature_matrix, positive_sets, benchmark_sets, scores, matrix, cfg, out_tables)
    sensitivity, threshold_sens = run_sensitivity(feature_matrix, matrix, positive_sets, benchmark_sets, scores, cfg, out_tables)
    write_results_summary(
        matrix, positive_sets, benchmark_sets, review_classes, comparator, leave_source,
        independent_validation, oof_summary, permutation, sensitivity, threshold_sens,
        evidence_cards, out_tables,
    )

