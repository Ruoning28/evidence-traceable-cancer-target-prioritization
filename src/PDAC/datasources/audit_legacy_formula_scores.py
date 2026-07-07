"""Audit the legacy PDAC fixed-formula scores included in the submission bundle.

The PDAC module ships an analysis-ready evidence matrix migrated from the
earlier ``pdac_target_prioritization`` project.  That matrix already contains
the legacy fixed-formula comparator columns used by the manuscript:
``clinical_actionability_score`` and ``discovery_potential_score``.  This
script recomputes those columns from their documented component formulas so
that reviewers can verify their provenance without needing the older project.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = [
    "patient_alteration_score",
    "proteomic_support_score",
    "external_reproducibility_score",
    "dependency_score",
    "druggability_score",
    "therapeutic_window_score",
    "target_tractability_score",
    "safety_risk_score",
]

AUDITED_COLUMNS = [
    "disease_relevance_score",
    "safety_factor",
    "safety_modifier",
    "clinical_actionability_core",
    "clinical_actionability_score",
    "discovery_potential_core",
    "discovery_potential_score",
]


def repository_root() -> Path:
    """Return the submission repository root.

    Parameters
    ----------
    None

    Returns
    -------
    pathlib.Path
        Absolute path to the ``evidence_traceable_submission`` directory.
    """
    return Path(__file__).resolve().parents[2]


def bounded_numeric(df: pd.DataFrame, column: str) -> pd.Series:
    """Read a score column as a bounded numeric series.

    Parameters
    ----------
    df : pandas.DataFrame
        Evidence matrix containing the requested column.
    column : str
        Column name to transform.

    Returns
    -------
    pandas.Series
        Numeric values with non-numeric or missing entries set to 0 and values
        clipped to the 0--1 score interval.
    """
    if column not in df.columns:
        raise KeyError(f"Required PDAC formula column is missing: {column}")
    return pd.to_numeric(df[column], errors="coerce").fillna(0.0).clip(0.0, 1.0)


def mean_frame(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    """Compute the row-wise arithmetic mean used by the legacy PDAC code.

    Parameters
    ----------
    frame : pandas.DataFrame
        Bounded component-score matrix.
    columns : list[str]
        Component columns included in the mean.

    Returns
    -------
    pandas.Series
        Row-wise arithmetic mean clipped to the 0--1 interval.
    """
    return frame[columns].fillna(0.0).clip(0.0, 1.0).mean(axis=1).clip(0.0, 1.0)


def geomean_frame(frame: pd.DataFrame, columns: list[str], epsilon: float = 1e-6) -> pd.Series:
    """Compute the row-wise geometric mean used by the legacy PDAC code.

    Parameters
    ----------
    frame : pandas.DataFrame
        Bounded component-score matrix.
    columns : list[str]
        Component columns included in the geometric mean.
    epsilon : float, default=1e-6
        Lower bound used before taking logarithms, matching the legacy code.

    Returns
    -------
    pandas.Series
        Row-wise geometric mean of the selected score columns.
    """
    arr = frame[columns].fillna(0.0).clip(epsilon, 1.0).to_numpy(float)
    return pd.Series(np.exp(np.mean(np.log(arr), axis=1)), index=frame.index)


def recompute_pdac_formula_scores(matrix: pd.DataFrame) -> pd.DataFrame:
    """Recompute the legacy PDAC fixed-formula score columns.

    Parameters
    ----------
    matrix : pandas.DataFrame
        PDAC evidence matrix containing the component-score columns listed in
        ``REQUIRED_COLUMNS``.

    Returns
    -------
    pandas.DataFrame
        Data frame indexed like ``matrix`` with recomputed disease relevance,
        safety, clinical-actionability, and discovery-potential formula scores.
    """
    components = pd.DataFrame({col: bounded_numeric(matrix, col) for col in REQUIRED_COLUMNS}, index=matrix.index)

    disease_components = pd.DataFrame(
        {
            "patient_alteration_score": components["patient_alteration_score"],
            "proteomic_support_score": components["proteomic_support_score"],
            "external_reproducibility_score": components["external_reproducibility_score"],
        },
        index=matrix.index,
    )
    disease_relevance = mean_frame(disease_components, list(disease_components.columns))

    safety_factor = (1.0 - components["safety_risk_score"]).clip(0.0, 1.0)

    clinical_components = pd.DataFrame(
        {
            "disease_relevance_score": disease_relevance,
            "dependency_score": components["dependency_score"],
            "druggability_score": components["druggability_score"],
        },
        index=matrix.index,
    )
    clinical_core = geomean_frame(clinical_components, list(clinical_components.columns))

    discovery_components = pd.DataFrame(
        {
            "disease_relevance_score": disease_relevance,
            "dependency_score": components["dependency_score"],
            "therapeutic_window_score": components["therapeutic_window_score"],
        },
        index=matrix.index,
    )
    discovery_core = geomean_frame(discovery_components, list(discovery_components.columns))

    return pd.DataFrame(
        {
            "disease_relevance_score": disease_relevance,
            "safety_factor": safety_factor,
            "safety_modifier": safety_factor,
            "clinical_actionability_core": clinical_core,
            "clinical_actionability_score": clinical_core.mul(safety_factor).clip(0.0, 1.0),
            "discovery_potential_core": discovery_core,
            "discovery_potential_score": discovery_core.mul(np.sqrt(components["target_tractability_score"])).clip(0.0, 1.0),
        },
        index=matrix.index,
    )


def audit_formula_scores(matrix: pd.DataFrame, tolerance: float = 1e-10) -> pd.DataFrame:
    """Compare shipped PDAC formula columns with recomputed values.

    Parameters
    ----------
    matrix : pandas.DataFrame
        PDAC evidence matrix containing shipped formula columns.
    tolerance : float, default=1e-10
        Maximum absolute difference allowed for a column to pass.

    Returns
    -------
    pandas.DataFrame
        One audit row per formula column, including maximum absolute difference
        and pass/fail status.
    """
    recomputed = recompute_pdac_formula_scores(matrix)
    rows = []
    for column in AUDITED_COLUMNS:
        if column not in matrix.columns:
            rows.append({"column": column, "status": "missing", "max_abs_difference": np.nan})
            continue
        observed = pd.to_numeric(matrix[column], errors="coerce").fillna(0.0)
        diff = (observed - recomputed[column]).abs()
        rows.append(
            {
                "column": column,
                "status": "pass" if float(diff.max()) <= tolerance else "fail",
                "max_abs_difference": float(diff.max()),
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the audit script.

    Parameters
    ----------
    None

    Returns
    -------
    argparse.Namespace
        Parsed command-line options.
    """
    root = repository_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix",
        type=Path,
        default=root / "PDAC" / "data" / "processed" / "pdac_target_evidence_matrix.tsv",
        help="PDAC evidence matrix to audit.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "PDAC" / "result" / "tables" / "pdac_legacy_formula_score_audit.tsv",
        help="Output TSV path for the audit summary.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-10,
        help="Maximum absolute difference allowed for a pass.",
    )
    return parser.parse_args()


def main() -> None:
    """Run the PDAC legacy fixed-formula score audit.

    Parameters
    ----------
    None

    Returns
    -------
    None
        Writes an audit TSV and raises ``SystemExit`` if any audited column
        fails or is missing.
    """
    args = parse_args()
    matrix = pd.read_csv(args.matrix, sep="\t")
    audit = audit_formula_scores(matrix, tolerance=args.tolerance)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output, sep="\t", index=False)
    print(audit.to_string(index=False))
    if not audit["status"].eq("pass").all():
        raise SystemExit(f"PDAC legacy formula score audit failed; see {args.output}")


if __name__ == "__main__":
    main()

