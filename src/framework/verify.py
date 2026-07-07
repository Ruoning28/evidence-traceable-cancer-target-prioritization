from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

EXPECTED = {
    "PDAC": {
        "module_path": PROJECT_ROOT / "PDAC",
        "candidate_universe": 40591,
        "translational_positives": 7,
        "functional_positives": 21,
        "class_I": 298,
        "class_II": 75,
        "class_III": 1884,
        "class_IV": 1184,
        "class_V": 2967,
        "class_VI": 34183,
        "translational_observed_top50": 0.363636,
        "translational_permutation_mean": 0.015273,
        "translational_empirical_p": 0.001996,
        "functional_observed_top50": 0.133333,
        "functional_permutation_mean": 0.014800,
        "functional_empirical_p": 0.005988,
    },
    "LUAD": {
        "module_path": PROJECT_ROOT / "LUDA",
        "candidate_universe": 35478,
        "translational_positives": 6,
        "functional_positives": 26,
        "class_I": 266,
        "class_II": 14,
        "class_III": 574,
        "class_IV": 1184,
        "class_V": 1392,
        "class_VI": 32048,
        "translational_observed_top50": 0.571429,
        "translational_permutation_mean": 0.015429,
        "translational_empirical_p": 0.001996,
        "functional_observed_top50": 0.055556,
        "functional_permutation_mean": 0.008167,
        "functional_empirical_p": 0.039920,
    },
}


def read_summary(module_path: Path) -> pd.Series:
    """Read the framework summary table for one module.

    Parameters
    ----------
    module_path:
        Path to ``PDAC`` or ``LUDA`` package directory.

    Returns
    -------
    pandas.Series
        Metric-to-value mapping from ``framework_results_summary.tsv``.
    """
    path = module_path / "result" / "tables" / "framework_results_summary.tsv"
    return pd.read_csv(path, sep="\t").set_index("metric")["value"].astype(float)


def read_permutation(module_path: Path) -> pd.DataFrame:
    """Read the label-permutation validation table.

    Parameters
    ----------
    module_path:
        Path to ``PDAC`` or ``LUDA`` package directory.

    Returns
    -------
    pandas.DataFrame
        Rows keyed by axis with observed recovery, permutation mean and
        empirical p-value columns.
    """
    path = module_path / "result" / "tables" / "framework_label_permutation_test.tsv"
    return pd.read_csv(path, sep="\t")


def close_enough(observed: float, expected: float, tolerance: float) -> bool:
    """Compare two numbers using an absolute tolerance."""
    return abs(float(observed) - float(expected)) <= tolerance


def verify_results() -> pd.DataFrame:
    """Compare result tables against the numbers reported in the manuscript.

    Parameters
    ----------
    None

    Returns
    -------
    pandas.DataFrame
        One row per checked value with PASS/FAIL status.
    """
    rows: list[dict[str, object]] = []
    for cancer, expected in EXPECTED.items():
        module_path = expected["module_path"]
        summary = read_summary(module_path)
        permutation = read_permutation(module_path).set_index("axis")

        for metric, exp_value in expected.items():
            if metric == "module_path":
                continue
            axis_metric_map = {
                "translational_observed_top50": ("translational", "observed_value"),
                "translational_permutation_mean": ("translational", "permutation_mean"),
                "translational_empirical_p": ("translational", "empirical_p_value"),
                "functional_observed_top50": ("functional", "observed_value"),
                "functional_permutation_mean": ("functional", "permutation_mean"),
                "functional_empirical_p": ("functional", "empirical_p_value"),
            }
            if metric in axis_metric_map:
                axis, column = axis_metric_map[metric]
                row = permutation.loc[axis]
                observed = float(row[column])
                tolerance = 5e-4
            else:
                observed = float(summary[metric])
                tolerance = 0.0
            rows.append({
                "cancer": cancer,
                "metric": metric,
                "observed": observed,
                "expected": exp_value,
                "tolerance": tolerance,
                "status": "PASS" if close_enough(observed, exp_value, tolerance) else "FAIL",
            })

        class_total = sum(float(summary[f"class_{roman}"]) for roman in ["I", "II", "III", "IV", "V", "VI"])
        universe = float(summary["candidate_universe"])
        rows.append({
            "cancer": cancer,
            "metric": "class_total_equals_universe",
            "observed": class_total,
            "expected": universe,
            "tolerance": 0.0,
            "status": "PASS" if class_total == universe else "FAIL",
        })
    return pd.DataFrame(rows)


def write_verification_report() -> Path:
    """Write the verification report to ``result/logs``.

    Parameters
    ----------
    None

    Returns
    -------
    Path
        Path to the generated TSV report.
    """
    report = verify_results()
    out_path = PROJECT_ROOT / "result" / "logs" / "result_verification.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out_path, sep="\t", index=False)
    return out_path

