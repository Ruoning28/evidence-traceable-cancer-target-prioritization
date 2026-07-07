from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from framework.pipeline import (
    AXES,
    axis_feature_columns,
    build_feature_matrix,
    build_validation_contract,
    copy_framework_inputs,
    read_config,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODULE_CONFIGS = {
    "PDAC": PROJECT_ROOT / "PDAC" / "config" / "framework_config.json",
    "LUAD": PROJECT_ROOT / "LUDA" / "config" / "framework_config.json",
}


class MethodContractTests(unittest.TestCase):
    """Check the predeclared split, embargo and feature-control contract."""

    def _build_contract(self, cancer: str):
        """Construct one cancer's role and feature tables from bundled inputs.

        Parameters
        ----------
        cancer:
            ``PDAC`` or ``LUAD``.

        Returns
        -------
        tuple
            Configuration, development labels, benchmark sets, role manifest,
            leakage audit, feature matrix and feature manifest.
        """
        cfg = read_config(PROJECT_ROOT, MODULE_CONFIGS[cancer])
        with tempfile.TemporaryDirectory() as directory:
            matrix, label_sources, _ = copy_framework_inputs(cfg, Path(directory), PROJECT_ROOT)
        development, benchmarks, roles, audit = build_validation_contract(label_sources, matrix, cfg)
        features, feature_manifest = build_feature_matrix(matrix, cfg)
        return cfg, development, benchmarks, roles, audit, features, feature_manifest

    def test_development_and_benchmark_genes_are_disjoint(self) -> None:
        """Require zero L/V gene overlap for both axes and cancers."""
        for cancer in MODULE_CONFIGS:
            with self.subTest(cancer=cancer):
                _, development, benchmarks, roles, audit, _, _ = self._build_contract(cancer)
                self.assertTrue(audit["status"].eq("PASS").all())
                self.assertTrue(audit["gene_overlap_count"].eq(0).all())
                for axis in AXES:
                    dev_genes = set(
                        development.loc[development[axis.label_column].eq(1), "gene_symbol"].astype(str)
                    )
                    self.assertFalse(dev_genes & benchmarks[axis.name])
                    benchmark_records = roles[
                        roles["axis"].eq(axis.name) & roles["record_role"].eq("V")
                    ]
                    self.assertEqual(set(benchmark_records["gene_symbol"]), benchmarks[axis.name])

    def test_axis_feature_exclusions_are_enforced(self) -> None:
        """Require every configured leakage feature and indicator to be absent."""
        for cancer in MODULE_CONFIGS:
            with self.subTest(cancer=cancer):
                cfg, _, _, _, _, features, _ = self._build_contract(cancer)
                for axis in AXES:
                    selected = set(axis_feature_columns(features, cfg, axis.name))
                    excluded = set(cfg["evaluation"]["axes"][axis.name].get("feature_exclusions", []))
                    for feature in excluded:
                        self.assertNotIn(feature, selected)
                        self.assertNotIn(f"{feature}_missing", selected)

    def test_formal_resampling_settings_are_symmetric(self) -> None:
        """Fix formal permutation count and observed/null repeat symmetry."""
        for cancer, config_path in MODULE_CONFIGS.items():
            with self.subTest(cancer=cancer):
                cfg = read_config(PROJECT_ROOT, config_path)
                self.assertEqual(cfg["validation"]["n_permutations"], 500)
                self.assertEqual(
                    cfg["validation"]["permutation_repeats"],
                    cfg["background_sampling"]["n_repeats"],
                )
                self.assertEqual(cfg["validation"]["leave_source_out_repeats"], 30)


if __name__ == "__main__":
    unittest.main()

