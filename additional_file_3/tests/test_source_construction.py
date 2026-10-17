"""Small-file checks for source-table saving and Excel ingestion."""
from contextlib import ExitStack
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scientific_source/reference'))
import v4_end_to_end as source
import build_source_matrix as builder


class SourceConstructionTests(unittest.TestCase):
    def test_save_gene_index_and_values(self):
        frame = pd.DataFrame({'support': [0.25, 0.75]}, index=['KRAS', 'EGFR'])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'nested/evidence.tsv'
            source._save_indexed_table(frame, path)
            restored = pd.read_csv(path, sep='\t').set_index('gene_symbol')
        pd.testing.assert_frame_equal(restored, frame.rename_axis('gene_symbol'))

    def test_pdac_entry_point_reaches_final_save(self):
        patient = pd.DataFrame({'patient_alteration_score': [0.2]}, index=['KRAS'])
        geo = pd.DataFrame({'geo_external_reproducibility_score': [0.5]}, index=['KRAS'])
        empty = pd.DataFrame()
        with tempfile.TemporaryDirectory() as tmp, ExitStack() as stack:
            for name, value in {
                'build_patient_and_protein_features': patient,
                'build_tcga_features': empty,
                'build_geo_validation_features': (geo, empty),
                'build_dependency_features': empty,
                'build_druggability_features': (empty, empty),
                'build_safety_features': empty,
                'load_known_targets': {},
                'extract_literature_benchmark': empty,
            }.items():
                stack.enter_context(patch.object(source, name, return_value=value))
            output = Path(tmp) / 'output'
            builder.pdac_matrix(Path(tmp), output)
            saved = pd.read_csv(output / 'pdac_target_evidence_matrix.tsv', sep='\t')
        self.assertEqual(saved.gene_symbol.tolist(), ['KRAS'])
        self.assertEqual(saved.external_reproducibility_score.iloc[0], 0.5)
        self.assertIn('clinical_actionability_score', saved)

    def make_ttd(self, root):
        folder = root / 'data/raw/ttd'
        folder.mkdir(parents=True)
        (folder / 'P1-01-TTD_target_download.txt').write_text(
            'x\tTARGETID\tT1\nx\tGENENAME\tEGFR\n', encoding='utf-8')
        (folder / 'P1-02-TTD_drug_download.txt').write_text(
            'x\tDRUG__ID\tD1\nx\tTRADNAME\tExample drug\n', encoding='utf-8')
        workbook = folder / 'P1-07-Drug-TargetMapping.xlsx'
        pd.DataFrame({'TargetID': ['T1'], 'DrugID': ['D1'],
                      'Highest_status': ['Approved']}).to_excel(
                          workbook, index=False, engine='openpyxl')
        return workbook

    def test_ttd_workbook_contributes_features_and_edges(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_ttd(root)
            features, edges = source.build_druggability_features(root)
        self.assertEqual(features.loc['EGFR', 'ttd_approved_mapping_count'], 1)
        self.assertEqual(features.loc['EGFR', 'ttd_clinical_mapping_count'], 1)
        self.assertEqual(edges.gene_symbol.tolist(), ['EGFR'])
        self.assertEqual(edges.source.tolist(), ['TTD'])

    def test_missing_excel_engine_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_ttd(root)
            with patch.object(source.pd, 'read_excel', side_effect=ImportError('openpyxl unavailable')):
                with self.assertRaisesRegex(RuntimeError, 'TTD mapping workbook') as caught:
                    source.build_druggability_features(root)
            self.assertIsInstance(caught.exception.__cause__, ImportError)

    def test_invalid_ttd_workbook_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workbook = self.make_ttd(root)
            workbook.write_bytes(b'invalid workbook')
            with self.assertRaisesRegex(RuntimeError, 'TTD mapping workbook'):
                source.build_druggability_features(root)

    def test_literature_excel_is_read_and_failures_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            folder = root / 'data/raw/literature_benchmark/nwosu2025'
            folder.mkdir(parents=True)
            workbook = folder / 'targets.xlsx'
            pd.DataFrame([['EGFR']]).to_excel(workbook, header=False, index=False, engine='openpyxl')
            result = source.extract_literature_benchmark(root, {'EGFR'})
            self.assertIn('EGFR', result.gene_symbol.tolist())
            workbook.write_bytes(b'invalid workbook')
            with self.assertRaisesRegex(RuntimeError, 'literature workbook'):
                source.extract_literature_benchmark(root, {'EGFR'})


if __name__ == '__main__':
    unittest.main(verbosity=2)
