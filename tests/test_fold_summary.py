import csv
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))


class FoldSummaryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.input = self.root / "scores ' $(touch NEVER).tsv"
        self.output = self.root / 'summary.tsv'
        self.fields = ['model_id', 'cell_type', 'chr', 'pos', 'allele1', 'allele2',
                       'variant_id', 'logfc', 'abs_logfc', 'jsd', 'active_allele_quantile',
                       'abs_logfc_x_active_allele_quantile', 'logfc.pval']
        self.rows = []
        for model, cell, logfc, quantile in [('a_fold0', 'CD4', 1, 0.2),
                                           ('a_fold1', 'CD4', -0.5, 0.8),
                                           ('b_fold0', 'CD4', 2, 0.9),
                                           ('nk', 'NK', 0, 0.5)]:
            for ident, pos in [('001', '50'), ('NA', '60')]:
                self.rows.append(dict(model_id=model, cell_type=cell, chr='chr1', pos=pos,
                    allele1='A', allele2='T', variant_id=ident, logfc=str(logfc),
                    abs_logfc=str(abs(logfc)), jsd='0.1', active_allele_quantile=str(quantile),
                    abs_logfc_x_active_allele_quantile=str(abs(logfc) * quantile),
                    **{'logfc.pval': '0.01'}))

    def write_input(self, rows=None):
        with self.input.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=self.fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(self.rows if rows is None else rows)

    def summarize(self, rows=None):
        from summarize_folds import summarize_folds
        self.write_input(rows)
        summarize_folds(self.input, self.output)
        with self.output.open() as stream:
            return list(csv.DictReader(stream, delimiter='\t'))

    def test_means_products_dispersion_and_direction_without_mixing_models(self):
        rows = self.summarize(list(reversed(self.rows)))
        self.assertEqual(len(rows), 6)
        self.assertEqual({r['model_group'] for r in rows}, {'a', 'b', 'nk'})
        a = next(r for r in rows if r['model_group'] == 'a' and r['variant_id'] == '001')
        self.assertEqual(a['n_folds'], '2')
        self.assertEqual(a['model_ids'], 'a_fold0,a_fold1')
        self.assertAlmostEqual(float(a['logfc.mean']), 0.25)
        self.assertAlmostEqual(float(a['abs_logfc.mean']), 0.75)
        # Mean of products is 0.3, not product of means (0.75 * 0.5 = 0.375).
        self.assertAlmostEqual(float(a['abs_logfc_x_active_allele_quantile.mean']), 0.3)
        self.assertAlmostEqual(float(a['logfc.sd']), math.sqrt(1.125))
        self.assertAlmostEqual(float(a['percent_change']), 18.920711500272102)
        self.assertEqual((a['n_positive'], a['n_negative'], a['n_zero']), ('1', '1', '0'))
        self.assertEqual(a['direction_agreement'], '0.5')
        self.assertFalse(any('pval' in field for field in a))

    def test_single_model_has_no_sample_sd_and_zero_direction_is_reported(self):
        rows = self.summarize()
        nk = next(r for r in rows if r['model_group'] == 'nk')
        self.assertEqual(nk['n_folds'], '1')
        self.assertEqual(nk['logfc.sd'], '')
        self.assertEqual(nk['n_zero'], '1')
        self.assertEqual(nk['direction_agreement'], '0.0')
        self.assertEqual(nk['variant_id'], '001')

    def test_cli_reads_local_input_and_writes_summary(self):
        self.write_input()
        result = subprocess.run([sys.executable, str(ROOT / 'scripts/summarize_folds.py'),
            '--scores', str(self.input), '--output', str(self.output)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(self.output.is_file())
        self.assertFalse((self.root / 'NEVER').exists())

    def test_rejects_invalid_or_incomplete_fold_tables(self):
        from summarize_folds import summarize_folds
        cases = [
            ('duplicate', self.rows + [self.rows[0]]),
            ('alleles', [dict(r, allele2='G') if r['model_id'] == 'a_fold1' else r for r in self.rows]),
            ('Variant set', self.rows[1:]),
            ('cell type', [dict(r, cell_type='NK') if r['model_id'] == 'a_fold1' else r for r in self.rows]),
            ('Nonfinite', [dict(r, abs_logfc='nan') if r is self.rows[0] else r for r in self.rows]),
            ('fold index', [dict(r, model_id='a_fold00') if r['model_id'] == 'a_fold1' else r for r in self.rows]),
            ('unsuffixed', [dict(r, model_id='a') if r['model_id'] == 'a_fold1' else r for r in self.rows]),
        ]
        for message, rows in cases:
            with self.subTest(message=message):
                self.write_input(rows)
                with self.assertRaisesRegex(ValueError, message):
                    summarize_folds(self.input, self.output)

    def test_cloud_uri_is_a_localization_error(self):
        from summarize_folds import summarize_folds
        with self.assertRaisesRegex(ValueError, 'localization'):
            summarize_folds('gs://bucket/merged.tsv', self.output)


if __name__ == '__main__':
    unittest.main()
