import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_cloud_uri_is_localization_error(self):
        from io_utils import readable
        for path in ['gs://bucket/model.h5', 's3://bucket/peaks.bed', 'https://example.org/x']:
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, 'localization'):
                readable(path)

    def test_localized_path_with_shell_characters(self):
        from io_utils import readable
        path = self.root / "model ' $(touch NEVER) file.h5"
        path.write_text('model')
        self.assertEqual(readable(str(path)), path)
        self.assertFalse((self.root / 'NEVER').exists())

    def test_variant_schema_and_duplicate_ids(self):
        from io_utils import read_variants
        path = self.root / 'variants.tsv'
        path.write_text('chr1\t50\tA\tT\tv1\n')
        self.assertEqual(read_variants(path)[0]['pos'], '50')
        path.write_text('chr\tpos\tallele1\tallele2\tvariant_id\n')
        with self.assertRaisesRegex(ValueError, 'headerless'):
            read_variants(path)
        path.write_text('chr1\t50\tA\tT\tv1\nchr1\t51\tA\tC\tv1\n')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            read_variants(path)

    def test_reference_mismatch_and_indel_context(self):
        from io_utils import validate_reference
        genome = {'chr1': 'A' * 100}
        row = dict(chr='chr1', pos='50', allele1='A', allele2='T', variant_id='v1')
        validate_reference([row], genome, {'chr1': 100}, 20)
        with self.assertRaisesRegex(ValueError, 'REF'):
            validate_reference([dict(row, allele1='C')], genome, {'chr1': 100}, 20)
        with self.assertRaisesRegex(ValueError, 'window'):
            validate_reference([dict(row, pos='89', allele1='AAAA', allele2='-')], genome, {'chr1': 100}, 20)

    def test_manifest_rejects_empty_duplicate_or_unsafe_ids(self):
        from validate_manifest import validate_manifest
        row = dict(model_id='CD4_fold0', cell_type='CD4 T cell', model='gs://b/m.h5', peaks='gs://b/p.bed')
        self.assertEqual(validate_manifest([row]), 1)
        for rows in [[], [row, row], [dict(row, model_id='../bad')]]:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                validate_manifest(rows)

    def write_scores(self, name, model, cell, variants=None):
        path = self.root / name
        variants = variants or [('v1', '50', 'A', 'T'), ('v2', '60', 'A', 'C')]
        with path.open('w') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['model_id', 'cell_type', 'chr', 'pos', 'allele1', 'allele2', 'variant_id', 'logfc', 'jsd', 'active_allele_quantile'])
            for ident, pos, ref, alt in variants:
                writer.writerow([model, cell, 'chr1', pos, ref, alt, ident, '0.25', '0.1', '0.9'])
        return path

    def test_merge_keeps_models_and_variant_identity(self):
        from merge_scores import merge_scores
        one = self.write_scores('one.tsv', 'cd4', 'CD4')
        two = self.write_scores('two.tsv', 'nk', 'NK', [('v2', '60', 'A', 'C'), ('v1', '50', 'A', 'T')])
        merge_scores([one, two], self.root / 'long.tsv', self.root / 'wide.tsv')
        with (self.root / 'long.tsv').open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual(len(rows), 4)
        self.assertEqual({row['cell_type'] for row in rows}, {'CD4', 'NK'})
        with (self.root / 'wide.tsv').open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['cd4.logfc'], '0.25')
        self.assertEqual(rows[0]['nk.logfc'], '0.25')

    def test_merge_rejects_missing_variant_and_changed_alleles(self):
        from merge_scores import merge_scores
        one = self.write_scores('one.tsv', 'cd4', 'CD4')
        for variants in [[('v1', '50', 'A', 'T')], [('v1', '50', 'A', 'G'), ('v2', '60', 'A', 'C')]]:
            two = self.write_scores('two.tsv', 'nk', 'NK', variants)
            with self.subTest(variants=variants), self.assertRaises(ValueError):
                merge_scores([one, two], self.root / 'long.tsv', self.root / 'wide.tsv')

    def test_reference_index_staging_does_not_modify_localized_input(self):
        import os
        from score_variants import stage_reference
        genome = self.root / 'input.fa'
        index = self.root / 'input.fa.fai'
        index.write_text('original-index')
        os.utime(index, (1, 1))
        genome.write_text('>chr1\nAAAA\n')
        output = self.root / 'effects'
        output.mkdir()
        local_genome = stage_reference(genome, index, output)
        local_index = Path(str(local_genome) + '.fai')
        self.assertGreaterEqual(local_index.stat().st_mtime, genome.stat().st_mtime)
        local_index.write_text('updated-index')
        self.assertEqual(index.read_text(), 'original-index')
        self.assertEqual(local_genome.read_text(), genome.read_text())

    def test_scorer_ids_preserve_null_and_numeric_identifiers(self):
        from io_utils import write_scorer_variants, restore_score_file
        original = [dict(chr='chr1', pos='50', allele1='A', allele2='T', variant_id=ident)
                    for ident in ['null', 'NA', '001']]
        mapping = write_scorer_variants(original, self.root / 'scorer.tsv')
        tokens = list(mapping)
        self.assertEqual(len(set(tokens)), 3)
        scored = self.root / 'scored.tsv'
        with scored.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=['chr', 'pos', 'allele1', 'allele2', 'variant_id', 'logfc'], delimiter='\t')
            writer.writeheader()
            for token in tokens:
                writer.writerow(dict(original[0], variant_id=token, logfc='0.5'))
        restore_score_file(scored, mapping, require_all=True)
        with scored.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual([row['variant_id'] for row in rows], ['null', 'NA', '001'])
        self.assertEqual([row['logfc'] for row in rows], ['0.5'] * 3)
        with self.assertRaisesRegex(ValueError, 'identity'):
            restore_score_file(scored, mapping, require_all=True)

    def test_score_cli_rejects_cloud_model_before_tensorflow(self):
        result = subprocess.run([sys.executable, str(ROOT / 'scripts/score_variants.py'),
            '--variants', 'gs://b/v.tsv', '--model', 'gs://b/m.h5', '--peaks', 'gs://b/p.bed',
            '--genome', 'gs://b/g.fa', '--genome-index', 'gs://b/g.fa.fai', '--chrom-sizes', 'gs://b/s.tsv',
            '--model-id', 'm1', '--cell-type', 'CD4', '--output-dir', str(self.root)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('localization', result.stderr)


if __name__ == '__main__':
    unittest.main()
