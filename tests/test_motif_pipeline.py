import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'motif_pipeline/scripts'))


class MotifPipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def fixture(self, fold, score=None, raw=None, region_text=None, head='counts'):
        path = self.root / f'fold{fold}.h5'
        raw = np.eye(4, dtype=np.int8)[np.arange(32) % 4].T[None] if raw is None else raw
        score = np.full(raw.shape, fold - 2, dtype=np.float32) if score is None else score
        with h5py.File(path, 'w') as handle:
            handle.create_dataset('raw/seq', data=raw)
            handle.create_dataset('shap/seq', data=score)
            handle.create_dataset('projected_shap/seq', data=score * raw)
            handle.attrs.update(head=head, cell_type='CD4', fold_index=fold)
        bed = self.root / f'fold{fold}.bed'
        bed.write_text(region_text or 'chr1\t40\t60\tp1\t0\t.\t0\t0\t0\t10\n')
        return path, bed

    def test_average_preserves_signed_hypothetical_scores_and_sequence(self):
        from average_contributions import average_contributions
        # Fold scores -2,-1,0,1,2 must cancel, not become mean(|score|)=1.2.
        files, beds = zip(*(self.fixture(i) for i in range(5)))
        out = self.root / 'average.h5'
        result = average_contributions(files, beds, out, self.root / 'regions.bed',
                                       expected_folds=5, chunk_rows=1)
        with h5py.File(out) as handle, h5py.File(files[0]) as first:
            np.testing.assert_array_equal(handle['shap/seq'][:], 0)
            np.testing.assert_array_equal(handle['projected_shap/seq'][:], 0)
            np.testing.assert_array_equal(handle['raw/seq'][:], first['raw/seq'][:])
        self.assertEqual(result['fold_count'], 5)
        self.assertEqual((self.root / 'regions.bed').read_text(), beds[0].read_text())

    def test_average_is_arithmetic_mean_with_float32_output(self):
        from average_contributions import average_contributions
        files, beds = zip(*(self.fixture(i, score=np.full((1, 4, 32), i / 10, dtype=np.float32)) for i in range(5)))
        out = self.root / 'average.h5'
        average_contributions(files, beds, out, self.root / 'regions.bed', chunk_rows=1)
        with h5py.File(out) as handle:
            np.testing.assert_allclose(handle['shap/seq'][:], 0.2)
            self.assertEqual(handle['shap/seq'].dtype, np.dtype('float32'))

    def test_average_rejects_changed_sequence_order(self):
        from average_contributions import average_contributions
        fixtures = [self.fixture(i) for i in range(4)]
        changed = np.eye(4, dtype=np.int8)[(np.arange(32) + 1) % 4].T[None]
        fixtures.append(self.fixture(4, raw=changed))
        files, beds = zip(*fixtures)
        with self.assertRaisesRegex(ValueError, 'sequence'):
            average_contributions(files, beds, self.root / 'average.h5', self.root / 'regions.bed')

    def test_average_rejects_changed_regions(self):
        from average_contributions import average_contributions
        fixtures = [self.fixture(i) for i in range(4)]
        fixtures.append(self.fixture(4, region_text='chr1\t41\t61\tp1\t0\t.\t0\t0\t0\t10\n'))
        files, beds = zip(*fixtures)
        with self.assertRaisesRegex(ValueError, 'region'):
            average_contributions(files, beds, self.root / 'average.h5', self.root / 'regions.bed')

    def test_average_rejects_missing_duplicate_or_mixed_folds(self):
        from average_contributions import average_contributions
        files, beds = zip(*(self.fixture(i) for i in range(5)))
        for fs, bs in [(files[:4], beds[:4]), (files[:4] + (files[0],), beds)]:
            with self.subTest(files=fs), self.assertRaisesRegex(ValueError, 'fold'):
                average_contributions(fs, bs, self.root / 'average.h5', self.root / 'regions.bed')
        with h5py.File(files[-1], 'a') as handle:
            handle.attrs['head'] = 'profile'
        with self.assertRaisesRegex(ValueError, 'head'):
            average_contributions(files, beds, self.root / 'average.h5', self.root / 'regions.bed')

    def test_average_rejects_nonfinite_or_invalid_onehot(self):
        from average_contributions import average_contributions
        files, beds = zip(*(self.fixture(i) for i in range(5)))
        with h5py.File(files[-1], 'a') as handle:
            handle['shap/seq'][0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, 'finite'):
            average_contributions(files, beds, self.root / 'average.h5', self.root / 'regions.bed')
        self.fixture(4)
        with h5py.File(files[-1], 'a') as handle:
            handle['raw/seq'][0, :, 0] = 1
        with self.assertRaisesRegex(ValueError, 'one-hot'):
            average_contributions(files, beds, self.root / 'average.h5', self.root / 'regions.bed')

    def test_peak_preparation_retains_all_valid_peaks_and_records_edge_exclusions(self):
        from prepare_motif_inputs import prepare_peaks
        bed = self.root / 'input.bed'
        bed.write_text('chr1\t0\t8\tedge\nchr1\t20\t40\ta\nchr1\t50\t70\tb\n')
        out = self.root / 'peaks.bed'
        meta = prepare_peaks(bed, {'chr1': 100}, 32, out, random_seed=7)
        self.assertEqual(meta['input_peaks'], 3)
        self.assertEqual(meta['retained_peaks'], 2)
        self.assertEqual(meta['excluded_edge_peaks'], 1)
        self.assertEqual({line.split('\t')[3] for line in out.read_text().splitlines()}, {'a', 'b'})
        self.assertTrue(all(len(line.split('\t')) == 10 for line in out.read_text().splitlines()))
        other = self.root / 'same.bed'
        prepare_peaks(bed, {'chr1': 100}, 32, other, random_seed=7)
        self.assertEqual(out.read_bytes(), other.read_bytes())

    def test_bad_summit_unknown_chromosome_and_duplicate_center_fail(self):
        from prepare_motif_inputs import prepare_peaks
        bad_rows = ['chr1\t20\t40\tp\t0\t.\t0\t0\t0\t-1\n',
                    'chr2\t20\t40\n', 'chr1\t20\t40\nchr1\t25\t35\n']
        for row in bad_rows:
            path = self.root / 'bad.bed'
            path.write_text(row)
            with self.subTest(row=row), self.assertRaises(ValueError):
                prepare_peaks(path, {'chr1': 100}, 32, self.root / 'out.bed')

    def test_discovery_explains_empty_threshold_distribution(self):
        from unittest.mock import patch
        from discover_motifs import main
        raw = np.eye(4, dtype=np.int8)[np.arange(64) % 4].T[None]
        source, _ = self.fixture(0, raw=raw)
        argv = ['discover_motifs.py', '--contributions', str(source), '--window', '40',
                '--output', str(self.root / 'motifs.h5'), '--metadata', str(self.root / 'metadata.json')]
        with patch.object(sys, 'argv', argv), patch('discover_motifs.shutil.which', return_value='/fake/modisco'), \
                patch('discover_motifs.runpy.run_path', side_effect=ValueError('Found array with 0 sample(s)')):
            with self.assertRaisesRegex(ValueError, 'empty threshold-fitting distribution'):
                main()

    def test_cloud_paths_fail_as_localization_errors(self):
        from motif_io import readable, read_file_list
        for value in ['gs://bucket/model.h5', 's3://bucket/scores.h5', 'https://example.org/x']:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'localization'):
                readable(value)
        file_list = self.root / 'files.list'
        file_list.write_text('gs://bucket/model.h5\n')
        with self.assertRaisesRegex(ValueError, 'localization'):
            read_file_list(file_list)

    def test_model_width_check_and_five_fold_validation(self):
        from prepare_motif_inputs import validate_models
        models = []
        for fold in range(5):
            path = self.root / f'model{fold}.h5'
            with h5py.File(path, 'w') as handle:
                handle.attrs['model_config'] = json.dumps({'config': {'layers': [
                    {'class_name': 'InputLayer', 'config': {'batch_input_shape': [None, 32, 4]}}]}})
                handle.create_dataset('test_weights', data=[fold])
            models.append(path)
        validate_models(models, 32, expected_folds=5)
        with self.assertRaisesRegex(ValueError, 'length'):
            validate_models(models, 64, expected_folds=5)
        with self.assertRaisesRegex(ValueError, 'five|5|fold'):
            validate_models(models[:4], 32, expected_folds=5)

    def test_reference_staging_handles_gzip_and_does_not_modify_supplied_index(self):
        import gzip
        from motif_io import stage_reference
        from pyfaidx import Fasta
        fasta = self.root / "reference ' $(touch NEVER).fa"
        fasta.write_text('>chr1\n' + 'ACGT' * 25 + '\n')
        with Fasta(str(fasta)):
            pass
        index = Path(str(fasta) + '.fai')
        original = index.read_bytes()
        zipped = self.root / 'reference.fa.gz'
        with gzip.open(zipped, 'wb') as stream:
            stream.write(fasta.read_bytes())
        target = self.root / 'stage'
        target.mkdir()
        staged = stage_reference(zipped, index, target)
        self.assertEqual(staged.read_bytes(), fasta.read_bytes())
        self.assertEqual(index.read_bytes(), original)
        with Fasta(str(staged), build_index=False) as reference:
            self.assertEqual(len(reference['chr1']), 100)


if __name__ == '__main__':
    unittest.main()
