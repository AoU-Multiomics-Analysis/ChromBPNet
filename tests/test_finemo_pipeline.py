import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'motif_pipeline/scripts'))
sys.path.insert(0, str(ROOT / 'finemo_pipeline/scripts'))


class AlleleTests(unittest.TestCase):
    def test_snv_and_indel_sequences_and_reference_maps(self):
        from variant_sequences import allele_windows
        genome = {'chr1': 'ACGT' * 50}
        for ref, alt in [('A', 'T'), ('A', 'AGG'), ('ACG', 'A')]:
            row = dict(chr='chr1', pos='81', allele1=ref, allele2=alt, variant_id='v')
            sequences, maps, spans = allele_windows(row, genome, 32)
            self.assertEqual(sequences[0][16:16+len(ref)], ref)
            self.assertEqual(sequences[1][16:16+len(alt)], alt)
            self.assertTrue(all(len(s) == 32 for s in sequences))
            self.assertEqual(maps.shape, (2, 32))
            if len(alt) > len(ref):
                self.assertEqual(maps[1, 17:19].tolist(), [-1, -1])
                self.assertEqual(maps[1, 19], 81)
            if len(ref) > len(alt):
                self.assertEqual(maps[1, 17], 83)
            self.assertEqual(spans[1], (16, 16 + len(alt)))

    def test_ref_mismatch_edge_and_identity_rejected(self):
        from variant_sequences import allele_windows
        genome = {'chr1': 'ACGT' * 50}
        for pos, ref, alt in [('81', 'C', 'T'), ('1', 'A', 'T'), ('81', 'A', 'A')]:
            with self.subTest(pos=pos, ref=ref, alt=alt), self.assertRaises(ValueError):
                allele_windows(dict(chr='chr1', pos=pos, allele1=ref, allele2=alt, variant_id='v'), genome, 32)

    def test_unanchored_empty_alleles_rejected(self):
        from variant_sequences import allele_windows
        with self.assertRaisesRegex(ValueError, 'anchored'):
            allele_windows(dict(chr='chr1', pos='81', allele1='-', allele2='A', variant_id='v'),
                           {'chr1': 'ACGT' * 50}, 32)

    def test_preparation_shards_have_exact_rows_and_checksums(self):
        from prepare_variant_motifs import prepare_variants
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            variants = root / 'variants.tsv'
            variants.write_text('chr1\t81\tA\tT\tv1\nchr1\t85\tA\tAGG\tv2\n')
            prepare_variants(variants, {'chr1': 'ACGT' * 50}, 32, 1, root, 'CD4', 'counts', ['hash'] * 5)
            for shard in range(2):
                with h5py.File(root / f'shard_{shard:06d}.h5') as handle:
                    self.assertEqual(handle['raw/seq'].shape, (2, 4, 32))
                    self.assertEqual(handle['reference_positions'].shape, (2, 32))
                    rows = list(csv.reader((root / f'shard_{shard:06d}.tsv').open(), delimiter='\t'))
                    self.assertEqual([r[2] for r in rows], ['REF', 'ALT'])
                    self.assertEqual(handle.attrs['row_count'], len(rows))


class SummaryTests(unittest.TestCase):
    def rows(self):
        return [dict(region_index=str(i), variant_id='v' if i < 2 else 'empty', allele='REF' if i % 2 == 0 else 'ALT',
                     chr='chr1', pos='101', ref='A', alt='T', edit_start='16', edit_end='17',
                     window_start='84', insertion_anchor='100') for i in range(4)]

    def hit(self, row, motif='pos_patterns.pattern_0', start=14, value=2):
        return dict(peak_id=str(row), start=str(start), end=str(start+6), motif_name=motif,
                    strand='+', hit_coefficient_global=str(value), hit_importance=str(value), hit_similarity='0.9')

    def test_gain_loss_retained_signed_effect_and_empty_variants(self):
        from summarize_finemo import compare_variants
        hits = [self.hit(0), self.hit(1, value=1), self.hit(0, 'neg_patterns.pattern_0'),
                self.hit(1, 'pos_patterns.pattern_1')]
        maps = np.tile(np.arange(84, 116), (4, 1))
        calls, changes, variants = compare_variants(self.rows(), hits, maps, {}, 3)
        self.assertEqual({r['status'] for r in changes}, {'retained', 'gained', 'lost'})
        retained = next(r for r in changes if r['status'] == 'retained')
        self.assertEqual(retained['delta_hit_coefficient_global'], -1)
        self.assertEqual(retained['overlaps_edit'], True)
        self.assertEqual(len(variants), 2)
        self.assertEqual(next(r for r in variants if r['variant_id'] == 'empty')['n_changes'], 0)
        self.assertEqual(len(calls), 4)

    def test_indel_shift_pairs_by_reference_map_and_preserves_inserted_bases(self):
        from summarize_finemo import compare_variants
        rows = self.rows()[:2]
        rows[1]['edit_end'] = '19'
        maps = np.tile(np.arange(84, 116), (2, 1))
        maps[1] = np.array(list(range(84, 101)) + [-1, -1] + list(range(101, 114)))
        calls, changes, _ = compare_variants(rows, [self.hit(0, start=19), self.hit(1, start=21)], maps, {}, 0)
        self.assertEqual(changes[0]['status'], 'retained')
        calls, _, _ = compare_variants(rows, [self.hit(1, start=16)], maps, {}, 3)
        self.assertEqual(calls[0]['inserted_base_count'], 2)
        self.assertEqual(calls[0]['insertion_anchor'], 100)

    def test_multiple_tf_matches_do_not_multiply_changes(self):
        from summarize_finemo import compare_variants
        matches = {'pos_patterns.pattern_0': ['TF1', 'TF2']}
        _, changes, _ = compare_variants(self.rows(), [self.hit(1)], np.tile(np.arange(32), (4, 1)), matches, 3)
        self.assertEqual(len(changes), 1)
        self.assertEqual(changes[0]['candidate_tfs'], 'TF1;TF2')

    def test_invalid_hit_or_missing_allele_fails(self):
        from summarize_finemo import compare_variants
        maps = np.tile(np.arange(32), (4, 1))
        for rows, hits in [(self.rows()[:1], []), (self.rows(), [self.hit(9)]),
                           (self.rows(), [self.hit(0, value=float('nan'))])]:
            with self.assertRaises(ValueError):
                compare_variants(rows, hits, maps, {}, 3)


class FinemoInputTests(unittest.TestCase):
    def test_existing_tfmodisco_match_schema(self):
        from finemo_io import tf_matches
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'matches.tsv'
            path.write_text('pattern\ttarget_motif\tcandidate_tf\tq_value\n'
                            'pos_patterns.pattern_0\tm1\tTF1\t0.01\n'
                            'pos_patterns.pattern_0\tm2\tTF2\t0.02\n')
            self.assertEqual(tf_matches(path)['pos_patterns.pattern_0'], ['TF1', 'TF2'])
    def test_sequence_and_nonfinite_validation(self):
        from prepare_finemo import validate_arrays
        raw = np.zeros((2, 4, 32), dtype=np.int8)
        raw[:, 0] = 1
        hyp = raw.astype(np.float32)
        validate_arrays(raw, hyp)
        for bad_raw, bad_hyp in [(raw * 2, hyp), (raw, hyp * np.nan), (raw, hyp[:, :, :-1])]:
            with self.assertRaises(ValueError):
                validate_arrays(bad_raw, bad_hyp)

    def test_unresolved_uri_rejected(self):
        from finemo_io import readable
        with self.assertRaisesRegex(ValueError, 'localization'):
            readable('gs://bucket/file.h5')


if __name__ == '__main__':
    unittest.main()
