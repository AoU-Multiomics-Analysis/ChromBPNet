import csv
import hashlib
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
sys.path.insert(0, str(ROOT / 'scripts'))


class TenK10KDownloadTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def archive(self, entries):
        path = self.root / 'peaks.zip'
        with zipfile.ZipFile(path, 'w') as archive:
            for name, content in entries:
                archive.writestr(name, content)
        return path

    def fixture(self):
        folder = 'TenK10K_ATAC_ct_specific_peaks/'
        archive = self.archive([
            (folder + 'ASDC.csv', ',chrom,start,end,name,score,strand,signal_value,p_value,q_value,peak\n'
             '0,chr1,100,180,.,15,.,2.5,3.5,4.5,17\n'),
            (folder + 'CD14_Mono.csv', 'chr1\t200\t280\t.\t0\t.\n'),
            (folder + 'ILC.csv', ',chrom,start,end,name,score,strand,signal_value,p_value,q_value,peak\n'
             '0,chr2,300,400,.,20,.,2,3,4,30\n'),
            ('__MACOSX/' + folder + '._ASDC.csv', b'ignored'),
        ])
        models = []
        for cell in ['ASDC', 'CD14_Mono']:
            for fold in range(5):
                path = self.root / f'{cell}_fold{fold}.h5'
                path.write_bytes(b'\x89HDF\r\n\x1a\n' + cell.encode() + bytes([fold]))
                models.append((f'ChromBPNet/{cell}/chrombpnet_model_full_bias_NKFull1_fold{fold}/models/chrombpnet_nobias.h5', path))
        combined = self.root / 'combined.bed'
        combined.write_text('chr1\t100\t180\tannotation\n')
        paths = models + [
            ('Miscellaneous/TenK10K_ATAC_ct_specific_peaks.zip', archive),
            ('Miscellaneous/TenK10K_ATAC_MACS3_Combined_Peaks_Annotated.bed', combined),
        ]
        catalog = {'repository': 'anglixue/TenK10K_multiome', 'revision': 'a' * 40,
                   'peak_cell_types': ['ASDC', 'CD14_Mono', 'ILC'],
                   'model_cell_types': ['ASDC', 'CD14_Mono'],
                   'files': [{'path': name, 'size': path.stat().st_size,
                              'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for name, path in paths]}
        index = self.root / 'index.json'
        index.write_text(json.dumps(catalog))
        urls = {name: path.as_uri() for name, path in paths}
        return index, urls, archive

    def test_converts_indexed_csv_and_bed6_without_losing_summits(self):
        # Catch a shifted CSV index, lost summit, or wrong cell-to-peak pairing.
        import prepare_tenk10k as tool
        _, _, archive = self.fixture()
        report = tool.extract_peaks(archive, self.root / 'out', ['ASDC', 'CD14_Mono', 'ILC'])
        self.assertEqual((self.root / 'out/peaks/ASDC.narrowPeak').read_text(),
                         'chr1\t100\t180\t.\t15\t.\t2.5\t3.5\t4.5\t17\n')
        self.assertEqual((self.root / 'out/peaks/CD14_Mono.bed').read_text(),
                         'chr1\t200\t280\t.\t0\t.\n')
        self.assertEqual(report['ASDC']['peak_count'], 1)
        self.assertEqual(report['CD14_Mono']['centering'], 'midpoint')
        self.assertFalse(list((self.root / 'out').rglob('._*')))

    def test_missing_fold_duplicate_fold_and_wrong_model_fail(self):
        # Catch silent partial downloads or selection of the full bias model.
        import prepare_tenk10k as tool
        index, _, _ = self.fixture()
        catalog = json.loads(index.read_text())
        with self.assertRaisesRegex(ValueError, 'fold'):
            tool.select_assets({**catalog, 'files': catalog['files'][1:]})
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            tool.select_assets({**catalog, 'files': catalog['files'] + [catalog['files'][0]]})
        bad = json.loads(index.read_text())
        bad['files'][0]['path'] = bad['files'][0]['path'].replace('chrombpnet_nobias', 'chrombpnet')
        with self.assertRaisesRegex(ValueError, 'fold'):
            tool.select_assets(bad)

    def test_missing_entire_model_cell_type_is_not_silently_skipped(self):
        import prepare_tenk10k as tool
        index, _, _ = self.fixture()
        catalog = json.loads(index.read_text())
        catalog['files'] = [entry for entry in catalog['files'] if '/ASDC/' not in entry['path']]
        with self.assertRaisesRegex(ValueError, 'missing.*ASDC'):
            tool.select_assets(catalog)

    def test_t7_output_requires_mounted_drive_before_creating_directories(self):
        import prepare_tenk10k as tool
        index, _, _ = self.fixture()
        with patch.object(Path, 'is_mount', return_value=False):
            with self.assertRaisesRegex(ValueError, 'T7.*mounted'):
                tool.prepare('/Volumes/T7/SHOULD_NOT_CREATE', 'gs://test-bucket/tenk10k',
                             skip_reference=True, source_index=index, plan_only=True)

    def test_rejects_bad_summit_missing_cell_and_unsafe_archive(self):
        import prepare_tenk10k as tool
        folder = 'TenK10K_ATAC_ct_specific_peaks/'
        header = 'chrom,start,end,name,score,strand,signal_value,p_value,q_value,peak\n'
        bad = self.archive([(folder + 'ASDC.csv', header + 'chr1,10,30,.,0,.,1,2,3,-1\n')])
        with self.assertRaisesRegex(ValueError, 'summit'):
            tool.extract_peaks(bad, self.root / 'bad', ['ASDC'])
        good = self.archive([(folder + 'ASDC.csv', 'chr1\t10\t30\n')])
        with self.assertRaisesRegex(ValueError, 'missing'):
            tool.extract_peaks(good, self.root / 'missing', ['ASDC', 'ILC'])
        unsafe = self.archive([('../ASDC.csv', 'chr1\t10\t30\n')])
        with self.assertRaisesRegex(ValueError, 'Unsafe'):
            tool.extract_peaks(unsafe, self.root / 'unsafe', ['ASDC'])

    def test_checksum_cache_is_verified_before_reuse(self):
        import prepare_tenk10k as tool
        source = self.root / 'source'
        source.write_bytes(b'verified')
        target = self.root / 'target'
        target.write_bytes(b'corrupt')
        checksum = hashlib.sha256(b'verified').hexdigest()
        tool.download(source.as_uri(), target, checksum, 8)
        self.assertEqual(target.read_bytes(), b'verified')
        source.unlink()
        tool.download(source.as_uri(), target, checksum, 8)
        with self.assertRaisesRegex(ValueError, 'checksum'):
            tool.download(target.as_uri(), self.root / 'bad', '0' * 64, 8)
        self.assertFalse((self.root / 'bad').exists())

    def test_manifest_upload_is_last_and_all_cloud_targets_exist(self):
        # Catch premature manifest publication and cloud paths that omit cell/fold.
        import prepare_tenk10k as tool
        from validate_manifest import read_manifest
        index, urls, _ = self.fixture()
        output = self.root / "local ' $(touch NEVER)"
        bucket = self.root / 'bucket'
        prefix = 'gs://test-bucket/tenk10k'
        uploads = []

        def upload(command, check):
            self.assertEqual(command[:2], ['gsutil', 'cp'])
            source = Path(command[2])
            relative = command[3].removeprefix(prefix + '/')
            if relative == 'models.tsv':
                for row in read_manifest(source):
                    for column in ['model', 'peaks']:
                        self.assertTrue((bucket / row[column].removeprefix(prefix + '/')).is_file())
            target = bucket / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read_bytes())
            uploads.append(relative)

        with patch.object(tool, 'asset_url', side_effect=lambda catalog, path: urls[path]), \
                patch.object(tool.subprocess, 'run', side_effect=upload):
            tool.prepare(output, prefix, skip_reference=True, source_index=index)
        rows = read_manifest(output / 'models.tsv')
        self.assertEqual(len(rows), 10)
        self.assertEqual(rows[0], {'model_id': 'TenK10K_ATAC_ASDC_fold0', 'cell_type': 'ASDC',
                                  'model': prefix + '/models/ASDC/fold_0/chrombpnet_nobias.h5',
                                  'peaks': prefix + '/peaks/ASDC.narrowPeak'})
        self.assertEqual(rows[5]['peaks'], prefix + '/peaks/CD14_Mono.bed')
        self.assertEqual(uploads[-1], 'models.tsv')
        self.assertTrue((bucket / 'peaks/ILC.narrowPeak').is_file())
        report = json.loads((output / 'sources.json').read_text())
        self.assertEqual(report['peak_cell_types_without_models'], ['ILC'])
        self.assertFalse((output / 'NEVER').exists())

    def test_failed_upload_does_not_create_ready_manifest(self):
        import prepare_tenk10k as tool
        index, urls, _ = self.fixture()
        output = self.root / 'failed'
        with patch.object(tool, 'asset_url', side_effect=lambda catalog, path: urls[path]), \
                patch.object(tool.subprocess, 'run', side_effect=RuntimeError('upload failed')):
            with self.assertRaisesRegex(RuntimeError, 'upload failed'):
                tool.prepare(output, 'gs://test-bucket/tenk10k', skip_reference=True, source_index=index)
        self.assertFalse((output / 'models.tsv').exists())

    def test_plan_needs_no_download_or_upload_and_is_not_ready_manifest(self):
        import prepare_tenk10k as tool
        from validate_manifest import read_manifest
        index, _, _ = self.fixture()
        output = self.root / 'plan'
        with patch.object(tool, 'download', side_effect=AssertionError('download during plan')), \
                patch.object(tool.subprocess, 'run', side_effect=AssertionError('upload during plan')):
            tool.prepare(output, 'gs://test-bucket/tenk10k', skip_reference=True,
                         source_index=index, plan_only=True)
        self.assertEqual(len(read_manifest(output / 'models.planned.tsv')), 10)
        self.assertFalse((output / 'models.tsv').exists())


if __name__ == '__main__':
    unittest.main()
