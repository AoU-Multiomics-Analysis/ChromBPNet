import csv
import gzip
import hashlib
import io
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
sys.path.insert(0, str(ROOT / 'scripts'))


class GM12878DownloadTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def archive(self, name, entries):
        path = self.root / name
        with tarfile.open(path, 'w:gz') as archive:
            for member_name, data in entries:
                member = tarfile.TarInfo(member_name)
                member.size = len(data)
                archive.addfile(member, io.BytesIO(data))
        return path

    def inputs(self):
        models = self.archive('models.tar.gz', [
            (f'./fold_{fold}/model.chrombpnet_nobias.fold_{fold}.ENCSR637XSC.h5',
             b'\x89HDF\r\n\x1a\n' + bytes([fold])) for fold in range(5)
        ] + [('./fold_0/model.bias_scaled.fold_0.ENCSR637XSC.h5', b'bias')])
        regions = self.archive('regions.tar.gz', [
            ('./peaks.all_input_regions.ENCSR637XSC.bed.gz',
             gzip.compress(b'chr1\t10\t30\tp1\t10\t.\t2\t3\t4\t8\n')),
            ('./fold_0/nonpeaks.all_input_regions.fold_0.ENCSR637XSC.bed.gz', b'unused')])
        return models, regions

    def test_extracts_five_nobias_models_and_full_uncompressed_peaks(self):
        import prepare_gm12878 as tool
        models, regions = self.inputs()
        output = self.root / "output ' $(literal)"
        paths = tool.extract_inputs(models, regions, output)
        self.assertEqual(len(paths), 6)
        for fold in range(5):
            self.assertEqual((output / 'models' / f'fold_{fold}' / 'chrombpnet_nobias.h5').read_bytes(),
                             b'\x89HDF\r\n\x1a\n' + bytes([fold]))
        self.assertEqual((output / 'peaks' / 'GM12878_ATAC.narrowPeak').read_text(),
                         'chr1\t10\t30\tp1\t10\t.\t2\t3\t4\t8\n')
        self.assertFalse(list(output.rglob('*bias_scaled*')))

    def test_rejects_missing_fold_and_bad_peak_summit(self):
        import prepare_gm12878 as tool
        models, regions = self.inputs()
        incomplete = self.archive('missing.tar.gz', [])
        with self.assertRaisesRegex(ValueError, 'fold_0'):
            tool.extract_inputs(incomplete, regions, self.root / 'missing')
        invalid_regions = self.archive('bad-regions.tar.gz', [
            ('peaks.all_input_regions.ENCSR637XSC.bed.gz',
             gzip.compress(b'chr1\t10\t30\tp1\t10\t.\t2\t3\t4\t-1\n'))])
        with self.assertRaisesRegex(ValueError, 'summit'):
            tool.extract_inputs(models, invalid_regions, self.root / 'bad')

    def test_required_archive_member_cannot_be_a_link_or_duplicate(self):
        import prepare_gm12878 as tool
        models, regions = self.inputs()
        name = 'fold_0/model.chrombpnet_nobias.fold_0.ENCSR637XSC.h5'
        duplicate = self.archive('duplicate.tar.gz', [(name, b'a'), (name, b'b')])
        with self.assertRaisesRegex(ValueError, 'one regular file'):
            tool.extract_inputs(duplicate, regions, self.root / 'duplicate')
        link_path = self.root / 'link.tar.gz'
        with tarfile.open(link_path, 'w:gz') as archive:
            member = tarfile.TarInfo(name)
            member.type = tarfile.SYMTYPE
            member.linkname = '/etc/passwd'
            archive.addfile(member)
        with self.assertRaisesRegex(ValueError, 'one regular file'):
            tool.extract_inputs(link_path, regions, self.root / 'link')

    def test_download_checks_checksum_and_replaces_corrupt_cache(self):
        import prepare_gm12878 as tool
        source = self.root / 'source'
        source.write_bytes(b'correct bytes')
        destination = self.root / 'download'
        destination.write_bytes(b'corrupt cache')
        md5 = hashlib.md5(b'correct bytes').hexdigest()
        tool.download(source.as_uri(), destination, md5)
        self.assertEqual(destination.read_bytes(), b'correct bytes')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            tool.download(source.as_uri(), self.root / 'bad', '0' * 32)
        self.assertFalse((self.root / 'bad').exists())

    def test_preparation_copies_files_then_publishes_workflow_manifest(self):
        import prepare_gm12878 as tool
        from validate_manifest import read_manifest
        models, regions = self.inputs()
        sources = [(models.as_uri(), 'models.tar.gz', hashlib.md5(models.read_bytes()).hexdigest()),
                   (regions.as_uri(), 'regions.tar.gz', hashlib.md5(regions.read_bytes()).hexdigest())]
        cloud = self.root / 'bucket'
        prefix = 'gs://test-bucket/GM12878_ATAC'
        output = self.root / "local ' $(touch NEVER)"

        # Google authentication is external. Model only the gsutil boundary;
        # downloads, archive selection, file writes, and manifest parsing stay real.
        def copy_to_test_bucket(command, check):
            self.assertEqual(command[:2], ['gsutil', 'cp'])
            source, target = Path(command[2]), command[3]
            relative = target.removeprefix(prefix + '/')
            destination = cloud / relative
            if relative == 'models.tsv':
                for row in read_manifest(source):
                    for field in ['model', 'peaks']:
                        self.assertTrue((cloud / row[field].removeprefix(prefix + '/')).is_file())
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(source.read_bytes())

        with patch.object(tool, 'SOURCES', sources), patch.object(tool.subprocess, 'run', copy_to_test_bucket):
            tool.prepare(output, prefix, skip_reference=True)
        rows = read_manifest(output / 'models.tsv')
        self.assertEqual(len(rows), 5)
        self.assertEqual(rows[0]['model_id'], 'GM12878_ATAC_ENCSR637XSC_fold0')
        self.assertEqual(rows[0]['cell_type'], 'GM12878')
        self.assertEqual(rows[0]['model'], prefix + '/models/fold_0/chrombpnet_nobias.h5')
        self.assertEqual(rows[0]['peaks'], prefix + '/peaks/GM12878_ATAC.narrowPeak')
        self.assertEqual((cloud / 'models.tsv').read_bytes(), (output / 'models.tsv').read_bytes())
        self.assertFalse((output / 'NEVER').exists())
        for row in rows:
            self.assertTrue((cloud / row['model'].removeprefix(prefix + '/')).is_file())
        self.assertTrue((cloud / 'peaks' / 'GM12878_ATAC.narrowPeak').is_file())

    def test_reference_sizes_match_the_downloaded_fasta(self):
        import prepare_gm12878 as tool
        source = self.root / 'reference-source'
        source.mkdir()
        data = gzip.compress(b'>chr1\nAAAA\nCCC\n>chr2 description\nTT\n')
        (source / 'hg38.fa.gz').write_bytes(data)
        with patch.object(tool, 'UCSC', source.as_uri()), patch.object(tool, 'FASTA_MD5', hashlib.md5(data).hexdigest()):
            paths = tool.prepare_reference(self.root / 'reference-output')
        self.assertEqual(paths[0].read_bytes(), data)
        self.assertEqual(paths[1].read_text(), 'chr1\t7\nchr2\t2\n')

    def test_failed_upload_does_not_write_ready_manifest(self):
        import prepare_gm12878 as tool
        models, regions = self.inputs()
        sources = [(path.as_uri(), name, hashlib.md5(path.read_bytes()).hexdigest())
                   for path, name in [(models, 'models.tar.gz'), (regions, 'regions.tar.gz')]]
        output = self.root / 'failed-upload'
        with patch.object(tool, 'SOURCES', sources), patch.object(tool.subprocess, 'run', side_effect=RuntimeError('upload failed')):
            with self.assertRaisesRegex(RuntimeError, 'upload failed'):
                tool.prepare(output, 'gs://test-bucket/GM12878', skip_reference=True)
        self.assertFalse((output / 'models.tsv').exists())

    def test_rejects_cloud_wildcards_and_control_characters(self):
        import prepare_gm12878 as tool
        self.assertEqual(tool.validate_prefix('gs://test-bucket/GM12878/'), 'gs://test-bucket/GM12878')
        for prefix in ['http://bucket/x', 'gs://', 'gs://bucket/x*', 'gs://bucket/x\n', 'gs://bucket/x?generation=1']:
            with self.subTest(prefix=prefix), self.assertRaises(ValueError):
                tool.validate_prefix(prefix)


if __name__ == '__main__':
    unittest.main()
