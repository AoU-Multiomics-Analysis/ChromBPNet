"""CPU integration test for the Linux image; no real models or cloud jobs."""
import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
os.environ['TF_NUM_INTRAOP_THREADS'] = '1'
os.environ['TF_NUM_INTEROP_THREADS'] = '1'
import numpy as np
import tensorflow as tf
from pyfaidx import Fasta


def main():
    tf.keras.utils.set_random_seed(1234)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        genome = root / "reference ' $(literal).fa"
        genome.write_text('>chr1\n' + 'ACGT' * 500 + '\n')
        with Fasta(str(genome)):
            pass
        compressed_genome = root / "reference ' $(literal).fa.gz"
        with genome.open('rb') as source, gzip.open(compressed_genome, 'wb') as target:
            target.write(source.read())
        sizes = root / 'chrom.sizes'
        sizes.write_text('chr1\t2000\n')
        variants = root / 'variants.tsv'
        variants.write_text('chr1\t801\tA\tC\t001\nchr1\t1001\tA\tA\tnull\nchr1\t1201\tAC\tA\tdeletion\n')
        peaks = root / 'peaks.bed'
        peaks.write_text('chr1\t750\t850\nchr1\t950\t1050\n')
        seq = tf.keras.Input(shape=(100, 4))
        profile = tf.keras.layers.Conv1D(1, 1, use_bias=False, kernel_initializer=tf.keras.initializers.Constant([[[1.], [2.], [3.], [8.]]]))(seq)
        profile = tf.keras.layers.Reshape((100,))(profile)
        counts = tf.keras.layers.GlobalAveragePooling1D()(seq)
        counts = tf.keras.layers.Dense(1, kernel_initializer=tf.keras.initializers.Constant([[1.], [2.], [3.], [8.]]), bias_initializer='zeros')(counts)
        model = tf.keras.Model(seq, [profile, counts])
        model_path = root / 'chrombpnet_nobias.h5'
        model.save(model_path, include_optimizer=False)
        score_paths = []
        scripts = Path('/opt/chrombpnet/scripts')
        manifest = root / 'models.tsv'
        manifest.write_text('model_id\tcell_type\tmodel\tpeaks\n' +
            f'cd4\tCD4 T cell\t{model_path}\t{peaks}\n' + f'nk\tNK\t{model_path}\t{peaks}\n')
        subprocess.run([sys.executable, str(scripts / 'validate_manifest.py'), '--manifest', str(manifest),
            '--output', str(root / 'validated.txt'), '--rows-output', str(root / 'normalized.tsv')], check=True)
        assert (root / 'validated.txt').read_text() == '2\n'
        reference_rows = None
        for ident, cell in [('cd4', 'CD4 T cell'), ('nk', 'NK')]:
            output = root / ident
            command = [sys.executable, str(scripts / 'score_variants.py'), '--allow-cpu',
                '--variants', str(variants), '--model', str(model_path), '--peaks', str(peaks),
                '--genome', str(genome if ident == 'cd4' else compressed_genome), '--chrom-sizes', str(sizes),
                '--model-id', ident, '--cell-type', cell, '--output-dir', str(output), '--threads', '1', '--batch-size', '2']
            if ident == 'cd4':
                command += ['--genome-index', str(genome) + '.fai']
            # Exercise the optional peak and shuffle arguments on both models.
            command += ['--max-peaks', '2', '--num-shuf', '2']
            subprocess.run(command, check=True)
            path = output / 'variant_effects.tsv'
            score_paths.append(path)
            with path.open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            assert len(rows) == 3
            if reference_rows is not None:
                assert [row['variant_id'] for row in rows] == [row['variant_id'] for row in reference_rows]
                for field in ['logfc', 'jsd', 'active_allele_quantile']:
                    np.testing.assert_allclose([float(row[field]) for row in rows],
                                               [float(row[field]) for row in reference_rows])
            reference_rows = rows
            null = next(row for row in rows if row['variant_id'] == 'null')
            assert abs(float(null['logfc'])) < 1e-6
            assert abs(float(null['jsd'])) < 1e-6
            assert float(rows[0]['jsd']) > 0
            assert json.loads((output / 'run_metadata.json').read_text())['gpu_count'] == 0
        files = root / 'score_files.txt'
        files.write_text('\n'.join(map(str, score_paths)) + '\n')
        subprocess.run([sys.executable, str(scripts / 'merge_scores.py'), '--score-files', str(files),
            '--long-output', str(root / 'long.tsv'), '--wide-output', str(root / 'wide.tsv')], check=True)
        with (root / 'long.tsv').open() as stream:
            assert len(list(csv.DictReader(stream, delimiter='\t'))) == 6
        with (root / 'wide.tsv').open() as stream:
            assert len(list(csv.DictReader(stream, delimiter='\t'))) == 3
        # A normal run must reject CPU-only hosts.
        command.remove('--allow-cpu')
        command[command.index('--output-dir') + 1] = str(root / 'gpu_required')
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode != 0 and 'GPU is required' in result.stderr
        print('PASS: TSV manifest, plain and gzip FASTA, optional index, real scorer, merge, and GPU requirement')


if __name__ == '__main__':
    main()
