"""Run real five-fold REF/ALT DeepSHAP, indel preparation, and signed averaging."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import h5py
import numpy as np
import tensorflow as tf

SCRIPTS = Path('/opt/motif_pipeline/scripts')
sys.path.insert(0, str(SCRIPTS))
from motif_io import sha256


def run(script, *args):
    subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], check=True)


def main():
    tf.config.threading.set_intra_op_parallelism_threads(2)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    os.environ['OMP_NUM_THREADS'] = '2'
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        motif, width = 'ACGTACGA', 64
        sequence = 'C' * 60 + motif + 'C' * 60
        fasta = root / 'genome.fa'
        fasta.write_text('>chr1\n' + sequence * 3 + '\n')
        sizes = root / 'sizes.tsv'
        sizes.write_text('chr1\t384\n')
        variants = root / "variants ' $(touch NEVER).tsv"
        variants.write_text('chr1\t65\tA\tT\tsnv\nchr1\t193\tA\tAGG\tins\nchr1\t321\tAC\tA\tdel\n')
        models = []
        for fold in range(5):
            inputs = tf.keras.Input(shape=(width, 4))
            conv = tf.keras.layers.Conv1D(1, len(motif), activation='relu', name='motif')(inputs)
            profile = tf.keras.layers.Flatten(name='profile')(conv)
            counts = tf.keras.layers.Dense(1, name='counts')(tf.keras.layers.GlobalAveragePooling1D()(conv))
            model = tf.keras.Model(inputs, [profile, counts])
            kernel = np.zeros((len(motif), 4, 1), dtype=np.float32)
            for i, base in enumerate(motif):
                kernel[i, 'ACGT'.index(base), 0] = 1
            model.get_layer('motif').set_weights([kernel, np.array([-7.5], dtype=np.float32)])
            model.get_layer('counts').set_weights([np.array([[1000 + 100 * fold]], dtype=np.float32), np.zeros(1, dtype=np.float32)])
            path = root / f'fold{fold}.h5'
            model.save(path, include_optimizer=False)
            models.append(path)
        hashes = [sha256(p) for p in models]
        provenance = root / 'motif_provenance.json'
        provenance.write_text(json.dumps(dict(cell_type='synthetic_CD4', head='counts', model_sha256=hashes)))
        models_list = root / 'models.list'
        models_list.write_text('\n'.join(map(str, models)) + '\n')
        prepared = root / 'prepared'
        run('prepare_variant_motifs.py', '--models-list', models_list, '--variants', variants,
            '--genome', fasta, '--chrom-sizes', sizes, '--motif-provenance', provenance,
            '--cell-type', 'synthetic_CD4', '--input-length', width, '--output-dir', prepared)
        common = ['--sequences', prepared / 'shard_000000.h5', '--regions', prepared / 'shard_000000.tsv',
                  '--preparation', prepared / 'shard_000000.json', '--cell-type', 'synthetic_CD4',
                  '--num-backgrounds', '2', '--threads', '2', '--batch-size', '2']
        denied = subprocess.run([sys.executable, str(SCRIPTS / 'variant_contributions.py'), *map(str, common),
                                '--model', str(models[0]), '--fold-index', '0', '--output-dir', str(root/'denied')],
                               capture_output=True, text=True)
        assert denied.returncode != 0 and 'GPU is required' in denied.stderr
        files, rows = [], []
        for fold in range(5):
            out = root / f'scores{fold}'
            run('variant_contributions.py', *common, '--model', models[fold], '--fold-index', fold,
                '--allow-cpu', '--output-dir', out)
            files.append(out/'scores.h5')
            rows.append(out/'regions.tsv')
        a, b = root/'scores.list', root/'rows.list'
        a.write_text('\n'.join(map(str, files))+'\n')
        b.write_text('\n'.join(map(str, rows))+'\n')
        run('average_contributions.py', '--scores-list', a, '--regions-list', b, '--output', root/'mean.h5',
            '--output-regions', root/'rows.tsv', '--metadata', root/'mean.json')
        values = []
        for path in files:
            with h5py.File(path) as handle:
                values.append(handle['shap/seq'][:])
        with h5py.File(root/'mean.h5') as handle:
            np.testing.assert_allclose(handle['shap/seq'][:], np.mean(values, axis=0), rtol=1e-5, atol=1e-6)
            assert handle['raw/seq'].shape == (6, 4, width)
            projected = handle['projected_shap/seq'][:].sum(axis=(1,2))
            assert projected[0] > projected[1], projected
        profile_prepared = root/'profile_prepared'
        provenance.write_text(json.dumps(dict(cell_type='synthetic_CD4', head='profile', model_sha256=hashes)))
        run('prepare_variant_motifs.py', '--models-list', models_list, '--variants', variants,
            '--genome', fasta, '--chrom-sizes', sizes, '--motif-provenance', provenance, '--head', 'profile',
            '--cell-type', 'synthetic_CD4', '--input-length', width, '--output-dir', profile_prepared)
        run('variant_contributions.py', '--sequences', profile_prepared/'shard_000000.h5',
            '--regions', profile_prepared/'shard_000000.tsv', '--preparation', profile_prepared/'shard_000000.json',
            '--model', models[0], '--fold-index', '0', '--cell-type', 'synthetic_CD4', '--head', 'profile',
            '--num-backgrounds', '2', '--threads', '2', '--allow-cpu', '--output-dir', root/'profile_scores')
        with h5py.File(root/'profile_scores/scores.h5') as handle:
            assert np.isfinite(handle['shap/seq'][:]).all()
            assert handle.attrs['head'] == 'profile'
        assert not (root/'NEVER').exists()
        print('PASS: real five-fold REF/ALT DeepSHAP, SNV/indel coordinates, profile head, and signed averaging')


if __name__ == '__main__':
    main()
