"""Exercise the real motif image on a planted motif and five synthetic models."""
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

import h5py
import numpy as np
import tensorflow as tf

SCRIPTS = Path('/opt/motif_pipeline/scripts')


def run(script, *args):
    subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], check=True)


def main():
    tf.config.threading.set_intra_op_parallelism_threads(2)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    os.environ['NUMBA_NUM_THREADS'] = '2'
    os.environ['OMP_NUM_THREADS'] = '2'
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        motif = 'ACGTCAGTGCAT'
        rng = np.random.default_rng(42)
        nregions, width = 128, 512
        regions = []
        peak_lines = []
        alphabet = np.array(list('ACGT'))
        for i in range(nregions):
            sequence = alphabet[rng.integers(0, 4, width)]
            for offset in [190, 290]:
                sequence[offset:offset + len(motif)] = list(motif)
            regions.append(''.join(sequence))
            center = i * width + width // 2
            peak_lines.append(f'chr1\t{center - 20}\t{center + 20}\tp{i}\t0\t.\t0\t0\t0\t20\n')
        fasta = root / 'genome.fa.gz'
        with gzip.open(fasta, 'wt') as stream:
            stream.write('>chr1\n' + ''.join(regions) + '\n')
        sizes = root / 'sizes.tsv'
        sizes.write_text(f'chr1\t{nregions * width}\n')
        peaks = root / "peaks ' $(touch NEVER).bed"
        peaks.write_text(''.join(peak_lines))
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
            model.get_layer('motif').set_weights([kernel, np.array([-len(motif) + 0.5], dtype=np.float32)])
            model.get_layer('counts').set_weights([np.array([[2000 + fold * 100]], dtype=np.float32), np.zeros(1, dtype=np.float32)])
            path = root / f"model ' $(touch NEVER) fold{fold}.h5"
            model.save(path, include_optimizer=False)
            models.append(path)
        model_list = root / 'models.list'
        model_list.write_text('\n'.join(map(str, models)) + '\n')
        prepared = root / 'prepared'
        run('prepare_motif_inputs.py', '--models-list', model_list, '--peaks', peaks,
            '--genome', fasta, '--chrom-sizes', sizes, '--cell-type', 'synthetic_CD4',
            '--input-length', width, '--output-dir', prepared)
        common = ['--model', models[0], '--peaks', prepared / 'peaks.bed',
                  '--genome', prepared / 'reference.fa', '--genome-index', prepared / 'reference.fa.fai',
                  '--preparation', prepared / 'preparation.json', '--cell-type', 'synthetic_CD4',
                  '--fold-index', '0', '--threads', '2', '--num-backgrounds', '2', '--batch-size', '16']
        denied = subprocess.run([sys.executable, str(SCRIPTS / 'fold_contributions.py'),
                                 *map(str, common), '--output-dir', str(root / 'denied')],
                                capture_output=True, text=True)
        assert denied.returncode != 0 and 'GPU is required' in denied.stderr
        scores, beds = [], []
        for fold in range(5):
            out = root / f'fold{fold}'
            args = common.copy()
            args[args.index('--model') + 1] = models[fold]
            args[args.index('--fold-index') + 1] = str(fold)
            run('fold_contributions.py', *args, '--allow-cpu', '--output-dir', out)
            scores.append(out / 'scores.h5')
            beds.append(out / 'regions.bed')
        # Exercise the separate profile graph with real DeepSHAP too.
        profile_prepared = root / 'prepared_profile'
        run('prepare_motif_inputs.py', '--models-list', model_list, '--peaks', peaks,
            '--genome', fasta, '--chrom-sizes', sizes, '--cell-type', 'synthetic_CD4',
            '--input-length', width, '--head', 'profile', '--max-peaks', '8',
            '--output-dir', profile_prepared)
        profile_args = common.copy()
        for flag, name in [('--peaks', 'peaks.bed'), ('--genome', 'reference.fa'),
                           ('--genome-index', 'reference.fa.fai'), ('--preparation', 'preparation.json')]:
            profile_args[profile_args.index(flag) + 1] = profile_prepared / name
        run('fold_contributions.py', *profile_args, '--head', 'profile', '--allow-cpu',
            '--output-dir', root / 'profile')
        with h5py.File(root / 'profile/scores.h5') as handle:
            assert handle.attrs['head'] == 'profile'
            values = handle['shap/seq'][:]
            assert values.shape == (8, 4, width) and np.isfinite(values).all()
            assert np.abs(values).max() > 0
        scores_list, beds_list = root / 'scores.list', root / 'beds.list'
        scores_list.write_text('\n'.join(map(str, scores)) + '\n')
        beds_list.write_text('\n'.join(map(str, beds)) + '\n')
        average = root / 'averaged.h5'
        run('average_contributions.py', '--scores-list', scores_list, '--regions-list', beds_list,
            '--output', average, '--output-regions', root / 'regions.bed', '--metadata', root / 'average.json')
        with h5py.File(average) as handle:
            fold_arrays = []
            for path in scores:
                with h5py.File(path) as fold_handle:
                    fold_arrays.append(fold_handle['shap/seq'][:])
            expected = np.mean(np.stack(fold_arrays), axis=0)
            np.testing.assert_allclose(handle['shap/seq'][:], expected, rtol=1e-5, atol=1e-7)
            assert handle['raw/seq'].shape == (nregions, 4, width)
        results = root / 'modisco_results.h5'
        run('discover_motifs.py', '--contributions', average, '--output', results,
            '--metadata', root / 'discovery.json', '--max-seqlets', '4096', '--n-leiden', '1')
        with h5py.File(results) as handle:
            patterns = [handle[group][name] for group in ['pos_patterns', 'neg_patterns']
                        if group in handle for name in handle[group]]
            assert patterns
            for pattern in patterns:
                assert all(key in pattern for key in ['sequence', 'contrib_scores', 'hypothetical_contribs', 'seqlets'])
        database = root / "known ' $(touch NEVER).meme"
        with database.open('w') as stream:
            stream.write('MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF SYNTHETIC Synthetic_TF\n')
            stream.write(f'letter-probability matrix: alength= 4 w= {len(motif)} nsites= 128\n')
            for base in motif:
                stream.write(' '.join('0.97' if item == base else '0.01' for item in 'ACGT') + '\n')
        report = root / 'annotation'
        run('report_motifs.py', '--motifs', results, '--motif-database', database, '--output-dir', report)
        with (report / 'candidate_tf_matches.tsv').open() as stream:
            matches = list(csv.DictReader(stream, delimiter='\t'))
        assert any(match['candidate_tf'] == 'Synthetic_TF' for match in matches), matches
        assert (report / 'report/report.html').stat().st_size > 0
        with zipfile.ZipFile(report / 'report_bundle.zip') as archive:
            assert archive.testzip() is None
            assert 'report/report.html' in archive.namelist()
        assert not (root / 'NEVER').exists()
        print('PASS: real five-fold DeepSHAP, signed averaging, motif discovery, Tomtom, portable report, and GPU requirement')


if __name__ == '__main__':
    main()
