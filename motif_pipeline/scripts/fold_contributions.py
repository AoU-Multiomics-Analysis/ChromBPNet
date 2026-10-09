"""Compute hypothetical DeepSHAP scores with pinned ChromBPNet functions."""
import argparse
import csv
import importlib.metadata
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np

from motif_io import create_contribution_datasets, readable, sha256, write_json


def main():
    parser = argparse.ArgumentParser()
    for name in ['model', 'peaks', 'genome', 'genome-index', 'preparation', 'cell-type', 'output-dir']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--fold-index', type=int, required=True)
    parser.add_argument('--head', choices=['counts', 'profile'], default='counts')
    parser.add_argument('--random-seed', type=int, default=1234)
    parser.add_argument('--num-backgrounds', type=int, default=20)
    parser.add_argument('--batch-size', type=int, default=64)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--allow-cpu', action='store_true', help='For the GitHub image smoke test only')
    args = parser.parse_args()
    for key in ['model', 'peaks', 'genome', 'genome_index', 'preparation']:
        setattr(args, key, readable(getattr(args, key)))
    if min(args.num_backgrounds, args.batch_size, args.threads) < 1 or args.fold_index < 0:
        raise ValueError('Background count, batch size, and threads must be positive; fold index cannot be negative')
    preparation = json.loads(args.preparation.read_text())
    if preparation['cell_type'] != args.cell_type or preparation['head'] != args.head:
        raise ValueError('Cell type or head differs from the prepared inputs')
    if args.fold_index >= preparation['fold_count'] or sha256(args.model) != preparation['model_sha256'][args.fold_index]:
        raise ValueError('Fold model does not match the prepared input list')
    if sha256(args.peaks) != preparation['prepared_peak_sha256']:
        raise ValueError('Prepared peak identity/order changed')
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    # Stage an index next to a task-local symlink; do not modify localized files.
    reference = output.resolve() / 'reference.fa'
    reference.symlink_to(args.genome)
    shutil.copyfile(args.genome_index, str(reference) + '.fai')
    import tensorflow as tf
    tf.compat.v1.disable_eager_execution()
    tf.config.threading.set_intra_op_parallelism_threads(args.threads)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    gpus = tf.config.list_physical_devices('GPU')
    if not gpus and not args.allow_cpu:
        raise ValueError('GPU is required, but TensorFlow cannot detect one')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    import shap
    from deeplift.dinuc_shuffle import dinuc_shuffle
    from chrombpnet.evaluation.interpret import input_utils, shap_utils
    from pyfaidx import Fasta
    model = input_utils.load_model_wrapper(SimpleNamespace(model_h5=str(args.model)))
    input_length = preparation['input_length']
    if model.input_shape != (None, input_length, 4) or len(model.outputs) != 2:
        raise ValueError('Use a standard bias-corrected ChromBPNet model with profile and counts outputs')
    rng = np.random.RandomState(args.random_seed)

    def background(inputs):
        if len(inputs) != 1:
            raise ValueError('Model must have one sequence input')
        return [np.stack([dinuc_shuffle(inputs[0], rng=rng) for _ in range(args.num_backgrounds)])]

    target = (tf.reduce_sum(model.outputs[1], axis=-1) if args.head == 'counts'
              else shap_utils.get_weightedsum_meannormed_logits(model))
    explainer = shap.explainers.deep.TFDeepExplainer(
        (model.input, target), background,
        combine_mult_and_diffref=shap_utils.combine_mult_and_diffref)
    with args.peaks.open() as stream:
        peaks = list(csv.reader(stream, delimiter='\t'))
    shape = (len(peaks), 4, input_length)
    with Fasta(str(reference), build_index=False, rebuild=False) as genome, h5py.File(output / 'scores.h5', 'w') as handle:
        create_contribution_datasets(handle, shape, args.batch_size)
        handle.attrs.update(head=args.head, cell_type=args.cell_type, fold_index=args.fold_index,
                            random_seed=args.random_seed, num_backgrounds=args.num_backgrounds,
                            model_sha256=preparation['model_sha256'][args.fold_index],
                            region_sha256=preparation['prepared_peak_sha256'], base_order='ACGT')
        for start in range(0, len(peaks), args.batch_size):
            sequences = []
            for row in peaks[start:start + args.batch_size]:
                center = int(row[1]) + int(row[9])
                sequence = str(genome[row[0]][center - input_length // 2:center + input_length // 2])
                if len(sequence) != input_length:
                    raise ValueError('Prepared peak has an incomplete sequence window')
                sequences.append(sequence)
            encoded = input_utils.one_hot.dna_to_one_hot(sequences).astype(np.float32)
            scores = explainer.shap_values(encoded, progress_message=100)
            if isinstance(scores, list):
                if len(scores) != 1:
                    raise ValueError('Unexpected number of attribution outputs')
                scores = scores[0]
            scores = np.asarray(scores, dtype=np.float32)
            if scores.shape != encoded.shape or not np.isfinite(scores).all():
                raise ValueError('Hypothetical contributions have invalid dimensions or nonfinite values')
            stop = start + len(encoded)
            raw = encoded.transpose(0, 2, 1).astype(np.int8)
            hyp = scores.transpose(0, 2, 1)
            handle['raw/seq'][start:stop] = raw
            handle['shap/seq'][start:stop] = hyp
            handle['projected_shap/seq'][start:stop] = raw * hyp
            print(f'[contributions] fold={args.fold_index} head={args.head} regions={stop}/{len(peaks)}', flush=True)
    shutil.copyfile(args.peaks, output / 'regions.bed')
    write_json(output / 'metadata.json', dict(
        cell_type=args.cell_type, head=args.head, fold_index=args.fold_index,
        region_count=len(peaks), input_length=input_length, num_backgrounds=args.num_backgrounds,
        random_seed=args.random_seed, gpu_count=len(gpus), allow_cpu=args.allow_cpu,
        model_sha256=preparation['model_sha256'][args.fold_index],
        region_sha256=preparation['prepared_peak_sha256'], tensorflow=tf.__version__,
        chrombpnet_commit='09938fdb4397ec0006510e5251e48920a505d4de',
        kundajelab_shap=importlib.metadata.version('kundajelab-shap')))


if __name__ == '__main__':
    main()
