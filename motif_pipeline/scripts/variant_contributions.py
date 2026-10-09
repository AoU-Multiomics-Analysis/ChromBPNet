"""Compute REF/ALT DeepSHAP on prepared sequences for one model fold."""
import argparse
import json
import shutil
from pathlib import Path

import h5py
import numpy as np

from motif_io import readable, sha256, create_contribution_datasets, write_json


def main():
    parser = argparse.ArgumentParser()
    for name in ['model', 'sequences', 'regions', 'preparation', 'cell-type', 'output-dir']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--fold-index', type=int, required=True)
    parser.add_argument('--head', choices=['counts', 'profile'], default='counts')
    parser.add_argument('--random-seed', type=int, default=1234)
    parser.add_argument('--num-backgrounds', type=int, default=20)
    parser.add_argument('--batch-size', type=int, default=64)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--allow-cpu', action='store_true', help='GitHub smoke test only')
    args = parser.parse_args()
    for key in ['model', 'sequences', 'regions', 'preparation']:
        setattr(args, key, readable(getattr(args, key)))
    meta = json.loads(args.preparation.read_text())
    if not 0 <= args.fold_index < 5 or meta['fold_count'] != 5 or sha256(args.model) != meta['model_sha256'][args.fold_index]:
        raise ValueError('Fold model does not match preparation')
    if meta['cell_type'] != args.cell_type or meta['head'] != args.head:
        raise ValueError('Cell type or head differs from preparation')
    if sha256(args.regions) != meta['region_sha256'] or sha256(args.sequences) != meta['sequences_sha256']:
        raise ValueError('Prepared sequence or row identity changed')
    if min(args.batch_size, args.num_backgrounds, args.threads) < 1:
        raise ValueError('Batch size, background count, and threads must be positive')
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
    from chrombpnet.evaluation.interpret import shap_utils
    from chrombpnet.training.utils import losses
    tf.keras.utils.get_custom_objects().update({'multinomial_nll': losses.multinomial_nll, 'tf': tf})
    model = tf.keras.models.load_model(str(args.model), compile=False)
    if model.input_shape != (None, meta['input_length'], 4) or len(model.outputs) != 2:
        raise ValueError('Use standard bias-corrected ChromBPNet models')
    # Reset the background RNG for each REF/ALT pair. Alleles share a seed;
    # their dinucleotide-shuffled sequences can differ because the alleles differ.
    def background(inputs):
        rng = np.random.RandomState(background.seed)
        return [np.stack([dinuc_shuffle(inputs[0], rng=rng) for _ in range(args.num_backgrounds)])]
    background.seed = args.random_seed
    target = (tf.reduce_sum(model.outputs[1], axis=-1) if args.head == 'counts'
              else shap_utils.get_weightedsum_meannormed_logits(model))
    explainer = shap.explainers.deep.TFDeepExplainer((model.input, target), background,
                combine_mult_and_diffref=shap_utils.combine_mult_and_diffref)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.sequences) as source, h5py.File(output / 'scores.h5', 'w') as scores:
        shape = source['raw/seq'].shape
        if shape != (meta['row_count'], 4, meta['input_length']) or shape[0] % 2:
            raise ValueError('Prepared allele row dimensions differ')
        create_contribution_datasets(scores, shape, args.batch_size)
        scores.attrs.update(cell_type=args.cell_type, head=args.head, fold_index=args.fold_index,
                            base_order='ACGT', region_sha256=meta['region_sha256'])
        for start in range(0, shape[0], args.batch_size):
            stop = min(shape[0], start + args.batch_size)
            raw = source['raw/seq'][start:stop]
            values = []
            for offset, sequence in enumerate(raw):
                background.seed = (args.random_seed + (start + offset) // 2) % (2 ** 32)
                value = explainer.shap_values(sequence.T[None].astype(np.float32), progress_message=100)
                if isinstance(value, list):
                    value = value[0]
                values.append(np.asarray(value)[0].T)
            hyp = np.asarray(values, dtype=np.float32)
            if hyp.shape != raw.shape or not np.isfinite(hyp).all():
                raise ValueError('Invalid hypothetical allele contributions')
            scores['raw/seq'][start:stop] = raw
            scores['shap/seq'][start:stop] = hyp
            scores['projected_shap/seq'][start:stop] = raw * hyp
            print(f'[variant contributions] fold={args.fold_index} rows={stop}/{shape[0]}', flush=True)
    shutil.copyfile(args.regions, output / 'regions.tsv')
    write_json(output / 'metadata.json', dict(meta, fold_index=args.fold_index,
               tensorflow=tf.__version__, num_backgrounds=args.num_backgrounds,
               random_seed=args.random_seed, gpu_count=len(gpus), allow_cpu=args.allow_cpu,
               chrombpnet_commit='09938fdb4397ec0006510e5251e48920a505d4de'))


if __name__ == '__main__':
    main()
