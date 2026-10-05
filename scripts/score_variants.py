"""Run the pinned official scorer with localized inputs and strict output checks."""
import argparse
import csv
import json
import os
import runpy
import sys
from pathlib import Path
from io_utils import (CORE_METRICS, VARIANT_FIELDS, check_label, read_chrom_sizes,
                      read_variants, readable, validate_reference, write_scorer_variants, restore_score_file)


def main():
    parser = argparse.ArgumentParser()
    for name in ['variants', 'model', 'peaks', 'genome', 'genome-index', 'chrom-sizes', 'model-id', 'cell-type', 'output-dir']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--batch-size', type=int, default=128)
    parser.add_argument('--num-shuf', type=int, default=0)
    parser.add_argument('--max-peaks', type=int)
    parser.add_argument('--random-seed', type=int, default=1234)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--allow-cpu', action='store_true', help='For the CI smoke test only')
    args = parser.parse_args()
    check_label(args.model_id, args.cell_type)
    for name in ['variants', 'model', 'peaks', 'genome', 'genome_index', 'chrom_sizes']:
        setattr(args, name, readable(getattr(args, name)).resolve())
    if args.batch_size < 1 or args.num_shuf < 0 or args.threads < 1 or (args.max_peaks is not None and args.max_peaks < 1):
        raise ValueError('Batch size, thread count, and max peaks must be positive; num shuf must be nonnegative')
    variants = read_variants(args.variants)
    sizes = read_chrom_sizes(args.chrom_sizes)
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    local_genome = output / 'reference.fa'
    local_genome.symlink_to(args.genome)
    (output / 'reference.fa.fai').symlink_to(args.genome_index)
    os.environ.setdefault('TF_FORCE_GPU_ALLOW_GROWTH', 'true')
    import tensorflow as tf
    from pyfaidx import Fasta
    tf.config.threading.set_intra_op_parallelism_threads(args.threads)
    tf.config.threading.set_inter_op_parallelism_threads(1)
    gpus = tf.config.list_physical_devices('GPU')
    if not gpus and not args.allow_cpu:
        raise ValueError('GPU is required, but TensorFlow cannot detect one')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    vendor = Path(os.environ.get('VARIANT_SCORER_DIR', '/opt/variant-scorer/src'))
    sys.path.insert(0, str(vendor))
    from utils.helpers import load_model_wrapper
    with tf.device('/CPU:0'):
        model = load_model_wrapper(str(args.model))
    shape = model.input_shape
    if not isinstance(shape, tuple) or len(shape) != 3 or shape[-1] != 4:
        raise ValueError('Use a standard bias-corrected ChromBPNet model with one sequence input')
    input_length = int(shape[1])
    with Fasta(str(local_genome), build_index=False, rebuild=False) as genome:
        validate_reference(variants, genome, sizes, input_length)
        # Check all peak records before the upstream scorer selects valid sequence windows.
        with args.peaks.open() as stream:
            count = 0
            valid_count = 0
            for row in csv.reader(stream, delimiter='\t'):
                if not 3 <= len(row) <= 10 or row[0] not in sizes or row[0] not in genome:
                    raise ValueError('Peaks must be BED3 through narrowPeak with known chromosomes')
                chrom, start, end = row[0], int(row[1]), int(row[2])
                if not 0 <= start < end <= sizes[chrom] or len(genome[chrom]) != sizes[chrom]:
                    raise ValueError('Invalid peak coordinates or reference length')
                summit = int(row[9]) if len(row) == 10 else (end - start) // 2
                if not 0 <= summit < end - start:
                    raise ValueError('Peak summit must be within the peak; narrowPeak -1 is not supported')
                count += 1
                half = input_length // 2
                valid_count += start + summit - half > 0 and start + summit + half <= sizes[chrom]
            if not count or not valid_count:
                raise ValueError('No peak has a valid model sequence window')
    del model
    tf.keras.backend.clear_session()
    prefix = str(output / 'scored')
    scorer_variants = output / 'scorer_input.tsv'
    id_mapping = write_scorer_variants(variants, scorer_variants)
    sys.argv = [str(vendor / 'variant_scoring.py'), '--list', str(scorer_variants),
                '--model', str(args.model), '--peaks', str(args.peaks), '--genome', str(local_genome),
                '--chrom_sizes', str(args.chrom_sizes), '--out_prefix', prefix,
                '--schema', 'chrombpnet', '--batch_size', str(args.batch_size),
                '--num_shuf', str(args.num_shuf), '--random_seed', str(args.random_seed), '--no_hdf5']
    if args.max_peaks is not None:
        sys.argv += ['--max_peaks', str(args.max_peaks)]
    print(f'[score] Start model={args.model_id}, cell_type={args.cell_type}, GPUs={len(gpus)}', flush=True)
    runpy.run_path(str(vendor / 'variant_scoring.py'), run_name='__main__')
    restore_score_file(prefix + '.variant_scores.tsv', id_mapping, require_all=True)
    if args.num_shuf > 0:
        restore_score_file(prefix + '.variant_scores.shuffled.tsv', id_mapping)
    with open(prefix + '.variant_scores.tsv') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames
        scores = list(reader)
    if len(scores) != len(variants) or any(tuple(row[f] for f in VARIANT_FIELDS) != tuple(original[f] for f in VARIANT_FIELDS)
                                          for row, original in zip(scores, variants)):
        raise ValueError('Scorer lost variants or changed their coordinates/alleles; check chromosome names and windows')
    labelled = output / 'variant_effects.tsv'
    with labelled.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=['model_id', 'cell_type'] + fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(dict(row, model_id=args.model_id, cell_type=args.cell_type) for row in scores)
    from io_utils import read_scores
    read_scores(labelled)
    (output / 'run_metadata.json').write_text(json.dumps(dict(model_id=args.model_id, cell_type=args.cell_type,
        variant_count=len(scores), tensorflow=tf.__version__, gpu_count=len(gpus), input_length=input_length,
        num_shuf=args.num_shuf, max_peaks=args.max_peaks, random_seed=args.random_seed,
        variant_scorer_commit='0e1e34199e63112aa618748bb79a206fc491300a'), indent=2) + '\n')
    print(f'[score] Complete: {len(scores)} variants', flush=True)


if __name__ == '__main__':
    main()
